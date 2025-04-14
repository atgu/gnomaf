__author__ = 'Lindo Nkambule & Sophie Parsa'

import argparse
import hail as hl

from hail.expr.functions import _cdf_combine, _quantile_from_cdf, _result_from_raw_cdf

AS_INFO_AGG_FIELDS_VQSR = {
    "sum_agg_fields": ["AS_QUALapprox", "AS_VarDP"],
    "median_agg_fields": ["AS_ReadPosRankSum_cdf", "AS_MQRankSum_cdf"],
    "array_sum_agg_fields": ["AS_SB_TABLE"],
    "array_extend_agg_fields": ["AS_MQ"],
}


def combine_two_cohorts(ht1: hl.Table, ht2: hl.Table, cdf_k: int = 200) -> hl.Table:
    """
    Given AS annotations from two cohorts, this function combines the annotations into one

    :param ht1: input Table from cohort 1 with AS annotations from compute_info_intermediate
    :param ht2: input Table from cohort 2 with AS annotations from compute_info_intermediate
    :param cdf_k: Parameter controlling the accuracy vs. memory usage tradeoff when retaining CDF
    :return: hl.Table with annotations from two cohorts combined
    """
    joined_ht = ht1.join(ht2, how='outer')

    # for AS_MQ, we need to keep the values in array with each element representing a cohort's AS_MQ
    # when AS_MQ is missing in a cohort, fill with 0
    joined_ht = joined_ht.annotate(AS_MQ=hl.if_else(hl.is_missing(joined_ht.AS_MQ),
                                                    hl.array([hl.float64(0)]), joined_ht.AS_MQ),
                                   AS_MQ_1=hl.if_else(hl.is_missing(joined_ht.AS_MQ_1),
                                                      hl.array([hl.float64(0)]), joined_ht.AS_MQ_1))

    annotations = {}
    fields_to_drop = []

    # each field from "left" cohort will be overwritten by the combined field
    for k, v in AS_INFO_AGG_FIELDS_VQSR.items():
        for field in v:
            fields_to_drop.append(f"{field}_1")
            ann1 = joined_ht[f"{field}"]
            ann2 = joined_ht[f"{field}_1"]

            if k == "median_agg_fields":
                combine_expr = hl.array([
                    _cdf_combine(k=cdf_k, left=ann1[0], right=ann2[0])
                ])
            elif k == "sum_agg_fields":
                combine_expr = hl.array([
                    ann1[0] + ann2[0]
                ])
            elif k == "array_sum_agg_fields":
                combine_expr = hl.map(lambda x: x[0] + x[1], hl.zip(ann1, ann2))

            else:  # collecting the annotation across cohorts into an array
                combine_expr = ann1.extend(ann2)

            annotations[field] = hl.if_else(
                hl.is_missing(ann1),
                ann2,
                hl.if_else(
                    hl.is_missing(ann2),
                    ann1,
                    combine_expr
                )
            )

    joined_ht = joined_ht.annotate(**annotations)

    # drop fields from the "right" cohort
    joined_ht = joined_ht.drop(*fields_to_drop)

    return joined_ht


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--intermediate-annotations', type=str, required=True)
    parser.add_argument('--out-dir', type=str, required=True)
    parser.add_argument('--output-prefix', type=str, required=True)
    parser.add_argument('--cdf_k', type=int, default=200)
    parser.add_argument('--target_table_size', type=int, default=200)
    parser.add_argument('--sample-sizes', type=str, required=False)
    # --sample-sizes order MUST correspond to --intermediate-annotations order !!!

    args = parser.parse_args()
    cohorts = args.intermediate_annotations.split(",")
    as_hts = [hl.read_table(f'{ht_p}') for ht_p in cohorts]

    # hierarchically combine annotations
    combined_ht = as_hts[0]
    for ht in as_hts[1:]:
        combined_ht = combine_two_cohorts(ht1=combined_ht, ht2=ht, cdf_k=args.cdf_k)

    # compute final (aggregated) annotations
    # 1. AS_ReadPosRankSum
    combined_ht = combined_ht.annotate(
        AS_ReadPosRankSum=hl.array([_quantile_from_cdf(_result_from_raw_cdf(combined_ht.AS_ReadPosRankSum_cdf[0]), 0.5)])
    )

    # 2. AS_MQRankSum
    combined_ht = combined_ht.annotate(
        AS_MQRankSum=hl.array([_quantile_from_cdf(_result_from_raw_cdf(combined_ht.AS_MQRankSum_cdf[0]), 0.5)])
    )

    # 3. AS_QD
    combined_ht = combined_ht.annotate(
        AS_QD=combined_ht.AS_QUALapprox / combined_ht.AS_VarDP
    )

    # 4. AS_SOR
    # for AS_SOR, we need to convert to float to avoid int32 overflow for large values in AS_SB_TABLE
    combined_ht = combined_ht.annotate(
        AS_SB_TABLE_tmp=hl.map(lambda AS_SB_TABLE: hl.map(lambda x: hl.float64(x + 1), AS_SB_TABLE), combined_ht.AS_SB_TABLE)
    )

    combined_ht = combined_ht.annotate(
        symmetricalRatio=((combined_ht.AS_SB_TABLE_tmp[0][0] * combined_ht.AS_SB_TABLE_tmp[1][1]) / (combined_ht.AS_SB_TABLE_tmp[0][1] * combined_ht.AS_SB_TABLE_tmp[1][0])) + ((combined_ht.AS_SB_TABLE_tmp[0][1] * combined_ht.AS_SB_TABLE_tmp[1][0]) / (combined_ht.AS_SB_TABLE_tmp[0][0] * combined_ht.AS_SB_TABLE_tmp[1][1])),
        refRatio=hl.min(combined_ht.AS_SB_TABLE_tmp[0][0], combined_ht.AS_SB_TABLE_tmp[0][1]) / hl.max(combined_ht.AS_SB_TABLE_tmp[0][0], combined_ht.AS_SB_TABLE_tmp[0][1]),
        altRatio=hl.min(combined_ht.AS_SB_TABLE_tmp[1][0], combined_ht.AS_SB_TABLE_tmp[1][1]) / hl.max(combined_ht.AS_SB_TABLE_tmp[1][0], combined_ht.AS_SB_TABLE_tmp[1][1])
    )

    combined_ht = combined_ht.annotate(
        AS_SOR=hl.array([hl.log(combined_ht.symmetricalRatio) + hl.log(combined_ht.refRatio) - hl.log(combined_ht.altRatio)])
    )

    # 5. AS_FS
    # For AS_FS, we need to normalize the SB_TABLE counts if their sum is greater than 2*TARGET_TABLE_SIZE
    target_table_size = args.target_table_size
    combined_ht = combined_ht.annotate(
        AS_SB_TABLE_sum=hl.sum(combined_ht.AS_SB_TABLE[0]) + hl.sum(combined_ht.AS_SB_TABLE[1])
    )

    combined_ht = combined_ht.annotate(
        AS_SB_TABLE_FS=hl.if_else(combined_ht.AS_SB_TABLE_sum > target_table_size * 2,
                                  hl.map(lambda AS_SB_TABLE: hl.map(lambda x: hl.int32(x / (combined_ht.AS_SB_TABLE_sum / target_table_size)), AS_SB_TABLE),combined_ht.AS_SB_TABLE),
                                  combined_ht.AS_SB_TABLE)
    )

    combined_ht = combined_ht.annotate(
        AS_FS=-10 * hl.log10(hl.fisher_exact_test(combined_ht.AS_SB_TABLE_FS[0][0],
                                                  combined_ht.AS_SB_TABLE_FS[0][1],
                                                  combined_ht.AS_SB_TABLE_FS[1][0],
                                                  combined_ht.AS_SB_TABLE_FS[1][1]).p_value
                             )
    )

    # https://hail.zulipchat.com/#narrow/channel/123010-Hail-Query-0.2E2-support/topic/.E2.9C.94.20Negative.20values.20after.20Phred-scaling
    combined_ht = combined_ht.annotate(
        AS_FS=hl.array([hl.if_else(combined_ht.AS_FS > 0, combined_ht.AS_FS, 0)])
    ) # to avoid -0.0 values

    # 6. AS_MQ
    # if sample size are given, use weighted mean, else use median
    if args.sample_sizes:
        sample_sizes = args.sample_sizes.split(",")
        sample_sizes = [int(i) for i in sample_sizes]
        assert len(sample_sizes) == len(as_hts), "--sample-sizes must have the same items as --intermediate-annotations"
        weights_expr = hl.literal(sample_sizes)

        combined_ht = combined_ht.annotate(
            AS_MQ=hl.sum(combined_ht.AS_MQ * weights_expr) / hl.sum(weights_expr)
        )
    else:
        combined_ht = combined_ht.annotate(
            AS_MQ=hl.median(combined_ht.AS_MQ)
        )

    # drop intermediate annotations
    fields_to_keep = {'locus', 'alleles', 'AS_MQ', 'AS_ReadPosRankSum', 'AS_MQRankSum', 'AS_QD', 'AS_SOR', 'AS_FS'}
    combined_ht = combined_ht.drop(*[k for k in combined_ht.row if k not in fields_to_keep])

    # nest annotations under info
    combined_ht = combined_ht.annotate(
        info=hl.struct(**{
            k: combined_ht[k] for k in combined_ht.row if k not in {'locus', 'alleles'}
        })
    )
    combined_ht = combined_ht.drop(*[k for k in combined_ht.row if k not in {'locus', 'alleles', 'info'}])

    # write HT with AS annotations
    combined_ht.write(f'{args.out_dir}/{args.output_prefix}_as_annotations_final.ht', overwrite=True)


if __name__ == '__main__':
    main()