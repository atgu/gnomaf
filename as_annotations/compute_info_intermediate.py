__author__ = 'Lindo Nkambule & Sophie Parsa'

import argparse
import hail as hl
import logging

from gnomad.utils.sparse_mt import default_compute_info
from typing import List

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def check_missing_annotations(fields: List[str]):
    """
    Check missing annotations that are required to generate AS annotations

    :param fields: list of fields present in the dataset
    :return: list of missing annotations
    """

    if 'RAW_MQandDP' in fields:
        required_fields = ['QUALapprox', 'VarDP', 'ReadPosRankSum', 'MQRankSum', 'SB', 'RAW_MQandDP']
    else:
        required_fields = ['QUALapprox', 'VarDP', 'ReadPosRankSum', 'MQRankSum', 'SB', 'MQ_DP', 'RAW_MQ']

    return list(set(required_fields).difference(fields))


def compute_missing_annotations(vds: hl.vds.variant_dataset.VariantDataset,
                                missing_annotations: List[str]) -> hl.vds.variant_dataset.VariantDataset:
    """
    Compute missing annotations

    :param vds: input VDS with missing annotations
    :param missing_annotations: list of missing annotations
    :return: VariantDataset
    """

    if ('VarDP' in missing_annotations) and ('LAD' in vds.variant_data.entry):
        logger.info("Computing `VarDP` as sum(LAD)")
        vds.variant_data = vds.variant_data.annotate_entries(gvcf_info=
                                                             vds.variant_data.gvcf_info.annotate(VarDP=
                                                                                                 hl.int32(hl.sum(
                                                                                                     vds.variant_data.LAD))))

    if ('QUALapprox' in missing_annotations) and ('LPL' in vds.variant_data.entry):
        logger.info("Computing `QUALapprox` as LPL[0]")
        vds.variant_data = vds.variant_data.annotate_entries(gvcf_info=
                                                             vds.variant_data.gvcf_info.annotate(QUALapprox=
                                                                                                 vds.variant_data.LPL[
                                                                                                     0]))
    if ('MQ_DP' in missing_annotations) and ('DP' in vds.variant_data.entry):
        logger.info("Computing `MQ_DP` as DP")
        vds.variant_data = vds.variant_data.annotate_entries(gvcf_info=
                                                             vds.variant_data.gvcf_info.annotate(MQ_DP=
                                                                                                 vds.variant_data.DP))

    return vds


def filter_format_ht(ht: hl.Table) -> hl.Table:
    """
    Format the input Table so it can be aggregated with other tables. These are the steps we apply:
        (1) split multi-allelics into bi-allelics and update the allele annotations
        (2) drop a variant (bi-allelic) if it's missing one of the RankSum annotations.
        For the other annotations, we set the annotation to 0 if it's missing
        (3) keep only AS_QUALapprox, AS_VarDP, AS_ReadPosRankSum_cdf, AS_MQRankSum_cdf, AS_MQ, AS_SB_TABLE as these
        are the only annotations we need to compute the final AS annotations

    :param ht: input Table with AS annotations from gnomAD's default_compute_info()
    :return: hl.Table
    """
    # Remove the nested structure in info (it will be added back in the end). Also drop fields we will not use
    ht = ht.annotate(**ht.info)
    ht = ht.select(*[f for f in ht.row_value if f != 'info'])
    fields_to_drop = ['AS_lowqual', 'AS_pab_max', 'AC_raw', 'AC', 'AS_QD', 'AS_FS', 'AS_SOR']
    ht = ht.drop(*fields_to_drop)

    # split multi-allelics
    bi = ht.filter(hl.len(ht.alleles) == 2)
    bi = bi.annotate(a_index=1, was_split=False, old_locus=bi.locus, old_alleles=bi.alleles)
    multi = ht.filter(hl.len(ht.alleles) > 2)
    split = hl.split_multi(multi)
    ht = split.union(bi)

    # split_multi will only split the alleles, not annotations. Update the annotations
    ht = ht.annotate(a_index=hl.array([ht.a_index])) # a_index starts at 1
    annotations = {}
    for field in ["AS_ReadPosRankSum", "AS_MQRankSum", "AS_QUALapprox", "AS_VarDP", "AS_ReadPosRankSum_cdf",
                  "AS_MQRankSum_cdf", "AS_MQ"]:
        annotations[field] = ht.a_index.map(lambda i: ht[field][i - 1])
    ht = ht.annotate(**annotations)

    ht = ht.annotate(
        AS_SB_TABLE=hl.array([ht.AS_SB_TABLE[0]]).extend(ht.a_index.map(lambda i: ht.AS_SB_TABLE[i]))
        # first element is REF allele
    )

    fields_to_drop = ['a_index', 'was_split', 'old_locus', 'old_alleles']
    ht = ht.drop(*fields_to_drop)

    # Filter out alleles missing either AS_ReadPosRankSum or AS_MQRankSum
    ht = ht.filter((hl.is_missing(ht.AS_ReadPosRankSum[0]) | hl.is_missing(ht.AS_MQRankSum[0])), keep=False)
    ht = ht.drop(*['AS_ReadPosRankSum', 'AS_MQRankSum'])

    return ht


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--inputs', type=str, required=True)
    parser.add_argument('--out-dir', type=str, required=True)
    parser.add_argument('--cdf_k', type=int, default=200)

    args = parser.parse_args()
    paths = args.inputs.split(",")

    for filepath in paths:
        print(f"---------- Processing {filepath} ----------")
        outname = f'{filepath.split("/")[-1].split(".")[0]}_as_annotations_intermediate'
        input_type = filepath.split("/")[-1].split(".")[1]

        if input_type == "vds":
            vds = hl.vds.read_vds(filepath)
            missing_anns = check_missing_annotations(list(vds.variant_data.gvcf_info) + list(vds.variant_data.entry))
            if len(missing_anns) > 0:
                vds = compute_missing_annotations(vds, missing_anns)

            # data to be passed to default_compute_info should be filtered to nonref sites
            mt = vds.variant_data.filter_entries(vds.variant_data.LGT.is_non_ref(), keep=True)

        else: # HGDP+1kGP file is a MT and has no missing annotations
            mt = hl.read_matrix_table(filepath)

            # data to be passed to default_compute_info should be filtered to nonref sites
            mt = mt.filter_entries(mt.LGT.is_non_ref(), keep=True)

            # We have RAW_MQ and MQ_DP as separate annotations but default_compute_info looks for RAW_MQandDP
            mt = mt.annotate_entries(gvcf_info=mt.gvcf_info.annotate(RAW_MQandDP=hl.array([mt.gvcf_info.RAW_MQ,
                                                                                           mt.gvcf_info.MQ_DP])))

        # compute AS annotations
        as_ht = default_compute_info(mt, retain_cdfs=True, cdf_k=args.cdf_k)
        # Filter out rows that have NA for all AS annotations
        # https://hail.zulipchat.com/#narrow/channel/123010-Hail-Query-0.2E2-support/topic/.E2.9C.94.20Filter.20out.20rows.20with.20NAs.20in.20a.20HT
        as_ht = as_ht.filter(hl.array([hl.is_defined(as_ht.info[f]) & as_ht.info[f].any(lambda x: hl.is_defined(x)) for f in as_ht.info]).any(lambda a: a))

        as_ht = filter_format_ht(as_ht)

        # write HT with AS annotations
        as_ht.write(f'{args.out_dir}/{outname}.ht', overwrite=True)


if __name__ == '__main__':
    main()
