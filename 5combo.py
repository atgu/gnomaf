import hail as hl
import os
import gnomad.utils.sparse_mt as m
import numpy as np

from gnomad.utils.annotations import (
    fs_from_sb,
    sor_from_sb,
)

from hail.vds.combiner.combine import merge_alleles

hl.init(log=hl.utils.timestamp_path(os.path.join('/tmp', 'hail'), suffix='.log'), default_reference='GRCh38', gcs_requester_pays_configuration='diverse-pop-seq-ref')

num_datasets = 5

# Fields that need to be combined as array element wise sums
array_sum_fields = ['AS_QUALapprox','AS_VarDP','AC','AS_ReadPosRankSum_n_smaller','AS_ReadPosRankSum_n_larger','AS_MQRankSum_n_smaller','AS_MQRankSum_n_larger','AS_pab_max','AS_RAW_MQ', 'AS_MQ_DP']

# Fields that need to be combined as element wise array of array sums
array_array_sum_fields = ['AS_SB_TABLE','AS_ReadPosRankSum_freq','AS_MQRankSum_freq']

# Prepares hail table to be combined with other tables
def clean_up_input_tables(ht):
    # Update bin edges
    ht = ht.annotate_globals(bin_edges = [hl.float64(i) for i in np.arange(-5, 5.01, .01)])

    # Explode the RankSum structs into their individual components
    ht = ht.annotate(info = ht.info.annotate(AS_ReadPosRankSum_freq = ht.info.AS_ReadPosRankSum.bin_freq,
                                     AS_ReadPosRankSum_n_smaller = ht.info.AS_ReadPosRankSum.n_smaller,
                                     AS_ReadPosRankSum_n_larger = ht.info.AS_ReadPosRankSum.n_larger))

    ht = ht.annotate(info = ht.info.annotate(AS_MQRankSum_freq = ht.info.AS_MQRankSum.bin_freq,
                                     AS_MQRankSum_n_smaller = ht.info.AS_MQRankSum.n_smaller,
                                     AS_MQRankSum_n_larger = ht.info.AS_MQRankSum.n_larger))

    ht = ht.annotate(info = ht.info.annotate(AS_RAW_MQ = ht.info.AS_RAW_MQandDP.map(lambda x: x[0]),
                                     AS_MQ_DP = ht.info.AS_RAW_MQandDP.map(lambda x: x[1])))

    # Clean up the info field by dropping unnecessary fields
    ht = ht.annotate(info = ht.info.drop("AS_ReadPosRankSum", "AS_MQRankSum", "AS_SB", "AS_QD", "AS_SOR","AS_FS","AC_raw","AS_RAW_MQandDP"))
    ht = ht.drop('AS_lowqual')

    # Key by locus as preparation for table join
    ht = ht.key_by('locus')

    return ht

# Adjusts AS annotation arrays in each dataset to be consistent with the list of merged alleles
def pad_array(alt_alleles, field, dataset_alleles):
    """
    param alt_alleles: list of merged alleles between all datasets excluding reference allele
    param field: an info field representing a particular AS annotation
    param dataset_alleles: the alleles found in a particular dataset
    """
    return alt_alleles.map(lambda x: hl.or_else(field[dataset_alleles.index(x)-1], 0))

# Adjusts AS annotation array<array> in each dataset to be consistent with the list of merged alleles
def pad_array_array(alleles, field, dataset_alleles, impute_len, is_sb_table):
    """
    param alt_alleles: list of merged alleles between all datasets excluding reference allele
    param field: an info field representing a particular AS annotation
    param dataset_alleles: the alleles found in a particular dataset
    param impute_len: length that the array should be in the case it is missing (2 for SB_TABLE, num_bins for RankSum)
    """
    if(is_sb_table):
        return alleles.map(lambda x: hl.or_else(field[dataset_alleles.index(x)], hl.literal([0 for i in range(impute_len)])))
    else:
        return alleles.map(lambda x: hl.or_else(field[dataset_alleles.index(x)-1], hl.literal([0 for i in range(impute_len)])))


# Combines arrays element wise from a list of arrays
def combine_arrays(to_be_combined):
    base_array = to_be_combined[0]
    for i in range(1,num_datasets):
        base_array += to_be_combined[i]
    return base_array


# Finds element wise max in a list of arrays
def max_arrays(to_be_maxed):
    new_arr = to_be_maxed[0].map(lambda x: hl.max(hl.array([to_be_maxed[i][to_be_maxed[0].index(x)] for i in range(len(to_be_maxed))])))
    return new_arr               
                     

# Calculates the median based on a histogram of bins, frequencies, n_smaller, and n_larger
def find_median(freqs, bins, n_smaller, n_larger):
    # Add n_smaller and n_larger as extra bins
    freqs = hl.array([n_smaller]).extend(freqs).append(n_larger)
    # Add -3.5 and +3.5 as new bin min and max
    bins = hl.literal([hl.float64(-5.5)]).extend(bins).append(hl.float64(5.5))
    # Calculate cumulative sum for each bin and identify the midpoint bin
    cumulative_sums = hl.cumulative_sum(freqs)
    total_sum = hl.sum(freqs)
    median_index = cumulative_sums.index(lambda x: x >= total_sum / 2)
    # Calculate median according to several cases
    median = hl.if_else(
        (total_sum % 2 == 0) & (cumulative_sums[median_index] == total_sum / 2),
        (bins[median_index + 1] + bins[freqs[median_index+1:].index(lambda x: x > 0) + median_index + 1])/2,
        (bins[median_index] + bins[median_index + 1]) / 2
    )
    median = hl.if_else(total_sum == 0, hl.missing(hl.tfloat), median)
    return median

def main():
    # Read in hail tables
    nigeria_ht = hl.read_table('gs://gnomaf/AS_annotations/nigeria_54gene_merged_gvcfs.ht')
    ggv_ht = hl.read_table('gs://gnomaf/AS_annotations/gambian_genomes_merged_gvcfs.ht')
    sa_ht = hl.read_table('gs://gnomaf/AS_annotations/neurogap_south_africa_genomes_merged_gvcfs.ht')
    neuro_extra_ht = hl.read_table('gs://gnomaf/AS_annotations/neurogap_highcov_all_sites_93_genomes_merged_gvcfs.ht')
    hgdp_ht = hl.read_table('gs://gnomaf/AS_annotations/gnomad.genomes.v3.1.2.hgdp_1kg_subset.ht')

    # Clean up the hail tables
    nigeria_ht = clean_up_input_tables(nigeria_ht)
    ggv_ht = clean_up_input_tables(ggv_ht)
    sa_ht = clean_up_input_tables(sa_ht)
    neuro_extra_ht = clean_up_input_tables(neuro_extra_ht)
    hgdp_ht = clean_up_input_tables(hgdp_ht)

    # Outer join all AS annotation tables
    jt = nigeria_ht.join(ggv_ht, how="outer")
    jt = jt.join(sa_ht, how = "outer")
    jt = jt.join(neuro_extra_ht, how = "outer")
    jt = jt.join(hgdp_ht, how = "outer")

    # Use merge_alleles to merge alleles from all datasets: alleles stores the merged alleles array and alleles_og stores a list of the original alleles in each dataset
    jt = jt.annotate(alleles_set_1=jt.alleles)
    jt = jt.annotate(alleles_set_2=jt.alleles_1)
    jt = jt.annotate(alleles_set_3=jt.alleles_2)
    jt = jt.annotate(alleles_set_4=jt.alleles_3)
    jt = jt.annotate(alleles_set_5=jt.alleles_4)
    jt = jt.annotate(alleles=merge_alleles(hl.array([jt.alleles_set_1, jt.alleles_set_2, jt.alleles_set_3, jt.alleles_set_4, jt.alleles_set_5])).globl)
    jt = jt.annotate(alleles_og=merge_alleles(hl.array([jt.alleles_set_1, jt.alleles_set_2, jt.alleles_set_3, jt.alleles_set_4, jt.alleles_set_5])).local)

    # Get the array sum fields from each dataset
    dic_set_1 = {f: jt.info[f] for f in array_sum_fields if f in jt.info}
    dic_set_2 = {f: jt.info_1[f] for f in array_sum_fields if f in jt.info_1}
    dic_set_3 = {f: jt.info_2[f] for f in array_sum_fields if f in jt.info_2}
    dic_set_4 = {f: jt.info_3[f] for f in array_sum_fields if f in jt.info_3}
    dic_set_5 = {f: jt.info_4[f] for f in array_sum_fields if f in jt.info_4}

    # Pad the array sum fields from each dataset
    jt = jt.annotate(**{field_name+'_set_1': pad_array(jt.alleles[1:],jt.info[field_name], jt.alleles_og[0]) for field_name, field in dic_set_1.items()})
    jt = jt.annotate(**{field_name+'_set_2': pad_array(jt.alleles[1:],jt.info_1[field_name], jt.alleles_og[1]) for field_name, field in dic_set_2.items()})
    jt = jt.annotate(**{field_name+'_set_3': pad_array(jt.alleles[1:],jt.info_2[field_name], jt.alleles_og[2]) for field_name, field in dic_set_3.items()})
    jt = jt.annotate(**{field_name+'_set_4': pad_array(jt.alleles[1:],jt.info_3[field_name], jt.alleles_og[3]) for field_name, field in dic_set_4.items()})
    jt = jt.annotate(**{field_name+'_set_5': pad_array(jt.alleles[1:],jt.info_4[field_name], jt.alleles_og[4]) for field_name, field in dic_set_5.items()})

    # Combine the padded arrays into a list
    dic = {f: [jt[f+'_set_'+str(i)] for i in range(1,num_datasets+1)] for f in array_sum_fields}

    # Combine array AS annotations
    jt = jt.annotate(**{field_name: max_arrays(value) if field_name == 'AS_pab_max' else combine_arrays(value) for field_name, value in dic.items()})
    

    # Get the array<array> sum fields from each dataset
    dic_set_1_aa = {f: jt.info[f] for f in array_array_sum_fields if f in jt.info}
    dic_set_2_aa = {f: jt.info_1[f] for f in array_array_sum_fields if f in jt.info_1}
    dic_set_3_aa = {f: jt.info_2[f] for f in array_array_sum_fields if f in jt.info_2}
    dic_set_4_aa = {f: jt.info_3[f] for f in array_array_sum_fields if f in jt.info_3}
    dic_set_5_aa = {f: jt.info_4[f] for f in array_array_sum_fields if f in jt.info_4}

    # Pad the array<array> fields in each dataset
    jt = jt.annotate(**{field_name+'_set_1': pad_array_array(jt.alleles if field_name == 'AS_SB_TABLE' else jt.alleles[1:],jt.info[field_name], jt.alleles_og[0], 2 if field_name == 'AS_SB_TABLE' else 1000, field_name == 'AS_SB_TABLE') for field_name, field in dic_set_1_aa.items()})
    jt = jt.annotate(**{field_name+'_set_2': pad_array_array(jt.alleles if field_name == 'AS_SB_TABLE' else jt.alleles[1:],jt.info_1[field_name], jt.alleles_og[1], 2 if field_name == 'AS_SB_TABLE' else 1000, field_name == 'AS_SB_TABLE') for field_name, field in dic_set_2_aa.items()})
    jt = jt.annotate(**{field_name+'_set_3': pad_array_array(jt.alleles if field_name == 'AS_SB_TABLE' else jt.alleles[1:],jt.info_2[field_name], jt.alleles_og[2], 2 if field_name == 'AS_SB_TABLE' else 1000, field_name == 'AS_SB_TABLE') for field_name, field in dic_set_3_aa.items()})
    jt = jt.annotate(**{field_name+'_set_4': pad_array_array(jt.alleles if field_name == 'AS_SB_TABLE' else jt.alleles[1:],jt.info_3[field_name], jt.alleles_og[3], 2 if field_name == 'AS_SB_TABLE' else 1000, field_name == 'AS_SB_TABLE') for field_name, field in dic_set_4_aa.items()})
    jt = jt.annotate(**{field_name+'_set_5': pad_array_array(jt.alleles if field_name == 'AS_SB_TABLE' else jt.alleles[1:],jt.info_4[field_name], jt.alleles_og[4], 2 if field_name == 'AS_SB_TABLE' else 1000, field_name == 'AS_SB_TABLE') for field_name, field in dic_set_5_aa.items()})

    # Combine the array<arrays> into a list
    dic_array_array = {f: [jt[f+'_set_'+str(i)] for i in range(1,num_datasets+1)] for f in array_array_sum_fields}

    # Combine array<array> annotations
    jt = jt.annotate(**{field_name: hl.zip(*value).map(lambda x: combine_arrays(list(x))) for field_name, value in dic_array_array.items()})

    # Calculate medians, QD, and MQ using the other annotations
    jt = jt.annotate(AS_QD = hl.map(lambda x, y: hl.if_else(y == 0, 0, x/y), jt.AS_QUALapprox, jt.AS_VarDP))

    jt = jt.annotate(AS_MQ = hl.map(lambda x, y: hl.if_else(y == 0, 0, hl.sqrt(x/y)), jt.AS_RAW_MQ, jt.AS_MQ_DP))

    jt = jt.annotate(AS_ReadPosRankSum = hl.map(lambda x, y, z: find_median(x, jt.bin_edges, y, z), jt.AS_ReadPosRankSum_freq, jt.AS_ReadPosRankSum_n_smaller, jt.AS_ReadPosRankSum_n_larger))

    jt = jt.annotate(AS_MQRankSum = hl.map(lambda x, y, z: find_median(x, jt.bin_edges, y, z), jt.AS_MQRankSum_freq, jt.AS_MQRankSum_n_smaller, jt.AS_MQRankSum_n_larger))


    # Use gnomad functions to calculate FS and SOR
    jt = jt.annotate(AS_FS = jt.AS_SB_TABLE[1:].map(lambda x: hl.if_else(hl.is_defined(x),fs_from_sb(hl.array([jt.AS_SB_TABLE[0], x])),hl.missing('float64'))))
    jt = jt.annotate(AS_SOR = jt.AS_SB_TABLE[1:].map(lambda x: hl.if_else(hl.is_defined(x),sor_from_sb(hl.array([jt.AS_SB_TABLE[0], x])),hl.missing('float64'))))

    # Re key table by allele and locus
    jt = jt.key_by('locus','alleles')


    # Export the final combined table
    condensed = jt.annotate(info = hl.struct(
        AS_QUALapprox = jt.AS_QUALapprox, 
        AS_VarDP = jt.AS_VarDP, 
        AC = jt.AC, 
        AS_pab_max = jt.AS_pab_max, 
        AS_SB_TABLE = jt.AS_SB_TABLE, 
        AS_QD = jt.AS_QD, 
        AS_ReadPosRankSum = jt.AS_ReadPosRankSum, 
        AS_ReadPosRankSum_n_smaller = jt.AS_ReadPosRankSum_n_smaller,
        AS_ReadPosRankSum_n_larger = jt.AS_ReadPosRankSum_n_larger,
        AS_ReadPosRankSum_freq = jt.AS_ReadPosRankSum_freq,
        AS_MQRankSum = jt.AS_MQRankSum,
        AS_MQRankSum_n_smaller = jt.AS_MQRankSum_n_smaller,
        AS_MQRankSum_n_larger = jt.AS_MQRankSum_n_larger,
        AS_MQRankSum_freq = jt.AS_MQRankSum_freq,
        AS_FS = jt.AS_FS, 
        AS_SOR = jt.AS_SOR, 
        AS_MQ = jt.AS_MQ,
        AS_RAW_MQ = jt.AS_RAW_MQ,
        AS_MQ_DP = jt.AS_MQ_DP))
    



    combo = condensed.select(condensed.info)
    combo.write('gs://sophie-bucket/gnomaf_combined_as_annotations.ht')

if __name__ == '__main__':
    main()













