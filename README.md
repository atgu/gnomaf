# gnomAF

# n + k annotations
We assume that a VDS (Hail Combiner from gVCFS) for each cohort has been generated. The as_annotations/ directory contains two scripts for merging AS annotations
1.[compute_info_intermediate.p](https://github.com/atgu/gnomaf/blob/main/as_annotations/compute_info_intermediate.py): run on each cohort VDS to generate per-cohort HailTable containing AS annotations. This script will also try to compute any missing annotations required downstream
2.[compute_info_final.py](https://github.com/atgu/gnomaf/blob/main/as_annotations/compute_info_final.py): takes the tables from the first script and combines them into one table

The idea: when we get a new cohort, the only thing we need to do is (1) generate AS annotations for that cohort using [compute_info_intermediate.p](https://github.com/atgu/gnomaf/blob/main/as_annotations/compute_info_intermediate.py), and then; (2) combine the annotations from this cohort with other already existing cohorts using [compute_info_final.py](https://github.com/atgu/gnomaf/blob/main/as_annotations/compute_info_final.py). This [document](https://docs.google.com/document/d/1iohLHNJOIV2n9Dr2AGnCNpCxi7jEGpvRSQyq-xXrGao/edit?tab=t.0) contains information about which annotations are used and how they are merged.
