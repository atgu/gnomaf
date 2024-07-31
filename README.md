# gnomaf

# n + k annotations
The script 5combo.py provides code to combine dataset allele specific annotations in an n+k manner

The process:

1. Merge gVCFS into one VDS
2. Convert a VDS into an AS annotation table (using Lindo's scripts to fill in missing annotations for each dataset)
3. Combine AS tables using the n + k combination code --> note: the RankSum combination is not currently exact, would need to do an update to the code with the new histogram CDF in hail 
