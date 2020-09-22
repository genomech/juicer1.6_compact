# ExoC

# What is sam_to_valid_pairs.sh ?

## Valid pairs

This sript created for processing Hi-C data. It supports all types of Hi-C protocols and only requires a .SAM input file.
This script resut in ValidPairs file. ValidPairs format is :


```
    chromosome1 position1 length1 strand1 chromosome2 position2 length2 strand2
```
## Satistics

Besides of ValidPairs file this script get some key features of Hi-C library quality:

- fracton of unique alignment in all pairs. (both reads alignes uniquely)
- fracton of multiple aligment in all pairs
- fracton of one side unmapped in all pairs
- fracton of both side unmapped in all pairs
- fracton of PCR duplicates in unique alignments
- fracton of dangling ends in unique alignments
- fracton of valid pairs in all pairs
- fracton of cis contacts in valid pairs
- fracton of trans contacts in valid pairs
- fracton of cis short range (distance less than 20kb) contacts in valid pairs
- fracton of cis long range (distance more than 20kb) contacts in valid pairs

## that is dangling ends ?

Dangling ends is a whole genome fragment. So it is not a ligation product. Dangling ends can appear in the Hi-C data in two ways: because of back ligations and from anligated genome fragments. Anyway this reads do not carry information about 3D conformation of chromatin.

# Run sam_to_valid_pairs.sh

Before running sam_to_valid_pairs.sh highly recommended to use BWA-mem for mapping Hi-C data.

To run sam_to_valid_pairs.sh use:

```
    ./sam_to_valid_pairs.sh sam_file out_dir
```

# Get hic file from valid pairs

If you want to get .hic file from valid pairs run:

```
    ./vaid_pairs_to_hic.sh path_to_juicer_tools validpairs_file out_dir
```
