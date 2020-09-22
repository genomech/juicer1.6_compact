# ExoC

# What is sam_to_valid_pairs.sh ?

This sript created for processing Hi-C data. It supports all types of Hi-C protocols and only requires a .SAM input file.
This script resut in ValidPairs file. ValidPairs format is :


```
    chromosome1 position1 length1 strand1 chromosome2 position2 length2 strand2
```

Besides of ValidPairs file this script get some key features of Hi-C library quality:

- fracton of multiple aligment in all pairs
- fracton of one side unmapped in all pairs
- fracton of both side unmapped in all pairs
- fracton of unique alignment in all pairs. (both reads alignes uniquely)
- fracton of PCR duplicates in unique alignments
- fracton of dangling end in unique alignments
- fracton of valid pairs in all pairs
- fracton of cis contacts in valid pairs
- fracton of trans contacts in valid pairs
- fracton of cis short range (distance less than 20kb) contacts in valid pairs
- fracton of cis long range (distance more than 20kb) contacts in valid pairs

# Run sam_to_valid_pairs.sh

```
    ./sam_to_valid_pairs.sh sam_file out_dir
```


