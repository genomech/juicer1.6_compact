# ExoC

## What is sam_to_valid_pairs.sh ?

### valid pairs

This sript created for processing Hi-C data. It supports all types of Hi-C protocols and only requires a .SAM input file.
This script resut in ValidPairs file. ValidPairs format is :


```
    read_name chromosome1 position1 length1 strand1 chromosome2 position2 length2 strand2
```
### satistics

Besides of ValidPairs file this script get some key features of quality of Hi-C libraries:

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

### what is dangling ends ?

Dangling ends is a whole genome fragment. Dangling ends can appear in the Hi-C data in two ways: because of back ligations and from anligated genome fragments. Anyway this reads do not carry information about 3D conformation of chromatin.

## Run sam_to_valid_pairs.sh

Before running sam_to_valid_pairs.sh highly recommended to use BWA-mem for mapping Hi-C data.

To run sam_to_valid_pairs.sh use:

```
./sam_to_valid_pairs.sh sam_file out_dir
```

## Get hic file from valid pairs

If you want to get .hic file from valid pairs run:

```
./vaid_pairs_to_hic.sh path_to_juicer_tools genome_version validpairs_file out_dir
```

For example, genome_version = "hg19" or "hg38"

## Optional requirements

- BWA-mem alinger https://github.com/lh3/bwa 
- Juicer Tools https://github.com/aidenlab/juicer/wiki/Juicer-Tools-Quick-Start
