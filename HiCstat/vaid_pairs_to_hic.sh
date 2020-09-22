
PATH_TO_JUICER_TOOLS=$1
GENOME_VERSION=$3
VALIDPAIRS_FILE=$4
OUT_DIR=$5

cd ${OUT_DIR}

time cat allValidPairs${SAMPLE} | awk '{ { print "0" " chr" $2 " " $3 " " "0" " " "1" " chr" $6 " " $7 " " "1"  } }' > allValidPairs${SAMPLE}_pre_hic
time java -Xmx2g -jar ${PATH_TO_JUICER_TOOLS} pre -r 5000000,2500000,1000000,500000,100000,50000,25000,10000,5000,1000 ${VALIDPAIRS_FILE} out.hic ${GENOME_VERSION}
rm allValidPairs${SAMPLE}_pre_hic

exit
