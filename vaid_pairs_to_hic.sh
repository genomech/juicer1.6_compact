
PATH_TO_JUICER_TOOLS=$1
GENOME_VERSION=$2
VALIDPAIRS_FILE=$3
OUT_DIR=$4
SAMPLE_NAME=$5

cd ${OUT_DIR}



time cat ${VALIDPAIRS_FILE} | awk '{ { print "0" " chr" $1 " " $2 " " "0" " " "1" " chr" $3 " " $4 " " "1"  } }' > allValidPairs${SAMPLE_NAME}_pre_hic
time java -Xmx2g -jar ${PATH_TO_JUICER_TOOLS} pre -r 5000000,2500000,1000000,500000,100000,50000,25000,10000,5000,1000 allValidPairs${SAMPLE_NAME}_pre_hic s${SAMPLE_NAME}.hic ${GENOME_VERSION}
rm allValidPairs${SAMPLE_NAME}_pre_hic

exit
