#!/bin/bash


## The script will result in the file : allValidPairs
# allValidPairs format: read_name chromosome1 position1 length1 strand1 chromosome2 position2 length2 strand2

# allValidPairs format:  chromosome1 position1 chromosome2 position2 strand1 strand2 read_name

SAM_FILE=$1
OUT_DIR=$2

cd ${OUT_DIR}



############################################
## Sam parsing
############################################
samtools view -O SAM ${SAM_FILE} | awk '{
if (NR>94) {
  print $0
}
}' > s1

samtools view -O SAM ${SAM_FILE} | grep -v "^@" > s1
###@@@@@@@@@@@@@@@@@@@@@@@@@@@@STAT
cat s1 | awk 'BEGIN{k=0;n=0}{
  if (n!=$1) {
    n=$1
    k=k+1
}
}END{ { print k } }' > stat




#extract nessesary columns and descript orientation from SAM flag 

cat s1 | awk '
function d2b5(d,  b) {
      k = 0
      while(k != 5) {
          k = k + 1
          b=d%2
          d=int(d/2)
      }
      return(b)
}
function d2b7(d,  b) {
      k = 0
      while(k != 7) {
          k = k + 1
          b=d%2
          d=int(d/2)
      }
      return(b)
}
{
print $1 " " d2b5($2) " " $3 " " $4 " " $5 " " d2b7($2)


}' > s2
########################################################################################

########################################################################################


############################################
## Identificate reported pairs, multiple alignment, unmapped
############################################
cat s1 | awk 'BEGIN{k=0; n=0; rep=0; mult=0; unmap=0; rep_count=0; mult_count=0; unmap_one_count=0; unmap_both_count=0;}{
if (n!=$1) {
  n=$1
  
  if(rep>1){
  rep_count=rep_count+1
  } else{
    if(mult>0){
      mult_count=mult_count+1
    } else{
      if(rep==1){
        unmap_one_count=unmap_one_count+1
      } else {
        unmap_both_count=unmap_both_count+1
      }
    }
  }
  
  rep=0;mult=0;unmap_one=0;unmap_both=0;
}

if($5>0 && length($3)<=5 && $3!="*"){
  rep=rep+1
}
if($5==0 && length($3)<=5 && $3!="*"){
  mult=mult+1
}
}END{
print rep_count " " mult_count " " unmap_one_count " " unmap_both_count-1
}' > is_rep_mult_unmap


rm s1
###@@@@@@@@@@@@@@@@@@@@@@@@@@@@STAT
#is_rep, mult, unmap_both, unmap_one 
cat is_rep_mult_unmap | awk '{ { print $1 } }'  >> stat
cat is_rep_mult_unmap | awk '{ { print $2 } }'  >> stat
cat is_rep_mult_unmap | awk '{ { print $3 } }'  >> stat
cat is_rep_mult_unmap | awk '{ { print $4 } }'  >> stat
rm is_rep_mult_unmap

############################################
## Extract reported reads only
############################################

cat s2 | awk 'BEGIN{ k=0;pos_in_read=0 }{ 
if($5>0 && length($3)<=5 && $3!="*")
print $1 " " $2 " " $3 " " $4 " " $6
}' > s3

rm s2
############################################
## Pairing
############################################

sed -i -r 's/chr//g' s3


#cat s3 > s4
sort -k1,1 -k4,4n s3 > s4

rm s3
# different chromosomes
cat s4 | awk 'BEGIN{n=0;k=0;is_dif_chr=0;}{
if (n!=$1) {
  if(k>1 && is_dif_chr==1){
    if(chr1<chr2){
      print chr1 " " pos1 " " chr2 " " pos2 " " str1 " " str2 " " r1 " " r2 " " n
    } else{
      print chr2 " " pos2 " " chr1 " " pos1 " " str2 " " str1 " " r2 " " r1 " " n
    }
  }
  n=$1
  k=0
  is_dif_chr=0
}

if(prev_chr!=$3 && k!=0){
  is_dif_chr=1
  chr1=prev_chr
  pos1=prev_pos
  str1=prev_str
  r1=prev_r
  chr2=$3
  pos2=$4
  str2=$2
  r2=$5
}

k=k+1
prev_chr=$3
prev_pos=$4
prev_str=$2
prev_r=$5
}' > allValidPairs_pre

# same chromosomes
cat s4 | awk 'BEGIN{n=0;k=0;is_dif_chr=0;prev_abs_val=0;}{
if (n!=$1) {
  if(k>1 && is_dif_chr==0){
    print chr1 " " pos1 " " chr2 " " pos2 " " str1 " " str2 " " r1 " " r2 " " n
    
  }
  n=$1
  k=0
  is_dif_chr=0
  chr1=$3
  pos1=$4
  str1=$2
  r1=$5
}


if(k!=0){
  chr2=$3
  pos2=$4
  str2=$2
  r2=$5
}

if(prev_chr!=$3 && k!=0){
  is_dif_chr=1
}

k=k+1
prev_chr=$3
}' >> allValidPairs_pre

rm s4

sed -i -r 's/X/2300000000/g' allValidPairs_pre
sed -i -r 's/Y/2400000000/g' allValidPairs_pre
cat allValidPairs_pre | awk '{
if($2!="M" && $3!="M"){
  print $0
}
}' > allValidPairs_pre1
rm allValidPairs_pre
sort -k1,1n -k3,3n -k2,2n -k4,4n allValidPairs_pre1 > allValidPairs_pre2
rm allValidPairs_pre1
sed -i -r 's/2300000000/X/g' allValidPairs_pre2
sed -i -r 's/2400000000/Y/g' allValidPairs_pre2



cat allValidPairs_pre2 | awk 'BEGIN{c1=0;c2=0;s1=0;s2=0}(c1!=$1 || c2!=$2 || s1!=$3 || s2!=$4){print;c1=$1;c2=$2;s1=$3;s2=$4}' > allValidPairs_pre3


###@@@@@@@@@@@@@@@@@@@@@@@@@@@@STAT
# report pairs before rm duplicates
cat allValidPairs_pre2 | wc -l >> stat
rm allValidPairs_pre2

###@@@@@@@@@@@@@@@@@@@@@@@@@@@@STAT
# report pairs after rm duplicates
cat allValidPairs_pre3 | wc -l >> stat




#cat allValidPairs_pre3 | awk '{
cat allValidPairs_pre3 | awk '{
if($7==$8){
  if($6==1){
    str2=0
  }
  if($6==0){
    str2=1
  }
  print $1 " " $2 " "$3 " " $4 " " $5 " " str2 " " $9
} else 
{

  print $1 " " $2 " "$3 " " $4 " " $5 " " $6 " " $9
}

}' > allValidPairs

rm allValidPairs_pre3
###@@@@@@@@@@@@@@@@@@@@@@@@@@@@STAT
# dangling ends
cat allValidPairs | awk 'BEGIN{same=0;dif=0;}{
if ($5==$6) {
  same=same+1
}
if ($5!=$6) {
  dif=dif+1
}
}END{ { print same " " dif } }' >> stat

###@@@@@@@@@@@@@@@@@@@@@@@@@@@@STAT
# cis trans
cat allValidPairs | awk 'BEGIN{cis=0;trans=0;cis_short=0;trans2=0;cis_long=0;}{
if ( $1==$3 && $5==$6) {
  cis=cis+1
}
if ( $1!=$3 && $5==$6) {
  trans=trans+1
}
if ( $1!=$3 && $5!=$6) {
  trans2=trans2+1
}
if ( $1==$3 &&  ($4-$2)<=20000 && $5==$6) {
  cis_short=cis_short+1
}
if ( $1==$3 &&  ($4-$2)>20000 && $5==$6) {
  cis_long=cis_long+1
}
}END{ { print cis " " trans " " cis_short " " cis_long " " trans2} }' >> stat

#############stat
cat stat | awk 'BEGIN{getline sample < "sample"}{

if (NR==1) {
  all_count=$1
}
if (NR==2) {
  rep=$1
}
if (NR==3) {
  mult=$1
}
if (NR==4) {
  unmap_one=$1
}
if (NR==5) {
  unmap_both=$1
}
if (NR==6) {
  rep_before_dup=$1
}
if (NR==7) {
  rep_after_dup=$1
}
if (NR==8) {
  same=$1
  dif=$2
}
if (NR==9) {
  cis=$1
  trans=$2
  cis_short=$3
  cis_long=$4
}
}END{
rep=100*rep/all_count
rep=rep - rep%1
mult=100*mult/all_count
mult=mult - mult%1
unmap_one=100*unmap_one/all_count
unmap_one=unmap_one - unmap_one%1
unmap_both=100*unmap_both/all_count
unmap_both=unmap_both - unmap_both%1


dups = 100*(1 - rep_after_dup/rep_before_dup)
dups = dups - dups%1

percent_count_of_valid_hic_interactions=percent_count_of_valid_hic_interactions - percent_count_of_valid_hic_interactions%1
de = 100*(dif - same)/rep_after_dup
de = de - de%1


count_of_valid_hic_interactions = (same)*2
percent_count_of_valid_hic_interactions = 100*count_of_valid_hic_interactions/all_count
percent_count_of_valid_hic_interactions=percent_count_of_valid_hic_interactions - percent_count_of_valid_hic_interactions%1

cis1 = 100*cis/(cis+trans)
cis1 = cis1 - cis1%1

trans1 = 100*trans/(cis+trans)
trans1 = trans1 - trans1%1

cis_short = 100*cis_short/(cis+trans)
cis_short = cis_short/1 - cis_short%1

cis_long = 100*cis_long/(cis+trans)
cis_long = cis_long/1 - cis_long%1

print "sample" "\t" "all_count_of_reads" "\t" "reported_pairs_fraction" "\t" "multiple_pairs_fraction" "\t" "unmap_one_pairs_fraction" "\t" "unmap_both_pairs_fraction" "\t" "duplicates_fraction" "\t" "valid_hic_interactions_count" "\t" "valid_hic_interactions_fraction" "\t" "dangling_ends" "\t" "cis_fraction" "\t" "trans_fraction" "\t" "cis_short_fraction" "\t" "cis_long_fraction"

print sample "\t" all_count "\t" rep "\t" mult "\t" unmap_one "\t" unmap_both "\t" dups "\t" count_of_valid_hic_interactions "\t" percent_count_of_valid_hic_interactions "\t" de "\t" cis1 "\t" trans1 "\t" cis_short "\t" cis_long
}'  > stat_res
rm stat

exit














