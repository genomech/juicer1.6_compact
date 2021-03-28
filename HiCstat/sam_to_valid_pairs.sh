#!/bin/bash


## The script will result in the file : allValidPairs
# allValidPairs fo#rmat: read_name chromosome1 position1 length1 strand1 chromosome2 position2 length2 strand2



SAM_FILE=$1
OUT_DIR=$2

cd ${OUT_DIR}



############################################
## Sam parsing
############################################
cat ${SAM_FILE} | awk '{
if (NR>94) {
  print $0
}
}' > s1


###@@@@@@@@@@@@@@@@@@@@@@@@@@@@STAT
cat s1 | awk 'BEGIN{k=0;n=0}{
  if (n!=$1) {
    n=$1
    k=k+1
}
}END{ { print k } }' >> stat


cat s1 | awk '{ print $1 " " $2 " " $3 " " $4 " " $5 " " $6 }' > s2
#rm s1

############################################
## Len of alingment
############################################

cat s2 | awk '{ print $6 }' > s2_cigar


sed -i -r 's/I/M/g' s2_cigar
sed -i -r 's/(D|S|H)/ /g' s2_cigar
sed -i -r 's/M/M /g' s2_cigar

cat s2_cigar | awk '{ m = ""; for (i = 1; i <= NF; ++i) { if ($i ~ "M") { m = m$i}};print m  }' | awk -F'M' '{ m = 0;for (i = 1; i <= NF; ++i) {  m = m+$i}; print m }' > len

#rm s2_cigar


############################################
## Orientation
############################################

## awk get 5 bit
cat s2 | awk '
function d2b(d,  b) {
      k = 0
      while(k != 5) {
          k = k + 1
          b=d%2
          d=int(d/2)
      }
      return(b)
}
{
    print d2b($2)
}' > is_reverse

## awk get 7 bit
cat s2 | awk '
function d2b(d,  b) {
      k = 0
      while(k != 7) {
          k = k + 1
          b=d%2
          d=int(d/2)
      }
      return(b)
}
{
    print 1 - d2b($2)
}' > is_r1

## awk get 12 bit
cat s2 | awk '
function d2b(d,  b) {
      k = 0
      while(k != 12) {
          k = k + 1
          b=d%2
          d=int(d/2)
      }
      return(b)
}
{
    print 1 - d2b($2)
}' > is_sapl

############################################
## Cbind all extracted data
############################################

cat s2 | awk 'BEGIN{ k=0;pos_in_read=0 }{ 
k=k+1
getline len < "len"
getline pos_in_read_left < "pos_in_read_left"
getline pos_in_read_right < "pos_in_read_right"
getline is_reverse < "is_reverse"
getline is_r1 < "is_r1"
getline is_sapl < "is_sapl"
print $1 " " $2 " " $3 " " $4 " " $5 " " $6 " " len " " k " " 0 " " is_reverse " " is_sapl " " is_r1

}' > s3
#rm s2



############################################
## Identificate reported pairs, multiple alignment, unmapped
############################################

cat s3 | awk 'BEGIN{k=0;prev_name=1;unmap_r1=0;unmap_r2=0;mult=0;mult_r1=0;mult_r2=0;unmap_both=0;unmap_one=0;is_rep=0;is_rep_r1=0;is_rep_r2=0;st=0}{
if (prev_name!=$1 && st!=0) {
  if ((unmap_r1*(1 - is_rep_r1) + unmap_r2*(1 - is_rep_r2))==1){
    unmap_one=1
  }
  if ((unmap_r1*(1 - is_rep_r1) + unmap_r2*(1 - is_rep_r2))==2){
    unmap_both=1
  }
  if ((is_rep_r1 + is_rep_r2)==2){
    is_rep=1
  }
  if (mult_r1*(1 - unmap_r1)*(1 - is_rep_r1)==1 || mult_r2*(1 - unmap_r2)*(1 - is_rep_r2)==1){
    mult=1
  }
  print is_rep " " mult " " unmap_both " " unmap_one
  unmap_r1=0
  unmap_r2=0
  unmap_both=0
  unmap_one=0
  mult=0
  mult_r1=0
  mult_r2=0
  is_rep=0
  is_rep_r1=0
  is_rep_r2=0
  k=0
}
st=1;
prev_name=$1
k=k+1
if ($12==0) {
  if ($5==0 || length($3)>5) {
    if ($3=="*" || length($3)>5) {
      unmap_r1=1
    } else {
      mult_r1=1
    }
  } else {
    is_rep_r1=1
  }
} else {
  if ($5==0 || length($3)>5) {
    if ($3=="*" || length($3)>5) {
      unmap_r2=1
    } else {
      mult_r2=1
    }
  } else {
    is_rep_r2=1
  }
}
}' > is_rep_mult_unmap


###@@@@@@@@@@@@@@@@@@@@@@@@@@@@STAT
#is_rep, mult, unmap_both, unmap_one 
cat is_rep_mult_unmap | awk '{ { print $1 } }' | grep "1" |  wc -l >> stat
cat is_rep_mult_unmap | awk '{ { print $2 } }' | grep "1" |  wc -l >> stat
cat is_rep_mult_unmap | awk '{ { print $3 } }' | grep "1" |  wc -l >> stat
cat is_rep_mult_unmap | awk '{ { print $4 } }' | grep "1" |  wc -l >> stat
#rm is_rep_mult_unmap

############################################
## Extract reported pairs only
############################################


cat s3 | awk 'BEGIN{k=0;prev_name=0;is_rep=0;is_rep_r1=0;is_rep_r2=0;st=0;}{
if (prev_name!=$1 && st!=0) {
  if ((is_rep_r1 + is_rep_r2)==2){
    is_rep=1
  }
  for (i = 1; i <= k; ++i) {
    print is_rep
  }
  is_rep=0
  is_rep_r1=0
  is_rep_r2=0
  k=0
}
st=1
getline is_r1 < "is_r1"
prev_name=$1
k=k+1
if (is_r1==1) {
  if ($5!=0 && length($3)<=5 && $3!="*") {
    is_rep_r1=1
  }
} else {
  if ($5!=0 && length($3)<=5 && $3!="*") {
    is_rep_r2=1
  }
}
}' > is_rep


cat s3 | awk '{
getline is_rep < "is_rep"
if (is_rep==1 && $5!=0 && length($3)<=5 && $3!="*") {
  print $0
}
}' > s4
#rm s3



############################################
## Pairing
############################################



cat s4 | awk 'BEGIN{k=1;prev_r=0;a=0;chr_r1="no_chr";st=0;chr_dif=0;}{
if (0==$12 && prev_r==1) {
  r1=$0
  chr_r1=$3
  a=1
}

chr_r2=$3
if (a==1 && chr_r1!=chr_r2) {
  print r1
  print $0
  a=0
}

prev_r=$12
}' > s5_1

## same_chr_all_name
cat s4 | awk 'BEGIN{k=0;prev_r=0;a=0;chr_r1="no_chr";st=0;chr_diff=0;}{
if (0==$12 && prev_r==1) {
  r1=$0
  chr_r1=$3
  if ( chr_diff==1 ){
    for (i = 1; i <= k; ++i) {
      print 0
    }
  } else {
    for (i = 1; i <= k; ++i) {
      print 1
    }
  }
  chr_diff=0
  a=1
  k=1
} else {
  k=k+1
}

chr_r2=$3
if (a==1 && chr_r1!=chr_r2) {
  chr_diff=1
  a=0
}

prev_r=$12
}' > same_chr_all_name

# get same_chr_all_name only
cat s4 | awk '{
getline same_chr_all_name < "same_chr_all_name"
if (same_chr_all_name==1) {
  print $0
}
}' > s5_2
#rm s4
#rm same_chr_all_name

# remaining
cat s5_2 | awk 'BEGIN{k=0;prev_r=0;st=0;a=0;d=-160;r1="";}{
st_r2=$4
en_r2=$4+$7

if (0==$12 && prev_r==1) {
  if ( st>1 ) {
    print r1
    print r2
  }
  st=st+1
  d=-160
  r1=$0
  st_r1=$4
  en_r1=$4+$7
  k=0
}
k=k+1
if ( r1!="" && k!=1) {
  if ( st_r1<st_r2 ){
    if ( d<(st_r2-en_r1) ) {
      r2=$0
      d=st_r2-en_r1
    }
  } else {
    if ( d<(st_r1-en_r2) ) {
      r2=$0
      d=st_r1-en_r2
    }
  }
}

prev_r=$12
}' >> s5_1
#rm s5_2



### #rm "chr"
sed -i -r 's/chr//g' s5_1



############################################
## allValidInteractions
############################################

# only neighbour in chimeras are pair

cat s5_1 | awk 'BEGIN{k=0;prev_name=0;is_rep=0;is_rep_r1=0;is_rep_r2=0;st=0;pos2=0;pos22=0;}{
if (prev_name==$1) {
  chr2=$3
  pos2=$4
  len2=$7
  str2=$10
  if ( 1==1 ) {
    if ( chr1==chr2 && pos1>pos2) {
      print $1 " " chr2 " " pos2 " " len2 " " str2 " " chr1 " " pos1 " " len1 " " str1
    }
    if ( chr1==chr2 && pos1<=pos2 ) {
      print $1 " " chr1 " " pos1 " " len1 " " str1 " " chr2 " " pos2 " " len2 " " str2
    }
    if ( chr1>chr2) {
      print $1 " " chr2 " " pos2 " " len2 " " str2 " " chr1 " " pos1 " " len1 " " str1
    }
    if ( chr1<chr2) {
      print $1 " " chr1 " " pos1 " " len1 " " str1 " " chr2 " " pos2 " " len2 " " str2
    }
  }
}
prev_name=$1
chr1=$3
pos1=$4
len1=$7
str1=$10
}' > allValidPairs_pre
#rm s5_1


############################################
## #rm duplicates
############################################


sed -i -r 's/X/2300000000/g' allValidPairs_pre
sed -i -r 's/Y/2400000000/g' allValidPairs_pre
grep -v "M" allValidPairs_pre > allValidPairs_pre1
#rm allValidPairs_pre
sort -k2,2n -k6,6n -k3,3n -k7,7n allValidPairs_pre1 > allValidPairs_pre2
#rm allValidPairs_pre1
sed -i -r 's/2300000000/X/g' allValidPairs_pre2
sed -i -r 's/2400000000/Y/g' allValidPairs_pre2



cat allValidPairs_pre2 | awk 'BEGIN{c1=0;c2=0;s1=0;s2=0}(c1!=$2 || c2!=$3 || s1!=$6 || s2!=$7){print;c1=$2;c2=$3;s1=$6;s2=$7}' > allValidPairs


###@@@@@@@@@@@@@@@@@@@@@@@@@@@@STAT
# report pairs before #rm duplicates
cat allValidPairs_pre2 | wc -l >> stat
#rm allValidPairs_pre2

###@@@@@@@@@@@@@@@@@@@@@@@@@@@@STAT
# report pairs after #rm duplicates
cat allValidPairs | wc -l >> stat

###@@@@@@@@@@@@@@@@@@@@@@@@@@@@STAT
# dangling ends
cat allValidPairs | awk 'BEGIN{FF=0;RR=0;RF=0;FR=0;}{
if ( $2==$6 &&  $5==1 && $9==1) {
  RR=RR+1
}
if ( $2==$6 &&  $5==0 && $9==0) {
  FF=FF+1
}
if ( $2==$6 &&  $5==1 && $9==0) {
  RF=RF+1
}
if ( $2==$6 &&  $5==0 && $9==1) {
  FR=FR+1
}
}END{ { print RR " " FF " " FR " " RF } }' >> stat

###@@@@@@@@@@@@@@@@@@@@@@@@@@@@STAT
# cis trans
cat allValidPairs | awk 'BEGIN{cis=0;trans=0;cis_short=0;trans2=0;cis_long=0;}{
if ( $2==$6 && $5==$9) {
  cis=cis+1
}
if ( $2!=$6 && $5==$9) {
  trans=trans+1
}
if ( $2!=$6 && $5!=$9) {
  trans2=trans2+1
}
if ( $2==$6 &&  ($7-$3)<=20000 && $5==$9) {
  cis_short=cis_short+1
}
if ( $2==$6 &&  ($7-$3)>20000 && $5==$9) {
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
if (NR==4) {
  unmap_both=$1
}
if (NR==5) {
  rep_before_dup=$1
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














