#get the sequence of a block
#$1 is the blocks.txt file
#$2 is the strain ID
#$3 is the start block ID for the left region
#$4 is the end block ID for the left region
#$5 is the start block ID for the right region
#$6 is the end block ID for the right region
#$7 is the faNamesPre.txt file
#$8 is the path for genomes
#$9 is the path for storing the reportA.txt
#$10 is the rearrangment type,eg. tp,itp,bi,ibi,hibi
export PATH="$PATH:$HOME/blast-2.2.26/bin"
#$3,4,5,6 is the number in the permutation which may be negative

#sequence 1
lpos1=( $(sh getBlkPos.sh $1 $2 ${3#-}) )
declare -i start1=${lpos1[1]}+1
rpos1=( $(sh getBlkPos.sh $1 $2 ${4#-}) )
declare -i end1=${rpos1[0]}-1
declare -i size1=$end1-$start1+1


#sequence 2
lpos2=( $(sh getBlkPos.sh $1 $2 ${5#-}) )
declare -i start2=${lpos2[1]}+1
rpos2=( $(sh getBlkPos.sh $1 $2 ${6#-}) )
declare -i end2=${rpos2[0]}-1
declare -i size2=$end2-$start2+1


if [ $size1 -gt 20 ] && [ $size2 -gt 20 ]; then
	sh extractSq.sh $2 $start1 $end1 $7 $8 >$9/$2_bt_$3_$4.fa
	sh extractSq.sh $2 $start2 $end2 $7 $8 >$9/$2_bt_$5_$6.fa
	bl2seq -i $9/$2_bt_$3_$4.fa -j $9/$2_bt_$5_$6.fa -p blastn -e 1e-50 -o $9/$2_${10}_$4_$5.txt
else
	echo "Too short sequences!" >$9/$2_${10}_$4_$5.txt					
fi
