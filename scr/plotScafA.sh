#$1 is the sequence for IR
#$2 is the faNamesPre.txt file
#$3 is the genome path
#$4 is the IR's letter representative
#$5 is the path .pos file.
#$6 is the blocks.txt path and file
#$7 is the length of IR
#$8 is the threshold for short IR

declare -i totGenome="$(grep -n '^' $2|tail -n 1|cut -d ':' -f1)"

sh mapA.sh $1 $2 $3 $5
sh getApos.sh $1 $4 $5
java getScafA $5 $6 $4 $totGenome $7 $8

rm $5/*.pos
rm $5/*.start
rm $5/*.end
rm $5/*.qstart
rm $5/*.qend


