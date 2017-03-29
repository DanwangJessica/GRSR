#$1 is the blocks.txt file
#$2 is the genome id
#$3 is the reversed left block id,can be positive or negative
#$4 is the search range
#$5 is the faNamesPre.txt file
#$6 is the genome absolute path

declare -i sqEnd=0
declare -i sqStart=0
declare -i RANGE=$4

if [ $3 -gt 0 ]; then
    lspos=( $(sh getBlkPos.sh $1 $2 $3) ) #lspos is an array, which stores the Left blk's start,end,strand on the certain genome
    if [ ${lspos[2]} -eq 1 ]; then
        sqEnd=${lspos[0]}+$RANGE
        sqStart=${lspos[0]}-$RANGE
        sh extractSq.sh $2 $sqStart $sqEnd $5 $6
    else
        sqStart=${lspos[1]}-$RANGE
        sqEnd=${lspos[1]}+$RANGE
        sh extractSq.sh $2 $sqStart $sqEnd $5 $6
    fi
else
    declare -i absLBlk="$((-1*$3))"
    lspos=( $(sh getBlkPos.sh $1 $2 $absLBlk) ) #lspos is an array, which stores the Left blk's start,end,strand on the certain genome
    if [ ${lspos[2]} -eq -1 ]; then
        sqEnd=${lspos[0]}+$RANGE
        sqStart=${lspos[0]}-$RANGE
        sh extractSq.sh $2 $sqStart $sqEnd $5 $6
    else
        sqStart=${lspos[1]}-$RANGE
        sqEnd=${lspos[1]}+$RANGE
        sh extractSq.sh $2 $sqStart $sqEnd $5 $6
    fi

fi
