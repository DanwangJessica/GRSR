#$1 is leftblock ID
#$2 is rightblock ID
#$3 is the block.txt file
#$4 is the strain ID

#$5 is the blast output file
#$6 is the blast length cutoff
#$7 is the blast similarity cutoff

declare -i LBlk=$1
declare -i RBlk=$2
declare -i lSign=0;
declare -i rSign=0;

declare -i absLBlk=$LBlk;
declare -i absRBlk=$RBlk;
if [ $LBlk -lt 0 ]; then
    absLBlk=-1*$LBlk
fi

if [ $RBlk -lt 0 ]; then
    absRBlk=-1*$RBlk
fi

lspos=( $(sh getBlkPos.sh $3 $4 $absLBlk) )
lSign=${lspos[2]}
rspos=( $(sh getBlkPos.sh $3 $4 $absRBlk) )
rSign=${rspos[2]}

declare -i time=$LBlk*$RBlk*$lSign*$rSign
if [ $time -gt 0 ]; then
    sh searchRevAA.sh $5 $6 $7
else
    sh searchAA.sh $5 $6 $7
fi                                                                       
