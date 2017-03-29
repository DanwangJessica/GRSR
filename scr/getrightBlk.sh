#$1 is current LBlk
#calculate the real left block being reversed
declare -i result=$(grep "=$1$" AftFilter.txt|cut -d '=' -f1)
result=$(grep "=$result$" BefFilter.txt|cut -d '=' -f1)
result=$(grep "=$result$" AftMerge.txt|cut -d ':' -f2|cut -d '=' -f1) #because it is get right blk id
result=$(grep "=$result$" BefMerge.txt|cut -d '=' -f1)
echo $result
