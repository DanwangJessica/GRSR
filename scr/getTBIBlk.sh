#$1 is current LBlk
#calculate the real left and right block involved in Transposition or Block Interchange
declare -i result=$(grep "=$1$" BefFilter.txt|cut -d '=' -f1)
declare -i s=$(grep "=$result$" AftMerge.txt|cut -d ':' -f1)
declare -i e=$(grep "=$result$" AftMerge.txt|cut -d ':' -f2|cut -d '=' -f1)
declare -i sB=$(grep "=$s$" BefMerge.txt|cut -d '=' -f1)
declare -i eB=$(grep "=$e$" BefMerge.txt|cut -d '=' -f1)
echo $sB $eB
