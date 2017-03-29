#$1 is the strain ID, start from 1
#$2 is the start
#$3 is the end
#$4 is the faNamesPre.txt file which store the genome name in correct order
#$5 is the path to store the genomes

declare -i start=$2
declare -i end=$3

namePre=$(sed -n ''$1'p' $4)
faName=$(find $5 -name "*$namePre*"|sed -n '1p')
echo ">$1 $faName $start $end"


#get the nucleotide number per Line in the genome file
nucLine=$(sed -n '2p' $faName)
declare -i nucPerLine=${#nucLine} 

#Extract Sequence Part
declare -i start_line=$start/$nucPerLine+1
declare -i end_line=$end/$nucPerLine+1
declare -i start_col=($start)%$nucPerLine
declare -i end_col=($end)%$nucPerLine
declare -i tot_line=$end_line-$start_line+1
declare -i end_pos=$end-$start+$start_col

if [ $start_col -ne 0 ] && [ $end_col -ne 0 ]; then
  grep -A$end_line '>' $faName|sed '1,'$start_line'd'|awk '{if(NR%'$tot_line'==0) ORS="\n";else ORS="";print}'|cut -c $start_col-$end_pos
elif [ $start_col -eq 0 ] && [ $end_col -ne 0 ]; then
  declare -i start_line1=$start_line-1
  declare -i tot_line1=$tot_line+1
  declare -i end_pos1=$end-$start+$nucPerLine
  grep -A$end_line '>' $faName|sed '1,'$start_line1'd'|awk '{if(NR%'$tot_line1'==0) ORS="\n";else ORS="";print}'|cut -c $nucPerLine-$end_pos1
elif [ $start_col -ne 0 ] && [ $end_col -eq 0 ]; then
  declare -i end_line2=$end_line-1
  declare -i tot_line2=$tot_line-1
  grep -A$end_line2 '>' $faName|sed '1,'$start_line'd'|awk '{if(NR%'$tot_line2'==0) ORS="\n";else ORS="";print}'|cut -c $start_col-$end_pos
else
  declare -i end_line3=$end_line-1
  declare -i start_line3=$start_line-1
  declare -i end_pos3=$end-$start+$nucPerLine
  grep -A$end_line3 '>' $faName|sed '1,'$start_line3'd'|awk '{if(NR%'$tot_line'==0) ORS="\n";else ORS="";print}'|cut -c $nucPerLine-$end_pos3
fi


