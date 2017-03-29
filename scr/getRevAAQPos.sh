#To get the position of repeat on the query
#$1 is the blast result of blBtBlkSq.sh
#$2 is the id of the query position file

grep -A5 'Strand' $1|grep 'Query'|cut -d ' ' -f2 >temp.qstart #A's start position in query
grep -B6 'Score'  $1|grep 'Query'|awk '{print $4}' >temp.qend    #A's end position in query
grep -B8 'Gapped' $1|grep 'Query'|awk '{print $4}' >>temp.qend
grep "Strand" $1|cut -d '=' -f2 >temp.strand
paste -d" " temp.qstart temp.qend temp.strand >temp"$2".poss #all the alignment start, end, strand
grep "Plus / Minus" temp"$2".poss|awk '{print $1,$2}' >temp"$2".pos #get the A/-A alignment start, end