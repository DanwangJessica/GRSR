#Usage: sh getApos.sh D_Sq.fa D
#To get the position of IR in each genome, usually use this after run the mapA.sh
#$1 is the name of the IR
#$2 is the letter rep of the IR
#$3 is the directory to store the output

declare -i totGenome="$(grep -n '^' faNamesPre.txt|tail -n 1|cut -d ':' -f1)"

for (( i=1; i<$totGenome+1; i=i+1 ))
do  
    declare -i result=$(grep -n 'No hits found' $3/"$i"_$1map.txt|cut -d ':' -f1)
    if [ $result -ne 0 ]; then
        echo 'No hits found' >$3/"$i"$2.pos
        continue
    fi

    grep -A5 'Strand' $3/"$i"_$1map.txt|grep 'Sbjct'|cut -d ' ' -f2 >$3/$2.start #A's start position in strain i 
    grep -B6 'Score' $3/"$i"_$1map.txt|grep 'Sbjct'|awk '{print $4}' >$3/$2.end #A's end position in strain i
    grep -B8 'Gapped' $3/"$i"_$1map.txt|grep 'Sbjct'|awk '{print $4}' >>$3/$2.end  
    
    grep -A5 'Strand' $3/"$i"_$1map.txt|grep 'Query'|cut -d ' ' -f2 >$3/$2.qstart #A's start position in itself seq 
    grep -B6 'Score' $3/"$i"_$1map.txt|grep 'Query'|awk '{print $4}' >$3/$2.qend    #A's end position in itself seq
    grep -B8 'Gapped' $3/"$i"_$1map.txt|grep 'Query'|awk '{print $4}' >>$3/$2.qend

    paste -d" " $3/$2.start $3/$2.end $3/$2.qstart $3/$2.qend >$3/"$i"$2.pos
done
