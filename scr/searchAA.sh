#$1 is the blast result
#$2 is the A/-A length cutoff
#$3 is the A/-A similarity cutoff
grep -B1 'Strand = Plus / Plus' $1 >tempA.txt
if [ -s tempA.txt ]; then
    ALength=( $(grep 'Identities' tempA.txt|cut -d '/' -f2|cut -d ' ' -f1) )
    ASim=( $(grep 'Identities' tempA.txt|cut -d '(' -f2|cut -d '%' -f1) )



    declare -i lengthCut=$2
    declare -i simCut=$3
    declare -i num=${#ALength[@]}
    declare -i maxLen=0
    declare -i maxSim=0

    for (( i=0; i<$num; i=i+1 ))
    do
        if [ ${ALength[$i]} -ge $lengthCut ] && [ ${ASim[$i]} -ge $simCut ]; then
            if [ ${ALength[$i]} -gt $maxLen ]; then
                maxLen=${ALength[$i]}
                maxSim=${ASim[$i]}
            fi
        fi
    done
    
    if [ $maxLen -ne 0 ]; then
        echo "Found,Length = $maxLen,Simlarity = $maxSim%"
    else
        echo "N"
    fi

else
    echo "N"
fi
