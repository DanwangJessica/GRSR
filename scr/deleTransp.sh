#$1 is the mgr_macro.txt file
#$2 is id for genome 0
#$3 is id for genome 1

g0=$(grep -A2 ">genome$2$" $1|sed '1d'|sed '1d')

g1=$(grep -A2 ">genome$3$" $1|sed '1d'|sed '1d')

blks=( $(grep -A2 ">genome$2$" $1|sed '1d'|sed '1d') )

declare -i blksNo=${#blks[@]}-1

java deleTransp $blksNo $g0 $g1
