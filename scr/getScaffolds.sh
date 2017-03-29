#$1 is the minimum block size
#$2 is the maximum gap threshold

#Delete strain-specific blocks
declare -i n=$(grep -n 'label=u0' ./Example/1.MSA_results/MSA_Result.maf|cut -d ':' -f1)
n=$n-1
sed ''$n',$d' ./Example/1.MSA_results/MSA_Result.maf >./Example/2.Scaffolds/cordis.maf

javac getCoreKwayf.java
java getCoreKwayf ./Example/2.Scaffolds/cordis.maf >./Example/2.Scaffolds/core_coords.txt

sh grimm_synt.sh $1 $2
