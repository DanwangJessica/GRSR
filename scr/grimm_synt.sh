#$1 is the minimum block size
#$2 is the maximum gap threshold
export PATH="$PATH:$HOME/GRIMM-synteny/GRIMM_SYNTENY-2.02"

#To filter out conflicting coordinates
#mkdir ./Example/2.Scaffolds/anchors
grimm_synt -A -f ./Example/2.Scaffolds/core_coords.txt -d ./Example/2.Scaffolds/anchors

#To get the scaffold map and k-way format blocks
#mkdir ./Example/2.Scaffolds/anchors_$1_$2_c
grimm_synt -f ./Example/2.Scaffolds/anchors/unique_coords.txt -d ./Example/2.Scaffolds/anchors_$1_$2_c -c -Q -m $1 -g $2


