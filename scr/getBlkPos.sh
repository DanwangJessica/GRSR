declare -i strainID=$2
declare -i blkID=$3
blockfile=$1
declare -i blkStartF=$strainID*4-1

declare -i blkSizeF=$blkStartF+1

declare -i blkStrandF=$blkSizeF+1

declare -i blkStart="$(grep '^'$blkID' ' $blockfile|cut -d ' ' -f$blkStartF)"

declare -i blkSize="$(grep '^'$blkID' ' $blockfile|cut -d ' ' -f$blkSizeF)"

declare -i blkStrand="$(grep '^'$blkID' ' $blockfile|cut -d ' ' -f$blkStrandF)"


declare -i blkEnd=$blkStart+$blkSize-1
echo $blkStart $blkEnd $blkStrand
