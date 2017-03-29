#!/bin/bash
#$1 is the scaffold_map file,like mgr_macro.txt
#$2 is the block.txt file
#$3 is the maf file
#$4 is the absolute genome path
#$5 is the absolute path for reportA
#$6 is the length of A threshold in bp
#$7 is the similarity threshold of A and -A %
#$8 is the search range of repeats in reversals
export PATH="$PATH:$HOME/GRIMM/GRIMM-2.01"
export PATH="$PATH:$HOME/blast-2.2.26/bin"

declare -i totGenomes="$(grep '^a' $3|sed '2,$d'|cut -d '=' -f4)"
declare -i RANGE=$8 #the search range of A, =5000 for all the result
declare -i revCut=5 #the threshold for reversd no. of blocks
declare -i ALCut=$6 #the length of A threshold
declare -i ASimCut=$7 #the similarity threshold of A and -A
#declare -i curDis=$6

declare -i distance=0
declare -i LBlk=0
declare -i RBlk=0

declare -i count=0
declare -i totRev=0
declare -i countBoth=0

declare -i totTp=0
declare -i cTp=0 #number of transposition found repeats
declare -i cTpBoth=0
declare -i totItp=0
declare -i cItp=0
declare -i cItpBoth=0
declare -i totBi=0
declare -i cBi=0
declare -i cBiBoth=0
declare -i totIbi=0
declare -i cIbi=0
declare -i cIbiBoth=0
declare -i totHib=0
declare -i cHib=0
declare -i cHibBoth=0

declare -i n=0
#grimm -f $1 -C >$5/reportA.txt
echo "Pairwise Genome Rearrangment Scenarios with repeats information:" >$5/reportA.txt


declare -i endLine=2+$totGenomes
grep -A$totGenomes 'mult='$totGenomes'' $3|sed ''$endLine',$d'|sed '1d'|cut -d ' ' -f2|cut -d '.' -f1 >faNamesPre.txt


for (( i=1; i<$totGenomes+1; i=i+1 ))
do
    for (( j=$i+1; j<$totGenomes+1; j=j+1 ))
    do
        #delete the transposition and block interchange cases in genome i and genome j
        sh deleTransp.sh $1 $i $j >tbi.txt
		declare -i tbiS="$(head -n 1 tbi.txt)"  #steps except inversion
		
        grimm -f mgr_macro2.txt -C -g 1,2 >$5/sn"$i"To"$j".sn #reversal senario from genome i to j
        distance="$(grep 'Reversal Distance' $5/sn"$i"To"$j".sn|cut -d$'\t' -f2)"  #reversal distance calculated by GRIMM
        
        declare -i totS=$tbiS+$distance
        echo ">From Genome $i to $j, total rearrangment step(s): $totS" >>$5/reportA.txt
		if [ $totS -eq 0 ]; then
			continue
		fi
		
		#The following deals with tranposition and block interchange
		if [ $tbiS -ne 0 ]; then
			gs=$(grep -A2 ">genome$i$" $1|sed '1d'|sed '1d') #genome for source permutation
			gd=$(grep -A2 ">genome$j$" $1|sed '1d'|sed '1d') #genome for destination permutation
			steps=( $(sed -n '2p' tbi.txt) )
			declare -i tp=${steps[0]}
			declare -i itp=${steps[1]}
			declare -i bi=${steps[2]}
			declare -i ibi=${steps[3]}
			declare -i hbi=${steps[4]}
			totTp=$totTp+$tp
			totItp=$totItp+$itp
			totBi=$totBi+$bi
			totIbi=$totIbi+$ibi
			totHib=$totHib+$hbi
			#Transpostion  steps
			if [ $tp -ne 0 ]; then
				echo "*Transposition: $tp step(s):" >>$5/reportA.txt
				#show transposition steps
				tpBlks=( $(grep '^T' tbi.txt) )
				n=${#tpBlks[@]}
				for (( k=1; k<$n; k=k+1 )) #check each block involved transposition
				do
					declare -i bID=${tpBlks[$k]}
					tpSE=( $(sh getTBIBlk.sh $bID) ) #a transposition region start and end block on s

					#search repeats on source genome
					declare -i nbLs=$(java getPrevBlk ${tpSE[0]} $gs) #left neighbor block on s
					declare -i nbRs=$(java getPostBlk ${tpSE[1]} $gs) #right neighbor block on s
					if [ $nbLs -ne 0 ] && [ $nbRs -ne 0 ]; then
						sh blBtBlkSq.sh $2 $i $nbLs ${tpSE[0]} ${tpSE[1]} $nbRs faNamesPre.txt $4 $5 tp #blast result now in $5/$i_tp_${tpSE[0]}_${tpSE[0]}.txt
						searchSA=$(sh searchAA.sh $5/"$i"_tp_${tpSE[0]}_${tpSE[1]}.txt $ALCut $ASimCut)
					else
						searchSA="N"
					fi
					
					#search repeats on destination genome
					nbd=( $(java getNbD ${tpSE[0]} ${tpSE[1]} $gd) )  #get the left neighbour, start of the block, end of the block, right neighour on the orders of d
					declare -i nbLd=${nbd[0]} #left neighbour of the rearrangement region on d
					declare -i sd=${nbd[1]} #start block of the rearrangement region on d
					declare -i ed=${nbd[2]} #end block of the rearrangement region on d
					declare -i nbRd=${nbd[3]} # right neighbour of the rearrangement region on d
					if [ $nbLd -ne 0 ]; then
						sh blBtBlkSq.sh $2 $j $nbLd $sd $ed $nbRd faNamesPre.txt $4 $5 tp #blast result now in $5/$j_tbi${tpSE[0]}_${tpSE[0]}.txt
						searchDA=$(sh searchAA.sh $5/"$j"_tp_"$sd"_$ed.txt $ALCut $ASimCut)
					else
						searchDA="N"
					fi
					#update $searchSA, when there exist a pair of repeats on s
					if [ "$searchSA" != "N" ]; then
						brks=( $(java getNbD $nbLd $nbRd $gs) ) #from the neighbors on d to get the insertion point on s ${brks[1]},${brks[2]}
						sh blBtBlkSq.sh $2 $i $nbLs ${tpSE[0]} ${brks[1]} ${brks[2]} faNamesPre.txt $4 $5 tp #blast result now in $5/$i_tp_$4_$5.txt
						searchSA1=$(sh searchAA.sh $5/"$i"_tp_${tpSE[0]}_${brks[1]}.txt $ALCut $ASimCut)
						if [ "$searchSA1" != "N" ]; then
							sh getAAQPos.sh $5/"$i"_tp_${tpSE[0]}_${tpSE[1]}.txt 1 #temp1.pos
							sh getAAQPos.sh $5/"$i"_tp_${tpSE[0]}_${brks[1]}.txt 2
							java getAInts temp1.pos temp2.pos $ALCut >temp1a2.pos
							maxL=$(java getMaxL temp1a2.pos)
							if [ $maxL -ne 0 ]; then
								searchSA="Found,Length=$maxL"
							else
								searchSA="N"
							fi
						else
							searchSA="N"
						fi
					fi

					#update $searchDA, when there exist a pair of repeats on d
					if [ "$searchDA" != "N" ]; then
						brks=( $(java getNbD $nbLs $nbRs $gd) )
						sh blBtBlkSq.sh $2 $j $nbLd $sd ${brks[1]} ${brks[2]} faNamesPre.txt $4 $5 tp #blast result now in $5/$i_tp_$4_$5.txt
						searchDA1=$(sh searchAA.sh $5/"$j"_tp_"$sd"_${brks[1]}.txt $ALCut $ASimCut)
						if [ "$searchDA1" != "N" ]; then
							sh getAAQPos.sh $5/"$j"_tp_"$sd"_$ed.txt 1 #temp1.pos
							sh getAAQPos.sh $5/"$j"_tp_"$sd"_${brks[1]}.txt 2
							java getAInts temp1.pos temp2.pos $ALCut >temp1a2.pos
							maxL=$(java getMaxL temp1a2.pos)
							if [ $maxL -ne 0 ]; then
								searchDA="Found,Length=$maxL"
							else
								searchDA="N"
							fi
						else
							searchDA="N"
						fi
					fi
					echo "Step $k: ${tpSE[0]} through ${tpSE[1]} Transposition:Source:$searchSA Destination:$searchDA" >>$5/reportA.txt

					if [[ $searchSA == *"Found"* ]] || [[ $searchDA == *"Found"* ]]; then
						cTp=$cTp+1          
					fi
					if [[ $searchSA == *"Found"* ]] && [[ $searchDA == *"Found"* ]]; then
						cTpBoth=$cTpBoth+1 
					fi
				done
				
			fi
			
			if [ $itp -ne 0 ]; then
				echo "*Inverted Transposition: $itp step(s):" >>$5/reportA.txt
				tpBlks=( $(grep '^IT' tbi.txt) )
				n=${#tpBlks[@]}
				for (( k=1; k<$n; k=k+1 )) #check each block involved transposition
				do
					declare -i bID=${tpBlks[$k]}
					tpSE=( $(sh getTBIBlk.sh $bID) )

					declare -i nbLs=$(java getPrevBlk ${tpSE[0]} $gs) #left neighbor block on s
					declare -i nbRs=$(java getPostBlk ${tpSE[1]} $gs) #right neighbor block on s
					if [ $nbLs -ne 0 ] && [ $nbRs -ne 0 ]; then
						sh blBtBlkSq.sh $2 $i $nbLs ${tpSE[0]} ${tpSE[1]} $nbRs faNamesPre.txt $4 $5 itp #blast result now in $5/$i_tp_${tpSE[0]}_${tpSE[0]}.txt
						searchSA=$(sh searchAA.sh $5/"$i"_itp_${tpSE[0]}_${tpSE[1]}.txt $ALCut $ASimCut)
					else
						searchSA="N"
					fi
					
					#search repeats on destination genome
					nbd=( $(java getNbD ${tpSE[0]} ${tpSE[1]} $gd) )  #get the left neighbour, start of the block, end of the block, right neighour on the orders of d
					declare -i nbLd=${nbd[0]} #left neighbour of the rearrangement region on d
					declare -i sd=${nbd[1]} #start block of the rearrangement region on d
					declare -i ed=${nbd[2]} #end block of the rearrangement region on d
					declare -i nbRd=${nbd[3]} # right neighbour of the rearrangement region on d
					if [ $nbLd -ne 0 ]; then
						sh blBtBlkSq.sh $2 $j $nbLd $sd $ed $nbRd faNamesPre.txt $4 $5 itp #blast result now in $5/$j_tbi${tpSE[0]}_${tpSE[0]}.txt
						searchDA=$(sh searchAA.sh $5/"$j"_itp_"$sd"_$ed.txt $ALCut $ASimCut)
					else
						searchDA="N"
					fi
					#update $searchSA, when there exist a pair of repeats on s
					if [ "$searchSA" != "N" ]; then
						brks=( $(java getNbD $nbLd $nbRd $gs) ) #from the neighbors on d to get the insertion point on s ${brks[1]},${brks[2]}
						sh blBtBlkSq.sh $2 $i $nbLs ${tpSE[0]} ${brks[1]} ${brks[2]} faNamesPre.txt $4 $5 itp #blast result now in $5/$i_tp_$4_$5.txt
						searchSA1=$(sh searchRevAA.sh $5/"$i"_itp_${tpSE[0]}_${brks[1]}.txt $ALCut $ASimCut)
						if [ "$searchSA1" != "N" ]; then
							sh getAAQPos.sh $5/"$i"_itp_${tpSE[0]}_${tpSE[1]}.txt 1 #temp1.pos A/A
							sh getRevAAQPos.sh $5/"$i"_itp_${tpSE[0]}_${brks[1]}.txt 2 #temp2.pos A/-A
							java getAInts temp1.pos temp2.pos $ALCut >temp1a2.pos
							maxL=$(java getMaxL temp1a2.pos)
							if [ $maxL -ne 0 ]; then
								searchSA="Found,Length=$maxL"
							else
								searchSA="N"
							fi
						else
							searchSA="N"
						fi
					fi

					#update $searchDA, when there exist a pair of repeats on d
					if [ "$searchDA" != "N" ]; then
						brks=( $(java getNbD $nbLs $nbRs $gd) )
						sh blBtBlkSq.sh $2 $j $nbLd $sd ${brks[1]} ${brks[2]} faNamesPre.txt $4 $5 itp #blast result now in $5/$i_tp_$4_$5.txt
						searchDA1=$(sh searchRevAA.sh $5/"$j"_itp_"$sd"_${brks[1]}.txt $ALCut $ASimCut)
						if [ "$searchDA1" != "N" ]; then
							sh getAAQPos.sh $5/"$j"_itp_"$sd"_$ed.txt 1 #temp1.pos
							sh getRevAAQPos.sh $5/"$j"_itp_"$sd"_${brks[1]}.txt 2
							java getAInts temp1.pos temp2.pos $ALCut >temp1a2.pos
							maxL=$(java getMaxL temp1a2.pos)
							if [ $maxL -ne 0 ]; then
								searchDA="Found,Length=$maxL"
							else
								searchDA="N"
							fi
						else
							searchDA="N"
						fi
					fi
					echo "Step $k: ${tpSE[0]} through ${tpSE[1]} Inverted Transposition: Source:$searchSA Destination:$searchDA" >>$5/reportA.txt
					if [[ $searchSA == *"Found"* ]] || [[ $searchDA == *"Found"* ]]; then
						cItp=$cItp+1          
					fi
					if [[ $searchSA == *"Found"* ]] && [[ $searchDA == *"Found"* ]]; then
						cItpBoth=$cItpBoth+1
					fi
				done
			fi
			
			if [ $bi -ne 0 ]; then #block interchange cases
				echo "*Block Interchange: $bi step(s):" >>$5/reportA.txt
				#show steps
				tpBlks=( $(grep '^B' tbi.txt) )
				n=${#tpBlks[@]}
				for (( k=1; k<$n-1; k=k+2 )) #check each block involved transposition
				do
					declare -i bID1=${tpBlks[$k]}
                    tpSE1=( $(sh getTBIBlk.sh $bID1) )
					declare -i bID2=${tpBlks[$k+1]}
					tpSE2=( $(sh getTBIBlk.sh $bID2) )

					#check A,B,A,B on s
					nbs1=( $(java getNbD ${tpSE1[0]} ${tpSE1[1]} $gs) )
					nbs2=( $(java getNbD ${tpSE2[0]} ${tpSE2[1]} $gs) )
					if [ ${nbs1[0]} -ne 0 ] && [ ${nbs2[0]} -ne 0 ]; then
						sh blBtBlkSq.sh $2 $i ${nbs1[0]} ${nbs1[1]} ${nbs2[0]} ${nbs2[1]} faNamesPre.txt $4 $5 bi #blast result now in $5/$j_tbi_${tpSE[0]}_${tpSE[0]}.txt
						sh blBtBlkSq.sh $2 $i ${nbs1[2]} ${nbs1[3]} ${nbs2[2]} ${nbs2[3]} faNamesPre.txt $4 $5 bi 
						searchSA=$(sh searchAA.sh $5/"$i"_bi_${nbs1[1]}_${nbs2[0]}.txt $ALCut $ASimCut) #search A A on s
						searchSB=$(sh searchAA.sh $5/"$i"_bi_${nbs1[3]}_${nbs2[2]}.txt $ALCut $ASimCut) #search B B on s
					else
						searchSA="N"
						searchSB="N"
					fi

					if [ "$searchSA" == "N" ] || [ "$searchSB" == "N" ]; then
						searchSAB="N"
					else
						searchSAB="A-$searchSA,B-$searchSB"
					fi

					#check A,B,A,B on d
					nbd1=( $(java getNbD ${tpSE1[0]} ${tpSE1[1]} $gd) )
					nbd2=( $(java getNbD ${tpSE2[0]} ${tpSE2[1]} $gd) )
					if [ ${nbd1[0]} -ne 0 ] && [ ${nbd2[0]} -ne 0 ]; then
						sh blBtBlkSq.sh $2 $j ${nbd1[0]} ${nbd1[1]} ${nbd2[0]} ${nbd2[1]} faNamesPre.txt $4 $5 bi #blast result now in $5/$j_tbi_${tpSE[0]}_${tpSE[0]}.txt
						sh blBtBlkSq.sh $2 $j ${nbd1[2]} ${nbd1[3]} ${nbd2[2]} ${nbd2[3]} faNamesPre.txt $4 $5 bi 
						searchDA=$(sh searchAA.sh $5/"$j"_bi_${nbd1[1]}_${nbd2[0]}.txt $ALCut $ASimCut) #search A A on d
						searchDB=$(sh searchAA.sh $5/"$j"_bi_${nbd1[3]}_${nbd2[2]}.txt $ALCut $ASimCut) #search B B on d
					else
						searchDA="N"
						searchDB="N"
					fi

					if [ "$searchDA" == "N" ] || [ "$searchDB" == "N" ]; then
						searchDAB="N"
					else
						searchDAB="A-$searchDA,B-$searchDB"
					fi
					echo "Step $k: ${tpSE1[0]} through ${tpSE1[1]} with ${tpSE2[0]} through ${tpSE2[1]} Block Interchange: Source:$searchSAB Destination:$searchDAB" >>$5/reportA.txt
					if [[ $searchSAB == *"Found"* ]] || [[ $searchDAB == *"Found"* ]]; then
						cBi=$cBi+1          
					fi
					if [[ $searchSA == *"Found"* ]] && [[ $searchDA == *"Found"* ]]; then
						cBiBoth=$cBiBoth+1
					fi
				done
			fi
			
			if [ $ibi -ne 0 ]; then
				echo "*Inverted Block Interchange: $ibi step(s):" >>$5/reportA.txt
				##############################show block interchange steps
				tpBlks=( $(grep '^IB' tbi.txt) )
				n=${#tpBlks[@]}
				for (( k=1; k<$n-1; k=k+2 )) #check each block involved transposition
				do
					declare -i bID1=${tpBlks[$k]}
                    tpSE1=( $(sh getTBIBlk.sh $bID1) )
					declare -i bID2=${tpBlks[$k+1]}
					tpSE2=( $(sh getTBIBlk.sh $bID2) )

					#check A,B,A,B on s
					nbs1=( $(java getNbD ${tpSE1[0]} ${tpSE1[1]} $gs) )
					nbs2=( $(java getNbD ${tpSE2[0]} ${tpSE2[1]} $gs) )
					if [ ${nbs1[0]} -ne 0 ] && [ ${nbs2[0]} -ne 0 ]; then
						sh blBtBlkSq.sh $2 $i ${nbs1[0]} ${nbs1[1]} ${nbs2[2]} ${nbs2[3]} faNamesPre.txt $4 $5 ibi #blast result now in $5/$j_tbi_${tpSE[0]}_${tpSE[0]}.txt
						sh blBtBlkSq.sh $2 $i ${nbs1[2]} ${nbs1[3]} ${nbs2[0]} ${nbs2[1]} faNamesPre.txt $4 $5 ibi 
						searchSA=$(sh searchRevAA.sh $5/"$i"_ibi_${nbs1[1]}_${nbs2[2]}.txt $ALCut $ASimCut) #search A -A on s
						searchSB=$(sh searchRevAA.sh $5/"$i"_ibi_${nbs1[3]}_${nbs2[0]}.txt $ALCut $ASimCut) #search B -B on s
					else
						searchSA="N"
						searchSB="N"
					fi

					if [ "$searchSA" == "N" ] || [ "$searchSB" == "N" ]; then
						searchSAB="N"
					else
						searchSAB="A-$searchSA,B-$searchSB"
					fi

					#check A,B,A,B on d
					nbd1=( $(java getNbD ${tpSE1[0]} ${tpSE1[1]} $gd) )
					nbd2=( $(java getNbD ${tpSE2[0]} ${tpSE2[1]} $gd) )
					if [ ${nbd1[0]} -ne 0 ] && [ ${nbd2[0]} -ne 0 ]; then
						sh blBtBlkSq.sh $2 $j ${nbd1[0]} ${nbd1[1]} ${nbd2[2]} ${nbd2[3]} faNamesPre.txt $4 $5 ibi #blast result now in $5/$j_tbi_${tpSE[0]}_${tpSE[0]}.txt
						sh blBtBlkSq.sh $2 $j ${nbd1[2]} ${nbd1[3]} ${nbd2[0]} ${nbd2[1]} faNamesPre.txt $4 $5 ibi 
						searchDA=$(sh searchRevAA.sh $5/"$j"_ibi_${nbd1[1]}_${nbd2[2]}.txt $ALCut $ASimCut) #search A -A on d
						searchDB=$(sh searchRevAA.sh $5/"$j"_ibi_${nbd1[3]}_${nbd2[0]}.txt $ALCut $ASimCut) #search B -B on d
					else
						searchDA="N"
						searchDB="N"
					fi

					if [ "$searchDA" == "N" ] || [ "$searchDB" == "N" ]; then
						searchDAB="N"
					else
						searchDAB="A-$searchDA,B-$searchDB"
					fi
					echo "Step $k: ${tpSE1[0]} through ${tpSE1[1]} with ${tpSE2[0]} through ${tpSE2[1]} Inverted Block Interchange: Source:$searchSAB Destination:$searchDAB" >>$5/reportA.txt
					if [[ $searchSAB == *"Found"* ]] || [[ $searchDAB == *"Found"* ]]; then
						cIbi=$cIbi+1          
					fi
					if [[ $searchSA == *"Found"* ]] && [[ $searchDA == *"Found"* ]]; then
						cIbiBoth=$cIbiBoth+1
					fi
				done
			fi
			
			if [ $hbi -ne 0 ]; then
				echo "*Half-Inverted Block Interchange: $hbi step(s):" >>$5/reportA.txt
				##############################show block interchange steps
				tpBlks=( $(grep '^HIB' tbi.txt) )
				n=${#tpBlks[@]}
				for (( k=1; k<$n-1; k=k+2 )) #check each block involved
				do
					declare -i bID1=${tpBlks[$k]}
                    tpSE1=( $(sh getTBIBlk.sh $bID1) )
					declare -i bID2=${tpBlks[$k+1]}
					tpSE2=( $(sh getTBIBlk.sh $bID2) )
					
					
					nbs1=( $(java getNbD ${tpSE1[0]} ${tpSE1[1]} $gs) )
					nbs2=( $(java getNbD ${tpSE2[0]} ${tpSE2[1]} $gs) )
					if [ ${nbs1[0]} -ne 0 ] && [ ${nbs2[0]} -ne 0 ]; then
						sh blBtBlkSq.sh $2 $i ${nbs1[0]} ${nbs1[1]} ${nbs1[2]} ${nbs1[3]} faNamesPre.txt $4 $5 hbi #blast result now in $5/$j_tbi_${tpSE[0]}_${tpSE[0]}.txt
						sh blBtBlkSq.sh $2 $i ${nbs1[0]} ${nbs1[1]} ${nbs2[0]} ${nbs2[1]} faNamesPre.txt $4 $5 hbi 
						sh blBtBlkSq.sh $2 $i ${nbs1[0]} ${nbs1[1]} ${nbs2[2]} ${nbs2[3]} faNamesPre.txt $4 $5 hbi
						searchSA1=$(sh searchRevAA.sh $5/"$i"_hbi_${nbs1[1]}_${nbs1[2]}.txt $ALCut $ASimCut) #search A -A on s
						searchSA2=$(sh searchAA.sh $5/"$i"_hbi_${nbs1[1]}_${nbs2[0]}.txt $ALCut $ASimCut)
						searchSA3=$(sh searchRevAA.sh $5/"$i"_hbi_${nbs1[1]}_${nbs2[2]}.txt $ALCut $ASimCut)
						if [ "$searchSA1" != "N" ] && [ "$searchSA2" != "N" ] && [ "$searchSA3" != "N" ]; then
							sh getRevAAQPos.sh $5/"$i"_hbi_${nbs1[1]}_${nbs1[2]}.txt 1
							sh getAAQPos.sh $5/"$i"_hbi_${nbs1[1]}_${nbs2[0]}.txt 2
							sh getRevAAQPos.sh $5/"$i"_hbi_${nbs1[1]}_${nbs2[2]}.txt 3
							java getAInts temp1.pos temp2.pos $ALCut >temp1a2.pos
							java getAInts temp1a2.pos temp3.pos $ALCut >tempall.pos
							maxL=$(java getMaxL tempall.pos)
							if [ $maxL -ne 0 ]; then
								searchSA="Found,Length=$maxL"
							else
								searchSA="N"
							fi
						else
							searchSA="N"
						fi

					else
						searchSA="N"
					fi


					nbd1=( $(java getNbD ${tpSE1[0]} ${tpSE1[1]} $gd) )
					nbd2=( $(java getNbD ${tpSE2[0]} ${tpSE2[1]} $gd) )
					if [ ${nbd1[0]} -ne 0 ] && [ ${nbd2[0]} -ne 0 ]; then
						sh blBtBlkSq.sh $2 $j ${nbd1[0]} ${nbd1[1]} ${nbd1[2]} ${nbd1[3]} faNamesPre.txt $4 $5 hbi #blast result now in $5/$j_tbi_${tpSE[0]}_${tpSE[0]}.txt
						sh blBtBlkSq.sh $2 $j ${nbd1[0]} ${nbd1[1]} ${nbd2[0]} ${nbd2[1]} faNamesPre.txt $4 $5 hbi 
						sh blBtBlkSq.sh $2 $j ${nbd1[0]} ${nbd1[1]} ${nbd2[2]} ${nbd2[3]} faNamesPre.txt $4 $5 hbi
						searchDA1=$(sh searchRevAA.sh $5/"$j"_hbi_${nbd1[1]}_${nbd1[2]}.txt $ALCut $ASimCut) #search A -A on s
						searchDA2=$(sh searchAA.sh $5/"$j"_hbi_${nbd1[1]}_${nbd2[0]}.txt $ALCut $ASimCut)
						searchDA3=$(sh searchRevAA.sh $5/"$j"_hbi_${nbd1[1]}_${nbd2[2]}.txt $ALCut $ASimCut)
						if [ "$searchDA1" != "N" ] && [ "$searchDA2" != "N" ] && [ "$searchDA3" != "N" ]; then
							sh getRevAAQPos.sh $5/"$j"_hbi_${nbd1[1]}_${nbd1[2]}.txt 1
							sh getAAQPos.sh $5/"$j"_hbi_${nbd1[1]}_${nbd2[0]}.txt 2
							sh getRevAAQPos.sh $5/"$j"_hbi_${nbd1[1]}_${nbd2[2]}.txt 3
							java getAInts temp1.pos temp2.pos $ALCut >temp1a2.pos
							java getAInts temp1a2.pos temp3.pos $ALCut >tempall.pos
							maxL=$(java getMaxL tempall.pos)
							if [ $maxL -ne 0 ]; then
								searchDA="Found,Length=$maxL"
							else
								searchDA="N"
							fi
						else
							searchDA="N"
						fi
					else
						searchDA="N"
					fi
					echo "Step $k: ${tpSE1[0]} through ${tpSE1[1]} with ${tpSE2[0]} through ${tpSE2[1]} Half-Inverted Block Interchange: Source:$searchSA Destination:$searchDA" >>$5/reportA.txt
					if [[ $searchSA == *"Found"* ]] || [[ $searchDA == *"Found"* ]]; then
						cHib=$cHib+1          
					fi
					if [[ $searchSA == *"Found"* ]] && [[ $searchDA == *"Found"* ]]; then
						cHibBoth=$cHibBoth+1 
					fi
				done
			fi
		fi
		
		#The following deals with reversal steps
        if [ $distance -ne 0 ]; then
			echo "*Reversal: $distance step(s):" >>$5/reportA.txt	                  
		fi
		for (( k=1; k<$distance+1; k=k+1 )) #check reversal in each steps
		do
			LBlk="$(grep 'Step '$k': ' $5/sn"$i"To"$j".sn|cut -d '[' -f2|sed 's/\']'.*$//g')"
			RBlk="$(grep 'Step '$k': ' $5/sn"$i"To"$j".sn|cut -d '[' -f3|sed 's/\']'.*$//g')"
			
			if [ $LBlk -gt 0 ]; then
				LBlk=$(sh getleftBlk.sh $LBlk)
			else
				LBlk=$LBlk*-1
				LBlk=$(sh getrightBlk.sh $LBlk)
				LBlk=$LBlk*-1
			fi

			if [ $RBlk -gt 0 ]; then
				RBlk=$(sh getrightBlk.sh $RBlk)
			else
				RBlk=$RBlk*-1
				RBlk=$(sh getleftBlk.sh $RBlk)
				RBlk=$RBlk*-1
			fi
							 
			
			totRev=$totRev+1
			sh getleftSq.sh $2 $i $LBlk $RANGE faNamesPre.txt $4 >$5/sLeft$LBlk.fa 
			sh getRightSq.sh $2 $i $RBlk $RANGE faNamesPre.txt $4 >$5/sRight$RBlk.fa
			bl2seq -i $5/sLeft$LBlk.fa -j $5/sRight$RBlk.fa -p blastn -e 1e-50 -o sBlast.txt
			searchSA=$(sh findA.sh $LBlk $RBlk $2 $i sBlast.txt $ALCut $ASimCut)
				
						   
			sh getleftSq.sh $2 $j $LBlk $RANGE faNamesPre.txt $4 >$5/dLeft$LBlk.fa
			sh getRightSq.sh $2 $j $RBlk $RANGE faNamesPre.txt $4 >$5/dRight$RBlk.fa
			bl2seq -i $5/dLeft$LBlk.fa -j $5/dRight$RBlk.fa -p blastn -e 1e-50 -o dBlast.txt
			searchDA=$(sh findA.sh $LBlk $RBlk $2 $j dBlast.txt $ALCut $ASimCut)
			

			echo "Step $k: $LBlk through $RBlk Reversal: Source:$searchSA Destination:$searchDA" >>$5/reportA.txt

			if [[ $searchSA == *"Found"* ]] || [[ $searchDA == *"Found"* ]]; then
				count=$count+1                   
			fi

			if [[ $searchSA == *"Found"* ]] && [[ $searchDA == *"Found"* ]]; then
				countBoth=$countBoth+1
			fi

		done
      

    done
done

if [ $totTp -ne 0 ]; then
	echo "Total $totTp transpositions, and $cTp of them have A/A/A. In $cTpBoth transpositions, A/A/A occurs in both genomes." >>$5/reportA.txt
fi
if [ $totItp -ne 0 ]; then
	echo "Total $totItp inverted transpositions, and $cItp of them have A/A/-A. In $cItpBoth inverted transpositions, A/A/-A occurs in both genomes." >>$5/reportA.txt
fi
if [ $totBi -ne 0 ]; then
	echo "Total $totBi block interchanges, and $cBi of them have A/B/A/B. In $cBiBoth block interchanges, A/B/A/B occurs in both genomes." >>$5/reportA.txt
fi
if [ $totIbi -ne 0 ]; then
	echo "Total $totIbi inverted block interchanges, and $cIbi of them have A/B/-B/-A. In $cIbiBoth inverted block interchanges, A/B/-B/-A occurs in both genomes." >>$5/reportA.txt
fi
if [ $totHib -ne 0 ]; then
	echo "Total $totHib half-inverted block interchanges, and $cHib of them have A/-A/A/-A. In $cHibBoth half-inverted block interchanges, A/-A/A/-A occurs in both genomes." >>$5/reportA.txt
fi
echo "Total $totRev reversals, and $count of them have A/-A. $countBoth of A/-A occurs in both genomes." >>$5/reportA.txt
#rm $5/*.fa
rm $5/*.sn