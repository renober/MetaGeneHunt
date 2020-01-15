#!/bin/bash
###################################################
clear
echo "MetaGeneHunt_V1 by Renaud Berlemont (Renoberlemont@gmail.com)"
echo "MetaGeneHunt for precise domain-specific GH identification in MG-RAST annotated metagenome"
echo "using custom domain-specific reannotation of the M5nr db using GeneHunt \(see"
echo "Nguyen SN, Flores A, Talamantes D, Dar F, Valdez A, Schwans J, Berlemont R."
echo "GeneHunt for rapid domain-specific annotation of glycoside hydrolases. Sci Rep."
echo "2019 Jul 12;9(1):10137. doi: 10.1038/s41598-019-46290-w. PubMed PMID: 31300677;"
echo "PubMed Central PMCID: PMC6626019.)"
echo "###################################################"
echo INPUT file \(e.g.\ Input_file.txt\)\:
read INPUT
echo Overlap between local alignment and predicted domain \(positive number\, e\.g\.\:20\)\:
read z

rm -f .temp.txt .temp2.txt .phit_superblat.txt .xx .sum_aa90.txt .RAT_id.txt .ANNOT
cut -f 3 RAT_GH.txt | sort -k 1 | uniq | grep -v "^$" > .RAT_id.txt   #RAT_GH is a custom Reference Annotation Table with GH domain specific annotation of the M5nr db.
timestamp=$(date +"%y_%m_%d_%Hh_%Mm")
echo -e MGMid"\t"QUERYid"\t"TARGETid"\t"ali_start"\t"ali_end"\t"Nseq"\t"Z_overlap_cutoff"\t"PFam >ANNOT_$timestamp.txt

##### for each MGM
#1.1#####downloading the mgm files!
while read mgmlist
do
mgm=$(echo $mgmlist | cut -f 1 | cut -d " " -f 1 )
echo $mgm
curl -X GET http://api.metagenomics.anl.gov/1/download/$mgm.3?file=550.1 -o $mgm.550  #downloading the files 550 (small)
curl -X GET http://api.metagenomics.anl.gov/1/download/$mgm.3?file=650.1 -o $mgm.650 #downloading the files 550 (big)

#1.2#summarizing the aa90.mapping file
awk '{print $1"\t"$NF}' < $mgm.550 | sed 's/\%/\t/g' | awk '{print $1"\t"NF}' |sort -k 1  >.sum_aa90.txt  
#1.3 extracting potential hits (phit)
echo summarizing $mgm.650 is going to take some time ... be patient
grep -h -F -f .RAT_id.txt < $mgm.650  |awk '!x[$1]++' >.xx  
echo xx is done
cut -f 1,2,9,10 <.xx | grep -h -F -f .RAT_id.txt | sort -k 1 > .phit_superblat.txt 

join .phit_superblat.txt .sum_aa90.txt -1 1 -2 1 -a 1 -o '1.1 1.2 1.3 1.4 2.2' -e 1 | tr ' ' '\t' >.temp.txt #Quick Merge (!): $1-SeqID, $2-md5id, $3-ali.start, $4-ali-end, $5-seq.count

#2#merging phit with aa90 cluster numbers
#awk 'NR==FNR {h[$1] = $2; next} {print $1,$2,$3,$4,h[$1]}'  .sum_aa90.txt .phit_superblat.txt | sed 's/ $/ 1/g'  | tr ' ' '\t' > .temp.txt 

echo Matching the alignments with domains ...
#####################################################
z=20
while read line
do
md5id=$(echo $line | cut -d " " -f 2)
blat_start=$(echo $line | cut -d " " -f 3)   
blat_end=$(echo $line | cut -d " " -f 4) 		
grep $md5id RAT_GH.txt >.temp2.txt # $1-PF, $2-PFLen, $3:md5id, $4:tar.Len, $5:Eval, $6:ali.start.dom, $7:ali.end.dom, $8:ali.start.targ, $9:ali.end.targ, $10:%cov.
while read HIT ; do #from .temp2.txt
rm .ANNOT
dom_start=$(echo $HIT | cut -d " " -f 8) 
dom_end=$(echo $HIT | cut -d " " -f 9) 
dom_PF=$(echo $HIT | cut -d " " -f 1)
if  [[ "$blat_end-$dom_start" -ge "z" &&  "$blat_end" -le "$dom_end" ]] #cond1: the alignment is shited toward the N-terminal end of the domain, but overlap gt z or completely overlap with the domain
then 
echo $mgm $line $z $dom_PF | tr ' ' '\t' >>.ANNOT
elif  [[ "$dom_end-$blat_start" -ge "z" && "$blat_start" -ge "$dom_start" ]] #cond2: the alignment is shited toward the C-terminal end of the domain, but overlap gt z or completely overlap with the domain
then
echo $mgm $line $z $dom_PF | tr ' ' '\t' >>.ANNOT
else 
echo $mgm $line $z "NO_DOMAIN_MATCH" | tr ' ' '\t' >>.ANNOT
fi 
done <.temp2.txt
sort -k 7 <.ANNOT | tail -n 1 | tee -a ANNOT_$timestamp.txt
done < .temp.txt
rm -f  $mgm.550 $mgm.650
#####################################################

done < $INPUT

rm -f .temp.txt .temp2.txt .phit_superblat.txt .xx .sum_aa90.txt .RAT_id.txt .ANNOT
