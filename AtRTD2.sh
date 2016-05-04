##################AtRTD2 analysis
##################April 19, 2016

##########################################################################
#### 1. merge two transcript annotations ###############################
#########################################################################

###this part of code is to merge Araport11 and the reference transcript dataset we generated. It could be adapted to merge two transcript annotations to remove redundancies and low quality transcripts as specified in the manuscript providing the same genome reference has been used.
 

###araport11
#Araport11_genes.20151202.gff3

###the reference transcript data we generated
##mergedatRTDClockYTranscript.gff

###attach "_P" to transcript names to distinguish it from transcripts from the other reference dataset 
cat Araport11_genes.20151202.gff3 | sed 's/%26#64257%3B/fi/g' | awk -F"\t" '{if($9~/\./){gsub(/\./, "_P", $9)};print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' > Araport11_genes.20151202_1.gff3

###add intron coordinates to both 
gt gff3 -sort yes -retainids yes -tidy yes -checkids yes -addintrons yes Araport11_genes.20151202_1.gff3 | grep -v "#" > Araport11_genes.20151202_intron.gff

gt gff3 -sort yes -retainids yes -tidy yes -checkids yes -addintrons yes mergedatRTDClockYTranscript.gff | grep -v "#" > mergedatRTDClockYTranscript_intron.gff

# #2. gene models
# ###33,635 genes
 cat mergedatRTDClockYTranscript.gff | awk '{if($3=="gene"){print}}' > mergedatRTDClockYTranscript_gene
# ###27,667 genes
# #cat /home/rz41242/Projects/2014AS/TranscriptAnnot/Marquez12/mergedTranscript3.gff |awk '{if($3=="gene"){split($9, a, "=");print a[4]}}'| sort | uniq > atRTD_gene
 cat Araport11_genes.20151202_intron.gff |awk '{if($3=="gene"){print}}'> Araport11_genes.20151202_intron_gene

# #1)compare the gene models with the same name, remove the smaller one being contained in a different gene model; 2) if the two gene models only partially overlap, keep the gene model in TAIR10. 
 cat Araport11_genes.20151202_intron_gene mergedatRTDClockYTranscript_gene| sort -k9  | awk 'BEGIN{crd="R"}{split($9, b, ";");split(b[1], c, "=");if(crd==c[2]&&($4>St&&$5<=Ed||$4>=St&&$5<Ed)){print id"\t"crd"\t"St"\t"Ed"\t"$2"\t"c[2]"\t"$4"\t"$5}else{if(crd==c[2]){print $2"\t"c[2]"\t"$4"\t"$5"\t"id"\t"crd"\t"St"\t"Ed}};crd=c[2];id=$2;St=$4;Ed=$5}'>geneToRemove

# ####if the longer gene overlap with its neighboring gene, use the shorter gene to replace the longer gene;if replacement of the longer with the shorter for the first gene is not enough, replace the shorter gene with longer for the second gene too.
###if two genes are of equal length, keep tair10
 cat geneToRemove | awk '{if(NR==1){a=$1; b=$2; l=$3;r=$4;c=$5; d=$6; e=$7;f=$8}else{if($3>r){print a"\t"b"\t"l"\t"r"\t"c"\t"d"\t"e"\t"f;a=$1; b=$2; l=$3;r=$4;c=$5; d=$6; e=$7;f=$8}else{if($3>=f){print c"\t"d"\t"e"\t"f"\t"a"\t"b"\t"l"\t"r;a=$1; b=$2; l=$3;r=$4;c=$5; d=$6; e=$7;f=$8}else{if($7>=r){print a"\t"b"\t"l"\t"r"\t"c"\t"d"\t"e"\t"f;a=$5; b=$6; l=$7;r=$8;c=$1; d=$2; e=$3;f=$4}else{if(f<r){print c"\t"d"\t"e"\t"f"\t"a"\t"b"\t"l"\t"r}else{print a"\t"b"\t"l"\t"r"\t"c"\t"d"\t"e"\t"f};if($7>$3){a=$5; b=$6; l=$7;r=$8;c=$1; d=$2; e=$3;f=$4}else{a=$1; b=$2; l=$3;r=$4;c=$5; d=$6; e=$7;f=$8}}}}}}END{if($3>r){print a"\t"b"\t"l"\t"r"\t"c"\t"d"\t"e"\t"f;}else{if($3>=f){print c"\t"d"\t"e"\t"f"\t"a"\t"b"\t"l"\t"r;}else{if($7>=r){print a"\t"b"\t"l"\t"r"\t"c"\t"d"\t"e"\t"f;}else{if(f<r){print c"\t"d"\t"e"\t"f"\t"a"\t"b"\t"l"\t"r}else{print a"\t"b"\t"l"\t"r"\t"c"\t"d"\t"e"\t"f}}}}}'> geneToRemove1

 cat geneToRemove1 | awk '$5!~/Araport11/{print $6}'> mergedatRTDClockYGene2Remove

 cat geneToRemove1 | awk '$5~/Araport11/{print $6}'> Araport11Gene2Remove

 ####remove the redundant gene models

 cat /home/rz41242/Projects/2014AS/Merge/RTDv1ClockYamile/mergedatRTDClockYTranscript.gff | awk 'BEGIN{while(getline<"mergedatRTDClockYGene2Remove"){a[$1]=5}}{split($9, b, ";");split(b[1], c, "="); if(a[c[2]]!=5){print}}' > temp1aa

 cat Araport11_genes.20151202_1.gff3 | awk 'BEGIN{while(getline<"Araport11Gene2Remove"){a[$1]=5}}{split($9, b, ";");split(b[1], c, "="); if(a[c[2]]!=5){print}}' > temp2aa

# ## 34,371 genes, 123,925 (48,389 + 75,536) transcripts
 cat temp1aa temp2aa > mergedGene.gff

# #####get the transcript proportions relative to the gene length
  cat mergedGene.gff | awk '$3~/gene/{split($9, a, ";"); split(a[1], b, "="); print b[2]"\t"$5-$4+1}' > GeneLength
  cat mergedGene.gff | awk '$3=="transcript"||$3=="mRNA"||$3~/miRNA/||$3~/mRNA_TE_gene/||$3~/ncRNA/||$3~/rRNA/||$3~/snoRNA/||$3~/snRNA/||$3~/tRNA/{split($9, a, ";"); split(a[1], b, "="); print b[2]"\t"$5-$4+1}'> TranscriptLength
  cat TranscriptLength | awk 'BEGIN{while(getline<"GeneLength"){a[$1]=$2}}{split($1, b, "."); split(b[1], c, "_"); print $1"\t"$2/a[c[1]]}' > TranscriptProportion

 ##take the intron coordinates for all RNAs with exons (excluding pseudogenes)
cat /home/rz41242/Projects/2014AS/Merge/RTDv1ClockYamile/mergedatRTDClockYTranscript_intron.gff| awk '{split($9, a, ";"); split(a[1], b, "="); if($3=="transcript"||$3=="mRNA"||$3~/miRNA/||$3~/mRNA_TE_gene/||$3~/ncRNA/||$3~/rRNA/||$3~/snoRNA/||$3~/snRNA/||$3~/tRNA/){printf "\n"; trpt=b[2]; printf trpt"\t"$3"\t"$4"\t"$5"\t"}; if($3~/intron/&&b[2]==trpt){printf $3"\t"$4"\t"$5"\t"}}END{printf "\n"}' | sed '/^$/d'>  mergedatRTDClockYTranscript_sorted_intron1
 
#4.take the intron coordinates for all RNAs with exons (excluding pseudogenes)
cat Araport11_genes.20151202_intron.gff | awk '{split($9, a, ";"); split(a[1], b, "="); if($3~/transcript/||$3~/mRNA/||$3~/miRNA/||$3~/mRNA_TE_gene/||$3~/ncRNA/||$3~/rRNA/||$3~/snoRNA/||$3~/snRNA/||$3~/tRNA/){printf "\n"; trpt=b[2]; printf trpt"\t"$3"\t"$4"\t"$5"\t"}; if($3~/intron/&&b[2]==trpt){printf $3"\t"$4"\t"$5"\t"}}END{printf "\n"}' | sed '/^$/d'> Araport11_genes20151202_sorted_intron1

###keep longer transcripts and remove the shorter transcripts 
##sort by transcript by intron coordinates
cat mergedatRTDClockYTranscript_sorted_intron1 Araport11_genes20151202_sorted_intron1 | sort -k3,3n | awk '{if(/intron/){a="intron";for(i=6; i<=NF; i++){a=a"-"$i};print $1"\t"$2"\t"$3"\t"$4"\t"a"\t"$4-$3}else{print $0"\t"$4-$3}}'| sort -k5> temp


#### remove transcripts longer than the gene length (118,642 transcripts)
cat temp | awk 'BEGIN{while(getline<"TranscriptProportion"){a[$1]=$2}}a[$1]<=1{print}'> temp0
###5,283 transcripts removed
cat temp | awk 'BEGIN{while(getline<"TranscriptProportion"){a[$1]=$2}}a[$1]>1{print}'> transcriptsToRemove0

# remove transcripts shorter than the 70% of the gene length (115,359 transcripts)
cat temp0 | awk 'BEGIN{while(getline<"TranscriptProportion"){a[$1]=$2}}a[$1]>0.7{print}'> temp01
###3,283 transcripts removed
cat temp0 | awk 'BEGIN{while(getline<"TranscriptProportion"){a[$1]=$2}}a[$1]<=0.7{print}'> transcriptsToRemove01

####remove redundant single exon transcripts 
##GeneName:gn
##Transcript coordinates Start End: TCS TCE
##Transcript Name: tn
##transcript length: tl

###4,916 (4801 non-redundant) transcripts removed
cat temp01 | grep -v intron | sort | awk '{split($1, a, "."); split(a[1], b, "_");if(NR==1){i=1; tl[i]=$5; tn[i]=$1; gn=b[1]; TCS[i]=$3; TCE[i]=$4}else{if(gn==b[1]){i++; tl[i]=$5; tn[i]=$1; TCS[i]=$3; TCE[i]=$4}else{for(j=2;j<=i;j++){for(k=1;k<j;k++){if((TCE[k]>TCS[j]&&TCE[k]<=TCE[j])||(TCE[j]>TCS[k]&&TCE[j]<=TCE[k])){if(tl[j]<=tl[k]){print tn[j]}else{print tn[k]}}}};delete TCS; delete TCE; i=1; tl[i]=$5; tn[i]=$1; gn=b[1]; TCS[i]=$3; TCE[i]=$4}}}END{for(j=2;j<=i;j++){for(k=1;k<j;k++){if((TCE[k]>TCS[j]&&TCE[k]<TCE[j])||(TCE[j]>TCS[k]&&TCE[j]<TCE[k])){if(tl[j]<=tl[k]){print tn[j]}else{print tn[k]}}}}}'> transcriptsToRemoveSingleExon

#### compare the transcripts with same introns and keep the longest transcript (70,757 transcripts)
cat temp01 | grep intron | awk '{if(NR==1){Int=$5;Lmax=$6;p=$0}else{if(Int==$5){if($6>Lmax){Lmax=$6;p=$0}}else{print p;Int=$5;Lmax=$6;p=$0}}}END{print p}'> temp1
###27,820 transcripts removed
cat temp01 | grep intron | awk 'BEGIN{while(getline<"temp1"){a[$0]=1}}a[$0]!=1{print}'> transcriptsToRemove1



###add in gene name
cat temp1 | awk '{split($1,a,"."); split(a[1], b, "_"); print b[1]"\t"$0}'| sort > temp2
####transcripts which whose intron coordinates are subset of other transcript intron coordinates are removed
##IC intron cooridnates
##ps transcript start position
##pe transcript end position
##g gene

cat temp2 | awk '{if(NR==1){i=1;IC[i]=$6;ps[i]=$4;pe[i]=$5;flag=0;g=$1}else{if(g==$1){for(j=1;j<=i;j++){if($6~IC[j]&&$6!=IC[j]){PosSt=match($6,IC[j]);PosEnd=PosSt+length(IC[j])-1; cd1=split(substr($6,1,PosSt-1), array1, "-"); split(substr($6,PosEnd,length($6)-PosEnd), array2, "-");if((PosSt<3||ps[j]>array1[cd1-1] )&&(PosEnd>=length($6)||pe[j]<array2[3])){IC[j]=$6;ps[j]=$4; pe[j]=$5;flag=1; delete array1; delete array2}};if(IC[j]~$6){PosSt=match(IC[j],$6);PosEnd=PosSt+length($6)-1; cd1=split(substr(IC[j],1,PosSt-1), array1, "-"); split(substr(IC[j],PosEnd,length(IC[j])-PosEnd), array2, "-");if((PosSt<3||$4>array1[cd1-1])&&(PosEnd>=length(IC[j])||$5<array2[3])){flag=1;break; delete array1; delete array2}}};if(flag==0){i=i+1;IC[i]=$6;ps[i]=$4; pe[i]=$5}else{flag=0}}else{for(j=1; j<=i;j++){print IC[j]"\t"ps[j]"\t"pe[j]};g=$1;delete IC; delete ps; delete pe; i=1; IC[i]=$6; ps[i]=$4; pe[i]=$5}}}END{for(j=1; j<=i;j++){print IC[j]"\t"ps[j]"\t"pe[j]}}'| sort | uniq > IntronToKeep1

####516 transcripts removed
cat temp1 | awk 'BEGIN{while(getline<"IntronToKeep1"){a[$1]=5;}}/intron/{if(a[$5]!=5){print}}'>transcriptsToRemove2
###70,241 transcripts
cat temp2 | awk 'BEGIN{while(getline<"IntronToKeep1"){a[$1]=5;}}/intron/{if(a[$6]==5){print}}'>temp3

#####remove transcript that skip 4 exons or retain 4 introns


cat temp3 | sed -e 's/intron-//g' | awk '{$6=$4"-"$6"-"$5; print}'> ForR.txt
#####
R

setwd("/your/directory")
 
rm(list=ls())  
##
x <- read.table(file="ForR.txt",sep = " ", header=F, as.is=T)

Dist<-9
Pa<-NA
for (i in unique(x[,1])){
##i<-"AT1G01010"
 Ind<-which(x[,1]==i);
 print(Ind);
	for (j in 2:length(Ind)){
		a=unlist(strsplit(x[Ind[j],6], "-"));
		for (k in 1:j){
		b=unlist(strsplit(x[Ind[k],6], "-"));
		d<-sort(c(a,b))
		e<-match(a,d)
		f<-e[2:length(e)]-e[1:length(e)-1];
		if(any(f>Dist)){g=c(x[Ind[j],2], x[Ind[k],2]);if(is.na(Pa)){Pa<-g}else{Pa<-rbind(g, Pa)}}
		e<-match(b,d)
		f<-e[2:length(e)]-e[1:length(e)-1];
		if(any(f>Dist)){g=c(x[Ind[k],2], x[Ind[j],2]);if(is.na(Pa)){Pa<-g}else{Pa<-rbind(g, Pa)}}
		}
	}
}

write.table(Pa, file = "Transcript2Remove.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep ="\t")
##################################################################################################################

##remove 32 transcripts
cat Transcript2Remove.txt | cut -f1 | sort | uniq > transcriptsToRemove3

 
### 41,735 (41,850 redundant) transcripts removed
cat transcriptsToRemoveSingleExon transcriptsToRemove0 transcriptsToRemove01 transcriptsToRemove1 transcriptsToRemove2 transcriptsToRemove3 > transcriptsToRemove



####remove all transcripts failed the filtering criteria 
cat mergedGene.gff|  awk '{split($9, b, ";"); split(b[2], e, "="); split(e[2], f, ","); for (i in f){a=$9; gsub(e[2],f[i],$9);print;$9=a}}'| awk 'BEGIN{while(getline<"transcriptsToRemove"){a[$1]=1}}{split($9, b, ";"); split(b[1], c, "="); split(c[2], d, ",");split(b[2], e, "="); if(a[d[1]]!=1&&a[e[2]]!=1){print}}'| awk '{if($3=="gene"){split($9, a, ";");$9=a[1]";"a[2]};print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}'| grep -v "#"> mergedatRTDClockYAraport.gff


##### sort and generate the final merged transcript annotation dataset (82,190 transcripts, 30,538 transcripts from Araport11, 33,674 (18,801 from TAIR10) from atTRD)
gt gff3 -sort yes -retainids yes -tidy yes -checkids yes -addids no  mergedatRTDClockYAraport.gff | grep -v "#" >  mergedatRTDClockYAraportTranscript.gff


##########################################################################
#### 2. Generate padded transcript reference in fasta and gtf format #####
##########################################################################

###This part of code is to generate padded transcript reference in fasta and gtf format

##mergedatRTDClockYAraportTranscript.gff contains the original Transcripts


####Get transcript start and end coordinates 
cat mergedatRTDClockYAraportTranscript.gff | awk '$3=="transcript"||$3=="mRNA"||$3~/miRNA/||$3~/mRNA_TE_gene/||$3~/ncRNA/||$3~/rRNA/||$3~/snoRNA/||$3~/snRNA/||$3~/tRNA/{split($9, a, ";"); split(a[1], b, "="); print b[2]"\t"$4"\t"$5"\t"$7}'> TranscriptCoordinates

####get the coordinates that covers the longest region in a gene
less TranscriptCoordinates | awk '{split($1, a, "[._]"); if(NR==1){GeneName=a[1];st=$2; ed=$3}else{if(a[1]==GeneName){if(st>$2){st=$2}; if(ed<$3){ed=$3}}else{print GeneName"\t"st"\t"ed;GeneName=a[1];st=$2;ed=$3}}}END{print GeneName"\t"st"\t"ed}'> GeneMaxCoordinates

####change start and end coordinates for every transcript to the padded coordinates
less TranscriptCoordinates | awk 'BEGIN{while(getline<"GeneMaxCoordinates"){a[$1]=$2;c[$1]=$3}}{split($1, b, "[._]");$3=c[b[1]];$2=a[b[1]];print $1"\t"$2"\t"$3"\t"$4}'> TranscriptCoordinates_Max

cat mergedatRTDClockYAraportTranscript.gff | awk 'BEGIN{while(getline<"TranscriptCoordinates_Max"){a[$1]=$2;b[$1]=$3}}{FS="\t";if($3~/transcript/||$3~/mRNA/||$3~/miRNA/||$3~/mRNA_TE_gene/||$3~/ncRNA/||$3~/rRNA/||$3~/snoRNA/||$3~/snRNA/||$3~/tRNA/){split($9, c, ";");split(c[1], d, "="); st=$4; ed=$5; $4=a[d[2]];$5=b[d[2]];print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9};if($3~/exon/){split($9, c, ";");split(c[2], d, "=");if($4>=a[d[2]]&&$4==st){$4=a[d[2]];print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9};if($4>st&&$5<ed){print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9};if($5==ed&&$5<=b[d[2]]){$5=b[d[2]];print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}}}'| uniq > atRTDv2_padding_19April2016.gff

####generated  fasta for padded transcripts
gffread -w atRTDv2_padding_19April2016.fa -g /home/rz41242/Projects/2014AS/TranscriptAnnot/Marquez12/genome/ArabidopsisGenome.fa atRTDv2_padding_19April2016.gff
########generated  gtf for padded transcripts
gffread -T -O -F atRTDv2_padding_19April2016.gff -o atRTDv2_padding_19April2016.gtf
