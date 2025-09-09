for i in *segments.csv;do echo grep -v "neutral" $i\>$i.noN;done> noneutral.sh
step1_mod.r
 for i in *.txt.mod.txt;do echo tail -n +2 $i \|cut -f2- \> $i.bed;done>step2.pbs
for i in *.noN;do echo tail -n +2 $i \|cut -f2- \> $i.bed;done>step2.pbs
 
for i in *.txt.bed;do echo intersectBed -a $i -b hg19.ref.bed -wa -wb -F 0.5 \>$i.ann.bed;done >step3.ann.50.pbs
for i in *.txt.bed;do echo intersectBed -a $i -b hg19.ref.bed -wa -wb -F 1 \>$i.ann.bed;done >step3.ann.100.pbs

 for i in *.noN.bed;do echo intersectBed -a $i -b hg19.ref.bed -wa -wb -F 1\>$i.ann.bed;done >step3.ann.100.pbs
