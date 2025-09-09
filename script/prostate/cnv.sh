
step1_mod.r
 for i in *.txt.mod.txt;do echo tail -n +2 $i \|cut -f2- \> $i.bed;done>step2.pbs
for i in *.txt.bed;do echo intersectBed -a $i -b hg19.ref.bed -wa -wb -F 0.5 \>$i.ann.bed;done >step3.ann.50.pbs
for i in *.txt.bed;do echo intersectBed -a $i -b hg19.ref.bed -wa -wb -F 1 \>$i.ann.bed;done >step3.ann.100.pbs
