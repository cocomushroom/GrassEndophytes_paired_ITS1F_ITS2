#############################
### ITS workflow Feb 2019 ###
#############################
## Endophytes Florencia dataset##
## Sequenced by ITS1F and ITS2##
### https://github.com/benjjneb/ITS-Workflow/blob/master/ITS_workflow.md ####
### https://astrobiomike.github.io/amplicon/workflow_ex ###
## 5' -> 3' ITS1F primer: CTTGGTCATTTAGAGGAAGTAA
## Reverse complement ITS1F: TTACTTCCTCTAAATGACCAAG
## 5'-> 3' ITS2 primer: GCTGCGTTCTTCATCGATGC
## Reverse complement ITS2 primer: GCATCGATGAAGAACGCAGC



mkdir cut
touch log
for file in `find . -name "*R1*fastq.gz"`; do fastaFile=${file}; cutadapt -g CTTGGTCATTTAGAGGAAGTAA -a GCATCGATGAAGAACGCAGC -o cut/${fastaFile}_trimmed.fastq.gz ${fastaFile} >>log; done

for file in `find . -name "*R2*fastq.gz"`; do fastaFile=${file}; cutadapt -g GCTGCGTTCTTCATCGATGC -a TTACTTCCTCTAAATGACCAAG -o cut/${fastaFile}_trimmed.fastq.gz ${fastaFile} >>log; done


for i in .; do /./Users/kc178/bin/bin/FastQC.app/Contents/MacOS/fastqc *; done



maxee 2,6 truc 11

      input filtered merged nonchim
SD100 34370    30185  28161   28161
SD101 34745    29997  28586   28586
SD99  39504    34570  25513   25513


maxee 3,6 truc 12
      input filtered merged nonchim
SD100 34370    30181  28164   28164
SD101 34745    29995  28585   28585
SD99  39504    34566  25514   25514


maxee 5,7 truc 13
      input filtered merged nonchim
SD100 34370    30181  28165   28165
SD101 34745    29995  28585   28585
SD99  39504    34566  25513   25513
SD99  39504    34566  25514   25514