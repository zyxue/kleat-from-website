#!/bin/bash
#$ -S /bin/bash
#$ -l mem_token=30G,mem_free=30G,h_vmem=30G
#$ -m a
#$ -M traymond@bcgsc.ca
#$ -j y
#$ -o log/$JOB_NAME.$JOB_ID.$TASK_ID

set -eux

lib=$JOB_NAME
sup=$SGE_TASK_ID
rpet=$1

mkdir -p $lib/

date

### cufflinks
./cssupport <(cat ../$lib/wgEncodeGisRnaPet*CellPapPlusRawRep*.bg) <(awk '{print "chr" $0 "\t+"}' ../../cufflinks/$lib/transcripts_+.bg) + $sup $rpet >$lib/cuff.1_${sup}_${rpet}+
./cssupport <(cat ../$lib/wgEncodeGisRnaPet*CellPapMinusRawRep*.bg) <(awk '{print "chr" $0 "\t-"}' ../../cufflinks/$lib/transcripts_-.bg) - $sup $rpet >$lib/cuff.1_${sup}_${rpet}-

cat $lib/cuff.1_${sup}_${rpet}[+-] >$lib/cuff.1_${sup}_${rpet}

FP=`grep -v TP= $lib/cuff.1_${sup}_${rpet} |awk '$6 == 0' \
   |perl mergecoordinates_long.pl \
   |wc -l`
TP=`grep -c TP=1 $lib/cuff.1_${sup}_${rpet}`
FN=`grep -c TP=0 $lib/cuff.1_${sup}_${rpet}`
TPR=`/home/traymond/bin/calc "$TP/($TP + $FN)"`
FDR=`/home/traymond/bin/calc "$FP/($TP + $FP)"`
echo $lib cuff $TP $FP $FN $TPR $FDR >$lib/cuff.1_${sup}_${rpet}.result.txt

date

## CLEAT
./cssupport <(cat ../$lib/wgEncodeGisRnaPet*CellPapPlusRawRep*.bg) $lib.CLEAT.polyA.tsv + $sup $rpet >$lib/CLEAT.1_${sup}_${rpet}+
./cssupport <(cat ../$lib/wgEncodeGisRnaPet*CellPapMinusRawRep*.bg) $lib.CLEAT.polyA.tsv - $sup $rpet >$lib/CLEAT.1_${sup}_${rpet}-


cat $lib/CLEAT.1_${sup}_${rpet}[+-] >$lib/CLEAT.1_${sup}_${rpet}

FP=`grep -v TP= $lib/CLEAT.1_${sup}_${rpet} |awk '$6 == 0' \
   |perl mergecoordinates_long.pl \
   |wc -l`
TP=`grep -c TP=1 $lib/CLEAT.1_${sup}_${rpet}`
FN=`grep -c TP=0 $lib/CLEAT.1_${sup}_${rpet}`
TPR=`/home/traymond/bin/calc "$TP/($TP + $FN)"`
FDR=`/home/traymond/bin/calc "$FP/($TP + $FP)"`
echo $lib cuff $TP $FP $FN $TPR $FDR >$lib/CLEAT.1_${sup}_${rpet}.result.txt

date

## Trinity + CLEAT
./cssupport <(cat ../$lib/wgEncodeGisRnaPet*CellPapPlusRawRep*.bg) $lib.trin.polyA.tsv + $sup $rpet >$lib/trin.1_${sup}_${rpet}+
./cssupport <(cat ../$lib/wgEncodeGisRnaPet*CellPapMinusRawRep*.bg) $lib.trin.polyA.tsv - $sup $rpet >$lib/trin.1_${sup}_${rpet}-

cat $lib/trin.1_${sup}_${rpet}[+-] >$lib/trin.1_${sup}_${rpet}

FP=`grep -v TP= $lib/trin.1_${sup}_${rpet} |awk '$6 == 0' \
   |perl mergecoordinates_long.pl \
   |wc -l`
TP=`grep -c TP=1 $lib/trin.1_${sup}_${rpet}`
FN=`grep -c TP=0 $lib/trin.1_${sup}_${rpet}`
TPR=`/home/traymond/bin/calc "$TP/($TP + $FN)"`
FDR=`/home/traymond/bin/calc "$FP/($TP + $FP)"`
echo $lib cuff $TP $FP $FN $TPR $FDR >$lib/trin.1_${sup}_${rpet}.result.txt

date
