#!/bin/sh

SEQUENCING_RUN="HU000"
I7_INDEX="B41 B42 B43 B44 B45 B46 B47 B48"
SAMPLE="YF000a"
SAMPLE_ALL="YF000"

mkdir -p ZUMIS/multx ZUMIS/merge
cd ZUMIS

## multx
cd multx
for name in ${I7_INDEX}
do
fastq-multx -m 0 -d 2 -x -b -B ~/${SEQUENCING_RUN}/${SEQUENCING_RUN}-${name}_index \
~/${SEQUENCING_RUN}/I7-${name}_1.fq.gz ~/${SEQUENCING_RUN}/I7-${name}_2.fq.gz -o %_R1.fq.gz %_R2.fq.gz \
> ${SEQUENCING_RUN}-${name}_multx.txt
done
cd ..

## merge reads
cd merge
for name in ${SAMPLE}
do
seqkit replace -p "\s.+" ~/${SEQUENCING_RUN}/ZUMIS/multx/${name}_R1.fq.gz > ./${name}_R1_001.fq
seqkit replace -p "\s.+" ~/${SEQUENCING_RUN}/ZUMIS/multx/${name}_R2.fq.gz > ./${name}_R2_002.fq
done
gzip *
Rscript ~/software/zUMIs/misc/merge_demultiplexed_fastq.R -d . -t 8
fastqc reads_for_zUMIs.R*.fastq.gz
cd ..

## zUMIs
bash ~/software/zUMIs/zUMIs.sh -d ~/software/zUMIs/ -y ~/${SEQUENCING_RUN}/${SAMPLE_ALL}.yaml

## multiqc
rm -r *._STARgenome *._STARpass1
multiqc .

rm ./multx/*.fq.gz ./merge/*.fq.gz
echo "\n*** SMART-zUMIs.sh FINISH! ***\n"
