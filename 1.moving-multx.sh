#!/bin/sh

SEQUENCING_RUN="HU000"
I7_INDEX="B41 B42 B43 B44 B45 B46 B47 B48"

mkdir rawReads multx

cd rawReads
for name in ${I7_INDEX}
do
cp /home/storage/${SEQUENCING_RUN}/I7-${name}_R1.fq.gz ./I7-${name}_1.fq.gz
cp /home/storage/${SEQUENCING_RUN}/I7-${name}_R2.fq.gz ./I7-${name}_2.fq.gz
done
cd ..

cd multx
for name in ${I7_INDEX}
do
fastq-multx -m 0 -d 2 -x -b -B ~/${SEQUENCING_RUN}/${SEQUENCING_RUN}-${name}_index \
~/${SEQUENCING_RUN}/rawReads/I7-${name}_1.fq.gz ~/${SEQUENCING_RUN}/rawReads/I7-${name}_2.fq.gz \
-o %_R1.fq.gz %_R2.fq.gz > ${SEQUENCING_RUN}-${name}_multx.txt
done
cd ..

echo "\n*** moving-multx.sh FINISH! ***\n"
