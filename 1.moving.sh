#!/bin/sh

SEQUENCING_RUN="HU000"
I7_INDEX="B41 B42 B43 B44 B45 B46 B47 B48"

#mkdir ${SEQUENCING_RUN}
cd ${SEQUENCING_RUN}
for name in ${I7_INDEX}
do
cp /home/storage/${SEQUENCING_RUN}/I7-${name}_R1.fq.gz ./I7-${name}_1.fq.gz
cp /home/storage/${SEQUENCING_RUN}/I7-${name}_R2.fq.gz ./I7-${name}_2.fq.gz
done

echo "\n*** moving.sh FINISH! ***\n"
