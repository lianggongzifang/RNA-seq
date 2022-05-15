#!/bin/sh

SEQUENCING_RUN="HU000"

mkdir FastMultx
cd FastMultx
cp ~/${SEQUENCING_RUN}/I7-B41_1.fq.gz ./YF000a_R1.fq.gz
cp ~/${SEQUENCING_RUN}/I7-B41_2.fq.gz ./YF000a_R2.fq.gz
cd ..

echo "\n*** rename.sh FINISH! ***\n"
