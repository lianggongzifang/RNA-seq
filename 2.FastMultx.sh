#!/bin/sh
# python2!

SEQUENCING_RUN="HU000"
I7_INDEX="B41 B42 B43 B44 B45 B46 B47 B48"

gunzip *
#mkdir FastMultx
cd FastMultx

# for all
for name in ${I7_INDEX}
do
python ~/software/FastMultx.py -i ./${SEQUENCING_RUN}-${name}_index \
-1 ~/${SEQUENCING_RUN}/I7-${name}_1.fq -2 ~/${SEQUENCING_RUN}/I7-${name}_2.fq -t 4 \
--multx_opt "-m 0 -d 2 -x -b"
done
rm -r trim/
cd ..

echo "\n*** FastMultx.sh FINISH! ***\n"
