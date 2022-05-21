#!/bin/sh

SEQUENCING_RUN="HU000"
I7_INDEX="B41 B42 B43 B44 B45 B46 B47 B48"
SAMPLE="YF000a"
SAMPLE_ALL="YF000"
PREFIX="SMART"

mkdir -p ${PREFIX}/cut-index ${PREFIX}/STAR ${PREFIX}/rsem ${PREFIX}/rseqc
cd ${PREFIX}

## cut ME and index
cd cut-index
for name in ${SAMPLE}
do
#seqkit -j 8 subseq -r 26:-1 -o ${name}_R1_MEcut.fq.gz ~/${SEQUENCING_RUN}/multx/${name}_R1.fq.gz
seqkit -j 8 subseq -r 20:-1 -o ${name}_R1_MEcut.fq.gz ~/${SEQUENCING_RUN}/multx/${name}_R1.fq.gz
seqkit -j 8 subseq -r 20:-1 -o ${name}_R2_MEcut.fq.gz ~/${SEQUENCING_RUN}/multx/${name}_R2.fq.gz
done
fastqc *.fq.gz
cd ..

## STAR
cd STAR
for name in ${SAMPLE}
do
STAR --runThreadN 8 --genomeDir ~/database/STAR/GENCODE_mm10/ \
--readFilesIn ~/${SEQUENCING_RUN}/${PREFIX}/cut-index/${name}_R1_MEcut.fq.gz ~/${SEQUENCING_RUN}/${PREFIX}/cut-index/${name}_R2_MEcut.fq.gz \
--readFilesCommand zcat \
--sjdbGTFfile ~/database/genomes/GENCODE/gencode.vM24.annotation.gtf \
--quantMode TranscriptomeSAM GeneCounts --twopassMode Basic \
--outSAMtype BAM SortedByCoordinate --outSAMattributes All --outFilterMultimapNmax 1 \
--outFileNamePrefix ${name}.
samtools index ${name}.Aligned.sortedByCoord.out.bam
igvtools count -w 1 -e 0 ${name}.Aligned.sortedByCoord.out.bam ${name}.STAR.tdf mm10
done
rm -r *._STARgenome *._STARpass1
cd ..

## rsem count
cd rsem
for name in ${SAMPLE}
do
rsem-calculate-expression --paired-end --no-bam-output --append-names -p 8 --forward-prob 0.5 \
-bam ~/${SEQUENCING_RUN}/${PREFIX}/STAR/${name}.Aligned.toTranscriptome.out.bam \
~/database/RSEM/GENCODE_mm10/mouse_gencode ${name}.rsem
done
cd ..

## rseqc
cd rseqc
for name in ${SAMPLE}
do
read_duplication.py -i ~/${SEQUENCING_RUN}/${PREFIX}/STAR/${name}.Aligned.sortedByCoord.out.bam -o ${name}
read_distribution.py -r ~/database/genomes/GENCODE/mm10_Gencode_VM18.bed \
-i ~/${SEQUENCING_RUN}/${PREFIX}/STAR/${name}.Aligned.sortedByCoord.out.bam > ${name}.readDistribution.txt
geneBody_coverage.py -r ~/database/genomes/GENCODE/mm10_Gencode_VM18.bed \
-i ~/${SEQUENCING_RUN}/${PREFIX}/STAR/${name}.Aligned.sortedByCoord.out.bam -o ${name}
junction_saturation.py -r ~/database/genomes/GENCODE/mm10_Gencode_VM18.bed \
-i ~/${SEQUENCING_RUN}/${PREFIX}/STAR/${name}.Aligned.sortedByCoord.out.bam -o ${name}
inner_distance.py -r ~/database/genomes/GENCODE/mm10_Gencode_VM18.bed \
-i ~/${SEQUENCING_RUN}/${PREFIX}/STAR/${name}.Aligned.sortedByCoord.out.bam -o ${name}
done
cd ..

multiqc .

list_csv=""
list_name="gene"
for i in  ./rsem/*genes.results
do
	echo $i
	list_csv=${list_csv}" "${i}
	i_name1=${i%/*}
	i_name2=${i_name1##*/}
	list_name=${list_name}"       "${i_name2}
	echo ${list_name}
done
echo $list_name > ./gene-count-matrix.txt
cp ~/software/rsem-count-extract.py .
python ./rsem-count-extract.py ${list_csv} >> ./${SAMPLE}.gene-count-matrix.txt

list_csv=""
list_name="gene"
for i in  ./rsem/*genes.results
do
	echo $i
	list_csv=${list_csv}" "${i}
	i_name1=${i%/*}
	i_name2=${i_name1##*/}
	list_name=${list_name}"       "${i_name2}
	echo ${list_name}
done
echo $list_name > ./gene-tpm-matrix.txt
cp ~/software/rsem-tpm-extract.py .
python ./rsem-tpm-extract.py ${list_csv} >> ./${SAMPLE}.gene-tpm-matrix.txt

list_csv=""
list_name="gene"
for i in  ./rsem/*genes.results
do
	echo $i
	list_csv=${list_csv}" "${i}
	i_name1=${i%/*}
	i_name2=${i_name1##*/}
	list_name=${list_name}"       "${i_name2}
	echo ${list_name}
done
echo $list_name > ./gene-fpkm-matrix.txt
cp ~/software/rsem-fpkm-extract.py .
python ./rsem-fpkm-extract.py ${list_csv} >> ./${SAMPLE}.gene-fpkm-matrix.txt

echo "\n*** SMART-STAR-rsem-count.sh FINISH! ***\n"
