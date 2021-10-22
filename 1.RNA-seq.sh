#!/bin/sh

SEQUENCING_RUN="HU000"
I7_INDEX="B41 B42 B43 B44 B45 B46 B47 B48"
SAMPLE="YF000"
SAMPLE_ALL="YF000"

for name in ${SAMPLE}
do
mkdir ${name}
cd ${name}

STAR --runThreadN 4 --genomeDir ~/database/STAR/GENCODE_mm10/ \
--readFilesIn ~/rna-seq/${SEQUENCING_RUN}/${name}_1.fq.gz ~/rna-seq/${SEQUENCING_RUN}/${name}_2.fq.gz --readFilesCommand zcat \
--outSAMtype BAM SortedByCoordinate --outSAMattributes All --outFilterMultimapNmax 1 \
--sjdbGTFfile ~/database/genomes/GENCODE/gencode.vM24.annotation.gtf \
--quantMode TranscriptomeSAM GeneCounts --twopassMode Basic --outFileNamePrefix ${name}.
samtools index ${name}.Aligned.sortedByCoord.out.bam
igvtools count -w 1 -e 0 ${name}.Aligned.sortedByCoord.out.bam ${name}.STAR.tdf mm10

rsem-calculate-expression --paired-end --no-bam-output --append-names -p 4 --forward-prob 0.5 \
-bam ${name}.Aligned.toTranscriptome.out.bam ~/database/RSEM/GENCODE_mm10/mouse_gencode ${name}.rsem

geneBody_coverage.py -r ~/database/genomes/GENCODE/mm10_Gencode_VM18.bed -i ${name}.Aligned.sortedByCoord.out.bam -o ${name}
read_distribution.py -r ~/database/genomes/GENCODE/mm10_Gencode_VM18.bed -i ${name}.Aligned.sortedByCoord.out.bam > ${name}.readDistribution.txt
read_duplication.py -i ${name}.Aligned.sortedByCoord.out.bam -o ${name}
junction_saturation.py -r ~/database/genomes/GENCODE/mm10_Gencode_VM18.bed -i ${name}.Aligned.sortedByCoord.out.bam -o ${name}
inner_distance.py -r ~/database/genomes/GENCODE/mm10_Gencode_VM18.bed -i ${name}.Aligned.sortedByCoord.out.bam -o ${name}

rm *._STARgenome
rm *._STARpass1

cd ..
done
