#!/bin/sh
for name in T01 T02 T03 T04 T05 T06 T07 T08 T09 T10 
do
mkdir X930-S01-S01-01-${name}
cd X930-S01-S01-01-${name}
STAR --runThreadN 4 --genomeDir /home2/yuefeng/database/STAR/GENCODE_mm10/ \
--readFilesIn ~/rna-seq/Yang20201221/Mouse_X930-S01-S01-01-${name}_good_1.fq ~/rna-seq/Yang20201221/Mouse_X930-S01-S01-01-${name}_good_2.fq \
--outSAMtype BAM SortedByCoordinate --outSAMattributes All --outFilterMultimapNmax 1 \
--sjdbGTFfile /home2/yuefeng/database/genomes/GENCODE/gencode.vM24.annotation.gtf \
--quantMode TranscriptomeSAM GeneCounts --twopassMode Basic --outFileNamePrefix ${name}.
samtools index ${name}.Aligned.sortedByCoord.out.bam
rsem-calculate-expression --paired-end --no-bam-output --append-names -p 4 --forward-prob 0.5 \
-bam ${name}.Aligned.toTranscriptome.out.bam ~/database/RSEM/GENCODE_mm10/mouse_gencode ${name}.rsem
geneBody_coverage.py -r ~/database/genomes/GENCODE/mm10_Gencode_VM18.bed -i ${name}.Aligned.sortedByCoord.out.bam -o ${name}
read_distribution.py -r ~/database/genomes/GENCODE/mm10_Gencode_VM18.bed -i ${name}.Aligned.sortedByCoord.out.bam > ${name}.readDistribution.txt
read_duplication.py -i ${name}.Aligned.sortedByCoord.out.bam -o ${name}
junction_saturation.py -r ~/database/genomes/GENCODE/mm10_Gencode_VM18.bed -i ${name}.Aligned.sortedByCoord.out.bam -o ${name}
inner_distance.py -r ~/database/genomes/GENCODE/mm10_Gencode_VM18.bed -i ${name}.Aligned.sortedByCoord.out.bam -o ${name}
cd ..
done