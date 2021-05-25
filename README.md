# RNA-seq

## From reads to gene expression level
Based on:  
`STAR` (https://github.com/alexdobin/STAR)  
`samtools` (https://www.htslib.org/)  
`rsem` (http://deweylab.github.io/RSEM/)  
`RSeQC` (http://rseqc.sourceforge.net/)  

1. RNA-seq.sh: from .FASTQ to .bam and .gene.results (count/fpkm/tpm)
2. rsem-extract.sh and the three .py file: count fkpm/tpm in different .gene.results file
  
![image](https://github.com/lianggongzifang/RNA-seq/blob/main/RNA-seq.jpg)  
  
## From gene expression level to differential expression gene analysis
Based on:  
`R` (https://www.r-project.org/)  
`DESeq2`(http://www.bioconductor.org/packages/release/bioc/html/DESeq2.html) 

3. DEG analysis by DESeq2.R: from .gene.results file to DEG results  
