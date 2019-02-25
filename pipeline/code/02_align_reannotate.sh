sample=$1
directory=$2

# Number of threads
nt=4

# File path to a bwa genome. Change aligner as so desired
bwa_index="/Volumes/dat/genomes/hg19_bwa/hg19.fa"

bwa mem -t $nt $bwa_index "${directory}/${sample}_1.fastq.gz" "${directory}/${sample}_2.fastq.gz" 2> "${directory}/${sample}.bwa.log" | samtools view -bS - |  samtools sort -@ $nt - -o "${directory}/${sample}.st.bam"

# Reannotate bam files with an annotation rather than having the ID in the name
python code/03_bamReannotate.py --input "${directory}/${sample}.st.bam" --output "${directory}/${sample}.m.bam"
samtools index "${directory}/${sample}.m.bam"

# Remove the original bam file
rm "${directory}/${sample}.st.bam"


