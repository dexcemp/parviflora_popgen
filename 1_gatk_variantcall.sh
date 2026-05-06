# Remember to make the dictionary and indices if needed from the reference
#gatk CreateSequenceDictionary -R H230-3.ref.fasta
#bwa index H230-3.ref.fasta
#samtools faidx H230-3.ref.fasta

# This is set up like a text sample list (see the bottom) because I wanted to match filtering steps from another analysis
# So I used an ASTRAL analysis for the sample list
# This could easily be rewritten to glob a list
while read p; do
echo $p

# Set variables
read1fn=`find /mnt/Botbot/sequencing_reads_processed/heuchera*/trimmed/ -name "${p}_P1.fastq" | tail -n 1`
# Resolve duplicates by taking last element
# Read folders are alphabetized so that they glob in chronological order
echo $read1fn

read2fn=`find /mnt/Botbot/sequencing_reads_processed/heuchera*/trimmed/ -name "${p}_P2.fastq" | tail -n 1`
echo $read2fn

#Align read files to reference sequence and map
bwa mem H230-3.ref.fasta $read1fn $read2fn | samtools view -bS - | samtools sort - -o "$p.sorted.bam"
gatk FastqToSam -F1 $read1fn -F2 $read2fn -O $p.unmapped.bam -SM $p.sorted.bam

#Replace read groups to mapped and unmapped bam files using library prep and sequencing information
gatk AddOrReplaceReadGroups -I  $p.sorted.bam -O $p.sorted-RG.bam -RGID 2 -RGLB lib1 -RGPL illumina -RGPU unit1 -RGSM $p
gatk AddOrReplaceReadGroups -I  $p.unmapped.bam -O $p.unmapped-RG.bam -RGID 2 -RGLB lib1 -RGPL illumina -RGPU unit1 -RGSM $p

#Combine mapped and unmapped BAM files
gatk MergeBamAlignment --ALIGNED_BAM $p.sorted-RG.bam --UNMAPPED_BAM $p.unmapped-RG.bam -O $p.merged.bam -R H230-3.ref.fasta

#Remove duplicate sequences
gatk MarkDuplicates -I $p.merged.bam -O $p.marked.bam -M $p.metrics.txt
samtools index $p.marked.bam

#Create GVCF 
gatk HaplotypeCaller -I $p.marked.bam -O $p-g.vcf -ERC GVCF -R H230-3.ref.fasta


done <sample_list.txt