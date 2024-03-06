# mRNA Seqencing Pre-Processing Pipelines
# Ronghao Zhang on 2023-01


# ---------- Script Information ------------------------------
# This script aims to pre-process the Nex-Generation-Sequencing (mRNA) Data.
# Before using this scripts, there are several pre-req need to be done. 
#   1. Download the following packages & tools. It is OK to download the following or newer version. 
#       a. if any package need to be download manually, put them into '~/Desktop/mRNA_rz_2022/seq_tools'
#   2. Download the referencing genome from GenBank or RefSeq of NCBI, rename as 'ref_genome.fna' -> '~/Desktop/mRNA_rz_2022/ref_genome'

# ---------- User Instruction ------------------------------
#   1. Double-check the directories of your work
#   2. Install the Tools listed in the Tools Information Section
#   3. Install bash and vim on user's PC
#   4. use vim to make any changes in the script 
#   5. use bash to run the script

# ---------- Tools Information ------------------------------
# package & tools used
#   1. FastQC High Throughput Sequence QC Report - version 0.11.9
#   2. Trimmomatic - version 0.39
#   3. bwa - version 0.7.17-r1188
#   4. samtools - version 1.13
#   5. bcftools - version 1.13
# referencing genome information: RefSeq NCBI Retreived 2022-07-08.
#   https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Rattus_norvegicus/
#     latest_assembly_versions/GCF_015227675.2_mRatBN7.2/GCF_015227675.2_mRatBN7.2_genomic.fna.gz
#   https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Rattus_norvegicus/
#     latest_assembly_versions/GCF_015227675.2_mRatBN7.2/GCF_015227675.2_mRatBN7.2_genomic.gbff.gz
# Adapter Sequence Information
#   Please Download NEBNext-PE.fa or other equivalent files if you are using different assay. 
#   Please put the .fa file inside your mRNA_data_raw directory together with your raw mRNA data. 


# ---------- Set-Up Directories ------------------------------
# /media/robb-e532/Data/mRNA_rz_2023 -> mRNA_data_raw (21 samples PE in .fastq format)
#                                    -> mRNA_data_processed -> fastqc_untirm
#                                                           -> mRNA_data_trim
#                                                           -> fastqc_trim
#                                                           -> sam
#                                                           -> bam
#                                                           -> docs -> qc_summary
#                                                                   -> bam_flagstat
# /media/robb-e532/Data/mRNA_rz_2022 -> mRNA_data_raw (36 samples PE in .fastq format)
#                                    -> ref_genome 
#                                    -> seq_tools (contains samtools & trimmomatic)
#                                    -> mRNA_data_processed -> fastqc_untirm
#                                                           -> mRNA_data_trim
#                                                           -> fastqc_trim
#                                                           -> sam
#                                                           -> bam
#                                                           -> bcf
#                                                           -> vcf
#                                                           -> docs -> qc_summary
#                                                                   -> bam_flagstat
echo "Setting-Up Directories ..."

cd /media/robb-e532/Data/mRNA_rz_2023/
mkdir mRNA_data_processed
cd /media/robb-e532/Data/mRNA_rz_2023/mRNA_data_processed
mkdir fastqc_untrim mRNA_data_trim fastqc_trim sam bam docs
mkdir /media/robb-e532/Data/mRNA_rz_2023/mRNA_data_processed/docs/qc_summary
mkdir /media/robb-e532/Data/mRNA_rz_2023/mRNA_data_processed/docs/bam_flagstat

echo "Finish Setting-Up Directories Successfully!"


# ---------- FASTQC for Untrimmed Data ------------------------------
echo "Running FASTQC for Untrimmed Data ..."

## perform quality check on all mRNA samples (both forward and reverse) and store in separate file
cd /media/robb-e532/Data/mRNA_rz_2023/mRNA_data_raw 
fastqc -t 16 ./*.fastq.gz 
mv ./*fastqc.html /media/robb-e532/Data/mRNA_rz_2023/mRNA_data_processed/fastqc_untrim
mv ./*fastqc.zip /media/robb-e532/Data/mRNA_rz_2023/mRNA_data_processed/fastqc_untrim 

## summarize qualities and export the fail ones
cd /media/robb-e532/Data/mRNA_rz_2023/mRNA_data_processed/fastqc_untrim
for filename in *.zip
do
  unzip $filename
done
cat */summary.txt > /media/robb-e532/Data/mRNA_rz_2023/mRNA_data_processed/docs/qc_summary/fastqc_summary_untrim.txt
# grep FAIL /media/robb-e532/Data/mRNA_rz_2022/mRNA_data_processed/docs/qc_summary/fastqc_summary_untrim.txt 

echo "FASTQC Summaries Exported!"


# ---------- Filter and Trim mRNA Data ------------------------------
echo "Trimming and Filtering Low-Quality mRNA Data ..."

# use trimmomatic to filter the low quality data
cd /media/robb-e532/Data/mRNA_rz_2023/mRNA_data_raw
for infile in *_1.fastq.gz
do
  base=$(basename ${infile} _1.fastq.gz)
  java -jar /media/robb-e532/Data/mRNA_rz_2022/seq_tools/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 16 \
            ${infile} ${base}_2.fastq.gz \
            ${base}_1_trim.fastq.gz ${base}_1_untrim.fastq.gz \
            ${base}_2_trim.fastq.gz ${base}_2_untrim.fastq.gz \
            SLIDINGWINDOW:4:20 MINLEN:25
done
mv *trim.fastq.gz /media/robb-e532/Data/mRNA_rz_2023/mRNA_data_processed/mRNA_data_trim

echo "Complete Filtering Process!"


# ---------- FASRQC for Trimmed Data ------------------------------
echo "Running FASTQC for Trimmed Data ..."

## perform quality check on trimmed mRNA samples
cd /media/robb-e532/Data/mRNA_rz_2023/mRNA_data_processed/mRNA_data_trim
fastqc -t 16 ./*_trim.fastq.gz 
mv ./*fastqc.html /media/robb-e532/Data/mRNA_rz_2023/mRNA_data_processed/fastqc_trim
mv ./*fastqc.zip /media/robb-e532/Data/mRNA_rz_2023/mRNA_data_processed/fastqc_trim 

## summarize qualities and export the fail ones
cd /media/robb-e532/Data/mRNA_rz_2023/mRNA_data_processed/fastqc_trim
for filename in *.zip
do
  unzip $filename
done
cat */summary.txt > /media/robb-e532/Data/mRNA_rz_2023/mRNA_data_processed/docs/qc_summary/fastqc_summary_trim.txt
# grep FAIL /media/robb-e532/Data/mRNA_rz_2023/mRNA_data_processed/docs/qc_summary/fastqc_summary_trim.txt
          
echo "FASTQC Summaries Exported!"


# ---------- Sequence Alignment w/z ref_genome ------------------------------
echo "Running Sequence Alignment ..."

cd /media/robb-e532/Data/mRNA_rz_2023/mRNA_data_processed/mRNA_data_trim

bwa index /media/robb-e532/Data/mRNA_rz_2022/ref_genome/ref_genome.fna ###index the referencing genome

for infile in *_1_trim.fastq.gz
do
  base=$(basename ${infile} _1_trim.fastq.gz)
  bwa mem -t 16 /media/robb-e532/Data/mRNA_rz_2022/ref_genome/ref_genome.fna ${infile} ${base}_2_trim.fastq.gz > \
          /media/robb-e532/Data/mRNA_rz_2023/mRNA_data_processed/sam/${base}_align.sam
done

echo "Sequence Alignment Completed!"


# ---------- Compress into BAM & Sorting ------------------------------
echo "Compressing SAM and Sorting ..."

cd /media/robb-e532/Data/mRNA_rz_2023/mRNA_data_processed/sam

for infile in *_align.sam
do
  base=$(basename ${infile} _align.sam)
  samtools view --threads 16 -S -b ./${infile} > /media/robb-e532/Data/mRNA_rz_2023/mRNA_data_processed/bam/${base}_align.bam
  samtools sort -@ 8 -m 4G -o /media/robb-e532/Data/mRNA_rz_2023/mRNA_data_processed/bam/${base}_align_sorted.bam \
                   /media/robb-e532/Data/mRNA_rz_2023/mRNA_data_processed/bam/${base}_align.bam 
  samtools flagstat /media/robb-e532/Data/mRNA_rz_2023/mRNA_data_processed/bam/${base}_align.bam > \
                    /media/robb-e532/Data/mRNA_rz_2023/mRNA_data_processed/docs/bam_flagstat/${base}_bam_flagstat.txt
done

echo "BAM Files and Flagstat Exported!"
