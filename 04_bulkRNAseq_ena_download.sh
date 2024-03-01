# Download Proximal Tubule Sequencing Data from ENA Database
cd /media/robb-e532/Data/mRNA_rz_2023/mRNA_data_raw
while read f; do
	echo "start downloading $f ..."
	~/Desktop/sratoolkit.3.0.2-ubuntu64/bin/fasterq-dump --split-files $f
	echo "start zipping $f ..."
	gzip *.fastq
	rm *.fastq
	echo "The download of $f completed successfully!"
done < filereport_PRJNA244440_tsv.txt
