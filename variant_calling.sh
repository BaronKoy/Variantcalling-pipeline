#This script was initially written by Petr Danecek, at Sanger Institute and later adapted by Baron Koylass for the EBI Bioinformatics Summer School 2022.

sample=samples.txt
output_dir=out
bams=$output_dir/bams
fasta_ref=Homo_sapiens.GRCh38.dna.chromosome.2_15_16.fa.gz
fastqs=fastq
vcfs=vcf

if [ ! -e $fasta_ref.sa ]; then
    echo "Creating bwa index..."
    bwa index $fasta_ref
fi

# Create fasta index
if [ ! -e $fasta_ref.fai ]; then
    echo "Creating fasta index..."
    samtools faidx $fasta_ref
fi

cat $sample | while read sample population; do

    # Do not run bwa for this sample if an indexed BAM already exists. Because
    # the indexing step is run only after the alignment step, it is sufficient
    # to check the presence of the index - when the index exists, also the BAM
    # must exist
    if [ -e $output_dir/bams/$sample.bam.bai ]; then continue; fi

    echo "Mapping $sample..."

    bwa mem $fasta_ref $fastqs/${sample}_1.fq $fastqs/${sample}_2.fq | samtools sort -o $bams/$sample.bam

    # Index the BAM file
    samtools index $bams/$sample.bam
done

for filename in $bams/*.bam; do
	if [ ! -e ${filename%.bam}_filtered.bcf.csi ]; then
		echo $filename
		bcftools mpileup -f $fasta_ref $filename -Ou | bcftools call -mv -Oz -o "${filename%.bam}.bcf"
		bcftools view -i '%QUAL>=20' "${filename%.bam}.bcf" -Ob -o "${filename%.bam}_filtered.bcf"
		bcftools index "${filename%.bam}_filtered.bcf"
		bcftools view "${filename%.bam}_filtered.bcf" > "${filename%.bam}_filtered.vcf"
	fi
done

#Intersect a series of files. 

firstfile=$bams/CHS00479_filtered.vcf
for file in CHS00501 CHS00537 JPT18942 JPT18947 JPT19065; do
#	echo $firstfile
#	echo $file
	bedtools intersect -a $firstfile -b $bams/${file}_filtered.vcf -wa -header > $bams/temp1.vcf
	cat $bams/temp1.vcf > $bams/temp2.vcf
	firstfile=$bams/temp2.vcf
	#exit
done

cat $firstfile > intersection_Asia.vcf

rm $bams/temp1.vcf
rm $bams/temp2.vcf

firstfile=$bams/FIN00272_filtered.vcf
for file in FIN00304 FIN00309 GBR00099 GBR00239 GBR00253; do
#       echo $firstfile
#       echo $file
        bedtools intersect -a $firstfile -b $bams/${file}_filtered.vcf -wa -header> $bams/temp1.vcf
        cat $bams/temp1.vcf > $bams/temp2.vcf
        firstfile=$bams/temp2.vcf
        #exit
done

cat $firstfile > intersection_Europe.vcf

bedtools intersect -a intersection_Asia.vcf -b intersection_Europe.vcf -v -header > Asia_only.vcf
bedtools intersect -b intersection_Asia.vcf -a intersection_Europe.vcf -v -header > Europe_only.vcf
