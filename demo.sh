Gffinfor GCF_000001405.38_GRCh38.p12_genomic.gff.gz


Gffinfor GCF_000001405.38_GRCh38.p12_genomic.gff.gz gene,mRNA


Gffinfor GCF_000001405.38_GRCh38.p12_genomic.gff.gz \
gene,pseudogene gene,gene_biotype,description,GeneID,HGNC \
gene:gene_name,description:gene_description|
awk '++x[$1]==1' > gene.infor.tsv

# pigz -dc GCF_000001405.38_GRCh38.p12_genomic.gff.gz |
# gffread -EFG - -T -o genomic.gtf

Gffinfor genomic.gtf exon \
transcript_id,gbkey,Name,product,gene_id,gene_biotype,gene_name,description  \
gbkey:transcript_biotype,Name:transcript_name,description:gene_description |
awk 'BEGIN{FS=OFS="\t"} ++x[$1]==1{print}' > transcript.infor.tsv


Gffinfor genomic.gtf exon 0,gbkey | awk 'BEGIN{FS=OFS="\t"}
$10=="mRNA" && NR>1{sub(" gbkey \"", " transcript_biotype \""); 
sub(" Name \"", " transcript_name \"");
sub(" description \"", " gene_description \""); print}' |
cut -f1-9 > mRNA.gtf
