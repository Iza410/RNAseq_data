#download reference genome sequence and associated annotation
#checking the downloaded files 
wget https://ftp.ensembl.org/pub/release-108/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
sum Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz #02415 860559



wget https://ftp.ensembl.org/pub/release-108/gtf/homo_sapiens/Homo_sapiens.GRCh38.108.gtf.gz
sum Homo_sapiens.GRCh38.108.gtf.gz #31848 52840

#unpack reference
gzip -d ${ref}/Homo_sapiens.GRCh38.108.gtf.gz
gzip -d ${ref}/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz