# GBS-SNP-CROP v.4.0 Tutorial 
> Arthur T. O. Melo, Radhika Bartaula and Iago Hale   
> Department of Agriculture, Nutrition and Food Systems, University of New Hampshire, Durham, NH, USA

Please, download the tutorial dataset using the link below:

https://www.dropbox.com/sh/oh8fvm3ekpuktkc/AAAKQ3ZEQFA4DcetKtZnkYqya?dl=0

Also, download the barcode ID file to the same dataset directory.

### Step 1: Parsing the raw reads (GBS-SNP-CROP-1.pl)
```bash
# Paired-end (PE) reads:
# All parameters:
perl GBS-SNP-CROP-1.pl -d PE -b barcodesID.txt -fq L001 -s 1 -e 2 -enz1 TGCA -enz2 CGG -t 10
# Required parameters:
perl GBS-SNP-CROP-1.pl -d PE -b barcodesID.txt -fq L001 -s 1 -e 2 -enz1 TGCA -enz2 CGG

# Single-end (SE) reads:
# All parameters:
perl GBS-SNP-CROP-1.pl -d SE -b barcodesID.txt -fq L001 -s 1 -e 2 -enz1 TGCA -enz2 CGG -t 10
# Required parameters:
perl GBS-SNP-CROP-1.pl -d SE -b barcodesID.txt -fq L001 -s 1 -e 2 -enz1 TGCA -enz2 CGG
```

### Step 2: Trim based on quality (GBS-SNP-CROP-2.pl)
```bash
# Paired-end (PE) reads:
# All parameters:
perl GBS-SNP-CROP-2.pl -tm /usr/local/bin/trimmomatic-0.33.jar -d PE -fq L001 -t 10 -ph 33 -ad TruSeq3-PE.fa:2:30:10 -l 30 -sl 4:30 -tr 30 -m 32
# Required parameters:
perl GBS-SNP-CROP/v.4.0/GBS-SNP-CROP-2.pl -d PE -fq L001

# Single-end (SE) reads:
# All parameters:
perl GBS-SNP-CROP-2.pl -tm /usr/local/bin/trimmomatic-0.33.jar -d SE -fq L001 -t 10 -ph 33 -ad TruSeq3-SE.fa:2:30:10 -l 30 -sl 4:30 -tr 30 -m 32
# Required parameters:
perl GBS-SNP-CROP-2.pl -d SE -fq L001 -ad TruSeq3-SE.fa:2:30:10
```

### Step 3: Demultiplex (GBS-SNP-CROP-3.pl)
```bash
# Paired-end (PE) reads:
# All and required parameters:
perl GBS-SNP-CROP-3.pl -d PE -b barcodesID.txt -fq L001

# Single-end (SE) reads:
# All and required parameters:
perl GBS-SNP-CROP-3.pl -d SE -b barcodesID.txt -fq L001
```

### Step 4: Cluster reads and assemble the Mock Reference (GBS-SNP-CROP-4.pl)
```bash
# Paired-end (PE) reads:
# All parameters:
perl GBS-SNP-CROP-4.pl -pr /usr/local/bin/pear -vs /usr/local/bin/vsearch -d PE -b barcodesID.txt -t 10 -cl consout -rl 150 -pl 32 -p 0.01 -id 0.93 -min 32 -MR MR
# Required parameters:
perl GBS-SNP-CROP-4.pl -d PE -b barcodesID.txt

# Single-end (SE) reads:
# All parameters:
perl GBS-SNP-CROP-4.pl -pr /usr/local/bin/pear -vs /usr/local/bin/vsearch -d SE -b barcodesID.txt -t 10 -cl consout -rl 150 -pl 32 -p 0.01 -id 0.93 -min 32 -MR MR
# Required parameters:
perl GBS-SNP-CROP-4.pl -d SE -b barcodesID.txt
```

### Step 5: Align with BWA-mem and process with SAMTools (GBS-SNP-CROP-5.pl)
```bash
# Paired-end (PE) reads:
# All parameters:
perl GBS-SNP-CROP-5.pl -bw /usr/local/bin/bwa -st /usr/local/bin/samtools -d PE -b barcodesID.txt -ref MR.Genome.fa -Q 30 -q 30 -F 2308 -f 2 -t 10 -Opt 0
# Required parameters:
perl GBS-SNP-CROP-5.pl -d PE -b barcodesID.txt

# Single-end (SE) reads:
# All parameters:
perl GBS-SNP-CROP-5.pl -bw /usr/local/bin/bwa -st /usr/local/bin/samtools -d SE -b barcodesID.txt -ref MR.Genome.fa -Q 30 -q 30 -F 2308 -f 0 -t 10 -Opt 0
# Required parameters:
perl GBS-SNP-CROP-5.pl -d SE -b barcodesID.txt -f 0
```

### Step 6: Parse mpileup output and produce the variants discovery master matrix (GBS-SNP-CROP-6.pl)
```bash
# Discovery only SNPs:
All parameters:
perl GBS-SNP-CROP-6.pl -b barcodesID.txt -out MasterMatrix.txt -p snp -t 10
# Required parameters:
perl GBS-SNP-CROP-6.pl -b barcodesID.txt -p snp

# Discovery SNPs and indels:
All parameters:
perl GBS-SNP-CROP-6.pl -b barcodesID.txt -out MasterMatrix.txt -p indel -t 10
# Required parameters:
perl GBS-SNP-CROP-6.pl -b barcodesID.txt -p indel
```

### Step 7: Filter the variants and call genotypes (GBS-SNP-CROP-7.pl)
```bash
# Call only SNPs:
# All parameters:
perl GBS-SNP-CROP-7.pl -in MasterMatrix.txt -out GenoMatrix.txt -p snp -mnHoDepth0 5 -mnHoDepth1 20 -mnHetDepth 3 -altStrength 0.8 -mnAlleleRatio 0.25 -mnCall 0.75 -mnAvgDepth 3 -mxAvgDepth 200
# Required parameters:
perl GBS-SNP-CROP-7.pl -p snp 

# Call both SNPs and indels:
# All parameters:
perl GBS-SNP-CROP-7.pl -in MasterMatrix.txt -out GenoMatrix.txt -p indel -mnHoDepth0 5 -mnHoDepth1 20 -mnHetDepth 3 -altStrength 0.8 -mnAlleleRatio 0.25 -mnCall 0.75 -mnAvgDepth 3 -mxAvgDepth 200
# Required parameters:
perl GBS-SNP-CROP-7.pl -p indel
```

### Downstream Tool 1: Creating input files to software packages R, TASSEL GUI and/or PLINK and create VCF output (GBS-SNP-CROP-8.pl)
```bash
perl GBS-SNP-CROP-8.pl -in GenoMatrix.txt -out output -b barcodeID.txt -formats R,Tassel,Plink,VCF,HetFreq 
```

### Downstream Tool 2: Provide the cluster/centroid ID and other descriptors for all called SNPs (GBS-SNP-CROP-9.pl)
```bash
perl GBS-SNP-CROP-9.pl -in GenotMatrix.txt -out output -ref MR.Clusters.fa 
```
