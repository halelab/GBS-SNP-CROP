# GBS-SNP-CROP v.4.1 Tutorial 
> Arthur T. O. Melo, Radhika Bartaula and Iago Hale   
> Department of Agriculture, Nutrition and Food Systems, University of New Hampshire, Durham, NH, USA

You can access the tutorial dataset using the following link:

https://www.dropbox.com/sh/oh8fvm3ekpuktkc/AAAKQ3ZEQFA4DcetKtZnkYqya?dl=0

Download the provided barcode ID file (barcodesID.txt) to the same directory and proceed with the following commands.

### Step 1: Parse the raw reads (GBS-SNP-CROP-1.pl)
```bash
# Paired-end (PE) reads:
perl GBS-SNP-CROP-1.pl -d PE -b barcodesID.txt -fq L001 -s 1 -e 2 -enz1 TGCA -enz2 CGG -t 10

# Single-end (SE) reads:
perl GBS-SNP-CROP-1.pl -d SE -b barcodesID.txt -fq L001 -s 1 -e 2 -enz1 TGCA -enz2 CGG -t 10
```

### Step 2: Trim based on quality and adaptors (GBS-SNP-CROP-2.pl)
```bash
# Paired-end (PE) reads:
perl GBS-SNP-CROP-2.pl -tm /usr/local/bin/trimmomatic-0.39.jar -d PE -fq L001 -t 10 -ph 33 -ad TruSeq3-PE.fa:2:30:10 -l 30 -sl 4:30 -tr 30 -m 32

# Single-end (SE) reads:
perl GBS-SNP-CROP-2.pl -tm /usr/local/bin/trimmomatic-0.39.jar -d SE -fq L001 -t 10 -ph 33 -ad TruSeq3-SE.fa:2:30:10 -l 30 -sl 4:30 -tr 30 -m 32
```

### Step 3: Demultiplex (GBS-SNP-CROP-3.pl)
```bash
# Paired-end (PE) reads:
perl GBS-SNP-CROP-3.pl -d PE -b barcodesID.txt -fq L001

# Single-end (SE) reads:
perl GBS-SNP-CROP-3.pl -d SE -b barcodesID.txt -fq L001
```

### Step 4: Cluster reads and assemble the Mock Reference (GBS-SNP-CROP-4.pl)
```bash
# Paired-end (PE) reads:
perl GBS-SNP-CROP-4.pl -pr /usr/local/bin/pear -vs /usr/local/bin/vsearch -d PE -b barcodesID.txt -rl 150 -p 0.01 -pl 32 -t 10 -cl consout -id 0.93 -min 32 -MR GSC.MR

# Single-end (SE) reads:
perl GBS-SNP-CROP-4.pl -vs /usr/local/bin/vsearch -d SE -b barcodesID.txt -rl 150 -t 10 -cl consout -id 0.93 -min 32 -MR GSC.MR
```

### Step 5: Align with BWA-mem and process with SAMtools (GBS-SNP-CROP-5.pl)
```bash
# Paired-end (PE) reads:
perl GBS-SNP-CROP-5.pl -bw /usr/local/bin/bwa -st /usr/local/bin/samtools -d PE -b barcodesID.txt -ref GSC.MR.Genome.fa -Q 30 -q 30 -F 2308 -f 2 -t 10 -opt 0

# Single-end (SE) reads:
perl GBS-SNP-CROP-5.pl -bw /usr/local/bin/bwa -st /usr/local/bin/samtools -d SE -b barcodesID.txt -ref GSC.MR.Genome.fa -Q 30 -q 30 -F 2308 -f 0 -t 10 -opt 0
```

### Step 6: Parse mpileup output and produce the variant discovery master matrix (GBS-SNP-CROP-6.pl)
```bash
# Identifying SNPs only:
perl GBS-SNP-CROP-6.pl -b barcodesID.txt -out GSC.MasterMatrix.txt -t 10

# Identifying both SNPs and indels:
perl GBS-SNP-CROP-6.pl -b barcodesID.txt -out GSC.MasterMatrix.txt -indel -t 10
```

### Step 7: Filter variants and call genotypes (GBS-SNP-CROP-7.pl)
```bash
# Calling SNPs only:
perl GBS-SNP-CROP-7.pl -in GSC.MasterMatrix.txt -out GSC.GenoMatrix.txt -mnHoDepth0 5 -mnHoDepth1 20 -mnHetDepth 3 -altStrength 0.8 -mnAlleleRatio 0.25 -mnCall 0.75 -mnAvgDepth 3 -mxAvgDepth 200

# Calling SNPs and indels:
perl GBS-SNP-CROP-7.pl -in GSC.MasterMatrix.txt -out GSC.GenoMatrix.txt -indel -mnHoDepth0 5 -mnHoDepth1 20 -mnHetDepth 3 -altStrength 0.8 -mnAlleleRatio 0.25 -mnCall 0.75 -mnAvgDepth 3 -mxAvgDepth 200
```

### Downstream Tool 1: Create input files for downstream analyses (see User Manual) (GBS-SNP-CROP-8.pl)
```bash
perl GBS-SNP-CROP-8.pl -in GSC.GenoMatrix.txt -out GSC -b barcodesID.txt -formats R,T,P,V,H 
```

### Downstream Tool 2: Provide descriptors for all called variants (GBS-SNP-CROP-9.pl)
```bash
perl GBS-SNP-CROP-9.pl -in GSC.GenoMatrix.txt -out GSC -ref GSC.MR.Clusters.fa 
```
