## GBS SNP Calling Reference Optional Pipeline (GBS-SNP-CROP v1.1) Tutorial

Authors: Arthur T. O. Melo, Radhika Bartaula, and Iago Hale

Department of Biological Sciences, College of Life Science and Agriculture, University of New Hampshire, Durham, NH, USA.

## Step 1: Parsing the raw reads (GBS-SNP-CROP-1.pl)
```bash
# Parsing paired-end (PE) reads:
perl /path-to-GBS-SNP-CROP/GBS-SNP-CROP-1.pl -d PE -b barcodesID.txt -fq L001 -s 1 -e 2 -enz1 TGCA -enz2 CGG 
# Parsing single-end (SE) reads:
perl /path-to-GBS-SNP-CROP/GBS-SNP-CROP-1.pl -d SE -b barcodesID.txt -fq L001 -s 1 -e 2 -enz1 TGCA -enz2 CGG 
```

## Step 2: Trim based on quality (GBS-SNP-CROP-2.pl)
```bash
# Trimming paired-end (PE) reads:
perl /path-to-GBS-SNP-CROP/GBS-SNP-CROP-2.pl -d PE -fq L001 -t 10 -ph 33 -l 30 -sl 4:30 -tr 30 -m 32
# Trimming single-end (SE) reads:
perl /path-to-GBS-SNP-CROP/GBS-SNP-CROP-2.pl -d SE -fq L001 -t 10 -ph 33 -l 30 -sl 4:30 -tr 30 -m 32
```

## Step 3: Demultiplex (GBS-SNP-CROP-3.pl)
```bash
# Demultiplexing paired-end (PE) reads:
perl /path-to-GBS-SNP-CROP/GBS-SNP-CROP-3.pl -d PE -b barcodesID.txt -fq  L001
# Demultiplexing single-end (SE) reads:
perl /path-to-GBS-SNP-CROP/GBS-SNP-CROP-3.pl -d SE -b barcodesID.txt -fq  L001
```

## Step 4: Cluster reads and assemble the Mock Reference (GBS-SNP-CROP-4.pl)
```bash
# Parsing paired-end (PE) reads:
perl /path-to-GBS-SNP-CROP/GBS-SNP-CROP-4.pl -d PE -b barcodeID.txt -rl 150 -pl 32 -p 0.01 -id 0.93 -t 10 -MR MockRefName
# Parsing single-end (SE) reads:
perl /path-to-GBS-SNP-CROP/GBS-SNP-CROP-4.pl -d SE -b barcodeID.txt -rl 150 -pl 32 -p 0.01 -id 0.93 -t 10 -MR MockRefName
```

## Step 5: Align with BWA-mem and process with SAMTools (GBS-SNP-CROP-5.pl)
```bash
# Mapping paired-end (PE) reads:
perl /path-to-GBS-SNP-CROP/GBS-SNP-CROP-5.pl -d PE-b barcodeID.txt -ref MockRefName.MockRef_Genome.fasta -Q 30 -q 0 -f 2 -F 2308 -t 10 -Opt 0 
# Mapping single-end (SE) reads:
perl /path-to-GBS-SNP-CROP/GBS-SNP-CROP-5.pl -d SE-b barcodeID.txt -ref MockRefName.MockRef_Genome.fasta -Q 30 -q 0 -f 0 -F 2308 -t 10 -Opt 0 
```

## Step 6: Parse mpileup output and produce the SNP discovery master matrix (GBS-SNP-CROP-6.pl)
```bash
perl /path-to-GBS-SNP-CROP/GBS-SNP-CROP-6.pl -b barcodeID.txt -out SNPs_master_matrix.txt
```

## Step 7: Filter SNPs and call genotypes (GBS-SNP-CROP-7.pl)
```bash
perl /path-to-GBS-SNP-CROP/GBS-SNP-CROP-7.pl -in SNPs_master_matrix.txt -out SNPs_genotyping_matrix.txt -mnHoDepth0 11 -mnHoDepth1 48 -mnHetDepth 3 -altStrength 0.9 -mnAlleleRatio 0.1 -mnCall 0.75 -mnAvgDepth 4 -mxAvgDepth 200  
```

## Downstream Tool (GBS-SNP-CROP-8.pl)
```bash
perl /path-to-GBS-SNP-CROP/GBS-SNP-CROP-8.pl -in SNPs_genotyping_matrix.txt -b barcodeID.txt -formats R,Tassel,Plink 
```
