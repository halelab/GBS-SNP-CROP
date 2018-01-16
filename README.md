# GBS SNP Calling Reference Optional Pipeline (GBS-SNP-CROP)

### Introduction
The GBS SNP Calling Reference Optional Pipeline (GBS-SNP-CROP) is executed via a sequence of [seven Perl scripts][4] that integrate custom parsing and filtering procedures with well-known, vetted bioinformatic tools, giving the user full access to all intermediate files. By employing a novel strategy of SNP calling based on the correspondence of within-individual to across-population patterns of polymorphism, the pipeline is able to identify and distinguish high-confidence SNPs from both sequencing and PCR errors. The pipeline adopts a clustering strategy to build a population-tailored "Mock Reference" using the same GBS data for downstream SNP calling and genotyping. Designed for libraries of either paired-end (PE) or single-end (SE) reads of arbitrary lengths, GBS-SNP-CROP maximizes data usage by eliminating unnecessary data culling due to imposed length uniformity requirements. GBS-SNP-CROP is a complete bioinformatics pipeline developed primarily to support curation, research, and breeding programs wishing to utilize GBS for the cost-effective genome-wide characterization of plant genetic resources, mainly in the absence of a reference genome. The pipeline, however, can also be used when a reference genome is available, either as a standalone analysis or as a complement to reference-based analyses via alternative pipelines (e.g. TASSEL-GBS) or indeed its own reference-independent analysis.

### Important Notes
**New Version v.2.1** (03/08/2017)

A new version 2.1 will be released within few weeks accomplishing an indel functionality, allowing both SNPs and Indels calls.

**Additional Trimmomatic flag recommended and error fix** (24/02/2017)

GBS-SNP-CROP users, please make note of the following important announcements:

ERROR FIX: We have discovered a minor genotyping error in Script 7 (Line 216), resulting in the incorrect genotyping of secondary allele homozygotes as primary allele homozygotes, specifically in the case where secondary allele read depth is high and primary allele read depth = 1.  This error affects <1% of genotyping calls in our test data. The error has now been corrected, and all users should replace their Script 7 with the version available as of this date (22/02/17). We sincerely apologize for this!

ADDITIONAL TRIMMOMATIC FLAG RECOMMENDED: It has come to our attention that Trimmomatic, by default, can discard high-quality R2 reads if they contain any adapter sequence (e.g. when a GBS fragment length is less than the Illumina read length).  To avoid this unnecessary creation of singletons, and thus data loss to downstream scripts, we recommend that users activate the "keepBothReads" option within the Trimmomatic ILLUMINACLIP parameter, when using GBS-SNP-CROP to analyze paired-end (PE) data.
For example, the current recommended Trimmomatic ILLUMINACLIP parameters for Script 2 are:

```bash
-ad TruSeq3-PE.fa:2:30:10:8:true
```

Please refer to the [Trimmomatic user manual][15] for more details. 

**Version 2.0** (5/11/2016) - Updated on 22/02/2017

The GBS-SNP-CROP v.2.0 was realized. Please access the [realized version][14] for more information.

**Version 1.1** (3/11/2016)

The GBS-SNP-CROP v.1.1 was realized. Please access the [realized version][13] for more information.

**Version 1.0** (1/12/2016)

This is the [original version][12] of the GBS-SNP-CROP pipeline.

### Pipeline workflow
* **Stage 1. Process the raw GBS data**

*Step 1: Parse the raw reads*  
*Step 2: Trim based on quality*   
*Step 3: Demultiplex*

* **Stage 2. Build the Mock Reference** 

*Step 4: Cluster reads and assemble the Mock Reference*

* **Stage 3. Map the processed reads and generate standardized alignment files**

*Step 5: Align with BWA-mem and process with SAMtools*  
*Step 6: Parse mpileup output and produce the SNP discovery master matrix*

* **Stage 4. Call SNPs and Genotypes**

*Step 7: Filter SNPs and call genotypes*

### User Manual
For more details, please see the [GBS-SNP-CROP User Manual][2].

### Discussion Forum
Follow this link to access the [GBS-SNP-CROP Google Group][5].

### Requirements
* [Java 7 or higher][6] - We used Java 8
* [Trimmomatic][7] v.0.33 (Bolger et al., 2014)
* [PEAR][8] v.0.96 (Zhang et al., 2014)
* [Usearch][9] v.8.0.1623 (Edgar, 2010)
* [BWA aligner][10] v.0.7.12 (Li & Durbin, 2009)
* [SAMTools][11] v.1.2 (Li et al., 2009)

### Citing GBS-SNP-CROP
[Melo et al. GBS-SNP-CROP: A reference-optional pipeline for SNP discovery and plant germplasm characterization using genotyping-by-sequencing data. BMC Bioinformatics. 2016. 17:29. DOI 10.1186/s12859-016-0879-y.][1]

[1]:https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-0879-y
[2]:https://github.com/halelab/GBS-SNP-CROP/wiki
[3]:http://www.halelab.org
[4]:https://github.com/halelab/GBS-SNP-CROP/tree/master/GBS-SNP-CROP-scripts
[5]:https://groups.google.com/forum/#!forum/gbs-snp-crop
[6]:https://www.java.com/en/
[7]:http://www.usadellab.org/cms/?page=trimmomatic
[8]:http://sco.h-its.org/exelixis/web/software/pear/
[9]: http://www.drive5.com/usearch/
[10]:http://bio-bwa.sourceforge.net
[11]:http://samtools.sourceforge.net
[12]:https://github.com/halelab/GBS-SNP-CROP/releases/tag/v.1.0
[13]:https://github.com/halelab/GBS-SNP-CROP/releases/tag/v.1.1
[14]:https://github.com/halelab/GBS-SNP-CROP/releases/tag/v.2.0
[15]:http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf
