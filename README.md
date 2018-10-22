# GBS-SNP-CROP

### Introduction
The **GBS** **SNP** **C**alling **R**eference **O**ptional **P**ipeline (GBS-SNP-CROP) is executed via a sequence of [seven Perl scripts][4] that integrate custom parsing and filtering procedures with well-known, vetted bioinformatic tools, giving the user full access to all intermediate files. By employing a novel strategy of variant (SNPs and indels) calling based on the correspondence of within-individual to across-population patterns of polymorphism, the pipeline is able to identify and distinguish high-confidence variants from both sequencing and PCR errors, whether or not a reference genome is available. In the latter case, the pipeline adopts a clustering strategy to build a population-tailored "Mock Reference" using the same GBS data for downstream calling and genotyping. Designed for libraries of either paired-end (PE) or single-end (SE) reads of arbitrary lengths, GBS-SNP-CROP maximizes data usage by eliminating unnecessary data culling due to imposed length uniformity requirements. GBS-SNP-CROP is a complete bioinformatics pipeline developed primarily to support curation, research, and breeding programs wishing to utilize GBS for the cost-effective genome-wide characterization of plant genetic resources.

### Pipeline workflow
**Stage 1. Process the raw GBS data**  
- *Step 1: Parse the raw reads*  
- *Step 2: Trim based on quality and adaptors*  
- *Step 3: Demultiplex*

**Stage 2. Build the Mock Reference**   
- *Step 4: Cluster reads and assemble the Mock Reference*

**Stage 3. Map the processed reads and generate standardized alignment files**  
- *Step 5: Align with BWA-mem and process with SAMtools*  
- *Step 6: Parse mpileup outputs and produce the variants discovery matrix*

**Stage 4. Call Variants and Genotypes**  
- *Step 7: Filter variants and call genotypes*

Below is a schematic of the workflow, with inputs and outputs (boxes) indicated for each step (arrows). 
![](https://github.com/halelab/GBS-SNP-CROP/blob/master/images/workflow.jpg)

### Released versions
[v.4.0][18]: Released on 10/22/2018  
[v.3.0][16]: Released on 2/8/2018  
[v.2.0][14]: Released on 2/22/2017  
[v.1.1][13]: Released on 3/11/2016  
[v.1.0][12]: Released on 1/12/2016

### Getting Help
Initially, go carefully through the [GBS-SNP-CROP User manual][2]. Before post a question or start a discussion, please check your barcode ID file for empty characters or blank spaces and verify that it was saved as a tab-delimited file and also see the [FAQ page][17]. If you're still facing an issue, please, submit it on our discussion [Google groups page][5].

### Requirements
* [Java 7 or higher][6] - We used Java 8
* [Trimmomatic][7] v.0.33 (Bolger et al., 2014)
* [PEAR][8] v.0.96 (Zhang et al., 2014)
* [Vsearch][9] v2.6.2 (Rognes et al., 2016)
* [BWA aligner][10] v.0.7.12 (Li & Durbin, 2009)
* [SAMTools][11] v.1.7 (Li et al., 2009)

### Citing GBS-SNP-CROP
[Melo et al. GBS-SNP-CROP: A reference-optional pipeline for SNP discovery and plant germplasm characterization using genotyping-by-sequencing data. BMC Bioinformatics. 2016. 17:29. DOI 10.1186/s12859-016-0879-y.][1]

[1]:https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-0879-y
[2]:https://github.com/halelab/GBS-SNP-CROP/wiki/GBS-SNP-CROP-User-Manual-(v.4.0)
[3]:http://www.halelab.org
[4]:https://github.com/halelab/GBS-SNP-CROP/tree/master/GBS-SNP-CROP-scripts
[5]:https://groups.google.com/forum/#!forum/gbs-snp-crop
[6]:https://www.java.com/en/
[7]:http://www.usadellab.org/cms/?page=trimmomatic
[8]:https://www.h-its.org/en/research/sco/software/#NextGenerationSequencingSequenceAnalysis
[9]: https://github.com/torognes/vsearch
[10]:http://bio-bwa.sourceforge.net
[11]:http://samtools.sourceforge.net
[12]:https://github.com/halelab/GBS-SNP-CROP/releases/tag/v.1.0
[13]:https://github.com/halelab/GBS-SNP-CROP/releases/tag/v.1.1
[14]:https://github.com/halelab/GBS-SNP-CROP/releases/tag/v.2.0
[15]:http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf
[16]:https://github.com/halelab/GBS-SNP-CROP/releases/tag/v.3.0
[17]:https://github.com/halelab/GBS-SNP-CROP/wiki/Frequently-Asked-Questions-(FAQs)
[18]:https://github.com/halelab/GBS-SNP-CROP/releases/tag/v.4.0
