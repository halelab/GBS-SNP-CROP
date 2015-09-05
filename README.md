# GBS SNP Calling Reference Optional Pipeline (GBS-SNP-CROP)

This is the official development repository for GBS-SNP-CROP

Please, see [GBS-SNP-CROP User manual][1].

See also Professor Iago Hale [Lab page][2].

###Introduction
The GBS SNP Calling Reference Optional Pipeline (GBS-SNP-CROP) is executed via a sequence of [seven PERL scripts][3] which integrate custom parsing and filtering procedures with well-known, vetted bioinformatic tools, giving the user full access to all intermediate files. By employing a novel strategy of SNP calling based on the correspondence of within-individual to across-population patterns of polymorphism, the pipeline is able to identify and distinguish high-confidence SNPs from both sequencing and PCR errors. The pipeline adopts a clustering strategy to build a population-tailored "Mock Reference" using the same GBS data for downstream SNP calling and genotyping. Designed for libraries of paired-end (PE) read of arbitrary lengths, GBS-SNP-CROP maximizes data usage by eliminating unnecessary data culling due to imposed length uniformity requirements. GBS-SNP-CROP is a complete bioinformatics pipeline developed primarily to support curation, research, and breeding programs wishing to utilize GBS for the cost-effective genome-wide characterization of plant genetic resources, mainly in the absence of a reference genome. The pipeline, however, can also be used when a reference genome is available, either as a standalone analysis or as a complement to reference-based analyses via alternative pipelines (e.g. TASSEL-GBS) or indeed its own reference-independent analysis.

### Pipeline workflow
* **Stage 1. Process the raw GBS data**

*Step 1: Parse the raw reads*

*Step 2: Trim based on quality* 

*Step 3: Demultiplex*

* **Stage 2. Build the Mock Reference** 

*Step 4: Cluster reads and assemble the Mock Reference*

* **Stage 3. Map the processed reads and generate standardized alignment files**

*Step 5: Align with BWA-mem and process with SAMTools*

*Step 6: Parse mpileup output and produce the SNP discovery master matrix*

* **Stage 4. Calling SNPs and Genotypes**

*Step 7: Filter SNPs and call genotypes*

### Citing GBS-SNP-CROP
Melo et al. (2015) GBS-SNP-CROP: A reference-optional pipeline for SNP discovery and plant germplasm characterization using genotyping-by-sequencing data. BMC Bioinformatics. DOI XX.

### Requirements
* [Java 7 or higher][4]
* [Trimmomatic][5] (Bolger et al., 2014)
* [PEAR][6] (Zhang et al., 2014)
* [Usearch][7] (Edgar, 2010)
* [BWA aligner][8] (Li & Durbin, 2009)
* [SAMTools][9](Li et al., 2009)

[1]:https://github.com/halelab/GBS-SNP-CROP/blob/master/GBS-SNP-CROP-User-Manual-v1.0.pdf
[2]:http://www.halelab.org
[3]:https://github.com/halelab/GBS-SNP-CROP/tree/master/GBS-SNP-CROP-scripts
[4]:https://www.java.com/en/
[5]:http://www.usadellab.org/cms/?page=trimmomatic
[6]:http://sco.h-its.org/exelixis/web/software/pear/
[7]: http://www.drive5.com/usearch/
[8]:http://bio-bwa.sourceforge.net
[9]:http://samtools.sourceforge.net


