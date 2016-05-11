# GBS SNP Calling Reference Optional Pipeline (GBS-SNP-CROP)

###Introduction
The GBS SNP Calling Reference Optional Pipeline (GBS-SNP-CROP) is executed via a sequence of [seven Perl scripts][4] that integrate custom parsing and filtering procedures with well-known, vetted bioinformatic tools, giving the user full access to all intermediate files. By employing a novel strategy of SNP calling based on the correspondence of within-individual to across-population patterns of polymorphism, the pipeline is able to identify and distinguish high-confidence SNPs from both sequencing and PCR errors. The pipeline adopts a clustering strategy to build a population-tailored "Mock Reference" using the same GBS data for downstream SNP calling and genotyping. Designed for libraries of either paired-end (PE) or single-end (SE) reads of arbitrary lengths, GBS-SNP-CROP maximizes data usage by eliminating unnecessary data culling due to imposed length uniformity requirements. GBS-SNP-CROP is a complete bioinformatics pipeline developed primarily to support curation, research, and breeding programs wishing to utilize GBS for the cost-effective genome-wide characterization of plant genetic resources, mainly in the absence of a reference genome. The pipeline, however, can also be used when a reference genome is available, either as a standalone analysis or as a complement to reference-based analyses via alternative pipelines (e.g. TASSEL-GBS) or indeed its own reference-independent analysis.

### Updates
**Version 1.1**

This updated version of GBS-SNP-CROP (v1.1) expands the functionality of the original pipeline by accommodating either Paired-End (PE) or Single-End (SE) GBS reads. To access this functionality, a new "data-type" flag (-d) was added so that users can specify their data type ("PE" or "SE") for Steps 1-5. Please refer to the [User Manual][2] for more details.

**Version 2.0**

This updated version of GBS-SNP-CROP (v.2.0) features the following improvements, by script:

*GBS-SNP-CROP-1.pl:*  
* Rather than loading entire FASTQ files into memory for parsing, the raw FASTQ files are now processed read-by-read. This drastically reduces memory requirements and allows the script to parse arbitrarily large datasets on even small computers.  
* Users are now presented a progress bar to let them know the status of this relatively time-intensive parsing step.  

*GBS-SNP-CROP-2.pl:*  
* A new flag **```-ad```** was added, giving users the option to trim Illumina adapters using Trimmomatic's ILLUMINACLIP option.  

*GBS-SNP-CROP-4.pl:*  
* To maximize both data usage and genotype representation in building the mock reference, the script now randomly samples FASTA reads from all designated genotypes for input into USEARCH, up to that program's maximum data limit of 4 Gb (Edgar, 2010). If the total data available for read clustering is less than 4 Gb, the entire dataset is used (i.e. no sampling).   
* When assembling the mock reference, adjacent centroids are now separated by a string of 20 high quality A’s, an adenine-based boundary found to enhance downstream alignment (BWA) performance.  
* A new file called “PosToMask.txt” is created that contains the coordinates of all adenine-based centroid boundaries.  Using this file, all erroneous SNPs called within those boundaries are identified and culled.  

*GBS-SNP-CROP-6.pl:*  
* SNP depth is now more rigorously based on the number of independently sequenced GBS fragments rather than on the number of reads.  For example, an overlapping R1 and R2 read pair from the same GBS fragment increases the depth tally by 1 (i.e. one GBS fragment) rather than 2 (i.e. 1 R1 read + 1 R2 read).  This improvement only affects PE analysis.  
* Monomorphic sites are identified and filtered more efficiently, and users are now presented a progress bar to let them know the status of this relatively time-intensive mpileup parsing step.  

*GBS-SNP-CROP-9.pl - **New Script** (Downstream Tools):*    
* To better support user decision-making in applying subsequent filters based on the distributions of SNPs within centroids/clusters, this new script was developed that extracts and presents the following information for each identified SNP:
(a) the ID of the centroid/cluster containing the SNP; 
(b) the start position (coordinate) of the cluster within the mock reference;
(c) the end position (coordinate) of the cluster within the mock reference;
(d) the length of the cluster;
(e) the SNP position (coordinate) within the cluster; 
(f) the total number of SNPs called within the cluster; 
(g) the distances (bp) between all adjacent SNPs within that same cluster;
(h) the minimum distance (bp) between adjacent SNPs in the same cluster; and
(i) the maximum distance (bp) between adjacent SNPs in the same cluster.
See [User Manual][2] for more details.

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
