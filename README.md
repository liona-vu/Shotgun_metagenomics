Workflow detailing the steps for shotgun sequencing comparing between vegans diet and omnivore diets.

## Introduction
The human gut microbiome has been suggested to play an important role in human health such as producing vital nutrients and vitamins and aiding with breaking down complex carbohydrates and fibres. Several bacterial phyla that include *Firmicutes, Bacteroidetes, Actinobacteria, Proteobacteria, Fusobacteria*, and *Verrucomicrobia* generally make up the gut microbiota (Hou et al., 2022). Notably, the gut microbiome has shown to be highly influenced by the type of diet consumed. In modern times, studies have suggested diets rich in fruits and vegetables (or a vegan diet) have different bacterial compositions as compared to diets rich in meat proteins or omnivore diets (Fackelmann et al., 2025). Specifically, vegan diets tend to be richer in fibre-based foods which have resulted in more abundance of bacterial taxa such as Bacteroidetes and Prevotella (Tomova et al., 2019). Recent advances in shotgun metagenomics sequencing have revolutionized the field of microbial studies and can help us gain insight into the different types of bacterial species or the microbial diversity present in the human gut microbiome. In fact, shotgun metagenomics has the advantage of characterizing bacteria that have been historically difficult to cultivate in laboratory. Therefore, the overall goal of shotgun metagenomics is to provide a comprehensive overview of the types of microbials organisms and composition that are present in a certain environment or sample. From this, we can infer their functional roles and their interactions with the host environment. Specifically for this analysis, the objective is to perform a comparison of the gut microbiome diversity, and the different types of bacteria present between vegan and omnivore diets. The current hypothesis that vegan diets harbouring more bacteria specialized in digesting fibre and plant matter and less bacteria that consume dairy or proteins compared to omnivores will be tested.

There are many bioinformatics tools specialized for analyzing shotgun metagenomics data. Some examples of metagenomic profilers include Kraken2 (Wood et al., 2019), CLAssifier based on Reduced K-mers (CLARK) (Ounit et al., 2015), Kaiju (Menzel et al., 2016), and MetaPhlan 4 (Blanco-Míguez et al., 2023). Each tool uses a different classification approach in which both Kraken2 and CLARK use a DNA-based classification that examines and matches exact k-mers, Kaiju is a protein-based alignment, and MetaPhlAn 4 is marker gene-based aligner. Recent benchmarking studies have compared the efficacy and accuracy of these methods. Overall, Kraken2 was determined to be the superior classifier as it has been reported to have the highest F1 score in terms of performance at 0.74 while MetaPhlAn 4 and Kaiju had F1 scores of 0.41 and 0.48, respectively (Edwin et al., 2024). In addition, protein-based classifiers have the highest misclassification rates, ranging from 5-15% compared to the k-mer-based classifiers at 1-5% (Ye et al., 2019). In addition, Kracken2 when combined with Bayesian Reestimation of Abundance with KrakEN (Bracken) (Lu et al., 2017) has a more accurate abundance profile than CLARK (Ye et al., 2019). A disadvantage of Kraken2 is its longer computation time due to its higher requirement of memory compared to Kaiju which uses less than 2GB of memory. A major disadvantage of Kraken2 is its false positive rates; with as high as 5%. KrakenUniq (Breitwieser et al., 2018) has been proposed as an alternative to Kraken2 to mitigate the false positive problem however, it is very memory-intensive and requires at least 200GB of RAM even with the usage of smaller databases such as the Standard database. To mitigate the RAM limitations, bioinformaticians have increased the confidence level to 0.15 instead of the default 0 for Kracken2. Notably, the confidence score is one of the special features of Kraken2. In fact, it has been reported that a confidence level of 0.15 is a reasonable threshold for classification (Nyström-Persson et al., 2025). Any higher confidence levels can also be used and leads to more accuracy but at the cost of reduced sensitivity and increased numbers of unclassified reads (Liu et al., 2024).

One of the biggest decisions when performing shotgun metagenomics is database selection. The completeness and accuracy of the database is critical to ensure accurate and comprehensive taxonomic classification of microbial metagenomic sequencing. Examples of databases include, but not limited to, Minikraken, the subsampled 8GB and 16GB Standard databases, the full Standard database, and the core_nt database. Benchmarking studies employing Kraken2 with simulated metagenomic datasets have compared the taxonomic performance of each database. In fact, the use of smaller databases such as Minikraken and the Standard-16GB database led to approximately 0% of classified reads when confidence scores of 0.40 were utilized (Liu et al., 2024). In contrast, the larger databases such as the core nt database and the GTDB r202 database had classification reads of 93% and 81%, respectively, and at the confidence score of 0.4 (Liu et al., 2024). Therefore, if space and computational resources permit, larger databases, especially the core nt database, are usually better performers for taxonomic classification compared to smaller or subsampled databases due to their coverage of all known current sequences.

Statistical software is also important to consider when performing shotgun metagenomics to determine if there are differentially abundant microbial species between different groups. Many programs exist such as edgeR (Robinson et al., 2010), the newer ANOVA-Like Differential Expression 3 (Aldex3) (Nixon et al., 2025), which is the successor to Aldex2, and the Analysis of Compositions of Microbiomes with Bias Correction 2 (ANCOM-BC2) (Lin & Peddada, 2020, 2024). These methods differ in their statistical methods in that edgeR uses a negative bionomical distribution model (Robinson et al., 2010), Aldex3 uses the Dirichlet–multinomial approach (Nixon et al., 2025) and ANCOM-BC2 uses a log-linear regression to perform differential abundance (Lin & Peddada, 2020, 2024). Benchmarking studies compared different types of differential abundance softwares and reported that Aldex2 (note that no studies have compared Aldex3’s performance yet) and ANCOM-BC2 are the most robust and agree the best when reporting results compared to other software such as edgeR where they reported more false positive discovery rates (FDR) (Nearing et al., 2022). In fact, edgeR was reported not to handle compositional data as well as ANCOM-BC2 which has the advantage of bias correction. In fact, ANCOM-BC2 was reported to be one of the most conservative tools as it identifies FDR of 0-5% compared to edgeR’s FDR ranging from 10 - 30% (Nearing et al., 2022).

## Metholodogy
Based on previously evaluated tools, Kraken2/Bracken, and ANCOM-BC2 will be used to perform analysis on shotgun metagenomics. 

### Obtaining dataset and taxonomic classification
Paired end fastq files of SRR8146968, SRR8146970- SRR8146978 samples were downloaded from the National Centre for Biotechnology Information (NCBI) (https://www.ncbi.nlm.nih.gov/sra/?term=SRP126540) with 5 samples each downloaded for 2 diet groups: vegan and omnivores for a total of 10 samples. Samples were taken from residents living in Turin, Italy. The data was previously generated using the Illumina HiSeq platform (De Filippis et al., 2019). All analysis and modules were loaded on the Digital Research Alliance of Canada with the usage of the Nibi high-performance computing cluster unless otherwise stated. The sra-toolkit (v3.0.9) (Leinonen et al., 2011) was used to download the data with the “prefetch” command employed to download the SRA files. Consequently, the “fasterq-dump” command was employed with the “--split-files” flag to convert the SRA files to fastq files, producing paired end reads (i.e. SRR8146972_1.fastq  SRR8146972_2.fastq). FastQC (v0.12.1) (Andrews, 2010) was used to assess for quality of the data files to determine adapter contamination. fastp (v1.0.1) (Chen et al., 2018) was employed to remove residual sequencing adapters. Subsequently, taxonomic classifications were  performed using  kraken2 (v2.1.6) (Wood et al., 2019) with the following flags “--confidence 0.15”, “--threads 32”, “--report” to generate a kraken report, and  “--paired data” with the core_nt database reference database (https://genome-idx.s3.amazonaws.com/kraken/k2_core_nt_20251015.tar.gz).  Bracken (v3.0) (Lu et al., 2017) was employed to estimate species abundance with the following flags employed: “-r 150” due to the 150bp short-read length, “-l S” for species classification, “-t 0”, and input from the kraken reports. kraken-biom (v 1.2.0) (Dabdoub, 2016) was downloaded on bioconda due to the module not being available on the Nibi cluster and was employed to convert the all species_report files to one combined biom file with the the flag “--fmt json" to convert to json format. 

### Alpha and beta diversity, and differentially abundant taxa
Biom files were read into RStudio (v4.5.1 2025-06-13) using the biomformat (v1.38.0) (P. McMurdie & Paulson, 2025) package. Host DNA was not removed prior to classification but taxa that were classified as non-microbial organisms (such as chordata) or viruses were filtered out prior to downstream analysis using the microviz (v0.13.0) package (Barnett et al., 2021) via the tax_select() function as bacterial reads can be misidentified as eukaryotic taxa reads in shotgun metagenomics (Lind & Pollard, 2021). The phyloseq package (v 1.54.2) (P. J. McMurdie & Holmes, 2013) was used to calculate relative abundance at the phylum level. Alpha diversity metrics included observed, chao1, Shannon’s index and Simpson’s index and were plotted as a scatter plot. For beta diversity, phyloseq was also used to perform ordination via the Principal Coordinate Analysis (PCoA) method using both Bray-Curtis’s and Jaccard’s distance where binary was set to TRUE for the latter. The Permutational Multivariate Analysis of Variance Using Distance Matrices (PERMANOVA) was calculated using the adonis2() function from the vegan (v 2.7-3) package (Oksanen et al., 2026). Statistical modeling was performed using ANCOM-BC2 (v2.12.0) package (Lin & Peddada, 2024) and was used to calculate differential abundant taxa with the formula ~diet_type with the Benjamini-Hochberg (BH) used for adjusted p values to control for false positives and struc_zero parameter set to TRUE. Plots were generated using the ggplot2 (v 4.0.2) package (Wickham, 2016).

## Results

|SRR Number | Sample Name | Spot Reads | Diet Type |
|-----------|--------|------------|-----------|
|SRR8146968|VOV58_metag| 42,163,237|Vegan|
|SRR8146970| VOV36_metag| 31,282,154| Omnivore|
|SRR8146971	|VOV28_metag|38,662,403|	Omnivore|
|SRR8146972 |VOV26_metag|27,306,095|Omnivore|
|SRR8146973	|VOV114_metag	|34,582,517	|Vegan|
|SRR8146974	|VOV77_metag	|35,418,839	|Vegan|
|SRR8146975	|VOV77_metag	|35,602,941	|Omnivore|
|SRR8146976	|VOV70_metag	|28,577,678 |Omnivore|
|SRR8146977	|VOV29_metag	|39,237,195	|Vegan|
|SRR8146978 |VOV20_metag	|41,028,738	|Vegan|

**Table 1. Metadata of samples used for this analysis including the SRR number, sample name, spot reads, and diet type.**

<br>

|SRR Number | Classified (%) | Unclassified (%)|
|-----------|--------|------------|
|SRR8146972 |79.86	|20.14|
|SRR8146970 |57.64	|42.36|
|SRR8146971	|59.84	|40.16|
|SRR8146972	|45.42	|54.58|
|SRR8146973	|45.42	|54.58|
|SRR8146974	|63.57	|36.43|
|SRR8146975	|55.27	|44.73|
|SRR8146976	|59.13	|40.87|
|SRR8146977	|39.42	|60.58|
|SRR8146978 |69.63	|30.37|

**Table 2. Percent of classified and unclassified reads for each SRR sample number.** Briefly, Kraken2 was performed and outputted a report with all classification of taxa and its corresponding classified and unclassified percentages.

<br>

<img width="1278" height="809" alt="20260321_relative_abundance_top_15_barplot" src="https://github.com/user-attachments/assets/427905b1-6473-41a6-935c-733bbc6f991d" />

**Figure 1. Relative abundance of gut microbiota bacterial at the phylum level.** Top 15 phyla from n=10 samples were plotted as determined by the number of raw counts while the other phyla were classified as Other. Samples are labeled on the x-axis as their corresponding SRR number. *Bacillota* and *Bacteroidata* make up most of the relative abundance across all 10 samples.

<br>

<img width="2291" height="1943" alt="20250317_alpha_diversity_labeled" src="https://github.com/user-attachments/assets/39cbbc7f-6233-46b3-bfe2-cf13b26791ec" />
<br>

**Figure 2. Alpha diversity index plots for all SRR samples.** Measures include **A.** observed, **B.** Chao1, **C.** Shannon and **D.** Simpson. Diet types for omnivores (salmon) and vegans (turquoise) are shown. No discernable patterns can be observed between the two diet types. 

|SRR Number	|Observed| Chao1   | Shannon | Simpson |
|-----------|--------|---------|---------|---------|
|SRR8146968	|829	| 1061.8750	|2.633250	| 0.8173645|
|SRR8146970 |1090	|1333.0794	|2.497788	| 0.7651380|
|SRR8146971	|1279	|1618.3000	|2.435393	| 0.6800527|
|SRR8146972	|972	|1145.3544	|3.471903	| 0.9288290|
|SRR8146973	|1116	|1354.7100	|3.762623	| 0.9573695|
|SRR8146974	|1112	|1339.1750	|3.319138	| 0.9322468|
|SRR8146975	|803	|940.0923	  |3.429546	| 0.9390154|
|SRR8146976	|797	|1015.3051	|3.490182	| 0.9180766|
|SRR8146977	|953	|1207.0806	|3.401755	| 0.9330893|
|SRR8146978 |1082	|1332.0303	|3.173553	| 0.9053566|

**Table 3. Alpha diversity indexes data for all samples.**

<img width="2350" height="899" alt="20260321_PCoA_bray_jaccard" src="https://github.com/user-attachments/assets/65cd2152-6d43-4548-87df-ba22397e3e85" />

**Figure 3. PCoA plots of dissimilarity across different types of beta diversity measures.** Principal Coordinate Analysis (PCoA) plots were generated using two different beta diversity measures for: A. Bray-Curtis where variance can be explained by 32.2% and 21.5% for each PCoA and B. Jaccard where 24.8% and 17.8% for each PCoA can be explained. PERMANOVA p values of 0.151 and 0.168 were observed for Bray-Curtis and Jaccard, respectively.

<br>

<img width="1348" height="1126" alt="20260320_diff_abundance" src="https://github.com/user-attachments/assets/43ef5a14-cd58-4f24-869d-a3cf3407f08e" />

**Figure 4. Differentially abundant bacterial species between vegan and omnivore diets.** Plot shows that *Gordonibacter urolithinfaciens, Ernoma phocaeenis,* and *Collinsella sterconis* are the most differentially abundant in vegans while Prevotella species are most differentially abundant in omnivore diets. The smallest p values were extracted.

## Discussion

## References
|
