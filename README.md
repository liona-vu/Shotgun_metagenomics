Workflow detailing the steps for shotgun sequencing comparing between vegans diet and omnivore diets.

## Introduction
The human gut microbiome has been suggested to play an important role in human health such as producing vital nutrients and vitamins and aiding with breaking down complex carbohydrates and fibres. Several bacterial phyla that include *Firmicutes, Bacteroidetes, Actinobacteria, Proteobacteria, Fusobacteria*, and *Verrucomicrobia* generally make up the gut microbiota (Hou et al., 2022). Notably, the gut microbiome has shown to be highly influenced by the type of diet consumed. In modern times, studies have suggested diets rich in fruits and vegetables (or a vegan diet) have different bacterial compositions as compared to diets rich in meat proteins or omnivore diets (Fackelmann et al., 2025). Recent advances in shotgun metagenomics sequencing have revolutionized the field of microbial studies and can help us gain insight into the different types of bacterial species or the microbial diversity present in the human gut microbiome. In fact, shotgun metagenomics has the advantage of characterizing bacteria that have been historically difficult to cultivate in laboratory. Therefore, the overall goal of shotgun metagenomics is to provide a comprehensive overview of the types of microbials organisms and composition that are present in a certain environment or sample. From this, we can infer their functional roles and their interactions with the host environment. Specifically for this analysis, the objective is to perform a comparison of the gut microbiome diversity, and the different types of bacteria present between vegan and omnivore diets. 

There are many bioinformatics tools specialized for analyzing shotgun metagenomics data. Some examples of metagenomic profilers include Kraken2 (Wood et al., 2019), CLAssifier based on Reduced K-mers (CLARK) (Ounit et al., 2015), Kaiju (Menzel et al., 2016), and MetaPhlan 4 (Blanco-Míguez et al., 2023). Each tool uses a different classification approach in which both Kraken2 and CLARK use a DNA-based classification that examines and matches exact k-mers, Kaiju is a protein-based alignment, and MetaPhlAn 4 is marker gene-based aligner. Recent benchmarking studies have compared the efficacy and accuracy of these methods. Overall, Kraken2 was determined to be the superior classifier as it has been reported to have the highest F1 score in terms of performance at 0.74 while MetaPhlAn 4 and Kaiju had F1 scores of 0.41 and 0.48, respectively (Edwin et al., 2024). In addition, protein-based classifiers have the highest misclassification rates, ranging from 5-15% compared to the k-mer-based classifiers at 1-5% (Ye et al., 2019). In addition, Kracken2 when combined with Bayesian Reestimation of Abundance with KrakEN (Bracken) (Lu et al., 2017) has a more accurate abundance profile than CLARK (Ye et al., 2019). A disadvantage of Kraken2 is its longer computation time due to its higher requirement of memory compared to Kaiju which uses less than 2GB of memory. A major disadvantage of Kraken2 is its false positive rates; with as high as 5%. KrakenUniq (Breitwieser et al., 2018) has been proposed as an alternative to Kraken2 to mitigate the false positive problem however, it is very memory-intensive and requires at least 200GB of RAM even with the usage of smaller databases such as the Standard database. To mitigate the RAM limitations, bioinformaticians have increased the confidence level to 0.15 instead of the default 0 for Kracken2. Notably, the confidence score is one of the special features of Kraken2. In fact, it has been reported that a confidence level of 0.15 is a reasonable threshold for classification (Nyström-Persson et al., 2025). Any higher confidence levels can also be used and leads to more accuracy but at the cost of reduced sensitivity and increased numbers of unclassified reads (Liu et al., 2024).

One of the biggest decisions when performing shotgun metagenomics is database selection. The completeness and accuracy of the database is critical to ensure accurate and comprehensive taxonomic classification of microbial metagenomic sequencing. Examples of databases include, but not limited to, Minikraken, the subsampled 8GB and 16GB Standard databases, the full Standard database, and the core_nt database. Benchmarking studies employing Kraken2 with simulated metagenomic datasets have compared the taxonomic performance of each database. In fact, the use of smaller databases such as Minikraken and the Standard-16GB database led to approximately 0% of classified reads when confidence scores of 0.40 were utilized (Liu et al., 2024). In contrast, the larger databases such as the core nt database and the GTDB r202 database had classification reads of 93% and 81%, respectively, and at the confidence score of 0.4 (Liu et al., 2024). Therefore, if space and computational resources permit, larger databases, especially the core nt database, are usually better performers for taxonomic classification compared to smaller or subsampled databases due to their coverage of all known current sequences.

## Metholodogy

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

**Figure 1. Relative abundance of gut microbiota bacterial at the phylum level.** Top 15 phyla from n=10 samples were plotted as determined by the number of raw counts while the other phyla were classified as Other. Samples are labeled on the x-axis as their corresponding SRR number. *Bacillota* and *Bacteroidata* make up most of the abundance phyla across all 10 samples.

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

<br>

<img width="1348" height="1126" alt="20260320_diff_abundance" src="https://github.com/user-attachments/assets/43ef5a14-cd58-4f24-869d-a3cf3407f08e" />

**Figure 4. Differentially abundant bacterial species between vegan and omnivore diets.** Plot shows that *Gordonibacter urolithinfaciens, Ernoma phocaeenis,* and *Collinsella sterconis* are the most differentially abundant in vegans while Prevotella species are most differentially abundant in omnivore diets. The smallest p values were extracted.

## Discussion

## References
|
