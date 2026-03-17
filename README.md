Workflow detailing the steps for shotgun sequencing comparing between vegans diet and omnivore diets.

## Introduction
The human gut microbiome has been suggested to play an important role in human health such as producing vital nutrients and vitamins and aiding with breaking down complex carbohydrates and fibre. Several phyla that include Firmicutes, Bacteroidetes, Actinobacteria, Proteobacteria, Fusobacteria, and Verrucomicrobia generally make up the gut microbiota (Hou et al., 2022). Notably, the gut microbiome has shown to be highly influenced by the type of diet consumed. In modern times, studies have suggested diets rich in fruits and vegetables (or a vegan diet) have different bacterial compositions as compared to diets rich in meat proteins or omnivore diets (Fackelmann et al., 2025). Recent advances in shotgun metagenomics sequencing have revolutionized the field of microbial studies and can help us gain insight into the different types of bacterial species or the microbial diversity present in the human gut microbiome. In fact, shotgun metagenomics has the advantage of characterizing bacteria that have been historically difficult to cultivate in laboratory. Therefore, the overall goal of shotgun metagenomics is to provide a comprehensive overview of the types of microbials organisms that are present in a certain environment or sample. From this, we can infer their functional roles and their interactions with the host. Specifically for this analysis, the objective is to perform a comparison of the gut microbiome diversity, and the different types of bacteria present between vegan and omnivore diets. 

There are many bioinformatics tools specialized for analyzing shotgun metagenomics data. Some examples of metagenomic profilers include Kraken2 (Wood et al., 2019), CLAssifier based on Reduced K-mers (CLARK) (Ounit et al., 2015), Kaiju (Menzel et al., 2016), and MetaPhlan 4 (Blanco-Míguez et al., 2023). Each tool uses a different classification approach in which both Kraken2 and CLARK use a DNA-based classification that matches exact kmers, Kaiju is a protein-based alignment, and MetaPhlAn 4 is marker gene-based aligner. Recent benchmarking studies have compared the efficacy and accuracy of these methods. Overall, Kraken2 was determined to be the superior classifier as it has been reported to have the highest F1 score in terms of performance at 0.74 while MetaPhlAn 4 and Kaiju had F1 scores of 0.41 and 0.48, respectively (Edwin et al., 2024). In addition, protein-based classifiers have the highest misclassification rates, ranging from 5-15% compared to the kmer-based classifiers at 1-5% (Ye et al., 2019). In addition, Kracken2 when combined with Bayesian Reestimation of Abundance with KrakEN (Bracken) (Lu et al., 2017) has a more accurate abundance profile than CLARK (Ye et al., 2019). A disadvantage of Kraken2 is its longer computation time due to its higher requirement of memory compared to Kaiju which uses less than 2GB of memory. A major disadvantage of Kraken2 is its false positive rates; with as high as 5%. KrakenUniq (Breitwieser et al., 2018) has been proposed as an alternative to Kraken2 to mitigate the false positive problem however, it is very memory-intensive and requires at least 200GB of RAM even with the usage of smaller databases such as the Standard database. To mitigate the RAM limitations, bioinformaticians have increased the confidence level to a reasonable level instead of the default of 0 for Kracken2. In fact, it has been reported that a confidence level of 0.15 is a reasonable threshold for classification (Nyström-Persson et al., 2025). Any higher confidence levels can be used and leads to more accuracy at the cost of reduced sensitivity and increased numbers of unclassified reads (Liu et al., 2024).

## Metholodogy

## Results

|SRR Number | Sample | Spot Reads | Diet Type |
|-----------|--------|------------|-----------|
|SRR8146972 |VOV26_metag|27,306,095|Omnivore|
|SRR8146973	|VOV114_metag	|34,582,517	|Vegan|
|SRR8146974	|VOV77_metag	|35,418,839	|Vegan|
|SRR8146975	|VOV77_metag	|35,602,941	|Omnivore|
|SRR8146976	|VOV70_metag	|28,577,678 |Omnivore|
|SRR8146977	|VOV29_metag	|39,237,195	|Vegan|

Table 1. Metadata of samples used for this analysis including the SRR number, sample name, spot reads, and diet type.

|SRR Number | Classified (%) | Unclassified (%)|
|-----------|--------|------------|
|SRR8146972 |||
|SRR8146973	|	|	|
|SRR8146974	|	|	|
|SRR8146975	|	|	|
|SRR8146976	|	| |
|SRR8146977	|	|	|

Table 2. Classification percentages of each sample.

## Discussion

## References
