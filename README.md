Workflow detailing the steps for shotgun sequencing comparing between vegans diet and omnivore diets.

## Introduction
The human gut microbiome has been suggested to play an important role in human health such as producing vital nutrients and vitamins and aiding with breaking down complex carbohydrates and fibres. Several bacterial phyla that include *Firmicutes, Bacteroidetes, Actinobacteria, Proteobacteria, Fusobacteria*, and *Verrucomicrobia* generally make up the gut microbiota (Hou et al., 2022). Notably, the gut microbiome has shown to be highly influenced by the type of diet consumed. In modern times, studies have suggested diets rich in fruits and vegetables or a vegan diet have different bacterial compositions as compared to diets rich in meat proteins or omnivore diets (Fackelmann et al., 2025). Specifically, vegan diets tend to be richer in fibre-based foods which have resulted in more abundance of bacterial taxa such as Bacteroidetes and Prevotella (Tomova et al., 2019). Recent advances in shotgun metagenomics sequencing have revolutionized the field of microbial studies and can help us gain insight into the different types of bacterial species or the microbial diversity present in the human gut microbiome. In fact, shotgun metagenomics has the advantage of characterizing bacteria that have been historically difficult to cultivate in laboratory. Therefore, the overall goal of shotgun metagenomics is to provide a comprehensive overview of the types of microbials organisms and composition that are present in a certain environment or sample. From this, we can infer their functional roles and their interactions with the host. Specifically for this analysis, the objective is to perform a comparison of the gut microbiome diversity, and the different types of bacteria present between vegan and omnivore diets. The current hypothesis of vegan diets harbouring more diversity and bacteria specialized in digesting fibre and plant matter compared to omnivores will be tested.

There are many bioinformatics tools specialized for analyzing shotgun metagenomics data. Some examples of metagenomic profilers include Kraken2 (Wood et al., 2019), CLAssifier based on Reduced K-mers (CLARK) (Ounit et al., 2015), Kaiju (Menzel et al., 2016), and MetaPhlan 4 (Blanco-Míguez et al., 2023). Each tool uses a different classification approach in which both Kraken2 and CLARK use a DNA-based classification that examines and matches exact k-mers, Kaiju is a protein-based alignment, and MetaPhlAn 4 is marker gene-based aligner. Recent benchmarking studies have compared the efficacy and accuracy of these methods. Overall, Kraken2 was determined to be the superior classifier as it has been reported to have the highest F1 score in terms of performance at 0.74 while MetaPhlAn 4 and Kaiju had F1 scores of 0.41 and 0.48, respectively (Edwin et al., 2024). In addition, protein-based classifiers have the highest misclassification rates, ranging from 5-15% compared to the k-mer-based classifiers at 1-5% (Ye et al., 2019). In addition, Kraken2 when combined with Bayesian Reestimation of Abundance with KrakEN (Bracken) (Lu et al., 2017) has a more accurate abundance profile than CLARK (Ye et al., 2019). A disadvantage of Kraken2 is its longer computation time due to its higher requirement of memory compared to Kaiju which uses less than 2GB of memory. A major disadvantage of Kraken2 is its false positive rates; with as high as 5%. KrakenUniq (Breitwieser et al., 2018) has been proposed as an alternative to Kraken2 to mitigate the false positive problem however, it is very memory-intensive and requires at least 200GB of RAM even with the usage of smaller databases such as the Standard database. To mitigate the RAM limitations, bioinformaticians have increased the confidence level to 0.15 instead of the default 0 for Kraken2. Notably, the confidence score is one of the special features of Kraken2. In fact, it has been reported that a confidence level of 0.15 is a reasonable threshold for classification (Nyström-Persson et al., 2025). Any higher confidence levels can also be used and leads to more accuracy but at the cost of reduced sensitivity and increased numbers of unclassified reads (Liu et al., 2024).

One of the biggest decisions when performing shotgun metagenomics is database selection. The completeness and accuracy of the database is critical to ensure accurate and comprehensive taxonomic classification of microbial metagenomic sequencing. Examples of databases include, but not limited to, Minikraken, the subsampled 8GB and 16GB Standard databases, the full Standard database, and the core_nt database. Benchmarking studies employing Kraken2 with simulated metagenomic datasets have compared the taxonomic performance of each database. In fact, the use of smaller databases such as Minikraken and the Standard-16GB database led to approximately 0% of classified reads when confidence scores of 0.40 were utilized (Liu et al., 2024). In contrast, the larger databases such as the core nt database and the GTDB r202 database had classification reads of 93% and 81%, respectively, and at the confidence score of 0.4 (Liu et al., 2024). Therefore, if space and resources permit, larger databases, especially the core nt database, are usually better performers for taxonomic classification compared to smaller or subsampled databases due to their better coverage of all current known sequences.

Statistical software is also important to consider when performing shotgun metagenomics to determine if there are differentially abundant microbial species between different groups. Many programs exist such as edgeR (Robinson et al., 2010), the newer ANOVA-Like Differential Expression 3 (Aldex3) (Nixon et al., 2025), which is the successor to Aldex2, and the Analysis of Compositions of Microbiomes with Bias Correction 2 (ANCOM-BC2) (Lin & Peddada, 2020, 2024). These methods differ in their statistical methods in that edgeR uses a negative bionomical distribution model (Robinson et al., 2010), Aldex3 uses the Dirichlet–multinomial approach (Nixon et al., 2025) and ANCOM-BC2 uses a log-linear regression to perform differential abundance (Lin & Peddada, 2020, 2024). Benchmarking studies compared different types of differential abundance softwares and reported that Aldex2 (note that no studies have compared Aldex3’s performance yet) and ANCOM-BC2 are the most robust and agree the best when reporting results compared to other software such as edgeR where they reported more false positive discovery rates (FDR) (Nearing et al., 2022). In fact, edgeR was reported not to handle compositional data as well as ANCOM-BC2 which has the advantage of bias correction. In fact, ANCOM-BC2 was reported to be one of the most conservative tools as it identifies FDR of 0-5% compared to edgeR’s FDR ranging from 10 - 30% (Nearing et al., 2022).

## Metholodogy

Based on previously evaluated tools, Kraken2/Bracken, and ANCOM-BC2 will be used to perform analysis on shotgun metagenomics. 

### Obtaining dataset and quality control of metegenomic reads
Paired end fastq files of SRR8146968, SRR8146970- SRR8146978 samples were downloaded from the National Centre for Biotechnology Information (NCBI) (https://www.ncbi.nlm.nih.gov/sra/?term=SRP126540) with 5 samples each downloaded for 2 diet groups: vegan and omnivores for a total of 10 samples. Samples were taken from residents living in Turin, Italy. The data was previously generated using the Illumina HiSeq platform (De Filippis et al., 2019). All data files and modules were loaded on the Digital Research Alliance of Canada with the usage of the Nibi high-performance computing cluster unless otherwise stated. The sra-toolkit (v3.0.9) (Leinonen et al., 2011) was used to download the data with the “prefetch” command employed to download the SRA files. Consequently, the “fasterq-dump” command was employed with the “--split-files” flag to convert the SRA files to fastq files, producing paired end reads (i.e. SRR8146972_1.fastq  SRR8146972_2.fastq). FastQC (v0.12.1) (Andrews, 2010) was used to assess for quality of the data files to determine adapter contamination. fastp (v1.0.1) (Chen et al., 2018) was employed to solely remove residual sequencing adapters. 

### Taxonomic abundance estimation
Taxonomic classifications were  performed using  Kraken2 (v2.1.6) (Wood et al., 2019) with the following flags: “--confidence 0.15”, “--threads 32”, “--report” to generate a kraken report, and  “--paired data” with the core_nt database reference database (https://genome-idx.s3.amazonaws.com/kraken/k2_core_nt_20251015.tar.gz).  Bracken (v3.0) (Lu et al., 2017) was employed to estimate species abundance with the following flags employed: “-r 150” due to the data having 150bp short-read length, “-l S” for species classification, “-t 0”, and input from the kraken reports. kraken-biom (v 1.2.0) (Dabdoub, 2016) was downloaded on bioconda due to the module not being available on the Nibi cluster and was employed to convert the species_report files to one combined biom file with the flag “--fmt json" to convert to json format. Biom files were read into RStudio (v4.5.1 2025-06-13) using the biomformat (v1.38.0) (P. McMurdie & Paulson, 2025) package. Host DNA was not removed prior to classification but taxa that were classified as non-microbial organisms (such as chordata) were filtered out prior to downstream analysis using the microviz (v0.13.0) package (Barnett et al., 2021) via the tax_select() function as bacterial reads can be misidentified as eukaryotic taxa reads in shotgun metagenomics (Lind & Pollard, 2021). The phyloseq package (v1.54.2) (P. J. McMurdie & Holmes, 2013) was used to calculate relative abundance at the phylum and species level. 

### Alpha and beta diversity
The aforementioned phyloseq package (P. J. McMurdie & Holmes, 2013) was used to perform alpha and beta diversity. Alpha diversity metrics employed included observed, Chao1, Shannon’s index and Simpson’s index and were plotted as a scatter plot. For beta diversity, ordination via the Principal Coordinate Analysis (PCoA) method was employed using both Bray-Curtis’s and Jaccard’s distance where binary was set to TRUE for the latter. The Permutational Multivariate Analysis of Variance Using Distance Matrices (PERMANOVA) was calculated using the adonis2() function from the vegan (v 2.7-3) package (Oksanen et al., 2026). 

### Differentially abundant taxa estimation
ANCOM-BC2 (v2.12.0) (Lin & Peddada, 2024) was used for differential abundant estimation. Statistical modeling was performed using the formula ~diet_type with the Benjamini-Hochberg (BH) employed for adjusted p values to control for false positives and both pseudo_sens and struc_zero parameters set to TRUE. The top differentially abundant species were extracted by the lowest adjusted p values and with absolute log fold change values > 1. Plots were generated using the ggplot2 (v4.0.2) package (Wickham, 2016).

### Code availability
The bash and R code used for data preprocessing and analysis used for this analysis can be accessed at (https://github.com/liona-vu/Shotgun_metagenomics).

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

<img width="1758" height="1211" alt="20260321_relative_abundance_barplot" src="https://github.com/user-attachments/assets/2158fe04-173f-4308-8fcb-de6cf840d302" />

**Figure 1. Relative abundance of gut microbiota bacterial at the phylum level.** Top 15 phyla from n=10 samples were plotted as determined by the number of raw counts. Each bar was grouped based on their diet type. Phyla not determined as the top 15 were classified as Other. Samples are labeled on the x-axis with their corresponding SRR sample number. Phyla *Bacillota* and *Bacteroidota* make up most of the relative abundance across all 10 samples. Bars in the legend are coloured by phylum.

<br>

<img width="1083" height="884" alt="20260321_boxplot_top_10_species" src="https://github.com/user-attachments/assets/9e569978-b41c-44bf-a2f7-b20dcd1ad4d8" />

**Figure 2. Boxplot of the top identified bacterial species between different diet types.** Top 10 abundance of bacterial species were determined between both omnivore and vegan diets and was plotted as relative abundance, with the bacterial species labeled on the x-axis. Each dot represents an individual sample outlier that was determined to have the corresponding abundant bacterial species.

<br>

### Certain bacterial taxa are prominent in either vegan or omnivores diets
After determining that Kraken2 has identified sufficient percent of classified reads ranging from 39% to 79% to perform analysis (Table 2), relative abundance was examined. *Bacillota* and *Bacteroidota* together accounted for the majority of relative abundance across all samples (Figure 1). Interestingly, one sample, SRR8146974 appeared as an outlier. In fact, it has an unusually higher relative abundance of the *Verrucomicrobiota* phylum compared to the rest of the samples which may reflect individual microbiome variability rather than diet (Figure 1). Notably, this was in the vegan diet sample. Past the phyla level, the bacterial taxa at the species level were investigated and to determine if there were differences between vegan and omnivore diets. Interestingly, *Segatella copri (S. copri)* is observed to be greater in abundance in omnivore diets than vegans (Figure 2). On the contrary, the *Oscillospiraceae bacterium* have higher relative abundance in vegans than omnivores (Figure 2). In fact, *Oscillospiraceae* has been reported to help decompose plant material(Vedel et al., 2023). Finally, *Faecalibacterium prausnitzii*, one of the most abundant bacterial species in the human gut, was determined to be one of the top abundant bacterial species in both diets (Lopez-Siles et al., 2017). Together, these findings suggests that there are certain types of bacterial taxa at both the phyla and species level that are influenced by vegan or omnivore diets.

<br>
 
<img width="2291" height="1943" alt="20250317_alpha_diversity_labeled" src="https://github.com/user-attachments/assets/39cbbc7f-6233-46b3-bfe2-cf13b26791ec" />

<br>

**Figure 3. Alpha diversity index plots for all SRR samples.** Measures include A. observed, B. Chao1, C. Shannon and D. Simpson. Diet types for omnivores (salmon) and vegans (turquoise) are shown. Alpha diversity metrics did not reveal clear patterns between the diet groups. 

<br>

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

**Table 3. Alpha diversity indexes measures for all samples.** The observed, Chao1, Shannon, and Simpson diversity indexes were determined for all samples in this analysis. A higher value for Observed, Chao1, Shannon, and Simpson generally indicate greater diversity.

<br>

<img width="2043" height="982" alt="20260322_PCoA_plot_bray_jaccard_legit" src="https://github.com/user-attachments/assets/4e9f7e40-ed7e-46a8-8bc6-21b9575bd639" />

**Figure 4. Principal Coordinate Analysis (PCoA) plots of dissimilarity across different types of diets.** PCoA plots were generated to display the spatial distribution of each sample coloured by diet type (vegan or omnivore) using two different beta diversity measures: A. Bray-Curtis where variance can be explained by 32.2% and 21.5% for PCoA 1 and PCoA 2, and B. Jaccard where 15.7% and 13.9% for each PCoA can be explained. PERMANOVA P values of 0.162 and 0.166 were observed for Bray-Curtis and Jaccard methods, respectively. Closer points indicate more similar microbial communities while farther points indicate greater dissimilarity.

<br>

### Alpha and beta diversity are not significant between vegan and omnivore diets
Next, alpha diversity measures were calculated for all samples. Observed represents how many distinct taxa or in this case, operational taxonomic units (OTUs) present in a sample. Chao1 is an indicator for species richness, with higher values indicating higher diversity (Chao, 1987). A higher index for Shannon indicates a more diverse community while the Simpson’s index ranges from 0 to 1. Since phyloseq was used to calculate the Simpson index, 0 indicate low diversity and 1 indicating higher diversity. Unlike the previous figures, results show no discernable pattern between the two diet types for all 4 metrics of Observed, Chao1, Shannon, and Simpson indexes (Figure 3). Notably, there was no big difference between vegan and omnivores for the Shannon’s and Simpson’s indexes as both ranged between 2.4 – 3.8, and 0.68 – 0.93, respectively (Table 3). Interestingly, the vegan diet sample SRR8146973 has the highest 3 metrics, with Observed, Shannon’s and Simpson’s indexes, having values of 1116, 3.762 and 0.9573695 which is consistent with what was observed in figure 1 with the greater number of different phylum types in this sample compared to omnivores samples.

<br>

<img width="1168" height="982" alt="20260322_differential_abundance_final" src="https://github.com/user-attachments/assets/060156d4-d720-489a-993e-74941821b2ee" />

**Figure 5. Differentially abundant bacterial species between vegan and omnivore diets.** Top 30 taxa with differential abundance difference in vegans and omnivore diets. Bacterial species were extracted with the smallest p values and filtered with absolute logfoldchange > 1. Negative log fold changes were associated with omnivore diets and positive log fold changes were associated with vegan diets. Error bars represent standard error (SE).

<br>

### Differentially abundant bacterial species between vegan and omnivore diets
To determine which species are differentially abundant, the ANCOM-BC2 (Lin & Peddada, 2024) statistical model was employed and determined 226 differentially abundant taxa. Unfortunately, none of the adjusted p values were considered statistically significant as the adjusted p values were all greater than 0.05 (Figure 5). Therefore, the smallest p values were chosen and analysis proceeded from there. There were certain taxa that were more present in vegan diets. In fact, *Gordonibacter urolithinfaciens (G. urolithinfaciens) * (q = 0.85) has the highest log fold change and *Gordonibacter massiliensis* showed both the lowest p value (q = 0.61) and one of the largest log fold changes suggesting as the most promising bacterium to investigate. On the contrary, the *Prevotella* species were some of the most differentially abundant species in omnivore diets, consistent with the finding from figure 2 where the most abundant bacterial species was also *S. copri* in omnivore diets. Overall, there were differentially abundant species between vegans and omnivores, but they were not statistically significant.


## Discussion
Shotgun metagenomics has revolutionized the field of microbiology and our understanding of the human gut microbiome without traditional laboratory culturing. Specifically, it allows scientists to determine whether factors such as diet play a significant role in shaping the diversity or microbial composition of the human gut microbiome. However, not all microbial organisms have been discovered. During quality control, the lower percentage of classified reads can be due several factors: rare or extremely low abundance of microbes and/or the completeness of the database in which undiscovered bacterial species have yet to be added. Thus, some of the reads remained unclassified. In fact, the percent of classified reads in this analysis is lower than what had been reported in literature where it has been as high as 90% in simulated datasets (Liu et al., 2024). 

Literature has reported that vegan diets often harbour higher diversity of bacterial species or composition in the human gut microbiome compared to omnivore diets. By solely employing visual assessment, there was a presence of bacterial phyla that were present in the vegan samples that were either not found or found in low abundance in the omnivore diet, supporting this notion. (Figure 1).  Unfortunately, none of the results in this analysis was reported to be significant in all aspects such as alpha or beta diversity, PERMANOVA analysis, or the ANCOM-BC2 statistical modeling of differential abundance between vegan and omnivore diets. Notably, all p values were greater than the 0.05 cutoff. One reasonable explanation is due to the low sample size of 10. Therefore, increasing the number of samples may yield higher statistical power to detect a significant difference. In fact, a study reported using 21,561 individuals to investigate the effects of vegan and omnivore in addition to vegetarian diets on the gut microbiome which is roughly 2100 times more samples (Fackelmann et al., 2025). 

Nonetheless, despite the statistical significance limitations, there are still some key results that can yield insights into the human gut microbiome and the influence of diet. For example, in Figure 1, the two major phylum that made up the gut microbiome included the *Bacillota* (also known as *Firmicutes*) and *Bacteroidota*, which is consistent with what has been reported in literature (Hou et al., 2022). Notably, the Firmicutes: Bacteroidota ratio is often used a key health marker metric where vegans usually exhibit lower *Firmicutes: Bacteroidota* ratio than omnivores (Tomova et al., 2019). Surprisingly, the findings contradict this notion as there are more Bacteroidota than *Firmicutes/Bacillota* in most of the vegan samples. It is likely that other factors such as environmental factors may play a role in microbiome composition.

When analyzed at the species level, a certain bacterium, *S. copri* (previously known as *Prevotella copri (P. copri*) was shown to have higher abundance in omnivore samples than vegan samples (Figure 2). In fact, *S. copri* is part of the *Bacteroidota* phylum. Notably, this was also the case when performing differential abundant analysis where *S. copri* was determined to be differentially abundant in omnivore diets (Figure 5). This was an unusual finding since *S. copri* has been reported to be associated with diets higher in fibre such as vegan diets (Ley, 2016). However, one study reported *S. copri* was not identified as a key characteristic in vegan diets, corroborating with the results in this analysis (Fackelmann et al., 2025). Perhaps this has to do with the diets prevalent in Turin, Italy, where the residents may follow a Mediterranean diet which is characterized by fresh vegetables (high fibre), olive oil, and pasta. Therefore, omnivores in this cohort may consume more plant-based foods compared to a Western diet and thus, no differences were determined.

On top of *S. copri* species, other species under the *Prevotella* genus such as uncultured *Prevotella* was found to be differentially abundant in vegans, and *Prevotella heparinolytica* in omnivores (Figure 5). In fact, the role of the genus *Prevotella* in the context of human health have been varied. Some studies have reported that certain *Prevotella* species such as *Prevotella copri* acts as a beneficial microbe where it can produce short-chain fatty acids which protects the mucosal barrier and may reduce inflammation. Consequently, inflammatory diseases such as rheumatoid arthritis and type I and II diabetes were associated with lower amounts of short-chain fatty acids and lower abundance of *P. copri* (Bedarf et al., 2017). In contrast, *P. copri* has been associated with insulin resistance and impaired glucose uptake (Pedersen et al., 2016). In vegan samples, *G. urolithinfaciens* was indicated as the most differentially abundant species (Figure 5). In fact, *G. urolithinfaciens* was determined to metabolize ellagic acid which is a phenolic compound that is abundantly found in fruits into urolithins. Notably, urolithins has implications in human health where it can exhibit anti-inflammatory and anticarcinogenic effects (Selma et al., 2014). This makes sense given that vegans often have diets high in fruits. 

In summary, the gut microbiome was compared between vegan and omnivore diets using key metrics such as alpha and beta diversity, differential abundance, and relative abundance. While statistical significance was not achieved, interesting biological trends were observed such as *S. copri* abundance in omnivore diets and the dominance of *Bacillota/Firmicutes* and *Bacteroidota* across all samples which is consistent with established literature. Future studies should include larger sample sizes to determine significant differences between different types of diets on the gut microbiome.


## References
Andrews, S. (2010). FastQC: A quality control tool for high throughput sequence data (0.12.1).

Barnett, D., Arts, I., & Penders, J. (2021). microViz: an R package for microbiome data visualization and statistics. Journal of Open Source Software, 6(63), 3201. https://doi.org/10.21105/joss.03201

Bedarf, J. R., Hildebrand, F., Coelho, L. P., Sunagawa, S., Bahram, M., Goeser, F., Bork, P., & Wüllner, U. (2017). Functional implications of microbial and viral gut metagenome changes in early stage L-DOPA-naïve Parkinson’s disease patients. Genome Medicine, 9(1), 39. https://doi.org/10.1186/s13073-017-0428-y

Blanco-Míguez, A., Beghini, F., Cumbo, F., McIver, L. J., Thompson, K. N., Zolfo, M., Manghi, P., Dubois, L., Huang, K. D., Thomas, A. M., Nickols, W. A., Piccinno, G., Piperni, E., Punčochář, M., Valles-Colomer, M., Tett, A., Giordano, F., Davies, R., Wolf, J., … Segata, N. (2023). Extending and improving metagenomic taxonomic profiling with uncharacterized species using MetaPhlAn 4. Nature Biotechnology, 41(11), 1633–1644. https://doi.org/10.1038/s41587-023-01688-w

Breitwieser, F. P., Baker, D. N., & Salzberg, S. L. (2018). KrakenUniq: confident and fast metagenomics classification using unique k-mer counts. Genome Biology, 19(1), 198. https://doi.org/10.1186/s13059-018-1568-0

Chao, A. (1987). Estimating the population size for capture-recapture data with unequal      catchability. Biometrics, 43(4), 783–791. https://doi.org/10.2307/2531532

Chen, S., Zhou, Y., Chen, Y., & Gu, J. (2018). fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics, 34(17), i884–i890. https://doi.org/10.1093/bioinformatics/bty560

Dabdoub, S. (2016). kraken-biom: Enabling interoperative format conversion for Kraken results (Version 1.2) [Software]. Available at https://github.com/smdabdoub/kraken-biom.

De Filippis, F., Pasolli, E., Tett, A., Tarallo, S., Naccarati, A., De Angelis, M., Neviani, E., Cocolin, L., Gobbetti, M., Segata, N., & Ercolini, D. (2019). Distinct Genetic and Functional Traits of Human Intestinal <em>Prevotella copri</em> Strains Are Associated with Different Habitual Diets. Cell Host & Microbe, 25(3), 444-453.e3. https://doi.org/10.1016/j.chom.2019.01.004

Edwin, N. R., Fitzpatrick, A. H., Brennan, F., Abram, F., & O’Sullivan, O. (2024). An in-depth evaluation of metagenomic classifiers for soil microbiomes. Environmental Microbiome, 19(1), 19. https://doi.org/10.1186/s40793-024-00561-w

Fackelmann, G., Manghi, P., Carlino, N., Heidrich, V., Piccinno, G., Ricci, L., Piperni, E., Arrè, A., Bakker, E., Creedon, A. C., Francis, L., Capdevila Pujol, J., Davies, R., Wolf, J., Bermingham, K. M., Berry, S. E., Spector, T. D., Asnicar, F., & Segata, N. (2025). Gut microbiome signatures of vegan, vegetarian and omnivore diets and associated health outcomes across 21,561 individuals. Nature Microbiology, 10(1), 41–52. https://doi.org/10.1038/s41564-024-01870-z

Hou, K., Wu, Z.-X., Chen, X.-Y., Wang, J.-Q., Zhang, D., Xiao, C., Zhu, D., Koya, J. B., Wei, L., Li, J., & Chen, Z.-S. (2022). Microbiota in health and diseases. Signal Transduction and Targeted Therapy, 7(1), 135. https://doi.org/10.1038/s41392-022-00974-4

Leinonen, R., Sugawara, H., Shumway, M., & Collaboration, on behalf of the I. N. S. D. (2011). The Sequence Read Archive. Nucleic Acids Research, 39(suppl_1), D19–D21. https://doi.org/10.1093/nar/gkq1019

Ley, R. E. (2016). Prevotella in the gut: choose carefully. Nature Reviews Gastroenterology & Hepatology, 13(2), 69–70. https://doi.org/10.1038/nrgastro.2016.4

Lin, H., & Peddada, S. Das. (2020). Analysis of compositions of microbiomes with bias correction. Nature Communications, 11(1), 3514. https://doi.org/10.1038/s41467-020-17041-7

Lin, H., & Peddada, S. Das. (2024). Multigroup analysis of compositions of microbiomes with covariate adjustments and repeated measures. Nature Methods, 21(1), 83–91. https://doi.org/10.1038/s41592-023-02092-7

Lind, A. L., & Pollard, K. S. (2021). Accurate and sensitive detection of microbial eukaryotes from whole metagenome shotgun sequencing. Microbiome, 9(1), 58. https://doi.org/10.1186/s40168-021-01015-y

Liu, Y., Ghaffari, M. H., Ma, T., & Tu, Y. (2024). Impact of database choice and confidence score on the performance of taxonomic classification using Kraken2. ABIOTECH, 5(4), 465–475. https://doi.org/https://doi.org/10.1007/s42994-024-00178-0

Lopez-Siles, M., Duncan, S. H., Garcia-Gil, L. J., & Martinez-Medina, M. (2017). Faecalibacterium prausnitzii: from microbiology to diagnostics and prognostics. The ISME Journal, 11(4), 841–852. https://doi.org/10.1038/ismej.2016.176

Lu, J., Breitwieser, F. P., Thielen, P., & Salzberg, S. L. (2017). Bracken: estimating species abundance in metagenomics data. PeerJ Computer Science, 3, e104. https://doi.org/10.7717/peerj-cs.104

McMurdie, P. J., & Holmes, S. (2013). phyloseq: An R Package for Reproducible Interactive Analysis and Graphics of Microbiome Census Data. PLOS ONE, 8(4), e61217-. https://doi.org/10.1371/journal.pone.0061217

McMurdie, P., & Paulson, J. (2025). biomformat: An interface package for the BIOM file format (R package version 1.38.0). Bioconductor. https://doi.org/10.18129/B9.BIOC.BIOMFORMAT

Menzel, P., Ng, K. L., & Krogh, A. (2016). Fast and sensitive taxonomic classification for metagenomics with Kaiju. Nature Communications, 7(1), 11257. https://doi.org/10.1038/ncomms11257

Nearing, J. T., Douglas, G. M., Hayes, M. G., MacDonald, J., Desai, D. K., Allward, N., Jones, C. M. A., Wright, R. J., Dhanani, A. S., Comeau, A. M., & Langille, M. G. I. (2022). Author Correction: Microbiome differential abundance methods produce different results across 38 datasets. Nature Communications, 13(1), 777. https://doi.org/10.1038/s41467-022-28401-w

Nixon, M. P., Gloor, G. B., & Silverman, J. D. (2025). Incorporating scale uncertainty in microbiome and gene expression analysis as an extension of normalization. Genome Biology, 26(1), 139. https://doi.org/10.1186/s13059-025-03609-3

Nyström-Persson, J., Bapatdhar, N., & Ghosh, S. (2025). Precise and scalable metagenomic profiling with sample-tailored minimizer libraries. NAR Genomics and Bioinformatics, 7(2), lqaf076. https://doi.org/10.1093/nargab/lqaf076

Oksanen, J., Simpson, G., Blanchet, G., Kindt, R., Legendre Pierre, Minchin Peter, O’Hara R.B., Solymos, P., Stevens, M., Szoecs, E., Wagner Helene, Barbour Matt, Bedward, M., Bolker, B., Borcard, D., Borman, T., Carvalho, G., Chirico, M., De Caceres, M., … Weedon, J. (2026). vegan: Community Ecology Package. (R package version 2.8-0). https://vegandevs.github.io/vegan/

Ounit, R., Wanamaker, S., Close, T. J., & Lonardi, S. (2015). CLARK: fast and accurate classification of metagenomic and genomic sequences using discriminative k-mers. BMC Genomics, 16(1), 236. https://doi.org/10.1186/s12864-015-1419-2

Pedersen, H. K., Gudmundsdottir, V., Nielsen, H. B., Hyotylainen, T., Nielsen, T., Jensen, B. A. H., Forslund, K., Hildebrand, F., Prifti, E., Falony, G., Le Chatelier, E., Levenez, F., Doré, J., Mattila, I., Plichta, D. R., Pöhö, P., Hellgren, L. I., Arumugam, M., Sunagawa, S., … Pedersen, O. (2016). Human gut microbes impact host serum metabolome and insulin sensitivity. Nature, 535(7612), 376–381. https://doi.org/10.1038/nature18646

Robinson, M. D., McCarthy, D. J., & Smyth, G. K. (2010). edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics, 26(1), 139–140. https://doi.org/10.1093/bioinformatics/btp616

Selma, M. V, Beltrán, D., García-Villalba, R., Espín, J. C., & Tomás-Barberán, F. A. (2014). Description of urolithin production capacity from ellagic acid of two human intestinal Gordonibacter species. Food & Function, 5(8), 1779–1784. https://doi.org/10.1039/C4FO00092G

Tomova, A., Bukovsky, I., Rembert, E., Yonas, W., Alwarith, J., Barnard, N. D., & Kahleova, H. (2019). The Effects of Vegetarian and Vegan Diets on Gut Microbiota. Frontiers in Nutrition, Volume 6-2019. https://doi.org/10.3389/fnut.2019.00047

Vedel, G., Triadó-Margarit, X., Linares, O., Moreno-Rojas, J. M., la Peña, E. de, García-Bocanegra, I., Jiménez-Martín, D., Carranza, J., & Casamayor, E. O. (2023). Exploring the potential links between gut microbiota composition and natural populations management in wild boar (Sus scrofa). Microbiological Research, 274, 127444. https://doi.org/https://doi.org/10.1016/j.micres.2023.127444

Wickham, H. (2016). ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York.

Wood, D. E., Lu, J., & Langmead, B. (2019). Improved metagenomic analysis with Kraken 2. Genome Biology, 20(1), 257. https://doi.org/10.1186/s13059-019-1891-0

Ye, S. H., Siddle, K. J., Park, D. J., & Sabeti, P. C. (2019). Benchmarking Metagenomics Tools for Taxonomic Classification. Cell, 178(4), 779–794. https://doi.org/10.1016/j.cell.2019.07.010
 

