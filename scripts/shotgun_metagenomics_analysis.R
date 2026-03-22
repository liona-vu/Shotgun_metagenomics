# BINF 6110 Assignment 03: Shotgun Metagenomics
#
# By Liona Vu
#
# March 23 2026

##############################
#Load in libraries
##############################

library(biomformat)
library(phyloseq)
library(ggplot2)
library(tidyverse)
library(microViz)
library(ANCOMBC)
library(Biostrings)
library(vegan)
library(gridExtra)  

##############################
#Load in files
##############################



biom_file <- import_biom("combined_kraken.biom")

#Check sample names
sample_names(biom_file)
physeq_rel <- transform_sample_counts(biom_file, function(x) x / sum(x))
physeq_phy <- tax_glom(physeq_rel, taxrank = "Rank2") #rank2 is phylum, rank 3 is class and so on

#melt to plot data
df <- psmelt(physeq_phy)

# Why are there chordata and mollusca? that does not seem right
# plot phylum
ggplot(df, aes(x = Sample, y = Abundance, fill = Rank2)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(y = "Relative Abundance", x = "Sample") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  theme_minimal() 

# Assess phylum names
table(tax_table(biom_file)[,2])

# Can see that there are phylum names that don't make sense in the context of human gut microbiome
# There are phyla such as  p__Apicomplexa, p__Arthropoda, p__Artverviricota, p__Brachiopoda, p__Cercozoa, p__Chlorophyta, p__Cercozoa, , p__Chordata, p__Euglenozoa which should be removed!

#Create vector to remove unwanted phyla
remove <- c("p__Apicomplexa", "p__Arthropoda", 
            "p__Artverviricota", "p__Brachiopoda", "p__Bryozoa",
            "p__Cercozoa", "p__Chlorophyta", "p__Bacillariophyta", 
            "p__Chordata", "p__Cnidaria", "p__Cressdnaviricota",
            "p__Discosea","p__Echinodermata", "p__Ciliophora",
            "p__Evosea", "p__Euglenozoa", "p__Kitrinoviricota",
            "p__Fornicata", "p__Haptophyta", "p__Negarnaviricota",
            "p__Hemichordata","p__Heterolobosea", "p__Nematoda", 
            "p__Hofneiviricota", "p__Lenarviricota", "p__Nemertea",
            "p__Mollusca", "p__Nucleocytoviricota", "p__Parabasalia",
            "p__Oomycota",  "p__Phixviricota", "p__Xenacoelomorpha",
            "p__Pisuviricota", "p__Platyhelminthes",
            "p__Porifera","p__Preplasmiviricota",
            "p__Produgelaviricota", "p__Rhodophyta",
            "p__Streptophyta","p__Uroviricota")
           
#Filter out unwanted phyla with microviz package
filtered_biom_file <- biom_file %>%
  tax_select(tax_list = remove , deselect = TRUE)

#CHeck to see if phyla is removed
table(tax_table(filtered_biom_file)[,2])

#Still there are unclassified p__, must remove
rank_names(filtered_biom_file) #Use Rank2 since it corresponds to phylum
filtered_biom_file_2 <- subset_taxa(filtered_biom_file, Rank2 !="p__")

table(tax_table(filtered_biom_file_2)[,2])
#p__ is gone

#Rename sample names to SRR number, makes graph look cleaner
SRR_names <- c("SRR8146968", "SRR8146970" ,"SRR8146971","SRR8146972", "SRR8146973", "SRR8146974", "SRR8146975", "SRR8146976", "SRR8146977", "SRR8146978")

sample_names(filtered_biom_file_2) <- SRR_names
sample_names(filtered_biom_file_2)

########################################
# Alpha diversity
########################################
#Calculate relative abundance
physeq_rel_2 <- transform_sample_counts(filtered_biom_file_2, function(x) x / sum(x))
physeq_phy_2 <- tax_glom(physeq_rel_2, taxrank = "Rank2") 
df_2 <- psmelt(physeq_phy_2)

#Plot filtered by phylum
ggplot(df_2, aes(x = Sample, y = Abundance, fill = Rank2)) +
  geom_bar(stat = "identity", 
                   position = "stack", 
                   color = "black", 
                   linewidth = 0.2) +
  labs(y = "Relative Abundance", x = "Samples") +
  scale_x_discrete(guide = guide_axis(angle = 62.5)) +
  theme_minimal() +
  guides(fill = guide_legend(title = "Phylum Level"))
#Above graph seems too busy with over 40+ phyla

#plot top 10 phyla instead
#get top phyla
top10_phyla <- sort(taxa_sums(physeq_phy_2), 
              decreasing=TRUE)
top10_phyla <- names(top10_phyla)[1:15]
top10_phyla_names <- tax_table(physeq_phy_2)[top10_phyla, "Rank2"]

#label other phyla as "Other"
df_3 <- df_2 %>%
  mutate(new_column = ifelse(Rank2 %in% top10_phyla_names, Rank2, "Other")) %>%
  mutate(clean_names = gsub("p__", "", new_column))

#plot
ggplot(df_3, aes(x = Sample, y = Abundance, fill = clean_names)) +
  geom_bar(stat = "identity", 
           position = "stack", 
           color = "black", 
           linewidth = 0.2) +
  labs(y = "Relative Abundance", x = "Samples") +
  scale_x_discrete(guide = guide_axis(angle = 62.5)) +
  theme_minimal() +
  guides(fill = guide_legend(title = "Phylum Level"))
#from cleaner plot above, seems that bacillota and bacteroidata are the most abundant phyla

#Create metadata table
metadata_df <- data.frame(diet_type = c("Vegan", "Omnivore", "Omnivore", "Omnivore", "Vegan", "Vegan",
                                        "Omnivore", "Omnivore", "Vegan", "Vegan"))

#ADD SRR # as row names
row.names(metadata_df) <- SRR_names

#Add metadata to phyloseq file
sample_data(filtered_biom_file_2) <- metadata_df

#plot alpha diversity
plot <- plot_richness(filtered_biom_file_2, color = "diet_type", measures = c("Observed", "Chao1", "Shannon", "Simpson"), nrow = 2) +
  labs(color = "Diet Type")
plot

#grabs actual alpha diversity values
estimates <- estimate_richness(filtered_biom_file_2)
estimates

#Box whisker plot for the top most abundant species
#Grabs counts using taxa sums, sort by decreasing order 
top10 <- sort(taxa_sums(filtered_biom_file_2), 
                    decreasing=TRUE)

#grabs top 10 rownames 
top10 <- names(top10)[1:10]

#Filter by the top 10 taxas
top10_physeq <- prune_taxa(top10, filtered_biom_file_2)

top10_rel <- transform_sample_counts(top10_physeq, 
                                     function(x) x/sum(x))
#clean names
top10_df <- psmelt(top10_rel)
top10_df_2 <- top10_df %>%
  mutate(genus = gsub("g__", "", Rank6)) %>%
  mutate(species = gsub("s__", "", Rank7)) %>%
  unite(col = "Genus_species", genus, species, sep = " ")

#plot box whisker plot of top 10 species
ggplot(top10_df_2, aes(x=Genus_species, y=Abundance, fill=diet_type)) +
  geom_boxplot() +
  scale_fill_manual(values=c("Omnivore"="salmon", 
                             "Vegan"="lightblue")) +
  labs(x="Species", y="Relative Abundance",
       fill="Diet Type") +
  scale_x_discrete(guide = guide_axis(angle = 62.5)) +
  theme_minimal()

########################################
# Beta diversity
########################################
#Calculate PCoA bray distance 
pcoa_bray <- ordinate(filtered_biom_file_2, method = "PCoA", distance = "bray")
#plot_ordination(filtered_biom_file_2, pcoa_bray, color = "diet_type") + 
 # geom_point(aes(color = "diet_type"), size = 4, inherit.aes = T)

#why is colour not working 
sample_data(filtered_biom_file_2)

#NEED TO PLOT BY HAND
#get the top variances 
pca <- round((pcoa_bray[["values"]]$Relative_eig),3) *100

#grabs the eigenvalues for each sample
pcoa_bray_df <- data.frame(pca = pcoa_bray$vectors)
pcoa_bray_df <- cbind(pcoa_bray_df, metadata_df) #combine metadata with dataframe

#perform permanova
permanova_bray <- adonis2(phyloseq::distance(filtered_biom_file_2, method = "bray") ~ diet_type,
                    data = metadata_df)  
p_value_bray <- permanova_bray$`Pr(>F)`[1]
class(p_value_bray)
#Not significant, p value of 0.163

#plot again 
plot_bray <- ggplot(data = pcoa_bray_df, aes(x = pca.Axis.1, y = pca.Axis.2, color = diet_type)) +
  geom_point(size = 3) +
  xlab(paste0("[",pca[1], " %]")) +
  ylab(paste0("[", pca[2], " %]")) +
  ggtitle(label = "PCoA plot of Bray Curtis dissimilarity", paste0("P value: ", p_value_bray)) +
  theme_minimal() +
  labs(color = "Diet Type") 

#Trying NMDS
nmds_bray <- ordinate(filtered_biom_file_2, method="NMDS", distance = "bray")
#plot_ordination(filtered_biom_file_2, nmds_bray, color = "diet_type") + geom_point()

#must plot by hand 
nmds_df <- data.frame(nmds_bray$points) #grabs mds1 values
nmds_df <- cbind(nmds_df, metadata_df)
nmds_df

#plot nmds in ggplot2
ggplot(data = nmds_df, aes(x = MDS1, y = MDS2, color = diet_type)) +
  geom_point(size = 3) +
  ggtitle(label = "NMDS plot of Bray Curtis dissimilarity", paste0("P value: ", p_value_bray)) +
  theme_minimal() +
  labs(color = "Diet Type")

#Plot Jaccard's distance
pcoa_jaccard <- ordinate(filtered_biom_file_2, method = "PCoA", distance = "jaccard", binary = TRUE)
plot_ordination(filtered_biom_file_2, pcoa_jaccard, color = "diet_type") 

pca_jaccard <- round(pcoa_jaccard$values$Relative_eig *100, 1)

pcoa_jaccard_df <- data.frame(pcoa_jaccard$vectors)
pcoa_jaccard_df <- cbind(pcoa_jaccard_df, metadata_df)

#Calculate permanova of jaccard
permanova_jaccard <- adonis2(phyloseq::distance(filtered_biom_file_2, method = "jaccard") ~ diet_type,
                          data = metadata_df, binary = TRUE)  
#grabs p value
p_value_jaccard <- permanova_jaccard$`Pr(>F)`[1]
p_value_jaccard

#plot jacccard distance 
plot_jaccard <- ggplot(data = pcoa_jaccard_df, aes(x = Axis.1, y = Axis.2, color = diet_type)) +
  geom_point(size = 3) +
  xlab(paste0("[",pca_jaccard[1], " %]")) +
  ylab(paste0("[", pca_jaccard[2], " %]")) +
  ggtitle(label = "PCoA plot of Jaccard dissimilarity", paste0("P value: ", p_value_jaccard)) +
  theme_minimal() +
  labs(color = "Diet Type")
#not significant in diversity p value of 0.168

#combine graphs together
plot_combined <- grid.arrange(plot_bray, plot_jaccard, ncol = 2)

########################################
# Differential abundant species
########################################
#Run ancombc using benjamini-hochberg correction
ancombc_out <- ancombc2(data = filtered_biom_file_2, 
                        tax_level = NULL, #want to keep the taxon ID number!
                        fix_formula = "diet_type", 
                        rand_formula = NULL,
                        p_adj_method = "BH", 
                        pseudo_sens = TRUE,
                        prv_cut = 0, 
                        lib_cut = 1000, 
                        s0_perc = 0.05,
                        group = "diet_type", 
                        struc_zero = TRUE,
                        neg_lb = TRUE)

View(ancombc_out$res)

#pull out results table
ancombc_res <- ancombc_out$res
rownames(ancombc_res) <- ancombc_res$taxon 

#clean genus and species names
taxa_names <- as.data.frame(filtered_biom_file_2@tax_table)
taxa_names$taxon <- rownames(tax_df)
taxa_names_clean <- taxa_names %>%
  mutate(Genus = gsub("g__", "", Rank6)) %>%
  mutate(Species = gsub("s__", "", Rank7)) %>%
  unite(col= "Genus_species", Genus, Species, sep = " ")

head(taxa_names_clean)

#merge by rownames
ancombc_res_taxon <- merge(ancombc_res, taxa_names_clean)
head(ancombc_res_taxon)

#Filter by < 0.05
ancombc_res_sig <- ancombc_res_taxon %>% 
  filter(q_diet_typeVegan < 0.05)

ancombc_res_sig
#Unfortunately, no results are significant from adjusted p values or values less than 0.05

#Therefore, will sort by the lowest q values available and pull the top "most significant", as in, which taxa has the smallest p values?
ancombc_res_sig <- ancombc_res_taxon %>%
  arrange(q_diet_typeVegan)

top_20_ancombc_res_sig <- head(ancombc_res_sig, n =20)
View(top_20_ancombc_res_sig)
#From here most "significant" are Gordonibacter massiliensis, Prevotellamassilia uncultured Prevotellamassilia sp

#plot top 20 different abundant species from vegan and omnivores
ggplot(top_20_ancombc_res_sig, aes(x = lfc_diet_typeVegan, y = reorder(Genus_species, lfc_diet_typeVegan))) +
  geom_point(aes(color = q_diet_typeVegan), size = 3) +
  geom_errorbar(aes(xmin = lfc_diet_typeVegan - se_diet_typeVegan, 
                    xmax = lfc_diet_typeVegan + se_diet_typeVegan)) +
  geom_vline(xintercept = 0, color = "red") +
  scale_colour_viridis_c(option = "viridis", name = "Adjusted p values") +
 # scale_color_gradient(low = "red", high = "blue", name = "p adj values") +
  labs(x = "Log Fold Change (diet type Vegan)", 
       y = "Species") +
theme_minimal()
