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
library(rhdf5)  

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

# Why are there chordata and mollusca? that does not seem right??
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

#Calculate relative abundance
physeq_rel_2 <- transform_sample_counts(filtered_biom_file_2, function(x) x / sum(x))
physeq_phy_2 <- tax_glom(physeq_rel_2, taxrank = "Rank2") 
df_2 <-psmelt(physeq_phy_2)

#Plot filtered by phylum
ggplot(df_2, aes(x = Sample, y = Abundance, fill = Rank2)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(y = "Relative Abundance", x = "Samples") +
  scale_x_discrete(guide = guide_axis(angle = 62.5)) +
  theme_minimal() +
  guides(fill = guide_legend(title = "Phylum Level"))

#Plot by class? graph is too busy
#ggplot(df_2, aes(x = Sample, y = Abundance, fill = Rank3)) +
 # geom_bar(stat = "identity", position = "stack") +
  #labs(y = "Relative Abundance", x = "Sample") +
  #scale_x_discrete(guide = guide_axis(angle = 62.5)) +
  #theme_minimal() 

#Create metadata table
metadata_df <- data.frame(diet_type = c("Vegan", "Omnivore", "Omnivore", "Omnivore", "Vegan", "Vegan",
                                        "Omnivore", "Omnivore", "Vegan", "Vegan"))
#ADD SRR # as row names
row.names(metadata_df) <- SRR_names

#Add metadata to phyloseq file
sample_data(filtered_biom_file_2) <- metadata_df

#plot alpha diversity
plot <- plot_richness(filtered_biom_file_2, color = "diet_type", measures = c("Observed", "Chao1", "Shannon", "Simpson")) +
  labs(color = "Diet Type")
plot
estimates <- estimate_richness(filtered_biom_file_2)
estimates
#beta diversity 
pcoa_bray <- ordinate(filtered_biom_file_2, method = "PCoA", distance = "bray")
#plot_ordination(filtered_biom_file_2, pcoa_bray, color = "diet_type") + 
 # geom_point(aes(color = "diet_type"), size = 4, inherit.aes = T)

#why is colour not working waaaaa
sample_data(filtered_biom_file_2)

#NEED TO PLOT BY HAND! :(
#get the top variances 
pca <- round((pcoa_bray[["values"]]$Relative_eig),3) *100

#grabs the eigenvalues for each sample
pcoa_bray_df <- data.frame(pca = pcoa_bray$vectors)
pcoa_bray_df <- cbind(pcoa_bray_df, metadata_df) #combine metadata with dataframe
  
#plot again 
ggplot(data = pcoa_bray_df, aes(x = pca.Axis.1, y = pca.Axis.2, color = diet_type)) +
  geom_point(size = 3) +
  xlab(paste0("[",pca[1], " %]")) +
  ylab(paste0("[", pca[2], " %]")) +
  ggtitle(label = "PCoA plot of Bray Curtis dissimilarity") +
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
  ggtitle(label = "NMDS plot of Bray Curtis dissimilarity") +
  theme_minimal() +
  labs(color = "Diet Type")

#Plot Jaccard's distance
pcoa_jaccard <- ordinate(filtered_biom_file_2, method = "PCoA", distance = "jaccard")
plot_ordination(filtered_biom_file_2, pcoa_jaccard, color = "diet_type") 

pca_jaccard <- round(pcoa_jaccard$values$Relative_eig *100, 1)

pcoa_jaccard_df <- data.frame(pcoa_jaccard$vectors)
pcoa_jaccard_df <- cbind(pcoa_jaccard_df, metadata_df)

#plot jacccard distance 
ggplot(data = pcoa_jaccard_df, aes(x = Axis.1, y = Axis.2, color = diet_type)) +
  geom_point(size = 3) +
  xlab(paste0("[",pca_jaccard[1], " %]")) +
  ylab(paste0("[", pca_jaccard[2], " %]")) +
  ggtitle(label = "PCoA plot of Jaccard dissimilarity") +
  theme_minimal() +
  labs(color = "Diet Type")

#perform permanova
permanoa <- adonis2(phyloseq::distance(filtered_biom_file_2, method = "bray") ~ diet_type,
        data = metadata_df)
#not significant in diversity p value of 0.14

#Species level
#ancombc_out <- ancombc2(data = filtered_biom_file_2, tax_level = "Rank7",
 #                       fix_formula = "diet_type", rand_formula = NULL,
  #                      p_adj_method = "holm", pseudo_sens = TRUE,
   #                     prv_cut = 0, lib_cut = 1000, s0_perc = 0.05,
    #                    group = "diet_type", struc_zero = TRUE, neg_lb = TRUE)

#Run ancombc using benjamini-hochberg correction
ancombc_out_bon <- ancombc2(data = filtered_biom_file_2, tax_level = "Rank7",
                        fix_formula = "diet_type", rand_formula = NULL,
                        p_adj_method = "BH", pseudo_sens = TRUE,
                        prv_cut = 0, lib_cut = 1000, s0_perc = 0.05,
                        group = "diet_type", struc_zero = TRUE, neg_lb = TRUE)

#Plot results
ggplot(ancombc_out$res, aes(x = lfc_diet_typeVegan, y = reorder(taxon, lfc_diet_typeVegan))) +
  geom_point(aes(color = p_diet_typeVegan < 0.05), size = 3) +
  geom_errorbar(aes(xmin = lfc_diet_typeVegan - se_diet_typeVegan, 
                    xmax = lfc_diet_typeVegan + se_diet_typeVegan)) +
  geom_vline(xintercept = 0, color = "red") +
  labs(x = "Log Fold Change (diet_typeVegan)", 
       y = "Genus")

#THERE are too many species, will probably need a more stringent definition
#pull out results table
ancombc_res <- ancombc_out_bon$res

#Filter by < 0.05
ancombc_res_sig <- ancombc_res %>% 
  filter(q_diet_typeVegan < 0.05)

ancombc_res_sig
#Unfortunately, no results are significant from adjusted p values or values less than 0.05

#Therefore, will sort by the lowest q values available and pull the top "most significant", as in, which have the smallest p values?
ancombc_res_sig <- ancombc_res %>%
  arrange(q_diet_typeVegan)

top_20_ancombc_res_sig <- head(ancombc_res_sig, n =20)
#From here most "significant" are s__massiliensis, Prevotella heparinolytica, and s__uncultured Prevotellamassilia

clean_taxon <- gsub("^k.*?g__", "",top_20_ancombc_res_sig$taxon)

top_20_ancombc_res_sig <- cbind(top_20_ancombc_res_sig, clean_taxon)

#plot top 20 different abundant species from vegan and omnivores
ggplot(top_20_ancombc_res_sig, aes(x = lfc_diet_typeVegan, y = reorder(clean_taxon, lfc_diet_typeVegan))) +
  geom_point(aes(color = q_diet_typeVegan), size = 3) +
  geom_errorbar(aes(xmin = lfc_diet_typeVegan - se_diet_typeVegan, 
                    xmax = lfc_diet_typeVegan + se_diet_typeVegan)) +
  geom_vline(xintercept = 0, color = "red") +
  scale_colour_viridis_c(option = "viridis", name = "Adjusted p values") +
 # scale_color_gradient(low = "red", high = "blue", name = "p adj values") +
  labs(x = "Log Fold Change (diet type Vegan)", 
       y = "Species") +
  theme_minimal()


ggplot(x, aes(x = lfc_diet_typeVegan, y = reorder(taxon, lfc_diet_typeVegan))) +
  geom_point(aes(color = p_diet_typeVegan), size = 3) +
  geom_errorbar(aes(xmin = lfc_diet_typeVegan - se_diet_typeVegan, 
                    xmax = lfc_diet_typeVegan + se_diet_typeVegan)) +
  geom_vline(xintercept = 0, color = "red") +
  labs(x = "Log Fold Change (diet_typeVegan)", 
       y = "Genus")
