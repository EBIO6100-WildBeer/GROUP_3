# Wild Beer Metagenomics- Data Clean Up and Exploratory Data Analysis
# Begun October 24, 2022
########################################

### FIRST, SET TO YOUR OWN WORKING DIRECTORY DIRECTORY WHERE project_data_and_tutorials IS ######
########## SET UP ##########

# Read in files and check them out (from github repo on local machine)
beerComm <- read.csv(file="project_data_and_tutorials/data/2022_10_20_WildBeer.phyloFlash.extractedSSUclassifications.csv") #ASV table (sort of)
beerMeta <- read.delim(file="project_data_and_tutorials/data/2022_10_12_WildBeer_metadata.tsv", sep="\t", header=TRUE) #metadata


head(beerComm)
head(beerMeta)
# View(beerComm)
# View(beerMeta)

# Libraries
library(vegan) #version 2.5-7
library(phyloseq)
library(tidyverse)
library(dendextend)

########## MAIN PART OF SCRIPT ##########

########################################
# MAKE PHYLOSEQ OBJECTS
########################################

# 1. CLEAN UP TAXONOMY
# Right now, the fungal and bacteria taxonomy is not ordered the same. For example, the fungal taxonomic information
# starts at Eukaryota, whereas the bacterial starts at bacteria (not prokaryote). And there is a whole lot of 
# extra taxonomic information between "eukaryote" and "fungi" in the fungal ones.

## i. first, take a look at the raw data
colnames(beerComm)
head(beerComm)
step1 <- beerComm[,c(5,6)] #get only dbhit and taxonomy columns; keep all rows
colnames(step1)

## ii. Remove extra taxonomic information between Eukaryote and Fungi, as well as subkingdom and subdivisions
# btw, here is a great "cheat sheet" for using regular expressions in R: https://evoldyn.gitlab.io/evomics-2018/ref-sheets/R_strings.pdf
step2 <- step1[,2] %>%  
  stringr::str_replace("^.*Nucletmycea;", "") %>% # "^" is start of string, "." is every character, and "*" is any length.
  # So, this line means replace any and everything between the start of the line and "Nucletmycea" with nothing. 
  stringr::str_replace("Dikarya;", "") %>% #replace "Dikarya" with "" (i.e. replace it with nothing!)
  stringr::str_replace("Saccharomycotina;", "") %>% 
  stringr::str_replace("Agaricomycotina;", "")
colnames(step2)
cleaned1 <- cbind(step1[,1], step2) #combine the cleaned taxonomy column with dbhit from step 1 
colnames(cleaned1) <- colnames(step1) #make dbHit and taxonomy the column names again
# View(cleaned1)

## iii. Divide out taxonomy into different columns by level
beerTax <- tidyr::separate(as.data.frame(cleaned1), col = taxonomy, into= c("Kingdom", "Phylum", "Class", 
                                                                  "Order", "Family", "Genus", "Species"),
                    sep = ";")
# View(beerTax) 
colnames(beerComm)
unique(beerTax$dbHit == beerComm$dbHit) # since the order and the names of the dbHits are the same for beerTax and beerCommLess, 
# we can just do a simple cbind to merge them
beerCommLess <- beerComm[,c(1:4,7:10)] #no need to keep dbhit or Taxonomy columns (since we fixed taxonomy)
dim(beerCommLess)
beerCommCleaned <- cbind(beerCommLess, beerTax)
# View(beerCommCleaned)

# 2. MAKE DATA "WIDER" (AS REQUIRED BY PHYLOSEQ)
commCompTable <- beerCommCleaned %>% #beerCommCleaned is on this line because it is the object we are manipulating
  select(read_cov, sample:Species) %>% #makes it so that we are manipulating the columns read_cov 
  # and all of the columns from sample to Species (i.e. sample, dbHit, and all the taxonomy info)
  filter(!(sample %in% c('NTC', 'ExtB1', 'CommStd'))) %>%  #remove these samples. Important because they are not in our metadata,
  # and to make a phyloseq object, samples must be the same across metadata and ASV table
  pivot_wider(names_from=sample, values_from=read_cov) %>% #make the samples new columns and the values in these new columns 
  # come from read_cov
  mutate_all(~replace(., is.na(.), 0)) #replace NAs with zeros
# View(commCompTable)

# 3. MAKE TAXONOMY TABLE FOR PHYLOSEQ
# ASV/OTU names (i.e. dbHit for us!) should be the rownames and then their should be a distinct column for each taxonomic level
taxTab <- commCompTable %>% 
  column_to_rownames("dbHit") %>% #make it so that dbHit is now the row names of this dataframe
  select(Kingdom:Species) #take only the columns corresponding to the taxonomy info from commCompTable
# View(taxTab) #looks good!

# 4. MAKE ASV/OTU TABLE FOR PHYLOSEQ
# Here, phyloseq wants each taxon/ASV/OTU (here our dbHit name) as the rowname, with the columns as the separate
# samples. The values in these columns are simply the number of counts of each ASV/OTU/taxon in each sample.
# This is now the rest of the info in "commCompTable" made above
ASVtab <- commCompTable %>% 
  column_to_rownames("dbHit") %>% #make it so that dbHit is now the row names of this dataframe
  select("HI4-3":"HI4-1")  #R is being tripped up by the fact that samples have "-3" and "-1" in their names, so I
# added quotation marks. This is usually not necessary when using the select function.
# View(ASVtab) #Yay! Looks as expected!

# 5. FORMAT METADATA FILE FOR PHYLOSEQ
# View(beerMeta)
metaDat <- beerMeta %>% column_to_rownames("Sample.ID") #make the rownames the sample IDs. These now match up with 
# ASVtab
# View(metaDat)

# 6. PUT TAXONOMY, ASV TABLE, AND METADATA TOGETHER TO MAKE PHYLOSEQ OBJECT
# Need to use built-in phyloseq functions (tax_table, otu_table, and sample_data) to build each component of the phyloseq object
TAX <- tax_table(as.matrix(taxTab))
OTU <- otu_table(ASVtab, taxa_are_rows=T) 
META <- sample_data(metaDat)
beerPhyloseq <- phyloseq(TAX, OTU, META)

# This is the original phyloseq object! To see any of the component parts, re-use the functions above
beerPhyloseq 
sample_data(beerPhyloseq) #wahoo! Metadata!
tax_table(beerPhyloseq) 
otu_table(beerPhyloseq)

########################################
# DATA TRANSFORMATIONS AND ORDINATION VISUALIZATIONS
########################################
 
# 1. Perform a Hellinger transformation on the data. This converts the data to proportions and 
# then takes the square root.  
beerHellinger_ps <- transform_sample_counts(beerPhyloseq, function(x) sqrt(x / sum(x)))
otu_table(beerHellinger_ps) #looks about as expected

# 2. Get dissimilarity matrices
# Bray-Curtis (abundance based)
beerBCdist <- distance(beerHellinger_ps, method= "bray")

# Binary Jaccard  (presence/absence)
beerJaccardPA_dist <- distance(beerHellinger_ps, method= "jaccard", binary= TRUE) #binary must be set to true, or this is quantitative Jaccard

# 3. Make (unconstrained) ordinations for Bray-Curtis
# PCoA
ordBC_pcoa <- ordinate(beerHellinger_ps, "PCoA", "bray")
BCordplot <- plot_ordination(beerHellinger_ps, ordBC_pcoa, type= "samples", color= "Hops", shape= "Location") +
  geom_point(size=5) +
  ggtitle("PCoA based On Bray-Curtis Dissimilarities")
BCordplot

# NMDS
ordBC_nmds <- ordinate(beerHellinger_ps, "NMDS", "bray")
ordBCNMDS_plot <- plot_ordination(beerHellinger_ps, ordBC_nmds, type= "samples", color= "Hops", shape= "Location") +
  geom_point(size=5) +
  ggtitle("NMDS based On Bray-Curtis Dissimilarities")
ordBCNMDS_plot

# 4 .Make (unconstrained) ordinations for binary Jaccard
# PCoA
ordPA_pcoa <- ordinate(beerHellinger_ps, "PCoA", "jaccard", binary= TRUE)
PAordplot <- plot_ordination(beerHellinger_ps, ordPA_pcoa, type= "samples", color= "Hops", shape= "Location") +
  geom_point(size=5) +
  ggtitle("PCoA based on Binary Jaccard Dissimilarities")
PAordplot

# NMDS
ordPA_nmds <- ordinate(beerHellinger_ps, "NMDS", "jaccard", binary= TRUE)
ordBCNMDS_plot <- plot_ordination(beerHellinger_ps, ordPA_nmds, type= "samples", color= "Hops", shape= "Location") +
  geom_point(size=5) +
  ggtitle("NMDS based on Binary Jaccard Dissimilarities")
ordBCNMDS_plot

########################################
# STACKED BARPLOT FOR RELATIVE ABUNDANCE
########################################
ps_rel_abund <- phyloseq::transform_sample_counts(beerHellinger_ps, function(x){x / sum(x)})

plot_bar(ps_rel_abund, fill = "Genus") +
  geom_bar(aes(color = Genus, fill = Genus), stat = "identity", position = "stack") +
  labs(x = "", y = "Relative Abundance\n") +
  facet_wrap(~ Hops, scales = "free") +
  theme(panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

# HEATMAP ---- 
plot_heatmap(ps_rel_abund, method = "PCoA", distance = "bray", 
             taxa.label = "Genus", taxa.order = "Genus", 
             low="beige", high="red", na.value="beige")

# CLUSTER ---- 

#Extract OTU table and compute BC
ps_rel_otu <- data.frame(phyloseq::otu_table(ps_rel_abund))
ps_rel_otu <- t(ps_rel_otu)
bc_dist <- vegan::vegdist(ps_rel_otu, method = "bray")
as.matrix(bc_dist)

#Save as dendrogram
ward <- as.dendrogram(hclust(bc_dist, method = "ward.D2"))

#Provide color codes (looking at hopped vs un-hopped)
meta <- data.frame(phyloseq::sample_data(ps_rel_abund))
colorCode <- c(H = "red", `U` = "blue")
labels_colors(ward) <- colorCode[meta$Hops][order.dendrogram(ward)]

#Plot
plot(ward)


###########################################
# Filter out controls and wonky samples
########################################
# Remove "control" location. This is now a phyloseq object with everything but the controls. For some
# reason, I couldn't filter out on location and "X" at the same time, so I did it in two steps
step1_ps <- subset_samples(beerPhyloseq, Location != "C") 

# Remove wonky samples (control are already removed). 
beerPS_cleaned <- subset_samples(step1_ps, X != "24") # I removed sample 24, i.e. the hopped indoors 5, week 1

# Make 2 phyloseq objects, one with only the hopped samples and the other with only unhopped ones
beerHopped_ps <- subset_samples(beerPS_cleaned, Hops == "H") #get only hopped samples

beerUnhopped_ps <- subset_samples(beerPS_cleaned, Hops == "U") #get only unhopped samples

# 






