# Code to make the phyloseq object in R rather than in excel 

# INIT ----

rm(list=ls()) 

# loading pacakges 
library(tidyverse)
library(vegan) 
library(phyloseq)
library(dendextend)


# loading data 
beerComm <- read.csv("~/Documents/PhD projects/beer metagenomics/2022_10_20_WildBeer.phyloFlash.extractedSSUclassifications.csv")
beerMeta <- read.csv("~/Documents/PhD projects/beer metagenomics/metadata.csv")


# CLEAN UP TAXONOMY ----
# Right now, the fungal and bacteria taxonomy is not ordered the same. For example, the fungal taxonomic information
# starts at Eukaryota, whereas the bacterial starts at bacteria (not prokaryote). And there is a whole lot of 
# extra taxonomic information between "eukaryote" and "fungi" in the fungal ones.

## first, take a look at the raw data
colnames(beerComm)
head(beerComm)
step1 <- beerComm[,c(5,6)] #get only dbhit and taxonomy columns; keep all rows
colnames(step1)

## ii. Remove extra taxonomic information between Eukaryote and Fungi, as well as subkingdom and subdivisions
step2 <- step1[,2] %>%  
    stringr::str_replace("^.*Nucletmycea;", "") %>% 
    stringr::str_replace("Dikarya;", "") %>% 
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

#MAKE DATA "WIDER" ----
# (AS REQUIRED BY PHYLOSEQ)
commCompTable <- beerCommCleaned %>% #beerCommCleaned is on this line because it is the object we are manipulating
    select(read_cov, sample:Species) %>% #makes it so that we are manipulating the columns read_cov 
    # and all of the columns from sample to Species (i.e. sample, dbHit, and all the taxonomy info)
    filter(!(sample %in% c('NTC', 'ExtB1', 'CommStd'))) %>%  #remove these samples. Important because they are not in our metadata,
    # and to make a phyloseq object, samples must be the same across metadata and ASV table
    pivot_wider(names_from=sample, values_from=read_cov) %>% #make the samples new columns and the values in these new columns 
    # come from read_cov
    mutate_all(~replace(., is.na(.), 0)) #replace NAs with zeros
# View(commCompTable)

# MAKE TAXONOMY TABLE FOR PHYLOSEQ ----
# ASV/OTU names (i.e. dbHit for us!) should be the rownames and then their should be a distinct column for each taxonomic level
taxTab <- commCompTable %>% 
    column_to_rownames("dbHit") %>% #make it so that dbHit is now the row names of this dataframe
    select(Kingdom:Species) #take only the columns corresponding to the taxonomy info from commCompTable
# View(taxTab) #looks good!

# MAKE ASV/OTU TABLE FOR PHYLOSEQ ----
# Here, phyloseq wants each taxon/ASV/OTU (here our dbHit name) as the rowname, with the columns as the separate
# samples. The values in these columns are simply the number of counts of each ASV/OTU/taxon in each sample.
# This is now the rest of the info in "commCompTable" made above
ASVtab <- commCompTable %>% 
    column_to_rownames("dbHit") %>% #make it so that dbHit is now the row names of this dataframe
    select("HI4-3":"HI4-1")  #R is being tripped up by the fact that samples have "-3" and "-1" in their names, so I
# added quotation marks. This is usually not necessary when using the select function.
# View(ASVtab) #Yay! Looks as expected!

# FORMAT METADATA FILE FOR PHYLOSEQ ----
# View(beerMeta)
metaDat <- beerMeta %>% column_to_rownames("Sample.ID") #make the rownames the sample IDs. These now match up with 
# View(metaDat)

# MAKE PHYLOSEQ OBJECT ----
# Need to use built-in phyloseq functions (tax_table, otu_table, and sample_data) to build each component of the phyloseq object
TAX <- tax_table(as.matrix(taxTab))
OTU <- otu_table(ASVtab, taxa_are_rows=T) 
META <- sample_data(metaDat)
beerPhyloseq <- phyloseq(TAX, OTU, META)

# This is the phyloseq object! To see any of the component parts, re-use the functions above
beerPhyloseq 
sample_data(beerPhyloseq) #wahoo! Metadata!
tax_table(beerPhyloseq) 
otu_table(beerPhyloseq)



# HELLINGER TRANSFORMATION ----
hell <- transform_sample_counts(beerPhyloseq, function(x) sqrt(x / sum(x)))

# 100% RELATIVE ABUNDANCE CHART ----
#Convert to relative abundance
ps_rel_abund <-  phyloseq::transform_sample_counts(hell, function(x){x / sum(x)})
phyloseq::otu_table(beerPhyloseq)
class(ps_rel_abund)

# plot the chart
# Genus grouped by timepoint
phyloseq::plot_bar(ps_rel_abund, fill = "Genus") +
    geom_bar(aes(color = Genus, fill = Genus), stat = "identity", position = "stack") +
    labs(x = "", y = "Relative Abundance\n") +
    facet_wrap(~ timepoint, scales = "free") +
    theme(panel.background = element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())


# Genus grouped by hops/no hops
phyloseq::plot_bar(ps_rel_abund, fill = "Genus") +
    geom_bar(aes(color = Genus, fill = Genus), stat = "identity", position = "stack") +
    labs(x = "", y = "Relative Abundance\n") +
    facet_wrap(~ Hops, scales = "free") +
    theme(panel.background = element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())


# HEATMAPS ---- 
plot_heatmap(ps_rel_abund, method = "PCoA", distance = "bray", 
             taxa.label = "Genus", taxa.order = "Genus", 
             low="beige", high="red", na.value="beige")


# BRAY-CURTIS DISTANCE  ----
# with non-transformed data
dist <- distance(beerPhyloseq, method="bray", type="samples")
ord <- ordinate(beerPhyloseq, "PCoA", distance="bray")
p <- plot_ordination(beerPhyloseq, ord, type = "samples")
p

# with transformed data
dist <- distance(ps_rel_abund, method="bray", type="samples")
ord <- ordinate(ps_rel_abund, "PCoA", distance="bray")
plot <- plot_ordination(ps_rel_abund, ord, type = "samples", color = "Hops") +
    geom_point(size = 5) +  theme_bw()
plot

# Add binary Jaccard for presence absence 

# PCA  ----
# Literally no idea what is happening in following steps... but there is a figure? 

#PCA via phyloseq
ord_2 <- ordinate(beerPhyloseq, "RDA", distance="bray")

#Plot scree plot
phyloseq::plot_scree(ord_2) + 
    geom_bar(stat = "identity", fill = "blue") +
    labs(x = "\nAxis", y = "Proportion of Variance\n")

#Scale axes and plot ordination
in1 <- ord_2$CA$eig[1] / sum(ord_2$CA$eig)
in2 <- ord_2$CA$eig[2] / sum(ord_2$CA$eig)

phyloseq::plot_ordination(beerPhyloseq, ord_2) + 
    geom_point(size = 2) +
    coord_fixed(in2 / in1) +
    stat_ellipse(aes(group = timepoint), linetype = 2)

# PCoA  ----
# Warning messages:
# 1: In plot_ordination(beerPhyloseq, ordu, color = "Phylum", shape = "taxa") :
#    Color variable was not found in the available data you provided.No color mapped.
# 2: In plot_ordination(beerPhyloseq, ordu, color = "Phylum", shape = "taxa") :
#    Shape variable was not found in the available data you provided.No shape mapped.

ordu = ordinate(beerPhyloseq, "PCoA", "bray", weighted=TRUE)
plot_ordination(beerPhyloseq, ordu, color="Class", shape="taxa")





# CLUSTER ANALYSIS ----

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




