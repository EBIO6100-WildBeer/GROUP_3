# Wild Beer Metagenomics- Group 3 community analysis
# Begun October 24, 2022

########################################
rm(list=ls())

### FIRST, SET TO YOUR OWN WORKING DIRECTORY DIRECTORY WHERE project_data_and_tutorials is ######
# for claire this should be setwd("~/Desktop/CU_Research/wildBeer")
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
library(indicspecies)

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
which(beerCommCleaned$dbHit== "KR952335.21166.22627") #find the row that has the mitochondrion in it, row 96
beerCommCleaned <- beerCommCleaned[-96,] #remove this mitochondrion row
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
# DATA TRANSFORMATIONS, SUBSETTING, AND DISSIMILARITIES
########################################
# 1. Perform a Hellinger transformation on the data. This converts the data to proportions and 
# then takes the square root.  
beerHellinger_ps <- transform_sample_counts(beerPhyloseq, function(x) sqrt(x / sum(x)))
otu_table(beerHellinger_ps) #looks about as expected

# 2. Remove "control" location. This is now a phyloseq object with everything but the controls. For some
# reason, I couldn't filter out on location and "X" at the same time, so I did it in two steps
step1_ps <- subset_samples(beerHellinger_ps, Location != "C") 

# Remove wonky samples (control are already removed). 
beerPS_cleaned <- subset_samples(step1_ps, Sample.Description != "Unhopped, Indoor, Jar 5, wk1, filter 1") #sample has no Saccharomyces and low # of reads

# Make 2 phyloseq objects, one with only the hopped samples and the other with only unhopped ones
beerHopped_ps <- subset_samples(beerPS_cleaned, Hops == "H") #get only hopped samples

beerUnhopped_ps <- subset_samples(beerPS_cleaned, Hops == "U") #get only unhopped samples

beerWeek1_ps <- subset_samples(beerPS_cleaned, timepoint == "1 week") #get only hopped samples

beerWeek3_ps <- subset_samples(beerPS_cleaned, timepoint == "3 week") #get only unhopped samples

# 2. Get dissimilarity matrices
# i. All samples
# Bray-Curtis (abundance based)
beerBCdist <- distance(beerPS_cleaned, method= "bray")
# Quantitative Jaccard
beerQJdist <- distance(beerPS_cleaned, method= "jaccard", binary = FALSE)
# Binary Jaccard  (presence/absence)
beerJaccardPA_dist <- distance(beerPS_cleaned, method= "jaccard", binary= TRUE) #binary must be set to true, or this is quantitative Jaccard
# ii. hops only
# Bray-Curtis
hoppedBCdist <- distance(beerHopped_ps, method= "bray")
# quantitative Jaccard
hoppedQJdist <- distance(beerHopped_ps, method= "jaccard", binary = FALSE)
# iii. unhopped only
unhoppedBCdist <- distance(beerUnhopped_ps, method= "bray")
# quantitative Jaccard
unhoppedQJdist <- distance(beerUnhopped_ps, method= "jaccard", binary = FALSE)
# iv. week 1 only
week1_BCdist <- distance(beerWeek1_ps, method= "bray")
# quantitative Jaccard
week1_QJdist <- distance(beerWeek1_ps, method= "jaccard", binary = FALSE)
# v. week 3 only
week3_BCdist <- distance(beerWeek3_ps, method= "bray")
# quantitative Jaccard
week3_QJdist <- distance(beerWeek3_ps, method= "jaccard", binary = FALSE)


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
ps_rel_abund <- phyloseq::transform_sample_counts(beerPS_cleaned, function(x){x / sum(x)})

plot_bar(ps_rel_abund, fill = "Genus") +
  geom_bar(aes(color = Genus, fill = Genus), stat = "identity", position = "stack") +
  labs(x = "", y = "Relative Abundance\n") +
  facet_wrap(~ Hops, scales = "free") +
  theme(panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

# HEATMAP ---- 
plot_heatmap(beerPS_cleaned, method = "PCoA", distance = "bray", 
             taxa.label = "Genus", taxa.order = "Genus", 
             low="beige", high="red", na.value="beige")

# Now, we want to make a heat map that is grouped by 


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


# Alpha Diversity (currently using shannon)
# Adding a col which is the shannon value 
sample_data(beerHellinger_ps)$shannon.physeq <- estimate_richness(beerHellinger_ps, measures = "Shannon")

# Plot 1 (all timepoints, sep by hops) 
plot_richness(beerPS_cleaned, "Hops", color = "treatment", measures = "Shannon") 
# Plot 2 (week 1 only, sep by hops)
plot_richness(beerWeek1_ps, "Hops", color = "treatment", measures = "Shannon")  
# Plot 3 (week 3 only, sep by hops)
plot_richness(beerWeek3_ps, "Hops", color = "treatment", measures = "Shannon")  
# Plot 4 (all timepoints, sep by location)
plot_richness(beerPS_cleaned, "Location", color = "treatment", measures = "Shannon") 
# Plot 5 (week 1 only, sep by location)
plot_richness(beerWeek1_ps, "Location", color = "treatment", measures = "Shannon")  
# Plot 6 (week 3 only, sep by hops)
plot_richness(beerWeek3_ps, "Location", color = "treatment", measures = "Shannon") 
# Plot 7 (all timepoints, sep by timepoints)
plot_richness(beerPS_cleaned, "timepoint", color = "treatment", measures = "Shannon") 
# Plot 8 (all timepoints, sep by treatment)
plot_richness(beerPS_cleaned, "treatment", color = "treatment", measures = "Shannon") + geom_boxplot() 
# Plot 9 (week 1 only, sep by treatment)
plot_richness(beerWeek1_ps, "treatment", color = "treatment", measures = "Shannon") + geom_boxplot() 
# Plot 10 (week 3 only, sep by treatment)
plot_richness(beerWeek3_ps, "treatment", color = "treatment", measures = "Shannon") + geom_boxplot() 



########################################
# HYPOTHESIS TESTING
########################################

# 1. QUESTION 1: Is hopped different than unhopped?
# Method 1: PERMANOVA-- this is a non-parametric method that tests the null hypothesis that the centroids and
# and dispersions of the groups as defined by measure space are equivalent for all groups. If the null is 
# rejected, then it means that EITHER the centroid or the spread of the objects is different between the groups.
# The PERMANOVA is done on the underlying distance matrix (i.e. NOT on the output of the ordination technique)
# Thus, Thus, if H0 were true, any observed differences among the centroids in a given set of data will be similar
# in size to what would be obtained under random allocation of individual sample units to the groups (i.e., under permuration)

# Look within whole dataset
# i. PERMANOVA: 
# set seed!!
set.seed(93)
adonis2(formula = beerBCdist ~ Hops, data= as.data.frame(as.matrix(sample_data(beerPS_cleaned))))
# This shows that hopped samples are different than unhopped samples (p < 0.001)

# ii. Mantel test. Although Mantel tests are usually 

# QUESTION 1.1: We know that hopped samples are significantly different than unhopped samples. Do hopped and unhopped 
# samples differ in how similar the samples within each group are to one another?
# Calculate multivariate dispersions
mod_disp1 <- betadisper(beerBCdist, group= factor(as.data.frame(as.matrix(sample_data(beerPS_cleaned)))$Hops), 
                        type = "centroid")
# test it
# set seed!!
set.seed(93)
anova(mod_disp1, permutations = 9999)
# This shows that there the differences among samples are the same hopped and unhopped, which is surprising

# Plot it!
plot(mod_disp1)
# 2. QUESTION 2: Is indoors different from outdoors?
# Since we know that hopped samples are different than unhopped samples, we will test this on hopped and unhopped samples
# separately.

#. i. PERMANOVA ON HOPPED:
# set seed!!
set.seed(93)
adonis2(formula = hoppedBCdist ~ Location, data= as.data.frame(as.matrix(sample_data(beerHopped_ps))))
# Hopped samples are not different indoors or outdoors, but it's somewhat close to significant

# ii. Is a difference in variation among hopped samples indoors or outdoors?
mod_disp2 <- betadisper(hoppedBCdist, group= factor(as.data.frame(as.matrix(sample_data(beerHopped_ps)))$Location), 
                        type = "centroid")
# test it
# set seed!!
set.seed(93)
anova(mod_disp2, permutations = 9999) #no, no difference (P = 0.93323)

#. iii. PERMANOVA ON UNHOPPED:
# set seed!!
set.seed(93)
adonis2(formula = unhoppedBCdist ~ Location, data= as.data.frame(as.matrix(sample_data(beerUnhopped_ps))))
# Unhopped samples are also not significantly different based on being indoors or outdoors

# iv. Is a difference in variation among unhopped samples indoors or outdoors?
mod_disp3 <- betadisper(unhoppedBCdist, group= factor(as.data.frame(as.matrix(sample_data(beerUnhopped_ps)))$Location), 
                        type = "centroid")
# test it
# set seed!!
set.seed(93)
anova(mod_disp3, permutations = 9999) #YES!!! P = 0.0217
# So this means that within unhopped samples, either indoor or outdoor samples are more different from each other
# Looking at the average distances to centroid shows that there is 0.28 average distance to centroid for indoor 
# samples and 0.1683 for outdoor. So, the indoor samples are actually a little more varied.
mod_disp3

plot(mod_disp3)

# PCoA
ordunhopped_pcoa <- ordinate(beerUnhopped_ps, "PCoA", "bray")
unhopped_pcoa_plot <- plot_ordination(beerUnhopped_ps, ordunhopped_pcoa, type= "samples", color= "timepoint", shape= "Location") +
  geom_point(size=5) +
  ggtitle("(Unhopped samples) PCoA based On Bray-Curtis Dissimilarities")
unhopped_pcoa_plot

# QUESTION 3: Timepoints: Is week 1 different from week 3? Can break down even further to say “is there an indoor outdoor difference
# at week 1 versus is there an indoor/outdoor difference during week 3? Which taxa are relatively more abundant in week 1 versus week 3?

# i. This asks, within hopped samples, is there a difference between week 1 and week 3 and does being inside or outside matter
# for this relationship?
# set seed!!
set.seed(93)
adonis2(formula = hoppedBCdist ~ timepoint + Location*timepoint, data= as.data.frame(as.matrix(sample_data(beerHopped_ps))))
# This tells us that that within hopped samples, there is a difference between week 1 and week 3 and between indoor and outdoor
# samples. However, there is no significant interaction effect between location and timepoint 

# ii. This asks, within unhopped samples, is there a difference between week 1 and week 3 and does being inside or outside matter
# for this relationship?
# set seed!!
set.seed(93)
adonis2(formula = unhoppedBCdist ~ timepoint + Location*timepoint, data= as.data.frame(as.matrix(sample_data(beerUnhopped_ps))))
# This shows that within unhopped samples, week 1 is different than week 3, and that it this relationship holds for both
# indoor and outdoor samples. Location does not seem to matter

# QUESTION 4: What taxa are driving the differences?
# For this, I perform an indicator species analysis, using a function (indSpecTable2) that I wrote previously,
# found here: https://github.com/clairecwinfrey/PhanBioMS_scripts/blob/master/R_scripts/analysis/IndicSpeciesAnalysisFInalJan2022.R

# 1. Code for function (found at link above)
# New function: indSpecTable2
##############################
# indSpecTable2 take the output of the indicator species analysis,
# the taxonomy table for the samples, and the ASV table for the 
# samples. It returns a dataframe that has:
# 1) the ASV TaxonID,
# 2) the taxonomic name of the ASV, 3) the biserial correlation coefficient
# ("stat"), 4) the p-value for the indicator ASV, 5) the
# group of samples that the ASV is an indicator for, 6) and the relative
# abundance of each indicator ASV in each group of samples!

indSpecTable2 <- function(myMultipatt, taxTable, ASVtable, signLevel=0.05) { 
  # Arguments are: 
  # 1) myMultipatt: an object produced by multipatt function;
  # 2) taxTable: a taxonomy table that has the taxonIDs as row names and 
  # 7 columns for taxonomic levels (Kingdom, Phylum, Class...Species);
  # 3) signLevel: significance level for p-value cut-off
  # 4) ASV table where samples/sites are rows and TaxonIDs are columns
  # (this will be inverted in script)
  require(indicspecies)
  groups <- unique(myMultipatt$cluster) #gets groups from myMultipatt
  # The lines below pull out only the indicator ASVs that have a
  # p-value less than or equal to the significance level specified in the 
  # function call. Indicator ASVs are represented by TaxonID
  ASVstoKeep <- rownames(myMultipatt$sign)[which(myMultipatt$sign$p.value <= signLevel)]
  ### replace ASVdfs with ASVsToKeep
  ASVtable_t <- t(ASVtable) #make rows TaxonID 
  ASVtableIndex <- which(rownames(ASVtable_t) %in% ASVstoKeep)
  ASV_sigs <- ASVtable_t[ASVtableIndex,] #new ASV table only has significant ASVs
  # Below makes it so that individual samples are labeled according to clustering group
  colnames(ASV_sigs) <- myMultipatt$cluster
  # Combine by clustering groups
  combByCluster <- t(rowsum(t(ASV_sigs), group=colnames(ASV_sigs), na.rm =T)) 
  # Below, pre-allocate a dataframe that we'll fill in with the relative abundances
  # of each ASV in each group
  relAbundDf <- data.frame(matrix(nrow=nrow(combByCluster), ncol=(length(groups)))) 
  # The for loop below makes column names for relAbundDf based on names of groups used for clustering
  relAbundCol <- rep(NA,length(groups)) #pre-allocate vector to hold relative abundance columns
  for (i in 1:length(groups)){ #make column names for the relative abundance of each group in analysis
    relAbundCol[i] <- paste("relAbund in", colnames(combByCluster)[i])
  }
  colnames(relAbundDf) <- relAbundCol
  totalCounts <- colSums(combByCluster) #get total ASV counts in each clustering group
  # For loop calculates the relative abundance of each significant ASV in each grouping:
  for(j in 1:nrow(combByCluster)){  #loop over all the significant ASVs
    relAbundDf[j,] <- (combByCluster[j,]/totalCounts)*100
  }
  rownames(relAbundDf) <- rownames(combByCluster) #make row names TaxonID
  # Merge taxonomy table and relAbundDF by rows that they share
  mergedDf_1 <- merge(taxTable, relAbundDf, by=0)
  rownames(mergedDf_1) <- mergedDf_1[,1] #make the column "row.names" the actual rownames again
  mergedDf_1[,1] <- NULL #remove "row.names" column
  # Get the last 3 rows of myMultipatt$sign, which are index, stat, and p-value:
  accessSign <- myMultipatt$sign[,c((ncol(myMultipatt$sign)),(ncol(myMultipatt$sign) -1),(ncol(myMultipatt$sign) -2))]
  # Merge mergedDf_1 with accessSign to add index, stat, and p-value
  mergedDf_2 <- merge(mergedDf_1,accessSign, by=0)
  results <- mergedDf_2[with(mergedDf_2, order(index)),]
  rownames(results) <- results[,1] #make row names taxon ID
  results[,1] <- NULL #remove "row.names" column
  return(results)
}

# 2. PERFORM INDICATOR SPECIES ANALYSIS
# Make a vector that lists timepoints for each sample (within hopped and unhopped data sets)
hoppedTimepoints <- sample_data(beerHopped_ps)$timepoint
hoppedTimepoints

unhoppedTimepoints <- sample_data(beerUnhopped_ps)$timepoint
unhoppedTimepoints

# Perform indicator species analysis:
# First, get ASV ans taxonomy tables out of phyloseq:
# A very handy little function from:
# https://jacobrprice.github.io/2017/08/26/phyloseq-to-vegan-and-back.html
psotu2veg <- function(physeq) {
  OTU <- otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU) #invert columns and rows
  }
  return(as(OTU, "matrix"))
}

# get ASV table for hopped
hoppedASVs <- psotu2veg(beerHopped_ps)

# get ASV table for unhopped
unhoppedASVs <- psotu2veg(beerUnhopped_ps)

# get tax table for hopped
# Remove phyloseq attributes from taxonomy file that was trimmed down automatically by phyloseq function:
hoppedTax <- as.data.frame(phyloseq::tax_table(beerHopped_ps), stringsAsFactors = F)
head(hoppedTax)

# get tax table for unhopped
unhoppedTax <- as.data.frame(phyloseq::tax_table(beerUnhopped_ps), stringsAsFactors = F)
head(unhoppedTax)

set.seed(9)
IndicHopsTime <- multipatt(x=hoppedASVs, cluster=hoppedTimepoints, func="r.g", control=how(nperm = 9999)) 
summary(IndicHopsTime) #shows all of the taxon IDs for ASVs associated with each group, as well as
# the  r.g. value (under stat, this is the correlation index that takes into account the
# variation within and among groups), and the p-value from a permutation test.
#View(IndicByRange)

set.seed(9)
IndicUnhopTime <- multipatt(x=unhoppedASVs, cluster=unhoppedTimepoints, func="r.g", control=how(nperm = 9999)) 
summary(IndicUnhopTime) #shows all of the taxon IDs for ASVs associated with each group, as well as
# the  r.g. value (under stat, this is the correlation index that takes into account the
# variation within and among groups), and the p-value from a permutation test.
#View(IndicByRange)

# 3. PERFORM INDICATOR SPECIES ANALYSIS
# i. hopped samples
# Apply function to the multipatt object created above:
HopTimeIStable <- indSpecTable2(myMultipatt=IndicHopsTime, taxTable= hoppedTax, ASVtable = hoppedASVs, signLevel=0.05)
#View(ByRangeISTable)
HopTimeIStable$index
# Index numbers correspond to the groups that were found above (see summary(IndicByRange)), in order.
summary(IndicHopsTime)
# Change these accordingly
HopTimeIStable$index[which(HopTimeIStable$index==1)] <- "1 week"
HopTimeIStable$index[which(HopTimeIStable$index==2)] <- "3 weeks"

#View(HopTimeIStable)

# ii. unhopped samples
# Apply function to the multipatt object created above:
unhopTimeIStable <- indSpecTable2(myMultipatt=IndicUnhopTime, taxTable= unhoppedTax, ASVtable = unhoppedASVs, signLevel=0.05)
#View(unhopTimeIStable)
unhopTimeIStable$index
# Index numbers correspond to the groups that were found above (see summary(IndicByRange)), in order.
unhopTimeIStable$index[which(HopTimeIStable$index==1)] <- "1 week"
unhopTimeIStable$index[which(HopTimeIStable$index==2)] <- "3 weeks"
#View(unhopTimeIStable)




# QUESTION 5: Does the fungal/bacterial ratio change across week?


# iv. Is a difference in variation among unhopped samples indoors or outdoors?
# week 1
mod_dispweek1 <- betadisper(week1_BCdist, group= factor(as.data.frame(as.matrix(sample_data(beerWeek1_ps)))$treatment), 
                            type = "centroid")
# test it
# set seed!!
set.seed(93)
anova(mod_dispweek1, permutations = 9999) #NO P = 0.1745
# So this means that within week 1 samples, there is no difference among treatments in terms of their dispersions
mod_dispweek1

plot(mod_dispweek1)

# week 3
mod_dispweek3 <- betadisper(week3_BCdist, group= factor(as.data.frame(as.matrix(sample_data(beerWeek3_ps)))$treatment), 
                            type = "centroid")
# test it
# set seed!!
set.seed(93)
anova(mod_dispweek3, permutations = 9999) #YES!!! P = 0.001788 **

mod_dispweek3

plot(mod_dispweek3)

# Trying the same thing as above with quantitative Jaccard as the dissimilarity index
mod_dispweek1_QJ <- betadisper(week1_QJdist, group= factor(as.data.frame(as.matrix(sample_data(beerWeek1_ps)))$treatment), 
                               type = "centroid")
# test it
# set seed!!
set.seed(93)
anova(mod_dispweek1_QJ, permutations = 9999) #NO!! P = 0.2585
# So this means that within week 1 samples, there is no difference among treatments in terms of their dispersions
mod_dispweek1_QJ

plot(mod_dispweek1_QJ)

# iv. Is a difference in variation among unhopped samples indoors or outdoors?
mod_dispweek3_QJ <- betadisper(week3_QJdist, group= factor(as.data.frame(as.matrix(sample_data(beerWeek3_ps)))$treatment), 
                               type = "centroid")
# test it
# set seed!!
set.seed(93)
anova(mod_dispweek3_QJ, permutations = 9999) #YES!!! P = 0.001308 **
# So this means that within week 3, there is a significant effect of treatment! 
mod_dispweek3_QJ

plot(mod_dispweek3_QJ)

centroids <- cbind(week1_BC_vec, week3_BC_vec , week1_qJ_vec, week3_qJ_vec)
rownames(centroids) <- c("HI", "HO", "UI", "UO")
