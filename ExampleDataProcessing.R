# 16S ONLY Exploratory Data Analysis and Data Clean Up
# Aug 17, 2021

# This script is to take a first look at 16S MiSeq data from the samples that I
# took in May 2021 to examine an edge effect across forest and savanna boun-
# daries at the Savannah River Site. Here, I look at sequence counts, rarefy,
# examine top classes and phyla, and look at ordinations. I repeat the ordinations
# and examination of top classes and phyla before and after removing rare species
# (here those that showed up LESS THAN 50 times across whole dataset).
# Then, I move forward with the rarefied, rare-removed data set and do preliminary
# investigations into why outlier samples might be outliers, how savannas/patches
# differ from forests, etc.

# ****IMPORTANT NOTE****: THE DIFFERENTIAL ABUNDANCE ANALYSIS PERFORMED HERE IS *NOT* THE ONE
# USED FOR LATER ANALYSES OF THE EDGE EFFECT TAXA.

# FILES SAVED IN THIS SCRIPT (both in Bioinformatics directory: (save code at end of script,
# but commented out so that new files are not saved every time in Bioinformatics folder):
# 1. "EDA16SAug2021"
# rarefied.ps, relabun.phylatop99.5, relabun.phylatop99, top_99.5p_phyla, 
# relabun.classtop95, ord, ordSoils, ASVsTrimmed, taxTrimmed.

# 2. "trimmedJustsoils.ps" 
# This is the phyloseq object that was made after rarefying, and keeping
# only ASVs that occurred at least 50 times across the (rarefied) dataset, including
# soils and controls (NOT biocrusts). No outliers were removed. 

# FUNCTIONS CREATED IN THIS SCRIPT:
# 1. Function to get the taxonomy table out of phyloseq:
# (Inspired by similar function at https://jacobrprice.github.io/2017/08/26/phyloseq-to-vegan-and-back.html)
taxtable_outta_ps <- function(physeq){ #input is a phyloseq object
  taxTable <- tax_table(physeq)
  return(as.data.frame(taxTable))
}

# 2. Function to get ASV table out of phyloseq so that we can view it better
# (Inspired by similar function at https://jacobrprice.github.io/2017/08/26/phyloseq-to-vegan-and-back.html)
ASVs_outta_ps <- function(physeq){ #input is a phyloseq object
  ASVTable <- otu_table(physeq)
  return(as.data.frame(ASVTable))
}

# 3. Function that gets OTU/ASV table out of phyloseq
# (from https://jacobrprice.github.io/2017/08/26/phyloseq-to-vegan-and-back.html)
psotu2veg <- function(physeq) {
  OTU <- otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

# 4. Function that gets average pH for each sample, based on the 2-3 pH readings taken
mean_pH_func <- function(pHs_mat) { #enter a dataframe that has 1-3 pHs per row
  # this function works for averaging 2-3 pHs (per row), but could easily be adjusted in
  # for loop to be more flexible.
  H_mat <- matrix(nrow=nrow(pHs_mat), ncol=ncol(pHs_mat))
  colnames(H_mat) <- colnames(pHs_mat)
  for(i in 1:nrow(pHs_mat)){ #get the H+ concentration for each pH reading
    H_mat[i,1] <- 10^(pHs_mat[i,1]*-1)
    H_mat[i,2] <- 10^(pHs_mat[i,2]*-1)
    H_mat[i,3] <- 10^(pHs_mat[i,3]*-1)
  }
  avgH <- rowMeans(H_mat, na.rm=TRUE) #get the average H+ concentration for each set of pHs
  avgpH <- -log10(avgH) #convert back to pH
  return(avgpH)
}

#######################################################################################

# Set working directory
setwd("~/Desktop/CU_Research/SoilEdgeEffectsResearch")

# Read in libraries
library("phyloseq")
library("readxl")       #necessary to import the data from Excel file
library("tidyverse")
library("vegan")
library("gridExtra")    #allows you to make multiple plots on the same page with ggplot
library("DESeq2") #for differential abundance analysis

##################################################################################
# I. SET-UP, DATA CLEANING, RAREFACTION, AND FIRST TAXONOMIC & ORDINATION PLOTS
##################################################################################

list.files()
seqtab_wTax_mctoolsr <- read.table("Bioinformatics/seqtab_wTax_mctoolsr.txt")
str(seqtab_wTax_mctoolsr)
head(seqtab_wTax_mctoolsr)
#View(seqtab_wTax_mctoolsr)
colnames(seqtab_wTax_mctoolsr)
rownames(seqtab_wTax_mctoolsr)

seqtab <- read.table("Bioinformatics/seqtab_final.txt", header=T)
dim(seqtab) #38,255 ASVs and 270 samples (soil samples, controls, etc.)
#View(seqtab)

#### LOAD FILES IN FOR PHYLOSEQ #####
# 1. OTU table: ASVs/OTUs as rows, samples as columns. This is "seqtab"
otu_mat <- seqtab
seqtab$X #X is the ASV names
otu_mat <- otu_mat %>%
  tibble::column_to_rownames("X") #make it so that "X"(which is column with ASV names is the rownames columnn!)
colnames(otu_mat)
rownames(otu_mat)
#View(otu_mat) #ASVs are rows; samples (with weird X in front of number) are columns

# 2. OTU taxonomy
# This is basically the same as the first and last columns of seqtab_wTax_mctoolsr
# So none of the guts that are essentially the ASV table
#View(seqtab_wTax_mctoolsr) #271 columns, labeled V1 through V271. V1 is bascially "rownames", ie. each ASV ID, V2:V270 are samples; V271 is taxonomy info.
colnames(seqtab_wTax_mctoolsr) #samples
head(seqtab_wTax_mctoolsr$V1) #THE ASV names!
head(seqtab_wTax_mctoolsr$V271) #the taxonomy!
step1 <- cbind(seqtab_wTax_mctoolsr$V1, seqtab_wTax_mctoolsr$V271) #this gets just the ASV names and the taxonomy associated with each
head(step1)
colnames(step1) <- step1[1,] #make ASV_ID and taxonomy the column names
step1 <- step1[-1,] #get rid of first row which is now same as column names
head(step1)
str(step1)
#View(step1)
# as.data.frame(step1) #This didn't work so I foreced it below

# Make new columns based on the ";" separator to separate out taxonomy
require(tidyr)
tax_sep <- separate(as.data.frame(step1), col = taxonomy, into= c("Kingdom", "Phylum", "Class", 
                                                                  "Order", "Family", "Genus", "Species"),
                    sep = ";")

# View(tax_sep)
str(tax_sep)
all.equal(tax_sep$`#ASV_ID`, tax_sep$Species) #Yep
colnames(tax_sep)
#tax_sep$Species <- NA
tax_sep <- tax_sep[,1:7] #get rid of species column, which is the same as ASV_ID anyway
#View(tax_sep)

# Format tax_sep for phyloseq by making ASV ID the rownames
tax_mat <- tax_sep %>%
  tibble::column_to_rownames("#ASV_ID")

# 3. Sample metadata
metadata <- read.csv("Bioinformatics/SRS_AllMetadataAug17.csv") #this new csv has ALL sample metadata, but is missing some mean pH values
# not just for biological samples (i.e. controls too)
# View(metadata)
colnames(metadata)
rownames(metadata)
metadata <- metadata[,-1] #get rid of weird X column that is the same as the rownames
samples_df <- metadata %>%
  tibble::column_to_rownames("MappingSampID") #rownames/sample ID should be the mapping sample IDs!
str(samples_df)
rownames(samples_df)
colnames(samples_df)
samples_df$mean_pH <- mean_pH_func(samples_df[,3:5]) #get mean pH for each sample
# View(samples_df) #several of these were hand checked to make sure that function was working correctly

# Need to rename samples_df to match those in otu_mat; because they don't match
# each other except for the blanks and controls, those are the only ones that got
# phyloseq'd
length(rownames(samples_df)) == length(colnames(otu_mat))
#Will remove x's from the colnames on otu_mat
newnames <- str_remove(colnames(otu_mat), "X") #use grep to remove Xs
colnames(otu_mat) <- newnames
all.equal(sort(colnames(otu_mat)), sort(rownames(samples_df))) #Cool, names are now the same

# Transform otu table and taxonomy tables into matrices
otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)
#View(otu_mat)

# Transform to phyloseq objects
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(samples_df)
SRS_16S_raw.ps <- phyloseq(OTU, TAX, samples)
SRS_16S_raw.ps 
colnames(otu_table(SRS_16S_raw.ps))

sample_names(SRS_16S_raw.ps) # IMPORTANT: THESE ARE THE NAMES BASED ON THE MAPPING
# FILE USED FOR DEMULITPLEXING; I.E. THEY CORRESPOND TO ORDER LOADED INTO ROWS
# OF PLATES (*NOT NOT NOT*) DNA TUBE LABELS 

rank_names(SRS_16S_raw.ps)
sample_variables(SRS_16S_raw.ps)

##################################################################
# FILTER OUT MITOCHONDRIA, CHLOROPLASTS,a and KINGDOM EUKARYOTA
# chloroplast is an order and mitochondria is a family in this version of SILVA
##################################################################

# Get taxonomy table out of phyloseq:
SRS_taxTable <- taxtable_outta_ps(SRS_16S_raw.ps)
#View(SRS_taxTable)

# These lines below (lines 199-207) are technically not needed, since I found a way to do it with phyloseq
# tax_noeuksorNAs is no longer used, but serves as a good check to make sure that phyloseq
# is doing what I want it to. 
# Remove chloroplasts, mitochondria, and all eukaryotes:
tax_noeuksorNAs <- SRS_taxTable %>% 
  filter(Order != "Chloroplast") %>% 
  filter(Family != "Mitochondria") %>% 
  filter(Kingdom != "Eukaryota") %>% 
  filter(Kingdom != "NA") %>% #Remove ASVs where kingdom is unknown ("NA" in first column) 
  filter(Phylum != "NA") #Remove ASVs where phylum is unknown ("NA" in first column) 
# View(tax_noeuksorNAs) 
dim(SRS_taxTable)[1] - dim(tax_noeuksorNAs)[1]
# We filtered out a total of 3255 taxa.

# We now have to trim the OTU table, because it still has the taxa in it that were filtered out
# above. It currently has 38255 taxa, whereas tax_noeuksorNAs has 38255 taxa. We'll want it this 
# way so that we re-make the phyloseq object and then rarefy correctly. Could have done this all in
# same step as above, but oh well!

# Find ASVs where the order is Chloroplast
chloros <- SRS_taxTable %>% 
  filter(Order == "Chloroplast")
dim(chloros) #149 chloroplast ASVs
chloro_names <- rownames(chloros)
# Find ASVs where the family is Mitochondria
mitos <- SRS_taxTable %>% 
  filter(Family == "Mitochondria")
dim(mitos) #959 mitochondria ASVs
mito_names <- rownames(mitos)
# Find ASVs where kingdom is Eukaryota
kingdomEuks <- SRS_taxTable %>% 
  filter(Kingdom == "Eukaryota")
dim(kingdomEuks) #482
euks_names <- rownames(kingdomEuks)
# Find ASVs where kingdom is NAS
kingdomNAs <- SRS_taxTable %>% 
  filter(Kingdom == "NA")
dim(kingdomNAs) #490
kNAs_names <- rownames(kingdomNAs)
# Find ASVs where phylum is NAS
PhylumNAs <- SRS_taxTable %>% 
  filter(Phylum == "NA")
dim(PhylumNAs) #2147
pNAs_names <- rownames(PhylumNAs)

149 + 959 + 482 + 490 + 2147 #4227, this is more than what we removed above because some kingdom that were NAs had NA phyla too
bad_ASVs <- c(chloro_names, mito_names, euks_names, kNAs_names, pNAs_names)
length(bad_ASVs) == 149 + 959 + 482 + 490 + 2147 

# Inspiration for using phyloseq here (Joey711's code: https://github.com/joey711/phyloseq/issues/652
# Remove the ASVs listed above:
all_Taxa <- taxa_names(SRS_16S_raw.ps) #get all tax names in original, uncleaned dataset
ASVstoKeep <- all_Taxa[!(all_Taxa %in% bad_ASVs)]
length(ASVstoKeep) #35000
noeuksorNAs_ps <- prune_taxa(ASVstoKeep, SRS_16S_raw.ps) #new phyloseq object with just the stuff we want!

# TAKING A LOOK AT THE DATA
# How are number of reads distributed across the samples?
otu_table(noeuksorNAs_ps)
seqspersample <- colSums(otu_table(noeuksorNAs_ps))
SeqNumberPlot <-barplot(seqspersample, main="All: Total Sequences Per Sample", ylab= "Number of Sequences", xlab= "Sample Number")

# Just controls:
seqsperctrl <- seqspersample[244:length(seqspersample)]
SeqNumberPlotCtrls <- barplot(seqsperctrl, main="Total Sequences Per Sample", ylab= "Number of Sequences", xlab= "Sample Number", cex.names=0.27)
SeqNumberPlotCtrls

# Controls with easier axis names
simple <- c("ExtB1", "ExtB2", "ExtB3", "ExtW1", "ExtW10", "ExtW11", "ExtW12",
            "ExtW13", "ExtW14", "ExtW15", "ExtW16", "ExtW17", "ExtW18", "ExtW19",
            "ExtW2", "ExtW20", "ExtW21","ExtW22", "ExtW3","ExtW4",
            "ExtW5","ExtW6","ExtW7","ExtW8", "ExtW9", "PCRNTC")
seqsperctrl_simple <- setNames(seqsperctrl, simple)

##quartz()
barplot(seqsperctrl_simple, main="Controls: Total Sequences Per Sample",
        ylab= "Number of Sequences", xlab= "Sample Number", cex.names=0.35)

######## RAREFACTION ##########
names(seqspersample) #shows the order of all the samples, and the ones to pull out bel
# What is the minimum for the non-control samples?
min(seqspersample[1:243]) #1224
max(seqspersample[1:243]) #157741
mean(seqspersample[1:243]) # 26730.58

# What does all the spread look like for ALL samples?
min(seqspersample) #77
max(seqspersample) #157741
mean(seqspersample) # 24435.77
sd(seqspersample) #13421.05

min(seqspersample[244:length(seqspersample)]) #77
max(seqspersample[244:length(seqspersample)]) #15769
mean(seqspersample[244:length(seqspersample)]) # 2988.115
sd(seqspersample[244:length(seqspersample)]) #5276.02

# Plot this:
##quartz()
seqsPerExSamp <- seqspersample[1:243]
SeqNumberPlotExSamp <- barplot(seqsPerExSamp, main="Soils: Total Sequences Per Sample",
                               ylab= "Number of Sequences", xlab= "Sample Number", cex.names=0.4)

sort(seqsPerExSamp) #take a look in order from least to greatest
sort(seqspersample) #take a look in order from least to greatest (all samples and controls)

# We'll rarefy at 14973, which only costs us 7 samples (plus all but two of the controls)
# The "random rarefaction is made without replacement so that the variance of rarefied
# communities is rather related to rarefaction proportion than to the size of the sample". 
# (quoted bit from vegan's manual page)

set.seed(19)
rarefied.ps <- rarefy_even_depth(noeuksorNAs_ps, sample.size = 14973, replace=FALSE, trimOTUs=TRUE)
# 1121OTUs were removed because they are no longer present in any sample after random subsampling

# Rarefaction curve:
samp.col <- c(rep("blue", 243), rep("grey", 26))

rare.plot <- rarecurve(t(otu_table(noeuksorNAs_ps)), step = 3000, cex = 0.5, col = samp.col, label = FALSE, xlab = "Number of Sequences", ylab = "Number of ASVs")
rarecurve(t(otu_table(noeuksorNAs_ps)), step = 3000, cex = 0.5, col = samp.col, label = FALSE, xlab = "Number of Sequences", ylab = "Number of ASVs")

# How many are left?
colSums(otu_table(rarefied.ps)) #All have 14973
length(colSums(otu_table(rarefied.ps))) #238 samples left (2 less than when I rarefied w/o removing phylum= NA)
# Which samples are these:
sampleNames <- names(colSums(otu_table(rarefied.ps)))
sampleNames

####################################################################
# The following lines through about line 513 is all done BEFORE removing rare taxa,
# (i.e., with rarefied dataset, but BEFORE removing those that do not occur at least 50 
# times in the dataset))

#####################
# TOP PHYLA
#####################

# Phylum level (adopted from my code at:
# https://github.com/clairecwinfrey/PhanBioMS_scripts/blob/master/R_scripts/figures/taxonomic_barplots.R)
# TURN ASVs INTO PHYLUM LEVEL
rarefied.phylum.glom <-  tax_glom(rarefied.ps, taxrank = "Phylum") 
tax_table(rarefied.phylum.glom) # 

# TRANSFORM SAMPLE COUNTS ON JUST GLOMMED SAMPLES (UNLIKE WHAT WE DID AT FIRST)
relabun.phyla.0 <- transform_sample_counts(rarefied.phylum.glom, function(x) x / sum(x) )
rownames(otu_table(relabun.phyla.0)) #weirdly, this is ASV_.... Is this right? 
# I think that the ASVs are just representative from each phylum

# MERGE SAMPLES so that we only have combined abundances for site and different kinds of controls
relabun.phyla.1 <- merge_samples(relabun.phyla.0, group = "EU")
sample_data(relabun.phyla.1) #shows that we still have samples from each EU, biocrust, and extcontrol (water)

# CONVERT TO PROPORTIONS AGAIN B/C TOTAL ABUNDANCE OF EACH SITE WILL EQUAL NUMBER OF SPECIES THAT WERE MERGED
relabun.phyla.2 <- transform_sample_counts(relabun.phyla.1, function(x) x / sum(x))
sample_data(relabun.phyla.2)

# NOW GET ONLY TAXA THAT COMPRISE AT LEAST 1% OF THE ABUNDANCE 
relabun.phyla.df <-psmelt(relabun.phyla.2)
dim(relabun.phyla.df) #
sum(relabun.phyla.df[,3]) 
colnames(relabun.phyla.df) 
relabun.phylatop99 <- relabun.phyla.df
relabun.phylatop99.5 <- relabun.phyla.df

relabun.phylatop99$Phylum[relabun.phylatop99$Abundance < 0.01] <- "< 1% abund."
relabun.phylatop99.5$Phylum[relabun.phylatop99.5$Abundance < 0.005] <- "< .5% abund."

top_99p_phyla <- unique(relabun.phylatop99$Phylum)
top_99p_phyla

top_99.5p_phyla <- unique(relabun.phylatop99.5$Phylum)
top_99.5p_phyla

# Surprised that Gemmatimonadetes not above >1%! But Gemmatimonadetes in top 16 and
# comprises at least 0.5% of the total abundance

# Phyla comprising at least 0.5% of total abundance
phylumPlot99.5percent <- ggplot(data=relabun.phylatop99.5, aes(x=Sample, y=Abundance, fill=Phylum)) + theme(axis.title.y = element_text(size = 14, face = "bold")) + theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(colour = "black", size = 12, face = "bold"))
#quartz()
phylumPlot99.5percent + geom_bar(aes(), stat="identity", position="fill") +
  theme(legend.position="bottom") +
  guides(fill=guide_legend(nrow=4)) + theme(legend.text = element_text(colour="black", size = 10))  + ggtitle("Phyla comprising at least 0.5% of total abundance")

# "Phyla comprising at least 1% of total abundance"
phylumPlot.99percent <- ggplot(data=relabun.phylatop99, aes(x=Sample, y=Abundance, fill=Phylum)) + theme(axis.title.y = element_text(size = 14, face = "bold")) + theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(colour = "black", size = 12, face = "bold"))
#quartz()
phylumPlot.99percent + geom_bar(aes(), stat="identity", position="fill") +
  theme(legend.position="bottom") +
  guides(fill=guide_legend(nrow=4)) + theme(legend.text = element_text(colour="black", size = 10))  + ggtitle("Phyla comprising at least 1% of total abundance")

## Below is junk from earlier script... keep for easier manipulation later!
# scale_fill_manual(values = c("#4575b4", "#d73027", "#fc8d59", "#fee090", "#91bfdb", "grey"), 
#   name= "Phylum", breaks= c("D_1__Firmicutes", "D_1__Proteobacteria", "D_1__Bacteroidetes", "D_1__Actinobacteria", "D_1__Fusobacteria", "< 1% abund."), 
#  labels =c("Firmicutes", "Proteobacteria", "Bacteroidetes", "Actinobacteria", "Fusobacteria", "< 1% abund.")) +

# Get exact abundances of each phyla (top 99.5%):
colnames(relabun.phylatop99.5)
relabun.phylatop99.5[,2] #This is EU... now "Sample" because of the glomming!

top_99.5p_phyla <- relabun.phylatop99.5 %>%
  group_by(Sample, Phylum) %>%
  summarize(Mean = mean(Abundance)) %>%
  arrange(-Mean) 
# View()

#####################
# TOP CLASSES:
#####################

# (adopted from my code at:
# https://github.com/clairecwinfrey/PhanBioMS_scripts/blob/master/R_scripts/figures/taxonomic_barplots.R)
# TURN ASVs INTO class LEVEL
rarefied.class.glom <-  tax_glom(rarefied.ps, taxrank = "Class") 
tax_table(rarefied.class.glom) # good, this is only class (42 different phyla!)

# TRANSFORM SAMPLE COUNTS ON JUST GLOMMED SAMPLES (UNLIKE WHAT WE DID AT FIRST)
relabun.class.0 <- transform_sample_counts(rarefied.class.glom, function(x) x / sum(x) )
rownames(otu_table(relabun.class.0)) # I think that the ASVs are just representative from each class

# MERGE SAMPLES so that we only have combined abundances for site and different kinds of controls
relabun.class.1 <- merge_samples(relabun.class.0, group = "EU")
sample_data(relabun.class.1) #shows that we still have samples from each EU, biocrust, and extcontrol (water)

# CONVERT TO PROPORTIONS AGAIN B/C TOTAL ABUNDANCE OF EACH SITE WILL EQUAL NUMBER OF SPECIES THAT WERE MERGED
relabun.class.2 <- transform_sample_counts(relabun.class.1, function(x) x / sum(x))
sample_data(relabun.class.2)

# NOW GET ONLY TAXA THAT COMPRISE AT LEAST 1% OF THE ABUNDANCE 
relabun.class.df <-psmelt(relabun.class.2)
dim(relabun.class.df) #
sum(relabun.class.df[,3]) 
colnames(relabun.class.df) 
relabun.classtop99 <- relabun.class.df
relabun.classtop95 <- relabun.class.df

relabun.classtop99$Class[relabun.classtop99$Abundance < 0.01] <- "< 1% abund."
relabun.classtop95$Class[relabun.classtop95$Abundance < 0.05] <- "< 5% abund."

top_99p_class <- unique(relabun.classtop99$Class)
top_99p_class

top_95p_class <- unique(relabun.classtop95$Class)
top_95p_class

# Phyla comprising at least 5% of total abundance
classPlot.95pt <- ggplot(data=relabun.classtop95, aes(x=Sample, y=Abundance, fill=Class)) + theme(axis.title.y = element_text(size = 14, face = "bold")) + theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(colour = "black", size = 12, face = "bold"))
#quartz()
classPlot.95pt + geom_bar(aes(), stat="identity", position="fill") +
  theme(legend.position="bottom") +
  guides(fill=guide_legend(nrow=4)) + theme(legend.text = element_text(colour="black", size = 5.5))  + ggtitle("Classes comprising at least 5% of total abundance")

# Ignore the lines below... not currently using them
# Bray-Curtis dissimilarities based on square-root transformed data
# Code from mctoolsR
#braydist <- vegan::vegdist(t(otu_table(rarefied.ps)), method = "bray")
#summary(braydist)
# bray.mat <- as.matrix(braydist)
# View(bray.mat) looks as expected
# Plot ordination:

# ord <- calc_ordination(braydist, 'nmds')

# Using a phyloseq tool for convenience!
# Bray-Curtis dissimilarities based on square-root transformed data
set.seed(19)
ord <- ordinate(rarefied.ps, method = "NMDS", distance = "bray", trymax = 100)

# With no black outline around points
rarefiedBrayNMDS <- phyloseq::plot_ordination(rarefied.ps, ord, type= "samples", color= "EU")
#quartz()
rarefiedBrayNMDS + geom_polygon(aes(fill=EU)) + geom_point(size=3) + ggtitle("NMDS based on Bray-Curtis Dissimilarities")

# With black outline around points:
rarefiedBrayNMDSoutlined <-rarefiedBrayNMDS + geom_polygon(aes(fill=EU)) + geom_point(aes(fill=EU),color="black",pch=21, size=3) + ggtitle("NMDS based on Bray-Curtis Dissimilarities")
#quartz()
rarefiedBrayNMDSoutlined

# With transect included as shape:
# This plot is not very informative!
rarefiedBrayNMDStran <- rrarefiedBrayNMDS <- phyloseq::plot_ordination(rarefied.ps, ord, type= "samples", color= "EU", shape= "Transect")
#quartz()
rarefiedBrayNMDStran + geom_polygon(aes(fill=EU)) + geom_point(size=3) + ggtitle("NMDS based on Bray-Curtis Dissimilarities")

# Now, remove biocrust and controls from ordination:
justsoils.ps <- subset_samples(rarefied.ps, Type != "BioCrust" & Type != "ExtContWater")
unique(sample_data(justsoils.ps)[,27]) #only soils
sample_names(justsoils.ps) #another check to show that we have only soils!

# Ordination plot of just soils:
set.seed(19)
ordSoils <- ordinate(justsoils.ps, method = "NMDS", distance = "bray", trymax = 100)
soilsrarefiedBrayNMDS <- phyloseq::plot_ordination(justsoils.ps, ordSoils, type= "samples", color= "EU")
#quartz()
soilsrarefiedBrayNMDS + geom_polygon(aes(fill=EU)) + geom_point(size=3) + ggtitle("NMDS based on Bray-Curtis Dissimilarities")

# Add in labels to figure out what weird samples are
soilsrarefiedBrayNMDS <- phyloseq::plot_ordination(justsoils.ps, ordSoils, type= "samples", color= "EU", shape= "Habitat", label = "Sample.ID")
#quartz()
soilsrarefiedBrayNMDS + geom_polygon(aes(fill=EU)) + geom_point(size=3) + ggtitle("NMDS based on Bray-Curtis Dissimilarities")

# Visible outliers: (going clockwise from the top left on the ordination plot above):
# 53ND_B_20, 10C_L_10, 53ND_R_40, 53SD_R_10, 52D_R_10, 52D_L_100, (Maybe!) 52D_R_90, 53ND_L_80
outliers <- c("53ND_B_20", " 10C_L_10", "53ND_R_40", "53SD_R_10", "52D_L_100", "52D_R_90", "53ND_L_80")

##################################################################################
# II. REMOVAL OF "RARE" TAXA, NEW TAXONOMIC PLOTS AND ORDINATIONS
# (tried out different numbers of taxa, but settled on removing taxa that did not
# occur at least 50 times across the whole dataset)
##################################################################################

###########################################################################
# REMOVE "RARE" TAXA AND THEN RE-RUN 1) sequences per sample, 
# 2) top phyla, 3) top classes, 4) ordinations
# Questions when comparing pre and post removal of rare taxa:
# 1) Do top phyla or classes change, or outliers, after this removal? If so, this is
# evidence of these shifts being driven by rare taxa:

rarefiedNoBC.ps <- subset_samples(rarefied.ps, Type != "BioCrust") #make a phyloseq object that does not have biocrusts, but
# that retains controls. 

rare_ASVtab <- ASVs_outta_ps(rarefiedNoBC.ps) 
# Add column for abundance of ASVs across all (pre-rarefied) samples
rare_ASVtab$Abundance <- rowSums(rare_ASVtab)
dim(rare_ASVtab) #33,878 ASVs
# View(rare_ASVtab)

# Most are rare! 
plot(rare_ASVtab$Abundance)
length(which(rare_ASVtab$Abundance <= 50)) #23,903 ASVs have 50 or fewer occurrences across the rarefied data set 
length(which(rare_ASVtab$Abundance <= 35)) #20,882 ASVs have 30 or fewer occurrences across the rarefied data set
length(which(rare_ASVtab$Abundance <= 10)) #9,946 ASVs have 10 or fewer occurrences across the rarefied data set

length(which(rare_ASVtab$Abundance >= 50)) #10137 ASVs have 50 or more occurrences across the rarefied data set 

# We'll get rid of the ASVs that do not occur AT LEAST 50 times across our dataset:
keptASVsindex <- which(rare_ASVtab$Abundance >= 50) #gives row numbers to keep in the dataset 
max(keptASVsindex) # last row to keep is 15719 (but not all before are kept probably b/c of rarfying)
rare_ASVtab$Abundance[15719] #has exactly fifty 

ASVsTrimmed <- rare_ASVtab[keptASVsindex,] #keep only rows of the index, and all columns (i.e. all samples)
#View(ASVsTrimmed) #10137 ASVs across 235 samples, as expected 
plot(ASVsTrimmed$Abundance) # tail is still long, because most ASVs are still rare
# relative to the most abundant ASVs

length(which(ASVsTrimmed$Abundance == 50)) # 161 ASVs were right at cut off!
range(ASVsTrimmed$Abundance) #50 74724
length(which(ASVsTrimmed$Abundance >= 1000)) #only 513 ASVs appear more than 1000 times
# across the data set, which is in average of about 4 times per sample!

# Get TAXONOMIC TABLE that matches the ASV table above using function we defined earlier
rare_taxTab <- taxtable_outta_ps(rarefied.ps) 
# (first check that the ASVs (i.e. rownames) are in the same order for the two data frames so that we can use the 
# same index)
unique(rownames(rare_taxTab) == rownames(rare_ASVtab)) # all are the same, son we're good to move to the next step 
taxTrimmed <- rare_taxTab[keptASVsindex,] #keep only rows of the index made above for tax appearing at least 50 times
unique(rownames(taxTrimmed) == rownames(ASVsTrimmed)) #ALl true so this is the same!

# How many ASVs occur in every sample?
# Said another way: how many rows (i.e. ASVs) occur >= 1 time across all columns (i.e. samples)
# Or: are there any rows that contain no zeros?
# Inspiration from this post: https://stackoverflow.com/questions/9977686/how-to-remove-rows-with-any-zero-value

# Let's you find rows that contain zeros and then filter them out
rownames(ASVsTrimmed[!(apply(ASVsTrimmed, 1, function(y) any(y == 0))),])
# "ASV_1"  "ASV_3"  "ASV_4"  "ASV_12" "ASV_17" are found in every sample 
# This corresponds to rows:
ASVsinAllindex <- c(1, 3, 4, 10, 15)
rownames(ASVsTrimmed)[ASVsinAllindex] # yes, this matches above
# What are these ASVs?
rare_taxTab[ASVsinAllindex,] #Candidatus_Udaeobacter ASV, Mycobacterium ASV, ASV within family Xanthobacteraceae, Bradyrhizobium ASV, another ASV within family Xanthobacteraceae
rowMeans(ASVsTrimmed)[ASVsinAllindex] #decent average amount in each sample
# ASV_1     ASV_3     ASV_4    ASV_12    ASV_17 
# 633.25424 279.06780 258.35593 107.32203  87.07627 

coreASVs <- cbind.data.frame(rare_taxTab[ASVsinAllindex,], rowMeans(ASVsTrimmed)[ASVsinAllindex])
colnames(coreASVs)[7] <- "meanAbundancePerSample"
#View(coreASVs)

# Merge trimmed ASV tables and taxonomy tables to get something that is formatted like
# "seqtab_wTax_mctoolsr"
seqtab_wTax_trimmed <- cbind.data.frame(ASVsTrimmed, taxTrimmed)
seqtab_wTax_trimmed[,239]

# Make Excel files of this! #last re-written Sept., 29, 2022
#write.csv(ASVsTrimmed[,-239], "SRSMay2021_16S_cleanedASVtable.csv")
#write.csv(taxTrimmed, "SRSMay2021_16S_cleanedTaxTable.csv") 
#write.csv(seqtab_wTax_trimmed[,-239], "SRSMay2021_16S_seqtab_wtax_trimmed.csv")

# Do these ASVs that we'll keep differ if we remove the controls?
rarefiedOnlyControls <- subset_samples(rarefied.ps, Type != "BioCrust" & Type != "ExtContWater")

############ RE-DO ANALYSES WITH DATASET WHERE THE RARE TAXA HAVE BEEN REMOVED (couting soils and controls
#### NEED TO REMAKE THESE INTO PHYLOSEQ OBJECTS FIRST (READ THROUGH THIS TOO)

class(taxTrimmed)
dim(ASVsTrimmed)
#View(ASVsTrimmed)

otu_mat <- as.matrix(ASVsTrimmed[,-239]) #remove "Abundance" column that is at the end
tax_mat <- as.matrix(taxTrimmed)
#View(otu_mat)
#View(tax_mat)

# Transform to phyloseq objects
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(samples_df)
TrimmedSRS_16S.ps <- phyloseq(OTU, TAX, samples)
TrimmedSRS_16S.ps #10137 taxa and 235 samples!
colnames(otu_table(TrimmedSRS_16S.ps))

#####################
# trimmedTop PHYLA
#####################

# Phylum level (adopted from my code at:
# https://github.com/clairecwinfrey/PhanBioMS_scripts/blob/master/R_scripts/figures/taxonomic_bartrimmedPlots.R)
# TURN ASVs INTO PHYLUM LEVEL
rarefiedTrimmed.phylum.glom <-  tax_glom(TrimmedSRS_16S.ps, taxrank = "Phylum") 
tax_table(rarefiedTrimmed.phylum.glom) # good, this is only phyla (42 different phyla!)

# TRANSFORM SAMPLE COUNTS ON JUST GLOMMED SAMPLES (UNLIKE WHAT WE DID AT FIRST)
relabunTrimmed.phyla.0 <- transform_sample_counts(rarefiedTrimmed.phylum.glom, function(x) x / sum(x) )
rownames(otu_table(relabunTrimmed.phyla.0)) #I think that the ASVs are just representative from each phylum

# MERGE SAMPLES so that we only have combined abundances for site and different kinds of controls
relabunTrimmed.phyla.1 <- merge_samples(relabunTrimmed.phyla.0, group = "EU")
sample_data(relabunTrimmed.phyla.1) 

# CONVERT TO PROPORTIONS AGAIN B/C TOTAL ABUNDANCE OF EACH SITE WILL EQUAL NUMBER OF SPECIES THAT WERE MERGED
relabunTrimmed.phyla.2 <- transform_sample_counts(relabunTrimmed.phyla.1, function(x) x / sum(x))
sample_data(relabunTrimmed.phyla.2)

# NOW GET ONLY TAXA THAT COMPRISE AT LEAST 1% OF THE ABUNDANCE 
relabunTrimmed.phyla.df <-psmelt(relabunTrimmed.phyla.2)
dim(relabunTrimmed.phyla.df) #
sum(relabunTrimmed.phyla.df[,3]) 
colnames(relabunTrimmed.phyla.df) 
relabunTrimmed.phylatrimmedTop99 <- relabunTrimmed.phyla.df
relabunTrimmed.phylatrimmedTop99.5 <- relabunTrimmed.phyla.df
relabunTrimmed.phylatrimmedTop95 <- relabunTrimmed.phyla.df

relabunTrimmed.phylatrimmedTop99$Phylum[relabunTrimmed.phylatrimmedTop99$Abundance < 0.01] <- "< 1% abund."
relabunTrimmed.phylatrimmedTop99.5$Phylum[relabunTrimmed.phylatrimmedTop99.5$Abundance < 0.005] <- "< .5% abund."
relabunTrimmed.phylatrimmedTop95$Phylum[relabunTrimmed.phylatrimmedTop95$Abundance < 0.05] <- "< 5% abund."

trimmedTop_99p_phyla <- unique(relabunTrimmed.phylatrimmedTop99$Phylum)
trimmedTop_99p_phyla

trimmedTop_99.5p_phyla <- unique(relabunTrimmed.phylatrimmedTop99.5$Phylum)
trimmedTop_99.5p_phyla

# Surprised that Gemmatimonadetes not above >1%! But Gemmatimonadetes in trimmedTop 16 and
# comprises at least 0.5% of the total abundance

# Phyla comprising at least 0.5% of total abundance
#quartz()
phylumtrimmedPlot99.5percent <- ggplot(data=relabunTrimmed.phylatrimmedTop99.5, aes(x=Sample, y=Abundance, fill=Phylum)) + theme(axis.title.y = element_text(size = 14, face = "bold")) + theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(colour = "black", size = 12, face = "bold"))
phylumtrimmedPlot99.5percent + geom_bar(aes(), stat="identity", position="fill") +
  theme(legend.position="bottom") +
  guides(fill=guide_legend(nrow=4)) + theme(legend.text = element_text(colour="black", size = 10))  + ggtitle("Phyla comprising at least 0.5% of total abundance")

# "Phyla comprising at least 1% of total abundance"
#quartz()
phylumtrimmedPlot.99percent <- ggplot(data=relabunTrimmed.phylatrimmedTop99, aes(x=Sample, y=Abundance, fill=Phylum)) + theme(axis.title.y = element_text(size = 14, face = "bold")) + theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(colour = "black", size = 12, face = "bold"))
phylumtrimmedPlot.99percent <- phylumtrimmedPlot.99percent + geom_bar(aes(), stat="identity", position="fill") +
  theme(legend.position="bottom") +
  guides(fill=guide_legend(nrow=4)) + theme(legend.text = element_text(colour="black", size = 10))  + ggtitle("Phyla comprising at least 1% of total abundance")
phylumtrimmedPlot.99percent

# Get exact abundances of each phyla (trimmedTop 99.5%):
colnames(relabunTrimmed.phylatrimmedTop99.5)
relabunTrimmed.phylatrimmedTop99.5[,2] #This is EU... now "Sample" because of the glomming!

trimmedTop_99.5p_phyla <- relabunTrimmed.phylatrimmedTop99.5 %>%
  group_by(Sample, Phylum) %>%
  summarize(Mean = mean(Abundance)) %>%
  arrange(-Mean) 

#####################
# trimmedTop CLASSES:
#####################

# (adopted from my code at:
# https://github.com/clairecwinfrey/PhanBioMS_scripts/blob/master/R_scripts/figures/taxonomic_bartrimmedPlots.R)
# TURN ASVs INTO PHYLUM LEVEL
rarefiedTrimmed.class.glom <-  tax_glom(TrimmedSRS_16S.ps, taxrank = "Class") 
tax_table(rarefiedTrimmed.class.glom) # good, this is only class 

# TRANSFORM SAMPLE COUNTS ON JUST GLOMMED SAMPLES (UNLIKE WHAT WE DID AT FIRST)
relabunTrimmed.class.0 <- transform_sample_counts(rarefiedTrimmed.class.glom, function(x) x / sum(x) )
rownames(otu_table(relabunTrimmed.class.0)) # I think that the ASVs are just representative from each class

# MERGE SAMPLES so that we only have combined abundances for site and different kinds of controls
relabunTrimmed.class.1 <- merge_samples(relabunTrimmed.class.0, group = "EU")
sample_data(relabunTrimmed.class.1) #shows that we still have samples from each EU, biocrust, and extcontrol (water)

# CONVERT TO PROPORTIONS AGAIN B/C TOTAL ABUNDANCE OF EACH SITE WILL EQUAL NUMBER OF SPECIES THAT WERE MERGED
relabunTrimmed.class.2 <- transform_sample_counts(relabunTrimmed.class.1, function(x) x / sum(x))
sample_data(relabunTrimmed.class.2)

# NOW GET ONLY TAXA THAT COMPRISE AT LEAST 1% OF THE ABUNDANCE 
relabunTrimmed.class.df <-psmelt(relabunTrimmed.class.2)
dim(relabunTrimmed.class.df) #
sum(relabunTrimmed.class.df[,3]) 
colnames(relabunTrimmed.class.df) 
relabunTrimmed.classtrimmedTop99 <- relabunTrimmed.class.df
relabunTrimmed.classtrimmedTop95 <- relabunTrimmed.class.df

relabunTrimmed.classtrimmedTop99$Class[relabunTrimmed.classtrimmedTop99$Abundance < 0.01] <- "< 1% abund."
relabunTrimmed.classtrimmedTop95$Class[relabunTrimmed.classtrimmedTop95$Abundance < 0.05] <- "< 5% abund."

trimmedTop_99p_class <- unique(relabunTrimmed.classtrimmedTop99$Class)
trimmedTop_99p_class

trimmedTop_95p_class <- unique(relabunTrimmed.classtrimmedTop95$Class)
trimmedTop_95p_class

# Phyla comprising at least 5% of total abundance
#quartz()
classtrimmedPlot.95pt <- ggplot(data=relabunTrimmed.classtrimmedTop95, aes(x=Sample, y=Abundance, fill=Class)) + theme(axis.title.y = element_text(size = 14, face = "bold")) + theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(colour = "black", size = 12, face = "bold"))
classtrimmedPlot.95pt + geom_bar(aes(), stat="identity", position="fill") +
  theme(legend.position="bottom") +
  guides(fill=guide_legend(nrow=4)) + theme(legend.text = element_text(colour="black", size = 5.5))  + ggtitle("Classes comprising at least 5% of total abundance")

# Make new phyloseq object by removing biocrust and controls before ordination:
trimmedJustsoils.ps <- phyloseq::subset_samples(TrimmedSRS_16S.ps, Type != "BioCrust" & Type != "ExtContWater")
unique(sample_data(trimmedJustsoils.ps)[,27]) #only soils
sample_names(trimmedJustsoils.ps) #another check to show that we have only soils!

# Plot ordination
# Bray-Curtis dissimilarities based on square-root transformed data
set.seed(19)
trimOrd <- ordinate(trimmedJustsoils.ps, method = "NMDS", distance = "bray", trymax = 100) 

# With no black outline around points
trimmedBrayNMDS <- phyloseq::plot_ordination(trimmedJustsoils.ps, trimOrd, type= "samples", color= "EU")
#quartz()
trimmedBrayNMDS + geom_polygon(aes(fill=EU)) + geom_point(size=3) + ggtitle("NMDS based on Bray-Curtis Dissimilarities (Trimmed)")

# Add in labels to figure out what weird samples are
labeledTrimmedBrayNMDS <- phyloseq::plot_ordination(trimmedJustsoils.ps, trimOrd, type= "samples", color= "EU", label = "Sample.ID")
#quartz()
labeledTrimmedBrayNMDS + geom_polygon(aes(fill=EU)) + geom_point(size=3) + ggtitle("NMDS based on Bray-Curtis Dissimilarities")


##################################################################################
# III. EXPLORATION OF "OUTLIER" TAXA-- taxonomic barcharts, differential abundance analysis
##################################################################################

# Visible outliers: (going clockwise from the top left on the ordination plot "labeledTrimmedBrayNMDS" above):
# 53ND_B_20, 10C_L_10, 53ND_R_40, 53SD_R_10, 52D_R_10, 52D_L_100, (Maybe!) 52D_R_90, 53ND_L_80
# 10C_R_60 is the one with a ton of samples!
outliersTrimmed <- c("53ND_B_20", "10C_L_10", "53ND_R_40", "53SD_R_10", "52D_R_10", "52D_L_100", "53ND_L_80", "10C_R_60")
outliersTrimmed <- sort(outliersTrimmed) #sorted this because next lines of code will automatically sort
outliersTrimmed

# 10C_R_60 is the one with a ton of samples!

# Do these outliers look weird?
# Get the rows for the outliers in sample_df
#View(samples_df)]
# Code below looks for the row numbers corresponding to each one of the outlier samples. It automatically sorts the data
outlier_index <- which(samples_df$Sample.ID=="53ND_B_20" | samples_df$Sample.ID=="10C_L_10" | samples_df$Sample.ID=="53ND_R_40"
                       | samples_df$Sample.ID=="53SD_R_10" | samples_df$Sample.ID=="52D_R_10" | samples_df$Sample.ID=="52D_L_100"
                       | samples_df$Sample.ID=="53ND_L_80" | samples_df$Sample.ID=="10C_R_60")
outlier_index 
samples_df[outlier_index,1]==outliersTrimmed # Names using outlier_index match the outliers defined above!
outlierSampNumb <- rownames(samples_df[outlier_index,]) # now these are the actual names used in bioinformatics and in my plate schemes.
# They differ from outlier_index because the row numbers are not the same as the row names.
# Can use "outlierSampNumb" with phyloseq's "prune_samples" function

# Just looking at ASVs
# Now, that we have the numbers corresponding to each sample, we can isolate the ASVs
outlierSamples <- get_taxa(trimmedJustsoils.ps, outlierSampNumb)
colnames(outlierSamples) == outlierSampNumb #Nice, this worked
str(outlierSamples)

# ALTERNATIVELY, CAN JUST USE PRUNE_SAMPLES IN PHYLOSEQ
sample_names(trimmedJustsoils.ps) #sample names are the numbers that we used in the PCR plate organization!
outlierSampNumb
# Also, these are the same as those in outlier_index
outliers.ps <- prune_samples(outlierSampNumb, trimmedJustsoils.ps)
sort(sample_data(outliers.ps)$Sample.ID) == sort(outliersTrimmed) #these are equal to the outliers we defined above!

outliersTax <- taxtable_outta_ps(outliers.ps)
outliersASVs <- ASVs_outta_ps(outliers.ps)
# View(outliersTax)

# Do the ASVs and the community abundances look weird for the outliers?
# Make an mctoolsr like table to view if we want!

outliersASVtax <- cbind.data.frame(outliersASVs, outliersTax)
#View(outliersASVtax)

# outliersASVtax2 <- merge(outliersASVs, outliersTax, by="row.names",all.x=TRUE) # this would have worked as well, 
# but then the rownames would have been a new column, not still a rowname!
outliersASVtax$Abundance <- rowSums(outliersASVtax[,1:7]) #What is the abundance of each ASV across all 7 samples? 
# Interestingly, the mycobacterium taxon that is the third most abundant across all samples has moved down a lot!

# What do the top phyla and classes look like?

# Top Phyla
# TURN ASVs INTO PHYLUM LEVEL
outliers.phylum.glom <-  tax_glom(outliers.ps, taxrank = "Phylum") 
tax_table(outliers.phylum.glom) # only phyla
sample_data(outliers.phylum.glom)

# TRANSFORM SAMPLE COUNTS ON JUST GLOMMED SAMPLES (UNLIKE WHAT WE DID AT FIRST)
outliers.phylum.0 <- transform_sample_counts(outliers.phylum.glom, function(x) x / sum(x) )
rownames(otu_table(outliers.phylum.0)) # I think that the ASVs are just representative from each class
sample_data(outliers.phylum.0)

# NOW GET ONLY TAXA THAT COMPRISE AT LEAST 1% OF THE ABUNDANCE 
outliers.phylum.df <-psmelt(outliers.phylum.0)
dim(outliers.phylum.df ) # 231 for the number of phyla
sum(outliers.phylum.df[,3]) #nice this is 7, the number of samples!
colnames(outliers.phylum.df) 
outliers.phylumTop99 <- outliers.phylum.df
outliers.phylumTop95 <- outliers.phylum.df

outliers.phylumTop99$Phylum[outliers.phylumTop99$Abundance < 0.01] <- "< 1% abund."
outliers.phylumTop95$Phylum[outliers.phylumTop95$Abundance < 0.05] <- "< 5% abund."

outliers.phylumTop_99p <- unique(outliers.phylumTop99$Phylum)
outliers.phylumTop_99p #"Acidobacteria" ,"Proteobacteria","Verrucomicrobia","Chloroflexi" ,"Planctomycetes","Actinobacteria","Firmicutes","Bacteroidetes","WPS-2", "Cyanobacteria", "Armatimonadetes", Rokubacteria" 

outliers.phylumTop_95p <- unique(outliers.phylumTop95$Phylum)
outliers.phylumTop_95p #"Acidobacteria"   "Proteobacteria"  "Verrucomicrobia" "Chloroflexi"     "Planctomycetes"  "Actinobacteria"  "Firmicutes"      "Bacteroidetes"

# Phyla comprising at least 5% of total abundance
outliers.phylumPlot.95pt <- ggplot(data=outliers.phylumTop95, aes(x=Sample, y=Abundance, fill=Phylum)) + theme(axis.title.y = element_text(size = 14, face = "bold")) + theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(colour = "black", size = 12, face = "bold"))
outliers.phylumPlot.95pt <- outliers.phylumPlot.95pt + geom_bar(aes(), stat="identity", position="fill") +
  scale_fill_manual(values = c("#999999", "#f781bf", "#a65628", "#ffff33", "#ff7f00", "#984ea3", "#4daf4a", "#377eb8", "#e41a1c")) +
  theme(legend.position="bottom") +
  guides(fill=guide_legend(nrow=4)) + theme(legend.text = element_text(colour="black", size = 6))  + ggtitle("Outliers: Phyla at least 5% of total abundance")
#quartz()
outliers.phylumPlot.95pt

# Compare this with soils aggregated by site/EU (trimmed, so aftr removing rare taxa):
# First, get the top phyla of just soil

# TURN ASVs INTO PHYLUM LEVEL
justsoilsphylum.glom <-  tax_glom(trimmedJustsoils.ps, taxrank = "Phylum") 
tax_table(justsoilsphylum.glom)

# TRANSFORM SAMPLE COUNTS ON JUST GLOMMED SAMPLES (UNLIKE WHAT WE DID AT FIRST)
justsoils.phyla.0 <- transform_sample_counts(justsoilsphylum.glom, function(x) x / sum(x) )
rownames(otu_table(justsoils.phyla.0))

# MERGE SAMPLES so that we only have combined abundances for site and different kinds of controls
justsoils.phyla.1 <- merge_samples(justsoils.phyla.0, group = "EU")
sample_data(justsoils.phyla.1) 

# CONVERT TO PROPORTIONS AGAIN B/C TOTAL ABUNDANCE OF EACH SITE WILL EQUAL NUMBER OF SPECIES THAT WERE MERGED
justsoils.phyla.2 <- transform_sample_counts(justsoils.phyla.1, function(x) x / sum(x))
sample_data(justsoils.phyla.2)

# NOW GET ONLY TAXA THAT COMPRISE AT LEAST 1% OF THE ABUNDANCE 
justsoils.phyla.df <-psmelt(justsoils.phyla.2)
dim(justsoils.phyla.df) 
justsoils.phyla.Top95 <- justsoils.phyla.df
justsoils.phyla.Top99 <- justsoils.phyla.df

justsoils.phyla.Top95$Phylum[justsoils.phyla.Top95$Abundance < 0.05] <- "< 5% abund."
justsoils.phyla.Top99$Phylum[justsoils.phyla.Top99$Abundance < 0.01] <- "< 1% abund."

# Plot 
# Phyla comprising at least 1% of total abundance 
# Colors are different for easier comparison with outliers in next section
justsoils.phylaPlot.95percent <- ggplot(data=justsoils.phyla.Top95, aes(x=Sample, y=Abundance, fill=Phylum)) + theme(axis.title.y = element_text(size = 14, face = "bold")) + theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(colour = "black", size = 12, face = "bold"))
justsoils.phylaPlot.95percent <- justsoils.phylaPlot.95percent + geom_bar(aes(), stat="identity", position="fill") + 
  scale_fill_manual(values = c("#999999", "#f781bf", "#a65628", "#ff7f00", "#4daf4a", "#377eb8", "#e41a1c")) +
  theme(legend.position="bottom") +
  guides(fill=guide_legend(nrow=4)) + theme(legend.text = element_text(colour="black", size = 6))  + ggtitle("Rarefied Soils: Phyla at least 5% of total abundance")
#quartz()
justsoils.phylaPlot.95percent 

#quartz()
grid.arrange(outliers.phylumPlot.95pt, justsoils.phylaPlot.95percent, nrow=1)
# You can see here that overall the outliers look pretty similar to the regular soil samples, with the exception of 
# Sample 11, which has a lot of Bacteroides and Firmicutes. 
# # Looking at 11 in the mctoolsR like file (View(outliersASVtax)), you can see that some of the top taxa are very different (with the
# exception of the Bradyrhizobium and another ASV in the same family).
# It's top ASV is in Firmicutes, and other Firmicutes are high on the list too

# 11  was sample 53ND_L_80, it was in Plate 1, well B1

#################
# Differential abundance analysis of outliers versus total data
# ***IMPORTANT NOTE****: THE DIFFERENTIAL ABUNDANCE ANALYSIS PERFORMED HERE IS *NOT* THE ONE
# USED FOR LATER ANALYSES OF THE EDGE EFFECT TAXA.
#################
# Paper for DESeq2 : Love, M.I., Huber, W. & Anders, S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. 
# Genome Biol 15, 550 (2014). https://doi.org/10.1186/s13059-014-0550-8
# This section works off example here: https://joey711.github.io/phyloseq-extensions/DESeq2.html

# First need to add a new metadata variable to trimmedJustsoils.ps so that we can compare
# outliers and rest of data
sample_data(trimmedJustsoils.ps)$Outlier <- rep("not_outlier", nrow(sample_data(trimmedJustsoils.ps)))
sample_data(trimmedJustsoils.ps)$Outlier[c(3, 12, 25, 37, 60, 221, 233)] <- "outlier"
sort(sample_data(trimmedJustsoils.ps)$Sample.ID[c(3, 12, 25, 37, 60, 221, 233)]) == sort(outliersTrimmed) #names are correct

outlierDeseq1 <- phyloseq_to_deseq2(trimmedJustsoils.ps, ~ Outlier)
# Note: uses default Benjamini-Hochberg correction 
outlierDeseqtested <- DESeq(outlierDeseq1, test="Wald", fitType = "parametric")

outlierDeSeq_res <- results(outlierDeseqtested, cooksCutoff = FALSE)
alpha <- 0.01 #note that this is less sensitive than 0.001 used for forests versus soils 
outlierSigtab <- outlierDeSeq_res[which(outlierDeSeq_res$padj < alpha), ]
outlierSigtab <- cbind(as(outlierSigtab, "data.frame"), as(tax_table(trimmedJustsoils.ps)[rownames(outlierSigtab), ], "matrix"))
head(outlierSigtab)
dim(outlierSigtab) 
#View(outlierSigtab)

# make new phyloseq object where samples are averaged together (ASV counts are summed) across Habitat type (here just forest and patch)
outlier.ps <- merge_samples(trimmedJustsoils.ps, "Outlier")

outlierASVs <- ASVs_outta_ps(outlier.ps)
outlierASVs <- t(outlierASVs) # make ASVs rows and outlier versus no outlier columns

# Merge dataframes so that abundance in outlier and non outlier is present, then rename columns 
outlierSigtab <- left_join(rownames_to_column(outlierSigtab), rownames_to_column(as.data.frame(outlierASVs)), by="rowname")
# View(outlierSigtab)
colnames(outlierSigtab)[14] <- "notOutlierAbundance"
colnames(outlierSigtab)[15] <- "outlierAbundance"
# Add in mean abundance in each category
outlierSigtab$meanPercentNotOutlier <- outlierSigtab[14]/sum(outlierASVs[,1])*100 #this is divided by total number of counts in the forest samples (after rarefying)
outlierSigtab$meanPercentOutlier <- outlierSigtab[15]/sum(outlierASVs[,2])*100 #this is divided by total number of counts in the patch samples (after rarefying)
#View(outlierSigtab)

# Plot 
#quartz()
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(outlierSigtab$log2FoldChange, outlierSigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
outlierSigtab$Phylum = factor(as.character(outlierSigtab$Phylum), levels=names(x))
# Genus order
x = tapply(outlierSigtab$log2FoldChange, outlierSigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
outlierSigtab$Genus = factor(as.character(outlierSigtab$Genus), levels=names(x))

# point size does not vary
#quartz()
ggplot(outlierSigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + ggtitle("Log2Fold Change between Outliers and non-outliers") + ylab("Log2FoldChange (Relative to Outliers)")

# point size varies based on baseMean (i.e. average of the normalized count values, dividing by size factors, taken over all samples)
#quartz()
ggplot(outlierSigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(aes(size = baseMean)) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + ggtitle("Log2Fold Change between Outliers and non-outliers") + ylab("Log2FoldChange (Relative to Outliers)")

# The fact that there are similar phyla and genera that are both higher and lower across these two groups suggests that these outliers might not
# be so weird. The tail of taxa (near the right hand side of the plot) is unsurprising; this just says that there are wuite a few taxa that
# are in all the rest of the samples but are not found in the outliers (which is far fewer samples)


##################################################################################
# IV. COMPARING FOREST AND PATCHES WITH ORDINATIONS AND DIFFERENTIAL ABUNDANCE ANALYSIS
##################################################################################

#################
# Ordination of forest versus patch
#################

# Added "habitat" column to metadata in Excel, re-ran whole script
#quartz()
HabitatBrayNMDS <- phyloseq::plot_ordination(trimmedJustsoils.ps, trimOrd, type= "samples", color= "Habitat")
HabitatBrayNMDS <- HabitatBrayNMDS + geom_polygon(aes(fill=Habitat)) + geom_point(size=3) + ggtitle("NMDS based on Bray-Curtis Dissimilarities")
HabitatBrayNMDS 
# Cool, you can see that the forest and the patch separate out, with edge somewhat in between!

#################
# Differential abundance analysis on forest v. patch
#################
sample_data(trimmedJustsoils.ps)$Habitat
# Remove edge for now, since I'm not sure how it fits into the dichotomy of forest v soil and it'll complicate
# differential abundance analysis 
soil_noedge.ps <- subset_samples(trimmedJustsoils.ps, Habitat != "edge") #210 samples remaining
colnames(sample_data(soil_noedge.ps))

# Note, if error with "Rcpp" package version here, but version is okay, re-install Rcpp
# and then restart R.
soilDeseq1 <- phyloseq_to_deseq2(soil_noedge.ps, ~ Habitat)
# Note: uses default Benjamini-Hochberg correction 
set.seed(19)
soilDeseqtested <- DESeq(soilDeseq1, test="Wald", fitType = "parametric")

# Good explanation of results is here: https://support.illumina.com/help/BS_App_RNASeq_DE_OLH_1000000071939/Content/Source/Informatics/Apps/DESeq2ResultFile_swBS.htm#:~:text=baseMean%E2%80%94The%20average%20of%20the,factors%2C%20taken%20over%20all%20samples.&text=log2FoldChange%E2%80%93The%20effect%20size%20estimate,the%20comparison%20and%20control%20groups.
DeSeq_res <- results(soilDeseqtested, cooksCutoff = FALSE)
alpha <- 0.001
sigtab <- DeSeq_res[which(DeSeq_res$padj < alpha), ]
sigtab <- cbind(as(sigtab, "data.frame"), as(tax_table(soil_noedge.ps)[rownames(sigtab), ], "matrix"))
head(sigtab)
dim(sigtab) #1,846 ASVs out of the 10,257 that we had, had a p value less than 0.001, forests and meadows are super different!
# BaseMean is the The average of the normalized count values, dividing by size factors, taken over all samples
sigtab <- rownames_to_column(sigtab) #make these rownames a column so that the merge below with habitatASVs can be done

# Add in abundance for forest or patch
# How many forest and patch samples?
# Forest Samples
length(which(sample_data(soil_noedge.ps)$Habitat=="forest")) #116
# Patch Samples
length(which(sample_data(soil_noedge.ps)$Habitat=="patch")) #94

# make new phyloseq object where samples are averaged together (ASV counts are summed) across Habitat type (here just forest and patch)
habitat.ps <- merge_samples(soil_noedge.ps, "Habitat") #10137 ASVs and 2 samples

habitatASVs <- ASVs_outta_ps(habitat.ps)
dim(habitatASVs) #rows are samples (i.e. forest and patch), columns are each ASV
habitatASVs <- t(habitatASVs) #now ASVs are rows
# class(habitatASVs)

# Merge dataframes so that abundance in forest and patch is present, then rename columns 
sigtab <- left_join(sigtab, rownames_to_column(as.data.frame(habitatASVs)), by="rowname")
#View(sigtab)
colnames(sigtab)[14] <- "forestAbundance"
colnames(sigtab)[15] <- "patchAbundance"
# Add in mean abundance in each category
sigtab$meanPercentForest <- sigtab[14]/sum(habitatASVs[,1])*100 #this is divided by total number of counts in the forest samples (after rarefying)
sigtab$meanPercentPatch <- sigtab[15]/sum(habitatASVs[,2])*100 #this is divided by total number of counts in the patch samples (after rarefying)
#View(sigtab)

# Plot 
#quartz()
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))

# point size does not vary
#quartz()
ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + ggtitle("Log2Fold Change between Forest and Patch Soils-Genera") + ylab("Log2FoldChange (Relative to Patch)")

# point size varies based on baseMean (i.e. average of the normalized count values, dividing by size factors, taken over all samples)
#quartz()
ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(aes(size = baseMean)) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + ggtitle("Log2Fold Change between Forest and Patch Soils-Genera") + ylab("Log2FoldChange (Relative to Patch)")

# Family instead
# Family order
# Family order
x2 = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x2) max(x2))
x2 = sort(x2, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x2))
# Family order
x2 = tapply(sigtab$log2FoldChange, sigtab$Family, function(x2) max(x2))
x2 = sort(x2, TRUE)
sigtab$Family = factor(as.character(sigtab$Family), levels=names(x2))

# point size varies based on baseMean (i.e. average of the normalized count values, dividing by size factors, taken over all samples)
#quartz()
ggplot(sigtab, aes(x=Family, y=log2FoldChange, color=Phylum)) + geom_point(aes(size = baseMean)) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + ggtitle("Log2Fold Change between Forest and Patch Soils-Families") + ylab("Log2FoldChange (Relative to Patch)")

########################
# Does the Bray-Curtis distance b/w samples tend to increase with distance?
########################
ASVtab <- psotu2veg(trimmedJustsoils.ps)
#View(ASVtab)

ASVbrayDist <- vegdist(ASVtab, method = "bray")
ASVBrayDist.mat <- as.matrix(ASVbrayDist)
diag(ASVBrayDist.mat) <- NA

# Make data frame for indexing matrix
rownames(ASVBrayDist.mat) == rownames(sample_data(trimmedJustsoils.ps))
# Since the thing above is true, we can get meter and replace the columns and rownames in the BC dist mat with it!
samps <- colnames(ASVBrayDist.mat)
meter <- sample_data(trimmedJustsoils.ps)$Meter
index.df <- data.frame(samps, meter)

# Get row and column numbers of different meters
m10 <- which(index.df$meter == 10)
m20 <- which(index.df$meter == 20)
m30<- which(index.df$meter == 30)
m40 <- which(index.df$meter == 40)
m50 <- which(index.df$meter == 50)
m60 <- which(index.df$meter == 60)
m70 <- which(index.df$meter == 70)
m80<- which(index.df$meter == 80)
m90 <- which(index.df$meter == 90)
m100 <- which(index.df$meter == 100)

# Distances between patch at 10m and each point along transect
m10_comp10 <- ASVBrayDist.mat[m10, m10]
m10_comp20 <- ASVBrayDist.mat[m10, m20]
m10_comp30 <- ASVBrayDist.mat[m10, m30]
m10_comp40 <- ASVBrayDist.mat[m10, m40]
m10_comp50 <- ASVBrayDist.mat[m10, m50]
m10_comp60 <- ASVBrayDist.mat[m10, m60]
m10_comp70 <- ASVBrayDist.mat[m10, m70]
m10_comp80 <- ASVBrayDist.mat[m10, m80]
m10_comp90 <- ASVBrayDist.mat[m10, m90]
m10_comp100 <- ASVBrayDist.mat[m10, m100]

# Distances between forest at 100m and each point along transect
m100_comp100 <- ASVBrayDist.mat[m100, m100]
m100_comp90 <- ASVBrayDist.mat[m100, m90]
m100_comp80 <- ASVBrayDist.mat[m100, m80]
m100_comp70 <- ASVBrayDist.mat[m100, m70]
m100_comp60 <- ASVBrayDist.mat[m100, m60]
m100_comp50 <- ASVBrayDist.mat[m100, m50]
m100_comp40 <- ASVBrayDist.mat[m100, m40]
m100_comp30 <- ASVBrayDist.mat[m100, m30]
m100_comp20 <- ASVBrayDist.mat[m100, m20]
m100_comp10 <- ASVBrayDist.mat[m100, m10]

######## MAKE BOXPLOTS ######## 
bold_a <- expression(bold("Dissimilarity Relative to 10m (patch)"))
bold_b <- expression(bold("Dissimilarity Relative to 100m (forest)"))
#quartz()
par(mfrow=c(1,2))
patch_box <- boxplot(list(m10_comp10, m10_comp20, m10_comp30, m10_comp40,
                          m10_comp50, m10_comp60, m10_comp70,
                          m10_comp80, m10_comp90, m10_comp100),
                     ylab = "Bray-Curtis Dissimilarity",
                     names = c("10", "20", "30", "40", "50",
                               "60", "70", "80", "90", "100"), cex.axis = 0.8,
                     xlab = "meter along transect",
                     cex.lab = 1,
                     ylim=c(0.0, 1.0))
mtext(text=bold_a, side=3, adj = -0.065, line = 2)
forest_box <- boxplot(list(m100_comp10, m100_comp20, m100_comp30, m100_comp40,
                           m100_comp50, m100_comp60, m100_comp70,
                           m100_comp80, m100_comp90, m100_comp100), 
                      ylab = "Bray-Curtis Dissimilarity",
                      names = c("10", "20", "30", "40", "50",
                                "60", "70", "80", "90", "100"), cex.axis = 0.8,
                      xlab = "meter along transect",
                      cex.lab = 1,
                      ylim=c(0.0, 1.0))
mtext(text=bold_b, side=3, adj = -0.065, line = 2)

##################################
# SAVE ALL OF THESE FOR EASY ACCESS
##################################

# all resaved Sept 29, 2022:
#save(rarefied.ps, samples_df, ASVsTrimmed, taxTrimmed, trimmedJustsoils.ps, trimOrd, outliersTrimmed, outliersASVtax, file = "RobjectsSaved/EDA16SAug2021")
#save(trimmedJustsoils.ps, file= "RobjectsSaved/trimmedJustSoils.ps") 
#soils_noEuks <- ASVs_outta_ps(noeuksorNAs_ps) #this is the dataset pre-rarefaction, but after removing eukaryotes and ASVs not assigned at least to phylum.
#soils_noEuksTaxTable <- taxtable_outta_ps(noeuksorNAs_ps)
#write.csv(soils_noEuks, file = "RobjectsSaved/SRSsoilsNoEuks.csv") 

# last saved April 25, 2022:
#write.csv(soils_noEuksTaxTable, file= "RobjectsSaved/SRSsoilsNoEuksTaxTable.csv") #saved April 25, 2022 and sent to Josep R.
#write.csv(samples_df, file= "RobjectsSaved/SRS_soilsFinalMetadata.csv") #saved April 25, 2022 and sent to Josep R.