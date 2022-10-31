# Wild Beer Metagenomics- Data Clean Up and Exploratory Data Analysis
# Begun October 24, 2022
########################################

########## SET UP ##########

# Read in files and check them out (from github repo on local machine)
beerComm <- read.csv(file="~/Desktop/CU_Research/wildBeer/project_data_and_tutorials/data/2022_10_20_WildBeer.phyloFlash.extractedSSUclassifications.csv") #ASV table (sort of)
beerMeta <- read.delim(file="~/Desktop/CU_Research/wildBeer/project_data_and_tutorials/data/2022_10_12_WildBeer_metadata.tsv", sep="\t", header=TRUE) #metadata
head(beerComm)
head(beerMeta)
#View(beerComm)

# Libraries
library(vegan) #version 2.5-7
library(phyloseq)
library(tidyverse)

########################################

########## MAIN PART OF SCRIPT ##########

# 1. Merge community data and metadata to make one file
#olnames(beerComm)
#colnames(beerMeta)
# Rename "sample" in beerComm to be "Sample.ID" to match beerMeta before merging based on this variable
#colnames(beerComm)[10] <- "Sample.ID"
# How do these dataframes overlap?
#sharedIDs <- intersect(unique(beerComm$Sample.ID), unique(beerMeta$Sample.ID))
#unique(beerComm$Sample.ID)[which(unique(beerComm$Sample.ID) %in% sharedIDs== FALSE)] #samples that occur in beerComm but not in beerMeta
#unique(beerMeta$Sample.ID)[which(unique(beerMeta$Sample.ID) %in% sharedIDs== FALSE)] #unhopped outdoors 3-1 is in beerMeta but not in beerComm
# since we lost the unhopped outdoor sample in the community data, we'll drop it from analysis 

#beerAllData <- merge(beerComm, beerMeta, by= "Sample.ID", all.x = TRUE)
#View(beerAllData)

# Taxonomy
#unique(beerComm$dbHit) #only 26 unique ASVs!!

# 1. CLEAN UP TAXONOMY 
colnames(beerComm)
head(beerComm)
step1 <- beerComm[,c(5,6)] #get only dbhit and taxonomy columns; keep all rows
colnames(step1)
# Right now, the fungal and bacteria taxonomy is not ordered the same. For example, the fungal taxonomic information
# starts at Eukaryota, whereas the bacterial starts at bacteria (not prokaryote). And there is a whole lot of 
# extra taxonomic information between "eukaryote" and "fungi" in the fungal ones.
step1[6,] #for example, here

# Remove extra taxonomic information between Eukaryote and Fungi, as well as subkingdom and subdivisions
step2 <- step1[,2] %>% 
  stringr::str_replace("^.*Nucletmycea;", "") %>% 
  stringr::str_replace("Dikarya;", "") %>% 
  stringr::str_replace("Saccharomycotina;", "") %>% 
  stringr::str_replace("Agaricomycotina;", "")
colnames(step2)
cleaned1 <- cbind(step1[,1], step2) #combine the cleaned taxonomy column with dbhit from step 1 
colnames(cleaned1) <- colnames(step1) #make dbHit and taxonomy the column names again
# View(cleaned1)

# divide out taxonomy into different columns by level
beerTax <- tidyr::separate(as.data.frame(cleaned1), col = taxonomy, into= c("Kingdom", "Phylum", "Class", 
                                                                  "Order", "Family", "Genus", "Species"),
                    sep = ";")

# View taxonomy file
# View(beerTax) 

colnames(beerComm)
unique(beerTax$dbHit == beerComm$dbHit) # since the order and the names of the dbHits are the same for beerTax and beerCommLess, 
beerCommLess <- beerComm[,c(1:4,7:10)] #no need to keep dbhit or Taxonomy columns (since we fixed taxonomy)
dim(beerCommLess)
# we can just do a simple cbind to merge them
beerCommCleaned <- cbind(beerCommLess, beerTax)
# View(beerCommCleaned)

# Make final table (As Noah asked for in class)
commCompTable <- beerCommCleaned %>% 
  select(dbHit, read_cov, sample, Species) %>% 
  pivot_wider(names_from=sample, values_from=read_cov) %>% 
  mutate_all(~replace(., is.na(.), 0))
# View(commCompTable)

################################
# 2. MAKE PHYLOSEQ OBJECT -Claire will do sometime before Monday :)
################################
# Make phyloseq object: Phyloseq wants 3-4 inputs: (1) a "taxonomy file" that has all the OTU names paired
# with their taxonomy (i.e. Kingdom, Phylum, Class, Order, etc...), (2) an OTU table that has these OTUs with
# information about how often these OTUs in the taxonomy file appear in each sample, and (3): metadata file, 
# which links each sample to their metadata
