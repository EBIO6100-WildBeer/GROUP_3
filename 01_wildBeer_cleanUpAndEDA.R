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
# View(beerMeta)

# Libraries
library(vegan) #version 2.5-7
library(phyloseq)
library(tidyverse)

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
# View(metaDat)

# 6. PUT TAXONOMY, ASV TABLE, AND METADATA TOGETHER TO MAKE PHYLOSEQ OBJECT
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




