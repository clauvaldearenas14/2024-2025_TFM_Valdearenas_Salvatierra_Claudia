##################################################################
##### CREATE A PHYLOSEQ OBJET WITH ALL MICROBIOTA INFORMATION ####
##################################################################

#### INSTALLING PACKAGES ####

#install.packages("vegan")    
#install.packages("ape")
#install.packages("ggplot2")


#### LOAD LIBRARIES ####

library(vegan)
library(ape)
library(ggplot2)

# Set the working directory

#setwd("~/xxx") #Navigate to the "Inputs" directory


#### IMPORT DATA INTO R ####

otumat = read.delim("CICHLIDS_SILVA/otumat_SILVA.txt", row.names = 1) #the table of microbial abundances 
taxmat = read.delim("CICHLIDS_SILVA/taxmat_SILVA.txt", row.names = 1) #the table with taxonomy information
#df <- replace(taxmat, taxmat=='', NA)
sampledata = read.delim("CICHLIDS_SILVA/Metadata_G1.txt", row.names = 1) # the table with sample metadata 

#Filter metadata to match otumat

rownames(sampledata)
sampledata=subset(sampledata, rownames(sampledata) %in% colnames(otumat))

#otumat and sampledata should have the same sample specimen order

index <- match(rownames(sampledata), colnames(otumat))
otumat  <- otumat[,index]
head(otumat)
  all(rownames(sampledata) == colnames(otumat)) #it should return TRUE. If FALSE the two tables do not match


#### IMPORT DATA INTO PHYLOSEQ ####

# see Phyloseq webpage https://joey711.github.io/phyloseq/index.html

### Phyloseq Installation

#source('http://bioconductor.org/biocLite.R')
#BiocManager::install('phyloseq')

#Load library
library(phyloseq)
  
#Create a phyloseq object, inckuding the three tables

SSV = otu_table(otumat, taxa_are_rows = TRUE) #Abundance table of microbial taxa per sample
TAX = tax_table(as.matrix(taxmat)) #Taxonomy table of SSVs

physeq = phyloseq(SSV, TAX)
physeq

#Prepare sample metadata and add it to physeq

sampledata = sample_data(data.frame(sampledata, 
                                    row.names=sample_names(physeq),
                                    stringsAsFactors=FALSE)
) 
sampledata

cichlids = merge_phyloseq(physeq, sampledata) # We merge the metadata table with the rest of the data
cichlids


### SAVE THE OBJECT AS A RDA ####

save(cichlids, file = "CICHLIDS_SILVA/cichlids_SILVA.rda")
cichlids


