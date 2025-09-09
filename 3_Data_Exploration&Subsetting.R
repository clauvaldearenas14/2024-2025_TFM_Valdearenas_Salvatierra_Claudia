##################################################################
##### DATA EXPLORATION ####
##################################################################


#### LOAD LIBRARIES ####

library(vegan)
library(ape)
library(phyloseq)
library(ggplot2)
library(MicrobiomeStat)
library(dplyr)


#### LOAD THE DATA ####

load("CICHLIDS_SILVA/cichlids_SILVA.rda") 
cichlids # See the content of the "cichlids" phyloseq object
#This object contains three tables and a phylogenetic tree: the abundance table (otu_table), the sample metadata table (sample_data), the taxonomy table (tax_table) and the phylogenetic tree (phy_tree)


#### EXPLORE THE DATA ####
#"samples" are the specimens, "taxa" are the bacteria 

ntaxa(cichlids) # Number of bacterial taxa (IDs), the distinct bacteria, here called "taxa".
taxa_names(cichlids)[1:10] # Name of  IDs for a subset of the data (first 10 taxa)
taxa_sums(cichlids) #Number of reads per each IDs 
nsamples(cichlids) # Total number of samples
sample_names(cichlids)[1:5] # Sample names for the subset
sample_sums(cichlids)#Total number of taxa reads per sample
rank_names(cichlids) # Taxonomic level of IDs
sample_variables(cichlids) # Variables associated to the sample
get_variable(cichlids, "Sex") # Visualizing an individual variable
otu_table(cichlids)[1:5, 1:5] # Visualizing a subset of the abundance table (first five rows and columns)
tax_table(cichlids)[1:5,c(1,7)] # Visualizing the taxonomic level for a subset of samples (first five rows and all seven columns)

data.frame(sample_sums(cichlids)) # number of reads per sample
data.frame(taxa_sums(cichlids))# number of reads per taxa

sum(taxa_sums(cichlids))#total number of reads

#### SUMMARIZE DATA  ####

metadata <- sample_data(cichlids)
otutable <-otu_table(cichlids)

#Summarize data
#check available variables
sample_variables(cichlids) 

# Frequency of variables: Sex
metadata %>% group_by(Sex) %>% summarise(n = n())
# Frequency of variables: Sex and Diet
metadata %>% group_by(Diet, Sex) %>% summarise(n = n())

### SUBSET THE DATA  ####

#The following are key phyloseq functions for filtering/subsetting
#prune_samples #it eliminates samples(individuals) based on some selection of variables
#subset_sample # it retains only some samples based on some selection of variables

#prune_taxa #it eliminates bacterial taxa based on read count (abundance), or taxonomy
#subset_taxa #it retains only some bacterial taxa based on some selection of variables
#filter_taxa #it eliminates bacterial taxa based on read count (abundance) and or frequency of occurrence across samples


### SUBSET SAMPLES ####

#BY SAMPLE ABUNDANCE
# keep only samples with total sequence count greater than 2000
prune_samples(sample_sums(cichlids) >= 2000, cichlids)

#BY A VARIABLE IN THE METADATA

#To filter according to some variables (e.g., Diet)

sample_variables(cichlids)#First check the available variables

get_variable(cichlids, "Diet") # check the content of the variable "Islet"
unique(get_variable(cichlids, "Diet"))
unique(get_variable(cichlids, "Sex"))# check only unique entries for Sample_type


cichlids_alg = subset_samples(cichlids, Diet == "Algae") #retain only Algae
cichlids_algfem = subset_samples(cichlids_alg, Sex == "Female") #retain only females from En Curt

cichlids_algfem


### SUBSET TAXA (microbial taxa) ####

###BY TAXA ABUNDANCE

# Remove singletons (sequences with only 1 count)
prune_taxa(taxa_sums(cichlids) >1, cichlids)

# Keep only taxa that with total sum of counts = or >2
prune_taxa(taxa_sums(cichlids) >= 2, cichlids)

# retain only  the first 1000 most abundant bacterial taxa
prune_taxa(names(sort(taxa_sums(cichlids),TRUE)[1:1000]), cichlids) 

#filter to just taxa that appeared in more than one sample
filter_taxa(cichlids, function(x){sum(x > 0) > 1}, prune = TRUE)

#filter to just taxa that appeared in more than one sample (prevalence >= 2) and at least 20% of the samples. 

filter_taxa(cichlids, function(x) sum(x > 2) > (0.2*length(x)), TRUE)


###BY TAXA NAME

#Keep Only EUKARYOTAS (there are some)
rank_names(cichlids)
get_taxa_unique(cichlids, "Domain") #get unique Domains
euk<-subset_taxa(cichlids, Domain == "Eukaryota")
euk


#Keep Only BACTERIA (there are some EUKARYOTAS)
rank_names(cichlids)
get_taxa_unique(cichlids, "Domain") #get unique Domains
bacteria<-subset_taxa(cichlids, Domain == "Bacteria")
bacteria
cichlids <- prune_taxa(taxa_sums(bacteria) > 0, bacteria)


#After subsetting, remove bacteria taxa with no counts

cichlids <- prune_taxa(taxa_sums(cichlids) > 0, cichlids)


####EXPORT THE NEW FILTERED DATA ####

#As R object (rda)

cichlids<-bacteria
save(cichlids, file = "cichlids_SILVA_bacteria.rda")


#As txt file
#the abundance table
cichlids
otumat = as(otu_table(cichlids), "matrix")
# Coerce to data.frame
OTUdf = as.data.frame(otumat)
write.table(OTUdf, "otumat.txt", sep="\t", quote=F)

#the taxonomy table
TAX = as(tax_table(cichlids), "matrix")
# Coerce to data.frame
TAXdf = as.data.frame(TAX)
write.table(TAXdf, "taxmat.txt", sep="\t", quote=F)

