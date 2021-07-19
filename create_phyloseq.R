#Loading R packages
library("phyloseq"); packageVersion("phyloseq")
library("ggplot2"); packageVersion("ggplot2")

#Defining a default theme for ggplot graphics
theme_set(theme_bw())

#Simulating a mock OTU table
OTU = matrix(sample(1:100, 1000, replace = TRUE), nrow = 100, ncol = 20)
OTU

#Defining the attributes (column and row)
rownames(OTU) <- paste0("OTU", 1:nrow(OTU))
colnames(OTU) <- paste0("Sample", 1:ncol(OTU))
OTU

#Simulating a mock taxonomy table
TAX = matrix(sample(letters, 70, replace = TRUE), nrow = nrow(OTU), ncol = 7)
rownames(TAX) <- rownames(OTU)
colnames(TAX) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
TAX

#Checking the type of objects (OTU and TAX)
class(OTU)
class(TAX)

#Constructing Operational Taxonomic Unit (OTU) abundance
OTU = otu_table(OTU, taxa_are_rows = TRUE) #to convert matrix data format
TAX = tax_table(TAX) #to convert matrix data format
OTU
TAX

#Creating phyloseq object (load the phyloseq package first)
physeq = phyloseq (OTU,TAX)
physeq

#Plotting the composition (just to check whether it works or not)
plot_bar(physeq, fill = "Phylum")

#Simulating metadata
#1
SAM1 = sample_data(data.frame(
  Location = sample(LETTERS[1:4], size=nsamples(physeq), replace=TRUE),
  Depth = sample(50:1000, size=nsamples(physeq), replace=TRUE),
  row.names=sample_names(physeq),
  stringsAsFactors=FALSE
))
SAM1

#2
SAM2 = sample_data(data.frame(
  Niche = c("Bulk soil", "Bulk soil", "Bulk soil", "Bulk soil", "Bulk soil", "Rhizosphere","Rhizosphere","Rhizosphere","Rhizosphere","Rhizosphere", "Bulk soil", "Bulk soil", "Bulk soil", "Bulk soil", "Bulk soil", "Rhizosphere","Rhizosphere","Rhizosphere","Rhizosphere","Rhizosphere"),
  Treatment = c("Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Fertilized",  "Fertilized", "Fertilized", "Fertilized", "Fertilized", "Fertilized", "Fertilized", "Fertilized", "Fertilized", "Fertilized"),
  row.names=sample_names(physeq),
  stringsAsFactors = FALSE
))
SAM2

#Simulating phylogenetic tree (only applies to futher analysis that requires to check for phylogenetic relationship)
library("ape")
random_tree = rtree(ntaxa(physeq), rooted=TRUE, tip.label=taxa_names(physeq))
plot(random_tree)

#4 input files have been created - OTU, TAX, SAM, random_tree
#Merging the input OR
physeq1 = merge_phyloseq(physeq, SAM2, random_tree)
physeq1

#Rebuld from scratch
physeq2 = phyloseq (OTU, SAM2, TAX, random_tree)
physeq2

# mean, max and min of sample read counts
smin <- min(sample_sums(physeq))
smean <- mean(sample_sums(physeq))
smax <- max(sample_sums(physeq))
smin
smean
smax

#Identical check
identical(physeq1, physeq2)
