# creatingphyloseq
Phyloseq is mainly a package to analyze data consisting of operational taxonomic unit (OTU), metadata, taxonomy, and phylogenetic tree (optional).

Objectives of the repository: 
1. To create a mock OTU table, taxonomy, metadata, and phylogenetic tree
2. To merge the input files as a phyloseq object for further metagenomic analysis

The scripts are adapted from https://joey711.github.io/phyloseq/

## Loading R Packages
```bash
library("phyloseq"); packageVersion("phyloseq")
library("ggplot2"); packageVersion("ggplot2")
```

## Defining a default theme for ggplot graphics
```
theme_set(theme_bw())
```

## 1. Creating a mock OTU table (OTU)
```
OTU = matrix(sample(1:100, 1000, replace = TRUE), nrow = 100, ncol = 20)
OTU
```

## Defining the attributes (column and row)
```
rownames(OTU) <- paste0("OTU", 1:nrow(OTU))
colnames(OTU) <- paste0("Sample", 1:ncol(OTU))
OTU
```

## 2. Creating a mock taxonomy table (TAX)
```
TAX = matrix(sample(letters, 70, replace = TRUE), nrow = nrow(OTU), ncol = 7)
rownames(TAX) <- rownames(OTU)
colnames(TAX) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
TAX
```

## Checking the type of objects (OTU and TAX)
```
class(OTU)
class(TAX)
```

##Constructing Operational Taxonomic Unit (OTU) abundance
```
OTU = otu_table(OTU, taxa_are_rows = TRUE) #to convert matrix data format
TAX = tax_table(TAX) #to convert matrix data format
OTU
TAX
```

## Creating phyloseq object
```
physeq = phyloseq (OTU,TAX)
physeq
```

## Plotting the composition (just to check whether it works or not)
```
plot_bar(physeq, fill = "Phylum")
```

## 3. Creating mock metadata
Option 1:
```
SAM1 = sample_data(data.frame(
  Location = sample(LETTERS[1:4], size=nsamples(physeq), replace=TRUE),
  Depth = sample(50:1000, size=nsamples(physeq), replace=TRUE),
  row.names=sample_names(physeq),
  stringsAsFactors=FALSE
))
SAM1
```

Option 2:
```
SAM2 = sample_data(data.frame(
  Niche = c("Bulk soil", "Bulk soil", "Bulk soil", "Bulk soil", "Bulk soil", "Rhizosphere","Rhizosphere","Rhizosphere","Rhizosphere","Rhizosphere", "Bulk soil", "Bulk soil", "Bulk soil", "Bulk soil", "Bulk soil", "Rhizosphere","Rhizosphere","Rhizosphere","Rhizosphere","Rhizosphere"),
  Treatment = c("Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Control", "Fertilized",  "Fertilized", "Fertilized", "Fertilized", "Fertilized", "Fertilized", "Fertilized", "Fertilized", "Fertilized", "Fertilized"),
  row.names=sample_names(physeq),
  stringsAsFactors = FALSE
))
SAM2
```

## 4. Creating mock phylogenetic tree (random_tree)
```
library("ape")
random_tree = rtree(ntaxa(physeq), rooted=TRUE, tip.label=taxa_names(physeq))
plot(random_tree)
```

## Creating phyloseq object
Input objects: OTU, TAX, SAM, random_tree

Option 1: Merging the input OR
```
physeq1 = merge_phyloseq(physeq, SAM2, random_tree)
physeq1
```

Option 2: Rebuld from scratch
```
physeq2 = phyloseq (OTU, SAM2, TAX, random_tree)
physeq2
```

## mean, max and min of sample read counts
```
smin <- min(sample_sums(physeq))
smean <- mean(sample_sums(physeq))
smax <- max(sample_sums(physeq))
smin
smean
smax
```

## Identical check
```
identical(physeq1, physeq2)
```
