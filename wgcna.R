# no longer need to install here because it was installed in our base image
#install.packages("BiocManager")
#install.packages(c('BiocManager'), repos='https://cloud.r-project.org/');BiocManager::install('WGCNA')
library(WGCNA)
 
install.packages("tidyverse")
library(tidyverse)

# read in the normalized expression
# expecting the output from DESeq2 where data are normalized
#data <- readr::read_delim("/sbgenomics/project-files/test_data_GenePhenotypeFile.csv",  
#                          delim = ",")
#take input from gene-median-splitter matrices - 2 one for HighMyc and one for LowMyc

matrix <- commandArgs(trailingOnly=TRUE)
data <- read.csv(matrix, sep="\t")

data[1:9,1:9]

de_input = as.matrix(data[,-1])
row.names(de_input) = data$GeneId
de_input[1:9,1:9]

meta_df <- data.frame( Sample = names(data[-1])) %>%
  mutate(
    Type = gsub("-.*","", Sample) %>% gsub("[.].*","", .)
  )
meta_df

#the deseq normlization step is executed outside of here
# input_mat <- t(expr_normalized) 
input_mat<- t(de_input)
input_mat[1:9,1:9]

names(data)[1] = "GeneId"
names(data)           # Look at the column names

input_mat[1:9,1:9]

allowWGCNAThreads()          # allow multi-threading (optional)
#> Allowing multi-threading with up to 4 threads.

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to = 20, by = 2))

# Call the network topology analysis function
sft = pickSoftThreshold(
  input_mat,             # <= Input data
  #blockSize = 30,
  powerVector = powers,
  verbose = 5
  )

col_sel = names(data)[-1]     # Get all but first column name
col_sel

# Optional step ---  order the groups in the plot.
# mdata$group = factor(mdata$group,
#                     levels = c("B", "B_L1", ....))  #<= fill the rest of this in


mdata <- data %>%
  tidyr::pivot_longer(
    .,                        # The dot is the the input data, magrittr tutorial
    col = all_of(col_sel)
    ) %>%  
  mutate(
    group = gsub("-.*","", name) %>% gsub("[.].*","", .)   # Get the shorter treatment names
  )


# ==== Plot groups (Sample Groups vs RNA Seq Counts) to identify outliers
(
 p <- mdata %>%
    ggplot(., aes(x = name, y = value)) +             # x = treatment, y = RNA Seq count
    geom_violin() +                                   # violin plot, show distribution
    geom_point(alpha = 0.2) +                         # scatter plot
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90)          # Rotate treatment text
    ) +
    labs(x = "Treatment Groups", y = "RNA Seq Counts") +
    facet_grid(cols = vars(group), drop = TRUE, scales = "free_x")      # Facet by hour
)

par(mfrow = c(1,2));
cex1 = 0.9;

plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     main = paste("Scale independence")
)
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red"
)
abline(h = 0.90, col = "red")
plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity")
)
text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers,
     cex = cex1, col = "red")

    picked_power = 9
temp_cor <- cor       
cor <- WGCNA::cor         
# Force it to use WGCNA cor function (fix a namespace conflict issue)
# <= input here

typeof(input_mat)


input_mat

netwk <- blockwiseModules(input_mat,               
    # == Adjacency Function ==
    power = picked_power,
    # <= power here
    networkType = "signed",
    # == Tree and Block Options ==
    deepSplit = 2,
    pamRespectsDendro = F,
    # detectCutHeight = 0.75,
    minModuleSize = 30,
    maxBlockSize = 4000,

    # == Module Adjustments ==
    reassignThreshold = 0,
    mergeCutHeight = 0.25,

    # == TOM == Archive the run results in TOM file (saves time)
    saveTOMs = T,
    saveTOMFileBase = "ER",

    # == Output Options
    numericLabels = T,
    verbose = 3)
cor <- temp_cor     # Return cor function to original namespace

# Convert labels to colors for plotting
mergedColors = labels2colors(netwk$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(
  netwk$dendrograms[[1]],
  mergedColors[netwk$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05 )

# netwk$colors[netwk$blockGenes[[1]]]
# table(netwk$colors)
module_df <- data.frame(
  gene_id = names(netwk$colors),
  colors = labels2colors(netwk$colors)
)

module_df[1:5,]
#>            gene_id    colors
#> 1 AC149818.2_FG001      blue
#> 2 AC149829.2_FG003      blue
#> 3 AC182617.3_FG001      blue
#> 4 AC186512.3_FG001 turquoise
#> 5 AC186512.3_FG007 turquoise

write_delim(module_df,
    file = "gene_modules.txt",
    delim = "\t")

# Get Module Eigengenes per cluster
MEs0 <- moduleEigengenes(input_mat, mergedColors)$eigengenes

# Reorder modules so similar modules are next to each other
MEs0 <- orderMEs(MEs0)
module_order = names(MEs0) %>% gsub("ME","", .)

# Add treatment names
MEs0$treatment = row.names(MEs0)

# tidy & plot data
mME = MEs0 %>%
  pivot_longer(-treatment) %>%
  mutate(
    name = gsub("ME", "", name),
    name = factor(name, levels = module_order)
  )

mME %>% ggplot(., aes(x=treatment, y=name, fill=value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1,1)) +
  theme(axis.text.x = element_text(angle=90)) +
  labs(title = "Module-trait Relationships", y = "Modules", fill="corr")




