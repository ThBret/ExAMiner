#!/usr/bin/env Rscript

##################################################
### ARGUMENTS
##################################################

suppressPackageStartupMessages(library(argparser))

banner <- r"{


 ██████╗ ██╗  ██╗██╗   ██╗███████╗███████╗██╗   ██╗██████╗ 
 ██╔══██╗██║  ██║╚██╗ ██╔╝██╔════╝██╔════╝██║   ██║██╔══██╗
 ██████╔╝███████║ ╚████╔╝ ███████╗███████╗██║   ██║██████╔╝
 ██╔═══╝ ██╔══██║  ╚██╔╝  ╚════██║╚════██║██║   ██║██╔══██╗
 ██║     ██║  ██║   ██║   ███████║███████║╚██████╔╝██║  ██║
 ╚═╝     ╚═╝  ╚═╝   ╚═╝   ╚══════╝╚══════╝ ╚═════╝ ╚═╝  ╚═╝

======== PHYlogram Statistics and SUbsampling in R ========

}"

p <- arg_parser(cat(banner), hide.opts = TRUE)

p <- add_argument(p, "--input", help="Full path to directory with input files. Example: /directory/with/input/starting/from/root/")
p <- add_argument(p, "--output", help="Name of output directory", default="physsur_output")
p <- add_argument(p, "--remove", help="Add flag to remove outliers", flag=TRUE)
p <- add_argument(p, "--boot", help="Threshold for average bootstrap support", default=0)
p <- add_argument(p, "--prop", help="Proportion of kept loci after filtering (loci are selected based on clocklikeness score)", default=1)

argv <- parse_args(p)

print(p)

##################################################
### CHECK USER INPUT
##################################################

stopifnot(

  "Please provide input directory" =
    (argv$input != "NA"),
  "Please provide a value between 0 and 100 for average bootstrap support threshold" =
    (argv$boot >= 0 & argv$boot <= 100),
  "Please provide a proportion of final loci set between 0 and 1" =
    (argv$prop >= 0 & argv$prop <= 1)

  )

##################################################
### LIBRARIES
##################################################

print("-----------------------------------------")
print("Loading packages")
print("-----------------------------------------")

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(phangorn))
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(seqinr))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(QSutils))
suppressPackageStartupMessages(library(TreeDist))
suppressPackageStartupMessages(library(rstatix))

##################################################
### SETUP
##################################################

print("-----------------------------------------")
print("Setting-up variables and functions")
print("-----------------------------------------")

### Set working directory

setwd(argv$input)

### Set output directories
dir.create(argv$output)
dir.create(file.path(argv$output, "tree_plots"))

tre_out_dir <- file.path(argv$output, "tree_plots")

### Set variables

input_dir <- file.path(argv$input)

trees_files <- dir(path=input_dir, pattern="*treefile$")
runs_files <- dir(path=input_dir, pattern="*runtrees$")
aln_files <- dir(path=input_dir, pattern="*all.fasta$")
ingroup_aln_files <- dir(path=input_dir, pattern="ingroup.fasta")

tree_regex <- "(\\S+).treefile"
runs_regex <- "(\\S+).runtrees"
aln_regex <- "(\\S+).all.fasta"
ingroup_aln_regex <- "(\\S+).ingroup.fasta"

##################################################
### DEFINE FUNCTIONS
##################################################

### Proportion of ingroup variable sites

ingroup_var <- function(file) {

  # get locus name from filename
  locus_no <- sub(ingroup_aln_regex, "\\1", perl=TRUE, x=file)
  # read in tree
  ingroup_aln <- readDNAStringSet(file.path(input_dir, file))
  # store proportion of variable sites in a vector
  prop_var <- (SegSites(ingroup_aln)/width(ingroup_aln[1]))
  return(c(locus_no,prop_var))

}

###  Average Clustering Information Distance between different runs of the same locus

ClustInfDist_mean <- function(file) {

  # get locus name from filename
  locus_no <- sub(runs_regex, "\\1", perl=TRUE, x=file)
  # read in tree
  runs <- read.tree(file.path(input_dir, file))
  # get the Clustering Information Distance
  dist_list <- unlist(as.list(TreeDistance(runs)))
  # get the percent of tree comparisons that are identical (i.e., clustering information distance = 0) 
  avg_ClustInfDist <- mean(dist_list)
  return(c(locus_no,avg_ClustInfDist))

}

### Average bootstrap support

Avg_support <- function(file) {

  # get UCE number from filename
  locus_no <- sub(tree_regex, "\\1", perl=TRUE, x=file)
  # read in tree
  tree <- read.tree(file.path(input_dir, file))
  # store support values in a vector
  support <- c(as.numeric(tree$node.label))
  # calculate average support
  avg_supp <- mean(support, na.rm=T)
  return(c(locus_no,avg_supp))

}

### Average branch lenghts

Br_length.trees <- function(file) {

  # reads the phylogenetic tree
  tree <- read.tree(file.path(input_dir, file))
  # gets number of tips
  no_tips <- length(tree$tip.label)
  # calculate avg branch length
  avg_br_length <- mean(tree$edge.length)
  # get UCE number from filename
  locus_no <- sub(tree_regex, "\\1", perl=TRUE, x=file)
  return(c(locus_no,avg_br_length))

}

### Clock-likeness

Clocklikeness <- function(file) {

  # get UCE number from filename
  locus_no <- sub(tree_regex, "\\1", perl=TRUE, x=file)
  # read in tree
  tree <- read.tree(file.path(input_dir, file))
  # record coefficient for all possible outgroups
  CV <- c()
  taxa <- tree$tip.label
  for (taxon in taxa) {
    # root tree
    rooted_tr <- root(phy=tree,outgroup=taxon,resolve.root=T)
    # get matrix diagonal of phylogenetic variance-covariance matrix
    # these are your distances from root
    root_dist <- diag(vcv.phylo(rooted_tr))
    std_dev_root_dist <- sd(root_dist)
    mean_root_dist <- mean(root_dist)
    CV <- c(CV, (std_dev_root_dist/mean_root_dist)*100)
  }
  # get lowest CV
  minCV <- min(CV)
  print(minCV)
  return(c(locus_no,minCV))

}

### Plot gene trees 

Plot_trees <- function(file) {

  # reads the phylogenetic tree
  tree <- read.tree(file.path(input_dir, file))
  # extracts plot name (locus) from file name 
  plot_name <- sub(tree_regex, "\\1", perl=TRUE, x=file)
  # open png file
  pdf(file=paste(tre_out_dir, plot_name, "-tree.pdf", sep=""), width=20, height=15)
  plot.phylo(ladderize(tree), show.node.label=TRUE)
  # give title as locus number and subtitle as tree length
  title(main=plot_name)
  # close png file
  dev.off()

}

##################################################
### CALCULATE STATISTICS
##################################################

print("-----------------------------------------")
print("Calculating statistics")
print("-----------------------------------------")

## Get the proportion of ingroup variable sites per locus

var_sites <- lapply(ingroup_aln_files, ingroup_var)
var_sites <- data.frame(matrix(unlist(var_sites), nrow=(length(var_sites)), byrow=T))
colnames(var_sites) <- c("Locus", "Ingroup_variation")

## Get average clustering information distance

average_ClustInfDist <- lapply(runs_files, ClustInfDist_mean)
average_ClustInfDist <- data.frame(matrix(unlist(average_ClustInfDist), nrow=(length(average_ClustInfDist)), byrow=T))
colnames(average_ClustInfDist) <- c("Locus", "Average_CID")

## Get average bootstrap support per tree

average_bootstrap <- lapply(trees_files, Avg_support)
average_bootstrap <- data.frame(matrix(unlist(average_bootstrap), nrow=(length(average_bootstrap)), byrow=T))
colnames(average_bootstrap) <- c("Locus", "Average_bootstrap")

## Get average branchlength per tree

br_lengths <- lapply(trees_files, Br_length.trees)
br_lengths <- data.frame(matrix(unlist(br_lengths), nrow=(length(br_lengths)), byrow=T))
colnames(br_lengths) <- c("Locus", "Average_branch_length")

## Get clock-likeness score per tree

cv_clocklikeness <- lapply(trees_files, Clocklikeness)
cv_clocklikeness <- data.frame(matrix(unlist(cv_clocklikeness), nrow=(length(cv_clocklikeness)), byrow=T))
colnames(cv_clocklikeness) <- c("Locus", "Clocklikeness")

## Plot gene trees

#loop over all files
lapply(trees_files, Plot_trees)

##################################################
### COLLECT AND WRITE PRE-FILTERING FILES
##################################################

full_table <- inner_join(
  average_ClustInfDist, inner_join(
    average_bootstrap, inner_join(
      cv_clocklikeness, inner_join(
        br_lengths, var_sites, by="Locus"), by="Locus"), by= "Locus"), by="Locus") %>%
  mutate_at(

      c("Average_CID",
      "Average_bootstrap",
      "Clocklikeness",
      "Average_branch_length",
      "Ingroup_variation"),

    as.numeric

    )

write_csv(full_table, file=paste(argv$output, "gene_tree_stats_prefilter.csv", sep="/"))
write_tsv(full_table, file=paste(argv$output, "gene_tree_stats_prefilter.tsv", sep="/"))

##################################################
### FILTERS
##################################################

if(argv$remove == TRUE) {

  print("-----------------------------------------")
  print("Removing dataset outliers")
  print("-----------------------------------------")

  loci_outliers <- full_table %>%
    select(Locus,Average_CID,Average_bootstrap,Ingroup_variation,Clocklikeness) %>% 
    pivot_longer(-Locus) %>% 
    group_by(name) %>% 
    identify_outliers(value) %>% 
    mutate(type = case_when(

      name == "Average_CID" & value < median(full_table$Average_CID, na.rm=TRUE) ~ "low",
      name == "Average_CID" & value > median(full_table$Average_CID, na.rm=TRUE) ~ "high",
      name == "Average_bootstrap" & value < median(full_table$Average_bootstrap, na.rm=TRUE) ~ "low",
      name == "Average_bootstrap" & value > median(full_table$Average_bootstrap, na.rm=TRUE) ~ "high",
      name == "Ingroup_variation" & value < median(full_table$Ingroup_variation, na.rm=TRUE) ~ "low",
      name == "Ingroup_variation" & value > median(full_table$Ingroup_variation, na.rm=TRUE) ~ "high",
      name == "Clocklikeness" & value < median(full_table$Clocklikeness, na.rm=TRUE) ~ "low",
      name == "Clocklikeness" & value > median(full_table$Clocklikeness, na.rm=TRUE) ~ "high"

      )

    ) %>%
    filter(

      name == "Clocklikeness" & is.extreme == "TRUE" & type == "high" | 
      name == "Ingroup_variation" & is.outlier == "TRUE" | 
      name == "Average_CID" & is.outlier == "TRUE" & type == "high" |
      name == "Average_bootstrap" & is.outlier == "TRUE" & type == "low")

  loci_outliers_list <- loci_outliers %>%
    select(Locus) %>%
    unique()

  filt_table <- anti_join(full_table,loci_outliers_list)

  print("-----------------------------------------")
  print("Applying extra filters")
  print("-----------------------------------------")

  filtered_data <- filt_table %>%
    filter(Average_bootstrap >= argv$boot) %>% 
    arrange(Clocklikeness) %>% 
    slice_head(prop=argv$prop)

  write_csv(loci_outliers, file=paste(argv$output, "removed_outliers.csv", sep="/"))
  write_tsv(loci_outliers, file=paste(argv$output, "removed_outliers.tsv", sep="/"))
  write_csv(filtered_data, file=paste(argv$output, "gene_tree_stats_postfilter.csv", sep="/"))
  write_tsv(filtered_data, file=paste(argv$output, "gene_tree_stats_postfilter.tsv", sep="/"))
  write_lines(filtered_data$Locus, file=paste(argv$output, "whitelist.txt", sep="/"))

} else {

  print("-----------------------------------------")
  print("Applying filters")
  print("-----------------------------------------")

  filtered_data <- full_table %>%
    filter(Average_bootstrap >= argv$boot) %>%
    arrange(Clocklikeness) %>%
    slice_head(prop=argv$prop)

  write_csv(filtered_data, file=paste(argv$output, "gene_tree_stats_postfilter.csv", sep="/"))
  write_tsv(filtered_data, file=paste(argv$output, "gene_tree_stats_postfilter.tsv", sep="/"))
  write_lines(filtered_data$Locus, file=paste(argv$output, "whitelist.txt", sep="/"))

}

##################################################
### PLOTS
##################################################

print("-----------------------------------------")
print("Plotting results")
print("-----------------------------------------")

a <- full_table %>%
  pivot_longer(-Locus) %>% 
  group_by(name) %>%
  ggboxplot(

    x="name",
    y="value",
    color = "#FFA319FF",
    linetype = "solid",
    size = 0.5,
    width = 0.5,
    outlier.shape = 1,
    xlab = FALSE,
    ylab = FALSE,
    title = "PRE-FILTERING DATASET STATISTICS",
    ggtheme = theme_bw()

    ) %>%
  facet(

    facet.by = "name",
    scales = "free",
    nrow = 1
    
    )

b <- filtered_data %>%
  pivot_longer(-Locus) %>% 
  group_by(name) %>%
  ggboxplot(

    x="name",
    y="value",
    color = "#155F83FF",
    linetype = "solid",
    size = 0.5,
    width = 0.5,
    outlier.shape = 1,
    xlab = FALSE,
    ylab = FALSE,
    title = "POST-FILTERING DATASET STATISTICS",
    ggtheme = theme_bw()

    ) %>%
  facet(

    facet.by = "name",
    scales = "free",
    nrow = 1

    )

final_plot <- ggarrange(a, b, labels = c("A", "B"), ncol = 1, nrow = 2)

ggsave(paste(argv$output, "physsur_stats.pdf", sep="/"), width = 30, height = 30, units = "cm")
