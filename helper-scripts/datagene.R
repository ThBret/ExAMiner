#!/usr/bin/env Rscript

##################################################
### ARGUMENTS
##################################################
suppressPackageStartupMessages(library(argparser))

banner <- r"{

██████╗  █████╗ ████████╗ █████╗  ██████╗ ███████╗███╗   ██╗███████╗██████╗
██╔══██╗██╔══██╗╚══██╔══╝██╔══██╗██╔════╝ ██╔════╝████╗  ██║██╔════╝██╔══██╗
██║  ██║███████║   ██║   ███████║██║  ███╗█████╗  ██╔██╗ ██║█████╗  ██████╔╝
██║  ██║██╔══██║   ██║   ██╔══██║██║   ██║██╔══╝  ██║╚██╗██║██╔══╝  ██╔══██╗
██████╔╝██║  ██║   ██║   ██║  ██║╚██████╔╝███████╗██║ ╚████║███████╗██║  ██║
╚═════╝ ╚═╝  ╚═╝   ╚═╝   ╚═╝  ╚═╝ ╚═════╝ ╚══════╝╚═╝  ╚═══╝╚══════╝╚═╝  ╚═╝

============================ DATAset GENERator ============================

}"

p <- arg_parser(cat(banner), hide.opts = TRUE)

p <- add_argument(p, "--input", help="Full path to directory with input files. Example: /path/from/root/to/directory/with/input/")
p <- add_argument(p, "--output", help="Name of output directory", default="filter_output")
p <- add_argument(p, "--amas", help="Full path to AMAS table")
p <- add_argument(p, "--ml", help="Add flag to consider ML trees for filtering", flag=TRUE)
p <- add_argument(p, "--bi", help="Add flag to consider BI trees for filtering", flag=TRUE)
p <- add_argument(p, "--remove", help="Add flag to remove outliers", flag=TRUE)
p <- add_argument(p, "--taxa", help="Threshold for number of taxa", default=0)
p <- add_argument(p, "--min_length", help="Low threshold for locus length", default=0)
p <- add_argument(p, "--max_length", help="Upper threshold for locus length", default=999999)
p <- add_argument(p, "--missing", help="Threshold for missing data percent", default=100)
p <- add_argument(p, "--min_pars", help="Lower percentile threshold for ingroup parsimony sites", default=0)
p <- add_argument(p, "--max_pars", help="Upper percentile threshold for ingroup parsimony sites", default=1)
p <- add_argument(p, "--cid", help="Threshold for normalised Clustering Information Distance between trees", default=1)
p <- add_argument(p, "--supp_ml", help="Threshold for low average bootstrap support", default=0)
p <- add_argument(p, "--supp_bi", help="Threshold for low average posterior probability", default=0)
p <- add_argument(p, "--clock_ml", help="Number of loci to keep after filtering (loci are selected based on clocklikeness score)", default=999999)
p <- add_argument(p, "--clock_bi", help="Number of loci to keep after filtering (loci are selected based on clocklikeness score)", default=999999)

argv <- parse_args(p)

print(p)

##################################################
### CHECK USER INPUT
##################################################

stopifnot(

  "Please provide input directory" =
    (argv$input != "NA"),
  "Please provide a proportion of final loci set between 0 and 1" =
    (argv$prop >= 0 & argv$prop <= 1)

  )

##################################################
### LIBRARIES
##################################################
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(phangorn))
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(TreeDist))
suppressPackageStartupMessages(library(rstatix))
suppressPackageStartupMessages(library(ggpubr))

##################################################
### SETUP
##################################################

### Set working directory

setwd(argv$input)

### Set output directories
dir.create(argv$output)

input_dir <- file.path(argv$input)

ml_tree_files <- dir(path=input_dir, pattern="*ML.treefile$")
bi_tree_files <- dir(path=input_dir, pattern="*BI.treefile$")
run_files <- dir(path=input_dir, pattern="*trees$")
#aln_files <- dir(path=input_dir, pattern="*all.fasta$")
#ingroup_aln_files <- dir(path=input_dir, pattern="ingroup.fasta")

ml_tree_regex <- "(\\S+).ML.treefile"
bi_tree_regex <- "(\\S+).BI.treefile"
runs_regex <- "(\\S+).trees"
#aln_regex <- "(\\S+).all.fasta"
#ingroup_aln_regex <- "(\\S+).ingroup.fasta"

##################################################
### DEFINE FUNCTIONS
##################################################

ClustInfDist_ml_bi <- function(file) {

  # get locus name from filename
  locus_no <- sub(runs_regex, "\\1", perl=TRUE, x=file)
  # read in tree
  runs <- read.tree(file.path(input_dir, file))
  # get the Clustering Information Distance
  dist_val <- unlist(as.list(TreeDistance(runs)))
  return(c(locus_no,dist_val))

}

### Average support

Avg_boot_support <- function(file) {

  # get UCE number from filename
  locus_no <- sub(ml_tree_regex, "\\1", perl=TRUE, x=file)
  # read in tree
  ml_tree <- read.tree(file.path(input_dir, file))
  # store support values in a vector
  boot_support <- c(as.numeric(ml_tree$node.label))
  # calculate average support
  avg_boot_supp <- mean(boot_support, na.rm=T)
  return(c(locus_no,avg_boot_supp))

}

Avg_prob_support <- function(file) {

  # get UCE number from filename
  locus_no <- sub(bi_tree_regex, "\\1", perl=TRUE, x=file)
  # read in tree
  bi_tree <- read.tree(file.path(input_dir, file))
  # store support values in a vector
  prob_support <- c(as.numeric(bi_tree$node.label))
  # calculate average support
  avg_prob_supp <- mean(prob_support, na.rm=T)
  return(c(locus_no,avg_prob_supp))

}

### Average branch lenghts

Br_length_ml_trees <- function(file) {

  # reads the phylogenetic tree
  ml_tree <- read.tree(file.path(input_dir, file))
  # gets number of tips
  no_tips <- length(ml_tree$tip.label)
  # calculate avg branch length
  avg_br_length <- mean(ml_tree$edge.length)
  # get UCE number from filename
  locus_no <- sub(ml_tree_regex, "\\1", perl=TRUE, x=file)
  return(c(locus_no,avg_br_length))

}

Br_length_bi_trees <- function(file) {

  # reads the phylogenetic tree
  bi_tree <- read.tree(file.path(input_dir, file))
  # gets number of tips
  no_tips <- length(bi_tree$tip.label)
  # calculate avg branch length
  avg_br_length <- mean(bi_tree$edge.length)
  # get UCE number from filename
  locus_no <- sub(bi_tree_regex, "\\1", perl=TRUE, x=file)
  return(c(locus_no,avg_br_length))

}

### Clock-likeness

Clocklikeness_ml <- function(file) {

  # get UCE number from filename
  locus_no <- sub(ml_tree_regex, "\\1", perl=TRUE, x=file)
  # read in tree
  ml_tree <- read.tree(file.path(input_dir, file))
  # record coefficient for all possible outgroups
  CV <- c()
  taxa <- ml_tree$tip.label
  for (taxon in taxa) {
    # root tree
    rooted_tr <- root(phy=ml_tree,outgroup=taxon,resolve.root=T)
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

Clocklikeness_bi <- function(file) {

  # get UCE number from filename
  locus_no <- sub(bi_tree_regex, "\\1", perl=TRUE, x=file)
  # read in tree
  bi_tree <- read.tree(file.path(input_dir, file))
  # record coefficient for all possible outgroups
  CV <- c()
  taxa <- bi_tree$tip.label
  for (taxon in taxa) {
    # root tree
    rooted_tr <- root(phy=bi_tree,outgroup=taxon,resolve.root=T)
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

##################################################
### MAIN BODY
##################################################

## Get the relevant columns from the AMAS table

table_amas <- read_tsv(argv$amas) %>% 
  mutate(Locus = gsub(".ingroup.fasta","", Alignment_name)) %>%
  select(Locus, No_of_taxa, Alignment_length, Missing_percent, Proportion_parsimony_informative, GC_content)

## Get average clustering information distance between ML and BI trees

tree_distance <- lapply(run_files, ClustInfDist_ml_bi)
tree_distance <- data.frame(matrix(unlist(tree_distance), nrow=(length(tree_distance)), byrow=T))
colnames(tree_distance) <- c("Locus", "CID")

## Get average branch support per ML tree

average_bootstrap <- lapply(ml_tree_files, Avg_boot_support)
average_bootstrap <- data.frame(matrix(unlist(average_bootstrap), nrow=(length(average_bootstrap)), byrow=T))
colnames(average_bootstrap) <- c("Locus", "Average_bootstrap_ML")

## Get average branch support per BI tree

average_prob <- lapply(bi_tree_files, Avg_prob_support)
average_prob <- data.frame(matrix(unlist(average_prob), nrow=(length(average_prob)), byrow=T))
colnames(average_prob) <- c("Locus", "Average_prob_BI")

## Get average branchlength per ML tree

br_lengths_ml <- lapply(ml_tree_files, Br_length_ml_trees)
br_lengths_ml <- data.frame(matrix(unlist(br_lengths_ml), nrow=(length(br_lengths_ml)), byrow=T))
colnames(br_lengths_ml) <- c("Locus", "Average_branch_length_ML")

## Get average branchlength per BI tree

br_lengths_bi <- lapply(bi_tree_files, Br_length_bi_trees)
br_lengths_bi <- data.frame(matrix(unlist(br_lengths_bi), nrow=(length(br_lengths_bi)), byrow=T))
colnames(br_lengths_bi) <- c("Locus", "Average_branch_length_BI")

## Get clock-likeness score per ML tree

cv_clocklikeness_ml <- lapply(ml_tree_files, Clocklikeness_ml)
cv_clocklikeness_ml <- data.frame(matrix(unlist(cv_clocklikeness_ml), nrow=(length(cv_clocklikeness_ml)), byrow=T))
colnames(cv_clocklikeness_ml) <- c("Locus", "Clocklikeness_ML")

## Get clock-likeness score per BI tree

cv_clocklikeness_bi <- lapply(bi_tree_files, Clocklikeness_bi)
cv_clocklikeness_bi <- data.frame(matrix(unlist(cv_clocklikeness_bi), nrow=(length(cv_clocklikeness_bi)), byrow=T))
colnames(cv_clocklikeness_bi) <- c("Locus", "Clocklikeness_BI")


full_table <- inner_join(
  table_amas, inner_join(
    tree_distance, inner_join(
      average_bootstrap, inner_join(
        average_prob, inner_join(
          br_lengths_ml, inner_join(
            br_lengths_bi, inner_join(
              cv_clocklikeness_ml, cv_clocklikeness_bi, by="Locus"), by="Locus"), by= "Locus"), by="Locus"), by="Locus"), by= "Locus"), by="Locus") %>%
  mutate(across(-Locus, as.numeric))

write_csv(full_table, file=file.path(argv$output, "crossval_stats_prefilter.csv"))
write_tsv(full_table, file=file.path(argv$output, "crossval_stats_prefilter.tsv"))

##################################################
### FILTERS
##################################################

if(argv$remove == TRUE) {

  print("-----------------------------------------")
  print("Removing dataset outliers")
  print("-----------------------------------------")

  loci_outliers <- full_table %>%
    pivot_longer(-Locus) %>%
    group_by(name) %>%
    identify_outliers(value) %>%
    mutate(type = case_when(

      name == "No_of_taxa" & value < median(full_table$No_of_taxa, na.rm=TRUE) ~ "low",
      name == "No_of_taxa" & value > median(full_table$No_of_taxa, na.rm=TRUE) ~ "high",
      name == "Alignment_length" & value < median(full_table$Alignment_length, na.rm=TRUE) ~ "low",
      name == "Alignment_length" & value > median(full_table$Alignment_length, na.rm=TRUE) ~ "high",
      name == "Missing_percent" & value < median(full_table$Missing_percent, na.rm=TRUE) ~ "low",
      name == "Missing_percent" & value > median(full_table$Missing_percent, na.rm=TRUE) ~ "high",
      name == "Proportion_parsimony_informative" & value < median(full_table$Proportion_parsimony_informative, na.rm=TRUE) ~ "low",
      name == "Proportion_parsimony_informative" & value > median(full_table$Proportion_parsimony_informative, na.rm=TRUE) ~ "high",
      name == "GC_content" & value < median(full_table$GC_content, na.rm=TRUE) ~ "low",
      name == "GC_content" & value > median(full_table$GC_content, na.rm=TRUE) ~ "high",
      name == "CID" & value < median(full_table$CID, na.rm=TRUE) ~ "low",
      name == "CID" & value > median(full_table$CID, na.rm=TRUE) ~ "high",
      name == "Average_bootstrap_ML" & value < median(full_table$Average_bootstrap_ML, na.rm=TRUE) ~ "low",
      name == "Average_bootstrap_ML" & value > median(full_table$Average_bootstrap_ML, na.rm=TRUE) ~ "high",
      name == "Average_prob_BI" & value < median(full_table$Average_prob_BI, na.rm=TRUE) ~ "low",
      name == "Average_prob_BI" & value > median(full_table$Average_prob_BI, na.rm=TRUE) ~ "high",
      name == "Average_branch_length_ML" & value < median(full_table$Average_branch_length_ML, na.rm=TRUE) ~ "low",
      name == "Average_branch_length_ML" & value > median(full_table$Average_branch_length_ML, na.rm=TRUE) ~ "high",
      name == "Average_branch_length_BI" & value < median(full_table$Average_branch_length_BI, na.rm=TRUE) ~ "low",
      name == "Average_branch_length_BI" & value > median(full_table$Average_branch_length_BI, na.rm=TRUE) ~ "high",
      name == "Clocklikeness_ML" & value < median(full_table$Clocklikeness_ML, na.rm=TRUE) ~ "low",
      name == "Clocklikeness_ML" & value > median(full_table$Clocklikeness_ML, na.rm=TRUE) ~ "high",
      name == "Clocklikeness_BI" & value < median(full_table$Clocklikeness_BI, na.rm=TRUE) ~ "low",
      name == "Clocklikeness_BI" & value > median(full_table$Clocklikeness_BI, na.rm=TRUE) ~ "high"

      )

    ) %>%
    filter(

      name == "Missing_percent" & is.extreme == "TRUE" & type == "high" |
      name == "CID" & is.outlier == "TRUE"

      )

  loci_outliers_list <- loci_outliers %>%
    select(Locus) %>%
    unique()

  filt_table <- anti_join(full_table,loci_outliers_list)

  print("-----------------------------------------")
  print("Applying extra filters")
  print("-----------------------------------------")

  if(argv$ml == TRUE & argv$bi == FALSE) {

    filtered_data <- filt_table %>%
      filter(

        No_of_taxa >= argv$taxa &
        between(Alignment_length, 
          argv$min_length, 
          argv$max_length) &
        Missing_percent <= argv$missing &
        between(Proportion_parsimony_informative, 
          quantile(Proportion_parsimony_informative, argv$min_pars),
          quantile(Proportion_parsimony_informative, argv$max_pars)) &
        CID <= argv$cid &
        Average_bootstrap_ML >= argv$supp_ml

        ) %>%
      
      arrange(Clocklikeness_ML) %>%
      slice_head(n = argv$clock_ml)
     
  } else if(argv$ml == FALSE & argv$bi == TRUE) {

    filtered_data <- filt_table %>%
      filter(

        No_of_taxa >= argv$taxa &
        between(Alignment_length, 
          argv$min_length, 
          argv$max_length) &
        Missing_percent <= argv$missing &
        between(Proportion_parsimony_informative, 
          quantile(Proportion_parsimony_informative, argv$min_pars),
          quantile(Proportion_parsimony_informative, argv$max_pars)) &
        CID <= argv$cid &
        Average_prob_BI >= argv$supp_bi

        ) %>%
      
      arrange(Clocklikeness_BI) %>%
      slice_head(n = argv$clock_bi)

  } else {

    print("Please choose the --ml or --bi flag to apply filtering")

  }

} else {

  print("-----------------------------------------")
  print("Applying filters")
  print("-----------------------------------------")

  if(argv$ml == TRUE & argv$bi == FALSE) {

    filtered_data <- full_table %>%
      filter(

        No_of_taxa >= argv$taxa &
        between(Alignment_length, 
          argv$min_length, 
          argv$max_length) &
        Missing_percent <= argv$missing &
        between(Proportion_parsimony_informative, 
          quantile(Proportion_parsimony_informative, argv$min_pars),
          quantile(Proportion_parsimony_informative, argv$max_pars)) &
        CID <= argv$cid &
        Average_bootstrap_ML >= argv$supp_ml

        ) %>%
      
      arrange(Clocklikeness_ML) %>%
      slice_head(n = argv$clock_ml)
     
  } else if(argv$ml == FALSE & argv$bi == TRUE) {

    filtered_data <- full_table %>%
      filter(

        No_of_taxa >= argv$taxa &
        between(Alignment_length, 
          argv$min_length, 
          argv$max_length) &
        Missing_percent <= argv$missing &
        between(Proportion_parsimony_informative, 
          quantile(Proportion_parsimony_informative, argv$min_pars),
          quantile(Proportion_parsimony_informative, argv$max_pars)) &
        CID <= argv$cid &
        Average_prob_BI >= argv$supp_bi

        ) %>%
      
      arrange(Clocklikeness_BI) %>%
      slice_head(n = argv$clock_bi)

  } else {

    print("Please choose the --ml or --bi flag to apply filtering")

  }

}

write_tsv(filtered_data, file=file.path(argv$output, paste0("crossval_stats_", argv$output, ".tsv")))
write_lines(filtered_data$Locus, file=file.path(argv$output, "crossval_whitelist.txt"))


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

ggsave(file.path(argv$output, "crossval_stats.pdf"), width = 70, height = 30, units = "cm")
