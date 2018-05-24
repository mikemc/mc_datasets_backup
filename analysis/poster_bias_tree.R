devtools::dev_mode(on = TRUE)
library(phyloseq)
library(tidyverse)
library(magrittr)
library(ggtree)
# library(ggthemes)

## Source bias estimation functions
devtools::load_all("~/active_research/src/metacal")

## Avoid confusing behavior from factor conversion
options(stringsAsFactors = FALSE)

## Load data
git_path <- file.path("~/active_research/metagenomics_calibration")
data_path <- file.path(git_path, "costea2017", "data")
ps <- readRDS(file.path(data_path, "2018-05-22_phase3_phyloseq_reads.Rds"))
ps.org <- readRDS(file.path(data_path, "phase3_phyloseq_with_tree.Rds"))
# Get rid of the samples of Individual = A or B b/c these were only extracted
# by one protocol
ps <- subset_samples(ps, !(Individual %in% c("A", "B")))
ps
ps.org <- subset_samples(ps.org, !(Individual %in% c("A", "B")))
ps.org
# Species in the mock community
mock_taxa <- c("Blautia_hansenii", "Clostridium_difficile",
    "Clostridium_perfringens", "Clostridium_saccharolyticum",
    "Fusobacterium_nucleatum", "Lactobacillus_plantarum",
    "Prevotella_melaninogenica", "Salmonella_enterica",
    "Vibrio_cholerae", "Yersinia_pseudotuberculosis")
# denom for clr
denom <- setdiff(mock_taxa, "Blautia_hansenii")

## Filtering
# For each Individual there are a large number of taxa that are found by all
# three protocols, and most taxa within an individual are found by all three.
# When working within a single Individual, I'll keep all taxa that meet this
# threshold by either of the three protocols.

## Get biases by pooling across Individuals
# Goal is to get bias estimates for as many taxa as possible.

# Compute bias posterior for a single observation from each protocol
bias_post <- function (obs1, obs2, n = 128, denom = "all") {
    p1.samples <- extraDistr::rdirichlet(n, 1 + obs1)
    p2.samples <- extraDistr::rdirichlet(n, 1 + obs2)
    colnames(p1.samples) <- names(obs1)
    colnames(p2.samples) <- names(obs2)
    clr1.samples <- p1.samples %>%
        apply(1, clr, denom = denom) %>%
        t
    clr2.samples <- p2.samples %>%
        apply(1, clr, denom = denom) %>%
        t
    clr1.samples - clr2.samples
}

# Function for looping over all individuals
bias_post_ind <- function(ind, protocols, ps, n = 128) {
    # Subset
    sns <- sample_data(ps) %>% as_tibble %>%
        filter(Individual == ind, Protocol %in% protocols) %>%
        arrange(Protocol) %$%
        Sample
    st <- otu_table(ps)[sns, ]
    # Filter
    st.ra <- transform_sample_counts(st, function (x) x / sum(x))
    taxa_pass <- filter_taxa(st.ra, function (x) any(x >= (500/1e6)))
    taxa_pass[denom] <- TRUE # Make sure the mock makes it (should always)
    st.filt <- st[, taxa_pass]
    # Get posterior
    obs1 <- st.filt[1,] %>% c
    obs2 <- st.filt[2,] %>% c
    names(obs1) <- taxa_names(st.filt)
    names(obs2) <- taxa_names(st.filt)
    bias.post <- bias_post(obs1, obs2, n = n, denom = denom)
    # Tidy
    tb <- bias.post %>%
        as_tibble %>%
        add_column(Individual = ind, .before = 1) %>%
        gather(key = "Taxon", value = "Bias", -Individual)
    tb
}

# Compute for each individual and them combine
individuals <- as.character(1:8)
biasHW.post <- individuals %>%
    map(bias_post_ind, c("H", "W"), ps) %>%
    bind_rows
biasQW.post <- individuals %>%
    map(bias_post_ind, c("Q", "W"), ps) %>%
    bind_rows
bias.post <- bind_rows(list(HW = biasHW.post, QW = biasQW.post), 
    .id = "Protocols")
bias <- bias.post %>%
    group_by(Protocols, Taxon) %>%
    summarize(Median = median(Bias), Mean = mean(Bias), SD = sd(Bias)) %>%
    mutate(Median = Median - mean(Median), Mean = Mean - mean(Mean)) %>%
    ungroup()

# Get tibbles for joining to tree
tax <- tax_table(ps) %>% as_tibble
bias <- bias %>%
    left_join(tax, by = "Taxon")
bias.HW <- bias %>% filter(Protocols == "HW") %>% select(-Protocols)
bias.QW <- bias %>% filter(Protocols == "QW") %>% select(-Protocols)


## Prune the tree
# MRCA node for the lower half:
two_taxa <- tax %>% 
    filter(Genus %in% c("Bacteroides", "Bifidobacterium")) %>%
    group_by(Genus) %>%
    sample_n(1) %$%
    Taxon
# Prune to the lower half and just the taxa with bias estimates
tree <- phy_tree(ps)
mrca_node <- ape::getMRCA(tree, two_taxa)
tree0 <- tree %>%
     ape::extract.clade(mrca_node) %>%
     prune_taxa(bias.QW$Taxon, .)
# Remove taxa with highly uncertain bias estimates
high_sd_taxa <- bias.QW %>%
    filter(SD > 1.2) %$%
    Taxon
tree1 <- ape::drop.tip(tree0, high_sd_taxa)

# Plot with ggtree
ttree <- tidytree::as_data_frame(tree1) %>%
    left_join(bias.QW, by = c("label" = "Taxon"))
td <- tidytree::as.treedata(ttree)
bias_mean <- mean(ttree$Mean, na.rm = T)
bias_range <- range(ttree$Mean, na.rm = T) - bias_mean

# Get Genus MRCAs
# taxa <- ttree %>% filter(Genus == "Bacteroides") %$% label
# bacteroides <- ape::getMRCA(td@phylo, taxa)
genera <- list("Bacteroides", "Prevotella", "Bifidobacterium", "Alistipes",
    "Parabacteroides", )
mrcas <- taxa %>%
    map(~ filter(ttree, Genus == .) %$% label %>% ape::getMRCA(td@phylo, .))
names(mrcas) <- taxa


# gt <- ggtree(td, layout = "slanted")
# gt <- ggtree(td, layout = "rectangular", ladderize = F)
# gt <- ggtree(td, layout = "fan", ladderize = T)
# gt <- ggtree(td, layout = "circular", ladderize = T)
line_size <- 1.0
gt <- ggtree(td, layout = "rectangular", ladderize = T, size = line_size) + 
    scale_x_reverse() + coord_flip()
# gt + geom_tiplab(aes(label = node))
gt + geom_tiplab(aes(label = label), angle = -90) + geom_nodelab(aes(label = node))

gt + geom_tippoint(aes(fill = Mean - bias_mean), 
    shape=21, color="black", size = 4) +
    scale_fill_distiller(type = "seq", palette = 1, 
        limits = bias_range %>% round(digits=1),
        # breaks = sort(c(bias_range, -2, 0, 2)) %>% round(digits=1)) +
        breaks = sort(c(bias_range, 0)) %>% round(digits=1)) +
    geom_strip(37, 36, "Bacteroides", align = F, offset = -0.08, 
        hjust = 0.5, offset.text = -0.15, barsize = line_size,
        family = "Cantarell", fontsize=8.4375) +
    geom_cladelabel(mrcas[["Prevotella"]], "Prevotella", 
        offset = -0.5, 
        hjust = 0.5, offset.text = -0.15, barsize = line_size,
        family = "Cantarell", fontsize=8.4375) +
    geom_cladelabel(mrcas[["Bifidobacterium"]], "Bifidobacterium", 
        offset = -0.25, 
        hjust = 0.5, offset.text = -0.15, barsize = line_size,
        family = "Cantarell", fontsize=8.4375) +
    # geom_cladelabel(mrcas[["Parabacteroides"]], "Parabacteroides", 
    #     offset = -0.25, 
    #     hjust = 0.5, offset.text = -0.04) +
    geom_cladelabel(mrcas[["Alistipes"]], "Alistipes", 
        offset = -0.74, 
        hjust = 0.5, offset.text = -0.15, barsize = line_size,
        family = "Cantarell", fontsize=8.4375) +
    geom_cladelabel(104, "Enterobacteriaceae", offset = -0.1, 
        hjust = 0.5, offset.text = -0.15, barsize = line_size,
        family = "Cantarell", fontsize=8.4375) +
    labs(fill = "Bias") +
    theme(text = element_text(family = "Cantarell"),
        legend.position=c(0.9, 0.75),
        legend.direction = "vertical", 
        legend.title = element_text(size=24*4/3),
        legend.text = element_text(size=18*4/3),
        legend.key.height = unit(3, "line"))
ggsave("/tmp/bias_tree_QW.svg", width = 13, height = 7, units = "in", scale = 4/3)
# ggsave("/tmp/bias_tree_QW0.svg", width = 13*4/3, height = 7*4/3, units = "in")

# fontsize specified to ggtree gets multiplied by 2.133333

# theme_update(text = element_text(family = "Cantarell", size=18))
my_theme <- theme(text = element_text(family = "Cantarell"),
    legend.position=c(0.9, 0.85),
    legend.direction = "vertical", 
    legend.title = element_text(size=24),
    legend.text = element_text(size=18))

###########
gt <- ggtree(td, layout = "rectangular")
g <- gt + geom_tippoint(aes(fill = Mean - mean(Mean, na.rm = T)), 
    shape=21, color="black", size = 3) +
    scale_fill_distiller(type = "seq", palette = 1) +
    theme(legend.position='right') +
    labs(fill = "Bias")
ggsave("/tmp/bias_tree_QW.pdf", g, width = 6, height = 8, units = "in")
ggsave("/tmp/bias_tree_QW0.svg", g, width = 6, height = 8, units = "in")
ggsave("/tmp/bias_tree_QW1.svg", g, width = 6, height = 8, units = "in", scale = 4/3)
ggsave("/tmp/bias_tree_QW2.svg", g, width = 6, height = 8, units = "in",
    device = grDevices::svg)
