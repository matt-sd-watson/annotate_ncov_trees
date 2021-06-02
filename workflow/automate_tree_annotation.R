library(tidyverse)
#install.packages("ggtree")
library(ggtree)
# install.packages("data.table")
library(data.table)
# install.packages("phytools")
library(phytools)
library(ape)
library(grid)
library(stringr)
library(Biostrings)
library(TraMineR)
library(ggplot2)
library(argparse)
library(lubridate)

## Inputs: 
## - a tree in Newick format generated through iqtree
## - A lineage ID corresponding to the Pango Lineage of the tree samples
## A qc_pass list with the names of samples that pass a QC threshold as "WGS_Id"

## Outputs: 
## - A customized annotated tree in pdf format for each the lineage samples

parser <- ArgumentParser(description='Annotate a phylogenetic tree generated through iqtree with custom labelling')
parser$add_argument('--input_tree', type = "character", 
                    help='path to the tree input')
parser$add_argument('--output_directory', type = "character", 
                    help='output directory for annotated tree in a PDF report')
parser$add_argument('--lineage_id', type = "character", 
                    help='Pango lineage identifier for file output annotation')
parser$add_argument('--qc_list', type = "character", 
                    help='Metadata for samples that pass a certain QC threshold')

args <- parser$parse_args()

tree <- read.tree(args$input_tree)
qc_pass <- read.table(args$qc_list, header = T, sep = ',', fill = TRUE)
qc_pass <- qc_pass[!is.na(qc_pass$Login_Date) && qc_pass$Login_Date != "",]
qc_pass <- qc_pass[nchar(as.character(qc_pass$Login_Date)) == 9,]

if (nrow(qc_pass) == 0) {
  warning("There are no rows for the metadata. Please double check the threshold
          values and metadata entries.")
}

tree <- drop.tip(tree, "MN908947")

## Create a dataframe of the tree
tr.df <- fortify(tree)
## Create a list of tree labels
tr.df.labs <- as.data.frame(tr.df) %>%
  filter(isTip == "TRUE") %>%
  select(label)

# if the sample is not on the qc list, drop it from the tree
# drop only the PHO sequences that do not pass thresholds
# only_pho <- data.frame(label = tr.df.labs[grepl("PHLON", tr.df.labs$label),])
# sample_no_pass <- as.data.frame(only_pho[!only_pho$label %in% qc_pass$WGS_Id,])
# colnames(sample_no_pass) <- "label"
# print(args$lineage_id)
# print(sample_no_pass)
# tree <- drop.tip(tree, as.vector(sample_no_pass$label))
# 
# ## Create a dataframe of the tree
# tr.df <- fortify(tree)
# ## Create a list of tree labels
# tr.df.labs <- tr.df %>%
#   filter(isTip == "TRUE") %>%
#   select(label)

## PARKCED CODE: Annotate the 5 most prevalent countries/regions
## IMPORTANT: Requires background Gisaid sequences that have normalized names

# meta.df <- tr.df.labs %>%
#   mutate(label.2 = str_split_fixed(tr.df.labs$label, "/", 4)[,1])
# 
# meta.df$all_label <- ifelse(grepl("PHLON", meta.df$label), "ON-PHO", meta.df$label.2)
# 
# # treat the countries as a factor for plotting
# meta.df$label.2 <- as.factor(meta.df$label.2)
# 
# # find the most abundant countries and regions for annotation
## exclude PHO sequences from being considered in regions
# regions <- as.data.frame(table(meta.df$label.2)) %>% filter(!grepl("PHL", Var1)) %>% arrange(-Freq)
# 
# # create a frame with the 4 most abundant regions and give each one a colour
# cols <- c("#C31130", "blue", "dark green", "light blue", "dark red", "light green")
# names(cols) <- c(as.character("ON-PHO"),as.character(regions[1:5,]$Var1))
# 
# meta.df$nml_lab <- paste("ON-PHL", str_split_fixed(meta.df$label, "PHLON|-SARS", 4)[,2],
#                          str_split_fixed(meta.df$label, "PHLON|-SARS", 4)[,3], sep = "-")

qc_pass$Date <- as.Date(qc_pass$Login_Date, format = "%d%B%Y")

qc_pass <- qc_pass[!is.na(qc_pass$Date),]

date_cat <- c(ymd(Sys.Date() - 7), ymd(Sys.Date()))

qc_pass$date_cat <- ifelse(as.Date(qc_pass$Date) > min(date_cat) &&
                             as.Date(qc_pass$Date) < max(date_cat), "Latest Week", "All Past Weeks")

qc_pass_keep <- subset(qc_pass, select = c(WGS_Id, Date, date_cat))

meta.df <- merge(tr.df.labs, qc_pass_keep, by.x = "label", by.y = "WGS_Id")

# treat the countries as a factor for plotting
meta.df$date_cat <- as.factor(meta.df$date_cat)

# if the lineage is B.1.1.7, set the background tips to a transparent grey
# otherwise, have a darker grey for lineage trees with fewer samples
grey_trans <- rgb(179, 179, 179, max = 255, alpha = 45)

cols <- c()
## Define colours for tree annotation depending on lineage identifier
if (args$lineage_id == "B.1.1.7") {
  cols <- c("Latest Week" = "dark blue",  "All Past Weeks" = grey_trans)
} else {
  cols <- c("Latest Week" = "dark blue",  "All Past Weeks" = "grey60")
}

meta.df$nml_lab <- paste("ON-PHL", str_split_fixed(meta.df$label, "PHLON|-SARS", 4)[,2],
                         str_split_fixed(meta.df$label, "PHLON|-SARS", 4)[,3], sep = "-")

pl.1 <- ggtree(tree, size = 0.25)

pl.2 <- pl.1 %<+% meta.df +
  geom_tippoint(aes(x=x+0.000001, subset=date_cat == "Latest Week", label = label, colour = date_cat),  size = 1.2, shape = 16) +
  geom_tippoint(aes(x=x+0.000001, subset=date_cat == "All Past Weeks", label = label, colour = date_cat),  size = 1.2, shape = 16) +  
  # geom_tiplab(aes(x=x+0.000001, subset=label.2 == "ON-PHO", label = nml_lab),  size = 1.5, offset = 0.000005) +
  # geom_tiplab(aes(x=x+0.000001, subset=label.2 == "non-PHO", label = label),  size = 1.5, offset = 0.000005) +
  theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),
        legend.text=element_text(size=6),
        legend.position = c(0.85, 0.6),
        legend.title=element_text(size=7, face = "bold"),
        legend.box = "vertical",
        legend.box.margin = margin(0.01,0.02,0.02,0.02,"cm"),
        legend.box.background = element_rect(colour = "grey50"),
        legend.spacing.y = unit(0.1, "cm"),
        legend.key.width = unit(0.3, "cm"),
        legend.key.height = unit(0.3, "cm"),
        legend.spacing.x = unit(0.3, 'cm'),
        legend.key.size = unit(3,"lines")) +
  scale_color_manual("Location", values=cols) +
  guides(color = guide_legend(override.aes = list(size = 1.75))) +
  ggtitle(paste("Phylo tree", args$lineage_id, format(Sys.Date(), "%B %e, %Y"),
                sep=", "))

setwd(args$output_directory)

pl.2
## SAVE annotated tree

plot_name <- paste(paste(args$lineage_id, "_phylo_tree", sep = ""), ".pdf", sep="")

ggsave(filename = plot_name, plot = pl.2, width = 11, height = 15)

file.remove("Rplots.pdf")

dev.off()
rm(list = ls())








