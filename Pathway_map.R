library(dplyr)
library(ggplot2)
library(readr)
library(stringr)
library(clusterProfiler)
library(mygene)
library(VennDiagram)
library(ggraph)
library(RColorBrewer)
library(igraph)
source("Tools/utilities.R")

##### Load in the data and correct the protein names #####

# Add in whtever protein names you wish to act as candidates for the protein maps.
sig_names = tibble(Original_name = c("APC", "A6r", "IGF1R", "Arg", "CREB1", "CaMK1d", "CaMK2a", "CTNNB1",
                                      "PRKD1", "ACLY", "MEK1", "Kit", "GSK3a", "WASP", "PFKFB3", "ACP1", 
                                      "Smad2", "ZAP70", "VIM", "VEGFR1","Tau", "SMG1", "SRF" ,"Paxillin 1"))

# Load in the corrected names for all kinexus antibodies. 
proper_names = read_csv("All_names/all_names.csv")

# Correct the names and add Entrez IDs 
df = get_entrez_ids(sig_names, proper_names)

# Load in the background Entrez IDs. 
All_names = read.csv("All_names/phospho_background.csv")

##### Run ORA on all KEGG pathways using cluster profiler #####

# Use clusterprofiler to run ORA on all pathways
kegg_res = enrichKEGG(as.character(df$entrezgene), universe = as.character(All_names$entrezgene), organism = "hsa",
                      keyType = "kegg",pvalueCutoff = 0.99, qvalueCutoff = 0.99, pAdjustMethod = "BH",
                      minGSSize = 1, maxGSSize = 50000)

# Load up the names of some pathways of interest. 
pathways = read_csv("Pathway_lists/relevant_paths_KEGG.csv", col_names = FALSE)
colnames(pathways) = c("Description", "ID", "type")

# Switch from entrez ids to readable names. 
kegg_res <- setReadable(kegg_res, 'org.Hs.eg.db', 'ENTREZID')

# Select the pathways of interest from the results. 
kegg_subset = kegg_res@result[kegg_res@result$Description %in% pathways$Description,]
kegg_subset = subset(kegg_subset, select = c(Description, geneID, Count))

##### Rearrange results for plotting ####

# Create one description for each count. 
Description = c()
for (i in 1:length(kegg_subset$Description)) {
  Description = c(Description, rep(kegg_subset$Description[i], kegg_subset$Count[i]))
}

# Get all of the gene IDs by splitting up the gene/gene into the individual components.
geneID = unlist(str_split(kegg_subset$geneID, "/"))

# Match the pathway description to the gene ID. 
to_plot = tibble(Description, geneID)

# Join the relevant paths to the Rearranged descriptions and genes. 
to_plot = inner_join(to_plot, pathways, by = c("Description"))

##### Create the venn diagram ##### 

# suppress log files 
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

# Extract sets for each dataframe. 
phago = to_plot[to_plot$type == "Phagocytosis related",]$geneID
apop = to_plot[to_plot$type == "Apoptosis related",]$geneID
gsk3 = to_plot[to_plot$type == "GSK3 Related",]$geneID

# Define the colors
myCol <- brewer.pal(3, "Pastel2")

# Create the diagram.
x = venn.diagram(
  x = list(phago, apop, gsk3),
  category.names = c("Phagocytosis" , "Apoptosis" , "GSK3"),
  filename = NULL,
  
  # Circles
  lwd = 2,
  fill = myCol,
)

# Draw the diagram
grid.newpage()
grid.draw(x)

# Initiate writing to PDF file
pdf("test.pdf", height = 4, width = 6, paper = "letter")

# Draw the graphical diagram into the PDF 
grid.draw(x) 

# Stop writing to the PDF file
dev.off()

##### Create the network diagram #####

# Create the edges
edges <- distinct(subset(to_plot, select = c("geneID", "type")))
edges <- transform(edges, x = as.numeric(factor(type)))
edges <- tibble(transform(edges, y = as.numeric(factor(geneID))+3))
edges$type = factor(edges$type)

# Create the nodes
temp1 = distinct(subset(edges, select = c(type, x)))
temp2 = distinct(subset(edges, select = c(geneID, y)))
colnames(temp1) = c("name", "id")
colnames(temp2) = c("name", "id")
nodes = rbind(temp1, temp2)
n = nrow(nodes)

# degree of nodes (number of ties for each dolphin)
tb = tibble(v = c(1:n, edges$x, edges$y))
d = count(tb, v)$n - 1
nodes = mutate(nodes, degree = d)
nodes$color = 1
nodes$color[nodes$degree > 5] = 2

# create graph from data frames
g = graph_from_data_frame(edges, directed = FALSE, nodes)

lay = create_layout(g, layout = "linear", circular = TRUE)

# Select your preferred pallete: 
cols_f <- colorRampPalette(RColorBrewer::brewer.pal(8, 'Pastel2'))
cols_y <- colorRampPalette(RColorBrewer::brewer.pal(9, 'Spectral'))

# plot with ggraph
ggraph(lay, circular = TRUE) + 
  geom_edge_arc(aes(color = as.factor(edges$type))) +  
  geom_node_point(aes(size = degree), colour = as.factor(nodes$color)) +
  geom_node_text(aes(label = name), repel=TRUE) +
  theme_graph()

# Save the plot. 
ggsave("network.svg", height = 6, width = 5)


