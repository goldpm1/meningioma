library (ggtree)
library (ggplot2)
library (ape)

args <- commandArgs(trailingOnly = TRUE)
INPUT_PATH = args[1]
OUTPUT_PATH = args[2]


tree = ape::read.tree ( INPUT_PATH )

#xlim(NA, 200) + ylim(NA, 10)
fig1 = ggtree( tree, size = 1  ) + 
  theme_tree2 ( ) +
  # geom_tippoint( colour = "black", size = 5, shape = 15 ) +
  geom_tiplab( size = 4, color="black", arrow = TRUE, offset = 1 ) +
  labs( caption = "#Mutation", size = 40, fontweight = "bold") 
print (fig1)
ggsave( OUTPUT_PATH, dpi = "retina", width = 24, height = 10, unit = "cm" )

#ggsave( "/Users/youngsoo/Downloads/190426_dendrogram.pdf", dpi = "retina", width = 24, height = 10, unit = "cm" )