# 1. Install and load the libraries
if(!require(ape)) install.packages("ape")
if(!require(phytools)) install.packages("phytools")
library(ape)
library(phytools)

# 2. Load your exact tree files
tree1 <- read.tree("./trees/IDXUT_Ace2_matching_mafft_cbrc_jp.txt")
tree2 <- read.tree("./trees/IDXUT_Cqm1_matching_mafft_cbrc_jp.txt")


# 3. Create a clean 1-to-1 association matrix
# This automatically maps matching names and leaves the single haplotypes visible but unlinked

# Strip the trailing "_H1" or "_H2" suffix to find the baseline mosquito sample IDs
base_names_1 <- gsub("_H[0-9]+$", "", tree1$tip.label) 
base_names_2 <- gsub("_H[0-9]+$", "", tree2$tip.label) 

smart_links <- list()

for (i in 1:length(tree1$tip.label)) {
  # Find which tips in Tree 2 share the same biological sample base name
  matching_targets <- tree2$tip.label[base_names_2 == base_names_1[i]]
  
  # Connect every match found (allows both H1 and H2 to connect to a single H1)
  for (target in matching_targets) {
    smart_links[[length(smart_links) + 1]] <- c(tree1$tip.label[i], target)
  }
}

# Combine the list into the structured matrix required by phytools
association <- do.call(rbind, smart_links)
colnames(association) <- c("Tree1_Tips", "Tree2_Tips")

# 4. Test without rotation for a quick check
cophylo_obj_test <- cophylo(tree1, tree2, assoc=association, rotate=FALSE)

plot(cophylo_obj_test, 
     link.type="straight", 
     link.col=lty <- "darkblue", 
     link.lwd=1, 
     fsize=0.3) # Tiny text size so 1,700 tips don't overlap


pdf("./figures/ace2-cqm1-tanglegram_test.pdf", width=14, height=35)
plot(cophylo_obj_test, 
     link.type="straight", 
     link.col=lty <- "darkblue", 
     link.lwd=1, 
     fsize=0.3) # Tiny text size so 1,700 tips don't overlap
dev.off()

# 5. Run with optimizing node rotations for your 1,700 lines automatically
cophylo_obj <- cophylo(tree1, tree2, assoc=association)

plot(cophylo_obj, 
     link.type="straight", 
     link.col=lty <- "darkblue", 
     link.lwd=1, 
     fsize=0.3) # Tiny text size so 1,700 tips don't overlap


pdf("./figures/ace2-cqm1-tanglegram.pdf", width=14, height=35)
plot(cophylo_obj, 
     link.type="straight", 
     link.col=lty <- "darkblue", 
     link.lwd=1, 
     fsize=0.3) # Tiny text size so 1,700 tips don't overlap
dev.off()


# 6. Run with optimizing node rotations after Rooting by a chosen outgroup sample name present in your files
tree1_rooted <- root(tree1, outgroup = "LSP_24_113_1", resolve.root = TRUE)
tree2_rooted <- root(tree2, outgroup = "LSP_24_113_1", resolve.root = TRUE)

cophylo_obj_rooted <- cophylo(tree1_rooted, tree2_rooted, assoc=association)

plot(cophylo_obj_rooted, 
     link.type="straight", 
     link.col=lty <- "darkblue", 
     link.lwd=1, 
     fsize=0.3) # Tiny text size so 1,700 tips don't overlap


pdf("./figures/ace2-cqm1-tanglegram_rooted.pdf", width=14, height=35)
plot(cophylo_obj_rooted, 
     link.type="straight", 
     link.col=lty <- "darkblue", 
     link.lwd=1, 
     fsize=0.3) # Tiny text size so 1,700 tips don't overlap
dev.off()