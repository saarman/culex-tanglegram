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
association <- cbind(tree1$tip.label, tree1$tip.label)

# 4. Run the high-performance cophylo alignment without rotation for a quick check
cophylo_obj_test <- cophylo(tree1, tree2, assoc=association, rotate=FALSE)

# This optimizes node rotations for your 1,700 lines automatically
cophylo_obj <- cophylo(tree1, tree2, assoc=association)

# 5. Plot
plot(cophylo_obj_test, 
     link.type="straight", 
     link.col=lty <- "darkblue", 
     link.lwd=1, 
     fsize=0.3) # Tiny text size so 1,700 tips don't overlap

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
