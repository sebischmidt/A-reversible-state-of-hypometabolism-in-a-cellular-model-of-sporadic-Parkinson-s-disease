######Fisher's exact test######
test=function(r,n,R,N)
{
  spalte1=c(r,n-r)
  spalte2=c(R-r,N-R-n+r)
  as.double(fisher.test(cbind(spalte1,spalte2))$p.value)
}
#N all proteins
#R DEPs
#n e.g. GLI3 target genes present in all proteins
#r overlap n and R


GLI3_targets <- test(r = 528,n = 2321,R = 1667,N = 7943); GLI3_targets


#####Venn Diagram#####
setwd("")

library(extrafont)
font_import()
library(VennDiagram)

##Function to add thousands separators
myvenn <- VennDiagram::draw.pairwise.venn
body(myvenn)[[46]][[3]] <- quote(function(x) {
  prettyNum(x ,big.mark=",",scientific=FALSE)
})

##Plot Venn diagram exemplary for the overlap of predicted GLI3 targets with DEPs
emf(file = "GLI3targets.emf")
pdf("GLI3targets.pdf")
myvenn(2321, #predicted GLI3 targets expressed in e.g. NSC1a (full gene list)
       1667, #DEGs
       528, #Overlap
       category = c("GLI3 targets", "DEPs"),
       scaled = TRUE, lwd = rep(6, 2), lty = rep("solid", 2), 
       col = c("orange", "red"), 
       fill = c("orange", "red"), alpha = rep(0.35, 2), 
       label.col = rep("black", 3), cex = rep(4, 3), 
       fontface = rep("bold", 3), 
       cat.pos = c(180, 180), cat.dist = rep(0.05, 2), 
       cat.cex =  rep(4, 2), cat.col = rep("black", 2), 
       cat.fontface = rep("bold", 2),
)
dev.off()
