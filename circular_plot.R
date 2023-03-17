library(ape)
library(ggplot2)
library(ggtree)
library(readxl)
library(ggnewscale)
library(tidytree)


fasta.file.name <- "core_genome_alignment.fasta"
newick.file.name <- "core_ML_tree.nhx"
metadata.file.name <- "tree_metadata.xlsx"

# METADATA ------------------------------------------------------

all.df <- read_excel(metadata.file.name)
#all.df <- as.data.frame(all.df)
colnames(all.df)[1] <- "id"
rownames(all.df) <- all.df$id

all.df$rUTI <- as.character(all.df$rUTI)
all.df$MDR <- as.character(all.df$MDR)

seqtyp.df <- data.frame('id' = all.df[,c('id')], 'ST' = all.df[,c('ST')])

stnode.df <- data.frame(
  node = c(386, 254, 228, 291, 367, 337, 352, 360),
  annot = c("ST95", "ST73", "ST69", "ST131", "ST12", "ST404", "ST1193", "ST127")
)

source.df <- data.frame("id" = all.df[,c("id")], "Source" = all.df[,c("Source")])
rownames(source.df) <- source.df$id
source.df <- subset(source.df, select = -c(id))

hospital.df <- data.frame("id" = all.df[,c("id")], "Hospital" = all.df[,c("Hospital")])
rownames(hospital.df) <- hospital.df$id
hospital.df <- subset(hospital.df, select = -c(id))

mdr.df <- data.frame("id" = all.df[,c("id")], "MDR" = all.df[,c("MDR")])
rownames(mdr.df) <- mdr.df$id
mdr.df <- subset(mdr.df, select = -c(id))

abx.df <- all.df[,c("id", "AMP", "AUG", "CLX", "CIP", "GEN", "NIT", "TRI")]
abx.df <- as.data.frame(abx.df)
rownames(abx.df) <- abx.df$id
abx.df <- subset(abx.df, select = -c(id))

genes.df <- all.df[,c("id", "aac(6')-Ib-cr5", "blaCTX-M-1", "blaCTX-M-15", "blaOXA-1")]
genes.df <- as.data.frame(genes.df)
rownames(genes.df) <- genes.df$id
genes.df <- subset(genes.df, select = -c(id))
genes.df[genes.df==0] <- NA
genes.df[] <- lapply(genes.df, factor)

bla.df <- data.frame("id" = all.df[,c("id")], "bla" = all.df[,c("bla")])
rownames(bla.df) <- bla.df$id
bla.df <- subset(bla.df, select = -c(id))

dfrA.df <- data.frame("id" = all.df[,c("id")], "dfrA" = all.df[,c("dfrA")])
rownames(dfrA.df) <- dfrA.df$id
dfrA.df <- subset(dfrA.df, select = -c(id))


phylo.dat <- data.frame(id=c(221,250,403,406,407,425,430,217,218), Phylogroup=c('D','B2','E','C','B1','A','F','G','G'))



# TREE PLOTTING -------------------------------------------------

iqtree <- ape::read.tree(newick.file.name)
gg = ggtree(iqtree, size=0.5, layout="fan", open.angle=30) %<+% all.df +
  #Colour phylogroups
  geom_hilight(
    data=phylo.dat,
    mapping=aes(node=id, fill=Phylogroup),
    alpha=.4) +
  scale_fill_discrete(
    guide=guide_legend(ncol=2, order=1)
  )+

  new_scale_fill()+
  new_scale_color() +

  # tip label rUTI cases
  geom_tippoint(
    aes(shape = rUTI, colour = rUTI),
    size=4,
    na.rm=TRUE,
    alpha=0.6) +
  scale_shape_manual(
    values=c("1" = 18, "2" = 17, "3" = 15, "4" = 19),
    na.translate=F,
    na.value=NA,
    name="rUTI\ncases",
    guide=guide_legend(order=2)) +
  scale_colour_manual(
    name="rUTI\ncases",
    values=c("1" = "deepskyblue3", "2" = "cornflowerblue", "3" = "firebrick", "4" = "darkorchid2"),
    na.value=NA,
    na.translate = F,
    guide=guide_legend(order=2)
  ) +

# Label reference
new_scale_color() + new_scale_fill() +
geom_point2(
  aes(subset=(node==165)),
  shape=16, size=4, color="darkgreen", alpha=0.7) +
geom_cladelab(node=165, label="Ref.", fontsize=4, offset.text=.01, angle=220) +
new_scale_color()

# Heatmap of source metadata
gg <- gheatmap(
    gg, 
    source.df, 
    offset=0.0005, 
    width=0.05, 
    colnames_position = 'top', 
    legend_title = "Source",
    colnames_offset_y = 0.5,
    colnames_angle = 60,
    font.size = 4,
    color=NA,
    hjust=0
  ) + 
  scale_fill_manual(
    values=c("Community" = "skyblue1", "Hospital" = "royalblue1", "NA" = "white"),
    na.value="white",
    name="Source",
    guide=guide_legend(order=3)
  ) + 
  new_scale_fill()

# Heatmap of hospital metadata
gg <- gheatmap(
  gg,
  hospital.df,
  offset=0.02,
  width=0.05,
  colnames_position = 'top', 
  legend_title = "Hospital",
  colnames_offset_y = 0.5,
  colnames_angle = 60,
  font.size = 4,
  color=NA,
  hjust=0
) + 
  scale_fill_discrete(
    name="Hospital",
    guide=guide_legend(ncol=2, order=4), 
    labels=c('JPUH', 'NCH', 'NNUH', 'NSFT', 'QEH'),
    na.value="white",
    na.translate=F
    ) +
  new_scale_fill()

# Heatmap of MDR
gg <- gheatmap(
  gg,
  mdr.df,
  offset=0.05,
  width=0.05,
  colnames_position = 'top',
  colnames_offset_y = 0.5,
  colnames_angle = 60,
  font.size = 4,
  color=NA,
  hjust=0
  ) +
  scale_fill_manual(
    values = c("1"="darkred"),
    name = "",
    labels = c("MDR"),
    na.value = "white",
    guide=guide_legend(title=NULL, order=5)
  ) + new_scale_fill()

# Heatmap for phenotypic resistance
gg <- gheatmap(
    gg, 
    abx.df, 
    offset=0.065, 
    width=0.5, 
    colnames_position = 'top',
    colnames_offset_y = 0.5,
    colnames_angle = 60,
    font.size = 4,
    color="gray75",
    hjust=0
  ) +
  scale_fill_manual(
    values=c(
      "R" = "darksalmon",
      "S" = "white"
    ),
    na.value = "white",
    name = "Antibiotic\nsusceptibility",
    labels = c('Resistant', 'Susceptible'),
    guide=guide_legend(order=6)
  ) + new_scale_fill()
  
# Heatmap for AMR determinants presences
gg <- gheatmap(
  gg,
  genes.df,
  offset=0.22,
  width=0.3,
  colnames_position = 'top',
  colnames_offset_y = 0.5,
  colnames_angle = 60,
  font.size = 4,
  color="gray75",
  hjust=0
  ) + 
  scale_fill_manual(
    values=c("1" = "tan2"),
    na.value="white",
    name = "AMR gene",
    labels = c("Present"),
    guide=guide_legend(order=7)
  )+
  new_scale_fill()

# Heatmap of Beta-lactamase presence
gg <- gheatmap(
  gg,
  bla.df,
  offset=0.32,
  width=0.05,
  colnames_position = 'top',
  colnames_offset_y = 0.5,
  colnames_angle = 60,
  font.size = 4,
  color=NA,
  hjust=0
  ) +
  scale_fill_discrete(
    name=expression(italic(bla)~genes),
    labels=c(
      expression(italic(bla)[CARB-2]),
      expression(italic(bla)[DHA-1]),
      expression(italic(bla)[SHV-1]),
      expression(italic(bla)[TEM]),
      expression(italic(bla)[TEM-1]),
      expression(italic(bla)[TEM-148]),
      expression(italic(bla)[TEM-3]),
      expression(italic(bla)[TEM-32]),
      expression(italic(bla)[TEM-34]),
      expression(italic(bla)[TEM-4])
    ),
    na.value="white",
    guide=guide_legend(label.hjust = 0, ncol=2, order=8)
  ) + new_scale_fill()

# Heatmap of dfrA presence
gg <- gheatmap(
  gg,
  dfrA.df,
  offset=0.341,
  width=0.05,
  colnames_position = 'top',
  colnames_offset_y = 0.5,
  colnames_angle = 60,
  font.size = 4,
  color=NA,
  hjust=0
  ) + 
  scale_fill_discrete(
    name=expression(italic(dfrA)~genes),
    labels=c(
      expression(italic(dfrA1)),
      expression(italic(dfrA1/dfrA12)),
      expression(italic(dfrA1/dfrA14)),
      expression(italic(dfrA12)),
      expression(italic(dfrA14)),
      expression(italic(dfrA15)),
      expression(italic(dfrA17)),
      expression(italic(dfrA36)),
      expression(italic(dfrA5)),
      expression(italic(dfrA7)),
      expression(italic(dfrA8))
    ),
    na.value="white",
    guide=guide_legend(label.hjust = 0, ncol=2, order=9)
    )+ new_scale_fill()

# Labelling sequence types
 gg <- gg + geom_cladelab(
    data = stnode.df,
    mapping = aes(
      node = node,
      label = annot
    ),
    align = TRUE,
    show.legend = FALSE,
    horizontal=FALSE,
    hjust="center",
    angle="auto",
    fontsize = 4,
    offset=0.37,
    offset.text=0.03,
    barsize=0.8
  ) + new_scale_fill() +
  
  geom_hilight(
    data=stnode.df, 
    mapping=aes(node=node), 
    extendto=0.66, 
    alpha=.2, 
    fill='grey', 
    color='grey50', 
    size=0.05
  ) + 
  theme(
    legend.text = element_text(size=9),
    legend.title = element_text(size=11),
    legend.box="vertical",
    legend.direction="vertical",
    legend.background = element_rect(),
    legend.key.size = unit(0.3, 'cm')
  ) 


gg

# ggsave(
#   "extended_circular_phylogenomic_tree.png",
#   last_plot(),
#   device="png",
#   width=14,
#   height=12,
# )
