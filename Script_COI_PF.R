# Script written by R. Morard and Ana Carolina Bercini Gusmao
# Load libraries
library(tidyverse)
library("ape")
library(openxlsx)
library(scales)
library(ggpubr)
library(ggpol)
library(ggpmisc)
library(patchwork)
library(broom)

# set working directory with the script and the supplement to the Article 
#setwd("C:/Working directory/Recherche/Writing/Manuscripts/2024/IN PROGRESS - COI PF/Code")

# Load data files
# Load supplementary material files from the article "Exploring the potential of the COI gene marker for DNA barcoding of planktonic foraminifera." available at Zenodo (doi: 10.5281/zenodo.14285737)
COI_Barcoding <- read.xlsx(xlsxFile = "Supplementary material 1.xlsx", sheet = 1)
qPCR_data     <- read.xlsx(xlsxFile = "Supplementary material 2.xlsx", sheet = 1)
LBF_data      <- read.xlsx(xlsxFile = "Supplementary material 3.xlsx", sheet = 1)

# Load the Table S1 from the article "" available at Zenodo (https://zenodo.org/records/10454474)
SSU_Barcoding <- read.xlsx(xlsxFile = "Table S1.xlsx", sheet = 1,startRow=2)

# Select unique sequences pattern among SSU and COI data
Selection_COI <- COI_Barcoding %>% 
  unite("Morphotaxa", Morphogroup, Genus, Species, sep = "|") %>% 
  select(Morphotaxa, Sequence_COI) %>%
  distinct() %>% 
  mutate(seq_ID = sprintf("seq%03d", row_number())) %>%
  unite("ID", seq_ID, Morphotaxa, sep = "|") %>% 
  rename("Sequence" = "Sequence_COI")
  
Selection_SSU <- SSU_Barcoding %>% 
  unite("Morphotaxa", Clade_curated, Genus_curated, species_curated, sep = "|") %>% 
  unite("Sequence", `37F_seq`:`49E_seq`, remove = FALSE,sep= "") %>% 
  filter(Quality == "BASETYPE") %>% 
  select(Morphotaxa, Sequence) %>% 
  distinct() %>% 
  mutate(seq_ID = sprintf("seq%03d", row_number())) %>%
  unite("ID", seq_ID, Morphotaxa, sep = "|") 

# Code for a function to write FASTA file
# Piece of code taken from: https://bootstrappers.umassmed.edu/guides/main/r_writeFasta.html
writeFasta<-function(data, filename){
  fastaLines = c()
  for (rowNum in 1:nrow(data)){
    fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,"ID"], sep = "")))
    fastaLines = c(fastaLines,as.character(data[rowNum,"Sequence"]))
  }
  fileConn<-file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}

# write fasta for DADA2
writeFasta(Selection_COI, "COI_selection.fasta")  
writeFasta(Selection_SSU, "SSU_selection.fasta")

# The two FASTA files generated must be align at MAFFT (https://mafft.cbrc.jp/alignment/server/) and named as below for the analyses on patristic distances to work

# Import ALN files after aligning them automatically with MAFFT.
COI_ALN <- read.dna("COI_ALN.fasta", "fasta")
SSU_ALN <- read.dna("SSU_ALN.fasta", "fasta")

# Produce a 3 column file with SEQ-A, SEQ-B and the distance
dist.SSU <- dist.dna(SSU_ALN, model = "K80", pairwise.deletion = FALSE) %>% 
            as.matrix() %>% 
            as.data.frame() %>% 
            rownames_to_column(var = "SEQA") %>% 
            pivot_longer(-SEQA, names_to = "SEQB", values_to = "distance") %>%
            filter(SEQA > SEQB)  

dist.COI <- dist.dna(COI_ALN, model = "K80", pairwise.deletion = FALSE) %>% 
            as.matrix() %>% 
            as.data.frame() %>% 
            rownames_to_column(var = "SEQA") %>% 
            pivot_longer(-SEQA, names_to = "SEQB", values_to = "distance") %>%
            filter(SEQA > SEQB)  

dist.SSU <- dist.gene(SSU_ALN, method = "percentage", pairwise.deletion = FALSE,
          variance = FALSE) %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "SEQA") %>% 
  pivot_longer(-SEQA, names_to = "SEQB", values_to = "distance") %>%
  filter(SEQA > SEQB)  

dist.COI <- dist.gene(COI_ALN, method = "percentage", pairwise.deletion = FALSE,
                      variance = FALSE) %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "SEQA") %>% 
  pivot_longer(-SEQA, names_to = "SEQB", values_to = "distance") %>%
  filter(SEQA > SEQB)  

# Now we produce discrete categories the distance values 
Dist.SSU.disc <- dist.SSU %>% 
  separate(col=SEQA, into = c("SeqID_A", "Clade_A", "Genus_A", "species_A"), sep = "\\|") %>% 
  separate(col=SEQB, into = c("SeqID_B", "Clade_B", "Genus_B", "species_B"), sep = "\\|") %>% 
  select(-SeqID_A, -SeqID_B) %>% 
  mutate("Intra Species" = if_else(Clade_A == Clade_B & Genus_A == Genus_B & species_A == species_B, "YES", "NO")) %>% 
  mutate("Inter Species" = if_else(Clade_A == Clade_B & Genus_A == Genus_B & species_A != species_B, "YES", "NO")) %>% 
  mutate("Inter Genus" = if_else(Clade_A == Clade_B & Genus_A != Genus_B, "YES", "NO")) %>% 
  mutate("Inter Clade" = if_else(Clade_A != Clade_B, "YES", "NO")) %>% 
  gather(key = Category, value=dist, 8:11) %>% 
  filter(dist == "YES") %>% 
  select(-dist) %>% 
  add_column(Marker = "SSU")

Dist.COI.disc <- dist.COI %>% 
  separate(col=SEQA, into = c("SeqID_A", "Clade_A", "Genus_A", "species_A"), sep = "\\|") %>% 
  separate(col=SEQB, into = c("SeqID_B", "Clade_B", "Genus_B", "species_B"), sep = "\\|") %>% 
  select(-SeqID_A, -SeqID_B) %>% 
  mutate("Intra Species" = if_else(Clade_A == Clade_B & Genus_A == Genus_B & species_A == species_B, "YES", "NO")) %>% 
  mutate("Inter Species" = if_else(Clade_A == Clade_B & Genus_A == Genus_B & species_A != species_B, "YES", "NO")) %>% 
  mutate("Inter Genus" = if_else(Clade_A == Clade_B & Genus_A != Genus_B, "YES", "NO")) %>% 
  mutate("Inter Clade" = if_else(Clade_A != Clade_B, "YES", "NO")) %>% 
  gather(key = Category, value=dist, 8:11) %>% 
  filter(dist == "YES") %>% 
  select(-dist) %>% 
  add_column(Marker = "COI")

#Merge DF before plotting the data
All_dist_formated <- rbind(Dist.SSU.disc, Dist.COI.disc)

# Produce sub selection by clade
Dist_Spinose         <- All_dist_formated %>% filter(Clade_A == "Spinose" | Clade_B == "Spinose")                 %>% add_column(Clade = "Spinose")
Dist_NonSpinose      <- All_dist_formated %>% filter(Clade_A == "Non-spinose" | Clade_B == "Non-spinose")         %>% add_column(Clade = "Non-spinose")
Dist_Microperforates <- All_dist_formated %>% filter(Clade_A == "Microperforates" | Clade_B == "Microperforates") %>% add_column(Clade = "Microperforates")

All_dist_formated <- rbind(Dist_Spinose, Dist_NonSpinose, Dist_Microperforates)

# And plot
p1 <- ggplot(All_dist_formated, aes(x=Category, y=distance, fill = Marker)) + 
  geom_boxplot() +
  scale_y_sqrt() +
  facet_grid(.~Clade)+
  scale_fill_brewer(palette = "Pastel1")+
  stat_compare_means(aes(group = Marker), label = "p.signif", label.y = 0.7)+ 
  labs(x = "Taxonomic level")+
  labs(y = "Percentage differing sites")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        aspect.ratio = 1,
        axis.text=element_text(size=10),
        axis.title=element_text(size=15),
        strip.text = element_text(size=15))

#Export the plot
ggsave("Figure_1.pdf", plot = p1)

# Select sequences for phylogenetic inference. We infer a topology with the same species in each tree.
# List the taxa present in each dataset
Species_COI <- COI_Barcoding %>% select(Morphogroup, Genus, Species) %>% distinct()
Species_SSU <- SSU_Barcoding %>% select(Clade_curated,Genus_curated, species_curated) %>% distinct() %>% 
  rename("Morphogroup" = "Clade_curated", "Genus" = "Genus_curated", "Species" = "species_curated")
Species_common <- inner_join(Species_COI, Species_SSU) %>% unite("Name", Morphogroup, Genus, Species, sep = "_") 

# Subset sequences for the inference
# Select the longest sequence as representative
List_sequence_phylogeny_COI <- COI_Barcoding %>%  mutate(Length_COI = str_length(Sequence_COI)) %>% 
  unite("Name", Morphogroup, Genus, Species, sep = "_") %>%
  arrange(desc(Length_COI)) %>% 
  arrange(desc(Name)) %>%
  select(Name, Sequence_COI) %>% 
  group_by(Name) %>% 
  slice(1:1) %>%
  ungroup()%>% 
  inner_join(Species_common) %>% 
  rename("ID" = "Name", "Sequence" = "Sequence_COI")

List_sequence_phylogeny_SSU <- SSU_Barcoding %>% 
  unite("Name", Clade_curated, Genus_curated, species_curated, sep = "_") %>%
  filter(Quality == "BASETYPE") %>%
  arrange(desc(Lenght)) %>% 
  arrange(desc(Name)) %>%
  select(Name, FullSequence) %>% 
  group_by(Name) %>% 
  slice(1:1) %>%
  ungroup() %>% 
  inner_join(Species_common) %>% 
  rename("ID" = "Name", "Sequence" = "FullSequence")

# Now write FASTA files for the phylogenetic inference.
writeFasta(List_sequence_phylogeny_COI, "COI_phylogeny_selection.fasta")  
writeFasta(List_sequence_phylogeny_SSU, "SSU_phylogeny_selection.fasta")

# Use the ALN generated, align them at MAFFT and infer a topology using PhyML (http://www.atgc-montpellier.fr/phyml/). Name the tree files as below for the code to function.

# To have a direct comparison of the rates of evolution, let's get the patristic distances from the tree
tree_SSU <- ape::read.tree("ssu_aln_phylo_phy_phyml_tree.txt")
tree_COI <- ape::read.tree("coi_aln_phylo_phy_phyml_tree.txt")

# Extract Patristic
Patristic_SSU <- cophenetic(tree_SSU) %>% as.data.frame() %>% 
  rownames_to_column(var = "SEQA") %>% 
  pivot_longer(-SEQA, names_to = "SEQB", values_to = "distance_SSU") %>%
  filter(SEQA > SEQB)  

Patristic_COI <- cophenetic(tree_COI) %>% as.data.frame() %>% 
  rownames_to_column(var = "SEQA") %>% 
  pivot_longer(-SEQA, names_to = "SEQB", values_to = "distance_COI") %>%
  filter(SEQA > SEQB)  

#Merge the tables
All_patristic <- inner_join(Patristic_SSU, Patristic_COI) %>% 
  separate(col=SEQA, into = c("Clade_A", "Genus_A", "species_A"), sep = "_") %>% 
  separate(col=SEQB, into = c("Clade_B", "Genus_B", "species_B"), sep = "_") %>% 
  mutate("Intra Species" = if_else(Clade_A == Clade_B & Genus_A == Genus_B & species_A == species_B, "YES", "NO")) %>% 
  mutate("Inter Species" = if_else(Clade_A == Clade_B & Genus_A == Genus_B & species_A != species_B, "YES", "NO")) %>% 
  mutate("Inter Genus" = if_else(Clade_A == Clade_B & Genus_A != Genus_B, "YES", "NO")) %>% 
  mutate("Inter Clade" = if_else(Clade_A != Clade_B, "YES", "NO")) %>% 
  gather(key = Category, value=dist, 9:12) %>% 
  filter(dist == "YES") %>% 
  select(-dist) %>% 
  mutate(Ratio_SSU_COI = distance_SSU/distance_COI )


# And plot
p2 <- ggplot(All_patristic %>% filter(distance_COI > 0.0001), aes(y=distance_SSU, x=distance_COI, fill=Category)) + 
  geom_point(aes(fill=Category), colour="black",pch=21, size=7) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +  
  annotation_logticks() + 
  coord_cartesian(xlim =c(0.0001 , 4), ylim = c(0.0001 , 4)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed",size = 1, color = "grey") + 
  coord_fixed(xlim =c(0.0001 , 4), ylim = c(0.0001 , 4))+
  scale_fill_brewer(palette = "BuPu") + 
  theme(legend.position = c(0.8, 0.2))+
  labs(x = "Patristic distances COI", size=6)+
  labs(y = "Patristic distances SSU", size=6)

ggsave("Figure_2C.pdf", plot = p2)


# Reformate before plotting
Patristic_sorted <- rbind(
  (All_patristic %>% filter(Clade_A == "Microperforates" |Clade_B == "Microperforates") %>% select(distance_SSU, distance_COI, Ratio_SSU_COI, Category) %>% gather(key = Distance, value=patristic, 1:2) %>% add_column(Clade = "Microperforates")),
  (All_patristic %>% filter(Clade_A == "Non-spinose" |Clade_B == "Non-spinose")         %>% select(distance_SSU, distance_COI, Ratio_SSU_COI, Category) %>% gather(key = Distance, value=patristic, 1:2) %>% add_column(Clade = "Non-spinose")),
  (All_patristic %>% filter(Clade_A == "Spinose" |Clade_B == "Spinose")                 %>% select(distance_SSU, distance_COI, Ratio_SSU_COI, Category) %>% gather(key = Distance, value=patristic, 1:2) %>% add_column(Clade = "Spinose"))
) %>% filter(Distance == "distance_COI")


# Plot the pairwise distance and provide tests between distribution
my_comparisons <- list( c("Microperforates", "Non-spinose"), c("Microperforates", "Spinose"), c("Non-spinose", "Spinose") )

ggplot(Patristic_sorted %>% filter(patristic > 0.0001), 
       aes(y=Ratio_SSU_COI, x=Clade, fill=Clade)) + 
  geom_boxjitter(outlier.color = NA, 
                 jitter.shape = 21, 
                 jitter.color = "black", 
                 jitter.size = 3,
#                jitter.height = 0.05, 
#                jitter.width = 0.050, 
                 errorbar.draw = TRUE) +
  scale_y_log10(breaks = c(10,50,100,500,1000,5000,10000),
                labels = c(10,50,100,500,1000,5000,10000)) +
  theme_bw() +  
  coord_fixed()+
  scale_fill_brewer(palette = "BuPu") + 
  labs(y = "Ratio SSU vs COI patristic distances", size=6)+
  labs(x = "") +
  facet_grid(.~Category)+
  stat_compare_means(comparisons = my_comparisons, label.y = c(3.4, 3.8, 3.6))

# Analyses qPCR data
# add column with absolute number of gene copies for both markers
qPCR_data <- qPCR_data %>% 
  mutate(Density_COI  =  COI.gene.copy.number/`Volume.foraminifera.(µm³)`) %>% 
  rename("Volume_foram" = `Volume.foraminifera.(µm³)`)%>% filter(Volume_foram != "0")



# Plot size vs COI gene copy
p3A <- ggplot(qPCR_data, aes(x = COI.gene.copy.number, y = Volume_foram, colour = Species, group =1)) +
  geom_point(aes(fill=Species), colour="black",pch=21, size=7)+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) + 
  stat_poly_eq(use_label(c("eq","R2","P")),label.y = 0.97, label.x = 0.05,size=5)+
  geom_smooth(method = "lm", colour= "black") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed",size = 1, colour="grey") + 
  theme_bw()+
  annotation_logticks(sides="trbl")+
  labs(x = "N COI copies per individual", size=6)+
  labs(y = "Volume individual µm³", size=6) +
  coord_fixed(y = c(100, 100000000),x = c(100, 100000000)) +
  scale_fill_brewer(palette = "Paired")+
  theme(axis.text=element_text(size=20),
         axis.title=element_text(size=20),
        legend.position="none",
        legend.text=element_text(size=12, face = "italic"))

# Plot size vs COI gene copy
p3B <- ggplot(qPCR_data, aes(x = COI.gene.copy.number, y = SSU.gene.copy.number, colour = Species, group =1)) +
  geom_point(aes(fill=Species), colour="black",pch=21, size=7)+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) + 
  stat_poly_eq(use_label(c("eq","R2","P")),label.y = 0.97, label.x = 0.05,size=5)+
  geom_smooth(method = "lm", colour= "black") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed",size = 1, colour="grey") + 
  theme_bw()+
  annotation_logticks(sides="trbl")+
  labs(x = "N COI copies per individual", size=6)+
  labs(y = "N SSU copies per individual", size=6) +
  coord_fixed(y = c(100, 10000000),x = c(100, 10000000)) +
  scale_fill_brewer(palette = "Paired")+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20),
        legend.position="none",
        legend.text=element_text(size=12, face = "italic"))

combined_plot <- p3A + p3B + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")

# Adjust the layout so the two plots take up the majority of vertical space
# and the legend uses the remaining vertical space at the bottom.
final_plot <- combined_plot + 
  plot_layout(ncol = 2, heights = c(10, 1))  # 10 parts for plots, 1 part for legend

print(final_plot)

# Plot size vs COI gene copy per species
ggplot(qPCR_data, aes(x = COI.gene.copy.number, y = Volume_foram, colour = Species)) +
  geom_point(aes(fill=Species), colour="black",pch=21, size=7)+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) + 
  stat_poly_eq(use_label(c("eq","R2","P")),label.y = 0.05, label.x = 0.95, colour="black",size=3)+
  geom_smooth(method = "lm", colour="black") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed",size = 1, colour="grey") + 
  theme_bw()+
  annotation_logticks(sides="trbl")+
  labs(x = "N COI copies per individual")+
  labs(y = "Volume individual µm³") +
  coord_fixed(y = c(100, 100000000),x = c(100, 100000000)) +
  scale_fill_brewer(palette = "Paired")+
  facet_wrap(.~Species,ncol=4)+
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=20),
        legend.position="none",
        strip.text = element_text(face = "italic"))


# Plot size vs COI gene copy per species
# We need to create a manual color scale to match other plots, because there are fewer species in this plot.
cols <- c("Globigerinella siphonifera" = "#a6cee3", 
          "Globigerinita glutinata" = "#1f78b4", 
          "Neogloboquadrina dutertrei" = "#cab2d6", 
          "Neogloboquadrina pachyderma" = "#6a3d9a", 
          "Trilobatus sacculifer" = "#b15928")

pS3 <- ggplot(qPCR_data %>% filter(SSU.gene.copy.number != "NA"), aes(x = COI.gene.copy.number, y = SSU.gene.copy.number, colour = Species)) +
  geom_point(aes(fill=Species), colour="black",pch=21, size=7)+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) + 
  stat_poly_eq(use_label(c("eq","R2","P")),label.y = 0.95, label.x = 0.95, colour="black",size=4)+
  geom_smooth(method = "lm", colour="black") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed",size = 1, colour="grey") + 
  theme_bw()+
  annotation_logticks(sides="trbl")+
  labs(x = "N COI copies per individual")+
  labs(y = "N SSU copies per individual") +
  coord_fixed(y = c(100, 100000000),x = c(100, 100000000)) +
  scale_fill_manual(values = cols)+
  facet_wrap(.~Species,ncol=3, drop = TRUE)+
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=20),
        legend.position="none",
        strip.text = element_text(face = "italic"))

# Look at density data, add statistics on the plot
# First let s see if our data are normally distributed
Res_Norm_test_COI_density <- qPCR_data %>% 
  mutate(log_density = log10(Density_COI)) %>% 
  group_by(Species) %>% 
  do(tidy(shapiro.test(.$log_density)))

# Not all species distribution are normal, even after logging the data. Do Non parametric test
ggplot(qPCR_data, aes(x = Species,y =Density_COI, fill=Species)) + 
  scale_y_log10()+
  scale_fill_brewer(palette = "Paired")+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = "black", jitter.size = 4,
                 jitter.height = 0.05, jitter.width = 0.050, errorbar.draw = TRUE) +
  theme_bw()  +
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        legend.text = element_text(size=15, face = "italic"),
        legend.position="bottom") +
  stat_compare_means(method = "kruskal", label.y = 0.0001) +      # Add global p-value
  stat_compare_means(label = "p.signif", method = "wilcox.test",
                     ref.group = ".all.")+
  labs(y = "N COI copies per µm³", size=20) 
  

# Compare PF and LBF data
LBF_PF_data <- rbind(
  (LBF_data %>% 
  select(Species, `Calculated.volume.chambers.(mm³)`,Total.gene.copies.in.specimen,`gene.copies.per.um³`) %>% 
  rename("COI.gene.copy.number"= "Total.gene.copies.in.specimen",
         "Density_COI" = "gene.copies.per.um³",
          "calculated_volume_chambers_mm3" = "Calculated.volume.chambers.(mm³)") %>% 
    add_column(Group = "Larger Benthic Foraminifera") %>% 
  mutate(Volume_foram = calculated_volume_chambers_mm3 * 1000 * 1000* 1000) %>% 
  select(-calculated_volume_chambers_mm3)),
  qPCR_data %>% select(Species,COI.gene.copy.number,Density_COI, Volume_foram) %>% 
    add_column(Group = "Planktonic Foraminifera")
)

LBF_PF_data$COI.gene.copy.number <- as.numeric(LBF_PF_data$COI.gene.copy.number)


Fig5A <- ggplot(LBF_PF_data, aes(x = Volume_foram, y = COI.gene.copy.number, colour = Group, group =1)) +
  geom_point(aes(fill=Group), colour="black",pch=21, size=7)+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) + 
  stat_poly_eq(use_label(c("eq","R2","P")),label.y = 0.95, label.x = 0.05,size=5)+
  geom_smooth(method = "lm", colour= "black") +
  theme_bw()+
  annotation_logticks(sides="trbl")+
  labs(x = "Volume cell µm³", size=6)+
  labs(y = "N COI copies per cell", size=6) +
  coord_fixed() +
  scale_fill_brewer(palette = "Spectral")+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20),
        legend.position = "none",
        legend.text=element_text(size=12),
        legend.title = element_blank())

Fig5B <- ggplot(LBF_PF_data, aes(x = Volume_foram, y = Density_COI, colour = Group, group =1)) +
  geom_point(aes(fill=Group), colour="black",pch=21, size=7)+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) + 
  stat_poly_eq(use_label(c("eq","R2","P")),label.y = 0.95, label.x = 0.95,size=5)+
  geom_smooth(method = "lm", colour= "black") +
  theme_bw()+
  annotation_logticks(sides="trbl")+
  labs(x = "Volume cell µm³", size=6)+
  labs(y = "N COI copies per µm³", size=6) +
  coord_fixed() +
  scale_fill_brewer(palette = "Spectral")+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20),
        legend.position = c(0.25, 0.1),
        legend.text=element_text(size=12),
        legend.title = element_blank())

(Fig5A | Fig5B)

