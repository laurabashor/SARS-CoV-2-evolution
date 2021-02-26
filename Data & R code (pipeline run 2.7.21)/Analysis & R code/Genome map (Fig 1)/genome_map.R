###mapping variants along genome in ggplot

#load packages
library(tidyverse)
library(ggforce)
library(RColorBrewer)
library(ggpubr)

#start with data frame of cleaned data, combined replicates
df <- read.csv("variant_summary_processed.csv")

#make data long and skinny 
df2 <- df %>%
  pivot_longer(!c(reference_sequence, position, gene, indel,variant, reference_base, variant_base, effect), names_to= "dataset_ID", values_to ="frequency")

#make nice axis labels
labels <- c("Passage 1", "Passage 2", "Passage 3", "Cat 1", "Cat 2", "Cat 3", "Cat 4", "Cat 5", "Cat 6", "Dog 1", "Dog 2","Dog 3", "Ferret 1", "Hamster 1", "Hamster 2", "Hamster 3")

#make a plot to show all variants that were above 25% frequency, shaded based on frequency
#for this plot get rid of the x axis so that we can align it to a genome schematic

df0.25 <- df2 %>%
  mutate(frequency=ifelse(frequency<0.25, NA, frequency)) %>% 
  mutate(category=cut(frequency,breaks=c(-Inf, 0.5, 0.75, Inf),labels=c("25-50%","50-75%","75-100%")))

p0.25 <- ggplot(df0.25) +
  geom_point(data=df0.25, 
             aes(y=factor(dataset_ID, levels = c("Passage_1","Passage_2", "Passage_3", 
                                                 "Cat_1", "Cat_2", "Cat_3", "Cat_4", "Cat_5", "Cat_6", 
                                                 "Dog_1", "Dog_2","Dog_3", "Ferret_1", 
                                                 "Hamster_1", "Hamster_2", "Hamster_3")), 
                 x=position), color="white") +
  geom_rect(aes(xmin=21563, xmax=25384,ymin=-Inf,ymax=Inf), fill="lavenderblush1", alpha=0.5)+
  geom_point(data=(df0.25 %>% filter(!is.na(frequency),indel=="FALSE")), 
             aes(y=dataset_ID,
                 x=position, fill=frequency), shape=21, size=4, alpha=0.8, col="transparent") +
  geom_point(data=(df0.25 %>% filter(!is.na(frequency), indel=="TRUE")), 
             aes(y=dataset_ID, x=position, fill=frequency), 
             shape=24,size=4,alpha=0.8, col="transparent") +
  scale_fill_gradient(low="skyblue2", high="#08306b", 
                      breaks=c(0.26,1),
                      labels=c("25%","100%")) +
  scale_x_continuous(limits=c(0,30000)) +
  theme_classic() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text=element_text(size=16), 
        axis.title.y=element_blank(), 
        axis.text.x = element_blank(),
        legend.text=element_text(size=9), 
        legend.title=element_text(size=12),
        legend.background=element_rect(color="black"),
        plot.margin = unit(c(0.7,0.7,0.7,0.7), "cm")) +
  scale_y_discrete(labels=labels) +
  labs (x="Position in genome (nt)", fill="Variant \nFrequency \n")

print(p0.25)


#now make variant effects plot
#first need to simplify the effects

df2 %>%
  na.omit(df2) %>%
  filter(indel==TRUE) %>%
  count(effect)
#indels can be seven things:
#conservative inframe deletion,
#conservative_inframe_insertion, 
# disruptive_inframe_deletion, 
#disruptive_inframe_insertion, 
# frameshift_variant
# frameshift_variant&stop_gained 
# frameshift_variant&stop_lost&splice_region_variant 
#we can shrink that to three things: inframe deletion, inframe insertion, frameshift

df2 %>%
  na.omit(df2) %>%
  filter(indel==FALSE) %>%
  count(effect)
#SNVs can be four things:
#missense
#stop_lost
#stop_retained
#synonymous
#we can shrink that to two things: nonsynonymous and synonymous (stop retained is synonymous)

df3 <- df2 %>%
  na.omit(df2) %>%
  mutate(effect2 = ifelse(effect %in% c("missense","stop_lost"), 
                  "nonsynonymous",
                          ifelse(effect %in% c("stop_retained_variant","synonymous"),
                                 "synonymous",
                                 ifelse(effect %in% c("conservative_inframe_deletion","disruptive_inframe_deletion"),
                                        "inframe_deletion",
                                        ifelse(effect %in% c("conservative_inframe_insertion", "disruptive_inframe_insertion"), 
                                               "inframe_insertion",
                                               ifelse(effect %in% c("frameshift","frameshift_variant&stop_gained","frameshift_variant&stop_lost&splice_region"),
                                                      "frameshift", NA))))))

df3 %>%
  count(effect2)

#now we can plot it

colors <- c("forestgreen","#FF7F00","gold1","skyblue2","maroon")

p<- ggplot(df3) +
  geom_point(data=(df3 %>% filter(!is.na(frequency))), 
             aes(y=factor(dataset_ID, 
                          levels = c("Passage_1","Passage_2", "Passage_3", 
                                     "Cat_1", "Cat_2", "Cat_3", "Cat_4", "Cat_5", "Cat_6",
                                     "Dog_1", "Dog_2","Dog_3", "Ferret_1",
                                     "Hamster_1", "Hamster_2", "Hamster_3")), 
                 x=position, col=factor(effect2), shape=indel), 
             size=2.5, alpha=0.7, stroke=FALSE)  +
  scale_color_manual(name="Variant Effect", 
                     labels=c("Frameshift", 
                              "Inframe deletion",
                              "Inframe insertion",
                              "Nonsynonymous SNV",
                              "Synonymous SNV"),
                     values=colors)+
  scale_shape(name="Variant Type", labels=c("Single nucleotide variant", "Structural variant"))+
  theme_bw()+
  labs(x="Position in genome (nt)") + 
  scale_y_discrete(labels=labels)+
  theme_classic()+
  theme(axis.title.x = element_text(size=14, family="sans"),
        legend.background=element_rect(color="black"),
        axis.title.y=element_blank(),
        axis.text=element_text(size=12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_text(size = 10), 
        legend.text  = element_text(size = 8),
        legend.key.size = unit(0.8, "lines"),
        plot.margin = unit(c(0.7,0.7,0.7,0.7), "cm")) +
  guides(shape = guide_legend(override.aes = list(size = 2)),
         color = guide_legend(override.aes = list(size = 2))) 

print(p)

p_combined <- ggarrange(p, p0.25, ncol=1, nrow=2, align="v", labels=c("A","B"))

print(p_combined)

pdf(file="genome_map.pdf", width=12, height=9)
p_combined
dev.off()
