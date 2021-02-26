#load libraries
library(tidyverse)
library(ggrepel)
library(ggpubr)
library(RColorBrewer)

#read data
df <- read.csv("variant_summary_processed.csv")
length <- read.csv("gene_lengths.csv")

df2 <- df %>%
  pivot_longer(!c(reference_sequence, position, gene, indel,variant, reference_base, variant_base, effect), names_to= "dataset_ID", values_to ="frequency")%>%  
  mutate(species = case_when(grepl("C", dataset_ID) ~ "Cats",
                            grepl("D", dataset_ID) ~"Dogs",
                            grepl("H", dataset_ID) ~"Hamsters",
                            grepl("F", dataset_ID) ~"Ferret",
                            grepl("P", dataset_ID) ~"Vero")) %>%
  na.omit(df2)

#(1)
#richness (#variants) by species

#count the number of distinct variants to get richness
richness <- df2 %>% 
  group_by(dataset_ID) %>% 
  summarise(n=n_distinct(variant)) %>%
  mutate(species = case_when(grepl("C", dataset_ID) ~ "Cats",
                             grepl("D", dataset_ID) ~"Dogs",
                             grepl("H", dataset_ID) ~"Hamsters",
                             grepl("F", dataset_ID) ~"Ferret",
                             grepl("P", dataset_ID) ~"Vero"))

p1 <- ggplot(richness, aes(x=species, y=n)) +
  geom_dotplot(aes(x=species, y=n, fill=species), show.legend=FALSE, 
               binaxis="y", stackdir="center", stroke=FALSE) + 
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
                 geom="errorbar", width=0.2) +
  stat_summary(fun=mean, geom="crossbar", width=0.4) +
  #stat_compare_means(method = "anova", label.x.npc = "middle", label.y.npc="top")+
  scale_fill_brewer(palette="Paired")+
  labs(x="Species",y="Number of Unique Variants") +
  theme_classic()+
  theme(axis.text=element_text(size=12), axis.title.x=element_blank(), 
        axis.title.y=element_text(size=14))

p1
# pdf("p1.pdf")
# p1
# dev.off()

compare_means(n~species, data=richness) #no significant difference between any of the animals doing t-tests
aov <- aov(n~species, data=richness)
summary(aov) #p=0.224


#now look at difference between two cat cohorts
cohorts <- richness %>%
  filter(species=="Cats") %>%
  mutate(cohort = ifelse(dataset_ID %in% c("Cat_1", "Cat_5","Cat_6"),2,1))

#paired, tailed (alternative: expect one to be bigger than the other)
t.test(cohorts$cohort, cohorts$n, paired=TRUE, alternative="less") #p=0.01162



#(2)
#now plot the proportion of variants that were detected in however many individuals of cats, dogs, hamsters

df.cats <- df %>% 
  select(reference_sequence,position, gene, indel,variant, reference_base, variant_base, effect, starts_with("C"))

df.dogs <- df %>% 
  select(reference_sequence,position, gene, indel,variant, reference_base, variant_base, effect, starts_with("D"))

df.hamsters <- df %>% 
  select(reference_sequence,position, gene, indel,variant, reference_base, variant_base, effect, starts_with("H"))

#use a function to delete any rows in which there are NAs for all cats
delete.na <- function(df.cats, n=0) {
  df.cats[rowSums(is.na(df.cats)) <= n,]
}
df.cats <- delete.na(df.cats,5) #5 is max NAs allowed for 6 cats

# #do it for dogs
delete.na.dogs <- function(df.dogs, n=0) {
  df.dogs[rowSums(is.na(df.dogs)) <= n,]
}
df.dogs <- delete.na.dogs(df.dogs, 2) #2 is max NAs allowed for 3 dogs

#do it for hamsters
delete.na.ham <- function(df.hamsters, n=0) {
  df.hamsters[rowSums(is.na(df.hamsters)) <= n,]
}
df.hamsters <- delete.na.ham(df.hamsters, 2) #2 is max NAs allowed for 3 hamsters

#now count variants and proportion of variants by #individuals of each species
cat.count <- df.cats %>% 
  mutate(number_of_cats=(6-rowSums(is.na(.)))) %>%
  count(number_of_cats) %>%
  mutate(percentage=(n/sum(n)*100))%>%
  mutate(number_of_variants=n) %>%
  select(-n)%>%
  mutate(species="cat") %>%
  mutate(number_of_individuals=number_of_cats) %>%
  select(-number_of_cats)

#dogs
dog.count <- df.dogs%>%
  mutate(number_of_dogs=(3-rowSums(is.na(.))))%>%
  count(number_of_dogs) %>%
  mutate(percentage=(n/sum(n)*100))%>%
  mutate(number_of_variants=n) %>%
  select(-n) %>%
  mutate(species="dog") %>%
  mutate(number_of_individuals=number_of_dogs)%>%
  select(-number_of_dogs)

#hamsters
ham.count <- df.hamsters%>%
  mutate(number_of_hams=(3-rowSums(is.na(.))))%>%
  count(number_of_hams) %>%
  mutate(percentage=(n/sum(n)*100)) %>%
  mutate(number_of_variants=n) %>%
  select(-n)%>% 
  mutate(species="hamster") %>%
  mutate(number_of_individuals=number_of_hams) %>%
  select(-number_of_hams)

#stick them together
all <- full_join(ham.count, cat.count)
all <- full_join(all,dog.count)

###do this all again, but with only variants that weren't detected in the inoculum
#first filter data frame so that you only keep rows where there was NA for all three passages of the inoculum
df.filtered <-
  df %>%
  filter_at(vars(starts_with("Passage")), all_vars(is.na(.)))

#start with breaking down variants by number of individual animals for each species
df.cats2 <- df.filtered %>% 
  select(reference_sequence,position, gene, indel,variant, reference_base, variant_base, effect, starts_with("C"))

df.dogs2 <- df.filtered %>% 
  select(reference_sequence,position, gene, indel,variant, reference_base, variant_base, effect, starts_with("D"))

df.hamsters2 <- df.filtered %>% 
  select(reference_sequence,position, gene, indel,variant, reference_base, variant_base, effect, starts_with("H"))

#use a function to delete any rows in which there are NAs for all cats
delete.na2 <- function(df.cats2, n=0) {
  df.cats2[rowSums(is.na(df.cats2)) <= n,]
}
df.cats2 <- delete.na2(df.cats2,5) #5 is max NAs allowed for 6 cats

# #do it for dogs
delete.na.dogs2 <- function(df.dogs2, n=0) {
  df.dogs2[rowSums(is.na(df.dogs2)) <= n,]
}
df.dogs2 <- delete.na.dogs2(df.dogs2, 2) #2 is max NAs allowed for 3 dogs

#do it for hamsters
delete.na.ham2 <- function(df.hamsters2, n=0) {
  df.hamsters2[rowSums(is.na(df.hamsters2)) <= n,]
}
df.hamsters2 <- delete.na.ham2(df.hamsters2, 2) #2 is max NAs allowed for 3 hamsters

#now count variants and proportion of variants by #individuals of each species
cat.count2 <- df.cats2 %>% 
  mutate(number_of_cats=(6-rowSums(is.na(.)))) %>%
  count(number_of_cats) %>%
  mutate(percentage=(n/sum(n)*100))%>%
  mutate(number_of_variants=n) %>%
  select(-n)%>%
  mutate(species="cat") %>%
  mutate(number_of_individuals=number_of_cats) %>%
  select(-number_of_cats)

#dogs
dog.count2 <- df.dogs2%>%
  mutate(number_of_dogs=(3-rowSums(is.na(.))))%>%
  count(number_of_dogs) %>%
  mutate(percentage=(n/sum(n)*100))%>%
  mutate(number_of_variants=n) %>%
  select(-n)%>%
  mutate(species="dog") %>%
  mutate(number_of_individuals=number_of_dogs)%>%
  select(-number_of_dogs)

#hamsters
ham.count2 <- df.hamsters2%>%
  mutate(number_of_hams=(3-rowSums(is.na(.))))%>%
  count(number_of_hams) %>%
  mutate(percentage=(n/sum(n)*100)) %>%
  mutate(number_of_variants=n) %>%
  select(-n)%>% 
  mutate(species="hamster") %>%
  mutate(number_of_individuals=number_of_hams) %>%
  select(-number_of_hams)

all2 <- full_join(ham.count2, cat.count2)
all2 <- full_join(all2,dog.count2)

pA <- ggplot(all, aes(x=species, y=percentage, fill=factor(number_of_individuals))) + 
  geom_bar(stat="identity")+
  scale_fill_brewer(palette="Paired")+
  theme_classic()+
  scale_x_discrete(labels=c("Cat", "Dog", "Hamster"))+
  coord_flip()+
  labs(y="Percentage of Variants", fill="Number of\nindividuals")+
  theme(text=element_text(size=14),
        axis.title=element_blank())
pA

pB <- ggplot(all2, aes(x=species, y=percentage, fill=factor(number_of_individuals))) + 
  geom_bar(stat="identity")+
  scale_fill_brewer(palette="Paired")+
  theme_classic()+
  scale_x_discrete(labels=c("Cat", "Dog", "Hamster"))+
  labs(y="Percentage of Variants", fill="Number of\nindividuals")+
  coord_flip()+
  theme(axis.title.x=element_text(size=14),axis.text = element_text(size=12),
        axis.title.y=element_blank(), legend.position="none")
pB

p2 <- ggarrange(pA, pB, ncol=1, nrow=2, labels=c("(a)", "(b)"),
                label.x=0.06,
                common.legend=TRUE, legend="right")

p2

# pdf("p2.pdf")
# p2
# dev.off()

#now put in distribution plot 

#count every single observation of all variants by CDS
df.long <- df %>%
  pivot_longer(!c(reference_sequence, position, gene, indel,variant, reference_base, variant_base, effect), names_to= "dataset_ID", values_to ="frequency")

count_all <- df.long %>%
  count(gene)

#add length of cds to data frame
count_all <- merge(count_all, length)

count_all <- count_all %>%
  mutate(proportion_of_genome=length/29903) %>%
  mutate(proportion_of_variants=(n/sum(n))) 

#now plot it
#geom smooth applies a linear model, shows line of best fit & confidence interval (0.95 by default)

p3 <- ggplot(count_all, aes(x=proportion_of_genome, y=proportion_of_variants, label=gene))+
  geom_point()+
  geom_text_repel(aes(label=gene))+
  geom_smooth(method="lm")+
  labs(x="Gene length/genome length", y="Proportion of variants")+
  theme_classic()+
  theme(legend.position = "none",axis.text = element_text(size=12), 
        axis.title.y = element_text(size=14), 
        axis.title.x = element_text(size=14))

p3

lm <- lm(proportion_of_variants~proportion_of_genome, count_all)
summary(lm) #R-squared = 0.6876, p=4.348e-05

#finally, look at whether variants were found in only one species vs in 2,3, or all 4
#only one variant that wasn't a TC variant was found in all 4 species
#do this for animal variants not found in the inoculum

df3 <- df.filtered %>%
  pivot_longer(!c(reference_sequence, position, gene, indel,variant, reference_base, variant_base, effect), names_to= "dataset_ID", values_to ="frequency")%>%
  mutate(Animal = case_when(grepl("C", dataset_ID) ~ "Cats",
                            grepl("D", dataset_ID) ~"Dogs",
                            grepl("H", dataset_ID) ~"Hamsters",
                            grepl("F", dataset_ID) ~"Ferret")) %>%
  filter(dataset_ID != "Passage_1") %>%
  filter(dataset_ID != "Passage_2") %>%
  filter(dataset_ID != "Passage_3") %>%
  na.omit(df3)

number_of_species <- df3 %>%
  group_by(variant)%>%
  count(Animal) %>%
  count(variant, name="number_of_species") %>% #now we have the count for how many variants are found in 1, 2, 3, or 4 animal species, we need to summarize it
  group_by(number_of_species) %>%
  tally() %>%
  mutate(percentage=(n/sum(n)*100))%>%
  mutate(type="not_in_inoc")
# 
# #NOTE: this plot is only for variants not detected in the inoculum
# p4 <-ggplot(number_of_species, aes(x=number_of_species, y=n, fill=factor(number_of_species))) +
#   geom_bar(stat="identity")+
#   theme_classic()+
#   labs(x= "Number of Species", y="Number of Variants")+
#   theme(axis.title=element_text(size=16),
#          axis.text=element_text(size=14), legend.position="none")

# pdf("p4.pdf")
# p4
# dev.off()


#finally put them all together

pdf("Variant_summary_plot.pdf", width=11, height=8)
ggarrange(p1,p2,p3, ncol=2, nrow=2, labels="AUTO")
dev.off()


# library(grid)
# # Move to a new page
# grid.newpage()
# # Create layout : nrow = 3, ncol = 2
# pushViewport(viewport(layout = grid.layout(nrow = 2, ncol = 2)))
# # A helper function to define a region on the layout
# define_region <- function(row, col){
#   viewport(layout.pos.row = row, layout.pos.col = col)
# } 
# 
# # Arrange the plots
# print(p1, vp = define_region(row = 1, col = 1))
# print(p2, vp= define_region(row = 1, col = 2))
# print(p3, vp = define_region(row = 2, col = 1:2))
# 
# library("cowplot")
# ggdraw() +
#   draw_plot(p1, x = 0, y = .5, width = .5, height = .5) +
#   draw_plot(p2, x = .5, y = .5, width = .5, height = .5) +
#   draw_plot(p3, x = 0, y = 0, width = 1, height = 0.5) +
#   draw_plot_label(label = c("A", "B", "C"), size = 15,
#                   x = c(0, 0.5, 0), y = c(1, 1, 0.5))


#wondering if number of variants would be significant if we just looked at higher frequency variants
# 
# richness0.5 <- df2 %>%
#   filter(frequency>0.5)%>%
#   group_by(dataset_ID) %>% 
#   summarise(n=n_distinct(variant)) %>%
#   mutate(species = case_when(grepl("C", dataset_ID) ~ "Cats",
#                              grepl("D", dataset_ID) ~"Dogs",
#                              grepl("H", dataset_ID) ~"Hamsters",
#                              grepl("F", dataset_ID) ~"Ferret",
#                              grepl("P", dataset_ID) ~"Vero"))
# 
# ggplot(richness0.5, aes(x=species, y=n)) +
#   geom_dotplot(aes(x=species, y=n, fill=species), show.legend=FALSE, 
#                binaxis="y", stackdir="center", stroke=FALSE) + 
#   stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
#                geom="errorbar", width=0.2) +
#   stat_summary(fun=mean, geom="crossbar", width=0.4) +
#   #stat_compare_means(method = "anova", label.x.npc = "middle", label.y.npc="top")+
#   scale_fill_brewer(palette="Paired")+
#   labs(x="Species",y="Number of Unique Variants") +
#   theme_classic()+
#   theme(axis.text=element_text(size=14), axis.title.x=element_blank(), 
#         axis.title.y=element_text(size=16))
# 
# compare_means(n~species, data=richness0.5) #no significant difference between any of the animals doing t-tests
# aov <- aov(n~species, data=richness0.5)
# summary(aov)
# TukeyHSD(aov)#looks like when you look at variants with >50% frequency, 
# #dogs do have significantly more variants than cats or hamsters

#wondering if it would be significant on the basis of structural vs. nonstructural vs. accessory

# richness_type <- df2 %>%
#   filter(frequency>0.5)%>%
#   group_by(gene) %>% 
#   summarise(n=n_distinct(variant)) %>%
#   mutate(gene_type = case_when(grepl("nsp", gene) ~ "nonstructural",
#                                grepl("ORF", gene) ~"accessory",
#                                grepl("S", gene) ~"structural",
#                                grepl("M", gene) ~"structural",
#                                grepl("N", gene) ~"structural",
#                                grepl("E", gene) ~"structural"))
# 
# aov_type <- aov(n~gene_type, data=richness_type)
# summary(aov_type)
# TukeyHSD(aov_type)
#no sig difference between number of variants based on type of protein whether you just look at >50% or not

