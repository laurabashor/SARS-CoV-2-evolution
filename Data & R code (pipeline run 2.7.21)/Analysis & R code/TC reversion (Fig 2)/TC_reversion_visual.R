##visualizing TC variant emergence and reversion

#load libraries
library(tidyverse)
library(ggplot2)

#load in cleaned data with replicates combined
df <- read.csv("variant_summary_processed.csv")

#NAs back to 0s for plotting
df[is.na(df)] <- 0

##want show the reversion of these variants in vivo for all the animal species
#pull out just TC mutations that went to fixation
df2 <- df %>%
  filter(variant %in% c("D135E", "D215H", "R685H", "T7I","S194T"))

#make data long and skinny for plotting, rename value as frequency
#get a new column for the species/passage
#get the order right for passages
#change frequency to a percentage
df2 <- df2 %>%
  pivot_longer(!c(reference_sequence, position, gene, indel,variant, reference_base, variant_base, effect), 
               names_to= "dataset_ID", values_to ="frequency") %>%  
  mutate(passage = case_when(grepl("C", dataset_ID) ~ "Cats",
                           grepl("D", dataset_ID) ~"Dogs",
                           grepl("H", dataset_ID) ~"Hamsters",
                           grepl("F", dataset_ID) ~"Ferret",
                           grepl("_1", dataset_ID) ~"P1",
                           grepl("_2", dataset_ID) ~"P2",
                           grepl("_3", dataset_ID) ~"P3")) %>%
  mutate(passage = fct_relevel(passage, "P1","P2","P3", "Cats", "Dogs","Hamsters","Ferret")) %>%
  mutate(percentage = (frequency*100))

#get the labels to say both the variant and cds
labels <- c("D135E (nsp12)", "D215H (Spike)", "R685H (Spike)", "S194T (Nucleocapsid)","T7I (Membrane)")
names(labels) <- c("D135E", "D215H", "R685H", "S194T", "T7I")

#plot it
p<- ggplot() +
  geom_point(data=df2, aes(x=passage, y=percentage, group=variant), show.legend=FALSE) +
  geom_line(data=df2 %>%
              filter(passage == "P1"|passage == "P2"|passage == "P3"), aes(x=passage, y=percentage, group=variant), show.legend=FALSE)+
  geom_boxplot(data=df2  %>%
                 filter(passage=="Cats"|passage=="Dogs"|passage=="Hamsters"|passage=="Ferret"), aes(x=passage, y=percentage))+
  geom_jitter(data=df2  %>%
                filter(passage=="Cats"|passage=="Dogs"|passage=="Hamsters"|passage=="Ferret"), aes(x=passage, y=percentage, group=variant),position=position_jitter(w=0.1))+
  #scale_colour_manual(values=colors)+
  facet_wrap(vars(variant), labeller=labeller(variant=labels), scales='free_x') +
  theme_bw()+
  theme(text = element_text(size=12, family="sans"), 
        strip.text.x=element_text(size=14), 
        axis.title=element_text(size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) + 
  xlab ("Passage") + 
  ylab ("Variant Frequency (%)") 

print(p)

pdf(file="TC_variant_reversion.pdf", width=12, height=10)
p
dev.off()
