##visualizing TC variant emergence and reversion

#load libraries
library(tidyverse)
library(data.table)
library(ggplot2)

#load in cleaned data with replicates combined
df <- as.data.frame(variant_table_cleaned_combined)

#NAs to 0s for plotting
df[is.na(df)] <- 0

###first graph just the inoculums
df1 <- df[, c("variant", "cds", "SARS2_P1","SARS2_P2","SARS2_P3")]

#rename passages as 1, 2, 3
names(df1) <- gsub("SARS2_P", "", names(df1))

#pull out just TC mutations that went to fixation
df1 <- df1 %>%
  filter(variant %in% c("D135E", "D215H", "R685H", "T7I","S194T"))

#make data long and skinny for plotting
df1.long <- df1 %>%
  pivot_longer(!c(variant, cds), names_to= "passage", values_to ="frequency")

#plot
pdf(file="TC_variants.pdf")
ggplot(df1.long, aes(x=passage, y=frequency, group=variant, color= variant)) +
  geom_point(show.legend=FALSE) +
  geom_line(show.legend=FALSE) +
  facet_wrap(vars(variant))+
  theme_bw()+
  theme(text = element_text(size=10, family="sans")) + 
  xlab ("Passage") + 
  ylab ("Variant Frequency") 
dev.off()

##want show the reversion of these variants in vivo for all the animal species

df3 <- df %>%
  select(!c("position", "N_or_S"))

#pull out just TC mutations that went to fixation
df3 <- df3 %>%
  filter(variant %in% c("D135E", "D215H", "R685H", "T7I","S194T"))

#make data long and skinny for plotting, rename value as frequency
df3.long <- df3 %>%
  pivot_longer(!c(variant,cds), names_to= "dataset", values_to ="Frequency")

df3.long <- df3.long %>%  
  mutate(Passage = case_when(grepl("C", dataset) ~ "Cats",
                             grepl("D", dataset) ~"Dogs",
                             grepl("H", dataset) ~"Hamsters",
                             grepl("F", dataset) ~"Ferret",
                             grepl("P1", dataset) ~"P1",
                             grepl("P2", dataset) ~"P2",
                             grepl("P3", dataset) ~"P3"))

#get the order right for passages
df3.long <- df3.long %>%
  mutate(Passage = fct_relevel(Passage, "P1","P2","P3", "Cats", "Dogs","Hamsters","Ferret"))
#df3.long$Passage <- factor(df3.long$Passage, levels = c("P1","P2","P3", "Cats", "Dogs","Hamsters","Ferret"))

#colors <- c("P1"="black","P2"="black","P3"="black", "Cats"="red", "Dogs"="blue","Hamsters"="green","Ferret"="purple")
# animals <- df3.long %>%
#   filter(Passage=="Cats"|Passage=="Dogs"|Passage=="Hamsters"|Passage=="Ferret")

#plot
pdf(file="TC_variants_reversion.pdf", width=10, height=10)
ggplot() +
  geom_point(data=df3.long, aes(x=Passage, y=Frequency, group=variant), show.legend=FALSE) +
  geom_line(data=df3.long %>%
              filter(Passage == "P1"|Passage == "P2"|Passage == "P3"), aes(x=Passage, y=Frequency, group=variant), show.legend=FALSE)+
  geom_boxplot(data=df3.long  %>%
                 filter(Passage=="Cats"|Passage=="Dogs"|Passage=="Hamsters"|Passage=="Ferret"), aes(x=Passage, y=Frequency))+
  geom_jitter(data=df3.long  %>%
                filter(Passage=="Cats"|Passage=="Dogs"|Passage=="Hamsters"|Passage=="Ferret"), aes(x=Passage, y=Frequency, group=variant),position=position_jitter(w=0.1))+
  #scale_colour_manual(values=colors)+
  facet_wrap(vars(variant, cds), scales='free_x') +
  theme_bw()+
  theme(text = element_text(size=12, family="sans")) + 
  xlab ("Passage") + 
  ylab ("Variant Frequency") 
dev.off()
