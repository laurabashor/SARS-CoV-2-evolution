#selection analysis using output from SNPGenie output 
#population_summary file and by product_results file
#files were generated individually for each dataset then combined in command line

#packages
library(tidyverse)
library(ggpubr)
library(openxlsx)

#read in data
df_pop <- read.csv("pop_summary.csv")
df_prod <- read.csv("prod_results.csv")

#clean up data
#rename file as dataset_ID
df_pop <- df_pop %>%
  mutate(dataset_ID=file)%>%
  select(-file)

df_prod <- df_prod %>%
  mutate(dataset_ID=file)%>%
  select(-file)

#get rid of file name strings
df_pop$dataset_ID <- gsub("vcfs_1.25.21/","",df_pop$dataset_ID)
df_pop$dataset_ID <- gsub("_R.vcf","",df_pop$dataset_ID)
df_pop$dataset_ID <- gsub(".vcf","",df_pop$dataset_ID)

df_prod$dataset_ID <- gsub("vcfs_1.25.21/","",df_prod$dataset_ID)
df_prod$dataset_ID <- gsub("_R.vcf","",df_prod$dataset_ID)
df_prod$dataset_ID <- gsub(".vcf","",df_prod$dataset_ID)

#make a species category, and get rid of the vero cell passaged samples for now
df_pop <- df_pop%>%
  mutate(species = case_when(grepl("C", dataset_ID) ~ "Cats",
                          grepl("D", dataset_ID) ~"Dogs",
                          grepl("H", dataset_ID) ~"Hamsters",
                          grepl("F", dataset_ID) ~"Ferret")) %>%
  filter(dataset_ID != "Passage_1") %>%
  filter(dataset_ID != "Passage_2") %>%
  filter(dataset_ID != "Passage_3") %>%
  mutate(species=fct_relevel(species,"Dogs", "Cats","Hamsters","Ferret"))

df_prod <- df_prod%>%
  mutate(species = case_when(grepl("C", dataset_ID) ~ "Cats",
                             grepl("D", dataset_ID) ~"Dogs",
                             grepl("H", dataset_ID) ~"Hamsters",
                             grepl("F", dataset_ID) ~"Ferret")) %>%
  filter(dataset_ID != "Passage_1") %>%
  filter(dataset_ID != "Passage_2") %>%
  filter(dataset_ID != "Passage_3")

#calculate piN/piS and add some things to compare piN and piS
df_pop <- df_pop %>%
  mutate(piN_div_piS=(piN/piS)) %>% #allows you to see if its >1, =1, or <1
  mutate(piN_greater=ifelse((piN>piS), TRUE, FALSE)) %>% #allows you to see if piN>piS, even if one of them is 0
  mutate(piS_greater=ifelse((piS>piN), TRUE, FALSE)) %>%
  select(dataset_ID, pi, species, piN, piS,piN_div_piS)

df_prod <- df_prod %>%
  mutate(piN_div_piS=(piN/piS)) %>% #allows you to see if its >1, =1, or <1
  mutate(piN_greater=ifelse((piN>piS), TRUE, FALSE)) %>% #allows you to see if piN>piS, even if one of them is 0
  mutate(piS_greater=ifelse((piS>piN), TRUE, FALSE)) %>%
  select(dataset_ID, gene_product, species, piN, piS,piN_div_piS)%>%
  arrange(gene_product)

#export these as an excel sheet for supplemental tables

table1 <- df_pop %>%
  select(-species)

table2 <- df_prod %>%
  select(-species)

wb <- createWorkbook("supplemental_selection_tables.xlsx")
addWorksheet(wb, "population_level")
addWorksheet(wb, "gene_product_level")
writeData(wb, "population_level", table1, borders="all")
writeData(wb, "gene_product_level", table2,borders="all")
saveWorkbook(wb, "supplemental_selection_tables.xlsx", overwrite = TRUE)

##now do t-tests to show that nonsynonymous is not equal to synonymous

#first, look at the difference between piN and piS at a population level
t.test(x=df_pop$piN, y=df_pop$piS,paired=TRUE)

#if you look at all the data together, significant difference between piN and piS
#paired t-test, piN overall was significantly greater p=0.001599
#if you include vero cells, p=0.0006363

#break it down by species 
#have to get rid of the ferret because there aren't multiple values of piN and piS
#this summarizes the data by species, then does a t-test between piN and piS for each of three species with multiple individuals
#dogs are significant, hamsters are significant
population_summary <- df_pop %>%
  filter(species!="Ferret") %>%
  pivot_longer(!c("dataset_ID","species","pi", "piN_div_piS"),
             names_to="variable", values_to="value") %>%
  group_by(species, variable) %>%
  summarise(value=list(value)) %>%
  spread(variable, value)  %>%
  group_by(species) %>%
  mutate(p_value = t.test(unlist(piN), unlist(piS),paired=TRUE)$p.value,
         t_value = t.test(unlist(piN), unlist(piS),paired=TRUE)$statistic)

#now look at product summary data 
#this summarizes the data by species and gene, and then does a ttest between piN and piS for each gene
#orf1ab, spike, membrane are significant
product_summary <- df_prod %>%
  pivot_longer(!c("dataset_ID","gene_product", "species", "piN_div_piS"),
               names_to="variable", values_to="value") %>%
  group_by(gene_product, variable) %>%
  summarise(value=list(value)) %>%
  spread(variable, value) %>% 
  group_by(gene_product) %>%
  mutate(p_value = t.test(unlist(piN), unlist(piS),paired=TRUE)$p.value,
         t_value = t.test(unlist(piN), unlist(piS),paired=TRUE)$statistic)

product_summary %>%
  filter(p_value<=0.05)
#when you take all of the data together, M, orf1ab, and S are significant

#now we can plot these results
#(1) piN and piS for all datasets at a population level, and piN vs piS by species at a population level
#(2) piN vs piS by protein across all datasets (at a gene level)

#make another data frame with all the same data, but change the species to "all"
dfpop <- df_pop

dfpop <- dfpop %>%
  mutate(species="All")

#then join the data frames together, now you have all the data, both broken down by species and as "all"
df_pop <- full_join(dfpop, df_pop)

#(1)
options(scipen=8) #this helped with the scientific notation situation

plot1 <- ggpaired(df_pop, cond1 = "piN", cond2 = "piS",
         color = "condition", line.color = "gray", line.size = 0.4,
         palette="npg",xlab=FALSE, ylab=FALSE, legend="none", ylim=c(0, 0.00032))+
  facet_grid(cols=vars(species))+
  stat_compare_means(method="t.test", paired=TRUE, size=4, 
                     vjust=-2, label="p.format")

plot1

#(2) 
#pull out the interesting genes
df2 <- df_prod %>%
  filter(gene_product %in% c("S","orf1ab", "M"))

plot2 <- ggpaired(df2, cond1 = "piN", cond2 = "piS",
         color = "condition", line.color = "gray", line.size = 0.4,
         palette="npg",xlab=FALSE, ylab=FALSE, legend="none", ylim=c(0,0.0012))+
  facet_grid(cols=vars(gene_product))+
  stat_compare_means(method="t.test", paired=TRUE, size=4, 
                     vjust=-2, label="p.format")

plot2

pdf("selection.pdf", onefile=F)
ggarrange(plot1, plot2, ncol=1, nrow=2, labels="AUTO")
dev.off()
