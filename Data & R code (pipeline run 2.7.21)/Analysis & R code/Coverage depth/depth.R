##mining more information from the depth file ouput by the pipeline
#this script is mostly the same, some manuscript-specific modifications to the depth script that is part of the pipeline

#load libraries
library(tidyverse)
library(openxlsx)
library(pdftools)

#load data, rename columns
depth_df <- read.delim("all.depth", sep="\t", header=FALSE)
colnames(depth_df) <- c("dataset", "virus", "position", "depth")

# rename datasets 
#source(paste0(r_bindir, "/process_dataset_names.R"))
depth_df <- process_dataset_names(depth_df)

# higlight coverage below a certain limit
# TODO: Parameterize this
min_depth_highlight <- 100
depth_df <- depth_df %>% mutate(above_highlight = if_else(depth > min_depth_highlight, TRUE, FALSE))

# max position for drawing a red rectangle showing low coverage
max_position = max(depth_df$position)

######### calculate average depth in windows
# %/% is the integer division operator
window_size = 10
depth_df <- depth_df %>% mutate (window = position %/% window_size)

# calculate average coverage depth in each window
df_windowed <- depth_df %>% 
  group_by(dataset, virus, window)  %>% 
  summarize(depth = mean(depth), .groups = "drop") %>% 
  mutate(position = (window*window_size) + 1) %>% 
  ungroup()

##now plot coverage data on multiple pdf pages

datasets <- depth_df %>% group_by(dataset) %>% summarize(.groups="drop") %>% pull()

plot_some_datasets <- function(datasets){

  datasets_per_page <- 12

  page_number <- 1

  # iterate through the viruses, doing so many per page
  for (i in seq(1, length(datasets), datasets_per_page)) {

    # plots_per_page at a time
    subset_datasets <- datasets[i:(i+(datasets_per_page-1))]

    # output to console which ones we're doing
    print(paste0(subset_datasets))

    pdf_name <- paste0("coverage_plot_page_", page_number, ".pdf")

    # generate & print plot
    plot_datasets(subset_datasets, pdf_name)

    page_number = page_number + 1
    
  }

}

plot_datasets <- function(dataset_names, pdf_name){
  
  # subset the main dataframes to get the data just for these viruses
  subset_df <- df_windowed %>% filter(dataset %in% dataset_names)
  
  # output the virus names to the console
  print(paste0(dataset_names))
  
  # p returned here is a ggplot object
  p <- ggplot(subset_df) + 
    geom_line(aes(x=position, y=depth), size=0.5) +
    geom_area(aes(x=position, y=depth), fill="lightgrey") +
    annotate("rect", xmin=0, xmax=max_position, ymin=1, ymax=min_depth_highlight, alpha=0.05, fill="red") +
    scale_color_manual(values = c("red", "black")) +
    theme_bw(base_size = 10) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank()) +
    scale_y_log10() +
    xlab("genome position (nt)") +
    ylab (paste0("coverage depth\n", "coverage below ", min_depth_highlight, " in red"))+
    facet_grid(dataset~virus, scales="free", space="free_x")
  
  ggsave(pdf_name, p, height=10.5, width=7.5, units="in")

    # print p will make the plots appear on viewer
    # print(p)

  }

plot_some_datasets(datasets)

pdf_combine(c(list.files(pattern="pdf$")), output="coverage_plot.pdf")

# calculated median depth of (total) coverage for each virus in each dataset and store it in a new df
median_depths <- depth_df %>% 
  group_by(dataset, virus) %>% 
  summarize(median_depth = median(depth))

median_depth_all_together <- depth_df %>%
  summarize(median_depth = median(depth))

mean_depths <- depth_df %>% 
  group_by(dataset, virus) %>% 
  summarize(mean_depth = mean(depth))

mean_depth_all_together <- depth_df %>%
  summarize(mean_depth = mean(depth))

wb <- createWorkbook("Median_depths.xlsx")
addWorksheet(wb, "median_depth")
addWorksheet(wb, "median_depth_all_together")
addWorksheet(wb, "mean_depth")
addWorksheet(wb, "mean_depth_all_together")
writeData(wb, "median_depth", median_depths)
writeData(wb, "median_depth_all_together", median_depth_all_together)
writeData(wb, "mean_depth", mean_depths)
writeData(wb, "mean_depth_all_together", mean_depth_all_together)
saveWorkbook(wb, "Median_depths.xlsx", overwrite = TRUE)



#look at important spike and other variants from Table 2
#make sure there was good depth of coverage, no chance that low coverage skewed the assessment high frequency variants

artic_v2 <- depth_df %>%
  filter(dataset %in% c("Cat_5", "Cat_4")) %>%
  filter(position %in% c(21768, 23403, 441, 3303, 4965, 8240,11083, 18763,21974, 22205, 23064, 23618, 28021, 28285)) %>%
  arrange(depth)
  
min(artic_v2$depth) #20x for N501T position, but that's below our cutoff and no animals other than ferret had a variant there

#the min for a table 2 variant actually detected in the v2 cats would be 95x at H69R for Cat 5 -- BUT this is validated by Cat 6 depth of coverage at this position
#the max is 13047x (Cat 4 at the N4N position)

#for the rest of the variants

artic_v3 <- depth_df %>%
  filter(!(dataset %in% c("Cat_5", "Cat_4"))) %>%
  filter(position %in% c(21768, 23403, 441, 3303, 4965, 8240,11083, 18763,21974, 22205, 23064, 23618, 28021, 28285)) %>%
  arrange(depth)

#the min for table 2 variants detected in v3 animals is 683x at D215 for Dog 2
#the max for v3 sequenced animals is 87274x for Cat 2 at N4N

#variant by variant summary

h69r <- depth_df %>%
  filter(position==21768) #Cat_5 is 95x, but Cat_6 & Cat 6_R were 3843 & 1722x

d614g <- depth_df %>%
  filter(position==23403) #Cat_5 is 622x, Cat_6/R were 2959x and 1681x

g59d <- depth_df %>%
  filter(position==441)

e195g <- depth_df %>%
  filter(position==3303)

t749i <- depth_df %>%
  filter(position==4965)

h1841y <- depth_df %>%
  filter(position==8240)

l37f <- depth_df %>%
  filter(position==11083)

i242v <- depth_df %>%
  filter(position==18763)

d138y <- depth_df %>%
  filter(position==21974)

d215n <- depth_df %>%
  filter(position==22205)

n501t <- depth_df %>%
  filter(position==23064)

s686g <- depth_df %>%
  filter(position==23618)

s43y <- depth_df %>%
  filter(position==28021)

n4n <- depth_df %>%
  filter(position==28285)

#h69r, d614g, g59d, e195g, t749i, h1841y, l37f, i242v, d138y, d215n, n501t, s686g, s43y, n4n


#TC variants
#already did D215

r685h <- depth_df %>%
  filter(position==23616)
