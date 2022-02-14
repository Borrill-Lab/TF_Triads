#'TF evolution analysis
#'7/9/21
#'Let's take coding_variants_chisquared_step6_c.R and make the plotting and chi squared into functions
#'which I can call with as many data frame + variables combos as I want
#'This script contains the functions
#'
#'21/9/21
#'Add optional coord_flip() to make.stacked.plot()
#'
#'30/01/22
#'no legend - add legend afterwards in Inkscape
#'
#'09/02/22
#'Add font sizes to base theme
#'
#'Built with R 4.0.5 and tidyverse

#PREP
require(tidyverse)
require(lubridate)
require(chisq.posthoc.test)
theme_set(theme_bw(base_size = 9))
#Get font sizes closer to that required for paper - but it's not possible to ask for Arial font
theme_update(axis.title = element_text(size = 10), axis.text = element_text(size = 9), legend.text = element_text(size = 9))

#FUNCTIONS
#'This should make a bar plot counting SNPs in each x category
make.barplot <- function(data, x = "Consequence", labels = labs(y= "Number of SNPs", x = "Predicted SNP effect"), description = "")
  {
  barplot <- ggplot(data, aes_string(x=x)) +
    geom_bar(stat="count") +
    facet_wrap(facets = vars(isTF), scale = "free_y") +
    theme(axis.text.x = element_text(angle = 90)) +
    ggtitle(description) +
    labels
  setwd("/rds/projects/b/borrillp-phd-students/Catherine/TF_homoeolog_SNPs/figures")
  ggsave(file = paste("barplot_of_", description, "_", x, "_", lubridate::today(),".pdf", sep = ""), plot = barplot, device = "pdf", width = 20, height = 20, units = "cm")
}

#'This makes a stacked plot where categories sit on top of each other
make.stacked.plot <- function(data, x = "isTF", fill = "Consequence", description = "", invert = FALSE)
{
  stacked_plot <- ggplot(data, aes_string(x=x, fill=fill)) +
    geom_bar(position = "fill") +
    ggtitle(description) +
    theme(panel.grid = element_blank(), legend.position = "right", plot.margin = unit(c(1,1,1,1), "cm"), line = element_line(colour = "black"), axis.text = element_text(colour = "black"))
  if(invert == TRUE){
    stacked_plot <- stacked_plot + coord_flip()
  }
  setwd("/rds/projects/b/borrillp-phd-students/Catherine/TF_homoeolog_SNPs/figures")
  ggsave(file = paste("stacked_plot_of_", description, "_", x, "_", lubridate::today(),".pdf", sep = ""), plot = stacked_plot, device = "pdf", width = 20, height = 20, units = "cm")
}

#'Return stacked plot - working function
return.stacked.plot <- function(data, x = "isTF", fill = "Consequence", description = "", invert = FALSE, legend.position = "none")
{
  stacked_plot <- ggplot(data, aes_string(x=x, fill=fill)) +
    geom_bar(position = "fill") +
    scale_y_continuous(expand = c(0.005,0.005)) + #remove whitespace at top and bottom of stacked bar
    ggtitle(description) +
    theme(panel.grid = element_blank(), legend.position = legend.position, line = element_line(colour = "black"), axis.text = element_text(colour = "black"))
  if(invert == TRUE){
    stacked_plot <- stacked_plot +
      geom_text(aes(label = ..count.., fill = NULL, group = NULL), stat = "count", position = position_fill(), hjust = 2) +
      coord_flip() #Set categories on vertical axis
  }else{
    stacked_plot <- stacked_plot +
      geom_text(aes(label = ..count.., fill = NULL, group = NULL), stat = "count", position = position_fill(), vjust = 2)
    #geom_text specified separately to use hjust or vjust based on axis angle
  }
  
  return(stacked_plot) #return whole plot so more things can be added later
}


#'Make a contingency table, run and print a chi squared test.
make.contingency.table.1 <- function(data, xvar, yvar)
  {
    contingency_table <- data %>%
      group_by(data[[xvar]], data[[yvar]]) %>% #the factors to group by
      count() %>%
      pivot_wider(names_from = `data[[xvar]]`, values_from = n) %>% #turn it into contingency table format
      ungroup %>%
      select(-(`data[[yvar]]`)) #not the information column
    return(contingency_table)
}

make.contingency.table.2 <- function(data, xvar1, xvar2, yvar)
{
  contingency_table <- data %>%
    group_by(data[[xvar1]], data[[xvar2]], data[[yvar]]) %>% #the factors to group by
    count() %>%
    pivot_wider(names_from = c(`data[[xvar1]]`, `data[[xvar2]]`), values_from = n) %>% #turn it into contingency table format
    ungroup %>%
    select(-(`data[[yvar]]`)) #not the information column
  return(contingency_table)
}

run.chisquared.test <- function(contingency_table)
{
  xsq <- chisq.test(contingency_table)
  print("CHI SQUARED TEST")
  print(xsq)
  print("Expected values")
  print(xsq$expected)
  print("Observed values")
  print(xsq$observed)
}




#'Calculate the proportion of synonymous variants in a vector, for use with fct_reorder()
prop.synonymous <- function(vector){
  return(sum(vector=="synonymous_variant.na")/length(vector))
}

#'Calculate the proportion of deleterious missense or stop variants in a vector, for use with fct_reorder()
prop.deleterious <- function(vector){
  return(sum(vector %in% c("stop_gained.na", "missense_variant.deleterious"))/length(vector))
}
