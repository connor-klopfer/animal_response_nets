#' Visualistions 
library(dplyr)
library(ggplot2) 
source("analysis_scripts/general_analysis.R", echo = F)


fix_timepoints <- function(d_f){
  cond_1 <- c("low_risk", "before", "Pre")
  cond_2 <- c("high_risk", "after", "Post")
  d_f$numerical_condition <- NA
  d_f$numerical_condition[d_f$condition %in% cond_1] <- 1
  d_f$numerical_condition[d_f$condition %in% cond_2] <- 2
  return(d_f)
}


get_difference <- function(d_f){
  #' Get the difference of values the 
  return(d_f %>% 
    fix_timepoints() %>% 
    select(-network, -condition) %>% 
    # group_by(animal, replicate_num) %>%
    arrange(animal, replicate_num, numerical_condition) %>% 
    mutate_at(vars(-animal, -replicate_num, -numerical_condition), difference_calc) %>% 
    filter(numerical_condition == 1))
}

difference_calc <- function(x){
  return(as.numeric(x) - as.numeric(lag(x)))
}

# For testing. 
# arrange(animal, replicate_num, numerical_condition) %>% View()

plot_comparison_plots <- function(d_f){
  
  final_plot <- d_f %>% 
    fix_timepoints() %>% 
    tidyr::pivot_longer(c("density", "eigenvector_centrality", 
                                "betweenness_centrality","n_nodes", 
                                "effective_information"), names_to = "metric", 
                              values_to = "value") %>% 
    ggplot(aes(x = as.factor(numerical_condition), y = as.numeric(value), fill = animal))+
    geom_boxplot()+
    facet_grid(metric~animal, scale = "free")+
    labs(x = "Pre (1) or Post (2) Stimulus", y = "Value", 
         title = "Pre and Post Stimuls in animal social networks.")
  
  ggsave(final_plot, path = "figures", height = 8, width = 6,
         filename = "comparing_pre_post_boxplot.png", device = "png")
}

plot_comparison_plots(final_df)

plot_all_change <- function(d_f){
  final_plot <- d_f %>% 
    get_difference() %>%
    tidyr::pivot_longer(c("density", "eigenvector_centrality", 
                          "betweenness_centrality","n_nodes", 
                          "effective_information"), names_to = "metric", 
                        values_to = "value") %>% 
    ggplot(aes(x = animal, y = as.numeric(value), fill = animal))+
    geom_boxplot()+
    theme(axis.text.x = element_text(angle = 20, vjust = 1, hjust = 1))+
    facet_wrap(.~metric, scales = "free") +
    labs(title = "Comparing Animal Social Networks", x = "", 
         y = "Change in Value", fill = "Species")

  ggsave(final_plot, path = "figures", filename = "all_change_boxplot.png", device = "png")
  
} 

plot_all_change(final_df)

all_individual_plots <- function(){
  
  final_df %>%
    ggplot(aes(x = condition, y = as.numeric(eigenvector_centrality), 
               fill = condition))+
    geom_violin(alpha = 0.4) + 
    geom_boxplot(width = 0.1, size = 1.5, fill = NA)+
    facet_grid(.~animal, scales = "free")+
    labs(title = "Comparing Eigenvector Centrality", 
         y = "Eigenvector Centrality")
  
  final_df %>%
    ggplot(aes(x = condition, 
               y = as.numeric(betweenness_centrality), 
               fill = condition))+
    geom_violin(alpha = 0.4) + 
    geom_boxplot(width = 0.1, size = 1.5, fill = NA)+
    facet_grid(.~animal, scales = "free")+
    labs(title = "Comparing Betweenness Centrality", 
         y = "Betweenness Centrality")
  
  final_df %>%
    ggplot(aes(x = condition, y = as.numeric(density), 
               fill = condition))+
    geom_boxplot(width = 0.1, size = 1.5, fill = NA)+
    facet_grid(.~animal, scales = "free")+
    labs(title = "Comparing Density", 
         y = "Density")
  
  final_df %>%
    ggplot(aes(x = condition, y =as.numeric(n_nodes), 
               fill = condition))+
    # geom_bar(alpha = 0.4, stat = "identity") + 
    geom_boxplot(width = 0.1, size = 1.5, fill = NA)+
    facet_grid(.~animal, scales = "free")+
    labs(title = "Comparing Population (Number of Nodes)", 
         y = "NUmber of Nodes")
  
  final_df %>% 
    ggplot(aes(x = condition, y =as.numeric(effective_information), 
               fill = condition))+
    # geom_bar(alpha = 0.4, stat = "identity") + 
    geom_violin(alpha = 0.4) + 
    geom_boxplot(width = 0.1, size = 1.5, fill = NA)+
    facet_grid(.~animal, scales = "free")+
    labs(title = "Comparing Effective Information Spread", 
         y = "Effective Information")
  
  ### COmparing the change in network metrics. 
  
  final_df %>% 
    get_difference() %>% 
    ggplot(aes(x = animal, y = as.numeric(eigenvector_centrality), 
               fill = animal))+
    geom_violin(alpha = 0.4) + 
    geom_boxplot(width = 0.1, size = 1.5, fill = NA)+
    # facet_grid(.~animal, scales = "free")+
    labs(title = "Change in Eigenvector Centrality", 
         y = "Eigenvector Centrality")
  
  final_df %>%
    get_difference() %>% 
    ggplot(aes(x = animal, 
               y = as.numeric(betweenness_centrality), 
               fill = animal))+
    geom_violin(alpha = 0.4) + 
    geom_boxplot(width = 0.1, size = 1.5, fill = NA)+
    # facet_grid(.~animal, scales = "free")+
    labs(title = "Change in Betweenness Centrality", 
         y = "Betweenness Centrality")
  
  final_df %>%  get_difference() %>% 
    ggplot(aes(x = animal, y = as.numeric(density), 
               fill = animal))+
    geom_boxplot(width = 0.1, size = 1.5, fill = NA)+
    # facet_grid(.~animal, scales = "free")+
    labs(title = "Change in Density", 
         y = "Density")
  
  final_df %>%  get_difference() %>%
    ggplot(aes(x = animal, y =as.numeric(n_nodes), 
               fill = animal))+
    # geom_bar(alpha = 0.4, stat = "identity") +
    geom_boxplot(width = 0.1, size = 1.5)+
    # facet_grid(.~animal, scales = "free")+
    labs(title = "Change in Population (Number of Nodes)", 
         y = "NUmber of Nodes")
  
  final_df %>% 
    get_difference() %>% 
    ggplot(aes(x = animal, y =as.numeric(effective_information), 
               fill = animal))+
    # geom_bar(alpha = 0.4, stat = "identity") + 
    geom_violin(alpha = 0.4) + 
    geom_boxplot(width = 0.1, size = 1.5, fill = NA)+
    # facet_grid(.~animal, scales = "free")+
    labs(title = "Change in Effective Infomation Spread", 
         y = "Effective Information")
}
