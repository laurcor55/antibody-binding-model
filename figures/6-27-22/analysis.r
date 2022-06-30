library(jsonlite)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)

get_shortest_distance <- function(point_1, point_2){
  n_points = length(point_1)
  distance = vector("list", n_points)
  for (ii in 1:n_points){
    distance[[ii]] = norm(point_1[[ii]] - point_2, type="2")
  }
  distance <- min(unlist(distance))
  return(distance)
}

get_substrate_spacing <- function(substrates_location){
  if (length(substrates_location) < 2){
    distance <- 0
  }
  else if (length(substrates_location) == 2){
    distance <- norm(substrates_location[[1]] - substrates_location[[2]], type="2")
  }
  else{
    distance <- get_shortest_distance(substrates_location[2:length(substrates_location)], substrates_location[[1]])
  }
  return(distance)
}


import_data <- function(file_name){
  df <- fromJSON(file_name, simplifyVector=TRUE)
  df <- unnest_wider(df,  names_sep='_',col = substrates ) #
  df <- unnest(df, col=c(reaction, ligand), names_sep='_')
  df <- df %>% 
    mutate(
      p_bind = reaction_binding_count/reaction_total_count
    ) %>% 
    rowwise() %>%
    mutate(
      min_substrate_ligand_distance = get_shortest_distance(substrates_location, ligand_location),
      n_substrates = length(substrates_n_docks),
      min_substrate_substrate_distance = get_substrate_spacing(substrates_location),
      horizontal_offset = sum(abs(substrates_location[[1]] - c(0, 0, 0))),
      vertical_offset = sum(abs(ligand_location - c(0, 0, 0))),
      ligand_rotation_offset_deg = 180 - round(ligand_rotation[[2]]*180/3.14) + round(ligand_rotation[[3]]*180/3.14)
    )
  df <- as.data.frame(df)
  return(df)
}
df_1 <- import_data('output.json')
df_2 <- import_data('output_test.json')
df <- rbind(df_1, df_2)
glimpse(df)



df_0deg <- df %>% filter(ligand_rotation_offset_deg==0)
ggplot(data=df_0deg) +
  stat_summary(
    mapping=aes(x=min_substrate_ligand_distance, y=p_bind),
    fun=mean,
    geom="errorbar",
    width=0.1,
    fun.min=function(x) mean(x)-sd(x),
    fun.max=function(x) mean(x)+sd(x)) +
  facet_grid(cols=vars(min_substrate_substrate_distance)) + 
  xlab('Start Distance') +
  ylab('P(binding)') +
  labs(subtitle='Substrate Spacing Distance')+
  theme(plot.subtitle = element_text(hjust=0.5))
ggsave('p_bind_0deg.png', width=7, height=3)


df_summary <- df_0deg %>% 
  group_by(min_substrate_ligand_distance, min_substrate_substrate_distance) %>%
  summarize(
    center_binding = mean(reaction_center_binding_proportion*p_bind),
    p_bind = mean(p_bind)
  )
ggplot(data=df_summary) + 
  geom_bar(aes(x=min_substrate_ligand_distance, y=p_bind), stat='identity') + 
  geom_bar(aes(x=min_substrate_ligand_distance, y=center_binding, fill='red'), stat='identity') + 
  facet_grid(cols=vars(min_substrate_substrate_distance)) + 
  xlab('Start Distance') +
  ylab('P(binding)') +
  labs(subtitle='Substrate Spacing Distance', fill='Center Binding Proportion')+
  theme(plot.subtitle = element_text(hjust=0.5))
ggsave('center_binding_0deg.png', width=9, height=3)

ggplot(data=df) +
  stat_summary(
    mapping=aes(x=min_substrate_ligand_distance, y=reaction_center_binding_proportion),
    fun=mean,
    geom="errorbar",
    width=0.1,
    fun.min=function(x) mean(x)-sd(x),
    fun.max=function(x) mean(x)+sd(x)) + 
  xlab('Start Distance') +
  ylab('P(binding)') +
  theme(plot.subtitle = element_text(hjust=0.5))
ggsave('center_binding_0deg_norm.png', width=4, height=3)


df_summary <- df %>% 
  filter(min_substrate_substrate_distance == 36) %>%
  group_by(min_substrate_ligand_distance, ligand_rotation_offset_deg) %>%
  summarize(
    center_binding = mean(reaction_center_binding_proportion*p_bind),
    p_bind = mean(p_bind)
  )
ggplot(data=df_summary) + 
  geom_bar(aes(x=min_substrate_ligand_distance, y=p_bind), stat='identity') + 
  geom_bar(aes(x=min_substrate_ligand_distance, y=center_binding, fill='red'), stat='identity') + 
  facet_wrap(~ligand_rotation_offset_deg) +
  #facet_grid(cols=vars(ligand_rotation_offset_deg)) + 
  xlab('Start Distance') +
  ylab('P(binding)') +
  labs(subtitle='Ligand Rotation (deg)', fill='Center Binding Proportion')+
  theme(plot.subtitle = element_text(hjust=0.5))
ggsave('center_binding.png', width=10, height=6)