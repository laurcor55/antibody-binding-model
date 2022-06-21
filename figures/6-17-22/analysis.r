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
      vertical_offset = sum(abs(ligand_location - c(0, 0, 0)))
    )
  return(df)
}
df_original <- import_data('output.json')
df_moveback <- import_data('output_moveback.json')

df_original <- df_original %>% mutate(overlap_strategy='original')
df_moveback <- df_moveback %>% mutate(overlap_strategy='moveback')
df <- rbind(df_original, df_moveback)
df <- as.data.frame(df)

df_3_docks <- df %>% filter(ligand_n_docks == 3 & n_substrates < 3)


ggplot(data=df_3_docks) + 
  geom_point(aes(x=min_substrate_ligand_distance, y=p_bind, color=overlap_strategy)) +
  facet_grid(cols=vars(n_substrates), labeller=label_both) 
ggsave('overlap_method.png', width=8, height=3)


df <- df %>% filter(ligand_n_docks == 4) 
df_0deg <- df %>% rowwise() %>% filter(sum(abs(ligand_rotation - c(0, 3.14, 0)))<0.01)
df_90deg <- df %>% rowwise() %>% filter(sum(abs(ligand_rotation - c(0, 3.14/2, 0)))<0.01)


df_plot <- df_0deg %>% 
  filter(min_substrate_substrate_distance==40 ||min_substrate_substrate_distance==0 ) %>% 
  rowwise() %>%
  filter(sum(abs(substrates_location[[1]]-c(0, 0, 0)))<0.01)
ggplot(data=df_plot) +
  geom_point(aes(x=min_substrate_ligand_distance, y=p_bind)) +
  facet_grid(cols=vars(n_substrates))
ggsave('n_substrates.png', height=3, width=12)


df_plot <- df_0deg %>% 
  filter(n_substrates==9) %>% 
  rowwise() %>%
  filter(sum(abs(substrates_location[[1]]-c(0, 0, 0)))<0.01)
ggplot(data=df_plot) +
  geom_point(aes(x=min_substrate_ligand_distance, y=p_bind))+
  facet_grid(cols=vars(min_substrate_substrate_distance))
ggsave('substrate_spacing.png', height=3, width=12)

ggplot(data=df_0deg) + 
  geom_point(aes(x=min_substrate_ligand_distance, y=p_bind)) +
  facet_grid(cols=vars(min_substrate_substrate_distance), rows=vars(n_substrates))

df_plot <- df_90deg %>% 
  rowwise() %>%
  filter(sum(abs(substrates_location[[1]]-c(0, 0, 0)))<0.01)
ggplot(data=df_plot) +
  geom_point(aes(x=min_substrate_ligand_distance, y=p_bind, color=factor(n_substrates)))
ggsave('rotated_ligand.png', height=3, width=5)

df_plot <- df_0deg %>% 
  filter(n_substrates == 1)
ggplot(data=df_plot) +
  geom_point(aes(x=min_substrate_ligand_distance, y=p_bind)) + 
  facet_grid(cols=vars(horizontal_offset))
ggsave('horizontal_offset.png', height=3, width=12)

df_plot <- df %>% rowwise() %>% 
  filter(sum(abs(ligand_rotation - c(0, 3.14, 3.14/4)))<0.01)
ggplot(data=df_plot) +
  geom_point(aes(x=horizontal_offset, y=p_bind)) +
  facet_grid(cols=vars(vertical_offset))
ggsave('vertical_offset_tilt.png', height=3, width=7)



glimpse(df)