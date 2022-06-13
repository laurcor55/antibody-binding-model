library(jsonlite)
library(tidyverse)
library(ggplot2)

get_shortest_distance <- function(point_1, point_2){
  n_points = length(point_1)
  distance = vector("list", n_points)
  for (ii in 1:n_points){
    distance[[ii]] = norm(point_1[[ii]] - point_2, type="2")
  }
  distance <- min(unlist(distance))
  return(distance)
}

add_value <- function(list_input, int_input){
  output <- list_input + int_input
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
      n_substrates = length(substrates_n_docks)
    )
  return(df)
}
df_original <- import_data('output.json')
df_moveback <- import_data('output_moveback.json')

df_original <- df_original %>% mutate(overlap_strategy='original')
df_moveback <- df_moveback %>% mutate(overlap_strategy='moveback')
df <- rbind(df_original, df_moveback)

ggplot(data=df) +
  geom_point(aes(x=min_substrate_ligand_distance, y=p_bind, color=overlap_strategy))

ggplot(data=df) +
  geom_point(aes(x=min_substrate_ligand_distance, y=p_bind, color=n_substrates))

ggplot(data=df) + 
  geom_point(aes(x=min_substrate_ligand_distance, y=p_bind, color=overlap_strategy)) +
  facet_grid(cols=vars(n_substrates))

ggplot(data=df) + 
  geom_point(aes(x=min_substrate_ligand_distance, y=p_bind, color=n_substrates)) +
  facet_grid(cols=vars(overlap_strategy))


glimpse(df)