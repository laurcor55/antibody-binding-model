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
  df <- filter(df, ligand_n_docks == 4)
  df <- as.data.frame(df)
  return(df)
}
df <- import_data('output.json')
glimpse(df)

ggplot(data=df) +
  geom_point(aes(x=min_substrate_ligand_distance, y=p_bind, color=factor(ligand_rotation_offset_deg)))
ggsave('ligand_rotation.png', width=7, height=3)

ggplot(data=df) +
  geom_point(aes(x=min_substrate_ligand_distance, y=reaction_center_binding_proportion, color=factor(ligand_rotation_offset_deg)))
ggsave('center_binding_proportion.png', width=7, height=3)
