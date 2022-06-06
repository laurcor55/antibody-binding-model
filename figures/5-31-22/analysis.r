library(readr)
library(dplyr)
library(ggplot2)
df <- read_csv('k_on_distances.csv')
df <- mutate(df, dock_distance = ligand_z - substrate_z, docks='1')

ggplot(data=df) + 
  geom_point(aes(x=dock_distance, y=p_on)) + 
  scale_y_log10()
ggsave('p_on.png', width=4, height=3)

ggplot(data=df) + 
  geom_point(aes(x=dock_distance, y=k_on)) + 
  scale_y_log10() 
ggsave('k_on.png', width=4, height=3)

df_two_docks <- read_csv('k_on_distances_two_docks.csv')
df_two_docks <- mutate(df_two_docks, dock_distance = ligand_z - substrate_z, docks='2')
df <- bind_rows(df, df_two_docks)

ggplot(data=df) + 
  geom_point(aes(x=dock_distance, y=p_on, color=docks)) + 
  scale_y_log10()
ggsave('p_on_docks.png', width=4, height=3)

ggplot(data=df) + 
  geom_point(aes(x=dock_distance, y=k_on, color=docks)) + 
  scale_y_log10() 
ggsave('k_on_docks.png', width=4, height=3)