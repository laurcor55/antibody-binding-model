library(jsonlite)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)

df <- read_csv('k_on_ligand_random_rotation.csv')

glimpse(df)
ggplot(data=df) + 
  geom_point(aes(x=n_substrates, y=p_binding)) 

ggplot(data=df) + 
  geom_point(aes(x=n_substrates, y=k_on))