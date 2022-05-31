library(dplyr)
library(readr)
library(ggplot2)

df <- read_csv('k_off.csv')

ggplot(data=df) +
  geom_point(aes(x=distance, y=k_off)) +
  ylab(expression('k'['off']*' (s'^-1*')')) + 
  xlab('Off distance (Å)')

#ggsave('k_off.png', width=4, height=3)

ggplot(data=df) +
  geom_point(aes(x=distance, y=`unbinding events`)) +
  ylab('Unbinding events') + 
  xlab('Off distance (Å)')

#ggsave('unbinding_events.png', width=4, height=3)


df <- read_csv('k_off_2.csv')

ggplot(data=df) +
  geom_point(aes(x=distance, y=k_off)) +
  ylab(expression('k'['off']*' (s'^-1*')')) + 
  xlab('Off distance (Å)')

ggsave('k_off_2.png', width=4, height=3)