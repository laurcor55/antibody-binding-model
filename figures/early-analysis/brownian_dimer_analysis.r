library(dplyr)
library(readr)

df <- read_csv('brownian_dimer.csv')
df_northrup <- read_csv('brownian_dimer_northrup.csv')
print(df)
print(df_northrup)