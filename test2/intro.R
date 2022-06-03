library(devtools)
library(tidyverse)

mpg %>%
  ggplot(aes(x = displ, y = hwy, color = class)) +
  geom_point()
