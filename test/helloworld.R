print("Hello world!")

install.packages("tidyverse")
library(tidyverse)

mpg %>%
  ggplot(aes(x = displ, y = hwy, color = class)) +
  geom_point()

install.packages("devtools")
library(devtools)
