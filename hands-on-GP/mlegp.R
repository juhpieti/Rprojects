library(mlegp)
library(tidyverse)
load("toJuho.Rdata")
load("whole_stack.Rdata")
load("metrics_matrix.Rdata")
source("model_metrics.R")

##### calculating model metrics via function model_metrics #####

model_metrics(site_no = 2, response = 1, n_knots = 25, n_test_points = 150) # example

### on Thursday run a couple more for site 1, largest no_knots
for (n in c(500, 600)) {
  for (resp in c(1,2)) {
    output <- model_metrics(site_no = 1, response = resp, n_knots = n, n_test_points = 150)
    metrics_matrix[nrow(metrics_matrix)+1, ] <- c(1, resp, n, 150, output[[1]], output[[2]], output[[3]], output[[4]])
  }
}

###for different sites with n_knots = 100
for (site in 7:12) {
  for (resp in c(1,2)) {
    output <- model_metrics(site_no = site, response = resp, n_knots = 100, n_test_points = 150)
    metrics_matrix[nrow(metrics_matrix)+1, ] <- c(site, resp, 100, 150, output[[1]], output[[2]], output[[3]], output[[4]])
  }
} 


save(metrics_matrix, file = "metrics_matrix.Rdata")

##### TESTING AREA #####


plot11 <- metrics_matrix %>%
  filter(knots == 100) %>%
  ggplot(mapping = aes(x = time_model+time_pred, y=RMSqE, color=site, shape=response)) +
  geom_point(size = 3) + labs(x = "time", title = "Site-to-site differences in prediction accuracies") +
  guides(color = guide_legend(ncol=2))

plot12 <- metrics_matrix %>%
  filter(knots == 100) %>%
  ggplot(mapping = aes(x = time_model+time_pred, y=MStdE, color=site, shape=response)) +
  geom_point(size = 3) + labs(x = "time") +
  guides(color = guide_legend(ncol=2))

plot21 <- metrics_matrix%>%
  filter(site == 1) %>%
  ggplot(mapping = aes(x = time_model+time_pred, y=RMSqE, color=knots, shape=response)) +
  geom_point(size = 3) + labs(x = "time", title = "Effects of increasing number of knots using site No.1")

plot22 <- metrics_matrix %>%
  filter(site == 1) %>%
  ggplot(mapping = aes(x = time_model+time_pred, y=MStdE, color=knots, shape=response)) +
  geom_point(size = 3) + labs(x = "time")


(plot1 <- ggpubr::ggarrange(plot11, plot12, nrow = 2, common.legend = TRUE, legend = "right"))
(plot2 <- ggpubr::ggarrange(plot21, plot22, nrow = 2, common.legend = TRUE, legend = "right"))