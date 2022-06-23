library(mlegp)
library(tidyverse)
load("toJuho.Rdata")
load("whole_stack.Rdata")
load("metrics_matrix.Rdata")
source("model_metrics.R")

##### calculating model metrics via function model_metrics #####

mod <- model_metrics(n_knots = 50)

mod <- model_metrics(site_no = 1, response = 1, n_knots = 25, n_test_points = 150, predictor_type = "quantile") # example

df <- as.data.frame(SS.stack[[1]][[1]])
colnames(df)[6] <- "likelihood"
df <- df[sample(1:nrow(df), size = nrow(df), replace = FALSE), ]

### testing if there's difference between using normal values for predictors or quantile values
rmse_normal <- rep(0, 10)
rmse_quantile <- rep(0,10)
mstde_normal <- rep(0,10)
mstde_quantile <- rep(0,10)
time_normal <- rep(0,10)
time_quantile <- rep(0,10)

for (i in 1:10) {
  df <- df[sample(1:nrow(df), size = nrow(df), replace = FALSE), ]
  metrics_normal <- model_metrics(site_no = 1, response = 1, n_knots = 75, n_test_points = 150, predictor_type = "normal", data = df)
  metrics_quantile <- model_metrics(site_no = 1, response = 1, n_knots = 75, n_test_points = 150, predictor_type = "quantile", data = df)
  rmse_normal[i] <- metrics_normal[[3]]
  rmse_quantile[i] <- metrics_quantile[[3]]
  mstde_normal[i] <- metrics_normal[[4]]
  mstde_quantile[i] <- metrics_quantile[[4]]
  time_normal[i] <- metrics_normal[[1]] + metrics_normal[[2]]
  time_quantile[i] <- metrics_quantile[[1]] + metrics_quantile[[2]]
}

df_metrics <- data.frame(rmse_normal = rmse_normal, rmse_quantile = rmse_quantile, mstde_normal = mstde_normal,
                         mstde_quantile = mstde_quantile, time_normal = time_normal, time_quantile = time_quantile)



### on Thursday run a couple more for site 1, largest no_knots
for (n in c(200, 300, 400, 500, 600)) {
  for (resp in c(1,2)) {
    output <- model_metrics(site_no = 12, response = resp, n_knots = n, n_test_points = 150)
    metrics_matrix[nrow(metrics_matrix)+1, ] <- c(12, resp, n, 150, output[[1]], output[[2]], output[[3]], output[[4]])
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
  filter(knots == 100, response == 1) %>%
  ggplot(mapping = aes(x = time_model+time_pred, y=RMSqE, color=site)) +
  geom_point(size = 3) + labs(x = "time", title = "Site-to-site differences in prediction accuracies",
                              subtitle = "100 knots used") +
  guides(color = guide_legend(ncol=2)) +
  scale_color_brewer(palette = "Paired")

plot12 <- metrics_matrix %>%
  filter(knots == 100, response == 1) %>%
  ggplot(mapping = aes(x = time_model+time_pred, y=MStdE, color=site)) +
  geom_point(size = 3) + labs(x = "time") +
  guides(color = guide_legend(ncol=2)) +
  scale_color_brewer(palette = "Paired")

plot21 <- metrics_matrix%>%
  filter(site == 1) %>%
  ggplot(mapping = aes(x = time_model+time_pred, y=RMSqE, color=knots, shape=response)) +
  geom_point(size = 3) + labs(x = "time", title = "Effects of increasing number of knots using site No.1") +
  ylim(0, 1500)

plot22 <- metrics_matrix %>%
  filter(site == 1) %>%
  ggplot(mapping = aes(x = time_model+time_pred, y=MStdE, color=knots, shape=response)) +
  geom_point(size = 3) + labs(x = "time") +
  ylim(0,600)


(plot1 <- ggpubr::ggarrange(plot11, plot12, nrow = 2, common.legend = TRUE, legend = "right"))
(plot2 <- ggpubr::ggarrange(plot21, plot22, nrow = 2, common.legend = TRUE, legend = "right"))


