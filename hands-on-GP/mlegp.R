library(mlegp)
library(tidyverse)
load("toJuho.Rdata")
load("whole_stack.Rdata")
load("metrics_matrix.Rdata")
source("model_metrics.R")

##### calculating model metrics via function model_metrics #####

mod <- model_metrics(n_knots = 50)

mod <- model_metrics(site_no = 1, response = 1, n_knots = 25, n_test_points = 150, predictor_type = "quantile") # example

model_metrics(1,1,50,150,"mlegp","normal",TRUE,sqrt(.Machine$double.eps),shuffled_df)
model_metrics(1,1,50,150,"mlegp","normal",FALSE,sqrt(.Machine$double.eps),shuffled_df)

df <- as.data.frame(SS.stack[[1]][[1]])
colnames(df)[6] <- "likelihood"
df <- df[sample(1:nrow(df), size = nrow(df), replace = FALSE), ]

model_metrics(1,1,50,150,"mlegp","normal",TRUE,sqrt(.Machine$double.eps),shuffled_df)
model_metrics(1,1,50,150,"mlegp","normal",TRUE,0,shuffled_df)
model_metrics(site_no = 1, response = 1, n_knots = 50, n_test_points = 150, package = "mlegp", predictor_type = "normal",
              nugget_known = TRUE, nugget = 0.1, data = shuffled_df)


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


rep <- 25

df_metrics <- data.frame(RMSqE = 1, MStdE = 1, time = 0.1, nan_var = 0.5, nugget = "est")
df_metrics <- df_metrics[-1, ]
df_metrics$nugget <- as.factor(df_metrics$nugget)
levels(df_metrics$nugget) <- c("eps", "0.1", "0.001", "est")

for (i in 1:rep) {
  df <- df[sample(1:nrow(df), size = nrow(df), replace = FALSE), ]
  metrics_eps <- model_metrics(site_no = 1, response = 1, n_knots = 50, n_test_points = 150, package = "mlegp",
                               predictor_type = "normal", nugget_known = TRUE, nugget = sqrt(.Machine$double.eps), data = df)
  metrics_01 <- model_metrics(site_no = 1, response = 1, n_knots = 50, n_test_points = 150, package = "mlegp",
                              predictor_type = "normal", nugget_known = TRUE, nugget = 0.1, data = df)
  metrics_0001 <- model_metrics(site_no = 1, response = 1, n_knots = 50, n_test_points = 150, package = "mlegp",
                                predictor_type = "normal", nugget_known = TRUE, nugget = 0.001, data = df)
  metrics_est <- model_metrics(site_no = 1, response = 1, n_knots = 50, n_test_points = 150, package = "mlegp",
                               predictor_type = "normal", nugget_known = FALSE, data = df)
  
  
  df_metrics[nrow(df_metrics) + 1 , ] <- c(metrics_eps[[3]], metrics_eps[[4]], metrics_eps[[1]]+metrics_eps[[2]],metrics_eps[[5]], "eps")
  df_metrics[nrow(df_metrics) + 1 , ] <- c(metrics_01[[3]], metrics_01[[4]], metrics_01[[1]]+metrics_01[[2]],metrics_01[[5]], "0.1")
  df_metrics[nrow(df_metrics) + 1 , ] <- c(metrics_0001[[3]], metrics_0001[[4]], metrics_0001[[1]]+metrics_0001[[2]],metrics_0001[[5]], "0.001")
  df_metrics[nrow(df_metrics) + 1 , ] <- c(metrics_est[[3]], metrics_est[[4]], metrics_est[[1]]+metrics_est[[2]],metrics_est[[5]], "est")
  
  
}
for (i in 1:4) {
  df_metrics[,i] <- as.numeric(df_metrics[,i])
}

df_metrics %>%
  group_by(nugget) %>%
  summarise(RMSqE = mean(RMSqE),
            MStdE = mean(MStdE, na.rm = TRUE),
            time = mean(time),
            'NaN%' = mean(nan_var)) %>% as.data.frame()


rmse_normal <- rep(0, 100)
rmse_quantile <- rep(0,100)
mstde_normal <- rep(0,100)
mstde_quantile <- rep(0,100)
time_normal <- rep(0,100)
time_quantile <- rep(0,100)

for (i in 1:100) {
  df <- df[sample(1:nrow(df), size = nrow(df), replace = FALSE), ]
  metrics_normal <- model_metrics(site_no = 1, response = 1, n_knots = 50, n_test_points = 150,
                                  package = "mlegp", predictor_type = "normal", data = df)
  metrics_quantile <- model_metrics(site_no = 1, response = 1, n_knots = 50, n_test_points = 150,
                                    package = "mlegp", predictor_type = "quantile", data = df)
  rmse_normal[i] <- metrics_normal[[3]]
  rmse_quantile[i] <- metrics_quantile[[3]]
  mstde_normal[i] <- metrics_normal[[4]]
  mstde_quantile[i] <- metrics_quantile[[4]]
  time_normal[i] <- metrics_normal[[1]] + metrics_normal[[2]]
  time_quantile[i] <- metrics_quantile[[1]] + metrics_quantile[[2]]
}

df_metrics <- data.frame(rmse_normal = rmse_normal, rmse_quantile = rmse_quantile, mstde_normal = mstde_normal,
                         mstde_quantile = mstde_quantile, time_normal = time_normal, time_quantile = time_quantile)

model_metrics(1,1,50,150,"mlegp","normal")
