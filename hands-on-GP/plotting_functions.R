library(tidyverse)
load("mlegp_matrix.Rdata")

mean_and_plot_rmsqe <- function(matrix, pred, site, resp, ymin = 0, ymax = NA) {
  h_line_value <- mlegp_matrix %>%
    filter(no_knots == 600, site_no == site, response == resp, package == "mlegp", pred_type == pred) %>%
    select(RMSqE) %>% as.numeric()
  
  mlegp_time <- mlegp_matrix %>%
    filter(no_knots == 600, site_no == site, response == resp, package == "mlegp", pred_type == pred) %>%
    select(time) %>% as.numeric %>% round(digits = 0)
  
  df <- matrix %>%
    group_by(site_no, response, package, pred_type, no_knots) %>%
    summarise(time = mean(time),
              sd_rmsqe = sd(RMSqE),
              sd_mstde = sd(MStdE),
              RMSqE = mean(RMSqE),
              MStdE = mean(MStdE)) %>% 
    filter(pred_type == pred,
           site_no == site,
           response == resp) %>% as.data.frame()
  
  mod_hetgp <- exp_decay_offset_mod_rmsqe(df %>% filter(package == "hetGP"))
  mod_lagp <- exp_decay_offset_mod_rmsqe(df %>% filter(package == "laGP"))
  mod_mgcv <- exp_decay_offset_mod_rmsqe(df %>% filter(package == "mgcv"))
  
  
  x_grid <- seq(min(df$time),max(df$time),.1)
  pred_hetgp <- predict(mod_hetgp, newdata = data.frame(time = x_grid))
  pred_lagp <- predict(mod_lagp, newdata = data.frame(time = x_grid))
  pred_mgcv <- predict(mod_mgcv, newdata = data.frame(time = x_grid))
  pred_data <- data.frame(x_grid = x_grid, pred_hetgp = pred_hetgp, pred_lagp = pred_lagp)#, pred_mgcv = pred_mgcv)
  
  plot <- df %>%
    ggplot(mapping = aes(x = time, y = RMSqE)) +
    geom_point(mapping = aes(shape = package, color = no_knots), size = 3) +
    geom_line(data = pred_data, mapping = aes(x = x_grid, y = pred_hetgp), alpha = .2) +
    geom_line(data = pred_data, mapping = aes(x = x_grid, y = pred_lagp), alpha = .2) +
   geom_line(data = pred_data, mapping = aes(x = x_grid, y = pred_mgcv), alpha = .2) +
    geom_hline(yintercept = h_line_value, linetype = 'dashed', color = 'red') + ylim(ymin, ymax) +
    geom_errorbar(aes(color = no_knots, ymin = RMSqE - sd_rmsqe, ymax = RMSqE + sd_rmsqe), alpha = 0.3, linetype = 5, size = .5, width = .05) +
    labs(title = "Comparison (prediction accuracy) between packages", subtitle = paste(pred, "predictors used at site",site,"with response",resp),
    caption = paste0("horizontal red line: RMSqE-value with mlegp using 600 knots (", mlegp_time, " seconds)"))
  return(plot)
}

# mean_and_plot_rmsqe(comp_mat, "quantile", 1, 1, 0, 385)

mean_and_plot_mstde <- function(matrix, pred, site, resp, ymin = 0, ymax = NA) {
  h_line_value <- mlegp_matrix %>%
    filter(no_knots == 600, site_no == site, response == resp, package == "mlegp", pred_type == pred) %>%
    select(MStdE) %>% as.numeric()
  
  mlegp_time <- mlegp_matrix %>%
    filter(no_knots == 600, site_no == site, response == resp, package == "mlegp", pred_type == pred) %>%
    select(time) %>% as.numeric %>% round(digits = 0)
  
  df <- matrix %>%
    group_by(site_no, response, package, pred_type, no_knots) %>%
    summarise(time = mean(time),
              sd_rmsqe = sd(RMSqE),
              sd_mstde = sd(MStdE),
              RMSqE = mean(RMSqE),
              MStdE = mean(MStdE)) %>% 
    filter(pred_type == pred,
           site_no == site,
           response == resp) %>% as.data.frame()
  
  mod_hetgp <- exp_decay_offset_mod_mstde(df %>% filter(package == "hetGP"))
  mod_lagp <- exp_decay_offset_mod_mstde(df %>% filter(package == "laGP"))
  #mod_mgcv <- exp_decay_offset_mod_mstde(df %>% filter(package == "mgcv"))
  
  x_grid <- seq(min(df$time),max(df$time),.1)
  pred_hetgp <- predict(mod_hetgp, newdata = data.frame(time = x_grid))
  pred_lagp <- predict(mod_lagp, newdata = data.frame(time = x_grid))
  #pred_mgcv <- predict(mod_mgcv, newdata = data.frame(time = x_grid))
  pred_data <- data.frame(x_grid = x_grid, pred_hetgp = pred_hetgp, pred_lagp = pred_lagp)#, pred_mgcv = pred_mgcv)
  
  plot <- df %>%
    ggplot(mapping = aes(x = time, y = MStdE)) +
    geom_point(mapping = aes(shape = package, color = no_knots), size = 3) +
    geom_line(data = pred_data, mapping = aes(x = x_grid, y = pred_hetgp), alpha = .2) +
    geom_line(data = pred_data, mapping = aes(x = x_grid, y = pred_lagp), alpha = .2) +
    #geom_line(data = pred_data, mapping = aes(x = x_grid, y = pred_mgcv), alpha = .2) +
    geom_hline(yintercept = h_line_value, linetype = 'dashed', color = 'red') + ylim(ymin, ymax) +
    geom_errorbar(aes(color = no_knots, ymin = MStdE - sd_mstde, ymax = MStdE + sd_mstde), alpha = 0.3, linetype = 5, size = .5, width = .05) +
    labs(title = "Comparison (prediction accuracy) between packages", subtitle = paste(pred, "predictors used at site",site,"with response",resp),
         caption = paste0("horizontal red line: RMSqE-value with mlegp using 600 knots (", mlegp_time, " seconds)"))
  return(plot)
}

# mean_and_plot_mstde(comp_mat, "normal", 1, 1, 0, 300)

facet_package_rmsqe <- function(matrix, site, resp) {
  plot <- matrix %>%
    group_by(site_no, response, package, pred_type, no_knots) %>%
    summarise(time = mean(time),
              sd_rmsqe = sd(RMSqE),
              sd_mstde = sd(MStdE),
              RMSqE = mean(RMSqE),
              MStdE = mean(MStdE)) %>% 
    filter(site_no == site_no,
           response == resp) %>%
    ggplot(mapping = aes(x = time, y = RMSqE, color = no_knots, shape = pred_type)) +
    geom_point(size = 2) + facet_wrap(~package, ncol = 1) +
    labs(title = "Differences (prediction accuracy) between predictor types",
         subtitle = paste0("site: ", site, ", response: ", resp))
  return(plot)
}

facet_package_mstde <- function(matrix, site, resp) {
  plot <- matrix %>%
    group_by(site_no, response, package, pred_type, no_knots) %>%
    summarise(time = mean(time),
              sd_rmsqe = sd(RMSqE),
              sd_mstde = sd(MStdE),
              RMSqE = mean(RMSqE),
              MStdE = mean(MStdE)) %>% 
    filter(site_no == site_no,
           response == resp) %>%
    ggplot(mapping = aes(x = time, y = MStdE, color = no_knots, shape = pred_type)) +
    geom_point(size = 2) + facet_wrap(~package, ncol = 1) +
    labs(title = "Differences (prediction certainty) between predictor types",
         subtitle = paste0("site: ", site, ", response: ", resp))
  return(plot)
}


### functions to fit exponential decay (now a bit specific data needs to be given, should work on that)
### assumes that data has columns named RMSqE and time, builds model with these

exp_decay_offset_mod_rmsqe <- function(data) {
  c0 <- 0.5*min(data$RMSqE)
  model0 <- lm(log(RMSqE - c0)~time, data = data)
  a_init <- exp(model0$coefficients[[1]])
  b_init <- model0$coefficients[[2]]
  mod <- nls(RMSqE ~a*exp(b*time) + c, data = data, start = list(a = a_init, b = b_init, c = c0), control = nls.control(maxiter = 100))
  return(mod)
}

exp_decay_offset_mod_mstde <- function(data) {
  c0 <- 0.5*min(data$MStdE)
  model0 <- lm(log(MStdE - c0)~time, data = data)
  a_init <- exp(model0$coefficients[[1]])
  b_init <- model0$coefficients[[2]]
  mod <- nls(MStdE ~a*exp(b*time) + c, data = data, start = list(a = a_init, b = b_init, c = c0), control = nls.control(maxiter = 100))
  return(mod)
}

