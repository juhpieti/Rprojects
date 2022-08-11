library(tidyverse)
source("helpers.R")
load("mlegp_matrix.Rdata")
load("comp_mat.Rdata")


### these functions assume you have a matrix of experiment results
### matrix should include columns for time, RMSqE, MStdE, pred_type, no_knots, site_no (from 1 to 12), response (1 or 2)
### example matrix called comp_mat is given at comp_mat.Rdata

### for mean_and_plot function you give the matrix of experiments, choose parameter space, study site & response, 
### as well as the metric to plot (basicly "RMSqE" or "MStdE")
### it will give you a plot of chosen metric against time, dots separated by shapes for different packages and by color
### for different amount of knots
### also horizontal line drawn for mlegp performance

mean_and_plot <- function(matrix, site = 1, resp = 1, pred, metric = "RMSqE", ymin = 0, ymax = NA) {
  
  ### prepare mlegp values for horizontal lines
  mlegp_matrix_mean <- mlegp_matrix %>%
    filter(no_knots == 600, site_no == site, response == resp, package == "mlegp", pred_type == pred) %>%
    group_by(pred_type) %>% 
    summarise(RMSqE = mean(RMSqE),
              MStdE = mean(MStdE),
              time = mean(time)) %>% as.data.frame()
  
  h_line_value <- mlegp_matrix_mean %>%
    select(ifelse(metric == "RMSqE", 'RMSqE', 'MStdE')) %>% as.numeric()
    #select(all_of(metric)) %>% as.numeric()
    
  mlegp_time <- mlegp_matrix_mean %>%
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
  
  df_het <- df %>% filter(package == "hetGP")
  x_het <- df_het$time
  if (metric == "RMSqE") {
    y_het <- df_het$RMSqE
  } else {y_het <- df_het$MStdE}
  
  df_lagp <- df %>% filter(package == "laGP")
  x_lagp <- df_lagp$time
  if (metric == "RMSqE") {
    y_lagp <- df_lagp$RMSqE
  } else {y_lagp <- df_lagp$MStdE}
  
  df_mgcv <- df %>% filter(package == "mgcv")
  x_mgcv <- df_mgcv$time
  if (metric == "RMSqE") {
    y_mgcv <- df_mgcv$RMSqE
  } else {y_mgcv <- df_mgcv$MStdE}
  
  
  #mod_hetgp <- exp_decay_offset_mod_rmsqe(df %>% filter(package == "hetGP"))
  mod_hetgp <- exp_decay_offset_mod(y_het, x_het)
  mod_lagp <- exp_decay_offset_mod(y_lagp, x_lagp)
  #mod_mgcv <- exp_decay_offset_mod(y_mgcv, x_mgcv)
  
  x_grid <- seq(min(df$time),max(df$time),.1)
  # pred_hetgp <- predict(mod_hetgp, newdata = data.frame(time = x_grid))
  pred_hetgp <- predict(mod_hetgp, newdata = data.frame(x = x_grid))
  pred_lagp <- predict(mod_lagp, newdata = data.frame(x = x_grid))
  #pred_mgcv <- predict(mod_mgcv, newdata = data.frame(x = x_grid))
  pred_data <- data.frame(x_grid = x_grid, pred_hetgp = pred_hetgp, pred_lagp = pred_lagp)#, pred_mgcv = pred_mgcv)
  
  if (metric == "RMSqE") {
    plot <- df %>%
      ggplot(mapping = aes(x = time, y = RMSqE)) +
      geom_point(mapping = aes(shape = package, color = no_knots), size = 3) +
      geom_line(data = pred_data, mapping = aes(x = x_grid, y = pred_hetgp), alpha = .2) +
      geom_line(data = pred_data, mapping = aes(x = x_grid, y = pred_lagp), alpha = .2) +
      #geom_line(data = pred_data, mapping = aes(x = x_grid, y = pred_mgcv), alpha = .2) +
      geom_hline(yintercept = h_line_value, linetype = 'dashed', color = 'red') + ylim(ymin, ymax) +
      geom_errorbar(aes(color = no_knots, ymin = RMSqE - sd_rmsqe, ymax = RMSqE + sd_rmsqe), alpha = 0.3, linetype = 5, size = .5, width = .05) +
      labs(title = "Comparison (prediction accuracy) between packages", subtitle = paste(pred, "predictors used at site",site,"with response",resp),
           caption = paste0("horizontal red line: RMSqE-value with mlegp using 600 knots (", mlegp_time, " seconds)")) + xlab("time (s)")
  } else {
    plot <- df %>%
      ggplot(mapping = aes(x = time, y = MStdE)) +
      geom_point(mapping = aes(shape = package, color = no_knots), size = 3) +
      geom_line(data = pred_data, mapping = aes(x = x_grid, y = pred_hetgp), alpha = .2) +
      geom_line(data = pred_data, mapping = aes(x = x_grid, y = pred_lagp), alpha = .2) +
      #geom_line(data = pred_data, mapping = aes(x = x_grid, y = pred_mgcv), alpha = .2) +
      geom_hline(yintercept = h_line_value, linetype = 'dashed', color = 'red') + ylim(ymin, ymax) +
      geom_errorbar(aes(color = no_knots, ymin = MStdE - sd_mstde, ymax = MStdE + sd_mstde), alpha = 0.3, linetype = 5, size = .5, width = .05) +
      labs(title = "Comparison (prediction certainty) between packages", subtitle = paste(pred, "predictors used at site",site,"with response",resp),
           caption = paste0("horizontal red line: MStdE-value with mlegp using 600 knots (", mlegp_time, " seconds)")) + xlab("time (s)")
  }

  return(plot)
}


### facet_by_package is giving you three (number of packages) plots, one for each package
### its plotting chosen metric ("RMSqE" or "MStdE") against time, predictor types separeted by shape, No. knots by color
### this is useful to notice how the parameter space used is affecting (e.g original space seems to perform worst wit GPs)

facet_by_package <- function(matrix, site = 1, resp = 1, metric = "RMSqE") {
  df <- matrix %>%
    group_by(site_no, response, package, pred_type, no_knots) %>%
    summarise(time = mean(time),
              sd_rmsqe = sd(RMSqE),
              sd_mstde = sd(MStdE),
              RMSqE = mean(RMSqE),
              MStdE = mean(MStdE)) %>% 
    filter(site_no == site_no,
           response == resp) %>% as.data.frame()
  
  if (metric == "RMSqE") {
    plot <- df %>%
      ggplot(mapping = aes(x = time, y = RMSqE, color = no_knots, shape = pred_type)) +
      geom_point(size = 2) + facet_wrap(~package, ncol = 1) +
      labs(title = "Differences (prediction accuracy) between predictor types by packages",
           subtitle = paste0("site: ", site, ", response: ", resp))
    return(plot)
  } else {
    plot <- df %>%
      ggplot(mapping = aes(x = time, y = MStdE, color = no_knots, shape = pred_type)) +
      geom_point(size = 2) + facet_wrap(~package, ncol = 1) +
      labs(title = "Differences (prediction certainty) between predictor types by packages",
           subtitle = paste0("site: ", site, ", response: ", resp))
    return(plot)
  }
}


##### tests #####

mean_and_plot(comp_mat, pred = "normal", metric = "RMSqE", ymin = 0, ymax = 400)
# facet_by_package(comp_mat, metric = "RMSqE")
# facet_by_package(comp_mat, metric = "MStdE")


