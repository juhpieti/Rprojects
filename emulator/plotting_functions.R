library(tidyverse)
source("helpers.R")
load("mlegp_matrix.Rdata") # to draw horizontal lines of mlegp performance

### these functions assume you have a matrix of experiment results (1st input)
### matrix should include columns for time, RMSqE, MStdE, pred_type, no_knots, site_no (from 1 to 12), response (1 or 2)
### example matrix called "experiments" is given at experiments.Rdata
### they RETURN different kinds of plots made from that matrix of experiments


### for compare_packages_plot() function you give the matrix of experiments, choose parameter space, study site & response,
### as well as the number of observations in the data (750 or 10000 in my case) and number of parameters (5,10,20,40)
### also give the metric to plot ("RMSqE" or "MStdE")
### first it will group your experiments and calculate mean performance
### e.g. mean rmse and time with hetGP using 2000 knots
### then it will give you a plot of chosen metric against time, dots separated by shapes for different packages and by color
### for different amount of knots
### also horizontal line drawn for mlegp performance, caption will give information about that

compare_packages_plot <- function(matrix, site = 1, resp = 1, pred = "original", metric = "RMSqE", n_obs = 750, n_params = 5,
                                  ymin = 0, ymax = NA) {
  
  ### prepare mlegp values for horizontal lines
  mlegp_matrix_mean <- mlegp_matrix %>%
    filter(no_knots == 600, site_no == site, response == resp, package == "mlegp", pred_type == pred,
           obs == n_obs, no_params == n_params) %>%
    group_by(pred_type) %>% 
    summarise(RMSqE = mean(RMSqE),
              MStdE = mean(MStdE),
              time = mean(time)) %>% as.data.frame()
  h_line_value <- mlegp_matrix_mean %>%
    select(ifelse(metric == "RMSqE", 'RMSqE', 'MStdE')) %>% as.numeric()
    
  mlegp_time <- mlegp_matrix_mean %>%
    select(time) %>% as.numeric %>% round(digits = 0)
  
  df <- matrix %>%
    filter(obs == n_obs,
           no_params == n_params) %>%
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
  
  ### uncommenting the rest of the commented lines will also fit you decreasing lines (modeling response value by time)
  ### I commented them because exponential decay model with nls() can cause errors with low number of experiments = weird trends in the data
  
  # mod_hetgp <- exp_decay_offset_mod_rmsqe(df %>% filter(package == "hetGP"))
  # mod_hetgp <- exp_decay_offset_mod(y_het, x_het)
  # mod_lagp <- exp_decay_offset_mod(y_lagp, x_lagp)
  # mod_mgcv <- exp_decay_offset_mod(y_mgcv, x_mgcv)
  # 
  # x_grid <- seq(min(df$time),max(df$time),.1)
  # pred_hetgp <- predict(mod_hetgp, newdata = data.frame(time = x_grid))
  # pred_hetgp <- predict(mod_hetgp, newdata = data.frame(x = x_grid))
  # pred_lagp <- predict(mod_lagp, newdata = data.frame(x = x_grid))
  # pred_mgcv <- predict(mod_mgcv, newdata = data.frame(x = x_grid))
  # pred_data <- data.frame(x_grid = x_grid, pred_hetgp = pred_hetgp, pred_lagp = pred_lagp)#, pred_mgcv = pred_mgcv)

  if (metric == "RMSqE") {
    plot <- df %>%
      ggplot(mapping = aes(x = time, y = RMSqE)) +
      geom_point(mapping = aes(shape = package, color = no_knots), size = 3) +
      #geom_line(data = pred_data, mapping = aes(x = x_grid, y = pred_hetgp), alpha = .2) +
      #geom_line(data = pred_data, mapping = aes(x = x_grid, y = pred_lagp), alpha = .2) +
      #geom_line(data = pred_data, mapping = aes(x = x_grid, y = pred_mgcv), alpha = .2) +
      geom_hline(yintercept = h_line_value, linetype = 'dashed', color = 'red') + ylim(ymin, ymax) +
      geom_errorbar(aes(color = no_knots, ymin = RMSqE - sd_rmsqe, ymax = RMSqE + sd_rmsqe), alpha = 0.3, linetype = 5, size = .5, width = .05) +
      labs(title = "Comparison (prediction accuracy) between packages",
           subtitle = paste0(n_params," ", pred, " predictors used at site ",site," with response ",resp,"  (n = ",n_obs,")"),
           caption = paste0("horizontal red line: RMSqE-value with mlegp using 600 knots (", mlegp_time, " seconds)")) + xlab("time (s)")
  } else {
    plot <- df %>%
      ggplot(mapping = aes(x = time, y = MStdE)) +
      geom_point(mapping = aes(shape = package, color = no_knots), size = 3) +
      #geom_line(data = pred_data, mapping = aes(x = x_grid, y = pred_hetgp), alpha = .2) +
      #geom_line(data = pred_data, mapping = aes(x = x_grid, y = pred_lagp), alpha = .2) +
      #geom_line(data = pred_data, mapping = aes(x = x_grid, y = pred_mgcv), alpha = .2) +
      geom_hline(yintercept = h_line_value, linetype = 'dashed', color = 'red') + ylim(ymin, ymax) +
      geom_errorbar(aes(color = no_knots, ymin = MStdE - sd_mstde, ymax = MStdE + sd_mstde), alpha = 0.3, linetype = 5, size = .5, width = .05) +
      labs(title = "Comparison (prediction certainty) between packages",
           subtitle = paste0(n_params," ", pred, " predictors used at site ",site," with response ",resp,"  (n = ",n_obs,")"),
           caption = paste0("horizontal red line: MStdE-value with mlegp using 600 knots (", mlegp_time, " seconds)")) + xlab("time (s)")
  }
  return(plot)
}


### compare_param_spaces_plot() RETURNS a plot with three (number of packages) rows/subplots, one for each package
### its plotting chosen metric ("RMSqE" or "MStdE") against time, predictor types separeted by shape, number of knots by color
### this is useful to notice how the parameter space used is affecting (e.g original space seems to perform worst wit GPs)

compare_param_spaces_plot <- function(matrix, site = 1, resp = 1, n_obs = 750, n_params = 5, metric = "RMSqE") {
  df <- matrix %>%
    filter(obs == n_obs,
           no_params == n_params) %>%
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
           subtitle = paste0("site: ", site, ", response: ", resp, ", parameters: ", n_params, ", n = ", n_obs))
    return(plot)
  } else {
    plot <- df %>%
      ggplot(mapping = aes(x = time, y = MStdE, color = no_knots, shape = pred_type)) +
      geom_point(size = 2) + facet_wrap(~package, ncol = 1) +
      labs(title = "Differences (prediction certainty) between predictor types by packages",
           subtitle = paste0("site: ", site, ", response: ", resp, ", parameters: ", n_params, ", n = ", n_obs))
    return(plot)
  }
}


### compare_no_params_plot() RETURNS a plot of chosen metric on y-axis and time on x-axis
### it's separating number of parameters used by shape, number of knots by color
### when calling, you need to specify metric, package used, parameter space and number of observations
### other inputs explained:
### no_param_list: list of number of parameters to include in the plot
###                I recommend using at maximum three (e.g c(5,10,20) or c(10,20,40)) to keep the plot clean
### ignore_knots_list: list of number of knots NOT to include in the plot. If runs have been done by a lot of different
###                    number of knots, you might want not to include them all, again to keep the plot kind of clean

compare_no_params_plot <- function(matrix, site = 1, resp = 1, pack = "hetGP", pred = "original", n_obs = 10000,
                                   metric = "RMSqE", no_param_list = c(5,10,20), ignore_knots_list = c(100,300,500,3500,4500),
                                   ymin = 0, ymax = NA) {
  df <- matrix %>%
    filter(obs == n_obs, package == pack, pred_type == pred, 
           no_params %in% no_param_list,
           !(no_knots %in% ignore_knots_list)) %>%
    group_by(site_no, response, package, no_knots, no_params, obs) %>%
    summarise(rmsqe = mean(RMSqE),
              mstde = mean(MStdE),
              time = mean(time),
              n = n()) %>% as.data.frame()
  
  if (metric == "RMSqE") {
    plot <- df %>%
      ggplot(mapping = aes(x = time, y = rmsqe, shape = no_params, color = no_knots)) +
      geom_point(size = 2.5) + ylim(0,ymax) +
      labs(title = "Effect (prediction accuracy) of increasing number of parameters to calibrate",
           subtitle = paste0(pack," with ",pred," predictors (n = ",n_obs, ")")) +
      guides(color = guide_legend(ncol = 2))
  } else {
    plot <- df %>%
      ggplot(mapping = aes(x = time, y = mstde, shape = no_params, color = no_knots)) +
      geom_point(size = 2.5) + ylim(0,ymax) +
      labs(title = "Effect (prediction certainty) of increasing number of parameters to calibrate",
           subtitle = paste0(pack," with ",pred," predictors (n = ",n_obs, ")")) +
      guides(color = guide_legend(ncol = 2))
  }
  return(plot)
}

##### tests / examples #####

# compare_no_params_plot(experiments, 1, 1, "hetGP", "original", no_param_list = c(10,20,40), ymax = NA)
# compare_packages_plot(comp_mat, pred = "original", metric = "RMSqE", n_obs = 750, ymin = 0, ymax = NA)
# compare_param_spaces_plot(experiments, metric = "RMSqE", n_obs = 10000, n_params = 10)

# facet_by_package(comp_mat, metric = "MStdE")

# plot3 <- compare_packages_plot(comp_mat, 1, 1, "normal", "MStdE", 750, 5, ymax = 300)


# compare_no_params_plot(experiments, 1, 1, "hetGP", "original", 10000,"RMSqE",c(10,20,40), c(100,300,500,1500,2500,3500,4500))
# compare_packages_plot(experiments, 1, 1, "quantile", "RMSqE", 10000, 10)
