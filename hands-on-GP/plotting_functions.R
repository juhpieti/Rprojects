library(tidyverse)
load("mlegp_matrix.Rdata")

mean_and_plot_rmsqe <- function(matrix, pred, site, resp, ymin = 0, ymax = NA) {
  h_line_value <- mlegp_matrix %>%
    filter(no_knots == 600, site_no == site, response == resp, package == "mlegp", pred_type == pred) %>%
    select(RMSqE) %>% as.numeric()
  
  mlegp_time <- mlegp_matrix %>%
    filter(no_knots == 600, site_no == site, response == resp, package == "mlegp", pred_type == pred) %>%
    select(time) %>% as.numeric %>% round(digits = 0)
  
  plot <- matrix %>%
    group_by(site_no, response, package, pred_type, no_knots) %>%
    summarise(time = mean(time),
              sd_rmsqe = sd(RMSqE),
              sd_mstde = sd(MStdE),
              RMSqE = mean(RMSqE),
              MStdE = mean(MStdE)) %>% 
    filter(pred_type == pred,
           site_no == site,
           response == resp) %>%
    ggplot(mapping = aes(x = time, y = RMSqE, color = no_knots, shape = package)) +
    geom_point(size = 3) +
    geom_hline(yintercept = h_line_value, linetype = 'dashed', color = 'red') + ylim(ymin, ymax) +
    geom_errorbar(aes(color = no_knots, ymin = RMSqE - sd_rmsqe, ymax = RMSqE + sd_rmsqe), alpha = 0.3, linetype = 5, size = .5, width = .05) +
    labs(title = "Comparison (prediction accuracy) between packages", subtitle = paste(pred, "predictors used at site",site,"with response",resp),
    caption = paste0("horizontal red line: RMSqE-value with mlegp using 600 knots (", mlegp_time, " seconds)"))
  return(plot)
}


mean_and_plot_mstde <- function(matrix, pred, site, resp, ymin = 0, ymax = NA) {
  h_line_value <- mlegp_matrix %>%
    filter(no_knots == 600, site_no == site, response == resp, package == "mlegp", pred_type == pred) %>%
    select(MStdE) %>% as.numeric()
  
  mlegp_time <- mlegp_matrix %>%
    filter(no_knots == 600, site_no == site, response == resp, package == "mlegp", pred_type == pred) %>%
    select(time) %>% as.numeric %>% round(digits = 0)
  
  plot <- matrix %>%
    group_by(site_no, response, package, pred_type, no_knots) %>%
    summarise(time = mean(time),
              sd_rmsqe = sd(RMSqE),
              sd_mstde = sd(MStdE),
              RMSqE = mean(RMSqE),
              MStdE = mean(MStdE)) %>% 
    filter(pred_type == pred,
           site_no == site,
           response == resp) %>%
    ggplot(mapping = aes(x = time, y = MStdE)) +
    geom_point(mapping = aes(color = no_knots, shape = package), size = 3) +
    geom_hline(yintercept = h_line_value, linetype = 'dashed', color = 'red') + ylim(ymin, ymax) +
    geom_errorbar(aes(color = no_knots, ymin = MStdE - sd_mstde, ymax = MStdE + sd_mstde), alpha = 0.3, linetype = 5, size = .5, width = .05) +
    labs(title = "Comparison (prediction certainty) between packages", subtitle = paste(pred, "predictors used at site",site,"with response",resp),
         caption = paste0("horizontal red line: MStdE-value with mlegp using 600 knots (", mlegp_time, " seconds)"))
  return(plot)
}

# mean_and_plot_rmsqe(comp_mat, "original", 1, 1)
# mean_and_plot_mstde(comp_mat, "original", 1, 1)

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

mean_and_plot_rmsqe(comp_mat, "original", 1, 1, 0, NA)

exp_decay_offset_mod <- function(data) {
  c0 <- 0.5*min(data$rmsqe)
  model0 <- lm(log(rmsqe - c0)~no_knots, data = data)
  a_init <- exp(model0$coefficients[[1]])
  b_init <- model0$coefficients[[2]]
  mod <- nls(rmsqe ~a*exp(b*no_knots) + c, data = data, start = list(a = a_init, b = b_init, c = c0))
  return(mod)
}