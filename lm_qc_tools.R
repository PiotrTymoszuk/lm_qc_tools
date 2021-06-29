# This script provides tools for summaries and QC of linear (lm) and generalized linear models

# libraries ------

  require(tidyverse)
  require(ggrepel)
  require(broom)

# auxiliary functions ----

  point_plot_ <- function(data, x_var, y_var, x_lab = x_var, y_lab = y_var, 
                          plot_title = NULL, smooth = T, silent = T, ...) {
    
    ## draws a simple point plot for diagnostic purposes, takes the output of get_qc_tbl() as data argument
    ## color-codes model missfits
    
    ## table for plotting 
    
    data <- data %>% 
      mutate(misslab = ifelse(.candidate_missfit == 'yes',
                              .rownames, 
                              NA))
    
    ## fill colors
    
    fill_colors <- c(no = 'cornflowerblue', 
                     yes = 'firebrick4')
    
    ## point plot
    
    point_plot <- data %>% 
      ggplot(aes(x = .data[[x_var]], 
                 y = .data[[y_var]], 
                 fill = .candidate_missfit)) + 
      geom_point(size = 2, 
                 shape = 21) + 
      geom_text_repel(aes(label = misslab), 
                      show.legend = F) + 
      scale_fill_manual(values = fill_colors, 
                        name = 'Candidate missfit') + 
      labs(x = x_lab, 
           y = y_lab, 
           title = plot_title)
    
    if(smooth) {
      
      if(silent) {
        
        suppressWarnings(point_plot <- point_plot + 
                           geom_smooth(show.legend = F, 
                                       color = 'black', 
                                       fill = 'dodgerblue2', ...))
        
      } else {
        
        point_plot <- point_plot + 
          geom_smooth(show.legend = F, 
                      color = 'black', 
                      fill = 'dodgerblue2', ...)
        
      }
      
    }
      
    return(point_plot)
    
  }
  
  calc_expected_ <- function(inp_data, observed) {
    
    ## calculates expected normal distribution of a variable observed
    ## credits to: https://stackoverflow.com/questions/43217104/coloring-points-in-a-geom-qq-plot
    
    
    inp_data <- inp_data[order(inp_data[[observed]]), ] %>% 
      mutate(.expect.norm = qnorm(ppoints(nrow(.))))
  
    return(inp_data)
    
  }

# model summary -----

  get_estimates <- function(linear_model, transf_fun = NULL, ...) {
    
    ## extract coefficients from a linear or generalized linear model (betas)
    ## together with confidence intervals. The argument trans_fun allows for 
    ## transformation of the coefficients and CI (e.g. to convert them to OR in logistic regression)
    ## ... specifies other arguments to the CI-calculating function confint()
    
    ## transforming function
    
    if(is.null(transf_fun)) {
      
      transf_fun <- function(x) x
      
    }
    
    ## model summary: to get se and p values
    
    mod_summary <- summary(linear_model)
    
    ## model estimates, error and CI, transforming
    
    model_coefs <- coefficients(linear_model)
    
    model_se <- mod_summary$coefficients[, 2]
    
    model_ci <- confint(linear_model, ...) %>% 
      as_tibble %>% 
      set_names(c('lower_ci', 
                  'upper_ci'))
    
    est_tibble <- model_ci %>% 
      mutate(estimate = unname(model_coefs)) %>% 
      map_dfc(transf_fun) %>% 
      mutate(parameter = names(model_coefs), 
             se = unname(model_se)) %>% 
      select(parameter, 
             estimate, 
             se, 
             lower_ci, 
             upper_ci)
    
    ## p values extracted from model summary
    ## n number of complete observations extracted from the model frame
    
    model_p <- mod_summary$coefficients[, 4]
    
    est_tibble <- est_tibble %>% 
      mutate(p_value = unname(model_p), 
             n_complete = nrow(model.frame(linear_model)))

    return(est_tibble)
    
  }
  
# model QC -----
  
  get_qc_tbl <- function(linear_model, ...) {
    
    ## generates a table with residuals of the model by calling 
    ## augment() from package broom. Finds potential model missfits
    ## by calculating normal distr. p for .std.resid (standardized residuals)
    ## 95% confidence region for the .fitted is estimated with normal distribution
    ## ... are further arguments to augment() function provided by broom
    
    qc_tbl <- augment(linear_model, 
                      se_fit = T, ...)
    
    qc_tbl <- qc_tbl %>% 
      mutate(.sq.std.resid = .std.resid^2, 
             .lower_ci.fit = .fitted + .se.fit * qnorm(0.025), 
             .upper_ci.fit = .fitted + .se.fit * qnorm(0.975), 
             .candidate_missfit = ifelse(abs(.std.resid) > qnorm(0.975), 'yes', 'no'))
    
    ## adding a variable holding expected normal distribution for the standardized residuals
    
    qc_tbl <- calc_expected_(qc_tbl, '.std.resid')
    
    return(qc_tbl)
    
  }
  
  get_qc_plots <- function(linear_model, ...) {
    
    ## draws standard model qc plots with missfits labeled by obs. numbers
    ## ... are additional arguments passed to get_qc_tbl()
    
    ## QC table
    
    qc_tbl <- get_qc_tbl(linear_model, ...)

    ## QC plots
    
    qc_plotting_lst <- list(x_var = c('.fitted', '.fitted', '.fitted', '.expect.norm', '.sigma'), 
                            y_var = c('.resid', '.std.resid', '.sq.std.resid', '.std.resid', '.cooksd'), 
                            plot_title = c('Residuals vs. fitted', 
                                           'Standardized residuals vs. fitted', 
                                           'Sqared residuals vs. fitted', 
                                           'QQ standardized residuals vs expected normal', 
                                           'Cook distance vs dropout sigma'),
                            method = c('loess', 'loess', 'loess', 'lm', 'loess'), 
                            smooth = c(T, T, T, T, F))
    
    qc_plots <- qc_plotting_lst %>% 
      pmap(point_plot_, 
           data = qc_tbl) %>% 
      set_names(c('resid_fitted', 
                  'std.resid_fitted', 
                  'sq.resid_fitted', 
                  'qq.std.resid', 
                  'cook_sigma'))
    

    
    return(qc_plots)
    
  }
  
# END ----