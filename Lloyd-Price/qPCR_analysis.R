#!/usr/bin/Rscript
library(readxl)
library(ggplot2)
library(ggpubr)
library(FSA)
library(tidyverse)
library(rstatix)

library(drc) # For 4PL dose-response curve
library(scales)

# Prevent scientific notation in plots
options(scipen = 999)

# Try to load Cairo for better PDF rendering, but continue if not available
tryCatch({
  library(Cairo)
  cairo_available <- TRUE
}, error = function(e) {
  message("Cairo package not available, using standard PDF device")
  cairo_available <- FALSE
})

# Set graphics device options to avoid font issues
options(device = function(...) {
  grDevices::png(..., type = "cairo")
})

# Try to set a safe font family
tryCatch({
  par(family = "sans")
}, error = function(e) {
  message("Could not set font family, using default")
})

library(dplyr)
library(stringr)

# Publication-ready theme
publication_theme <- function(base_size = 14, base_family = "sans") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      # Title and labels
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 10)),
      plot.subtitle = element_text(size = 12, hjust = 0.5, margin = margin(b = 10)),
      axis.title = element_text(size = 14, face = "bold"),
      axis.title.x = element_text(margin = margin(t = 10)),
      axis.title.y = element_text(margin = margin(r = 10)),

      # Axis text
      axis.text = element_text(size = 12, color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),

      # Legend
      legend.position = "top",
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 11),
      legend.background = element_rect(fill = "white", color = "grey80"),
      legend.key = element_rect(fill = "white"),

      # Panel and background
      panel.background = element_rect(fill = "white", color = NA),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
      panel.grid.minor = element_line(color = "grey95", linewidth = 0.1),

      # Plot background
      plot.background = element_rect(fill = "white", color = NA),

      # Margins
      plot.margin = margin(t = 20, r = 20, b = 20, l = 20),

      # Facet labels
      strip.text = element_text(size = 12, face = "bold"),
      strip.background = element_rect(fill = "grey90", color = NA)
    )
}

# Professional color palette
publication_colors <- c(
  "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
  "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"
)


# Load the data
setwd('/home/borisvdm/Documents/PhD/gut_brain/IBD/Lloyd-Price2019/Exp_validation')
# Read Excel with

file <- 'Analyse_qPCR.xlsx'

# Create output directories if they don't exist
dir.create('./Boxplots', showWarnings = FALSE, recursive = TRUE)
dir.create('./Dotplots', showWarnings = FALSE, recursive = TRUE)
dir.create('./Dose-response', showWarnings = FALSE, recursive = TRUE)
dir.create('./linear_model', showWarnings = FALSE, recursive = TRUE)
dir.create('./Statistical_summaries', showWarnings = FALSE, recursive = TRUE)

# Initialize statistical summary files
stats_summary_file <- './Statistical_summaries/dose_response_statistics.txt'
cat("DOSE-RESPONSE CURVE STATISTICAL SUMMARY\n", 
    "Generated on:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n",
    strrep("=", 80), "\n\n", file = stats_summary_file)

linear_stats_file <- './Statistical_summaries/linear_model_statistics.txt'
cat("LINEAR MODEL STATISTICAL SUMMARY\n", 
    "Generated on:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n",
    strrep("=", 80), "\n\n", file = linear_stats_file)

# PDF saving configuration - set to FALSE to disable PDF saving temporarily
enable_pdf_saving <- FALSE

# Statistical annotations configuration - set to FALSE to disable temporarily
enable_stat_annotations <- TRUE

# Inform user about configurations
if (cairo_available && enable_pdf_saving) {
  message("✓ Cairo PDF support available - high-quality PDFs will be generated")
} else if (!enable_pdf_saving) {
  message("ℹ PDF saving disabled - generating PNG files only")
  message("💡 To re-enable PDF saving, change 'enable_pdf_saving <- FALSE' to 'enable_pdf_saving <- TRUE'")
} else {
  message("ℹ Cairo PDF support not available - using standard PDF device")
}

if (!enable_stat_annotations) {
  message("ℹ Statistical annotations temporarily disabled to avoid font issues")
  message("💡 To re-enable statistical annotations, change 'enable_stat_annotations <- FALSE' to 'enable_stat_annotations <- TRUE'")
}

tabs <- excel_sheets(file)

get_significance <- function(p_val) {
  if (p_val < 0.001) {
    return("***")
  } else if (p_val < 0.01) {
    return("**")
  } else if (p_val < 0.05) {
    return("*")
  } else {
    return("ns")
  }
}

for (tab in tabs[2:28]) {
  # read 1 tab
  #tab <- tabs[2]
  dat <- read_excel(file, sheet = tab)
  print(paste0('Analyzing data for ', tab))
  # row 2 are column names
  #dat <- dat[-1,]
  colnames(dat) <- dat[1,]; dat <- dat[-1,]
  # replace ' ' by _' in colnames
  colnames(dat) <- gsub(' ', '_', colnames(dat))
  # set all columns excep 1:3 to numeric
  dat[,4:ncol(dat)] <- sapply(dat[,4:ncol(dat)], as.numeric)
  # Reverse the order of Group levels for all plots
  dat$Group <- factor(dat$Group, levels = rev(unique(dat$Group)))
  
  # set reference genes: HPRT and HMBS if group is 'aGPC', GAPDH and HMBS if group contains 'Car' or 'Trig'
  if (grepl('aGPC', dat$Group[1])) {
    treatment_type <- 'aGPC'
    print(treatment_type)
    # set reference gene columns
    ref1 <- 'HPRT'; ref2 <- 'HMBS'
  } else if (grepl('Car', dat$Group[1])) {
    #print('Car or Trig')
    treatment_type <- 'Car'
    print(treatment_type)
    # set reference gene columns
    ref1 <- 'GAPDH'; ref2 <- 'HMBS'
  } else if (grepl('Trig', dat$Group[1])) {
    #print('Car or Trig')
    treatment_type <- 'Trig'
    print(treatment_type)
    # set reference gene columns
    ref1 <- 'GAPDH'; ref2 <- 'HMBS'
  } 
  
  
  # in column geom_mean, calculate the geometric mean of the reference genes
  dat$Geom_mean <- apply(dat[,c(ref1, ref2)], 1, function(x) prod(as.numeric(x))^(1/length(x)))
  Cts <- dat$Ct_gene[1:48]
  # calculate the average of Ct1 and Ct2, Ct3 and Ct4, Ct5 and Ct6
  averages <- numeric()
  # Loop with step size 2
  for (i in seq(1, length(Cts) - 1, by = 2)) {
    avg <- mean(Cts[i:(i+1)])  # Calculate the average of the current and next value
    averages <- c(averages, avg)  # Store the result
  }
  # Drop rows 25:48
  dat <- dat[1:24,]
  dat$Average_Ct <- averages
  # Drop samples that we don't want to include (... 0 4 was pos control for LDH assay, ... 100uM 1 also seems off in comparison to similar samples)
  # That's row 24 and 1
  dat <- dat[-c(24, 1),]
  # calculate delta_Ct
  dat$delta_Ct <- abs(dat$Geom_mean - dat$Average_Ct)
  delta_Ct_control_avg <- mean(dat$delta_Ct[dat$Group == paste0(treatment_type, ' 0')])
  dat$delta_delta_Ct <- dat$delta_Ct - delta_Ct_control_avg
  dat$Fold_gene_expression <- 2^(-dat$delta_delta_Ct)
  
  
  # Test if the data is normally distributed
  # Shapiro-Wilk normality test
  # Perform t-tests comparing each condition to baseline (0 dose group)
  baseline_group <- paste0(treatment_type, ' 0')
  print(paste('Performing t-tests comparing each condition to baseline:', baseline_group))
  
  # Get baseline data
  baseline_data <- dat$Fold_gene_expression[dat$Group == baseline_group]
  
  # Perform t-tests for each group vs baseline
  t_test_results <- dat %>%
    group_by(Group) %>%
    filter(Group != baseline_group) %>%  # Exclude baseline from comparison
    summarise(
      n = n(),
      mean_expr = mean(Fold_gene_expression, na.rm = TRUE),
      sd_expr = sd(Fold_gene_expression, na.rm = TRUE),
      p_value = if(n() > 0) {
        group_data <- Fold_gene_expression
        if(length(baseline_data) > 0 && length(group_data) > 0) {
          t.test(group_data, baseline_data)$p.value
        } else NA
      } else NA,
      .groups = 'drop'
    ) %>%
    mutate(
      p.adj = p.adjust(p_value, method = "BH"),
      significance = case_when(
        p.adj < 0.001 ~ "***",
        p.adj < 0.01 ~ "**", 
        p.adj < 0.1 ~ "*",
        TRUE ~ "ns"
      )
    ) %>%
    # Add baseline group for visualization
    bind_rows(
      dat %>%
        filter(Group == baseline_group) %>%
        summarise(
          Group = baseline_group,
          n = n(),
          mean_expr = mean(Fold_gene_expression, na.rm = TRUE),
          sd_expr = sd(Fold_gene_expression, na.rm = TRUE),
          p_value = NA,
          p.adj = NA,
          significance = "baseline",
          .groups = 'drop'
        )
    )
  
  # Print results
  print("T-test results vs baseline:")
  print(t_test_results)
  
  # Prepare data for stat_pvalue_manual (only significant comparisons)
  sig_comparisons <- t_test_results %>%
    filter(significance %in% c("*", "**", "***")) %>%
    mutate(
      group1 = baseline_group,
      group2 = Group
    )
  
  # Calculate y.position for each significant comparison
  if (nrow(sig_comparisons) > 0) {
    max_y <- max(dat$Fold_gene_expression, na.rm = TRUE)
    min_y <- min(dat$Fold_gene_expression, na.rm = TRUE)
    y_range <- max_y - min_y
    
    # Position bars starting at max_y + 8% of range, with 6% spacing between bars for text
    sig_comparisons <- sig_comparisons %>%
      mutate(
        y.position = max_y + y_range * 0.08 + (row_number() - 1) * y_range * 0.06
      )
    
    # Calculate the new upper limit for the plot to accommodate bars and text
    max_bar_y <- max(sig_comparisons$y.position) + y_range * 0.05
  } else {
    max_bar_y <- max(dat$Fold_gene_expression, na.rm = TRUE)
  }

  # Create a summary for subtitle
  n_comparisons <- nrow(t_test_results) - 1  # Exclude baseline
  n_significant <- sum(t_test_results$significance %in% c("*", "**", "***"), na.rm = TRUE)
  test_summary <- paste0("T-tests vs baseline (", baseline_group, "): ", 
                        n_significant, "/", n_comparisons, " significant (p.adj < 0.1)")

  # Debug: Print sig_comparisons structure
  if (nrow(sig_comparisons) > 0) {
    print("Debug - sig_comparisons structure:")
    print(sig_comparisons)
    print("Debug - column names:")
    print(colnames(sig_comparisons))
  } else {
    print("Debug - No significant comparisons found")
  }

  # Plot the boxplot with significance annotations
  boxplot <- ggboxplot(dat,
            x = "Group",
            y = "Fold_gene_expression",
            fill = "Group",
            palette = publication_colors,
            add = "jitter",
            add.params = list(size = 1.5, alpha = 0.7))

  # Add manual significance annotations if enabled and there are significant comparisons
  if (enable_stat_annotations && nrow(sig_comparisons) > 0) {
    print(paste("Debug - Adding", nrow(sig_comparisons), "significance bars manually"))
    
    # Add significance bars and labels manually
    for (i in 1:nrow(sig_comparisons)) {
      comparison <- sig_comparisons[i, ]
      
      # Get x positions for baseline and treatment groups using factor levels
      all_groups <- levels(dat$Group)  # Use factor levels to get correct order
      baseline_x <- which(all_groups == comparison$group1)
      treatment_x <- which(all_groups == comparison$group2)
      
      # Debug: print group positions
      print(paste("Debug - Groups:", paste(all_groups, collapse = ", ")))
      print(paste("Debug - Connecting", comparison$group1, "(position", baseline_x, ") to", comparison$group2, "(position", treatment_x, ")"))
      
      # Add horizontal line (significance bar)
      boxplot <- boxplot + 
        annotate("segment", 
                x = baseline_x, xend = treatment_x,
                y = comparison$y.position, yend = comparison$y.position,
                linewidth = 0.5, color = "black") +
        # Add vertical connectors
        annotate("segment", 
                x = baseline_x, xend = baseline_x,
                y = comparison$y.position - max(dat$Fold_gene_expression, na.rm = TRUE) * 0.01, 
                yend = comparison$y.position,
                linewidth = 0.5, color = "black") +
        annotate("segment", 
                x = treatment_x, xend = treatment_x,
                y = comparison$y.position - max(dat$Fold_gene_expression, na.rm = TRUE) * 0.01, 
                yend = comparison$y.position,
                linewidth = 0.5, color = "black") +
        # Add significance symbol above the bar
        annotate("text", 
                x = (baseline_x + treatment_x) / 2,  # Center the text between the two groups
                y = comparison$y.position + max(dat$Fold_gene_expression, na.rm = TRUE) * 0.01,  # Closer to the bar
                label = comparison$significance, 
                size = 4, color = "black", fontface = "bold")
      
      # Debug: print positioning info
      print(paste("Debug - Bar", i, "for", comparison$group2, "vs", comparison$group1))
    }
  } else {
    print("Debug - No significance bars to add")
  }

  boxplot <- boxplot +
    labs(subtitle = test_summary,
         title = paste0("Gene Expression: ", tab),
         x = "Treatment Group",
         y = "Relative expression") +
    scale_y_continuous(labels = function(x) format(x, scientific = FALSE),
                      limits = c(NA, max_bar_y)) +  # Expand upper limit for significance bars
    publication_theme() +
    theme(legend.position = "none")  # Remove legend since colors are redundant with x-axis

  print(boxplot)
  ggsave(paste0('./Boxplots/', tab, '_boxplot.png'), plot = boxplot,
         width = 8, height = 6, units = 'in', dpi = 600, bg = "white")

  # Save PDF only if enabled
  if (enable_pdf_saving) {
    tryCatch({
      if (cairo_available) {
        ggsave(paste0('./Boxplots/', tab, '_boxplot.pdf'), plot = boxplot,
               width = 8, height = 6, units = 'in', dpi = 600, bg = "white",
               device = cairo_pdf)
      } else {
        ggsave(paste0('./Boxplots/', tab, '_boxplot.pdf'), plot = boxplot,
               width = 8, height = 6, units = 'in', dpi = 600, bg = "white")
      }
    }, error = function(e) {
      message("PDF saving failed for boxplot, continuing with PNG only: ", e$message)
    })
  }

  # Enhanced dot plot with error bars
  dotplot <- ggplot(dat, aes(x = Group, y = Fold_gene_expression, color = Group)) +
    geom_jitter(width = 0.2, size = 2, alpha = 0.8) +
    stat_summary(fun = mean, geom = "point", shape = 18, size = 4, color = "black", stroke = 1.5) +
    stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1),
                 geom = "errorbar", width = 0.15, color = "black", linewidth = 1)

  # Add manual significance annotations if enabled and there are significant comparisons
  if (enable_stat_annotations && nrow(sig_comparisons) > 0) {
    print(paste("Debug - Adding", nrow(sig_comparisons), "significance bars to dotplot manually"))
    
    # Add significance bars and labels manually
    for (i in 1:nrow(sig_comparisons)) {
      comparison <- sig_comparisons[i, ]
      
      # Get x positions for baseline and treatment groups using factor levels
      all_groups <- levels(dat$Group)  # Use factor levels to get correct order
      baseline_x <- which(all_groups == comparison$group1)
      treatment_x <- which(all_groups == comparison$group2)
      
      # Add horizontal line (significance bar)
      dotplot <- dotplot + 
        annotate("segment", 
                x = baseline_x, xend = treatment_x,
                y = comparison$y.position, yend = comparison$y.position,
                linewidth = 0.5, color = "black") +
        # Add vertical connectors
        annotate("segment", 
                x = baseline_x, xend = baseline_x,
                y = comparison$y.position - max(dat$Fold_gene_expression, na.rm = TRUE) * 0.01, 
                yend = comparison$y.position,
                linewidth = 0.5, color = "black") +
        annotate("segment", 
                x = treatment_x, xend = treatment_x,
                y = comparison$y.position - max(dat$Fold_gene_expression, na.rm = TRUE) * 0.01, 
                yend = comparison$y.position,
                linewidth = 0.5, color = "black") +
        # Add significance symbol above the bar
        annotate("text", 
                x = (baseline_x + treatment_x) / 2,  # Center the text between the two groups
                y = comparison$y.position + max(dat$Fold_gene_expression, na.rm = TRUE) * 0.01,  # Closer to the bar
                label = comparison$significance, 
                size = 4, color = "black", fontface = "bold")
    }
  }

  dotplot <- dotplot +
    labs(subtitle = test_summary,
         title = paste0("Gene Expression: ", tab),
         x = "Treatment Group",
         y = "Relative expression") +
    scale_y_continuous(labels = function(x) format(x, scientific = FALSE),
                      limits = c(NA, max_bar_y)) +  # Expand upper limit for significance bars
    scale_color_manual(values = publication_colors) +
    publication_theme() +
    theme(legend.position = "none")  # Remove legend since colors are redundant

  print(dotplot)
  ggsave(paste0('./Dotplots/', tab, '_dotplot.png'), plot = dotplot,
         width = 8, height = 6, units = 'in', dpi = 600, bg = "white")

  # Save PDF only if enabled
  if (enable_pdf_saving) {
    tryCatch({
      if (cairo_available) {
        ggsave(paste0('./Dotplots/', tab, '_dotplot.pdf'), plot = dotplot,
               width = 8, height = 6, units = 'in', dpi = 600, bg = "white",
               device = cairo_pdf)
      } else {
        ggsave(paste0('./Dotplots/', tab, '_dotplot.pdf'), plot = dotplot,
               width = 8, height = 6, units = 'in', dpi = 600, bg = "white")
      }
    }, error = function(e) {
      message("PDF saving failed for dotplot, continuing with PNG only: ", e$message)
    })
  }
  
  # Dose-response only if meaningful dose gradient is present
  if (grepl('Trig|Car|aGPC', treatment_type)) {
    print('Attempting to fit dose-response curve with ggplot2...')
    
    
    # --- Data Preparation ---
    dat <- dat %>%
      mutate(
        # Extract dose and units - improved regex to capture 0
        Dose_raw = str_extract(Group, "\\d+\\.?\\d*|\\b0\\b"),
        Unit = str_extract(Group, "(nM|uM)"),
        
        # Convert to nM
        Dose_nM = case_when(
          Unit == "uM" ~ as.numeric(Dose_raw) * 1000,
          Unit == "nM" ~ as.numeric(Dose_raw),
          # Handle the 0 dose case (control group)
          str_detect(Group, "\\b0\\b") ~ 0,
          TRUE ~ NA_real_
        ),
        
        # Convert to uM (for more intuitive log scale)
        Dose_uM = Dose_nM / 1000,
        
        # Add pseudovalue for log(0)
        Dose_uM_pseudo = ifelse(Dose_uM == 0, 0.001, Dose_uM),
        
        # Log10-transformed dose for plotting
        Dose_log10 = log10(Dose_uM_pseudo)
      ) %>%
      filter(!is.na(Dose_nM), !is.na(Fold_gene_expression))  # Remove NAs
    
    # --- Model Fitting ---
    drc_model <- tryCatch({
      drm(Fold_gene_expression ~ Dose_nM,
          data = dat,
          fct = LL.4())
    }, error = function(e) {
      warning(paste("Dose-response model failed for", tab))
      return(NULL)
    })
    
    # --- Plotting ---
    if (!is.null(drc_model)) {
      
      # Prediction data - include the baseline (0 dose) by using the pseudovalue
      min_dose_for_seq <- min(dat$Dose_nM[dat$Dose_nM > 0]) * 0.1  # Start slightly below lowest non-zero dose
      dose_seq <- data.frame(Dose_nM = c(0.001, exp(seq(log(min_dose_for_seq),
                                                        log(max(dat$Dose_nM)),
                                                        length.out = 199))))
      predictions <- predict(drc_model, newdata = dose_seq, interval = "confidence")
      dose_seq$fit <- predictions[,1]
      dose_seq$lower <- predictions[,2]
      dose_seq$upper <- predictions[,3]
      
      # Convert prediction x-axis to log10(uM) scale for plotting
      dose_seq$Dose_log10 <- log10(dose_seq$Dose_nM / 1000)
      
      # Extract model summary and p-value more robustly
      model_summary <- summary(drc_model)

      # Calculate multiple statistical measures
      slope_pval <- tryCatch({
        # Method 1: Use model comparison with null model (no dose effect)
        null_model <- lm(Fold_gene_expression ~ 1, data = dat)  # Null model: just intercept
        dose_model <- lm(Fold_gene_expression ~ Dose_nM, data = dat)  # Simple linear model

        # Compare models using ANOVA
        model_comparison <- anova(null_model, dose_model)
        model_comparison[2, "Pr(>F)"]  # p-value from F-test
      }, error = function(e) {
        print(paste("Error in model comparison:", e$message))
        NA
      })

      # Alternative: Use correlation test between dose and response
      if (is.na(slope_pval)) {
        slope_pval <- tryCatch({
          cor_test <- cor.test(dat$Dose_nM, dat$Fold_gene_expression, method = "spearman")
          cor_test$p.value
        }, error = function(e) {
          print(paste("Error in correlation test:", e$message))
          NA
        })
      }
      
      # Calculate R-squared for goodness of fit
      r_squared <- tryCatch({
        # Calculate R-squared for DRC model
        observed <- dat$Fold_gene_expression
        predicted <- predict(drc_model)
        ss_res <- sum((observed - predicted)^2)
        ss_tot <- sum((observed - mean(observed))^2)
        1 - (ss_res / ss_tot)
      }, error = function(e) {
        print(paste("Error calculating R-squared:", e$message))
        NA
      })
      
      # Calculate AIC for model quality
      model_aic <- tryCatch({
        AIC(drc_model)
      }, error = function(e) {
        print(paste("Error calculating AIC:", e$message))
        NA
      })

      # Calculate EC50
      ec50 <- tryCatch({
        ED(drc_model, 50, type = "absolute", display = FALSE)[1]
      }, error = function(e) {
        NA
      })

      # Format comprehensive statistical summary
      stat_lines <- c()
      
      # P-value
      if (!is.na(slope_pval)) {
        if (slope_pval < 0.001) {
          stat_lines <- c(stat_lines, "p < 0.001 ***")
        } else if (slope_pval < 0.01) {
          stat_lines <- c(stat_lines, sprintf("p = %.3f **", slope_pval))
        } else if (slope_pval < 0.05) {
          stat_lines <- c(stat_lines, sprintf("p = %.3f *", slope_pval))
        } else {
          stat_lines <- c(stat_lines, sprintf("p = %.3f (ns)", slope_pval))
        }
      } else {
        stat_lines <- c(stat_lines, "p-value: N/A")
      }
      
      # R-squared
      if (!is.na(r_squared)) {
        stat_lines <- c(stat_lines, sprintf("R² = %.3f", r_squared))
      }
      
      # EC50
      if (!is.na(ec50)) {
        ec50_uM <- ec50 / 1000
        if (ec50_uM < 1) {
          stat_lines <- c(stat_lines, sprintf("EC50 = %.1f nM", ec50))
        } else {
          stat_lines <- c(stat_lines, sprintf("EC50 = %.2f µM", ec50_uM))
        }
      }
      
      # Model quality
      if (!is.na(model_aic)) {
        stat_lines <- c(stat_lines, sprintf("AIC = %.1f", model_aic))
      }
      
      # Combine all statistics
      stat_label <- paste(stat_lines, collapse = "\n")
      
      # Add sample size information
      n_points <- nrow(dat)
      stat_label <- paste0("n = ", n_points, "\n", stat_label)
      
      p <- ggplot(dat, aes(x = Dose_log10, y = Fold_gene_expression)) +
        # Data points with error bars
        geom_point(size = 3, alpha = 0.8, color = publication_colors[1], stroke = 0.5) +
        # Confidence interval ribbon
        geom_ribbon(data = dose_seq,
                   aes(x = Dose_log10, ymin = lower, ymax = upper),
                   alpha = 0.2, fill = publication_colors[2],
                   inherit.aes = FALSE) +
        # Fitted curve
        geom_line(data = dose_seq, aes(x = Dose_log10, y = fit),
                 color = publication_colors[2], linewidth = 1.5,
                 inherit.aes = FALSE) +
        # Formatting - dynamic x-axis based on actual data range
        scale_x_continuous(
          breaks = log10(c(0.001, 0.01, 0.1, 1, 10, 100)),
          labels = c("0", "10 nM", "100 nM", "1 µM", "10 µM", "100 µM"),
          limits = c(min(dat$Dose_log10) - 0.1, max(dat$Dose_log10) + 0.1)
        ) +
        scale_y_continuous(labels = function(x) format(x, scientific = FALSE)) +
        labs(
          title = paste0("Dose-Response Curve: ", tab),
          subtitle = paste0("Treatment: ", treatment_type, " | ", 
                           if(!is.na(slope_pval)) paste0("p = ", sprintf("%.3f", slope_pval), 
                           if(slope_pval < 0.001) " ***" else if(slope_pval < 0.01) " **" else if(slope_pval < 0.05) " *" else " (ns)") else "",
                           if(!is.na(r_squared)) paste0(" | R² = ", sprintf("%.3f", r_squared)) else "",
                           if(!is.na(ec50)) paste0(" | EC50 = ", sprintf("%.1f", if(ec50/1000 < 1) ec50 else ec50/1000), if(ec50/1000 < 1) " nM" else " µM") else ""),
          x = "Concentration",
          y = "Relative expression"
        ) +
        publication_theme()      # Print & save with high quality
      print(p)
      ggsave(paste0('./Dose-response/', tab, '_dose_response.png'), plot = p,
             width = 8, height = 6, units = 'in', dpi = 600, bg = "white")

      # Save PDF only if enabled
      if (enable_pdf_saving) {
        tryCatch({
          if (cairo_available) {
            ggsave(paste0('./Dose-response/', tab, '_dose_response.pdf'), plot = p,
                   width = 8, height = 6, units = 'in', dpi = 600, bg = "white",
                   device = cairo_pdf)
          } else {
            ggsave(paste0('./Dose-response/', tab, '_dose_response.pdf'), plot = p,
                   width = 8, height = 6, units = 'in', dpi = 600, bg = "white")
          }
        }, error = function(e) {
          message("PDF saving failed for dose-response curve, continuing with PNG only: ", e$message)
        })
      }
      
      # Print comprehensive model summary and statistics
      cat("\n", strrep("=", 60), "\n")
      cat("DOSE-RESPONSE ANALYSIS RESULTS FOR:", tab, "\n")
      cat(strrep("=", 60), "\n")
      cat("Treatment type:", treatment_type, "\n")
      cat("Sample size:", n_points, "data points\n")
      cat("Model: 4-Parameter Logistic (LL.4)\n\n")
      
      cat("STATISTICAL MEASURES:\n")
      if (!is.na(slope_pval)) {
        cat("• P-value (dose effect):", sprintf("%.4f", slope_pval))
        if (slope_pval < 0.001) cat(" (***)")
        else if (slope_pval < 0.01) cat(" (**)")
        else if (slope_pval < 0.05) cat(" (*)")
        else cat(" (ns)")
        cat("\n")
      }
      if (!is.na(r_squared)) cat("• R-squared:", sprintf("%.4f", r_squared), "\n")
      if (!is.na(ec50)) {
        ec50_uM <- ec50 / 1000
        cat("• EC50:", sprintf("%.2f µM (%.1f nM)", ec50_uM, ec50), "\n")
      }
      if (!is.na(model_aic)) cat("• AIC:", sprintf("%.2f", model_aic), "\n")
      
      cat("\nDOSE RANGE:\n")
      cat("• Minimum dose:", min(dat$Dose_uM[dat$Dose_uM > 0]), "µM\n")
      cat("• Maximum dose:", max(dat$Dose_uM), "µM\n")
      cat("• Dose levels:", length(unique(dat$Dose_uM)), "\n")
      
      cat("\nRESPONSE RANGE:\n")
      cat("• Minimum response:", sprintf("%.3f", min(dat$Fold_gene_expression)), "\n")
      cat("• Maximum response:", sprintf("%.3f", max(dat$Fold_gene_expression)), "\n")
      cat("• Response range:", sprintf("%.3f", max(dat$Fold_gene_expression) - min(dat$Fold_gene_expression)), "\n")
      cat(strrep("=", 60), "\n")
      
      # Save the same information to the summary file
      cat("\nGENE:", tab, "\n",
          "Treatment:", treatment_type, "\n",
          "Sample size:", n_points, "\n",
          "Model: 4-Parameter Logistic\n",
          file = stats_summary_file, append = TRUE)
      
      if (!is.na(slope_pval)) {
        significance <- if (slope_pval < 0.001) "***" else if (slope_pval < 0.01) "**" else if (slope_pval < 0.05) "*" else "ns"
        cat("P-value:", sprintf("%.4f", slope_pval), significance, "\n",
            file = stats_summary_file, append = TRUE)
      }
      if (!is.na(r_squared)) cat("R-squared:", sprintf("%.4f", r_squared), "\n", 
                                 file = stats_summary_file, append = TRUE)
      if (!is.na(ec50)) {
        ec50_uM <- ec50 / 1000
        cat("EC50:", sprintf("%.2f µM (%.1f nM)", ec50_uM, ec50), "\n",
            file = stats_summary_file, append = TRUE)
      }
      if (!is.na(model_aic)) cat("AIC:", sprintf("%.2f", model_aic), "\n", 
                                 file = stats_summary_file, append = TRUE)
      
      cat("Dose range:", min(dat$Dose_uM[dat$Dose_uM > 0]), "-", max(dat$Dose_uM), "µM\n",
          "Response range:", sprintf("%.3f", min(dat$Fold_gene_expression)), "-", 
          sprintf("%.3f", max(dat$Fold_gene_expression)), "\n",
          strrep("=", 50), "\n\n", file = stats_summary_file, append = TRUE)
      
      # Print model summary for debugging
      print(paste("Detailed model summary for", tab, ":"))
      print(model_summary)
      
    } else {
      warning(paste("Could not fit dose-response model for", tab))
    }
  }
  
  # ============================================================================
  # LINEAR MODEL ANALYSIS
  # ============================================================================
  
  # Linear model analysis (for dose-response data)
  if (grepl('Trig|Car|aGPC', treatment_type)) {
    print(paste('Fitting linear model for', tab))
    
    # Prepare data for linear model (log-transform dose)
    dat_linear <- dat %>%
      filter(Dose_uM > 0) %>%  # Exclude control (dose = 0) for log transformation
      mutate(
        log_dose = log10(Dose_uM),
        log_dose_centered = log_dose - mean(log_dose, na.rm = TRUE)
      )
    
    # Fit linear model
    tryCatch({
      # Simple linear regression: expression ~ log(dose)
      linear_model <- lm(Fold_gene_expression ~ log_dose, data = dat_linear)
      model_summary <- summary(linear_model)
      
      # Extract statistics
      slope <- coef(linear_model)[2]
      slope_pval <- model_summary$coefficients[2, 4]
      r_squared <- model_summary$r.squared
      adj_r_squared <- model_summary$adj.r.squared
      intercept <- coef(linear_model)[1]
      n_points <- nrow(dat_linear)
      
      # Calculate confidence intervals
      conf_int <- confint(linear_model)
      slope_ci_lower <- conf_int[2, 1]
      slope_ci_upper <- conf_int[2, 2]
      
      # Create prediction data for plotting
      dose_range <- seq(min(dat_linear$Dose_uM), max(dat_linear$Dose_uM), length.out = 100)
      log_dose_range <- log10(dose_range)
      pred_data <- data.frame(log_dose = log_dose_range)
      predictions <- predict(linear_model, pred_data, interval = "confidence")
      plot_data <- data.frame(
        Dose_uM = dose_range,
        fit = predictions[, 1],
        lwr = predictions[, 2],
        upr = predictions[, 3]
      )
      
      # Create the linear model plot
      linear_plot <- ggplot() +
        # Data points (non-zero doses only)
        geom_point(data = dat_linear, aes(x = Dose_uM, y = Fold_gene_expression), 
                   size = 3, alpha = 0.7, color = "#2E86AB") +
        # Fitted line
        geom_line(data = plot_data, aes(x = Dose_uM, y = fit), 
                  color = "#A23B72", linewidth = 1.2) +
        # Confidence interval
        geom_ribbon(data = plot_data, aes(x = Dose_uM, ymin = lwr, ymax = upr), 
                    alpha = 0.2, fill = "#A23B72") +
        # Scale and labels
        scale_x_log10(
          breaks = c(0.01, 0.1, 1, 10, 100),
          labels = c("10nM", "100nM", "1µM", "10µM", "100µM"),
          limits = c(0.01, 100)
        ) +
        labs(
          title = paste0("Linear Model: ", tab),
          subtitle = paste0("log(Expression) ~ log(Dose) | R² = ", sprintf("%.3f", r_squared), 
                           " | p = ", sprintf("%.4f", slope_pval),
                           if(slope_pval < 0.001) " (***)" else if(slope_pval < 0.01) " (**)" else if(slope_pval < 0.05) " (*)" else " (ns)"),
          x = paste0(treatment_type, " Concentration"),
          y = "Relative Gene Expression"
        ) +
        publication_theme() +
        theme(
          plot.title = element_text(size = 14, face = "bold"),
          plot.subtitle = element_text(size = 11)
        )
      
      print(linear_plot)
      
      # Save the plot
      ggsave(paste0('./linear_model/', tab, '_linear_model.png'), plot = linear_plot,
             width = 8, height = 6, units = 'in', dpi = 600, bg = "white")
      
      # Save PDF if enabled
      if (enable_pdf_saving) {
        tryCatch({
          if (cairo_available) {
            ggsave(paste0('./linear_model/', tab, '_linear_model.pdf'), plot = linear_plot,
                   width = 8, height = 6, units = 'in', dpi = 600, bg = "white",
                   device = cairo_pdf)
          } else {
            ggsave(paste0('./linear_model/', tab, '_linear_model.pdf'), plot = linear_plot,
                   width = 8, height = 6, units = 'in', dpi = 600, bg = "white")
          }
        }, error = function(e) {
          message("PDF saving failed for linear model, continuing with PNG only: ", e$message)
        })
      }
      
      # Print comprehensive model summary
      cat("\n", strrep("=", 60), "\n")
      cat("LINEAR MODEL ANALYSIS RESULTS FOR:", tab, "\n")
      cat(strrep("=", 60), "\n")
      cat("Treatment type:", treatment_type, "\n")
      cat("Sample size:", n_points, "data points (excluding baseline)\n")
      cat("Model: Linear regression on log-transformed dose\n\n")
      
      cat("STATISTICAL MEASURES:\n")
      cat("• Slope:", sprintf("%.4f", slope))
      if (slope_pval < 0.001) cat(" (***)")
      else if (slope_pval < 0.01) cat(" (**)")
      else if (slope_pval < 0.05) cat(" (*)")
      else cat(" (ns)")
      cat("\n")
      cat("• 95% Confidence Interval: [", sprintf("%.4f", slope_ci_lower), ", ", sprintf("%.4f", slope_ci_upper), "]\n")
      cat("• P-value (slope ≠ 0):", sprintf("%.4f", slope_pval), "\n")
      cat("• R-squared:", sprintf("%.4f", r_squared), "\n")
      cat("• Adjusted R-squared:", sprintf("%.4f", adj_r_squared), "\n")
      
      cat("\nMODEL PARAMETERS:\n")
      cat("• Intercept:", sprintf("%.4f", intercept), "\n")
      cat("• Dose range:", min(dat_linear$Dose_uM), "-", max(dat_linear$Dose_uM), "µM\n")
      cat("• Response range:", sprintf("%.3f", min(dat_linear$Fold_gene_expression)), "-", 
          sprintf("%.3f", max(dat_linear$Fold_gene_expression)), "\n")
      cat(strrep("=", 60), "\n")
      
      # Print detailed model summary
      print(paste("Detailed linear model summary for", tab, ":"))
      print(model_summary)
      
    }, error = function(e) {
      warning(paste("Could not fit linear model for", tab, ":", e$message))
    })
  }
}

# Final summary message
cat("\n", strrep("=", 80), "\n")
cat("ANALYSIS COMPLETE!\n")
cat(strrep("=", 80), "\n")
cat("STATISTICAL INFORMATION AVAILABLE:\n")
cat("• Console output: Detailed statistics printed above for each gene\n")
cat("• Dose-response summary: './Statistical_summaries/dose_response_statistics.txt'\n")
cat("• Linear model summary: './Statistical_summaries/linear_model_statistics.txt'\n")
cat("• Dose-response plots: Enhanced with p-values, R², EC50, and AIC\n")
cat("• Linear model plots: With slope, confidence intervals, and R²\n")
cat("• Plot annotations: Statistical info displayed prominently on each curve\n\n")
cat("KEY STATISTICAL MEASURES INCLUDED:\n")
cat("DOSE-RESPONSE ANALYSIS:\n")
cat("• P-value: Significance of dose-response relationship\n")
cat("• R²: Goodness of fit (0-1, higher = better fit)\n")
cat("• EC50: Effective concentration for 50% response\n")
cat("• AIC: Akaike Information Criterion (model quality)\n\n")
cat("LINEAR MODEL ANALYSIS:\n")
cat("• Slope: Rate of change per log unit dose\n")
cat("• 95% Confidence Interval: Uncertainty range for slope\n")
cat("• R²: Proportion of variance explained by linear relationship\n")
cat("• P-value: Significance of linear dose-response relationship\n\n")
cat("• Sample size and dose/response ranges for both analyses\n")
cat(strrep("=", 80), "\n")

  



