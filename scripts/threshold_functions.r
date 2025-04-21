#===========================================================================
#                            Visual Inspection on Conditional Scores
#===========================================================================
vics_method <- function(response_matrix, rt_matrix, min_threshold = 1.0) {
  # Check input dimensions
  if (!is.matrix(response_matrix) || !is.matrix(rt_matrix)) {
    stop("Both inputs must be matrices")
  }
  
  if (nrow(response_matrix) != nrow(rt_matrix) || ncol(response_matrix) != ncol(rt_matrix)) {
    stop("Response matrix and RT matrix must have the same dimensions")
  }
  
  # Number of items
  J <- ncol(response_matrix)
  
  # Create a vector to store thresholds for each item
  thresholds <- numeric(J)
  
  # Process each item
  for (j in 1:J) {
    # Extract response and RT data for current item
    responses <- response_matrix[, j]
    rts <- rt_matrix[, j]
    
    # Remove NAs (missing responses)
    valid_indices <- !is.na(responses) & !is.na(rts)
    item_responses <- responses[valid_indices]
    item_rts <- rts[valid_indices]
    
    # Sort by response time
    sorted_indices <- order(item_rts)
    sorted_responses <- item_responses[sorted_indices]
    sorted_rts <- item_rts[sorted_indices]
    
    # Find threshold using VICS method
    # Calculate proportion correct in sliding windows
    window_size <- max(20, round(length(sorted_responses) * 0.05)) # Minimum window size
    
    # Calculate the chance level (could be modified based on item type)
    chance_level <- 0.25  # Default for 4-option multiple choice
    if (all(item_responses %in% c(0, 1))) {
      # For binary items, chance level is typically 0.5 for 2-option items or lower for open responses
      chance_level <- 0.2  # Can be adjusted based on item format
    }
    
    # Find the threshold
    threshold <-  min_threshold  # Default to minimum RT if no threshold is found
    
    for (i in window_size:length(sorted_rts)) {
      window_start <- max(1, i - window_size + 1)
      window_end <- i
      
      # Calculate proportion correct in current window
      prop_correct <- mean(sorted_responses[window_start:window_end])
      
      # If proportion correct rises above chance level
      if (prop_correct > chance_level) {
        threshold <- sorted_rts[window_start]
        break
      }
    }
    
    thresholds[j] <- threshold
  }
  
  return(thresholds)
}

#===========================================================================
#                            Visual Inspection on Item Information 
#===========================================================================

vii_thresholds <- function(response_matrix, response_time_matrix, 
                          reference_corr = 0.2, 
                          show_plots = FALSE,
                          bins = 10) {
  
  # Check input dimensions
  if(nrow(response_matrix) != nrow(response_time_matrix) || 
     ncol(response_matrix) != ncol(response_time_matrix)) {
    stop("Dimensions of response matrix and response time matrix must match")
  }
  
  n_students <- nrow(response_matrix)
  n_items <- ncol(response_matrix)
  
  # Calculate provisional proficiency estimates (simple sum score)
  provisional_scores <- rowSums(response_matrix, na.rm = TRUE)
  
  # Initialize results storage
  thresholds <- numeric(n_items)
  
  # Process each item
  for(j in 1:n_items) {
    # Get response times for this item
    item_times <- response_time_matrix[, j]
    # Get responses for this item
    item_responses <- response_matrix[, j]
    
    # Remove NAs (if any)
    valid_indices <- !is.na(item_times) & !is.na(item_responses)
    valid_times <- item_times[valid_indices]
    valid_responses <- item_responses[valid_indices]
    valid_scores <- provisional_scores[valid_indices]
    
    if(length(valid_times) == 0) {
      thresholds[j] <- NA
      next
    }
    
    # Sort by response time
    sorted_indices <- order(valid_times)
    sorted_times <- valid_times[sorted_indices]
    sorted_responses <- valid_responses[sorted_indices]
    sorted_scores <- valid_scores[sorted_indices]
    
    # Create time bins 
    time_breaks <- unique(quantile(sorted_times, probs = seq(0, 1, length.out = bins + 1)))
    bin_indices <- cut(sorted_times, breaks = time_breaks, include.lowest = TRUE)
    
    # Calculate correlation within each bin
    bin_correlations <- numeric(length(time_breaks) - 1)
    bin_median_times <- numeric(length(time_breaks) - 1)
    
    for(b in 1:length(bin_correlations)) {
      bin_mask <- bin_indices == levels(bin_indices)[b]
      if(sum(bin_mask) > 1 && var(sorted_responses[bin_mask]) > 0) {
        bin_correlations[b] <- cor(sorted_responses[bin_mask], sorted_scores[bin_mask])
      } else {
        bin_correlations[b] <- NA
      }
      bin_median_times[b] <- median(sorted_times[bin_mask])
    }
    
    # Find where correlation exceeds reference value
    threshold_bin <- which(bin_correlations >= reference_corr)[1]
    
    if(is.na(threshold_bin)) {
      # If no bin exceeds reference correlation, use minimum time
      thresholds[j] <- min(sorted_times)
    } else {
      # Use the lower bound of the bin where correlation exceeds reference
      thresholds[j] <- time_breaks[threshold_bin]
    }
    
    # Optional plotting
    if(show_plots) {
      par(mfrow = c(1, 1))
      plot(bin_median_times, bin_correlations, type = "b", 
           main = paste("Item", j, "- VII Method"),
           xlab = "Response Time", ylab = "Correlation with Proficiency",
           ylim = c(-0.1, 1))
      abline(h = reference_corr, col = "red", lty = 2)
      abline(v = thresholds[j], col = "blue", lty = 2)
      legend("bottomright", 
             legend = c("Correlation", "Reference Corr", "Threshold"),
             col = c("black", "red", "blue"),
             lty = c(1, 2, 2))
    }
  }
  
  return(thresholds)
}
