consecutive_CpG_finder <- function(data, include_non_consecutive = FALSE) {
  # Ensure the data is sorted by chr and start
  data <- data[order(data$chr, data$start), ]
  
  # Initialize variables to store results
  chr_result <- c()
  start_result <- c()
  end_result <- c()
  consecutive_result <- c()
  
  i <- 1
  while (i < nrow(data)) {
    # Start of a potential sequence
    start_chr <- data$chr[i]
    start_val <- data$start[i]
    
    # Initialize count for consecutive rows
    consecutive_count <- 1
    end_val <- start_val
    
    # Check for consecutive rows
    while (i < nrow(data) && data$chr[i + 1] == data$chr[i] && data$start[i + 1] == data$start[i] + 2) {
      consecutive_count <- consecutive_count + 1
      end_val <- data$start[i + 1]
      i <- i + 1
    }
    
    # Store the result, depending on the value of include_non_consecutive
    if (consecutive_count > 1 || include_non_consecutive) {
      chr_result <- c(chr_result, start_chr)
      start_result <- c(start_result, start_val)
      end_result <- c(end_result, end_val)
      consecutive_result <- c(consecutive_result, consecutive_count)
    }
    
    # Move to the next row
    i <- i + 1
  }
  
  # Create a result data frame with 'start', 'end', and 'consecutive'
  result <- data.frame(chr = chr_result, start = start_result, end = end_result, consecutive = consecutive_result)
  
  return(result)
}
