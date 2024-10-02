# Function to find the closest site with matching character column and return a table
Assigning_Closest_TSS <- function(TSS_Info, input_dataframe) {
  # Initialize an empty data frame to store the results
  results <- data.frame(
    start = numeric(),             # start for input dataframe
    chr = character(),             # chr column for both input and closest site
    closest_TSS = numeric(),       # TSS for TSS_Info table
    distance = numeric(),          # Distance between input start and closest TSS
    gene = character(),            # Gene column from TSS_Info
    stringsAsFactors = FALSE
  )
  
  # Loop through each target entry
  for (i in 1:nrow(input_dataframe)) {
    # Extract the current row from input_dataframe
    current_row <- input_dataframe[i, ]
    
    target_start <- current_row$start  # Use 'start' for input dataframe values
    target_chr <- current_row$chr      # Use 'chr' for character column
    
    # Filter the TSS_Info table to only those with matching character column
    filtered_sites <- TSS_Info[TSS_Info$chr == target_chr, ]
    
    # If no matching character rows are found in filtered_sites, skip this iteration
    if (nrow(filtered_sites) == 0) {
      next  # Skip if no matching chr is found
    }
    
    # Calculate the absolute difference between the input start and filtered_sites TSS
    differences <- abs(filtered_sites$TSS - target_start)
    
    # If there are no valid differences, skip to the next iteration
    if (length(differences) == 0) {
      next
    }
    
    # Find the index of the minimum difference
    closest_idx <- which.min(differences)
    
    # Ensure closest_idx is valid
    if (length(closest_idx) == 0) {
      next  # Skip if no valid closest index is found
    }
    
    # Extract the closest TSS value and the corresponding gene
    closest_TSS <- filtered_sites$TSS[closest_idx]
    gene <- filtered_sites$gene[closest_idx]  # Get the corresponding gene
    
    # Calculate the distance
    distance <- differences[closest_idx]
    
    # Create a new row with all columns from the current input row and additional results
    new_row <- cbind(current_row, closest_TSS = closest_TSS, distance = distance, gene = gene)
    
    # Append the new row to the results
    results <- rbind(results, new_row)
  }
  
  return(results)
}
