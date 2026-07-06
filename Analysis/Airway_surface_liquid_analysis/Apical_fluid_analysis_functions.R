# Script for functions to support apical fluid analysis
# Sam Siljee
# 27th June 2025

# Function to get matrix offset values
get_alignment_offsets <- function(image_matrix, reference_row, min_offset, max_offset) {
  # Set some variables
  n_rows <- nrow(image_matrix)

  # Initialise vector for offset results
  best_offsets <- numeric()

  # Loop through rows
  for (x in 1:n_rows) {
    # Initialise offset score for this row
    best_score <- Inf

    # Set image row
    image_row <- image_matrix[x, ]

    # Loop through offset values
    for (o in min_offset:max_offset) {
      # Calculate offset score
      score <- score_alignment(image_row, reference_row, o)
      
      # Compare calculated score to best score, and update if better
      if (score < best_score) {
        best_score <- score
        minimum_index <- o
      }
    }

    # Add minimum score in the minimum offsets vector
    best_offsets <- c(best_offsets, minimum_index)
  }
  
  # return the results
  return(best_offsets)
  
}

# Function to score row alignment
score_alignment <- function(image_row, reference_row, x_offset) {
  # Get original row length
  row_length <- length(image_row)

  # Create a row using the supplied offset
  if (x_offset > 0) { # Positive offset
    # Trim the end of the row
    offset_row <- image_row[1:(row_length - x_offset)]
  } else if (x_offset < 0) { # Negative offset
    # Trim the start of the row
    offset_row <- image_row[(-x_offset + 1):row_length]
  } else { # No offset
    offset_row <- image_row
  }

  # Create a reference row of the same length
  if (x_offset > 0) { # Positive offset
    # Trim the start of the row
    trimmed_reference_row <- reference_row[(x_offset + 1):row_length]
  } else if (x_offset < 0) { # Negative offset
    # Trim the end of the row
    trimmed_reference_row <- reference_row[1:(row_length + x_offset)]
  } else { # No offset
    trimmed_reference_row <- reference_row
  }

  # Return the mean absolute differences
  return(mean(abs(trimmed_reference_row - offset_row)))
}

# Function to align a matrix given an image matrix and set of offset values
align_matrix <- function(image_matrix, offsets) {
  # Set some values
  min_offset <- min(offsets)
  max_offset <- max(offsets)
  
  # Calculate new width given offset values
  new_width <- ncol(image_matrix) - (max_offset - min_offset)
  
  # Intialise new matrix
  aligned_matrix <- matrix(nrow = nrow(image_matrix), ncol = new_width)
  
  # Loop through rows
  for(x in 1:nrow(image_matrix)) {
    # Calculate start position in original matrix for this row
    start_col <- offsets[x] - min_offset + 1
    end_col <- start_col + new_width - 1
    
    # Extract the section from original matrix
    aligned_matrix[x,] <- image_matrix[x, start_col:end_col]
  }
  
  # Return the aligned matrix
  return(aligned_matrix)
}

# Function to normalise the image for .png export
normalise_image <- function(input_matrix) {
  normalised_matrix <- (input_matrix - min(input_matrix)) / (max(input_matrix) - min(input_matrix))
  return(normalised_matrix)
}
