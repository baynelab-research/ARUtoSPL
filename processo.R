# Load required libraries
library(tuneR)
library(seewave)
library(dplyr)
library(beepr)
library(lubridate)
load("wip.Rdata")
create_new_pure_tone_starts <- function(spl_interval_df) {
  # spl_interval_df: a data frame for a given continuous interval.
  # It must have a 'start_time' column (POSIXct) and a number of rows that is an exact multiple of 3.
  
  n <- nrow(spl_interval_df)
  if(n %% 3 != 0) {
    stop("The number of rows in the interval is not a multiple of 3.")
  }
  
  n_triads <- n / 3
  
  # Create a new data frame with one row per triad
  new_pure_tone_starts <- data.frame(
    strip = 1:n_triads,         # This assigns strip numbers sequentially (1,2,...)
    global_start_row = NA,      # first row index of each triad
    global_end_row   = NA,      # last row index of each triad
    pure_tone_time   = as.POSIXct(rep(NA, n_triads), origin="1970-01-01", tz="UTC"),
    stringsAsFactors = FALSE
  )
  
  for(i in 1:n_triads) {
    idx_start <- (i - 1) * 3 + 1
    idx_end   <- i * 3
    # Use the first row's start_time of each triad as the pure tone start time
    pure_time <- spl_interval_df$start_time[idx_start]
    new_pure_tone_starts$pure_tone_time[i] <- pure_time
    new_pure_tone_starts$global_start_row[i] <- idx_start
    new_pure_tone_starts$global_end_row[i] <- idx_end
  }
  
  return(new_pure_tone_starts)
}

# Example usage:
# Use the continuous_intervals rows for instance four as your spl_interval_df:
example_SPL <- continuous_intervals[continuous_intervals$instance == 4, ]
new_pure_tone_starts <- create_new_pure_tone_starts(example_SPL)
print(new_pure_tone_starts)


# Function to extract file start time from the filename
get_file_start_time <- function(wav_file) {
  base_name <- basename(wav_file)
  time_str <- sub(".*_(\\d{8}_\\d{6})\\.wav$", "\\1", base_name)
  as.POSIXct(time_str, format = "%Y%m%d_%H%M%S", tz = "UTC")
}

# Function to get the 1kHz band values from a WAV file
get_1khz_values <- function(wav_file, wl = 1024, ovlp = 50) {
  wav <- readWave(wav_file)
  if(wav@stereo) wav <- mono(wav, "left")
  spec <- spectro(wav, f = wav@samp.rate, wl = wl, ovlp = ovlp, plot = FALSE)
  target_bin <- which.min(abs(spec$freq - 1))  # spec$freq in kHz
  amp <- spec$amp[target_bin, ]
  amp_db <- 20 * log10(amp + 1e-6)
  data.frame(Time = spec$time, Amplitude = amp, Amplitude_dB = amp_db)
}

# Function to detect the 1kHz spike time (in seconds from file start)
get_spike_absolute_time <- function(wav_file, wl = 1024, ovlp = 50) {
  file_start <- get_file_start_time(wav_file)
  one_khz_data <- get_1khz_values(wav_file, wl, ovlp)
  diff_amp <- diff(one_khz_data$Amplitude)
  spike_index <- which.max(diff_amp)
  spike_time <- one_khz_data$Time[spike_index]
  file_start + seconds(spike_time)
}

# Function to update new_pure_tone_starts with alignment times.
# For each row (strip), we read the single file from each corresponding folder 
# and store its absolute spike time in new columns: ARU1_align_time, ..., ARU4_align_time.
update_alignment_offsets <- function(new_pure_tone_starts, base_align_dir = "../arus/align/") {
  for (r in 1:nrow(new_pure_tone_starts)) {
    strip <- new_pure_tone_starts$strip[r]
    for (aru in 1:4) {
      folder_num <- strip * 10 + aru  # e.g. strip 1, ARU1 -> 11
      folder_path <- file.path(base_align_dir, as.character(folder_num))
      wav_files <- list.files(path = folder_path, pattern = "\\.wav$", full.names = TRUE)
      if (length(wav_files) == 0) {
        cat("No WAV file found in folder", folder_path, "\n")
        next
      }
      # Use the only file in the folder
      wav_file <- wav_files[1]
      spike_time_abs <- get_spike_absolute_time(wav_file)
      col_name <- paste0("ARU", aru, "_align_time")
      new_pure_tone_starts[r, col_name] <- spike_time_abs
      cat("Strip", strip, "ARU", aru, "from folder", folder_num, "spike time:", spike_time_abs, "\n")
    }
  }
  return(new_pure_tone_starts)
}

tone_starts <- update_alignment_offsets(new_pure_tone_starts)

# Create a simple offset table by computing the differences
offset_table <- tone_starts %>%
  mutate(
    ARU1_offset = as.numeric(difftime(ARU1_align_time, pure_tone_time, units = "secs")),
    ARU2_offset = as.numeric(difftime(ARU2_align_time, pure_tone_time, units = "secs")),
    ARU3_offset = as.numeric(difftime(ARU3_align_time, pure_tone_time, units = "secs")),
    ARU4_offset = as.numeric(difftime(ARU4_align_time, pure_tone_time, units = "secs"))
  ) %>%
  select(strip, pure_tone_time, ARU1_offset, ARU2_offset, ARU3_offset, ARU4_offset)

print(offset_table)


###############
#############

#############################
#############

library(tuneR)
library(seewave)
library(dplyr)
library(lubridate)

# Parameter: number of seconds to remove from the beginning and end of the interval
shrink_sec <- 10

# Function to extract the file’s start time from its filename.
# Assumes filename contains a segment _YYYYMMDD_HHMMSS.wav
get_file_start_time <- function(wav_file) {
  base_name <- basename(wav_file)
  time_str <- sub(".*_(\\d{8}_\\d{6})\\.wav$", "\\1", base_name)
  as.POSIXct(time_str, format = "%Y%m%d_%H%M%S", tz = "UTC")
}

# We'll assume your offset table is stored in a data frame 'offset_table'
# with columns: strip, pure_tone_time, ARU1_offset, ARU2_offset, ARU3_offset, ARU4_offset.

# Now loop over each continuous interval (from continuous_intervals)
# and check, for each strip in offset_table, whether all ARU files in "../arus/{strip}{aru}/Data/"
# cover the (shrunken) effective interval after applying their offset.
# When one strip is found, record that interval & strip.
results <- data.frame(
  IntervalStart = as.POSIXct(character()),
  IntervalEnd = as.POSIXct(character()),
  EffectiveStart = as.POSIXct(character()),
  EffectiveEnd = as.POSIXct(character()),
  Strip = integer(),
  stringsAsFactors = FALSE
)

for (i in 1:nrow(continuous_intervals)) {
  interval <- continuous_intervals[i, ]
  interval_start <- interval$start_time
  interval_end <- interval$end_time
  cat("Processing interval:", interval_start, "to", interval_end, "\n")
  
  # Define the effective (shrunken) interval
  effective_start <- interval_start + seconds(shrink_sec)
  effective_end <- interval_end - seconds(shrink_sec)
  
  strip_found <- FALSE
  
  # Loop over each strip available in your offset table
  for (r in 1:nrow(offset_table)) {
    strip <- offset_table$strip[r]
    cat("  Checking strip:", strip, "\n")
    all_ARU_ok <- TRUE
    
    # Process each ARU (1 to 4) for this strip
    for (aru in 1:4) {
      # Folder path follows the pattern: "../arus/{strip}{aru}/Data/"
      folder_path <- sprintf("../arus/%d%d/Data/", strip, aru)
      cat("    Checking ARU", aru, "folder:", folder_path, "\n")
      
      wav_files <- list.files(path = folder_path, pattern = "\\.wav$", full.names = TRUE)
      if (length(wav_files) == 0) {
        cat("      No WAV files found in", folder_path, "\n")
        all_ARU_ok <- FALSE
        break
      }
      # Use the only file in the folder
      wav_file <- wav_files[1]
      
      # Get file start and compute file end
      file_start <- get_file_start_time(wav_file)
      wav_obj <- readWave(wav_file)
      duration <- length(wav_obj@left) / wav_obj@samp.rate
      file_end <- file_start + seconds(duration)
      
      # Retrieve ARU-specific offset from offset_table for this strip and ARU
      offset_col <- paste0("ARU", aru, "_offset")
      offset_val <- offset_table[[offset_col]][r]
      
      # Adjust the effective (shrunken) interval by the ARU offset
      adjusted_start <- effective_start + seconds(offset_val)
      adjusted_end <- effective_end + seconds(offset_val)
      cat("      ARU", aru, "adjusted effective interval:", adjusted_start, "to", adjusted_end, "\n")
      
      # Check if the file covers the entire adjusted effective interval
      if (!(file_start <= adjusted_start && file_end >= adjusted_end)) {
        cat("      File", basename(wav_file), "does NOT cover the adjusted effective interval.\n")
        all_ARU_ok <- FALSE
        break
      } else {
        cat("      File", basename(wav_file), "covers the adjusted effective interval.\n")
      }
    } # end ARU loop
    
    if (all_ARU_ok) {
      cat("  Found complete ARU recordings for strip", strip, "\n")
      results <- rbind(results, data.frame(
        IntervalStart = interval_start,
        IntervalEnd = interval_end,
        EffectiveStart = effective_start,
        EffectiveEnd = effective_end,
        Strip = strip,
        stringsAsFactors = FALSE
      ))
      strip_found <- TRUE
      break  # use only one strip per interval
    } else {
      cat("  Incomplete data for strip", strip, "; trying next strip.\n")
    }
  } # end strip loop
  
  if (!strip_found) {
    cat("No complete ARU data found for interval:", interval_start, "to", interval_end, "\n")
  }
}

print(results)


##############

