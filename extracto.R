library(tuneR)
library(seewave)
library(dplyr)
library(lubridate)

# Define frequency bands (in Hz) up to 16 kHz
frequency_bands <- c(12.5, 16, 20, 25, 31.5, 40, 50, 63, 80, 100, 125, 160,
                     200, 250, 315, 400, 500, 630, 800, 1000, 1250, 1600, 
                     2000, 2500, 3150, 4000, 5000, 6300, 8000, 10000, 12500, 16000)

# --- Helper Function: extract file start time from filename ---
get_file_start_time <- function(wav_file) {
  base_name <- basename(wav_file)
  # Expected pattern: anything then _YYYYMMDD_HHMMSS.wav
  time_str <- sub(".*_(\\d{8}_\\d{6})\\.wav$", "\\1", base_name)
  as.POSIXct(time_str, format = "%Y%m%d_%H%M%S", tz = "UTC")
}

# --- Helper Function: get candidate WAV file that covers a given target time ---
get_candidate_wav <- function(folder_path, target_time) {
  wav_files <- list.files(path = folder_path, pattern = "\\.wav$", full.names = TRUE)
  if(length(wav_files) == 0) return(NULL)
  
  info_list <- lapply(wav_files, function(wav) {
    start <- get_file_start_time(wav)
    wav_obj <- readWave(wav)
    duration <- length(wav_obj@left) / wav_obj@samp.rate
    list(file = wav, start = start, end = start + seconds(duration))
  })
  
  candidate <- NULL
  for(info in info_list) {
    if(info$start <= target_time && info$end >= target_time) {
      candidate <- info
      break
    }
  }
  return(candidate)
}

# --- Helper Function: extract a spectral point at a given relative time ---
# target_time here is relative (seconds from file start)
extract_spectral_point <- function(wav_file, target_time, frequency_bands, wl = 1024, ovlp = 50) {
  wav <- readWave(wav_file)
  if(wav@stereo) wav <- mono(wav, "left")
  spec <- spectro(wav, f = wav@samp.rate, wl = wl, ovlp = ovlp, plot = FALSE)
  
  # Convert frequency bands from Hz to kHz for comparison with spec$freq
  target_bins <- sapply(frequency_bands, function(fb) {
    which.min(abs(spec$freq - (fb/1000)))
  })
  
  # Interpolate amplitude at target_time for each frequency bin
  spectral_values <- sapply(target_bins, function(idx) {
    approx(x = spec$time, y = spec$amp[idx, ], xout = target_time, rule = 2)$y
  })
  
  df <- as.data.frame(t(spectral_values))
  colnames(df) <- paste0("Freq_", frequency_bands, "_Hz")
  return(df)
}

# --- Main Loop ---
# final_results will store the combined data.
final_results <- data.frame()

# Loop over each row of total_table (each SPL reading, assumed one per second)
for (r in 1:nrow(total_table)) {
  spl_time <- total_table$Timestamp[r]  # SPL meter timestamp (POSIXct)
  cat("Processing SPL time:", spl_time, "\n")
  
  found_strip <- FALSE
  
  # Try each strip from 1 to 5 (adjust if needed)
  for (strip in 1:5) {
    cat("  Checking strip:", strip, "\n")
    candidate_list <- list()
    all_ARU_found <- TRUE
    
    # For each ARU in the strip (ARU 1 to 4)
    for (aru in 1:4) {
      # Folder for ARU: "../arus/{strip}{aru}/Data"
      folder_num <- as.character(strip * 10 + aru)
      folder_path <- file.path("../arus", folder_num, "Data")
      
      candidate <- get_candidate_wav(folder_path, spl_time)
      if (is.null(candidate)) {
        cat("    No WAV file in folder", folder_path, "covers SPL time", spl_time, "\n")
        all_ARU_found <- FALSE
        break
      }
      # Compute relative time within the candidate file
      candidate$rel_time <- as.numeric(difftime(spl_time, candidate$start, units = "secs"))
      candidate_list[[aru]] <- candidate
      cat("    Folder", folder_path, "covers SPL time. Relative time:", candidate$rel_time, "sec\n")
    } # end ARU loop
    
    # If all four ARUs have candidate files covering spl_time, extract spectral data.
    if (all_ARU_found) {
      spectral_list <- list()
      for (aru in 1:4) {
        cand <- candidate_list[[aru]]
        sp <- extract_spectral_point(cand$file, cand$rel_time, frequency_bands)
        # Rename spectral columns with ARU prefix (keep the value, no Timestamp inside sp)
        spec_cols <- colnames(sp)
        colnames(sp) <- paste0("ARU", aru, "_", spec_cols)
        spectral_list[[aru]] <- sp
      }
      
      # Merge spectral data for all ARUs side-by-side (they all have one row)
      merged_spec <- do.call(cbind, spectral_list)
      
      # Combine with the SPL row and add the strip number
      spl_row <- total_table[r, , drop=FALSE]
      spl_row$Strip <- strip
      final_row <- cbind(spl_row, merged_spec)
      final_results <- bind_rows(final_results, final_row)
      cat("  Using data from strip", strip, "for SPL time", spl_time, "\n")
      found_strip <- TRUE
      break  # use only one strip per SPL row
    }
  } # end strip loop
  
  if (!found_strip) {
    cat("  No complete ARU data found for SPL time", spl_time, "\n")
  }
}

# Display final results
print(final_results)