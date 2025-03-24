library(tuneR)
library(seewave)
library(dplyr)
library(tidyr)
library(ggplot2)

total_table <- read.csv("total_table.csv")

# Convert Timestamp column back to POSIXct datetime
total_table$Timestamp <- as.POSIXct(total_table$Timestamp, format="%Y-%m-%d %H:%M:%S")

# Ensure the table is sorted by Timestamp
total_table <- total_table %>% arrange(Timestamp)

# Calculate differences between consecutive timestamps (in seconds)
diff_seconds <- as.numeric(diff(total_table$Timestamp), units = "secs")

# Identify indices where the gap is more than one second
gap_indices <- which(diff_seconds > 1)

# The start of each continuous segment: first row and row following each gap
segment_starts <- c(1, gap_indices + 1)
# The end of each continuous segment: row at each gap and the final row
segment_ends <- c(gap_indices, nrow(total_table))

# Build a summary table with the start and end timestamps and sample count per segment
continuous_intervals <- data.frame(
  start_time = total_table$Timestamp[segment_starts],
  end_time   = total_table$Timestamp[segment_ends],
  num_samples = segment_ends - segment_starts + 1
)

# View the continuous intervals table
print(continuous_intervals)


##33
##33
##33 make heatmaps

convert_freq <-
function(freq_str) {
  # remove any spaces
  freq_str <- gsub(" ", "", freq_str)
  if (grepl("kHz$", freq_str)) {
    # Remove "kHz", convert to numeric and multiply by 1000
    return(as.numeric(sub("kHz$", "", freq_str)) * 1000)
  } else if (grepl("Hz$", freq_str)) {
    # Remove "Hz" and convert directly
    return(as.numeric(sub("Hz$", "", freq_str)))
  } else {
    return(NA)
  }
}


for(q in 1:length(continuous_intervals[,1]))
{
	start_row<-continuous_intervals[q,1]
	end_row<-continuous_intervals[q,2]
	heat_data <- total_table %>%
	  slice(start_row:end_row) %>%
	  select(-Timestamp, -TickCount)

	# Create a time index column (representing the order/time)
	heat_data <- heat_data %>% 
	  mutate(TimeIndex = row_number())

	# Pivot to long format: each row now represents a measurement for one frequency band at a given time index.
	heat_long <- heat_data %>%
	  pivot_longer(cols = -TimeIndex, names_to = "Frequency", values_to = "LogValue")

	# Create a numeric frequency column (in Hz) using our conversion function
	heat_long <- heat_long %>%
	  mutate(Freq_Hz = sapply(Frequency, convert_freq))

	# Reorder the Frequency factor levels based on the numeric frequency values
	heat_long$Frequency <- factor(heat_long$Frequency, 
		                      levels = unique(heat_long$Frequency[order(heat_long$Freq_Hz)]))

	# Plot the rotated heatmap: x-axis = time, y-axis = frequency bands in proper order
pe<-ggplot(heat_long, aes(x = TimeIndex, y = Frequency, fill = LogValue)) +
	  geom_tile() +
	  scale_fill_viridis_c() +
	  theme_minimal() +
	  labs(x = "Time Index", y = "Frequency Band", fill = "Log Value") +
	  theme(axis.text.x = element_text(angle = 90, hjust = 1))
plot(pe)
}

##33
##33
##33

detect_1kHz_tone <- function(file_path, threshold = 0.9, nfft = 1024, overlap = 512) {
  # Read the WAV file
  wav <- readWave(file_path)
  
  # Convert to mono if stereo
  if (wav@stereo) {
    wav <- mono(wav, "left")
  }
  
  # Compute the spectrogram; note: set noisereduction = 1 to avoid the error
  spec <- spectro(wav, f = wav@samp.rate, n = nfft, 
                  ovlp = (overlap/nfft)*100, plot = FALSE, norm = FALSE, 
                  noisereduction = 1)
  
  # Since spec$freq is in kHz (e.g., 1.000 represents 1 kHz), target is 1
  freq_bin <- which.min(abs(spec$freq - 1))
  
  # Extract the amplitude over time for the 1 kHz bin
  amp_1kHz <- spec$amp[freq_bin, ]
  
  # Normalize amplitude by its maximum
  norm_amp <- amp_1kHz / max(amp_1kHz, na.rm = TRUE)
  
  # Detect the first time point where the normalized amplitude exceeds the threshold
  onset_index <- which(norm_amp > threshold)[1]
  
  if (is.na(onset_index)) {
    warning("No 1 kHz onset detected.")
    return(NA)
  }
  
  # Retrieve the corresponding time from the spectrogram's time vector (in seconds)
  onset_time <- spec$time[onset_index]
  
  return(onset_time)
}

####### Example usage:
file_path <- "../arus/test.wav"
onset_time <- detect_1kHz_tone(file_path)
print(onset_time)


# Function to process WAV files for a given strip and ARU
process_aru_files <- function(strip, aru, pure_tone_starts) {
  # Define the folder path for the ARU
  folder_path <- sprintf("./%d%d/Data/", strip, aru)
  
  # List all WAV files in the folder
  wav_files <- list.files(folder_path, pattern = "\\.wav$", full.names = TRUE)
  
  # Filter the pure_tone_starts data frame for the current strip
  tone_times <- pure_tone_starts %>% filter(strip == !!strip)
  
  # Initialize a vector to store offsets
  offsets <- numeric(nrow(tone_times))
  
  # Loop over each expected pure tone time
  for (i in seq_len(nrow(tone_times))) {
    spl_time <- as.POSIXct(tone_times$pure_tone_time[i])
    
    # Find the WAV file that likely contains the SPL time
    wav_file <- wav_files[which.min(abs(spl_time - as.POSIXct(sub(".*_(\\d{8}_\\d{6})\\.wav$", "\\1", basename(wav_files)), format = "%Y%m%d_%H%M%S")))]
    
    # Detect the 1 kHz tone onset in the selected WAV file
    onset_time <- detect_1khz_tone(wav_file)
    
    if (!is.na(onset_time)) {
      # Extract the start time from the WAV file name
      wav_start_time <- as.POSIXct(sub(".*_(\\d{8}_\\d{6})\\.wav$", "\\1", basename(wav_file)), format = "%Y%m%d_%H%M%S")
      
      # Calculate the ARU's recorded tone time
      aru_tone_time <- wav_start_time + onset_time
      
      # Calculate the offset between the ARU and SPL meter
      offsets[i] <- as.numeric(difftime(aru_tone_time, spl_time, units = "secs"))
    } else {
      offsets[i] <- NA
    }
  }
  
  # Add the offsets to the pure_tone_starts data frame
  pure_tone_starts[[paste0("ARU", aru, "_offset")]] <- offsets
  
  return(pure_tone_starts)
}

