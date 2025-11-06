#
#
# Read in Pupillometry Data -----------------------------------------------
#
#

# Get the Pupil file names ------------------------------------------------

files = list.files(
  path = 'data/Experiment 2/'
  , pattern = '.edf'
  , full.names = T
  , recursive = T
)

# Read in the data --------------------------------------------------------


overwrite_rds = FALSE # If set to TRUE will always re-process data, if not, loads if able
e2_study_pd = NULL # Where study data is stored
for(edf_path in files)
{
  sid = str_split(edf_path, "/", simplify = TRUE)[6] # SID from file path
  
  # Check if the file exists; if so, read it, if not create and save it
  rds_file = str_replace(edf_path, 'edf', 'rds')
  if(file.exists(rds_file) & !overwrite_rds)
  {
    message(sprintf("File Exists. Reading file %s", rds_file))
    pupil_data = readRDS(str_replace(edf_path, 'edf', 'rds'))
  } else {
    
    # Read in the data
    message(sprintf("Reading file %s", edf_path))
    edf_data <- eyelinkReader::read_edf(edf_path, import_samples=TRUE)
    
    # Store current sampling rate
    sampling_rate = edf_data$recordings$sample_rate[1]
    
    # Store eye measured
    sampled_eye = edf_data$recordings$eye[1]
    message(sprintf("    Measuring %s Eye", sampled_eye))
    
    # Store Events
    events = edf_data$events %>%
      select(trial, time, sttime, entime, sttime_rel, entime_rel, message)
    
    # Create trial conditions
    message("    Reading Trial Meta Data...")
    trial_cond = extract_trial_metadata_delayed(events) %>% filter(phase == 'study')
    
    # Grab the samples
    message("    Accessing samples...")
    pupil_data <- edf_data$samples %>%
      mutate(paR = ifelse(eye=='RIGHT', paR, paL),
             gxR = ifelse(eye=='RIGHT', gxR, gxL),
             gyR = ifelse(eye=='RIGHT', gyR, gyL)) %>%
      select(trial, time, eye, paR, gxR, gyR) %>%
      rename(rps = paR,
             x = gxR, y = gyR) %>%
      mutate(sid = sid)
    
    # Run the basic processing steps
    pupil_data = pupil_data %>%
      { message("    Marking blinks..."); identity(mark_blinks(., edf_data$blinks)) } %>%
      { message("    Detecting dropouts..."); identity(detect_dropouts(., sampling_rate=sampling_rate)) } %>%
      { message("    Detecting artifacts..."); identity(detect_artifacts(., sampling_rate=sampling_rate)) } %>%
      { message("    Marking data missing due to gaze..."); identity(missing_gaze(.)) } %>%
      { message("    Marking unreliable data..."); identity(mark_unreliable_data(.)) } %>%
      { message("    Marking bad trials (missing data)..."); identity(remove_excessive_missing(.)) } %>%
      { message("    Interpolating missing data..."); identity(interpolate_pupil(.)) } %>% 
      { message("    Apply a 4 hz low-pass filter..."); identity(lowpass_filter(., sampling_rate=sampling_rate)) } %>%
      { message("    Combining meta data..."); identity(left_join(., trial_cond)) } %>%
      { message("    Epoching data..."); identity(epoch_data(., lower_end=-200, upper_end=7000)) } %>%
      { message("    Baselining data..."); identity(baseline_pupil(.)) } %>%
      { message("    Downsampling to 50Hz..."); identity(downsample_pupil(.)) } %>%
      { message("    Flag bad trials..."); identity(bad_trials_by_missing_data(., missing_threshold_per_bin=.8)) } %>%
      { message("    Creating normalized values..."); identity(normalize_pupil_zscore(.)) }
    
    message("    Plotting study phase trials...")
    
    # Create a plot of the trials without condition labels for bad trial inspection
    g1 = pupil_data %>% filter(phase=='study') %>% ggplot(aes(x = bin, y = pupil_z)) +
      facet_wrap(~ trial) +
      geom_point(aes(color = pmissing)) +
      scale_color_gradient(low = "green", high = "red") +
      theme_minimal()
    
    ggsave(g1, file=sprintf("data/Experiment 2/Pupil_Plots/%s_study.pdf", sid), width=12, height=12)
    
    # Save the file for easy read-in
    saveRDS(pupil_data, str_replace(edf_path, 'edf', 'rds'))
  }
  
  # Save only the test phase data
  e2_study_pd = bind_rows(e2_study_pd, pupil_data %>% filter(phase=='study'))
}

# Check Bad Trials --------------------------------------------------------

# Calculate the % bad_trials (flagged with the earlier threshold)
# also store the % pmissing
e2_study_pd %>%
  group_by(sid, trial) %>%
  summarize(bad_trials = mean(bad_trial), pmissing = mean(pmissing)) %>%
  group_by(sid) %>%
  summarize(bad_trials = mean(bad_trials), pmissing = mean(pmissing)) -> e2_bad_trial_summary


# Exclude Bad Participants ------------------------------------------------

# Identify all participants with bad trials greater than the threshold
e2_too_many_bad_trials = e2_bad_trial_summary %>% filter(bad_trials > bad_trial_threshold_per_part) %>% pull(sid)

# sub_4_AB reported being unfocused due to sleep deprivation
# sub_19 didn't finish the experiment (withdrew and not counted)
# sub_21 blinked constantly due to an eye condition
# sub_22 technical issues with the data file producing strange wavy data
e2_flagged_participants = c(
  'sub_4_AB', 'sub_19', 'sub_21', 'sub_22'
)

e2_bad_part = unique(c(e2_flagged_participants))

# Exclude Bad Trials ------------------------------------------------------

# Remove all trials flagged as bad, also convert some variables to factors
e2_study_pd = e2_study_pd %>%
  filter(!bad_trial, !(sid %in% e2_bad_part)) %>%
  mutate(condition = factor(condition), sid = factor(sid), word = factor(word))

# Summary of bad trials
e2_bad_trial_summary %>% filter(!(sid %in% e2_bad_part)) %>% summarize_if(is.numeric, list(mean=mean, sd=sd))

e2_study_pd %>%
  group_by(sid, trial) %>%
  summarize(bad_trials = mean(bad_trial), pmissing = mean(pmissing)) %>%
  group_by(sid) %>%
  summarize(bad_trials = mean(bad_trials), pmissing = mean(pmissing)) %>%
  summarize(m=mean(pmissing), sd=sd(pmissing))

#
#
# Read in Behavioural Data ------------------------------------------------
#
#

# Get the Behavioural file names ------------------------------------------

e2_test_dat_dir = 'data/Experiment 2/Behavioural'

# Read-in data
e2_test_dat = NULL
for(l in list.files(e2_test_dat_dir, pattern='csv', recursive=TRUE, full.names = TRUE))
{
  print(l)
  dat = read.csv(l, header=TRUE) %>%
    mutate(sid=extract_participant_id(l), exp='e2') %>%
    select(sid, exp, word, condition, soa=SOA, key=response, testtrial = count_new_sketchpad_2) %>%
    mutate(said_yes = as.numeric(key > 3),
           condition = recode(condition, 'Plus'='aloud', 'Exes'='silent', 'foil'='foil'),
           word = tolower(word))
  
  e2_test_dat = rbind(e2_test_dat, dat)
}

# Exclude the same test phase participants as pupil participants
e2_test_dat %>%
  filter(!(sid %in% e2_bad_part), !is.na(testtrial)) %>%
  mutate(soa = ifelse(is.na(soa), 'foil', soa), soa = recode(soa, '495'='early', '1495'='late', '3495'='catch', 'foil'='foil'), key=as.numeric(key)) %>%
  mutate(soa = factor(soa, levels=c('foil', 'early', 'late', 'catch')),
         condition=factor(condition, levels=c('foil', 'silent', 'aloud'))) -> e2_test_dat

e2_test_dat %>%
  group_by(sid) %>%
  summarize(m=mean(key > 3), n=n())-> e2_test_sum

#
#
# Combine Pupil and Behavioural -------------------------------------------
#
#

# Study phase and test phase data use different SIDs, this makes them comparable
test_to_study_id_convert = c('PS001.csv'='PS001', 'PS002.csv'='PS002', 'PS003.csv'='PS003', 'PS004.csv'='PS004',
                             'PS006.csv'='PS006', 'PS007.csv'='PS007', 'PS008.csv'='PS008', 'PS009.csv'='PS009',
                             'PS010.csv'='PS010', 'PS011.csv'='PS011', 'PS012.csv'='PS012', 'PS013.csv'='PS013',
                             'PS014.csv'='PS014', 'PS015.csv'='PS015', 'PS016.csv'='PS016', 'PS017.csv'='PS017',
                             'subject-4.csv'='subject_4', 'subject-1'='sub_1_AA', 'subject-10.csv'='sub_10', 'subject-11.csv'='sub_11',
                             'subject-12.csv'='sub_12', 'subject-13.csv'='sub_13', 'subject-14.csv'='sub_14', 'subject-15.csv'='sub_15',
                             'subject-17.csv'='sub_17', 'subject-18.csv'='sub_18', 'subject-19.csv'='sub_19', 'subject-20.csv'='sub_20',
                             'subject-21.csv'='sub_22', 'subject-22.csv'='sub_22', 'subject-23.csv'='sub_23', 'subject-24.csv'='sub_24',
                             'subject-27.csv'='sub_27', 'subject-28.csv'='sub_28', 'subject-29.csv'='sub_29', 'subject-30.csv'='sub_30',
                             'subject-31.csv'='sub_31', 'subject-4'='sub_4_AB', 'subject-5'='sub_5_AC', 'subject-6'='sub_6_BA',
                             'subject-7'='sub_7_BB', 'subject-8'='sub_8_BC', 'subject-9.csv'='sub_9')

e2_test_dat = e2_test_dat %>% mutate(sid=recode(sid, !!!test_to_study_id_convert)) %>% filter(!(sid %in% e2_bad_part))

# Join study and test phase data
e2_study_pd = e2_study_pd %>%
  mutate(sid=as.character(sid), condition=as.character(condition), word=as.character(word)) %>%
  left_join(e2_test_dat) %>%
  mutate(mem = ifelse(said_yes==0, 'no', 'yes'), condmem = paste(condition, mem))

# Create a summary data frame
e2_study_pd_sum = e2_study_pd %>%
  group_by(sid, bin, condition, soa) %>%
  summarize(pupil_z = mean(pupil_z))

