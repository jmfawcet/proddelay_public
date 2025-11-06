
# Exclude Bad Participants ------------------------------------------------

# There were none
e1_flagged_participants = c(
  ''
)

e1_bad_part = unique(c(e1_flagged_participants))


#
#
# Read in Behavioural Data ------------------------------------------------
#
#

# Get the Behavioural file names ------------------------------------------

e1_test_dat_dir = 'data/Experiment 1'

# Read-in data
e1_test_dat = NULL
for(l in list.files(e1_test_dat_dir, pattern='csv', recursive=TRUE, full.names = TRUE))
{
  print(l)
  dat = read.csv(l, header=TRUE) %>%
    mutate(sid=extract_participant_id(l), exp='e1') %>%
    select(sid, exp, word, condition, soa=SOA, key=response, testtrial = count_new_sketchpad_2) %>%
    mutate(said_yes = as.numeric(key > 3),
           condition = recode(condition, 'Plus'='aloud', 'Exes'='silent', 'foil'='foil'),
           word = tolower(word))
  
  e1_test_dat = rbind(e1_test_dat, dat)
}

# Exclude bad participants (there were none)
e1_test_dat %>%
  filter(!(sid %in% e1_bad_part), !is.na(testtrial)) %>%
  mutate(soa = ifelse(is.na(soa), 'foil', soa), soa = recode(soa, '495'='early', '1495'='late', '3495'='catch', 'foil'='foil'), key=as.numeric(key)) %>%
  mutate(soa = factor(soa, levels=c('foil', 'early', 'late', 'catch')),
         condition=factor(condition, levels=c('foil', 'silent', 'aloud'))) -> e1_test_dat
