# Plot changes in effect sizes of LPS on Belize vampire bat social networks
# Gerry Carter and Simon Ripperger

# clear workspace
rm(list=ls())

# set directory
setwd("~/Dropbox/Dropbox/_working/_ACTIVE/belize_LPS/2018_Belize_analysis")
setwd("C:/Users/simon.ripperger/Dropbox/2018_Belize_analysis")
.libPaths("C:/R libraries")

# load packages
library(tidyverse)
library(boot)
library(lubridate)
library(lme4)
library(lmerTest) 
library(igraph)
library(cowplot)

# run simulations to get new data or plot existing data??
GET_NEW_DATA <- F
GET_NEW_DATA <- T

# or load simulation data or run simulations----
if(GET_NEW_DATA==F){
  load('simulation_results_2020-03-06_08_23.Rdata')
  }else{

  # load data from main R script
    load('results_2020-09-15_21_02.Rdata')      
    
  # EFFECT OF NUMBER OF MEETINGS ON LPS EFFECT SIZE ----- 
  # define association using MANY RSSI values ----
  thresholds <- 
    quantile(Belize$RSSI, probs = seq(from = 76, to = 98, by = 2)/100)
  thresholds
  
  # create df to store effect sizes
  effect_sizes <- tibble(RSSI= names(thresholds), effect_size_degree= NA, n.meetings= NA, p=NA)
  
  
  # for loop ----
  start <- Sys.time()
  options(dplyr.show_progress = T)
  for (i in 1:length(thresholds)) {
    RSSI_threshold = thresholds[i] #set RSSI threshold 
    #threshold for Current Biology paper was -26dmb at 90%; here -27dbm at 85%
    
    # clean data 
    df <- 
      Belize %>%
      filter(RSSI > RSSI_threshold) %>%
      filter(SenderID %in% BatIDs ) %>%
      filter(EncounteredID %in% BatIDs ) %>% 
      select(- PacketID, -ChunkID) %>%
      mutate(dyad= if_else(SenderID<EncounteredID, 
                           paste(SenderID,EncounteredID, sep="_"), 
                           paste(EncounteredID,SenderID, sep="_"))) %>%
      mutate(EndOfMeeting = StartOfMeeting + MeetingDuration) %>%
      group_by(dyad) %>%
      arrange(StartOfMeeting) %>% 
      mutate(indx = c(0, cumsum(as.numeric(lead(StartOfMeeting)) >
                                  cummax(as.numeric(EndOfMeeting)))[-n()])) %>%
      group_by(dyad, indx) %>%
      summarise(StartOfMeeting = min(StartOfMeeting), 
                EndOfMeeting = max(EndOfMeeting),
                duration = difftime(EndOfMeeting, StartOfMeeting,unit = "secs"), 
                RSSI = max(RSSI)) %>%
      mutate(bat1 = sub( "_.*$", "", dyad ), bat2 = sub('.*_', '', dyad)) %>%
      mutate(duration = as.numeric(duration, unit = "secs")) 
    
    # count meetings
    n.meetings <- nrow(df)
    
    # insert hour breaks
    # convert start and end times to interval
    df$interval <- interval(df$StartOfMeeting, df$EndOfMeeting)
    #insert "event" column with series of numbers
    df$event <- c(1:nrow(df))
    
    # create function to get hours within a time interval
    get_hours <- function(event, StartOfMeeting, EndOfMeeting){
      hours <-  seq(StartOfMeeting-minute(StartOfMeeting)*60-second(StartOfMeeting),
                    EndOfMeeting-minute(EndOfMeeting)*60-second(EndOfMeeting),
                    "hour")
      dateseq <- hours
      dateseq[1] <-  StartOfMeeting
      r <-  c(dateseq, EndOfMeeting)
      dur <-  as.numeric(difftime(r[-1], r[-length(r)], unit = 'secs'))
      data.frame(event, hour = hours, duration = dur)}
    
    # create new events with event durations within each hour
    df2 <- 
      df %>%
      rowwise %>%
      do(get_hours(.$event, .$StartOfMeeting, .$EndOfMeeting)) %>%
      ungroup() %>%
      group_by(event, hour) %>%
      summarize(duration = sum(duration)) %>% 
      as.data.frame()
    
    # match original start time back into new events
    df2$StartOfMeeting <- df$StartOfMeeting[match(df2$event, df$event)]
    # if start time is past the hour use that start time, otherwise use the hour slot as the start time
    df2$StartOfMeeting <- if_else(df2$StartOfMeeting>df2$hour, df2$StartOfMeeting, df2$hour)
    # match original end time back into new events
    df2$EndOfMeeting <- df$EndOfMeeting[match(df2$event, df$event)]
    # if end time is before the next hour (start hour+ 1 hour), use that end time, otherwise use the next hour
    df2$EndOfMeeting <- if_else(df2$EndOfMeeting<(df2$hour+3600), df2$EndOfMeeting, df2$hour+3600)
    # match other data back in
    df2$dyad <- df$dyad[match(df2$event, df$event)]
    df2$RSSI <- df$RSSI[match(df2$event, df$event)]
    
    # set end of meeting
    # set timezone to BelizeTime
    # set start of study to 3pm on 25th
    df <- 
      df2 %>%
      mutate(bat1 = sub( "_.*$", "", dyad ), bat2 = sub('.*_', '', dyad)) %>%
      mutate(StartBelizeTime = force_tz(StartOfMeeting, tzone = "CST6CDT")) %>%
      mutate(EndBelizeTime = force_tz(EndOfMeeting, tzone = "CST6CDT")) %>%
      mutate(StartBelizeTime= StartBelizeTime - hours(8), EndBelizeTime=EndBelizeTime - hours(8)) %>%
      filter(StartBelizeTime >= as.POSIXct("2018-04-25 15:00:00", tz = "CST6CDT")) %>% 
      select(StartBelizeTime,bat1,bat2,duration, RSSI)
    
    # match treatment from bats table into df
    df$treatment_bat1 <- bats$treatment[match(df$bat1, bats$sensor_node)]
    df$treatment_bat2 <- bats$treatment[match(df$bat2, bats$sensor_node)]
    
    # remove other dataframe
    rm(df2)
    
    # define treatment and post-treatment periods
    d <- 
      df %>% 
      mutate(datetime= as.POSIXct(StartBelizeTime,tz = "CST6CDT")) %>% 
      mutate(treatment_period= datetime >= treatment_start & datetime < treatment_stop) %>% 
      mutate(post24_period= datetime >= post24_start & datetime < post24_stop) %>% 
      mutate(post48_period= datetime >= post48_start & datetime < post48_stop) %>% 
      mutate(hour= substring(datetime, 1,13))
    
    # get treated bats
    treated.bats <- 
      d %>% 
      mutate(treated= ifelse(treatment_bat1=="LPS", bat1, 
                             ifelse(treatment_bat2=="LPS", bat2, NA))) %>% 
      filter(!is.na(treated)) %>% 
      pull(treated) %>% 
      unique()
    
    # make six-hour network for each period----
    network.treatment <-
      d %>% 
      filter(treatment_period) %>% 
      mutate(dyad= paste(bat1,bat2, sep="_")) %>%
      select(dyad, hour, duration) %>% 
      group_by(dyad, hour) %>% 
      summarize(duration= sum(duration)) %>% 
      group_by(dyad) %>% 
      summarize(duration= sum(duration)) %>% 
      ungroup() %>% 
      a_b_edgelist_to_matrix(directed=F) %>%
      graph_from_adjacency_matrix("undirected", weighted=T, diag=F) 
    
    network.post24 <- 
      d %>% 
      filter(post24_period) %>% 
      mutate(dyad= paste(bat1,bat2, sep="_")) %>%
      select(dyad, hour, duration) %>% 
      group_by(dyad, hour) %>% 
      summarize(duration= sum(duration)) %>% 
      group_by(dyad) %>% 
      summarize(duration= sum(duration)) %>% 
      ungroup() %>% 
      a_b_edgelist_to_matrix(directed=F) %>%
      graph_from_adjacency_matrix("undirected", weighted=T, diag=F) ###
    
    network.post48 <- 
      d %>% 
      filter(post48_period) %>% 
      mutate(dyad= paste(bat1,bat2, sep="_")) %>%
      select(dyad, hour, duration) %>% 
      group_by(dyad, hour) %>% 
      summarize(duration= sum(duration)) %>% 
      group_by(dyad) %>% 
      summarize(duration= sum(duration)) %>% 
      ungroup() %>% 
      a_b_edgelist_to_matrix(directed=F) %>%
      graph_from_adjacency_matrix("undirected", weighted=T, diag=F) ###
    
    # compare network centrality by treatment period
    centrality.treatment <- 
      tibble(bat= names(degree(network.treatment)), 
             degree= degree(network.treatment),
             strength= strength(network.treatment),
             eigenvector= eigen_centrality(network.treatment)$vector,
             period= "during treatment")
    centrality.post24 <- 
      tibble(bat= names(degree(network.post24)), 
             degree= degree(network.post24),
             strength= strength(network.post24),
             eigenvector= eigen_centrality(network.post24)$vector,
             period= "post-treatment") %>% 
      # define two post-treatment periods or one 
      mutate(period= paste(period, "24 hours later")) %>% 
      ungroup()
    centrality.post48 <- 
      tibble(bat= names(degree(network.post48)), 
             degree= degree(network.post48),
             strength= strength(network.post48),
             eigenvector= eigen_centrality(network.post48)$vector,
             period= "post-treatment") %>% 
      # define two post periods or one
      mutate(period= paste(period, "48 hours later")) %>% 
      ungroup()
    
    # get centrality data per period and bat
    data <- 
      rbind(centrality.treatment, centrality.post24, centrality.post48) %>% 
      mutate(treated= bat %in% treated.bats) %>% 
      mutate(treatment= ifelse(treated, "LPS", "saline"))
    
    # get effect sizes (slope)
    effect_sizes$effect_size_degree[i] <- 
      summary(lm(degree~ treated, data= data[which(data$period=="during treatment"),]))$coefficients[2,1]  #effect during treatment
    
    # p-values
    effect_sizes$p[i] <- 
      summary(lm(degree~ treated, data= data[which(data$period=="during treatment"),]))$coefficients[2,4] 
    
    # store number of meetings
    effect_sizes$n.meetings[i] <- n.meetings
    
    # print progress
    print(paste(i, "of", length(thresholds)))
  }
  
  stop <- Sys.time()
  loop1time <- stop-start
  loop1time
  
  
  # EFFECT OF MEETING DURATION ON LPS EFFECT SIZE-----
  # define association using one RSSI value ----
  thresholds <- quantile(Belize$RSSI, probs = c(85)/100)
  
  # define minimum duration using many values----
  min.durations <- seq(from=0, to= 1200, by= 60)
  
  # create df to store effect sizes
  effect_sizes2 <- tibble(duration= min.durations, effect_size_degree= NA, n.meetings= NA, p=NA)
  
  # for loop  ----
  start <- Sys.time()
  
  
  for (i in 1:length(min.durations)) {
    RSSI_threshold = thresholds[1] #set RSSI threshold 
    
    # clean data and filter by duration
    df <- 
      Belize %>%
      filter(RSSI > RSSI_threshold) %>%
      filter(SenderID %in% BatIDs ) %>%
      filter(EncounteredID %in% BatIDs ) %>% 
      select(- PacketID, -ChunkID) %>%
      mutate(dyad= if_else(SenderID<EncounteredID, 
                           paste(SenderID,EncounteredID, sep="_"), 
                           paste(EncounteredID,SenderID, sep="_"))) %>%
      mutate(EndOfMeeting = StartOfMeeting + MeetingDuration) %>%
      group_by(dyad) %>%
      arrange(StartOfMeeting) %>% 
      mutate(indx = c(0, cumsum(as.numeric(lead(StartOfMeeting)) >
                                  cummax(as.numeric(EndOfMeeting)))[-n()])) %>%
      group_by(dyad, indx) %>%
      summarise(StartOfMeeting = min(StartOfMeeting), 
                EndOfMeeting = max(EndOfMeeting),
                duration = difftime(EndOfMeeting, StartOfMeeting,unit = "secs"), 
                RSSI = max(RSSI)) %>%
      mutate(bat1 = sub( "_.*$", "", dyad ), bat2 = sub('.*_', '', dyad)) %>%
      mutate(duration = as.numeric(duration, unit = "secs")) %>% 
      filter(duration > min.durations[i]) %>% 
      ungroup()
    
    # count meetings
    n.meetings <- nrow(df)
    
    # insert hour breaks
    # convert start and end times to interval
    df$interval <- interval(df$StartOfMeeting, df$EndOfMeeting)
    #insert "event" column with series of numbers
    df$event <- c(1:nrow(df))
    
    # create function to get hours within a time interval
    get_hours <- function(event, StartOfMeeting, EndOfMeeting){
      hours <-  seq(StartOfMeeting-minute(StartOfMeeting)*60-second(StartOfMeeting),
                    EndOfMeeting-minute(EndOfMeeting)*60-second(EndOfMeeting),
                    "hour")
      
      dateseq <- hours
      dateseq[1] <-  StartOfMeeting
      r <-  c(dateseq, EndOfMeeting)
      dur <-  as.numeric(difftime(r[-1], r[-length(r)], unit = 'secs'))
      data.frame(event, hour = hours, duration = dur)}
    
    # create new events with event durations within each hour
    df2 <- 
      df %>%
      rowwise %>%
      do(get_hours(.$event, .$StartOfMeeting, .$EndOfMeeting)) %>%
      ungroup() %>%
      group_by(event, hour) %>%
      summarize(duration = sum(duration)) %>% 
      as.data.frame()
    
    # match original start time back into new events
    df2$StartOfMeeting <- df$StartOfMeeting[match(df2$event, df$event)]
    # if start time is past the hour use that start time, otherwise use the hour slot as the start time
    df2$StartOfMeeting <- if_else(df2$StartOfMeeting>df2$hour, df2$StartOfMeeting, df2$hour)
    # match original end time back into new events
    df2$EndOfMeeting <- df$EndOfMeeting[match(df2$event, df$event)]
    # if end time is before the next hour (start hour+ 1 hour), use that end time, otherwise use the next hour
    df2$EndOfMeeting <- if_else(df2$EndOfMeeting<(df2$hour+3600), df2$EndOfMeeting, df2$hour+3600)
    # match other data back in
    df2$dyad <- df$dyad[match(df2$event, df$event)]
    df2$RSSI <- df$RSSI[match(df2$event, df$event)]
    
    # set end of meeting
    # set timezone to BelizeTime
    # set start of study to 3pm on 25th
    df <- 
      df2 %>%
      mutate(bat1 = sub( "_.*$", "", dyad ), bat2 = sub('.*_', '', dyad)) %>%
      mutate(StartBelizeTime = force_tz(StartOfMeeting, tzone = "CST6CDT")) %>%
      mutate(EndBelizeTime = force_tz(EndOfMeeting, tzone = "CST6CDT")) %>%
      mutate(StartBelizeTime= StartBelizeTime - hours(8), EndBelizeTime=EndBelizeTime - hours(8)) %>%
      filter(StartBelizeTime >= as.POSIXct("2018-04-25 15:00:00", tz = "CST6CDT")) %>% 
      select(StartBelizeTime,bat1,bat2,duration, RSSI)
    
    # match treatment from bats table into df
    df$treatment_bat1 <- bats$treatment[match(df$bat1, bats$sensor_node)]
    df$treatment_bat2 <- bats$treatment[match(df$bat2, bats$sensor_node)]
    
    # remove other dataframe
    rm(df2)
    
    # define treatment and post-treatment periods
    d <- 
      df %>% 
      mutate(datetime= as.POSIXct(StartBelizeTime,tz = "CST6CDT")) %>% 
      mutate(treatment_period= datetime >= treatment_start & datetime < treatment_stop) %>% 
      mutate(post24_period= datetime >= post24_start & datetime < post24_stop) %>% 
      mutate(post48_period= datetime >= post48_start & datetime < post48_stop) %>% 
      mutate(hour= substring(datetime, 1,13))
    
    # get treated bats
    treated.bats <- 
      d %>% 
      mutate(treated= ifelse(treatment_bat1=="LPS", bat1, 
                             ifelse(treatment_bat2=="LPS", bat2, NA))) %>% 
      filter(!is.na(treated)) %>% 
      pull(treated) %>% 
      unique()
    
    # make six-hour network for each period----
    network.treatment <-
      d %>% 
      filter(treatment_period) %>% 
      mutate(dyad= paste(bat1,bat2, sep="_")) %>%
      select(dyad, hour, duration) %>% 
      group_by(dyad, hour) %>% 
      summarize(duration= sum(duration)) %>% 
      group_by(dyad) %>% 
      summarize(duration= sum(duration)) %>% 
      ungroup() %>% 
      a_b_edgelist_to_matrix(directed=F) %>%
      graph_from_adjacency_matrix("undirected", weighted=T, diag=F) 
    
    network.post24 <- 
      d %>% 
      filter(post24_period) %>% 
      mutate(dyad= paste(bat1,bat2, sep="_")) %>%
      select(dyad, hour, duration) %>% 
      group_by(dyad, hour) %>% 
      summarize(duration= sum(duration)) %>% 
      group_by(dyad) %>% 
      summarize(duration= sum(duration)) %>% 
      ungroup() %>% 
      a_b_edgelist_to_matrix(directed=F) %>%
      graph_from_adjacency_matrix("undirected", weighted=T, diag=F) ###
    
    network.post48 <- 
      d %>% 
      filter(post48_period) %>% 
      mutate(dyad= paste(bat1,bat2, sep="_")) %>%
      select(dyad, hour, duration) %>% 
      group_by(dyad, hour) %>% 
      summarize(duration= sum(duration)) %>% 
      group_by(dyad) %>% 
      summarize(duration= sum(duration)) %>% 
      ungroup() %>% 
      a_b_edgelist_to_matrix(directed=F) %>%
      graph_from_adjacency_matrix("undirected", weighted=T, diag=F) ###
    
    # compare network centrality by treatment period
    centrality.treatment <- 
      tibble(bat= names(degree(network.treatment)), 
             degree= degree(network.treatment),
             strength= strength(network.treatment),
             eigenvector= eigen_centrality(network.treatment)$vector,
             period= "during treatment")
    centrality.post24 <- 
      tibble(bat= names(degree(network.post24)), 
             degree= degree(network.post24),
             strength= strength(network.post24),
             eigenvector= eigen_centrality(network.post24)$vector,
             period= "post-treatment") %>% 
      # define two post-treatment periods or one 
      mutate(period= paste(period, "24 hours later")) %>% 
      ungroup()
    centrality.post48 <- 
      tibble(bat= names(degree(network.post48)), 
             degree= degree(network.post48),
             strength= strength(network.post48),
             eigenvector= eigen_centrality(network.post48)$vector,
             period= "post-treatment") %>% 
      # define two post periods or one
      mutate(period= paste(period, "48 hours later")) %>% 
      ungroup()
    
    # get centrality data per period and bat
    data <- 
      rbind(centrality.post24, centrality.post48, centrality.treatment) %>% 
      mutate(treated= bat %in% treated.bats) %>% 
      mutate(treatment= ifelse(treated, "LPS", "saline")) %>% 
      # combine post periods
      mutate(period= ifelse(period=="during treatment", period, "post-treatment"))
    
    # get effect sizes (slope)
    effect_sizes2$effect_size_degree[i] <- 
      summary(lm(degree~ treated, data= data[which(data$period=="during treatment"),]))$coefficients[2,1]  #effect during treatment
    
    # p-values
    effect_sizes2$p[i] <- summary(lm(degree~ treated, data= data[which(data$period=="during treatment"),]))$coefficients[2,4] 
    
    # store number of meetings
    effect_sizes2$n.meetings[i] <- n.meetings
    
    # print progress
    print(paste(i, "of", length(min.durations)))
  }
  stop <- Sys.time()
  loop2time <- stop-start
  loop2time
  
  
  
  # EFFECT OF RSSI ON LPS EFFECT SIZE (AT CONSTANT SAMPLE SIZE)-----
  
  # choose number of effect sizes to sample at each RSSI level
  replicates <- 200
  
  # define association using MANY RSSI values ----
  thresholds <- 
    quantile(Belize$RSSI, probs = rep(seq(from = 76, to = 94, by = 2)/100, each=replicates) )
  thresholds
  unique(thresholds)
  
  # choose number of meetings to randomly sample 
  # set number to be 95% of the number of meetings at the largest RSSI quantile
  random_N <- 
    Belize %>%
    filter(RSSI > thresholds[length(thresholds)]) %>%
    filter(SenderID %in% BatIDs ) %>%
    filter(EncounteredID %in% BatIDs ) %>% 
    select(- PacketID, -ChunkID) %>%
    mutate(dyad= if_else(SenderID<EncounteredID, 
                         paste(SenderID,EncounteredID, sep="_"), 
                         paste(EncounteredID,SenderID, sep="_"))) %>%
    mutate(EndOfMeeting = StartOfMeeting + MeetingDuration) %>%
    group_by(dyad) %>%
    arrange(StartOfMeeting) %>% 
    mutate(indx = c(0, cumsum(as.numeric(lead(StartOfMeeting)) >
                                cummax(as.numeric(EndOfMeeting)))[-n()])) %>%
    group_by(dyad, indx) %>%
    summarise(StartOfMeeting = min(StartOfMeeting), 
              EndOfMeeting = max(EndOfMeeting),
              duration = difftime(EndOfMeeting, StartOfMeeting,unit = "secs"), 
              RSSI = max(RSSI)) %>%
    mutate(bat1 = sub( "_.*$", "", dyad ), bat2 = sub('.*_', '', dyad)) %>%
    mutate(duration = as.numeric(duration, unit = "secs")) %>% 
    ungroup() %>% 
    sample_frac(size=0.95) %>% 
    nrow() 
  random_N
  
  # define minimum duration
  min.durations <- 0
  
  # create df to store effect sizes
  effect_sizes3 <- tibble(RSSI= names(thresholds), effect_size_degree= NA, n.meetings= NA, p=NA)
  
  # for loop ----
  start <- Sys.time()
  options(dplyr.show_progress = F)
  
  for (i in 1:length(thresholds)) {
    RSSI_threshold = thresholds[i] #set RSSI threshold 
    
    # clean data and filter by duration
    df <- 
      Belize %>%
      filter(RSSI > RSSI_threshold) %>% ###
      filter(SenderID %in% BatIDs ) %>%
      filter(EncounteredID %in% BatIDs ) %>% 
      select(- PacketID, -ChunkID) %>%
      mutate(dyad= if_else(SenderID<EncounteredID, 
                           paste(SenderID,EncounteredID, sep="_"), 
                           paste(EncounteredID,SenderID, sep="_"))) %>%
      mutate(EndOfMeeting = StartOfMeeting + MeetingDuration) %>%
      group_by(dyad) %>%
      arrange(StartOfMeeting) %>% 
      mutate(indx = c(0, cumsum(as.numeric(lead(StartOfMeeting)) >
                                  cummax(as.numeric(EndOfMeeting)))[-n()])) %>%
      group_by(dyad, indx) %>%
      summarise(StartOfMeeting = min(StartOfMeeting), 
                EndOfMeeting = max(EndOfMeeting),
                duration = difftime(EndOfMeeting, StartOfMeeting,unit = "secs"), 
                RSSI = max(RSSI)) %>%
      mutate(bat1 = sub( "_.*$", "", dyad ), bat2 = sub('.*_', '', dyad)) %>%
      mutate(duration = as.numeric(duration, unit = "secs")) %>% 
      ungroup() %>% 
      sample_n(size= random_N, replace = F)
    
    # count meetings
    n.meetings <- nrow(df)
    
    # insert hour breaks
    # convert start and end times to interval
    df$interval <- interval(df$StartOfMeeting, df$EndOfMeeting)
    #insert "event" column with series of numbers
    df$event <- c(1:nrow(df))
    
    # create function to get hours within a time interval
    get_hours <- function(event, StartOfMeeting, EndOfMeeting){
      hours <-  seq(StartOfMeeting-minute(StartOfMeeting)*60-second(StartOfMeeting),
                    EndOfMeeting-minute(EndOfMeeting)*60-second(EndOfMeeting),
                    "hour")
      
      dateseq <- hours
      dateseq[1] <-  StartOfMeeting
      r <-  c(dateseq, EndOfMeeting)
      dur <-  as.numeric(difftime(r[-1], r[-length(r)], unit = 'secs'))
      data.frame(event, hour = hours, duration = dur)}
    
    # create new events with event durations within each hour
    df2 <- 
      df %>%
      rowwise %>%
      do(get_hours(.$event, .$StartOfMeeting, .$EndOfMeeting)) %>%
      ungroup() %>%
      group_by(event, hour) %>%
      summarize(duration = sum(duration)) %>% 
      as.data.frame()
    
    # match original start time back into new events
    df2$StartOfMeeting <- df$StartOfMeeting[match(df2$event, df$event)]
    # if start time is past the hour use that start time, otherwise use the hour slot as the start time
    df2$StartOfMeeting <- if_else(df2$StartOfMeeting>df2$hour, df2$StartOfMeeting, df2$hour)
    # match original end time back into new events
    df2$EndOfMeeting <- df$EndOfMeeting[match(df2$event, df$event)]
    # if end time is before the next hour (start hour+ 1 hour), use that end time, otherwise use the next hour
    df2$EndOfMeeting <- if_else(df2$EndOfMeeting<(df2$hour+3600), df2$EndOfMeeting, df2$hour+3600)
    # match other data back in
    df2$dyad <- df$dyad[match(df2$event, df$event)]
    df2$RSSI <- df$RSSI[match(df2$event, df$event)]
    
    # set end of meeting
    # set timezone to BelizeTime
    # set start of study to 3pm on 25th
    df <- 
      df2 %>%
      mutate(bat1 = sub( "_.*$", "", dyad ), bat2 = sub('.*_', '', dyad)) %>%
      mutate(StartBelizeTime = force_tz(StartOfMeeting, tzone = "CST6CDT")) %>%
      mutate(EndBelizeTime = force_tz(EndOfMeeting, tzone = "CST6CDT")) %>%
      mutate(StartBelizeTime= StartBelizeTime - hours(8), EndBelizeTime=EndBelizeTime - hours(8)) %>%
      filter(StartBelizeTime >= as.POSIXct("2018-04-25 15:00:00", tz = "CST6CDT")) %>% 
      select(StartBelizeTime,bat1,bat2,duration, RSSI)
    
    # match treatment from bats table into df
    df$treatment_bat1 <- bats$treatment[match(df$bat1, bats$sensor_node)]
    df$treatment_bat2 <- bats$treatment[match(df$bat2, bats$sensor_node)]
    
    # remove other dataframe
    rm(df2)
    
    # define treatment and post-treatment periods
    d <- 
      df %>% 
      mutate(datetime= as.POSIXct(StartBelizeTime,tz = "CST6CDT")) %>% 
      mutate(treatment_period= datetime >= treatment_start & datetime < treatment_stop) %>% 
      mutate(post24_period= datetime >= post24_start & datetime < post24_stop) %>% 
      mutate(post48_period= datetime >= post48_start & datetime < post48_stop) %>% 
      mutate(hour= substring(datetime, 1,13))
    
    # get treated bats
    treated.bats <- 
      d %>% 
      mutate(treated= ifelse(treatment_bat1=="LPS", bat1, 
                             ifelse(treatment_bat2=="LPS", bat2, NA))) %>% 
      filter(!is.na(treated)) %>% 
      pull(treated) %>% 
      unique()
    
    # make six-hour network for each period----
    network.treatment <-
      d %>% 
      filter(treatment_period) %>% 
      mutate(dyad= paste(bat1,bat2, sep="_")) %>%
      select(dyad, hour, duration) %>% 
      group_by(dyad, hour) %>% 
      summarize(duration= sum(duration)) %>% 
      group_by(dyad) %>% 
      summarize(duration= sum(duration)) %>% 
      ungroup() %>% 
      a_b_edgelist_to_matrix(directed=F) %>%
      graph_from_adjacency_matrix("undirected", weighted=T, diag=F) 
    
    network.post24 <- 
      d %>% 
      filter(post24_period) %>% 
      mutate(dyad= paste(bat1,bat2, sep="_")) %>%
      select(dyad, hour, duration) %>% 
      group_by(dyad, hour) %>% 
      summarize(duration= sum(duration)) %>% 
      group_by(dyad) %>% 
      summarize(duration= sum(duration)) %>% 
      ungroup() %>% 
      a_b_edgelist_to_matrix(directed=F) %>%
      graph_from_adjacency_matrix("undirected", weighted=T, diag=F) ###
    
    network.post48 <- 
      d %>% 
      filter(post48_period) %>% 
      mutate(dyad= paste(bat1,bat2, sep="_")) %>%
      select(dyad, hour, duration) %>% 
      group_by(dyad, hour) %>% 
      summarize(duration= sum(duration)) %>% 
      group_by(dyad) %>% 
      summarize(duration= sum(duration)) %>% 
      ungroup() %>% 
      a_b_edgelist_to_matrix(directed=F) %>%
      graph_from_adjacency_matrix("undirected", weighted=T, diag=F) ###
    
    # compare network centrality by treatment period
    centrality.treatment <- 
      tibble(bat= names(degree(network.treatment)), 
             degree= degree(network.treatment),
             strength= strength(network.treatment),
             eigenvector= eigen_centrality(network.treatment)$vector,
             period= "during treatment")
    centrality.post24 <- 
      tibble(bat= names(degree(network.post24)), 
             degree= degree(network.post24),
             strength= strength(network.post24),
             eigenvector= eigen_centrality(network.post24)$vector,
             period= "post-treatment") %>% 
      # define two post-treatment periods or one 
      mutate(period= paste(period, "24 hours later")) %>% 
      ungroup()
    centrality.post48 <- 
      tibble(bat= names(degree(network.post48)), 
             degree= degree(network.post48),
             strength= strength(network.post48),
             eigenvector= eigen_centrality(network.post48)$vector,
             period= "post-treatment") %>% 
      # define two post periods or one
      mutate(period= paste(period, "48 hours later")) %>% 
      ungroup()
    
    # get centrality data per period and bat
    data <- 
      rbind(centrality.post24, centrality.post48, centrality.treatment) %>% 
      mutate(treated= bat %in% treated.bats) %>% 
      mutate(treatment= ifelse(treated, "LPS", "saline")) %>% 
      # combine post periods
      mutate(period= ifelse(period=="during treatment", period, "post-treatment"))
    
    # get effect sizes (slope)
    effect_sizes3$effect_size_degree[i] <- 
      summary(lm(degree~ treated, data= data[which(data$period=="during treatment"),]))$coefficients[2,1]  #effect during treatment
    
    # p-values
    effect_sizes3$p[i] <- summary(lm(degree~ treated, data= data[which(data$period=="during treatment"),]))$coefficients[2,4] 
    
    # store number of meetings
    effect_sizes3$n.meetings[i] <- n.meetings
    
    # print progress
    print(paste(i, "of", length(thresholds)))
    
  }
  stop <- Sys.time()
  loop3time <- stop-start
  loop3time
  
  # save output with timestamp-----
  timestamp <- substr(gsub(x=gsub(":","_",Sys.time()), 
                           pattern=" ", replace="_"), start=1, stop=16)
  save.image(file= paste("sim_results_", timestamp, ".Rdata", sep=""))
  
  # reset these options
  options(dplyr.show_progress = T)
  GET_NEW_DATA <- F #this is important
}

# time to run for loop code
loop1time # 8 min
loop2time # 13 min
loop3time # 11 hours

# plot LPS effect sizes with number of meetings----
N_meetings_plot <- 
  effect_sizes %>% 
  mutate(effect_size_degree= effect_size_degree*-1) %>% 
  mutate(pvalue= ifelse(p>0.05, 0.05, p)) %>% 
  mutate(significance= ifelse(p>0.05, "NS", "p<0.05")) %>% 
  select(-p) %>% 
  pivot_longer(cols= c(effect_size_degree, pvalue, n.meetings), names_to = "measure", values_to = "value") %>%
  mutate(measure= case_when(
    measure == "effect_size_degree" ~ "LPS treatment coefficient",
    measure == "n.meetings" ~ "number of dyadic encounters",
    measure == "pvalue" ~ "p-value (if < 0.05)")) %>% 
  ggplot(aes(x=RSSI, y=value, group= 1))+
  facet_wrap(~measure, scales= "free", ncol=1, 
             strip.position = "left")+
  geom_point(aes(shape=significance, color= significance), size=3) +
  geom_line()+
  geom_vline(xintercept = 5.5, color= 'blue', linetype='dashed', size=1)+
  xlab("minimum proximity index to define an encounter")+
  ylab("")+
  scale_color_manual(values= c("grey", "black"))+ 
  scale_shape_manual(values= c("triangle","circle"))+
  theme_cowplot()+
  theme(strip.background = element_blank(),
        strip.placement = "outside")
N_meetings_plot  


# plot effect sizes with meeting duration----
duration_plot <- 
  effect_sizes2 %>% 
  mutate(effect_size_degree= effect_size_degree*-1) %>% 
  mutate(pvalue= ifelse(p>0.05, 0.05, p)) %>% 
  mutate(significance= ifelse(p>0.05, "NS", "p<0.05")) %>% 
  select(-p) %>% 
  pivot_longer(cols= c(effect_size_degree, pvalue, n.meetings), names_to = "measure", values_to = "value") %>% 
  mutate(measure= case_when(
    measure == "effect_size_degree" ~ "LPS treatment coefficient",
    measure == "n.meetings" ~ "number of dyadic encounters",
    measure == "pvalue" ~ "p-value (if < 0.05)")) %>% 
  ggplot(aes(x=duration, y=value, group= 1))+
  facet_wrap(~measure, scales= "free", ncol=1, 
             strip.position = "left")+
  geom_point(aes(shape=significance, color= significance), size=3) +
  geom_line()+
  geom_vline(xintercept = 1, color= 'blue', linetype='dashed', size=1)+
  xlab("minimum seconds to define an encounter")+
  ylab("")+
  scale_color_manual(values= c("grey", "black"))+ 
  scale_shape_manual(values= c("triangle","circle"))+
  theme_cowplot()+
  theme(strip.background = element_blank(),
        strip.placement = "outside")
duration_plot

# plot effect sizes with RSSI at smallest sample size bin----

# get absolute effect size
effect_sizes3 <- 
  effect_sizes3 %>% 
  mutate(effect_size= effect_size_degree*-1) 

# get means and 95% CI
means_ci <- 
  effect_sizes3 %>% 
  select(RSSI, effect_size) %>% 
  boot_ci2(y= .$effect_size, x=.$RSSI)

RSSI_plot <- 
  effect_sizes3 %>% 
  ggplot(aes(x=RSSI, y=effect_size))+
  geom_jitter(alpha=0.2, width=0.1, color= "grey")+
  geom_violin(width=1, fill=NA)+
  geom_point(data= means_ci, aes(x=effect, y=mean), size=1)+
  geom_errorbar(data= means_ci, aes(ymin=low, ymax=high, width=.2, x=effect, y=mean), size=1)+
  ylim(c(0,5))+
  ylab("LPS treatment coefficient")+
  xlab("minimum proximity index to define an encounter")+
  ggtitle(paste("LPS effect size with", random_N, "encounters"))+
  theme_cowplot()
RSSI_plot

# combine plots
plot1 <-
  duration_plot+
  theme(legend.position= 'none')

plot2 <-
  N_meetings_plot+
  theme(legend.position= 'none')

plot3 <- 
  RSSI_plot+
  theme(plot.title = element_blank())

# make figure 3  
(sim1 <- plot_grid(plot1,plot2, ncol=2, labels= 'AUTO'))
ggsave("simulation1_results.pdf", width= 14, height= 8, units="in", dpi=1200)

(sim2 <- plot3)
ggsave("simulation2_results.pdf", width= 14, height= 7, units="in", dpi=1200)
