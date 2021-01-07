# Analyze effects of LPS on Belize vampire bat social networks
# Gerry Carter and Simon Ripperger

# # set directory
# setwd(dirname(file.choose()))
# setwd("~/Dropbox/Dropbox/_working/_ACTIVE/belize_LPS/2018_Belize_analysis")
# setwd("C:/Users/simon.ripperger/Dropbox/2018_Belize_analysis")
# .libPaths("C:/R libraries")

# clear workspace
rm(list=ls())

# load packages
library(tidyverse)
library(boot)
library(lubridate)
library(lme4)
library(lmerTest) 
library(igraph)
library(cowplot)
library(patchwork)
library(scales)

# choose number of permutations
perms <- 100
perms <- 10000

# make random numbers consistent 
set.seed(123)

# functions----

# get the mean and 95% CI of a vector by bootstrapping
boot_ci <- function(x, perms=5000) {
  mean.w=function(x,w) sum(x*w)
  numNA <- sum(is.na(x))
  x <- as.vector(na.omit(x))
  mean <- mean(x)
  boot <- boot.ci(boot(data=x, statistic=mean.w, R=perms, stype="w", parallel = "multicore", ncpus = 4), type="bca")
  low <- boot$bca[1,4]
  high <- boot$bca[1,5]
  c(low=low,mean=mean,high=high, N=round(length(x)))
}

# get the mean and 95% CI within groups of a dataframe by bootstrapping
boot_ci2 <- function(d=d, y=d$y, x=d$x, perms=5000){
  df <- data.frame(effect=unique(x))
  df$low <- NA
  df$mean <- NA
  df$high <- NA
  df$n.obs <- NA
  for (i in 1:nrow(df)) {
    ys <- y[which(x==df$effect[i])]
    if (length(ys)>1){
      b <- boot_ci(y[which(x==df$effect[i])], perms=perms) 
      df$low[i] <- b[1]
      df$mean[i] <- b[2]
      df$high[i] <- b[3]
      df$n.obs[i] <- b[4]
    }else{
      df$low[i] <- NA
      df$mean[i] <- ys
      df$high[i] <- NA
      df$n.obs[i] <- 1
    }
  }
  df
}

# convert a list of dyadic interactions to a sociomatrix
# (requires igraph)
a_b_edgelist_to_matrix <- function(el=el, symbol="_", directed= T, make.NA.zero=T){
  a <- str_split(as.data.frame(el)[,1],symbol, simplify = TRUE)[,1]
  r <- str_split(as.data.frame(el)[,1],symbol, simplify = TRUE)[,2]
  y <- as.data.frame(el)[,2]
  e <- data.frame(a,r,y, stringsAsFactors = F)
  if (make.NA.zero){
    g <- graph_from_data_frame(e, directed=directed)
    m <- get.adjacency(g, attr='y', sparse=FALSE)
    m
  }else{
    e$y <- e$y+1 # temporarily add one to distinguish between 0 and NA
    g <- graph_from_data_frame(e, directed=directed)
    m <- get.adjacency(g, attr='y', sparse=FALSE)
    m[m==0] <- NA # relabel missing values as NA
    m <- m-1 # subtract one to adjust values back
    m
  }
}

# convert interactions to rates (duration within time bin)
# this function adds zeros for possible associations/interactions within each time bin
# it assumes that every individual ("ids") is present in every time bin ("bin")
# it requires a df with dyad ("a_b"), bin, and duration

events_to_rates <- function(df= df, bin= df$bin, ids= ids, directed= T){
  if(directed){print("Assuming interactions are directed")}
  if(!directed){print("Assuming interactions are undirected")}
  if (directed){
    actual.rates <- 
      df %>% 
      mutate(bin= bin) %>% 
      group_by(bin, dyad) %>% 
      summarize(duration= sum(duration)) %>% 
      separate(dyad, into=c("actor", "receiver"), sep="_")
    possible.rates <- 
      expand_grid(actor= ids, receiver= ids, bin= unique(bin)) %>% 
      filter(actor!=receiver) %>% 
      mutate(duration=0)
    rates <- 
      rbind.data.frame(actual.rates, possible.rates) %>% 
      group_by(bin, actor, receiver) %>% 
      summarize(duration= sum(duration)) %>% 
      ungroup() %>% 
      mutate(dyad= paste(actor, receiver, sep="_"))
    rates
  }else{
    if (!directed){
      actual.rates <- 
        df %>% 
        mutate(bin= bin) %>% 
        separate(dyad, into=c("id1", "id2"), sep="_") %>% 
        mutate(dyad= if_else(id1<id2, paste(id1,id2, sep="_"), paste(id2,id1, sep="_"))) %>% 
        group_by(bin, dyad) %>% 
        summarize(duration= sum(duration)) %>% 
        separate(dyad, into=c("id1", "id2"), sep="_")
      possible.rates <- 
        expand_grid(id1= ids, id2= ids, bin= unique(bin)) %>% 
        filter(id1!=id2) %>%
        filter(id1<id2) %>% 
        mutate(duration=0)
      rates <- 
        rbind.data.frame(actual.rates, possible.rates) %>% 
        group_by(bin, id1, id2) %>% 
        summarize(duration= sum(duration)) %>% 
        ungroup() %>% 
        mutate(dyad= paste(id1, id2, sep="_"))
      rates
    }
  }
  return(rates)
}

# plot permutation test results
hist_perm <- function(exp=exp, obs=obs, perms=perms){
  ggplot()+
    geom_histogram(aes(x=exp), color="black",fill="light blue")+
    xlim(min= min(c(exp,obs)), max= max(c(exp,obs)))+
    geom_vline(aes(xintercept=obs), color="red", size=1)+
    xlab("expected values from null model")+
    labs(subtitle= paste('obs = ',round(obs, digits=2),", one-sided p = ", mean(exp>=obs),", permutations=",perms, sep=""))
}

# get standard error of the mean
se <- function(x=x){sd(x, na.rm=T)/sqrt(sum(!is.na(x)))}

# define treatment time period-----
# all bats were injected by 2 pm (1400 h) on April 25
# LPS effect begins 3 hours later (1700 h)
# captive data shows LPS effects at 3 hours and 6 hours post-injection
# LPS effects could last longer

# start sampling behavior at 5 pm on April 25 
treatment_start <- as.POSIXct("2018-04-25 17:00:00", tz = "CST6CDT")

# stop sampling before midnight (before bats are likely to forage)
LPS_duration_hours <- 6
treatment_stop <- treatment_start + LPS_duration_hours*60*60

# get first post-treatment time period  (same time of day 24 hours later)
post24_start <- treatment_start + 60*60*24
post24_stop <- treatment_stop + 60*60*24

# get second post-treatment time period  (same time of day 48 hours later)
post48_start <- treatment_start + 60*60*48
post48_stop <- treatment_stop + 60*60*48


# clean data----

# get raw meeting data (unfused)
Belize <- 
  read_delim("../data/Belize_pre-fusion.csv", delim = ";") 

# exclude missing bats and dropped sensors ----
# bat ID 23 also seems to have left 1 day early
excluded_bats <- c(10,14,26,41) 
BatIDs <- unique(c(unique(Belize$SenderID),unique(Belize$EncounteredID)))
BatIDs <- BatIDs[! BatIDs %in% excluded_bats]

# get bat attributes
bats <- 
  read.csv("../data/Belize_tracked_bats02.csv", stringsAsFactors = F) %>% 
  mutate(sex = "female") %>%
  filter(sensor_node %in% BatIDs)

# define association using RSSI  ----
thresholds <- quantile(Belize$RSSI, probs = c(0, 25, 50, 75, 80, 85, 90, 95, 97.5, 99)/100)
thresholds
RSSI_threshold = thresholds[6] #set RSSI threshold 
RSSI_threshold 
#threshold for 2019 Current Biology paper was -26dmb at 90%; here -27dbm at 85%

# plot and print RSSI threshold-----
rssi.plot <- 
  Belize %>% 
  ggplot(aes(x=RSSI))+
  geom_histogram(binwidth=1, fill= "grey", color= 'black')+
  geom_vline(xintercept= RSSI_threshold)+
  ggtitle("proximity sensor signal strength threshold", 
          subtitle= paste(names(RSSI_threshold), "threshold at", RSSI_threshold, "dBm"))+
  xlab("Received Signal Strength Indicator (RSSI)")+
  theme_cowplot()
ggsave("RSSI_plot.pdf", width= 6, height= 3, units="in", dpi=1200)

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

# insert hour breaks
# convert start and end times to interval
df$interval <- interval(df$StartOfMeeting, df$EndOfMeeting)
# label "events" sequentially
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
  data.frame(event, hour = hours, duration = dur)
}

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
# set start of study to 3pm on April 25th
df <- 
  df2 %>%
  mutate(bat1 = sub( "_.*$", "", dyad ), bat2 = sub('.*_', '', dyad)) %>%
  mutate(StartBelizeTime = force_tz(StartOfMeeting, tzone = "CST6CDT")) %>%
  mutate(EndBelizeTime = force_tz(EndOfMeeting, tzone = "CST6CDT")) %>%
  mutate(StartBelizeTime= StartBelizeTime - hours(8), EndBelizeTime=EndBelizeTime - hours(8)) %>%
  filter(StartBelizeTime >= as.POSIXct("2018-04-25 15:00:00", tz = "CST6CDT")) %>% 
  select(StartBelizeTime,bat1,bat2,duration, RSSI)

# assign treatment for each bat into df
df$treatment_bat1 <- bats$treatment[match(df$bat1, bats$sensor_node)]
df$treatment_bat2 <- bats$treatment[match(df$bat2, bats$sensor_node)]

# remove other dataframe
rm(df2)

# label treatment and post-treatment periods 
# label dyad types
d <- 
  df %>% 
  mutate(datetime= as.POSIXct(StartBelizeTime,tz = "CST6CDT")) %>% 
  mutate(treatment_period= datetime >= treatment_start & datetime < treatment_stop) %>% 
  mutate(post24_period= datetime >= post24_start & datetime < post24_stop) %>% 
  mutate(post48_period= datetime >= post48_start & datetime < post48_stop) %>% 
  mutate(dyad_type= case_when(
    treatment_bat1== 'LPS' & treatment_bat2== 'LPS' ~ 'sick-sick',
    treatment_bat1== 'PBS' & treatment_bat2== 'PBS' ~ 'control-control',
    treatment_bat1!= treatment_bat2 ~ 'sick-control',
    TRUE ~ 'NA')) %>% 
  mutate(hour= substring(datetime, 1,13))

# get treated bats
treated.bats <- 
  d %>% 
  mutate(treated= ifelse(treatment_bat1=="LPS", bat1, 
                         ifelse(treatment_bat2=="LPS", bat2, NA))) %>% 
  filter(!is.na(treated)) %>% 
  pull(treated) %>% 
  unique()

# number of treated bats = 16
length(treated.bats)

# number of untreated bats = 15
d %>% 
  mutate(untreated= ifelse(treatment_bat1!="LPS", bat1, 
                         ifelse(treatment_bat2!="LPS", bat2, NA))) %>% 
  filter(!is.na(untreated)) %>% 
  pull(untreated) %>% 
  unique() %>% 
  length()

# get hourly associations----
h <- 
  d %>% 
  mutate(dyad= paste(bat1, bat2, sep="_")) %>%
  group_by(hour, dyad) %>% 
  summarize(duration= sum(duration)) %>% 
  ungroup() %>% 
  events_to_rates(df=., bin= .$hour, ids= BatIDs, directed=F) %>% 
  select(hour= bin, dyad, duration)

# make one-hour networks----------
# make multilayer network where every layer is an hour
# net.list is a list of networks (one per hour)
all.hours <- unique(h$hour)
net.list <- list() 
for (i in 1:length(all.hours)){
  layer <-  sort(all.hours)[i]
  net <- 
    h %>% 
    filter(hour==layer) %>%
    select(dyad, duration) %>% 
    a_b_edgelist_to_matrix() %>% ###
    graph_from_adjacency_matrix("undirected", weighted=T, diag=F) ###
  net.list[[i]] <- net
}

# get mean centrality per treatment within each hour (as a list)
d.list <- list()
for (i in 1:length(all.hours)){
 net <- net.list[[i]]
 d.list[[i]] <-
   tibble(hour= as.POSIXct(all.hours[i], format= "%Y-%m-%d %H", tz="CST6CDT"),
          bat= names(degree(net)),
          degree= degree(net),
          strength= strength(net),
          eigenvector= eigen_centrality(net)$vector) %>%
   mutate(treated= bat %in% treated.bats)
}

# d2 is the mean centrality per treatment group and per hour ----
d2 <- 
  bind_rows(d.list) %>% 
  group_by(hour, treated) %>% 
  summarize(degree= mean(degree),
            strength= mean(strength),
            eigenvector= mean(eigenvector)) %>% 
  pivot_longer(cols= degree:eigenvector, 
               names_to = "centrality", values_to = "value") %>% 
  filter(hour <= as.POSIXct("2018-04-28 01:00:00", tz = "CST6CDT")) %>% 
  mutate(period= case_when(
    hour >= treatment_start & hour < treatment_stop ~ 'treatment',
    hour >= post24_start & hour < post24_stop ~ '24 hours later',
    hour >= post48_start & hour < post48_stop ~ '48 hours later',
    TRUE ~"none")) %>% 
  ungroup()

# get the corresponding standard errors
d2.se <- 
  bind_rows(d.list) %>% 
  group_by(hour, treated) %>% 
  summarize(degree= se(degree),
            strength= se(strength),
            eigenvector= se(eigenvector)) %>% 
  pivot_longer(cols= degree:eigenvector, 
               names_to = "centrality", values_to = "se") %>% 
  filter(hour <= as.POSIXct("2018-04-28 01:00:00", tz = "CST6CDT")) %>% 
  ungroup()

# add standard errors 
d2 <- full_join(d2, d2.se)  
rm(d2.se)

# d3 is the mean centrality per bat and per hour----
d3 <- 
  bind_rows(d.list) %>% 
  mutate(period= case_when(
    hour >= treatment_start & hour < treatment_stop ~ 'treatment',
    hour >= post24_start & hour < post24_stop ~ '24 hours later',
    hour >= post48_start & hour < post48_stop ~ '48 hours later',
    TRUE ~"none")) %>% 
  filter(period!="none") %>% 
  mutate(treatment= ifelse(treated, "sick", "control")) 

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

# function to plot networks
network_plot <- function(net, title=""){
  set.seed(123)
  V(net)$treated <- ifelse(V(net)$name %in% treated.bats, 'dark blue', "light blue")
  layout <- layout_with_gem(net)
  png(filename=paste(title,".png", sep=""), width = 1200, height = 1200)
  plot(net, 
       edge.color="slate grey",
       vertex.shape= 'sphere',
       vertex.color=V(net)$treated, 
       vertex.label="",
       vertex.size=8,
       edge.width=log10(E(net)$weight),
       layout=layout)
  title(title,cex.main=5,col.main="black")
  dev.off()
}

network_plot(network.treatment, 'social network during treatment period')
network_plot(network.post24, "social network after 24 hours")
network_plot(network.post48, 'social network after 48 hours')

# get network centrality by treatment period-----
centrality.treatment <- 
  tibble(bat= names(degree(network.treatment)), 
         degree= degree(network.treatment),
         strength= strength(network.treatment),
         eigenvector= eigen_centrality(network.treatment)$vector,
         period= "treatment")
centrality.post24 <- 
  tibble(bat= names(degree(network.post24)), 
         degree= degree(network.post24),
         strength= strength(network.post24),
         eigenvector= eigen_centrality(network.post24)$vector,
         period= "24 hours later") 
centrality.post48 <- 
  tibble(bat= names(degree(network.post48)), 
         degree= degree(network.post48),
         strength= strength(network.post48),
         eigenvector= eigen_centrality(network.post48)$vector,
         period= "48 hours later") 
d.period <- 
  rbind(centrality.treatment, centrality.post24, centrality.post48) %>% 
  mutate(treated= bat %in% treated.bats) %>% 
  mutate(treatment= ifelse(treated, "sick", "control"))
  
# get mean centrality values
means.degree <- 
  d.period %>% 
  mutate(group= paste(treated, period, sep="_")) %>% 
  boot_ci2(y=.$degree, x=.$group) %>% 
  separate(effect, into=c('treated', 'period'), sep="_")
means.str<- 
  d.period %>% 
  mutate(group= paste(treated, period, sep="_")) %>% 
  boot_ci2(y=.$strength/3600, x=.$group) %>% 
  separate(effect, into=c('treated', 'period'), sep="_")
means.ec<- 
  d.period %>% 
  mutate(group= paste(treated, period, sep="_")) %>% 
  boot_ci2(y=.$eigenvector, x=.$group) %>% 
  separate(effect, into=c('treated', 'period'), sep="_")

# create function to plot effect sizes
plot_effect <- function (df){
  df %>% 
    mutate(period2= factor(period, levels= c('treatment', '24 hours later', '48 hours later'))) %>% 
    mutate(treatment= ifelse(treated, 'sick', 'control')) %>% 
    ggplot(aes(x=treatment, y=mean, color=treatment, shape=treatment))+
    facet_wrap(~period2)+
    geom_point(size=3)+
    geom_errorbar(aes(ymin=low, ymax=high, width=.1), size=1)+
    theme_bw()+
    theme(legend.position= 'none',
          axis.title.y=element_blank(),
          axis.title.x=element_blank())+
    scale_color_manual(values=c("light blue", "dark blue"))+
    scale_shape_manual(values=c("circle", 'triangle'))
}

# plot centrality by hour-----
ch.plot <- 
  d2 %>% 
  mutate(value= ifelse(centrality=='strength', value/3600, value)) %>% 
  mutate(se= ifelse(centrality=='strength', se/3600, se)) %>% 
  mutate(centrality= paste(centrality, "centrality")) %>% 
  mutate(centrality= factor(centrality, 
                            levels= c("degree centrality", 'strength centrality', 'eigenvector centrality'))) %>% 
  mutate(injection= ifelse(treated, "sick", "control")) %>%
  mutate(high= value+se, low= value-se) %>% 
  ggplot(aes(x=hour, y=value, color=injection, shape=injection))+
  #facet_wrap(~centrality, scales = "free_y", nrow = 3, labeller = "label_value")+
  facet_wrap(~centrality, scales = "free_y", nrow = 3, 
             strip.position = "left",
             labeller = as_labeller(c(`degree centrality` = "degree centrality",
                                      `strength centrality` = "strength centrality",
                                      `eigenvector centrality` = "eigenvector centrality")))+
  geom_point()+
  geom_errorbar(aes(ymin=low, ymax= high), width=0.1)+
  geom_line()+
  geom_rect(aes(xmin=treatment_stop, xmax=post24_start, ymin=0, ymax=Inf), fill='white', color=NA, alpha=0.02)+
  geom_rect(aes(xmin=post24_stop, xmax=post48_start, ymin=0, ymax=Inf), fill='white', color=NA, alpha=0.02)+
  geom_rect(aes(xmin=post48_stop, xmax=max(d2$hour)+1500, ymin=0, ymax=Inf), fill='white', color=NA, alpha=0.02)+
  geom_vline(xintercept = treatment_start, color= "black")+
  geom_vline(xintercept=treatment_stop, color= "black")+
  geom_vline(xintercept = post24_start, linetype= "solid")+
  geom_vline(xintercept=post24_stop, linetype= "solid")+
  geom_vline(xintercept = post48_start, linetype= "solid")+
  geom_vline(xintercept=post48_stop, linetype= "solid")+
  theme_bw()+
  ylab(NULL) +
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        axis.title.x=element_blank(),
        legend.position = c(0.18, 0.72),
        legend.margin=margin(c(0.5,0.5,0.5,0.5)),
        legend.justification = "left",
        legend.title = element_blank(),
        legend.direction = "horizontal",
        legend.background = element_rect(fill='white', 
                                         size=0.5, linetype="solid"))+
  scale_x_datetime(breaks=date_breaks("6 hour"), labels=date_format("%H:%M", tz = "CST6CDT")) + 
  scale_color_manual(values=c("light blue", "dark blue"))+
  scale_shape_manual(values=c("circle", 'triangle')) 
ch.plot

# add panels of effect sizes
p1 <- means.degree %>% plot_effect()  + theme(axis.title.x=element_blank())
p2 <- means.str %>% plot_effect()  + theme(axis.title.x=element_blank())
p3 <- means.ec  %>% plot_effect()  + theme(axis.title.x=element_blank())
p4 <- plot_grid(p1,p2,p3, ncol=1, align = 'hv', rel_heights = c(1,1,1))

# plot and print to PDF
plot_grid(ch.plot, p4, ncol = 2, axis= 't',align = "v", rel_widths = c(2,1))
ggsave("centrality_by_hour.pdf", width= 11, height= 5, units="in", dpi=1200)

# fit models----
# general linear mixed effect model of LPS effect on degree
# response = degree centrality (one network per time period)
# fixed effect = time period, treatment, and their interaction
# random effect = bat

# first fit parametric GLM model

# add day as variable
data <- 
  d.period %>% 
  mutate(day= case_when(
    period=="treatment" ~ 1,
    period=="24 hours later" ~ 2,
    period=="48 hours later" ~ 3))

# get standardized effect sizes
# get degree effects
fit1 <- summary(lmer(scale(degree)~ day*treated+(1|bat), data= data))
fit2 <- summary(lm(scale(degree)~ treated, data= data[which(data$day==1),])) #effect during treatment
fit2b <- summary(lm(scale(degree)~ treated, data= data[which(data$day==2),])) 
fit3 <- summary(lm(scale(degree)~ treated, data= data[which(data$day==3),])) #effect outside treatment

# get strength effects
fit4 <- summary(lmer(scale(strength)~ day*treated+(1|bat), data= data))
fit5 <- summary(lm(scale(strength)~ treated, data= data[which(data$day==1),])) 
fit5b <- summary(lm(scale(strength)~ treated, data= data[which(data$day==2),])) 
fit6 <- summary(lm(scale(strength)~ treated, data= data[which(data$day==3),])) 

# get eigenvector effects
fit7 <- summary(lmer(scale(eigenvector)~ day*treated+(1|bat), data= data))
fit8 <- summary(lm(scale(eigenvector)~ treated, data= data[which(data$day==1),])) 
fit8b <- summary(lm(scale(eigenvector)~ treated, data= data[which(data$day==2),])) 
fit9 <- summary(lm(scale(eigenvector)~ treated, data= data[which(data$day==3),])) 

# get biologically meaningful effect sizes

# get change in number of bats
summary(lm((degree)~ treated, data= data[which(data$day==1),]))
# 4 fewer bats

# get change in time spent per bat
summary(lm((strength/30)~ treated, data= data[which(data$day==1),]))
# 1510 fewer seconds/six hour

# get observed slopes 
# degree
obs1 <- fit1$coefficients[4,1] # interaction effect treatedTrue:period
obs1.1 <- fit1$coefficients[3,1] # LPS effect controlling for period 
obs2 <- fit2$coefficients[2,1] # LPS effect during treatment
obs2.1 <- fit2b$coefficients[2,1] # LPS effect on day 2
obs3 <- fit3$coefficients[2,1] # LPS effect post treatment
# strength
obs4 <- fit4$coefficients[4,1] # interaction effect treatedTrue:period
obs4.1 <- fit4$coefficients[3,1] # LPS effect controlling for period 
obs5 <- fit5$coefficients[2,1] # LPS effect during treatment
obs5.1 <- fit5b$coefficients[2,1] # LPS effect on day 2
obs6 <- fit6$coefficients[2,1] # LPS effect post treatment
# eigenvector
obs7 <- fit7$coefficients[4,1] # interaction effect treatedTrue:period
obs7.1 <- fit7$coefficients[3,1] # LPS effect controlling for period 
obs8 <- fit8$coefficients[2,1] # LPS effect during treatment
obs8.1 <- fit8b$coefficients[2,1] # LPS effect on day 2
obs9 <- fit9$coefficients[2,1] # LPS effect post treatment

# permutation test to obtain non-parametric p-values

if (TRUE){
  # get observed coefficients
  
  # get expected
  exp1 <- rep(NA, perms)
  exp1.1 <- rep(NA, perms)
  exp2 <- rep(NA, perms)
  exp2.1 <- rep(NA, perms)
  exp3 <- rep(NA, perms)
  exp4 <- rep(NA, perms)
  exp4.1 <- rep(NA, perms)
  exp5 <- rep(NA, perms)
  exp5.1 <- rep(NA, perms)
  exp6 <- rep(NA, perms)
  exp7 <- rep(NA, perms)
  exp7.1 <- rep(NA, perms)
  exp8 <- rep(NA, perms)
  exp8.1 <- rep(NA, perms)
  exp9 <- rep(NA, perms)
  
  start <- Sys.time()
  for (i in 1:perms){
    # swap which bats are treated
    random.treated.bats <- 
      d.period %>% 
      group_by(bat) %>% 
      summarize(treated=sum(treated)==3) %>% 
      mutate(random.treated= sample(treated)) %>% 
      filter(random.treated) %>% 
      pull(bat)
    
    # refit models with random treated bats
    rdata <- 
      d.period %>% 
      mutate(day= case_when(
        period=="treatment" ~ 1,
        period=="24 hours later" ~ 2,
        period=="48 hours later" ~ 3)) %>% 
      # relabel treated bats
      mutate(treated= bat %in% random.treated.bats) 
    
    # get degree effects
    rfit1 <- summary(lmer(scale(degree)~ day*treated+(1|bat), data= rdata))
    rfit2 <- summary(lm(scale(degree)~ treated, data= rdata[which(rdata$day==1),])) #effect during treatment
    rfit2b <- summary(lm(scale(degree)~ treated, data= rdata[which(rdata$day==2),])) # effect on day 2
    rfit3 <- summary(lm(scale(degree)~ treated, data= rdata[which(rdata$day==3),])) #effect outside treatment
    
    # get strength effects
    rfit4 <- summary(lmer(scale(strength)~ day*treated+(1|bat), data= rdata))
    rfit5 <- summary(lm(scale(strength)~ treated, data= rdata[which(rdata$day==1),])) 
    rfit5b <- summary(lm(scale(strength)~ treated, data= rdata[which(rdata$day==2),])) # effect on day 2
    rfit6 <- summary(lm(scale(strength)~ treated, data= rdata[which(rdata$day==3),])) 
    
    # get eigenvector effects
    rfit7 <- summary(lmer(scale(eigenvector)~ day*treated+(1|bat), data= rdata))
    rfit8 <- summary(lm(scale(eigenvector)~ treated, data= rdata[which(rdata$day==1),])) 
    rfit8b <- summary(lm(scale(eigenvector)~ treated, data= rdata[which(rdata$day==2),])) # effect on day 2
    rfit9 <- summary(lm(scale(eigenvector)~ treated, data= rdata[which(rdata$day==3),])) 
    
    # save coefficients
    exp1[i]<- rfit1$coefficients[4,1]    
    exp1.1[i]<- rfit1$coefficients[3,1] 
    exp2[i] <- rfit2$coefficients[2,1]
    exp2.1[i] <- rfit2b$coefficients[2,1]
    exp3[i] <- rfit3$coefficients[2,1]
    exp4[i]<- rfit4$coefficients[4,1]    
    exp4.1[i]<- rfit4$coefficients[3,1] 
    exp5[i] <- rfit5$coefficients[2,1]
    exp5.1[i] <- rfit5b$coefficients[2,1]
    exp6[i] <- rfit6$coefficients[2,1]
    exp7[i]<- rfit7$coefficients[4,1]    
    exp7.1[i]<- rfit7$coefficients[3,1] 
    exp8[i] <- rfit8$coefficients[2,1]
    exp8.1[i] <- rfit8b$coefficients[2,1]
    exp9[i] <- rfit9$coefficients[2,1]
    if(i%%10==0) print(paste(i, "of", perms))
  }
  save(list= c('exp1', 'exp1.1','exp2','exp2.1', 'exp3','exp4', 'exp4.1','exp5','exp5.1', 'exp6','exp7', 'exp7.1','exp8', 'exp8.1', 'exp9'), 
       file= "perm_results.Rdata") ###
}else{ load('perm_results.Rdata') }
# end permutation test
speed <- Sys.time()- start
speed


# get perm test results----
# model coefficients (slopes) and p-values
# degree
# interaction effect (does LPS effect differ by treatment period?)
t1 <- hist_perm(exp1, obs1, perms)+labs(title='does LPS effect on degree differ by treatment period?')
# does LPS have an effect (controlling for period)?
t2 <- hist_perm(exp1.1, obs1.1, perms)+labs(title='does LPS affect degree?')
# does LPS have effect during treatment period?
t3 <- hist_perm(exp2, obs2, perms)+labs(title='does LPS affect degree in treatment period?')
t3.1 <- hist_perm(exp2.1, obs2.1, perms)+labs(title='does LPS affect degree on day 2?')
# does LPS have effect outside the treatment period?
t4 <- hist_perm(exp3, obs3, perms)+labs(title='does LPS affect degree in post-treatment period?')

# strength
# interaction effect (does LPS effect differ by treatment period?)
t5 <- hist_perm(exp4, obs4, perms)+labs(title='does LPS effect on strength differ by treatment period?')
# does LPS have an effect (controlling for period)?
t6 <- hist_perm(exp4.1, obs4.1, perms)+labs(title='does LPS affect strength?')
# does LPS have effect during treatment period?
t7 <- hist_perm(exp5, obs5, perms)+labs(title='does LPS affect strength in treatment period?')
t7.1 <- hist_perm(exp5.1, obs5.1, perms)+labs(title='does LPS affect strength on day 2?')
# does LPS have effect outside the treatment period?
t8 <- hist_perm(exp6, obs6, perms)+labs(title='does LPS affect strength in post-treatment period?')

# eigenvector
# interaction effect (does LPS effect differ by treatment period?)
t9 <- hist_perm(exp7, obs7, perms)+labs(title='does LPS effect on eigenvector centrality differ by treatment period?')
# does LPS have an effect (controlling for period)?
t10 <- hist_perm(exp7.1, obs7.1, perms)+labs(title='does LPS affect eigenvector centrality?')
# does LPS have effect during treatment period?
t11 <- hist_perm(exp8, obs8, perms)+labs(title='does LPS affect eigenvector centrality in treatment period?')
t11.1 <- hist_perm(exp8.1, obs8.1, perms)+labs(title='does LPS affect eigenvector centrality in treatment period?')
# does LPS have effect outside the treatment period?
t12 <- hist_perm(exp9, obs9, perms)+labs(title='does LPS affect eigenvector centrality in post-treatment period?')

# combine plots
(perm.test <- t1+t2+t3+t3.1+t4+t5+t6+t7+t7.1+t8+t9+t10+t11+t11.1+t12+plot_layout(ncol=3))
ggsave("perm_tests.pdf", width=18, height= 15, units="in", dpi=1200)

# make results table-----
results <- 
  data.frame(response= rep(c('degree', 'strength', 'eigenvector'), times=1, each=5),
             fixed_effect= rep(c('interaction', 'treatment', 'treatment.within.period','treatment.day2', 'treatment.within.post'), times=3),
             coefficient= c(obs1, obs1.1, obs2, obs2.1, obs3, obs4, obs4.1, obs5, obs5.1, obs6, obs7, obs7.1, obs8, obs8.1, obs9),
             pvalue= c(mean(exp1>=obs1),
                       mean(exp1.1>=obs1.1),
                       mean(exp2>=obs2),
                       mean(exp2.1>=obs2.1),
                       mean(exp3>=obs3),
                       mean(exp4>=obs4),
                       mean(exp4.1>=obs4.1),
                       mean(exp5>=obs5),
                       mean(exp5.1>=obs5.1),
                       mean(exp6>=obs6),
                       mean(exp7>=obs7),
                       mean(exp7.1>=obs7.1),
                       mean(exp8>=obs8),
                       mean(exp8.1>=obs8.1),
                       mean(exp9>=obs9))) %>% 
  mutate(pvalue= if_else(coefficient<0, (1-pvalue),pvalue)) %>% 
  mutate(pvalue2= c(mean(abs(exp1) >= abs(obs1)),
                      mean(abs(exp1.1) >= abs(obs1.1)),
                      mean(abs(exp2) >= abs(obs2)),
                      mean(abs(exp2.1) >= abs(obs2.1)),
                      mean(abs(exp3) >= abs(obs3)),
                      mean(abs(exp4) >= abs(obs4)),
                      mean(abs(exp4.1) >= abs(obs4.1)),
                      mean(abs(exp5) >= abs(obs5)),
                      mean(abs(exp5.1) >= abs(obs5.1)),
                      mean(abs(exp6) >= abs(obs6)),
                      mean(abs(exp7) >= abs(obs7)),
                      mean(abs(exp7.1) >= abs(obs7.1)),
                      mean(abs(exp8) >= abs(obs8)),
                      mean(abs(exp8.1) >= abs(obs8.1)),
                      mean(abs(exp9) >= abs(obs9))))
results


# Alternate analysis----
# this is an alternate analysis requested by reviewers
if(TRUE){
  # fit models on HOURLY DATA----
  # general linear mixed effect model of LPS effect on degree 
  # response = degree centrality (one network per hour)
  # fixed effect = time period, treatment, and their interaction
  # random effect = bat, hour
  time1=Sys.time()
  perms2=1000
  
  # first fit parametric GLM model
  
  # add day as variable
  data <- 
    d3 %>% 
    mutate(day= case_when(
      period=="treatment" ~ 1,
      period=="24 hours later" ~ 2,
      period=="48 hours later" ~ 3)) %>% 
    # convert datetime to hours (1-6)
    mutate(hour= as.numeric(substring(hour, 12,13)) -16 ) %>%  
    # convert to string
    mutate(hour=as.character(hour))
  
  # get degree effects
  fit1 <- summary(lmer(scale(degree)~ day*treated+(1|hour) +(1|bat), data= data))
  fit2 <- summary(lmer(scale(degree)~ treated+(1|hour), data= data[which(data$day==1),])) #effect during treatment
  fit2b <- summary(lmer(scale(degree)~ treated+(1|hour), data= data[which(data$day==2),])) 
  fit3 <- summary(lmer(scale(degree)~ treated+(1|hour), data= data[which(data$day==3),])) #effect outside treatment
  
  # get strength effects
  fit4 <- summary(lmer(scale(strength)~ day*treated+(1|hour)+(1|bat), data= data))
  fit5 <- summary(lmer(scale(strength)~ treated+(1|hour), data= data[which(data$day==1),])) 
  fit5b <- summary(lmer(scale(strength)~ treated+(1|hour), data= data[which(data$day==2),])) 
  fit6 <- summary(lmer(scale(strength)~ treated+(1|hour), data= data[which(data$day==3),])) 
  
  # get eigenvector effects
  fit7 <- summary(lmer(scale(eigenvector)~ day*treated+(1|hour)+(1|bat), data= data))
  fit8 <- summary(lmer(scale(eigenvector)~ treated+(1|hour), data= data[which(data$day==1),])) 
  fit8b <- summary(lmer(scale(eigenvector)~ treated+(1|hour), data= data[which(data$day==2),])) 
  fit9 <- summary(lmer(scale(eigenvector)~ treated+(1|hour), data= data[which(data$day==3),])) 
  
  # get observed slopes 
  # degree
  obs1 <- fit1$coefficients[4,1] # interaction effect treatedTrue:period
  obs1.1 <- fit1$coefficients[3,1] # LPS effect controlling for period 
  obs2 <- fit2$coefficients[2,1] # LPS effect during treatment
  obs2.1 <- fit2b$coefficients[2,1] # LPS effect on day 2
  obs3 <- fit3$coefficients[2,1] # LPS effect post treatment
  # strength
  obs4 <- fit4$coefficients[4,1] # interaction effect treatedTrue:period
  obs4.1 <- fit4$coefficients[3,1] # LPS effect controlling for period 
  obs5 <- fit5$coefficients[2,1] # LPS effect during treatment
  obs5.1 <- fit5b$coefficients[2,1] # LPS effect on day 2
  obs6 <- fit6$coefficients[2,1] # LPS effect post treatment
  # eigenvector
  obs7 <- fit7$coefficients[4,1] # interaction effect treatedTrue:period
  obs7.1 <- fit7$coefficients[3,1] # LPS effect controlling for period 
  obs8 <- fit8$coefficients[2,1] # LPS effect during treatment
  obs8.1 <- fit8b$coefficients[2,1] # LPS effect on day 2
  obs9 <- fit9$coefficients[2,1] # LPS effect post treatment
  
  # permutation test to obtain non-parametric p-values
  
  if (TRUE){
    # get observed coefficients
    
    # get expected
    exp1 <- rep(NA, perms2)
    exp1.1 <- rep(NA, perms2)
    exp2 <- rep(NA, perms2)
    exp2.1 <- rep(NA, perms2)
    exp3 <- rep(NA, perms2)
    exp4 <- rep(NA, perms2)
    exp4.1 <- rep(NA, perms2)
    exp5 <- rep(NA, perms2)
    exp5.1 <- rep(NA, perms2)
    exp6 <- rep(NA, perms2)
    exp7 <- rep(NA, perms2)
    exp7.1 <- rep(NA, perms2)
    exp8 <- rep(NA, perms2)
    exp8.1 <- rep(NA, perms2)
    exp9 <- rep(NA, perms2)
    
    start <- Sys.time()
    for (i in 1:perms2){
      # swap which bats are treated
      random.treated.bats <- 
        d.period %>% 
        group_by(bat) %>% 
        summarize(treated=sum(treated)==3) %>% 
        mutate(random.treated= sample(treated)) %>% 
        filter(random.treated) %>% 
        pull(bat)
      
      # refit models with random treated bats
      rdata <- 
        d3 %>% 
        mutate(day= case_when(
          period=="treatment" ~ 1,
          period=="24 hours later" ~ 2,
          period=="48 hours later" ~ 3)) %>% 
        # convert datetime to hours (1-6)
        mutate(hour= as.numeric(substring(hour, 12,13)) -16 ) %>%  
        # convert to string
        mutate(hour=as.character(hour)) %>% 
        # relabel treated bats
        mutate(treated= bat %in% random.treated.bats) 
      
      # get degree effects
      rfit1 <- summary(lmer(scale(degree)~ day*treated+(1|hour), data= rdata))
      rfit2 <- summary(lmer(scale(degree)~ treated+(1|hour), data= rdata[which(rdata$day==1),])) #effect during treatment
      rfit2b <- summary(lmer(scale(degree)~ treated+(1|hour), data= rdata[which(rdata$day==2),])) # effect on day 2
      rfit3 <- summary(lmer(scale(degree)~ treated+(1|hour), data= rdata[which(rdata$day==3),])) #effect outside treatment
      
      # get strength effects
      rfit4 <- summary(lmer(scale(strength)~ day*treated+(1|hour), data= rdata))
      rfit5 <- summary(lmer(scale(strength)~ treated +(1|hour), data= rdata[which(rdata$day==1),])) 
      rfit5b <- summary(lmer(scale(strength)~ treated +(1|hour), data= rdata[which(rdata$day==2),])) # effect on day 2
      rfit6 <- summary(lmer(scale(strength)~ treated+(1|hour), data= rdata[which(rdata$day==3),])) 
      
      # get eigenvector effects
      rfit7 <- summary(lmer(scale(eigenvector)~ day*treated+(1|hour), data= rdata))
      rfit8 <- summary(lmer(scale(eigenvector)~ treated+(1|hour), data= rdata[which(rdata$day==1),])) 
      rfit8b <- summary(lmer(scale(eigenvector)~ treated+(1|hour), data= rdata[which(rdata$day==2),])) # effect on day 2
      rfit9 <- summary(lmer(scale(eigenvector)~ treated+(1|hour), data= rdata[which(rdata$day==3),])) 
      
      # save coefficients
      exp1[i]<- rfit1$coefficients[4,1]    
      exp1.1[i]<- rfit1$coefficients[3,1] 
      exp2[i] <- rfit2$coefficients[2,1]
      exp2.1[i] <- rfit2b$coefficients[2,1]
      exp3[i] <- rfit3$coefficients[2,1]
      exp4[i]<- rfit4$coefficients[4,1]    
      exp4.1[i]<- rfit4$coefficients[3,1] 
      exp5[i] <- rfit5$coefficients[2,1]
      exp5.1[i] <- rfit5b$coefficients[2,1]
      exp6[i] <- rfit6$coefficients[2,1]
      exp7[i]<- rfit7$coefficients[4,1]    
      exp7.1[i]<- rfit7$coefficients[3,1] 
      exp8[i] <- rfit8$coefficients[2,1]
      exp8.1[i] <- rfit8b$coefficients[2,1]
      exp9[i] <- rfit9$coefficients[2,1]
      if(i%%10==0) print(paste(i, "of", perms2))
    }
    save(list= c('exp1', 'exp1.1','exp2','exp2.1', 'exp3','exp4', 'exp4.1','exp5','exp5.1', 'exp6','exp7', 'exp7.1','exp8', 'exp8.1', 'exp9'), 
         file= "perm_results02.Rdata") ###
  }else{ load('perm_results02.Rdata') }
  # end permutation test
  
  
  # get perm test results
  # model coefficients (slopes) and p-values
  # degree
  # interaction effect (does LPS effect differ by treatment period?)
  t1 <- hist_perm(exp1, obs1, perms2)+labs(title='does LPS effect on degree differ by treatment period?')
  # does LPS have an effect (controlling for period)?
  t2 <- hist_perm(exp1.1, obs1.1, perms2)+labs(title='does LPS affect degree?')
  # does LPS have effect during treatment period?
  t3 <- hist_perm(exp2, obs2, perms2)+labs(title='does LPS affect degree in treatment period?')
  t3.1 <- hist_perm(exp2.1, obs2.1, perms2)+labs(title='does LPS affect degree on day 2?')
  # does LPS have effect outside the treatment period?
  t4 <- hist_perm(exp3, obs3, perms2)+labs(title='does LPS affect degree in post-treatment period?')
  
  # strength
  # interaction effect (does LPS effect differ by treatment period?)
  t5 <- hist_perm(exp4, obs4, perms2)+labs(title='does LPS effect on strength differ by treatment period?')
  # does LPS have an effect (controlling for period)?
  t6 <- hist_perm(exp4.1, obs4.1, perms2)+labs(title='does LPS affect strength?')
  # does LPS have effect during treatment period?
  t7 <- hist_perm(exp5, obs5, perms2)+labs(title='does LPS affect strength in treatment period?')
  t7.1 <- hist_perm(exp5.1, obs5.1, perms2)+labs(title='does LPS affect strength on day 2?')
  # does LPS have effect outside the treatment period?
  t8 <- hist_perm(exp6, obs6, perms2)+labs(title='does LPS affect strength in post-treatment period?')
  
  # eigenvector
  # interaction effect (does LPS effect differ by treatment period?)
  t9 <- hist_perm(exp7, obs7, perms2)+labs(title='does LPS effect on eigenvector centrality differ by treatment period?')
  # does LPS have an effect (controlling for period)?
  t10 <- hist_perm(exp7.1, obs7.1, perms2)+labs(title='does LPS affect eigenvector centrality?')
  # does LPS have effect during treatment period?
  t11 <- hist_perm(exp8, obs8, perms2)+labs(title='does LPS affect eigenvector centrality in treatment period?')
  t11.1 <- hist_perm(exp8.1, obs8.1, perms2)+labs(title='does LPS affect eigenvector centrality in treatment period?')
  # does LPS have effect outside the treatment period?
  t12 <- hist_perm(exp9, obs9, perms2)+labs(title='does LPS affect eigenvector centrality in post-treatment period?')
  
  # combine plots
  (perm.test2 <- t1+t2+t3+t3.1+t4+t5+t6+t7+t7.1+t8+t9+t10+t11+t11.1+t12+plot_layout(ncol=3))
  ggsave("perm_tests2.pdf", width=18, height= 15, units="in", dpi=1200)
  
  # make results table
  results2 <- 
    data.frame(response= rep(c('degree', 'strength', 'eigenvector'), times=1, each=5),
               fixed_effect= rep(c('interaction', 'treatment', 'treatment.within.period','treatment.day2', 'treatment.within.post'), times=3),
               coefficient= c(obs1, obs1.1, obs2, obs2.1, obs3, obs4, obs4.1, obs5, obs5.1, obs6, obs7, obs7.1, obs8, obs8.1, obs9),
               pvalue= c(mean(exp1>=obs1),
                         mean(exp1.1>=obs1.1),
                         mean(exp2>=obs2),
                         mean(exp2.1>=obs2.1),
                         mean(exp3>=obs3),
                         mean(exp4>=obs4),
                         mean(exp4.1>=obs4.1),
                         mean(exp5>=obs5),
                         mean(exp5.1>=obs5.1),
                         mean(exp6>=obs6),
                         mean(exp7>=obs7),
                         mean(exp7.1>=obs7.1),
                         mean(exp8>=obs8),
                         mean(exp8.1>=obs8.1),
                         mean(exp9>=obs9))) %>% 
    mutate(pvalue= if_else(coefficient<0, (1-pvalue), pvalue)) %>% 
    mutate(pvalue2= c(mean(exp1>=obs1),
                      mean(exp1.1>=obs1.1),
                      mean(exp2>=obs2),
                      mean(exp2.1>=obs2.1),
                      mean(exp3>=obs3),
                      mean(exp4>=obs4),
                      mean(exp4.1>=obs4.1),
                      mean(exp5>=obs5),
                      mean(exp5.1>=obs5.1),
                      mean(exp6>=obs6),
                      mean(exp7>=obs7),
                      mean(exp7.1>=obs7.1),
                      mean(exp8>=obs8),
                      mean(exp8.1>=obs8.1),
                      mean(exp9>=obs9)),
           pvalue2= c(mean(abs(exp1) >= abs(obs1)),
                      mean(abs(exp1.1) >= abs(obs1.1)),
                      mean(abs(exp2) >= abs(obs2)),
                      mean(abs(exp2.1) >= abs(obs2.1)),
                      mean(abs(exp3) >= abs(obs3)),
                      mean(abs(exp4) >= abs(obs4)),
                      mean(abs(exp4.1) >= abs(obs4.1)),
                      mean(abs(exp5) >= abs(obs5)),
                      mean(abs(exp5.1) >= abs(obs5.1)),
                      mean(abs(exp6) >= abs(obs6)),
                      mean(abs(exp7) >= abs(obs7)),
                      mean(abs(exp7.1) >= abs(obs7.1)),
                      mean(abs(exp8) >= abs(obs8)),
                      mean(abs(exp8.1) >= abs(obs8.1)),
                      mean(abs(exp9) >= abs(obs9))))
  results2
  time2= Sys.time() - time1
  time2
}



# get associations by dyad type----

# label edge types in hourly data
dt <- 
  d %>%
  mutate(dyad= paste(bat1, bat2, sep='_')) %>% 
  group_by(dyad, dyad_type) %>% 
  summarize(n=n())
h$dyad_type <- dt$dyad_type[match(h$dyad, dt$dyad)]

# get mean association time and prob per edge type per hour
e <- 
  h %>% 
  mutate(hour= as.POSIXct(hour, format= "%Y-%m-%d %H", tz="CST6CDT")) %>% 
  group_by(hour, dyad_type) %>% 
  summarize(`mean association duration`= mean(duration), 
            `mean association probability`= mean(duration>0)) %>% 
  ungroup() %>% 
  filter(hour <= as.POSIXct("2018-04-28 01:00:00", tz = "CST6CDT")) %>% 
  pivot_longer(cols= 'mean association duration' : 'mean association probability', 
               values_to = 'value', names_to = 'measure') %>% 
  mutate(period= case_when(
    hour >= treatment_start & hour < treatment_stop ~ 'treatment',
    hour >= post24_start & hour < post24_stop ~ '24 hours later',
    hour >= post48_start & hour < post48_stop ~ '48 hours later',
    TRUE ~"none"))

# plot dyad-type associations per hour-----
(eh.plot <- 
   e %>% 
   ggplot(aes(x=hour, y=value, color= dyad_type, shape=dyad_type, linetype= dyad_type))+
   facet_wrap(~measure, ncol=1, scales="free_y")+  
   geom_point()+
   geom_line()+
   geom_rect(aes(xmin=treatment_stop, xmax=post24_start, ymin=0, ymax=Inf), fill='white', color=NA, alpha=0.01)+
   geom_rect(aes(xmin=post24_stop, xmax=post48_start, ymin=0, ymax=Inf), fill='white', color=NA, alpha=0.01)+
   geom_rect(aes(xmin=post48_stop, xmax=max(d2$hour)+1500, ymin=0, ymax=Inf), fill='white', color=NA, alpha=0.01)+
   geom_vline(xintercept = treatment_start, color= "black")+
   geom_vline(xintercept=treatment_stop, color= "black")+
   geom_vline(xintercept = post24_start, linetype= "solid")+
   geom_vline(xintercept=post24_stop, linetype= "solid")+
   geom_vline(xintercept = post48_start, linetype= "solid")+
   geom_vline(xintercept=post48_stop, linetype= "solid")+
   theme_bw()+
   ylab("")+
   theme(legend.position = c(0.18, 0.58),
         legend.margin=margin(c(0.5,0.5,0.5,0.5)),
         legend.justification = "left",
         legend.title = element_blank(),
         legend.direction = "horizontal",
         legend.background = element_rect(fill='white', 
                                          size=0.5, linetype="solid"))+
   scale_x_datetime(breaks=date_breaks("6 hour"), labels=date_format("%H:%M", tz = "CST6CDT")) + 
  scale_color_manual(values=c("black", "red", "springgreen"))+
  scale_linetype_manual(values= c('solid', 'dashed', 'dotted'))+
  scale_shape_manual(values=c("circle", "triangle", "square")))

# get association times and probs by dyad type during treatment and post-treatment periods----
d4 <- 
  h %>% 
  mutate(hour= as.POSIXct(hour, format= "%Y-%m-%d %H", tz="CST6CDT")) %>% 
  filter(hour <= as.POSIXct("2018-04-28 01:00:00", tz = "CST6CDT")) %>% 
  mutate(period= case_when(
    hour >= treatment_start & hour < treatment_stop ~ 'treatment',
    hour >= post24_start & hour < post24_stop ~ '24 hours later',
    hour >= post48_start & hour < post48_stop ~ '48 hours later',
    TRUE ~"none")) %>% 
  filter(period!='none') %>% 
  group_by(hour, dyad, dyad_type, period) %>% 
  summarize(time= mean(duration), 
            prob= mean(duration>0)) %>% 
  group_by(dyad, dyad_type, period) %>% 
  summarize(time= mean(time), 
            prob= mean(prob),
            n= n()) %>% 
  ungroup()

# get means and 95% CI
assoc.untreated <- 
  d4 %>% 
  filter(dyad_type=='control-control') %>% 
  select(period, time, prob) %>% 
  boot_ci2(y=.$time, x=.$period) %>% 
  mutate(type= 'control-control')

assoc.mixed <- 
  d4 %>% 
  filter(dyad_type=='sick-control') %>% 
  select(period, time, prob) %>% 
  boot_ci2(y=.$time, x=.$period) %>% 
  mutate(type= 'sick-control')

prob.untreated <- 
  d4 %>% 
  filter(dyad_type=='control-control') %>% 
  select(period, time, prob) %>% 
  boot_ci2(y=.$prob, x=.$period) %>% 
  mutate(type= 'control-control')

prob.mixed <- 
  d4 %>% 
  filter(period!='none') %>% 
  filter(dyad_type=='sick-control') %>% 
  select(period, time, prob) %>% 
  boot_ci2(y=.$prob, x=.$period) %>% 
  mutate(type= 'sick-control')

# plot association time
p1 <- 
    rbind(assoc.untreated, assoc.mixed) %>% 
    mutate(period= factor(effect, levels= c('treatment', '24 hours later', '48 hours later'))) %>% 
    ggplot(aes(x=type, y=mean, color=type, shape=type))+
    facet_wrap(~period)+
    geom_point(size=3)+
    geom_errorbar(aes(ymin=low, ymax=high, width=.1), size=1)+
    scale_color_manual(values=c("black", 'red'))+
    scale_shape_manual(values=c('circle', 'triangle'))+
    theme_bw()+
  theme(legend.position= 'none',
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())
    

# plot association prob
p2 <- 
  rbind(prob.untreated, prob.mixed) %>% 
  mutate(period= factor(effect, levels= c('treatment', '24 hours later', '48 hours later'))) %>% 
  mutate(type= ifelse(type=="control-control", "c-c", 's-c')) %>% 
  ggplot(aes(x=type, y=mean, color=type, shape=type))+
  facet_wrap(~period)+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=low, ymax=high, width=.1), size=1)+
  scale_color_manual(values=c("black", 'red'))+
  scale_shape_manual(values=c('circle', 'triangle'))+
  theme_bw()+
  theme(legend.position= 'none',
        axis.title.y=element_blank())+
  xlab("dyad type")


# plot and make PDF
(eh.plot2 <- (eh.plot|(p1+p2+plot_layout(ncol=1)))+plot_layout(widths= c(2,1)))
ggsave("dyad_types_by_hour.pdf", width=11, height= 5, units="in", dpi=1200)
  

# save output with timestamp-----
timestamp <- substr(gsub(x=gsub(":","_",Sys.time()), 
                         pattern=" ", replace="_"), start=1, stop=16)
workspace <- paste("results_", timestamp, ".Rdata", sep="")
save.image(file= workspace)

# print results table
results
results2

# print current results file (paste and copy below)
workspace 

# code to load old data
if(FALSE){
  load("results_2020-09-16_17_12.Rdata")    
}



