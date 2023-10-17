# analysis code in R

# 1) for calculating allele frequencies from filtered VCF file 
# 2) for calculating pairwise FST
# 3) for analysing patterns of genomic diversity
# 4) for PCA and t-SNE analysis

# the VCF file used in these analyses was generated using a snakemake workflow available from https://github.com/pascalangst/Angst_etal_2022_MBE and filtered using filter_VCF.sh

# load dependencies
library(poolfstat)
library(data.table)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(nls2)
library(ggprism)
library(lmerTest)
library(plyr)
library(ggsignif)
library(pcadapt)
library(tibble)
library(Rtsne)

## 1) calculating allele frequencies from filtered VCF file 

# get poolsizes (while keeping order of samples in vcf)
# this file was generated as descibed in https://doi.org/10.1101/2020.07.10.196758 and https://doi.org/10.1093/molbev/msac264
covariates <- read.csv("data/metadata_subpopulations/all.covariates.csv", sep=",")
covariates$season <- ifelse(covariates$sample == "1","spr", ifelse(covariates$sample == "2","smr",NA))
covariates$samplename <- paste0(covariates$poolname, "_", covariates$season, covariates$year)
samplenames <- system("bcftools query -l Dmagna_2014-2018_filtered.vcf.gz", intern = TRUE) # poolnames extracted from vcf
samplenames.frame <- data.frame(samplename=samplenames, row=1:length(samplenames))
covariates.merge <- merge(covariates, samplenames.frame, by="samplename")
covariates.merge <- covariates.merge[order(covariates.merge$row), ]
covariates.merge$animal_number[is.na(covariates.merge$animal_number)] <- 50

# covert vcf to pooldata:
pooldata <- vcf2pooldata("Dmagna_2014-2018_filtered.vcf.gz", poolsizes=2*covariates.merge$animal_number, poolnames=samplenames)

# get allele frequencies
data <- pooldata@refallele.readcount/pooldata@readcoverage
data[is.na(data)] <- 0

# get folded sfs
data_t <- transpose(as.data.frame(data))
rownames(data_t) <- pooldata@poolnames
m <- as.matrix(data_t)
m[m > 0.5] <- 1-m[m > 0.5] # this is much faster as matrix
data_t <- as.data.frame(m)

# get median MAFs
medians <- vector(mode = "list")
for (i in 1:nrow(data_t)) {
  medians[[i]] <- median(data_t[i,][data_t[i,] > 0])
}

medians.frame <- as.data.frame(cbind(samplename=rownames(data_t), medians))

# get info about immigration events
timeseries_info <- read.csv("data/2014-2018_MetapopData_timeseries_all.csv", stringsAsFactors = FALSE)

# get island (island_group)
medians.frame$island <- gsub(pattern = "-[0-9A]*_s.r201.", "", medians.frame$samplename)
medians.frame$island_group <- gsub(pattern = "^FS*", "FS*", gsub(pattern = "^SK.*", "SK*", gsub(pattern = "^LON.*", "LON*", medians.frame$island)))

# load grenedalf script (grenedalf command follows below)
pi.frame <- read.csv("grenedalf.csv", header=F, row.names = 1, sep = "\t")
pi.frame$samplename <- rownames(pi.frame)
colnames(pi.frame) <- c("snp_count", "cov_fract", "pi_abs", "pi", "tetha_abs", "tetha_rel", "tajimaD", "NA", "samplename")
covariates <- merge(covariates, pi.frame, by= "samplename")

# check if extinction-recolonization happened
cov_ponds <- vector(mode = "list")
for (i in ponds) cov_ponds[[i]] <- covariates[grep(paste0("^",gsub(".", "-", i, fixed=TRUE),"_"), covariates$samplename),][c("samplename", "medians", "tetha_rel", "age", "mean_dist_to_closest_2")]

# get relative age
for (i in 1:length(cov_ponds)) {
  if (is.data.frame(cov_ponds[[i]])){
    cov_ponds[[i]]["year"] <- gsub("^.*_s.r", "", cov_ponds[[i]]$samplename)
    cov_ponds[[i]]["season"] <- gsub("201.", "", gsub("^.*_", "", cov_ponds[[i]]$samplename))
    cov_ponds[[i]]["year"] <- as.numeric(as.data.frame(cov_ponds[[i]])$year) - 2014
    cov_ponds[[i]]["season"] <- ifelse(cov_ponds[[i]]["season"]=="spr", 0, 0.2)
  }
}

# vectors for pools with recol events
recol.vector <- vector()
last.sample <- vector()
last.sample.age <- vector()
first.sample <- vector()
first.sample.age <- vector()

cov_ponds_2 <- vector(mode = "list")

# get recol events
for (i in 1:length(cov_ponds)) {
  if (is.data.frame(cov_ponds[[i]])){
    
    # order samples relative age
    df.temp <- cov_ponds[[i]][order(cov_ponds[[i]]["year"], cov_ponds[[i]]["season"]),]
    
    # get discontinuous ages
    if(any(diff(df.temp$age) <= 0) & !all(is.na(df.temp$age))){
      recol.vector <- c(recol.vector, names(cov_ponds[i]))
      if(nrow(df.temp) > 2){
        last.sample <- c(last.sample, df.temp[diff(df.temp$age) <= 0,]$samplename)
        last.sample.age <- c(last.sample.age, df.temp[diff(df.temp$age) <= 0,]$age)
        first.sample <- c(first.sample, df.temp[c(F, diff(df.temp$age) <= 0),]$samplename)
        first.sample.age <- c(first.sample.age, df.temp[c(F, diff(df.temp$age) <= 0),]$age)
      } 
      else{
        last.sample <- c(last.sample, df.temp[1,]$samplename)
        last.sample.age <- c(last.sample.age, df.temp[1,]$age)
        first.sample <- c(first.sample, df.temp[2,]$samplename)
        first.sample.age <- c(first.sample.age, df.temp[2,]$age)
      }
      cov_ponds_2[[names(cov_ponds[i])]] <- df.temp[1:(which(diff(df.temp$age) <= 0)),]
      cov_ponds_2[[paste0(names(cov_ponds[i]), "B")]] <- df.temp[(1+which(diff(df.temp$age) <= 0)):nrow(df.temp),]
      
      cov_ponds_2[[names(cov_ponds[i])]]["samplingtime"] <- cov_ponds_2[[names(cov_ponds[i])]]["year"] +cov_ponds_2[[names(cov_ponds[i])]]["season"]
      cov_ponds_2[[names(cov_ponds[i])]]["relage"] <- cov_ponds_2[[names(cov_ponds[i])]]["samplingtime"] - min(cov_ponds_2[[names(cov_ponds[i])]]["samplingtime"])
      
      cov_ponds_2[[paste0(names(cov_ponds[i]), "B")]]["samplingtime"] <- cov_ponds_2[[paste0(names(cov_ponds[i]), "B")]]["year"] +cov_ponds_2[[paste0(names(cov_ponds[i]), "B")]]["season"]
      cov_ponds_2[[paste0(names(cov_ponds[i]), "B")]]["relage"] <- cov_ponds_2[[paste0(names(cov_ponds[i]), "B")]]["samplingtime"] - min(cov_ponds_2[[paste0(names(cov_ponds[i]), "B")]]["samplingtime"])
      
    }
    else{cov_ponds_2[[names(cov_ponds[i])]] <- df.temp
      cov_ponds_2[[names(cov_ponds[i])]]["samplingtime"] <- cov_ponds_2[[names(cov_ponds[i])]]["year"] +cov_ponds_2[[names(cov_ponds[i])]]["season"]
      cov_ponds_2[[names(cov_ponds[i])]]["relage"] <- cov_ponds_2[[names(cov_ponds[i])]]["samplingtime"] - min(cov_ponds_2[[names(cov_ponds[i])]]["samplingtime"])
    }
  }
}

# remove also samples for false positives (only consisting of 1 or 2 samples)
false_positives <- c("FS.27B", "FS.31B", "N.71B")

# this will be used below to visualise recolonization events
true_pos_indices <- grep(invert = T, pattern = paste(gsub("\\.","-", gsub("B","",false_positives)), collapse = "|"), x = first.sample)
last.sample <- last.sample[true_pos_indices]
last.sample.age <- last.sample.age[true_pos_indices]
first.sample <- first.sample[true_pos_indices] 
first.sample.age <- first.sample.age[true_pos_indices]
cov_ponds_2 <- cov_ponds_2[!names(cov_ponds_2) %in% false_positives]

# investigate relationships to covariates ...

# ... for example to age
cov_ponds_2_age <- as.data.frame(do.call("rbind",cov_ponds_2))

# ... relative age
summary(lmerTest::lmer(medians ~ relage + (1|island_group) + (1|island_group:pop), cov_ponds_2_age %>% group_by(pop) %>% mutate(long_series = relage >= 2) %>% filter(any(long_series))))
anova(lmerTest::lmer(medians ~ relage + (1|island_group) + (1|island_group:pop), cov_ponds_2_age %>% group_by(pop) %>% mutate(long_series = relage >= 2) %>% filter(any(long_series))))

# ... age class
cov_ponds_2_age_class <- as.data.frame(cov_ponds_2_age %>%
  group_by(pop) %>%
  summarize(medians_mean = mean(medians, na.rm=TRUE),
            tetha_mean = mean(tetha_rel, na.rm=TRUE),
            age_mean = mean(age, na.rm=TRUE)) %>%
  mutate(Age_class = case_when(
    age_mean <= 2 ~ "Young",
    age_mean > 2 & age_mean <= 15 ~ "Intermediate",
    age_mean > 15 ~ "Old")))

t.test(cov_ponds_2_age_class[cov_ponds_2_age_class$Age_class=="Young",]$medians_mean, cov_ponds_2_age_class[cov_ponds_2_age_class$Age_class=="Intermediate",]$medians_mean)
t.test(cov_ponds_2_age_class[cov_ponds_2_age_class$Age_class=="Young",]$medians_mean, cov_ponds_2_age_class[cov_ponds_2_age_class$Age_class=="Old",]$medians_mean)
t.test(cov_ponds_2_age_class[cov_ponds_2_age_class$Age_class=="Intermediate",]$medians_mean, cov_ponds_2_age_class[cov_ponds_2_age_class$Age_class=="Old",]$medians_mean)

young.vector.relage <- unique(as.data.frame(cov_ponds_2_age %>% group_by(pop) %>% filter(row_number() %in% c(1, n())) %>% mutate(long_series = relage >= 2 & lag(age <= 4, n=1)) %>% filter(any(long_series))) %>% pull(pop))
old.frame <- as.data.frame(cov_ponds_2_age %>% group_by(pop) %>% mutate(long_series = relage >= 2) %>% filter(age > 4, any(long_series)))
age_adj_names <- cov_ponds_2_age %>% group_by(pop) %>% mutate(long_series = relage >= 2) %>% filter(age <= 2, any(long_series)) %>% pull(pop)

old.frame[old.frame$pop %in% age_adj_names,]$relage <- old.frame[old.frame$pop %in% age_adj_names,]$relage - 2
old.frame <- old.frame %>% group_by(pop) %>% mutate(long_series = relage >= 2) %>% filter(any(long_series))

summary(lmerTest::lmer(medians ~ relage + (1|island_group) + (1|island_group:pop), cov_ponds_2_age[cov_ponds_2_age$pop %in% young.vector.relage,]))
anova(lmerTest::lmer(medians ~ relage + (1|island_group) + (1|island_group:pop), cov_ponds_2_age[cov_ponds_2_age$pop %in% young.vector.relage,]))

summary(lmerTest::lmer(medians ~ relage + (1|island_group) + (1|island_group:pop), old.frame))
anova(lmerTest::lmer(medians ~ relage + (1|island_group) + (1|island_group:pop), old.frame))


# ... hybrid vigor yes/no
# vector for samples from before and after hybrid vigor
FSTdata_imm_adjacent <- vector()
FSTdata_imm_adjacent_sample <- vector()
gain_vs_loss_imm_adjacent <- vector()
names_imm_adjacent <- vector()
# vector for ponds without hybrid vigor
FSTdata_noimm_adjacent <- vector()
FSTdata_noimm_adjacent_sample <- vector()
gain_vs_loss_noimm_adjacent <- vector()
names_noimm_adjacent <- vector()
# vector for stable ponds
FSTdata_stable_adjacent <- vector()
gain_vs_loss_stable_adjacent <- vector()
names_stable_adjacent <- vector()
FSTdata_rest_adjacent <- vector()
gain_vs_loss_rest_adjacent <- vector()
names_rest_adjacent <- vector()

samples.not_after_recol <- cov_ponds_2_age[grep("B",cov_ponds_2_age$pop, invert = T),]$samplename
set.seed(seed = 100)

stable.char <- c("stable" ,
"false_pos_recolonization" , # until (not-)event
"stable; potentially founded by two",
"stable; additional peak",
"recolonization" # until event; new subpops all have just one sample, i.e., not considered here
)

unstable.char <- c(
"potential recolonization - then recolonization",
"potential recolonization - then hybrid vigor - then recolonization",
"potential hybrid vigor",
"potential recolonization or hybrid vigor with bottleneck"  ,
"potential recolonization" ,
"hybrid vigor",
"potential recolonization or hybrid vigor with bottleneck; strange because all late samples show high FST",
"potential recolonization or hybrid vigor with bottleneck; strange because it gets more similar again; plus additional peak",
"recolonization based on multiple founders or false positive recolonization but hybrid vigor"             ,
"potential recolonization or hybrid vigor",
"potential recolonization or hybrid vigor with bottleneck; strange because it gets more similar again"
)

potImmig.char <- c(
  "potential hybrid vigor",
  "hybrid vigor",
  "potential recolonization or hybrid vigor with bottleneck"  ,
  "potential recolonization or hybrid vigor with bottleneck; strange because it gets more similar again",
  "potential recolonization or hybrid vigor with bottleneck; strange because all late samples show high FST",
  "potential recolonization or hybrid vigor with bottleneck; strange because it gets more similar again; plus additional peak",
  "potential recolonization or hybrid vigor",
  "potential recolonization - then hybrid vigor - then recolonization"
)

for (i in 1:length(FSTdata_ponds)) {
  if(any(names(FSTdata_ponds) != names(gain_vs_loss_ponds))) {next}
  
  if (is.data.frame(FSTdata_ponds[[i]])){
    FSTdata_pond_temp <- FSTdata_ponds[[i]]
    FSTdata_pond_temp$year <- gsub("^.*_s.r", "", rownames(FSTdata_pond_temp)) 
    FSTdata_pond_temp$season <- gsub("201.", "", gsub("^.*_", "", rownames(FSTdata_pond_temp)))
    FSTdata_pond_temp$year <- as.numeric(FSTdata_pond_temp$year) - 2014
    FSTdata_pond_temp$season <- ifelse(FSTdata_pond_temp$season=="spr", 0, 0.2)
    FSTdata_pond_temp <- FSTdata_pond_temp[order(FSTdata_pond_temp$year, FSTdata_pond_temp$season),order(FSTdata_pond_temp$year, FSTdata_pond_temp$season)]
    
    gain_vs_loss_pond_temp <- gain_vs_loss_ponds[[i]]
    gain_vs_loss_pond_temp$year <- gsub("^.*_s.r", "", rownames(gain_vs_loss_pond_temp)) 
    gain_vs_loss_pond_temp$season <- gsub("201.", "", gsub("^.*_", "", rownames(gain_vs_loss_pond_temp)))
    gain_vs_loss_pond_temp$year <- as.numeric(gain_vs_loss_pond_temp$year) - 2014
    gain_vs_loss_pond_temp$season <- ifelse(gain_vs_loss_pond_temp$season=="spr", 0, 0.2)
    gain_vs_loss_pond_temp <- gain_vs_loss_pond_temp[order(gain_vs_loss_pond_temp$year, gain_vs_loss_pond_temp$season),order(gain_vs_loss_pond_temp$year, gain_vs_loss_pond_temp$season)]
    
    FSTdata_pond_temp <- FSTdata_pond_temp[colnames(FSTdata_pond_temp) %in% samples.not_after_recol, colnames(FSTdata_pond_temp) %in% samples.not_after_recol]
    gain_vs_loss_pond_temp <- gain_vs_loss_pond_temp[colnames(gain_vs_loss_pond_temp) %in% samples.not_after_recol, colnames(gain_vs_loss_pond_temp) %in% samples.not_after_recol]
    
    if (is.data.frame(FSTdata_pond_temp)){
      if (any(colnames(FSTdata_pond_temp) %in% unlist(timeseries_info[timeseries_info$comment %in% potImmig.char,c("X.1", "X.2")], use.names = F))) {
        print(names(FSTdata_ponds[i]))
        names.vector <- unlist(timeseries_info[timeseries_info$comment %in% potImmig.char,c("X.1", "X.2")], use.names = F)[unlist(timeseries_info[timeseries_info$comment %in% potImmig.char,c("X.1", "X.2")], use.names = F) %in% names(FSTdata_pond_temp)]
        
        FSTdata_imm_adjacent <- c(FSTdata_imm_adjacent, FSTdata_pond_temp[names.vector[1], names.vector[2]])
        #gain_vs_loss_imm_adjacent <- c(gain_vs_loss_imm_adjacent, gain_vs_loss_pond_temp[names.vector[1], names.vector[2]])
        gain_vs_loss_imm_adjacent <- c(gain_vs_loss_imm_adjacent, mean(gain_vs_loss_pond_temp[col(gain_vs_loss_pond_temp)+1==row(gain_vs_loss_pond_temp)]))
        names_imm_adjacent <- c(names_imm_adjacent, names(FSTdata_ponds[i]))
      }else{
        FSTdata_noimm_adjacent <- c(FSTdata_noimm_adjacent, mean(FSTdata_pond_temp[col(FSTdata_pond_temp)+1==row(FSTdata_pond_temp)]))
        FSTdata_noimm_adjacent_sample <- c(FSTdata_noimm_adjacent_sample, sample(FSTdata_pond_temp[col(FSTdata_pond_temp)+1==row(FSTdata_pond_temp)], 1))
        gain_vs_loss_noimm_adjacent <- c(gain_vs_loss_noimm_adjacent, mean(gain_vs_loss_pond_temp[col(gain_vs_loss_pond_temp)+1==row(gain_vs_loss_pond_temp)]))
        names_noimm_adjacent <- c(names_noimm_adjacent, names(FSTdata_ponds[i]))
      }
      if(names(FSTdata_ponds[i]) %in% gsub("\\.", "-", timeseries_info[timeseries_info$comment %in% stable.char,"X"])){
        #print(names(FSTdata_ponds[i]))
        
        FSTdata_stable_adjacent <- c(FSTdata_stable_adjacent, mean(FSTdata_pond_temp[col(FSTdata_pond_temp)+1==row(FSTdata_pond_temp)]))
        gain_vs_loss_stable_adjacent <- c(gain_vs_loss_stable_adjacent, mean(gain_vs_loss_pond_temp[col(gain_vs_loss_pond_temp)+1==row(gain_vs_loss_pond_temp)]))
        names_stable_adjacent <- c(names_stable_adjacent, names(FSTdata_ponds[i]))
      }else{
        if(!any(colnames(FSTdata_pond_temp) %in% unlist(timeseries_info[timeseries_info$comment %in% potImmig.char,c("X.1", "X.2")], use.names = F))){
          FSTdata_rest_adjacent <- c(FSTdata_rest_adjacent, max(FSTdata_pond_temp[col(FSTdata_pond_temp)+1==row(FSTdata_pond_temp)]))
          gain_vs_loss_rest_adjacent <- c(gain_vs_loss_rest_adjacent, mean(gain_vs_loss_pond_temp[col(gain_vs_loss_pond_temp)+1==row(gain_vs_loss_pond_temp)]))
          names_rest_adjacent <- c(names_rest_adjacent, names(FSTdata_ponds[i]))
        }
      }
    }
  }
}  

cov_ponds_2_age <- merge(cov_ponds_2_age, timeseries_info, by.x = "pop", by.y = "X", all = T)

# stable
summary(lmerTest::lmer(medians ~ relage + (1|island_group) + (1|island_group:pop), cov_ponds_2_age[cov_ponds_2_age$comment %in% stable.char,]))
anova(lmerTest::lmer(medians ~ relage + (1|island_group) + (1|island_group:pop), cov_ponds_2_age[cov_ponds_2_age$comment %in% stable.char,]))

long_series.frame <- as.data.frame(cov_ponds_2_age %>% group_by(pop) %>% mutate(long_series = relage >= 2) %>% filter(any(long_series)))
summary(lmerTest::lmer(medians ~ relage + (1|island_group) + (1|island_group:pop), long_series.frame))
anova(lmerTest::lmer(medians ~ relage + (1|island_group) + (1|island_group:pop), long_series.frame))

summary(lmerTest::lmer(medians ~ relage + (1|island_group) + (1|island_group:pop), long_series.frame[long_series.frame$pop %in% young.vector.relage,]))
anova(lmerTest::lmer(medians ~ relage + (1|island_group) + (1|island_group:pop), long_series.frame[long_series.frame$pop %in% young.vector.relage,]))

long_series.frame.no_immigration <- as.data.frame(long_series.frame %>%
                                            mutate(Stable_class = case_when(
                                              comment %in% stable.char ~ "Stable",
                                              comment %in% unstable.char ~ "Unstable"),
                                              Immigration = case_when(
                                                comment %in% potImmig.char ~ "Yes",
                                                TRUE ~ "No"),
                                              `comb` = paste0(Stable_class, "_", Immigration)))

summary(lmerTest::lmer(medians ~ relage + (1|island_group) + (1|island_group:pop), long_series.frame.no_immigration[long_series.frame.no_immigration$Immigration == "No",]))
anova(lmerTest::lmer(medians ~ relage + (1|island_group) + (1|island_group:pop), long_series.frame.no_immigration[long_series.frame.no_immigration$Immigration == "No",]))

long_series.frame.no_immigration <- merge(long_series.frame.no_immigration, long_series.frame.no_immigration[long_series.frame.no_immigration$X.1 == long_series.frame.no_immigration$samplename & long_series.frame.no_immigration$comb == "Unstable_Yes", c("pop", "relage")], by = "pop", all = T, suffixes = c("", "_X.1"))
long_series.frame.no_immigration$relage_adj <- long_series.frame.no_immigration$relage - long_series.frame.no_immigration$relage_X.1

old.frame <- as.data.frame(cov_ponds_2_age %>% group_by(pop) %>% mutate(long_series = relage >= 2) %>% filter(age > 2, any(long_series)))
age_adj_names <- cov_ponds_2_age %>% group_by(pop) %>% mutate(long_series = relage >= 2) %>% filter(age <= 2, any(long_series)) %>% pull (pop)

old.frame[old.frame$pop %in% age_adj_names,]$relage <- old.frame[old.frame$pop %in% age_adj_names,]$relage - 2
old.frame <- old.frame %>% group_by(pop) %>% mutate(long_series = relage >= 2) %>% filter(any(long_series))

old.frame.no_immigration <- as.data.frame(old.frame %>%
  mutate(Stable_class = case_when(
    comment %in% stable.char ~ "Stable",
    comment %in% unstable.char ~ "Unstable"),
    Immigration = case_when(
      comment %in% potImmig.char ~ "Yes",
      TRUE ~ "No"),
    `comb` = paste0(Stable_class, "_", Immigration)))

summary(lmerTest::lmer(medians ~ relage + (1|island_group) + (1|island_group:pop), old.frame.no_immigration[old.frame.no_immigration$comb == "Unstable_Yes",]))
anova(lmerTest::lmer(medians ~ relage + (1|island_group) + (1|island_group:pop), old.frame.no_immigration[old.frame.no_immigration$comb == "Unstable_Yes",]))


## Metapop. overall SFS per timepoint

# Define a function to extract the timepoint ID
get_timepoints <- function(x) {
  sub("^[A-Za-z]+[^_]*_(s.r201.)$", "\\1", x)
}

# Split the dataframe by timepoint
data_t_timepoints <- split(data_t, get_timepoints(rownames(data_t)))

# average per SNP frequency
colMeans_dTimepoints <- data.frame()

for (i in 1:length(data_t_timepoints)) {
  colMeans_dTimepoints <- rbind(colMeans_dTimepoints,cbind(average_freq=colMeans(data_t_timepoints[[i]]), timepoint=names(data_t_timepoints[i])))
}

colMeans_dTimepoints$average_freq <- as.numeric(colMeans_dTimepoints$average_freq)
colMeans_dTimepoints$timepoint <- as.factor(colMeans_dTimepoints$timepoint)

# median overall MAF
med <- colMeans_dTimepoints[colMeans_dTimepoints$average_freq > 0,] %>%
  group_by(timepoint) %>%
  summarize(median = median(average_freq))


## 2) calculating pairwise FST

PairwiseFST<-compute.pairwiseFST(pooldata)
FSTdata_all<- PairwiseFST@PairwiseFSTmatrix
rownames(FSTdata_all)<-gsub("-", ".", rownames(FSTdata_all))
FSTdata_all[FSTdata_all < 0 & !is.na(FSTdata_all)] <- 0

# investigate relationships to covariates ...
covariates <- merge(covariates, as.data.frame(cbind(samplename=gsub("\\.", "-", colnames(FSTdata_all)))), by="samplename")

# ... per pond
# get pond names
ponds <- unique(gsub("_s.r201.", "", colnames(FSTdata_all)))

# aggregate samples from each pond
FSTdata_ponds <- vector(mode = "list")
for (i in ponds) FSTdata_ponds[[i]] <- FSTdata_all[grep(paste0("^",i,"_"), colnames(FSTdata_all)),grep(paste0("^",i,"_"), colnames(FSTdata_all))]

# relative age from 2014
for (i in 1:length(FSTdata_ponds)) {
  if (is.data.frame(FSTdata_ponds[[i]])){
	  FSTdata_ponds[[i]]["year"] <- gsub("^.*_s.r", "", rownames(FSTdata_ponds[[i]])) 
	  FSTdata_ponds[[i]]["season"] <- gsub("201.", "", gsub("^.*_", "", rownames(FSTdata_ponds[[i]])))
	  FSTdata_ponds[[i]]["year"] <- as.numeric(as.data.frame(FSTdata_ponds[[i]])$year) - 2014
	  FSTdata_ponds[[i]]["season"] <- ifelse(FSTdata_ponds[[i]]["season"]=="spr", 0, 0.2)
  }
}

dict <- list("0" = "spr2014", "0.2" = "smr2014",
             "1" = "spr2015", "1.2" = "smr2015",
             "2" = "spr2016", "2.2" = "smr2016",
             "3" = "spr2017", "3.2" = "smr2017",
             "4" = "spr2018", "4.2" = "smr2018")

# function to get absolute differences between samples
custom_fun <- function(x, y) {
    z <- abs(x - y)
    return(z)
}

# check if extinction-recolonization happened
# get name and age
cov_ponds <- vector(mode = "list")
for (i in ponds) cov_ponds[[i]] <- covariates[grep(paste0("^",gsub(".", "-", i, fixed=TRUE),"_s"), covariates$samplename),][c("samplename","age", "mean_dist_to_closest_2")]

# get relative age
for (i in 1:length(cov_ponds)) {
  if (is.data.frame(cov_ponds[[i]])){
    cov_ponds[[i]]["year"] <- gsub("^.*_s.r", "", cov_ponds[[i]]$samplename) 
    cov_ponds[[i]]["season"] <- gsub("201.", "", gsub("^.*_", "", cov_ponds[[i]]$samplename))
    cov_ponds[[i]]["year"] <- as.numeric(as.data.frame(cov_ponds[[i]])$year) - 2014
    cov_ponds[[i]]["season"] <- ifelse(cov_ponds[[i]]["season"]=="spr", 0, 0.2)
  }
}

# vector for pools with recol events and list for such without
recol.vector <- vector()
FSTdata_cov_ponds <- vector(mode = "list")

# get recol events
for (i in names(cov_ponds)) {
  if (is.data.frame(cov_ponds[[i]])){
    
    # order samples relative age
    df.temp <- cov_ponds[[i]][order(cov_ponds[[i]]["year"], cov_ponds[[i]]["season"]),]
    
    if(nrow(df.temp)>1){
      df.temp2 <- FSTdata_ponds[[i]][order(FSTdata_ponds[[i]]["year"], FSTdata_ponds[[i]]["season"]),order(FSTdata_ponds[[i]]["year"], FSTdata_ponds[[i]]["season"])]
      
      
      df.temp2$year <- gsub("^.*_s.r", "", rownames(df.temp2)) 
      df.temp2$season <- gsub("201.", "", gsub("^.*_", "", rownames(df.temp2)))
      df.temp2$year <- as.numeric(df.temp2$year) - 2014
      df.temp2$season <- ifelse(df.temp2$season=="spr", 0, 0.2)
      df.temp2$age <- df.temp$age
      df.temp2$mean_dist_to_closest_2 <- df.temp$mean_dist_to_closest_2
    
      # get discontinuous ages
      if(any(diff(df.temp$age) <= 0) & !all(is.na(df.temp$age))){
        recol.vector <- c(recol.vector, i)
    
        FSTdata_cov_ponds[[i]] <- df.temp2[1:(which(diff(df.temp2$age) <= 0)),c(1:(which(diff(df.temp2$age) <= 0)), (ncol(df.temp2)-3):ncol(df.temp2))]
        FSTdata_cov_ponds[[paste0(i, "B")]] <- df.temp2[(1+which(diff(df.temp2$age) <= 0)):nrow(df.temp2),(1+which(diff(df.temp2$age) <= 0)):(4+nrow(df.temp2))]
      }
      else{
        FSTdata_cov_ponds[[i]] <- df.temp2
      }
    }
  }
}

# remove also samples for false positives (only consisting of 1 or 2 samples)
false_positives <- c("FS.27B", "FS.31B", "N.71B")
FSTdata_cov_ponds <-  FSTdata_cov_ponds[!names( FSTdata_cov_ponds) %in% false_positives]

# data frame for absolute differences between samples
delta.frame <- data.frame()
delta.frame.connected <- data.frame()
delta.frame.isolated <- data.frame()
mean_dist_to_closest_2.check <- vector()

# calculate delta between all samples
for (i in 1:length(FSTdata_cov_ponds)) {
  if (is.data.frame(FSTdata_cov_ponds[[i]])){
    # delta age calculation
    age_matrix <- outer(as.numeric(unlist(FSTdata_cov_ponds[[i]]["year"])) + as.numeric(unlist(FSTdata_cov_ponds[[i]]["season"])),as.numeric(unlist(FSTdata_cov_ponds[[i]]["year"])) + as.numeric(unlist(FSTdata_cov_ponds[[i]]["season"])),custom_fun)
    
    # delta generations calculation
    gen_matrix <- outer(as.numeric(unlist(FSTdata_cov_ponds[[i]]["year"]))*8 + ifelse(as.numeric(unlist(FSTdata_cov_ponds[[i]]["season"]))==0, 0, 4),as.numeric(unlist(FSTdata_cov_ponds[[i]]["year"]))*8 + ifelse(as.numeric(unlist(FSTdata_cov_ponds[[i]]["season"]))==0, 0, 4),custom_fun)
    
    mean_dist_to_closest_2.check <- c(mean_dist_to_closest_2.check, mean(FSTdata_cov_ponds[[i]]$mean_dist_to_closest_2, na.rm = T))
    
    # add to data frame
    delta.frame <- rbind(delta.frame, cbind(age_matrix[lower.tri(age_matrix)], FSTdata_cov_ponds[[i]][1:(length(FSTdata_cov_ponds[[i]])-2)][lower.tri(FSTdata_cov_ponds[[i]][1:(length(FSTdata_cov_ponds[[i]])-2)])], gen_matrix[lower.tri(gen_matrix)], rep(mean(FSTdata_cov_ponds[[i]]$mean_dist_to_closest_2, na.rm = T), length(gen_matrix[lower.tri(gen_matrix)]))))
    
    if(all(is.na(FSTdata_cov_ponds[[i]]$mean_dist_to_closest_2))) next
    if(mean(FSTdata_cov_ponds[[i]]$mean_dist_to_closest_2, na.rm = T) >= 50){
      delta.frame.isolated <- rbind(delta.frame.isolated, cbind(age_matrix[lower.tri(age_matrix)], FSTdata_cov_ponds[[i]][1:(length(FSTdata_cov_ponds[[i]])-2)][lower.tri(FSTdata_cov_ponds[[i]][1:(length(FSTdata_cov_ponds[[i]])-2)])], gen_matrix[lower.tri(gen_matrix)]))
    } 
    if(mean(FSTdata_cov_ponds[[i]]$mean_dist_to_closest_2, na.rm = T) < 50){
      delta.frame.connected <- rbind(delta.frame.connected, cbind(age_matrix[lower.tri(age_matrix)], FSTdata_cov_ponds[[i]][1:(length(FSTdata_cov_ponds[[i]])-2)][lower.tri(FSTdata_cov_ponds[[i]][1:(length(FSTdata_cov_ponds[[i]])-2)])], gen_matrix[lower.tri(gen_matrix)]))
    }
  }
}


# plot pairwise FST vs. delta years
ggplot(delta.frame, aes(V1, V2)) + geom_point() + theme_bw() +
  stat_smooth(method='loess',se=T) + ylab("pairwise FST") + xlab("delta years")

# explore non-linearity (formula from Bergland 2014)
non.linear.model <- nls(V2 ~ a*b^(1/V1),
                        data = delta.frame,
                        start = list(a = 1, b = 1))

# 95CI for plotting
delta.frame4 <- delta.frame %>%
  mutate(group = ifelse(V4 > median(V4, na.rm=T), "High", "Low")) %>% drop_na(V4)

nsmodel1<-nls(V2 ~ a*b^(1/V1),data=subset(delta.frame4, group=="High"),start = list(a = 1, b = 1))
nsmodel2<-nls(V2 ~ a*b^(1/V1),data=subset(delta.frame4, group=="Low"), start = list(a = 1, b = 1))

#http://www.leg.ufpr.br/~walmes/cursoR/ciaeear/as.lm.R
# if object is an "nls" object then its often used like this:
# predict(as.lm(object), ...) where ... are any predict.lm args

# as.lm.nls effectively just does this:
# lm(lhs ~ gradient - 1, offset = fitted(object),
#   list(gradient = object$m$gradient(), lhs = object$m$lhs()))
# so most of the code is just to get the names right.

as.lm <- function(object, ...) UseMethod("as.lm")

as.lm.nls <- function(object, ...) {
  if (!inherits(object, "nls")) {
    w <- paste("expected object of class nls but got object of class:", 
               paste(class(object), collapse = " "))
    warning(w)
  }
  
  gradient <- object$m$gradient()
  if (is.null(colnames(gradient))) {
    colnames(gradient) <- names(object$m$getPars())
  }
  
  response.name <- if (length(formula(object)) == 2) "0" else 
    as.character(formula(object)[[2]])
  
  lhs <- object$m$lhs()
  L <- data.frame(lhs, gradient)
  names(L)[1] <- response.name
  
  fo <- sprintf("%s ~ %s - 1", response.name, 
                paste(colnames(gradient), collapse = "+"))
  fo <- as.formula(fo, env = proto:::as.proto.list(L))
  
  do.call("lm", list(fo, offset = substitute(fitted(object))))
  
}

fit1<-predict(as.lm.nls(nsmodel1), interval = "confidence") 
fit2<-predict(as.lm.nls(nsmodel2), interval = "confidence") 
delta.frame4$lowerfit[delta.frame4$group=="High"]<-fit1[,2]
delta.frame4$upperfit[delta.frame4$group=="High"]<-fit1[,3]
delta.frame4$lowerfit[delta.frame4$group=="Low"]<-fit2[,2]
delta.frame4$upperfit[delta.frame4$group=="Low"]<-fit2[,3]

# Figure with overall and subset
Fig4A <- delta.frame4 %>%
  mutate(group = ifelse(V4 > median(V4, na.rm=T), "High", "Low")) %>% drop_na(V4) %>%
  ggplot(aes(x = V1, y = V2, color = group)) +
  geom_point(alpha = 0.3) + theme_bw() +  coord_cartesian(ylim=c(0, 0.6)) + 
  labs(x=expression("delta years (="~paste(delta)~"t)"), y= expression("pairwise"~paste(italic(F) [ST])))+
  theme(axis.text=element_text(size=15), axis.title=element_text(size=15), legend.text = element_text(size=15), legend.title = element_text(size=15), legend.position = "bottom") +
  scale_color_manual(name="mean NN2:", values=c("#332288", "#DDCC77"), labels=c("top half (more isolated)", "bottom half (less isolated)"))+
  geom_ribbon(aes(x=V1,ymin=lowerfit,ymax=upperfit),data=subset(delta.frame4, group=="High"),alpha=0.3, colour = NA)+
  geom_ribbon(aes(x=V1,ymin=lowerfit,ymax=upperfit),data=subset(delta.frame4, group=="Low"),alpha=0.3, colour = NA)+
  stat_smooth(data=delta.frame, aes(V1, V2), method='nls', formula = y ~ a*b^(1/x), method.args = list(start = list(a = 1, b = 1)), se = FALSE, color = "black") + 
  annotate("text", x=0.75, y=0.56, label=expression(atop(paste(italic(F) [ST])~" = 0.13*0.51"^{"1/"~paste(delta)~"t"}, paste(italic("R")^2," = 0.089, ", italic("p")," = 6.43e-18"))), size = 5)+
  geom_smooth(method='nls', formula = y ~ a*b^(1/x), method.args = list(start = list(a = 1, b = 1)), se = F)



# ... and (delta) per year or season
set.seed(123)

# vector for delta spring
spring.vector <- vector()

for (i in 1:length(FSTdata_cov_ponds)) {
  if (is.data.frame(FSTdata_cov_ponds[[i]])){
    
    # same season
    FSTdata_cov_ponds.tmp <- FSTdata_cov_ponds[[i]][grep("spr",rownames(FSTdata_cov_ponds[[i]])), grep(paste(c("spr","year","season"), collapse = "|"),colnames(FSTdata_cov_ponds[[i]]))]
    
    # delta age calculation
    age_matrix <- outer(as.numeric(unlist(FSTdata_cov_ponds.tmp["year"])) + as.numeric(unlist(FSTdata_cov_ponds.tmp["season"])),as.numeric(unlist(FSTdata_cov_ponds.tmp["year"])) + as.numeric(unlist(FSTdata_cov_ponds.tmp["season"])),custom_fun)
    
    # diff 1 year
    FSTdata_cov_ponds.season <- FSTdata_cov_ponds.tmp[-c(ncol(FSTdata_cov_ponds.tmp)-1, ncol(FSTdata_cov_ponds.tmp))][age_matrix == 1 & lower.tri(age_matrix)]
    
    # add to vector
    #spring.vector <- c(spring.vector, mean(FSTdata_cov_ponds.season))
    #spring.vector <- c(spring.vector, FSTdata_cov_ponds.season)
    if (length(FSTdata_cov_ponds.season)!=0) {
      spring.vector <- c(spring.vector, sample(FSTdata_cov_ponds.season, 1))
    }else {spring.vector <- c(spring.vector,NA)}
  }
}

# vector for delta summer
summer.vector <- vector()

for (i in 1:length(FSTdata_cov_ponds)) {
  if (is.data.frame(FSTdata_cov_ponds[[i]])){
    
    # same season
    FSTdata_cov_ponds.tmp <- FSTdata_cov_ponds[[i]][grep("smr",rownames(FSTdata_cov_ponds[[i]])), grep(paste(c("smr","year","season"), collapse = "|"),colnames(FSTdata_cov_ponds[[i]]))]
    
    # delta age calculation
    age_matrix <- outer(as.numeric(unlist(FSTdata_cov_ponds.tmp["year"])) + as.numeric(unlist(FSTdata_cov_ponds.tmp["season"])),as.numeric(unlist(FSTdata_cov_ponds.tmp["year"])) + as.numeric(unlist(FSTdata_cov_ponds.tmp["season"])),custom_fun)
    
    # diff 1 year
    FSTdata_cov_ponds.season <- FSTdata_cov_ponds.tmp[-c(ncol(FSTdata_cov_ponds.tmp)-1, ncol(FSTdata_cov_ponds.tmp))][age_matrix == 1 & lower.tri(age_matrix)]
    
    # add to vector
    #summer.vector <- c(summer.vector, mean(FSTdata_cov_ponds.season))
    #summer.vector <- c(summer.vector, FSTdata_cov_ponds.season)
    if (length(FSTdata_cov_ponds.season)!=0) {
      summer.vector <- c(summer.vector, sample(FSTdata_cov_ponds.season, 1))
    }else {summer.vector <- c(summer.vector,NA)}
  }
}

# vector for delta active season
sprsmr.vector <- vector()

for (i in 1:length(FSTdata_cov_ponds)) {
  if (is.data.frame(FSTdata_cov_ponds[[i]])){
    
    # delta age calculation
    age_matrix <- outer(as.numeric(unlist(FSTdata_cov_ponds[[i]]["year"])) + as.numeric(unlist(FSTdata_cov_ponds[[i]]["season"])),as.numeric(unlist(FSTdata_cov_ponds[[i]]["year"])) + as.numeric(unlist(FSTdata_cov_ponds[[i]]["season"])),custom_fun)
    
    # diff 0.2 year
    FSTdata_cov_ponds.season <- FSTdata_cov_ponds[[i]][-seq(ncol(FSTdata_cov_ponds[[i]])-3, ncol(FSTdata_cov_ponds[[i]]))][age_matrix == "0.2" & lower.tri(age_matrix)]
    
    # add to vector
    #sprsmr.vector <- c(sprsmr.vector, mean(FSTdata_cov_ponds.season))
    #sprsmr.vector <- c(sprsmr.vector, FSTdata_cov_ponds.season)
    if (length(FSTdata_cov_ponds.season)!=0) {
      sprsmr.vector <- c(sprsmr.vector, sample(FSTdata_cov_ponds.season, 1))
    }else {sprsmr.vector <- c(sprsmr.vector,NA)}
  }
}

# vector for delta dormant season
smrspr.vector <- vector()

for (i in 1:length(FSTdata_cov_ponds)) {
  if (is.data.frame(FSTdata_cov_ponds[[i]])){
    # delta age calculation
    age_matrix <- outer(as.numeric(unlist(FSTdata_cov_ponds[[i]]["year"])) + as.numeric(unlist(FSTdata_cov_ponds[[i]]["season"])),as.numeric(unlist(FSTdata_cov_ponds[[i]]["year"])) + as.numeric(unlist(FSTdata_cov_ponds[[i]]["season"])),custom_fun)
    
    # diff 0.8 year
    FSTdata_cov_ponds.season <- FSTdata_cov_ponds[[i]][-seq(ncol(FSTdata_cov_ponds[[i]])-3, ncol(FSTdata_cov_ponds[[i]]))][age_matrix == "0.8" & lower.tri(age_matrix)]
    
    # add to vector
    #smrspr.vector <- c(smrspr.vector, mean(FSTdata_cov_ponds.season))
    #smrspr.vector <- c(smrspr.vector, FSTdata_cov_ponds.season)
    if (length(FSTdata_cov_ponds.season)!=0) {
      smrspr.vector <- c(smrspr.vector, sample(FSTdata_cov_ponds.season, 1))
    }else {smrspr.vector <- c(smrspr.vector,NA)}
  }
}

# pond names
names.vector <- vector()
for (i in 1:length(FSTdata_cov_ponds)) {
  if (is.data.frame(FSTdata_cov_ponds[[i]])){
    names.vector <- c(names.vector, names(FSTdata_cov_ponds[i]))
  }
}

# all delta
plot.frame <- rbind(data.frame(FST = spring.vector, comp = "sprspr", row.names = paste0(names.vector,"_sprspr")), 
                    data.frame(FST = summer.vector, comp = "smrsmr", row.names = paste0(names.vector,"_smrsmr")), 
                    data.frame(FST = sprsmr.vector, comp = "sprsmr", row.names = paste0(names.vector,"_sprsmr")), 
                    data.frame(FST = smrspr.vector, comp = "smrspr", row.names = paste0(names.vector,"_smrspr")))

# calculate significance for:
my_comparisons <- list( c("smrsmr", "sprsmr"), c("sprspr", "sprsmr"), c("smrspr", "sprsmr") )

# trafo
means_untrafo <- plot.frame %>%
  group_by(comp) %>%
  dplyr::summarize(Mean = mean(FST, na.rm=TRUE))
plot.frame$FST <- log((0.001 + plot.frame$FST)/(1 + 0.001 + plot.frame$FST))

Fig4B <- ggboxplot(plot.frame[plot.frame$comp %in% c("sprsmr", "smrspr"),], x = "comp", y = "FST", notch = F, xlab = "", ylab = "pairwise FST / (1 - pairwise FST) (ln)", main = "") + 
  add_pvalue(data.frame(group1 = "sprsmr",group2 = "smrspr",label = "0.0047",y.position = -0.9), xmin = "group1", xmax = "group2", label = "label", y.position = "y.position", label.size=5, bracket.size = 0.3) + # add pairwise comparisons p-value
  stat_summary(fun.data = function(x) data.frame(y=-7.75, label = round(mean(x), digits = 3)), geom="text", size = 5) +# add mean values
  geom_text(data=means_untrafo[means_untrafo$comp %in% c("sprsmr", "smrspr"),], aes(x=comp, y=-8.25, label=round(Mean, digits = 3)), size = 5, hjust = 0.38) +
  theme_bw() +
  theme(axis.text=element_text(size=15), axis.title=element_text(size=15), legend.text = element_text(size=15)) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2), labels = c("spring to summer", "summer to spring")) +
  ylim(c(-8.4, -0.6))


# 3) analysing genomic diversity (greendalf output) 

# grenedalf command (data loaded above)
#grenedalf diversity --window-width 100000 --pool-sizes 50 --min-coverage 10 --vcf-path Dmagna_2014-2018_filtered.vcf.gz

# models to explain genomic diversity
anova(lmerTest::lmer(tetha_rel ~ relage + mean_dist_to_closest_2 + (1|island_group) + (1|island_group:pop), long_series.frame))
anova(lmerTest::lmer(tetha_rel ~ relage + mean_dist_to_closest_2 + relage*mean_dist_to_closest_2 + (1|island_group) + (1|island_group:pop), long_series.frame))


# 4) PCA and t-SNE analysis

# transpose data for pcadapt
data_t <- transpose(as.data.frame(data))

# convert to pcadapt object
data_pcadapt <- read.pcadapt(data_t, type = "pool")

# PCA
x<-pcadapt(
  input=data_pcadapt,
  K = (nrow(data_pcadapt) - 1),
  method = "mahalanobis",
  min.maf = 0.05,
  ploidy = NULL,
  LD.clumping = NULL,
  pca.only = FALSE,
  tol
)

# t-SNE
tsne_out <- Rtsne(as.matrix(data_t),verbose=TRUE, max_iter = 5000, initial_dims= 9, check_duplicates = FALSE)
