# readRDS(file = "rime.data.rds"); readRDS("rime.data.metadata.rds")


dat.tmp <- dat.tmp[ ,grep("Abundances",colnames(dat.tmp))]
# identify all missing
index.tmp <- rowSums(is.na(dat.tmp[ , ])) == ncol(dat.tmp)
# remove all missing
dat.tmp  <- dat.tmp[ !index.tmp , ]
dat.tmp$Accession <- row.names(dat.tmp)
###
dat <- dat.tmp


################################################################################
#####################################################################
##################################################
################################################## create back up of non.normalised dataset
##################################################
#####################################################################
################################################################################
dat.no.norm <-  dat
row.names(dat.no.norm) <- paste( dat.no.norm$Accession )
dat.no.norm$Accession <- NULL
colnames(dat.no.norm) <- sub("Abundances..Grouped...F1..", "", colnames(dat.no.norm))
colnames(dat.no.norm) <-  sub(metadata$File.name[1], "", colnames(dat.no.norm))
colnames(dat.no.norm) <-  sub("_", "", colnames(dat.no.norm))


################################################################################
#####################################################################
##################################################
################################################## START - Q plex rime analysis
##################################################
#####################################################################
################################################################################
# Q plex rime analysis 

library(qPLEXanalyzer); library(gridExtra)
# prep expression and metadata files
dat$Accession <- NULL
row.names(metadata) <- metadata$Column.name
# make sure order of samples is the same in metadata file and expression matrix
idx <- match(colnames(dat), row.names(metadata))
idx
metadata <- metadata[ idx, ]
all(colnames(dat) == row.names(metadata)) # !!! should be TRUE
colnames(dat) <- make.unique(paste(metadata$Sample.ID, metadata$Biological.replicate, sep = "."), sep = ".")
row.names(metadata) <- make.unique(paste(metadata$Sample.ID, metadata$Biological.replicate, sep = ".") )


# rename metadata columns for convenience
dat$Accession <- row.names(dat)
metadata$SampleName <- row.names(metadata)
metadata$SampleGroup <- metadata$Sample.ID
metadata$BioRep <- metadata$Biological.replicate
metadata$TechRep <- metadata$Technical.replicate
metadata$Grp <- metadata$Antibody


# convert to MSnset, required by qPLEXanalyzer
MSnset_data <- convertToMSnset(dat,
                               metadata = metadata,
                               indExpData = c(1:(ncol(dat)-1) ), 
                               Sequences = ncol(dat), 
                               Accessions = ncol(dat))
# check if all went OK
MSnset_data

# diagnostic and QC plots 
# standard plots we normally look at
intensityPlot(MSnset_data, title = "Protein intensity distribution")
intensityBoxplot(MSnset_data, title = "Protein intensity distribution")
rliPlot(MSnset_data, title = "Relative Protein intensity")
corrPlot(MSnset_data)
exprs(MSnset_data) <- exprs(MSnset_data) + 0.01
hierarchicalPlot(MSnset_data)
pcaPlot(MSnset_data, labelColumn = "BioRep", pointsize = 3)


# compare normalization methods
p1 <- intensityPlot(MSnset_data, title = "No normalization")
p1
MSnset_norm_q <- normalizeQuantiles(MSnset_data)
p2 <- intensityPlot(MSnset_norm_q, title = "Quantile")
p2
MSnset_norm_ns <- normalizeScaling(MSnset_data, scalingFunction = median)
p3 <- intensityPlot(MSnset_norm_ns, title = "Scaling")
p3
MSnset_norm_gs <- groupScaling(MSnset_data, 
                               scalingFunction = median, 
                               groupingColumn = "SampleGroup")
p4 <- intensityPlot(MSnset_norm_gs, title = "Within Group Scaling")
p4 

grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2)



# run diff expression stats steps - this uses limma in the background
contrasts <- c(CellLine.1 = "CellLine.1_KO - CellLine.1_WT") 
diffstats <- computeDiffStats(MSnset_norm_gs, contrasts = contrasts)
diffexp <- getContrastResults(diffstats, 
                              contrast = contrasts,
                              controlGroup = "IgG.CellLine.1_WT")

maVolPlot(diffstats, contrast = contrasts, plotType = "MA", title = contrasts)
maVolPlot(diffstats, contrast = contrasts, plotType = "Volcano", title = contrasts)

################################################################################
#####################################################################
##################################################
################################################## END Q plex rime analysis
##################################################
#####################################################################
################################################################################




################################################################################
#####################################################################
##################################################
################################################## START custom analysis using code modules from Pxanalytics/PASTA package
##################################################
#####################################################################
################################################################################
set.seed(5)

# Transform the data to the long format
dat.long <- dat.tmp %>%
  pivot_longer(cols = starts_with("Abundance"),
               names_to = "TMT.channel",
               values_to =  "Abundances") %>%
  left_join(metadata, by = c("TMT.channel" = "Column.name")) %>%
  mutate(TMT.channel = gsub("Abundance.", "", gsub("Abundances..Grouped...", "", TMT.channel))) %>%
  separate(TMT.channel, into = c("TMT", "Channel", NA),remove = FALSE) %>%
  #  mutate(TMT = paste("F", Batch, sep = "")) %>%
  mutate(TMT.channel = paste(TMT, Channel, sep = "."))

# look at missing data pattern in more detail
d.na <- dat.long %>%
  group_by(TMT.channel) %>%
  summarise(No.obs = n(), No.missing_values = length(which(is.na(Abundances))), Sample.ID = unique(Sample.ID))
datatable(d.na)
d.na$No.obs_values <- d.na$No.obs - d.na$No.missing_values
# 
ggplot(data=d.na, aes(x=No.missing_values, y=paste(TMT.channel, Sample.ID ) )) +
  geom_bar(stat="identity") + 
  labs(title = "") +
  theme_bw()

d.na$TMT.channel <- sub("FNA.", "", d.na$TMT.channel)
ggplot(data=d.na, aes(x=No.missing_values, y=paste(TMT.channel ) )) +
  geom_bar(stat="identity") + 
  labs(title = "") +
  theme_bw()

ggplot(data=d.na, aes(x=No.missing_values, y=Sample.ID )) +
  geom_bar(stat="identity") + 
  labs(title = "") +
  theme_bw()

Amelia::missmap(dat.tmp[ , grep("Accession", colnames(dat.tmp), invert = T)], rank.order = FALSE, main = "", xlab = "TMT.Channel", ylab = "Protein.index")



##### visualize distributions 
# this can be done in many different ways, you need to know what are the factors in the exp. etc.
dat.long <- dat.long[ !is.na(dat.long$Cell.line) , ]
ggplot(dat.long, aes(x=log2(Abundances), color=Sample.ID , linetype  = as.factor(Crosslink)))+ 
  #  facet_wrap(~TMT) + 
  stat_density(geom = "line",position = "identity", size=1) + 
  labs(title = "Density plot for raw data") + 
  theme_classic()

#### additional plots
ggplot(data = dat.long, aes(x = TMT, y = log2(Abundances))) +
  facet_wrap(~Channel, nrow = 2) + 
  geom_boxplot() +
  labs(title = "Box plot for raw data", y = "log2(Abundances)") +
  theme_bw()+ 
  theme(axis.text.x = element_text(angle = 75, hjust = 0.1, vjust = -0.01)) 


ggplot(dat.long, aes(x=log2(Abundances),  color = Sample.ID, linetype  = as.factor(Biological.replicate) )   )+ 
  facet_wrap(~Crosslink) + 
  stat_density(geom = "line",position = "identity" ) + 
  labs(title = "Density plot for raw data") + 
  theme_classic()
###
ggplot(dat.long, aes(x=log2(Abundances), color=TMT))+ 
  facet_wrap(~Channel) + 
  stat_density(geom = "line",position = "identity") + 
  labs(title = "Density plot for raw data") + 
  theme_classic()

###########



###
# Normalization routines
##
#####
## Sample loading normalization
#####
# This assumes the same amount of protein labelled in each sample and the total signal in each channel summing to the same value within the batch. 
# The normalization can be done across the entire experiment, within batches, or within selected conditions (e.g. separately for Crosslinked and non-Crosslinked data) - user has to decide what is appropriate
dat.signal.stats <- dat.long %>%
  group_by(Crosslink, Antibody, Batch, Sample.ID, Biological.replicate) %>%
  summarise (SumEachChannel = sum(Abundances, na.rm = TRUE)) %>%
  group_by(Crosslink, Antibody, Batch) %>%
  mutate(MeanofSumWithinBatch = mean(SumEachChannel, na.rm = TRUE),
         MedianofSumWithinBatch = median(SumEachChannel, na.rm = TRUE),
         correct.factor.sl = MeanofSumWithinBatch/SumEachChannel)

# inspect signal across samples, can be useful to find issues with samples/sample groups
barplot(dat.signal.stats$correct.factor.sl, names.arg = dat.signal.stats$Sample.ID, las = 2)

# apply SL normalization 
dat.sl <- dat.long %>%
  group_by(Crosslink, Antibody, Batch, Sample.ID, Biological.replicate) %>%
  summarise (SumEachChannel = sum(Abundances, na.rm = TRUE))%>%
  group_by(Crosslink, Antibody, Batch) %>%
  mutate(MeanofSumWithinBatch = median(SumEachChannel, na.rm = TRUE),
         correct.factor.sl = MeanofSumWithinBatch/SumEachChannel) %>%
  right_join(dat.long, by = c("Crosslink", "Antibody", "Batch", "Sample.ID", "Biological.replicate")) %>%
  mutate(Abundances.sl = Abundances * correct.factor.sl) 
###


ggplot(dat.sl, aes(x=log2(Abundances.sl), color=Sample.ID , linetype  = as.factor(Crosslink)))+ 
  #  facet_wrap(~TMT) + 
  stat_density(geom = "line",position = "identity", size=1) + 
  labs(title = "Density plot for SL normalized data") + 
  theme_classic()


# plot output - additional plots 
ggplot(data = dat.sl, aes(x = Channel, y = log2(Abundances.sl))) +
  facet_wrap(~TMT, nrow = 2) + 
  geom_boxplot() +
  labs(title = "Box plot for SL data", y = "log2(Abundances)") +
  theme_bw()+ 
  theme(axis.text.x = element_text(angle = 75, hjust = 0.1, vjust = -0.01)) 


ggplot(data = dat.sl, aes(x = TMT, y = log2(Abundances.sl))) +
  facet_wrap(~Channel, nrow = 2) + 
  geom_boxplot() +
  labs(title = "Box plot for SL data", y = "log2(Abundances)") +
  theme_bw()+ 
  theme(axis.text.x = element_text(angle = 75, hjust = 0.1, vjust = -0.01)) 
###
ggplot(dat.sl, aes(x=log2(Abundances.sl), color=Channel))+ 
  facet_wrap(~TMT) + 
  stat_density(geom = "line",position = "identity") + 
  labs(title = "Density plot for SL data", x = "log2(Abundances)") + 
  theme_classic()
######

#####
ggplot(dat.sl, aes(x=log2(Abundances.sl), color=TMT))+ 
  facet_wrap(~Channel) + 
  stat_density(geom = "line",position = "identity") + 
  labs(title = "Density plot for SL data", x = "log2(Abundances)") + 
  theme_classic()
####

#########
#####
## Batch effect correction
#####
# The median correction of common control - "Median.CC" (default) or the internal reference scaling - "IRS" can be selected. 
Batch.correction <- function(data, method = "Median.CC", CC.channel = "126"){
  if (method == "Median.CC"){
    dat.correct <- data %>%
      filter(Channel == CC.channel)%>%
      group_by(Accession) %>%
      na_if(0) %>%
      mutate(correct.across.batch = 2^(median(log2(Abundances.sl), na.rm = TRUE)),
             correct.factor.batch = correct.across.batch/Abundances.sl) %>%
      
      dplyr::select(Accession, TMT, correct.factor.batch) %>%
      left_join(data, by = c("Accession", "TMT"))%>%
      mutate(Abundances.correct = Abundances.sl * correct.factor.batch)
  } else {
    dat.correct <- data %>%
      group_by(Accession, TMT) %>%
      summarise(sum.channel = sum(Abundances.sl, na.rm = TRUE)) %>%
      group_by(Accession) %>%
      na_if(0) %>%
      mutate(average.of.sum = 2^(mean(log2(sum.channel), na.rm = TRUE)),
             correct.factor.batch = average.of.sum/sum.channel) %>%
      left_join(data, by = c("Accession", "TMT"))%>%
      mutate(Abundances.correct = Abundances.sl * correct.factor.batch)
  }
  return(dat.correct)
}
######
# 
# # Apply batch correction function
dat.correct <- Batch.correction(data = dat.sl, method = "Median.CC", CC.channel = "126")


# # plot output
ggplot(dat.correct, aes(x=log2(Abundances.correct),  color = Sample.ID, linetype  = as.factor(Biological.replicate) )   )+ 
  facet_wrap(~Crosslink) + 
  stat_density(geom = "line",position = "identity" ) + 
  labs(title = "Density plot for SL/Batch correction (Median.CC) data") + 
  theme_classic()

ggplot(data = dat.correct, aes(x = Channel, y = log2(Abundances.correct))) +
  facet_wrap(~TMT, nrow = 2) +
  geom_boxplot() +
  labs(title = "Box plot for SL/Batch correction (Median.CC) data", y = "log2(Abundances)") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 75, hjust = 0.1, vjust = -0.01))

ggplot(dat.correct, aes(x=log2(Abundances.correct), color=Channel))+
  facet_wrap(~TMT) +
  stat_density(geom = "line",position = "identity") +
  labs(title = "Density plot for SL/Batch correction (Median.CC) data", x = "log2(Abundances)") +
  theme_classic()
#######
ggplot(dat.correct, aes(x=log2(Abundances.correct), color=TMT))+
  facet_wrap(~Channel) +
  stat_density(geom = "line",position = "identity") +
  labs(title = "Density plot for SL/Batch correction (Median.CC) data", x = "log2(Abundances)") +
  theme_classic()
#######

#####
## TMM (trimmed mean of M values) normalization final step
#####
# To ensure the median across TMM channels are as similar as possible 
dat.tmm <- dat.correct %>%
  pivot_wider(names_from = Accession, id_cols = c("TMT", "Channel"), 
              values_from = Abundances.correct)

tmm_raw <- list()# c()

for (i in unique(dat.tmm$TMT)){
  d.tmm <- filter(dat.tmm, TMT == i) %>%
    select(-TMT)%>%
    column_to_rownames(var = "Channel") 
  
  temp <- as.data.frame(calcNormFactors(na.omit(t(d.tmm)), method = "TMM"))# , refColumn = "126" 
  # Change refColum to the channel for common control
  colnames(temp) <- i
  
  tmm_raw[[i]] <- as.data.frame(t(temp))   
  
  # tmm_raw <- rbind(tmm_raw, t(temp))
}

dat.tmm <- plyr::rbind.fill(tmm_raw) %>%
  `rownames<-`(unlist(lapply(tmm_raw, row.names))) %>% 
  rownames_to_column(var = "TMT") %>%
  reshape2::melt(id.var = "TMT", variable.name = "Channel", value.name = "tmm_factor") %>%
  right_join(dat.correct, by = c("TMT", "Channel")) %>%
  mutate(Abundances.tmm= Abundances.correct/tmm_factor)
##
####
######
ggplot(dat.tmm, aes(x=log2(Abundances.tmm), color=Sample.ID , linetype  = as.factor(Crosslink)))+ 
  #  facet_wrap(~TMT) + 
  stat_density(geom = "line",position = "identity", size=1) + 
  labs(title = "Density plot for SL/TMM normalized data") + 
  theme_classic()


## additional plots
ggplot(data = dat.tmm, aes(x = Channel, y = log2(Abundances.tmm))) +
  facet_wrap(~TMT, nrow = 2) + 
  geom_boxplot() +
  labs(title = "Boxplot for SL/Batch correction (Median.CC)/TMM data", y = "log2(Abundances)") +
  theme_bw()+ 
  theme(axis.text.x = element_text(angle = 75, hjust = 0.1, vjust = -0.01)) 


ggplot(dat.tmm, aes(x=log2(Abundances.tmm), color=TMT))+ 
  facet_wrap(~Channel) + 
  stat_density(geom = "line",position = "identity") + 
  labs(title = "Density plot for SL/Batch correction (Median.CC)/TMM data", x = "log2(Abundances)") + 
  theme_classic()

ggplot(dat.tmm, aes(x=log2(Abundances.tmm), color=Channel))+ 
  facet_wrap(~TMT) + 
  stat_density(geom = "line",position = "identity") + 
  labs(title = "Density plot for SL/Batch correction (Median.CC)/TMM data", x = "log2(Abundances)") + 
  theme_classic()

ggplot(dat.tmm, aes(x=log2(Abundances.tmm), color=Sample.ID))+ 
  facet_wrap(~TMT) + 
  stat_density(geom = "line",position = "identity") + 
  labs(title = "Density plot for SL/Batch correction (Median.CC)/TMM data", x = "log2(Abundances)") + 
  theme_classic()


ggplot(dat.tmm, aes(x=log2(Abundances.tmm), color=Sample.ID))+ 
  facet_wrap(~Crosslink) + 
  stat_density(geom = "line",position = "identity") + 
  labs(title = "Density plot for SL/Batch correction (Median.CC)/TMM data", x = "log2(Abundances)") + 
  theme_classic()






# convert to data frame in wide format
dat.normalise <- dat.tmm %>%
  select(Accession,  TMT.channel, Abundances.tmm) %>%
  pivot_wider(id_cols = c("Accession"), names_from = "TMT.channel", values_from = "Abundances.tmm")

dat <- as.data.frame(dat.normalise)
row.names(dat) <- paste( dat$Accession )
dat$Accession <- NULL
# dat.normalise.sl.tmm <- dat



metadata <- metadata %>% 
  mutate(TMT.channel = gsub("Abundance.", "", gsub("Abundances..Grouped...", "", Column.name))) %>% 
  mutate(TMT.channel = gsub("_PE.*", "", TMT.channel)) %>% 
  mutate(TMT.channel = gsub("\\.\\.", ".", TMT.channel))
row.names(metadata) <- metadata$TMT.channel
####
colnames(dat)
row.names(dat)
# dat$F1.Accession <- NULL
idx <- match(colnames(dat), row.names(metadata))
idx
metadata <- metadata[ idx, ]
####
all(colnames(dat) == row.names(metadata))
######
# dat.main <- dat
# metadata.main <-metadata
# 
#########################
#########################
##
## limma analysis
###############################
metadata$Sample.ID <- factor(metadata$Sample.ID)

# construct model matrix
mod.matrix1 <- model.matrix( ~ 0 + Sample.ID,
                             data = metadata) # + Factor.5_concentration.ug.kg.day
colnames(mod.matrix1) <- sub("Sample.ID", "", colnames(mod.matrix1))
#########################
#########################
##
## run model
##
#########################
#########################
# Fit the model for the treatment effects to obtain the pooled variance
fit <- lmFit(log2(dat), mod.matrix1)
quantile(fit$sigma, c(.32, .57, 0.8, 0.9,.98)) 
fit.main <- as.data.frame(fit)
fit.main$Accession <- row.names(fit)

contr.mat <- makeContrasts(
  CellLine.1_KO - CellLine.1_WT, 
  CellLine.2_KO - CellLine.2_WT,
  (CellLine.2_KO - CellLine.2_WT) - (CellLine.1_KO - CellLine.1_WT),
  
  IgG.CellLine.1_WT - CellLine.1_WT,
  IgG.CellLine.1_KO - CellLine.1_KO,
  IgG.CellLine.2_WT - CellLine.2_WT,
  IgG.CellLine.2_KO -  CellLine.2_KO, 
  
  (IgG.CellLine.1_WT+IgG.CellLine.1_KO) - (CellLine.1_WT+CellLine.1_KO),
  
  (IgG.CellLine.2_WT+IgG.CellLine.2_KO) - (CellLine.2_WT+CellLine.2_KO),
  levels = mod.matrix1 )

colnames(mod.matrix1)
fit2 <-contrasts.fit(fit, contr.mat)		
fit2 <- eBayes(fit2)
fit2.main.2 <- as.data.frame(fit2)
fit2.main.2$Accession <- row.names(fit2)

### subset data for samples of interest
dat <- dat.main[ , which(metadata.main$Crosslink == "Crosslink:yes" & metadata.main$Antibody != "IgG")  ]
metadata <- metadata.main[ which(metadata.main$Crosslink == "Crosslink:yes" & metadata.main$Antibody != "IgG") , ]
all(colnames(dat) == row.names(metadata))

# construct model matrix
mod.matrix1 <- model.matrix( ~ 0 + Sample.ID,
                             data = metadata) # + Factor.5_concentration.ug.kg.day
colnames(mod.matrix1) <- sub("Sample.ID", "", colnames(mod.matrix1))
#########################
#########################
##
## run model
##
#########################
#########################
# Fit the model for the treatment effects to obtain the pooled variance
fit <- lmFit(log2(dat), mod.matrix1)
quantile(fit$sigma, c(.32, .57, 0.8, 0.9,.98)) 

contr.mat <- makeContrasts(
  CellLine.1_KO - CellLine.1_WT, 
  CellLine.2_KO - CellLine.2_WT,
  (CellLine.2_KO - CellLine.2_WT) - (CellLine.1_KO - CellLine.1_WT),
  
  
  levels = mod.matrix1 )

fit2 <-contrasts.fit(fit, contr.mat)		
fit2 <- eBayes(fit2)
hist((fit2$sigma) )
mean((fit2$sigma) )
############
result <- decideTests(fit2, method="separate",adjust.method="BH",p.value=0.05,lfc=0)
summary(decideTests(fit2, method="separate",adjust.method="BH",p.value=0.05,lfc=0))

# write.fit(fit2, results=NULL, "fit2.tmt3.txt", digits=3, adjust="BH", method="separate", F.adjust="none", sep="\t")


de.common <- which(row.names(fit2) %in% c("O96028_NSD2", "P51124_GZMM", "P09012_SNRPA"))
pheatmap::pheatmap(fit2$coefficients[ (de.common), ], scale = "none", cluster_rows = T, cluster_cols = F )


# explore individual coef./comparisons

n.coef <- 1
tab <- topTable(fit2, coef = n.coef, adjust="BH", number=Inf, sort.by = NULL, resort.by = NULL)
tab <- data.frame(tab) %>%  rownames_to_column(var = "Accession")

p1 <- ggplot(data=tab, aes(x=logFC, y=-log10(adj.P.Val))) +
  geom_point(alpha=0.5, col = AZcolor$Graphite, pch = 21, cex = 1.8) +
  
  geom_point(data= filter(tab, adj.P.Val < 0.05 & logFC < 0), #& Accession %in% peptides.filter
             aes(x=logFC, y=-log10(adj.P.Val)), col="black", pch = 23, fill = AZcolor$Mulberry, alpha=0.8, cex = 2)+
 
 
  geom_hline(aes(yintercept = -log10(0.05)), linetype = 6, col = "black") +
  
  labs(x= expression( ~ log[2] ~ FC), ylab = expression( ~ -log[10] ~ (Adjusted.p.value)),        title = paste(colnames(fit2$coefficients)[n.coef], "")  ) +
  # label some of the points
  geom_text_repel(aes(label=ifelse(adj.P.Val < 0.0000000000005 & logFC   < 0  &  Accession  != "O96028_NSD2" ,as.character(Accession),'')),hjust=0,vjust=0, cex = 3.5, box.padding = 0.7, max.overlaps = Inf, fontface = "bold")  +
  geom_text_repel(aes(label= ifelse(Accession  == "O96028_NSD2",  as.character(Accession),'')),hjust=0,vjust=0, cex = 3.5, col = AZcolor$Gold, box.padding = 0.7, max.overlaps = Inf, fontface = "bold")

p1