#Open R (v. 4.2.1) and set working directory 
set.seed(1000)

getwd()
setwd()

###### Part one: process the 16S rRNA amplicon sequences 

# Load dada2 package 

library(dada2) # (v.1.16.0)

# 1. Quality check, filter and trim the fastqs

# Define the path containing the fastq files 		

path <- ("/home/Final/Fastqs/BALFs")

# Define the samples´ names 

fnFs <- sort(list.files("/home/Final/Fastqs/BALFs", pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files("/home/Final/Fastqs/BALFs", pattern="_R2_001.fastq" , full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#inspect read quality profiles to define best suited filtering parameters

pdf("read_Quality.pdf", height=7, width=24)
par(mfrow=c(1,2))
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])
dev.off()

pdf("read_Quality2.pdf", height=7, width=24)
par(mfrow=c(1,2))
plotQualityProfile(fnFs[3:4])
plotQualityProfile(fnRs[3:4])
dev.off()

pdf("read_Quality3.pdf", height=7, width=30)
par(mfrow=c(1,2))
plotQualityProfile(fnFs[5:7])
plotQualityProfile(fnRs[5:7])
dev.off()


# Filter and trim and forward and reverse read 
# Place filtered fastqs in a sub-directory named filtered

filtFs <- file.path("/home/Final/Fastqs/BALFs", "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path("/home/Final/Fastqs/BALFs", "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(270,250),maxN=0, maxEE=c(3,3), rm.phix=TRUE,compress=TRUE, multithread=T)
head(out)
write.csv(data.frame(out), "filtered.csv")

file.exists(filtFs) 
file.exists(filtRs)

# Learn the Error Rates
errF <- learnErrors(filtFs, MAX_CONSIST=100,  randomize=TRUE,  multithread=T)
errR <- learnErrors(filtRs,MAX_CONSIST=100,  randomize=TRUE,  multithread=T)


## visualize the estimated error rates: 
pdf("Estimated_error_rates.pdf", height=18, width=20)
par(mfrow=c(1,2))
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)
dev.off()


pdf("Estimated_error_rates2.pdf", height=18, width=20)
par(mfrow=c(1,2))
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)
dev.off()


# Dereplication

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

#Denoising sample inference 

dadaFs <- dada(derepFs, err=errF, pool=FALSE, selfConsist=T,multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, pool=FALSE, selfConsist=T,multithread=TRUE)

dadaFs[[1]]
dadaRs[[1]]

# merging forward and reverse reads 

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE, minOverlap = 20) 
seqtab <- makeSequenceTable(mergers)
dim(seqtab) 


## remove chimeras with the option which best suits the dataset 

seqtab.nochim <- removeBimeraDenovo(seqtab, "per-sample", multithread=T, verbose=TRUE)

dim(seqtab.nochim) 
sum(seqtab.nochim)*100/sum(seqtab)   # Percenatge of non-chimeras in the merged sequences reads. 

# track the number of reads/sample through the steps 

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))

colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
write.csv(track, "track.csv")


#2. Assign taxonomy up to the genus level to "seq.nochim" chimera free merged using RDP classifier version 18

taxa.minboot80 <- assignTaxonomy(seqtab.nochim, minBoot=80, tryRC=T,  "./rdp_train_set_18.fa.gz", multithread=F, verbose=T, outputBootstraps = T)

#3. define ASVs which did not merge properly, and ASVs which did not classify as Bacteria 

# wrong length +/- 10% of 465 bp 

summary(nchar(getSequences(seqtab.nochim)))

wrong.length = which(nchar(getSequences(seqtab.nochim)) <= 418) 

# rename seqtab.nochim as "asvs.table", extract ASVs sequences, assign new ASVs names, 

asvs.table = seqtab.nochim

asvs.names <- vector()

for (i in 1:ncol(asvs.table)){
  a <- paste ("ASV_",i, sep="")
  asvs.names <- c(asvs.names,a)
}

asvs.names.data <- data.frame(cbind(asvs.names),colnames(asvs.table))
colnames(asvs.names.data)  <- c("ASV_Name", "ASV_Seq")
write.csv(asvs.names.data, "asvs.names.data.csv")

colnames(asvs.table) = asvs.names


# define a "taxonomy" data.frame as for asvs.table data.frame

taxonomy = taxa.minboot80$tax 
rownames(taxonomy) = asvs.names

taxonomy = data.frame(taxonomy) 


# define ASVs which did not classified as Bacteria 

wrong.tax = which(taxonomy$Kingdom != "Bacteria")

# remove wrongly merged ASVs, and non Bacteria ASVs 

asvs.table = asvs.table[,-c(wrong.length,wrong.tax)]

taxonomy = taxonomy[-c(wrong.length,wrong.tax),]

write.csv(asvs.table, "asvs.table.csv")
write.csv(taxonomy, "taxonomy.csv")


# generate fasta file of the ASVs sequences 


final.asvs.fasta <- data.frame()

for (i in 1:nrow(asvs.names.data)){
  a <- paste(">", asvs.names.data[i,]$ASV_Name, sep="")
  b <- rbind(a,asvs.names.data[i,]$ASV_Seq)
  final.asvs.fasta <- rbind(final.asvs.fasta,b)
}

write.table(final.asvs.fasta, "asvs.fasta", sep="\t")


#4. summary statistics of the defined ASVs 

# explore the summary statistics ASVs abundances 

summary.stat.asvs.table <- data.frame()

for (i in 1:nrow(asvs.table)){
  vec.summary <- vector()
  #detected.asvs <- which(colSums(asvs.table[i,]) > 0)
  detected.asvs <- names(which(asvs.table[i,] > 0))
  min_i <-  min(as.numeric(asvs.table[i,detected.asvs]))
  qt5_i  <- quantile(as.numeric(asvs.table[i,detected.asvs]), probs=0.05)
  qt10_i <- quantile(as.numeric(asvs.table[i,detected.asvs]), probs=0.10)
  qt25_i <- quantile(as.numeric(asvs.table[i,detected.asvs]))[2]
  qt50_i <- quantile(as.numeric(asvs.table[i,detected.asvs]))[3] 
  qt75_i <- quantile(as.numeric(asvs.table[i,detected.asvs]))[4]
  max_i  <- max(as.numeric(asvs.table[i,detected.asvs]))
  median_i <- median(as.numeric(asvs.table[i,detected.asvs]))
  mean_i <- mean(as.numeric(asvs.table[i,detected.asvs]))
  sd_i   <- sd(as.numeric(asvs.table[i,detected.asvs]))
  sum_i <- sum(as.numeric(asvs.table[i,detected.asvs]))
  max_ASV <- which(as.numeric(asvs.table[i,])==max_i)
  vec.summary <- cbind(length(detected.asvs),min_i,qt5_i, qt10_i, qt25_i,qt50_i,qt75_i,max_i, median_i,mean_i,sd_i, sum_i,colnames(asvs.table)[max_ASV][1])
  summary.stat.asvs.table <- rbind(summary.stat.asvs.table,vec.summary)
}


str(summary.stat.asvs.table)
rownames(summary.stat.asvs.table) <- rownames(asvs.table)
colnames(summary.stat.asvs.table) <- c("Detected_ASVs","Min", "5_quantile", "10_quantile", "25_quantile", "50_quantile", "75_quantile", "Max","Median","Mean","STD", "ASVs_Sum", "max_ASV")
write.csv(summary.stat.asvs.table, "summary.stat.asvs.table.csv")

# inspect distribution of ASVs across samples 

asvs.distribution = data.frame()

for(i in 1:ncol(asvs.table)){
  vec = vector()
  a = length(which(asvs.table[,i]>0))
  b = a *100/nrow(asvs.table)
  vec = c(a,b)
  asvs.distribution = rbind(asvs.distribution,vec)
}

colnames(asvs.distribution) = c("ASV_Freq", "ASV_Freq_Per")
rownames(asvs.distribution) = colnames(asvs.table)

write.csv(asvs.distribution, "asvs.distribution.csv")


### Part two: Adjusting ASVS clustering 

## Step 1 : Alignement of the ASVs sequences, and calculation of a pairwise distance matrix 

library(DECIPHER)

seqs <- asvs.names.data$ASV_Seq # ASVs´sequences 
length(seqs)
names(seqs) = asvs.names.data$ASV_Name

# align ASVs sequences, and derive a distance matrix 
asvs.alignment <- DECIPHER::AlignSeqs(DNAStringSet(seqs), verbose = T,processors = 4,iterations=10000,refinements=1000)

dm.asvs = DECIPHER::DistanceMatrix(asvs.alignment,type="matrix", correction = "none", processors = 4, verbose = T,includeTerminalGaps=T)


## make a diagonal distance matrix, to reduce unnecessary computing time

asvs = colnames(asvs.table)
dm.asvs.final = data.frame(dm.asvs[asvs,asvs])

dm.asvs.diag = data.frame(nrow=nrow(dm.asvs.final), ncol=ncol(dm.asvs.final))
vec = vector()

for(i in 1:ncol(dm.asvs.final)){
  
  if (i !=1){
    
    a = dm.asvs.final[,i][i:length(dm.asvs.final[,i])]
    
  } else {
    a = dm.asvs.final[,i]
  }
  vec = c(rep(NA, ncol(dm.asvs.final) - length(a)),a) 
  dm.asvs.diag = cbind(dm.asvs.diag,vec)
}


dim(dm.asvs.diag)
dm.asvs.diag = dm.asvs.diag[,-c(1,2)]
colnames(dm.asvs.diag) = colnames(dm.asvs.final)
rownames(dm.asvs.diag) = rownames(dm.asvs.final)
rownames(dm.asvs.diag)[1] = "ASV_1"

## Step 2: Define the range of ASVs pairwise-distances detected within pure culture samples (positive controls sequenced along): pure cultures, and their technical replicates (n=21), assessed distinctively 

## load the metadata file  "metadata.all.final"

## inspect the summary statistics of the pairwise-distances of true ASVs detected across the pure culture samples 

pure.cultures = rownames(metadata.all.final[which(metadata.all.final$Category=="Positive_control"),])[1:21]

pure.cultures.stats = vector()

for (i in 1:length(pure.cultures)){ ## assess each pure culture separately 
  vec = vector()
  a = names(which(asvs.table[pure.cultures[i],] > 0))  ## asvs detected within the pure culture i 
  
  if (length(unique(taxonomy[a,]$Genus)) > 1 ){
    #print(pure.cultures[i])
    #print(pure.cultures[i])
    
    a = a[which(taxonomy[a,]$Genus == names(sort(table(taxonomy[a,]$Genus), decreasing=T)[1]))]  ## excluding contaminants  so that only true ASV residents 
    
  }
  
  for (j in 1:length(a)){
    
    if (a[j] != a[length(a)]) { ## not for the last asv ib the vector a 
      #b = dm.asvs.final[a[j+1],a[j]]  ## why not with a diag 
      b = dm.asvs.diag[a[j+1],a[j]]  ## pairwise distances 
    } else {
      b= NA
    }
    vec = c(vec,b)
    
  }
  
  pure.cultures.stats = c(pure.cultures.stats, vec)
}

length(pure.cultures.stats) ## 478 

# the summary statistics of the pairwise-distances within the pure culture samples 

summary(pure.cultures.stats) 

sd(pure.cultures.stats) 

## Define the best-fitted distance to cluster ASVs into one adjusted ASV (referred to as OTU in what follows) 

# for each tested pairwise distance (media,mean, 3rd quantile, 3rd quantile+standard deviation, and the maximum), apply the following, distinctively: 

# Cluster the ASVs to a new adjusted ASV (here below referred to as OTU): ASVs which pairwise distance equals of is inferior than the tested distance (e.g. the mean), are assigned to the same new adjusted ASV, and their abundances are summed up 

# double check that asvs are indeed sorted by abundance 
colSums(asvs.table)


### Define the first adjusted ASV / OTU 

a_1 = which(dm.asvs.diag[,1] <= 0.0168)
OTU_1 = rownames(dm.asvs.diag)[a_1]
length(OTU_1) ## 320 
b_1 = setdiff(rownames(dm.asvs.diag),OTU_1)  # b are the ASVs which are NOT yet clustered into an OTU


otus.data = data.frame(nrow=nrow(dm.asvs.diag))
otus.asvs = OTU_1 ## vector list of all asvs already clustered into an OTU 
b_i = b_1 

for (i in 2:ncol(dm.asvs.diag)){
  
  OTU_i = vector()
  
  if (length(intersect(colnames(dm.asvs.diag)[i],b_i)) == 0){ # if the length of the intersection is 0, means that this ASV has already been clustered into an OTU
    
    OTU_i = NULL
    otus.asvs = otus.asvs
    otus.data = otus.data
    b_i = b_i    # b_i is the vector of the ASVs not yet assigned/clustered to an OTU 
    
  } else { 
    
    #a_i = which(dm.asvs.diag[,b_i[1]] <= 0.0168)
    a_i = which(dm.asvs.diag[,i] <= 0.0168) ## here assign the distance to be tested, here the rounded maximum pairwise-distance is used 
    OTU_i = rownames(dm.asvs.diag)[a_i]
    OTU_i= setdiff(OTU_i,otus.asvs) ### some ASVs could belong to several OTUs which lead/centroid ASVs are distant more than 0.0168 They are clustered to the most abundant ASV, when several matches occur
    otus.asvs = c(otus.asvs, OTU_i)
    OTU_i = c(OTU_i, rep("none", nrow(dm.asvs.diag) - length(OTU_i))) 
    otus.data = cbind(otus.data, OTU_i)
    b_i = setdiff(rownames(dm.asvs.diag), otus.asvs) # b_i is the vector of the ASVs not yet assigned/clustered to an OTU 
    
  } 
  #print(length(otus.asvs))
  # when length b_i = 0 all OTUs have been defined 
  if (length(b_i) ==0) {
    print(i)
  }
  
  if(sort(table(otus.asvs), decreasing = T)[1] != 1){
    print(i)
  }
  
}

length(otus.asvs) 
length(b_i) # if 0 all asvs have been assigned to an OTU cluster 

otus.data = otus.data[,-1]

otus.data = data.frame(cbind(c(OTU_1, rep("none", nrow(dm.asvs.diag) - length(OTU_1))),otus.data))
ncol(otus.data) ## there 2084 adjusted ASVs 
colnames(otus.data)= paste(rep("OTU_", ncol(otus.data)),1:ncol(otus.data), sep="")


# the resulting otus.data is a data.frame which columns are the defined OTUs (adjusted ASVs), and rows are the ASVs which belong to that given OTU (adjusted ASV)

# Step 3: check the accuracy of the adjusted ASVs clustering based on the ratio of shared and non-shared, and unique adjusted ASVs, sample-wise within each pure culture samples (n=4), and across its technical replicates (between 2 to 6 )

# make an OTU abundance table (adjusted ASVs abundance table)

otus.table = data.frame(nrow=nrow(asvs.table), ncol=ncol(otus.data))

for (i in 1:ncol(otus.data)){
  #vec = vector()
  a = which(otus.data[,i] != "none")
  a = otus.data[,i][a]
  
  if(length(a) > 1){
    vec = apply(asvs.table[,a],1,sum)
  } else {
    vec = asvs.table[,a]
  }
  otus.table = cbind(otus.table,vec)
}

dim(otus.table)
otus.table = otus.table[,-c(1,2)]

colnames(otus.table) = colnames(otus.data)
rownames(otus.table) = rownames(asvs.table)

apply(asvs.table,1,sum) == apply(otus.table,1,sum)
table(apply(asvs.table,1,sum) == apply(otus.table,1,sum))

# define the number of detected ASVs, and OTUs (adjusted ASVs) for the pure cultures, and their technical replicates, distinctively (n=21)
# number of detected ASVs 
length(pure.cultures)

vec_asvs = vector()

for (i in 1:length(pure.cultures)){
  a = which(rownames(metadata.all.final) == pure.cultures[i])
  b = length(which(asvs.table[a,] > 0))
  vec_asvs = c(vec_asvs,b)
}

summary(vec_asvs)


# number of detected OTUs (adjusted ASVs)
vec_otus = vector()

for (i in 1:length(pure.cultures)){
  a = which(rownames(metadata.all.final) == pure.cultures[i])
  b = length(which(otus.table[a,] > 0))
  vec_otus = c(vec_otus,b)
}

summary(vec_otus)



pure.cultures.asvs.otus = cbind(pure.cultures,vec_asvs, vec_otus)
write.csv(pure.cultures.asvs.otus, "pure.cultures.asvs.otus.csv")

# Step 4: compare the technical replicates of the pure culture samples 

# define the fraction of shared OTUs within a sample, across its technical replicates 

pure.cultures.replicated = names(which(table(metadata.all.final[pure.cultures,]$Sample_name) > 1)) ## 4 samples, replicated 2 to 6 times 

shared.otus = vector()

for(i in 1:length(pure.cultures.replicated)){
  vec = vector()
  a= rownames(metadata.all.final)[which(metadata.all.final$Sample_name==pure.cultures.replicated[i])]
  
  for (j in 1:length(a)){
    
    b = colnames(otus.table)[which(otus.table[a[j],] > 0)]
    vec = c(vec,b)
  }
  
  shared.otus = c(shared.otus, length(which(table(vec) > 1)))
  
}

shared.otus = data.frame(shared.otus)
rownames(shared.otus) =  pure.cultures.replicated
colnames(shared.otus) = "Shared_OTUs"


# define the fraction of non shared otus within a sample, across its technical replicates 

non.shared.otus = vector()

for(i in 1:length(pure.cultures.replicated)){
  vec = vector()
  a= rownames(metadata.all.final)[which(metadata.all.final$Sample_name==pure.cultures.replicated[i])] ## the technical replicate samples of pure culture i 
  
  for (j in 1:length(a)){
    
    b = colnames(otus.table)[which(otus.table[a[j],] > 0)] ## detected otus in technical replicate J 
    vec = c(vec,b) ## otus names detected in the technical replicates
  }
  
  non.shared.otus = c(non.shared.otus, length(which(table(vec) == 1)))
  
}

non.shared.otus = data.frame(non.shared.otus)
rownames(non.shared.otus) =  pure.cultures.replicated
colnames(non.shared.otus) = "Non_Shared_OTUs"

# check the taxonomy of the OTUs (adjusted ASVs): Do the ASVs assigned to the same OTU have the same Taxonomy? 

test = data.frame()

for (i in 1:ncol(otus.data)){
  vec= vector()
  
  a = which(otus.data[,i] != "none")
  a = otus.data[,i][a] ## a is the list of ASVs which belong to OTU i  # take the taxonomy of the lead ASV / most abundant one 
  
  for (j in 1:length(a)){
    b = taxonomy[a[j],]
    vec = rbind(vec,b)
  }
  
  #if (length(unique(vec$Genus)) > 1){
  if (length(unique(vec$Genus)) > 2){
    print(i)
    #print(unique(vec$Genus))
    
    test = rbind(test, c(unique(vec$Genus), rep("none", 5 - length(unique(vec$Genus)))))
    
  }
  
}

colnames(test) = c("Genus1","Genus2","Genus3","Genus4","Genus5")

### ! there are 15 OTUs which have mis-classified ASVs, occuring mostly in the mock samples. Double check again after all decontamination/curations steps. Possibly the mis-classified asvs are contaminants within a given sample, and will be excluded 


# define the taxonomy table of the OTUs (adjusted ASVs): To each OTU bin, assign the taxonomy of the lead (most abundant) ASV

taxonomy.otus = data.frame()

for (i in 1:ncol(otus.data)){
  
  vec= vector()
  
  a = which(otus.data[,i] != "none")
  a = otus.data[,i][a] ## a is the list of ASVs which belong to OTU i 
  
  b = taxonomy[a[1],]
  
  taxonomy.otus = rbind(taxonomy.otus, b)
  
}

rownames(taxonomy.otus) = colnames(otus.data)


#compare the technical replicates of all replicated samples: postive controls, and BALFs specimens 

all.replicates = names(which(table(metadata.all.final$Sample_name)> 1))  ## 16 samples 

#  define the fraction of shared otus within a sample, across its technical replicates 

shared.otus.all = vector()

for (i in 1:length(all.replicates)){
  vec = vector()
  a = which(metadata.all.final$Sample_name== all.replicates[i])
  b = rownames(metadata.all.final)[a]
  
  for (j in 1:length(b)){
    
    g = colnames(otus.table)[which(otus.table[b[j],] > 0)]
    vec = c(vec,g)
  }
  shared.otus.all = c(shared.otus.all, length(which(table(vec) > 1)))
}

shared.otus.all = data.frame(shared.otus.all)
rownames(shared.otus.all) =  all.replicates
colnames(shared.otus.all) = "Shared_OTUs"


##  define the fraction of non shared otus within a sample, across its technical replicates 

non.shared.otus.all = vector()

for(i in 1:length(all.replicates)){
  vec = vector()
  a= rownames(metadata.all.final)[which(metadata.all.final$Sample_name==all.replicates[i])] ## the technical replicate samples of pure culture i 
  
  for (j in 1:length(a)){
    
    b = colnames(otus.table)[which(otus.table[a[j],] > 0)] ## detected otus in technical replicate J 
    vec = c(vec,b) ## otus names detected in the technical replicates
  }
  
  non.shared.otus.all = c(non.shared.otus.all, length(which(table(vec) == 1)))
  
}

non.shared.otus.all = data.frame(non.shared.otus.all)
rownames(non.shared.otus.all) =  all.replicates
colnames(non.shared.otus.all) = "Non_Shared_OTUs"


# define the fraction of shared ASVs within a sample, across its technical replicates

shared.asvs.all = vector()

for (i in 1:length(all.replicates)){
  vec = vector()
  a = which(metadata.all.final$Sample_name== all.replicates[i])
  b = rownames(metadata.all.final)[a]
  
  for (j in 1:length(b)){
    
    g = colnames(asvs.table)[which(asvs.table[b[j],] > 0)]
    vec = c(vec,g)
  }
  shared.asvs.all = c(shared.asvs.all, length(which(table(vec) > 1)))
}

shared.asvs.all = data.frame(shared.asvs.all)
rownames(shared.asvs.all) =  all.replicates
colnames(shared.asvs.all) = "Shared_ASVs"


# define the fraction of non shared otus within a sample, across its technical replicates

non.shared.asvs.all = vector()

for(i in 1:length(all.replicates)){
  vec = vector()
  a= rownames(metadata.all.final)[which(metadata.all.final$Sample_name==all.replicates[i])] ## the technical replicate samples of pure culture i 
  
  for (j in 1:length(a)){
    
    b = colnames(asvs.table)[which(asvs.table[a[j],] > 0)] ## detected asvs in technical replicate J 
    vec = c(vec,b) ## asvs names detected in the technical replicates
  }
  
  non.shared.asvs.all = c(non.shared.asvs.all, length(which(table(vec) == 1)))
  
}

non.shared.asvs.all = data.frame(non.shared.asvs.all)
rownames(non.shared.asvs.all) =  all.replicates
colnames(non.shared.asvs.all) = "Non_Shared_ASVs"

# define the number of technical replicates for each sample 

replication.numbers = vector()

for(i in 1:length(all.replicates)){
  vec = vector()
  a= rownames(metadata.all.final)[which(metadata.all.final$Sample_name==all.replicates[i])] ## the technical replicate samples of pure culture i 
  
  b = length(a)
  
  replication.numbers = c(replication.numbers, b)
  
}

replication.numbers = data.frame(cbind(all.replicates, replication.numbers))


# define summary data.frame of : i) sample replicated, number of technical replicates, number of shared ASVs across technical replicates, number of non shared ASVs across technical replicates, number of shared OTUs across technical replicates,number of non shared OTUs across technical replicates,

data.summary.replicates = data.frame(cbind(replication.numbers,shared.asvs.all,non.shared.asvs.all, shared.otus.all, non.shared.otus.all))


# define the number of deteced ASVs, and OTUs (adjusted ASVs) for the mock samples 

# number of detected ASVs
vec.asvs = vector()

for (i in 1:nrow(asvs.table[c(129:136),])){
  
  a= length(which(asvs.table[c(129:136)[i],] > 0))
  
  vec.asvs = c(vec.asvs,a)
}


# number of detected OTUs (adjusted ASVs)

vec.otus = vector()


for (i in 1:nrow(otus.table[c(129:136),])){
  
  a= length(which(otus.table[c(129:136)[i],] > 0))
  
  vec.otus = c(vec.otus,a)
}


data.summary.mocks = data.frame(cbind(rownames(asvs.table[c(129:136),]), vec.asvs, vec.otus))
colnames(data.summary.mocks) = c("Mock_Sample", "detected_asvs", "detected_otus")


##NOTE : The accuracy of the adjusted ASVs clustering is repeated, and evaluated for each tested distance to cluster the ASVs into OTUs (adjusted ASVs) namely median, mean, 3rd quantile, 3rd quantile + sd, and max.
# After evaluating all tested clustering distances, maximum distance (0.0168) yielded the highest fraction of shared OTUs (adjusted ASVs) across the technical replicates, and the best estimation of the true number of expected species in the mock samples.
# Thus, the final clustering of ASVs into OTUs (adjusted ASVs) was performed using the maximum distance defined in step 2.


####Part three :  Detection of contaminant, and spurious ASVs 

# build a reference of true resident, and true contaminant ASVs for the mock samples 

## make a list of all ASVs detected across the mock samples (stock, and dilutions)



which(rownames(asvs.table)=="Ultra-clean-Mock-D0-DNA_S142") #129 
which(rownames(asvs.table)=="Ultra-clean-Mock-D0-DNA_S44") # 130 
which(rownames(asvs.table)=="Ultra-clean-Mock-D1-DNA_S45") # 131 
which(rownames(asvs.table)=="Ultra-clean-Mock-D2-DNA_S46") # 132 
which(rownames(asvs.table)=="Ultra-clean-Mock-D3-DNA_S47") # 133 
which(rownames(asvs.table)=="Ultra-clean-Mock-D4-DNA_S48") # 134 
which(rownames(asvs.table)=="Ultra-clean-Mock-D8-DNA_S49") # 135


## 
mock.uc.asvs = vector()

for(i in 1:nrow(asvs.table[c(129:135),])){
  
  a = names(which(asvs.table[c(129:135)[i],] > 0))
  
  mock.uc.asvs = c(mock.uc.asvs, a)
  
}

length(mock.uc.asvs) ## 1015
length(unique(mock.uc.asvs)) ## 934

mock.uc.asvs = unique(mock.uc.asvs) ## there are 934 detected unique ASVs sequences in the mock samples 

#5.2. make a fasta file of only the detected asvs in the mock samples  

asvs.names.data$ASV_Name
asvs.names.data$ASV_Seq


asvs.fasta.mocks <- data.frame()

for (i in 1:length(mock.uc.asvs)){
  
  g <- which(mock.uc.asvs[i] == asvs.names.data$ASV_Name)
  a <- paste(">", asvs.names.data[g,]$ASV_Name, sep="")
  b <- rbind(a,asvs.names.data[g,]$ASV_Seq)
  asvs.fasta.mocks <- rbind(asvs.fasta.mocks,b)
  
}
write.table(asvs.fasta.mocks, "asvs.fasta.mocks.fasta", sep="\t")


## in Geneious Prime program, do the following: 

#downolead the 16S rRNA sequences of the species present in the mock samples from the manufacturer Zymo Research
#align the mock reference sequences, with the all the ASVs sequences (n=934)  detected in the mocks DNA extracted with the ultra clean kit (n=7), D0 x2, D1, D2, D3, D4, D8: missing D5, D6, D7  
# using Clustaw algorithm with the following parameters Clustal Omega 1.2.2, Number of refinement iterations = 100, Cluster size for mBed Trees n= 100
# export the resulting distance matrix calculated based on the alignment trimed to equal lengths across all aligned sequences 

#n R session, load the distance matrix caculated in Geneious 
dis.matrix = read.csv("Mock_ASVs_database_Clustaw_final.csv", header=T, row.names=1) 
str(dis.matrix)


colnames(dis.matrix)[935:944] #  species originally present in the mock samples 
rownames(dis.matrix)[935:944] # species originally present in the mock samples 

# for each detected ASV in the mock sample, define a status "MEMBER" or "CONTAMINANT", based on 97% similarity threshold:if and ASV does not reach 97% similarity threshold to one (or more for some Mock strains), it is considered a contaminant

status = vector() 
MOCK_ASV = vector()

for (i in 1:ncol(dis.matrix[,1:934])){
  
  #print(colnames(dis.matrix)[i])
  
  b =  which(dis.matrix[c(935:944),i] >= 97.00) 
  
  
  if (length(b) >  0) {
    
    g = "MEMBER"
  } else { 
    
    g = "CONTAMINANT"
    
  }
  
  if (g == "MEMBER"){
    m = rownames(dis.matrix)[935:944][b[1]]
    
  } else {
    m = "NONE"
  }
  
  status = c(status, g)
  MOCK_ASV = c(MOCK_ASV,m)
  
}

table(status)
# CONTAMINANT      MEMBER 
#       185          749

# define data1 as a two.column data.frame: first column vector names of all detected ASVs in the mock samples, second column the defined status of the ASV (contaminant or member)
data1 = data.frame(colnames(dis.matrix[,1:934]), cbind(status, MOCK_ASV))
colnames(data1)[1] = "Mock_ASV"

# Extract the abundance of each contaminant ASV in the mock samples, among the detected contaminant ASVs within each sample separately
data1.contaminants = data1[which(data1$status=="CONTAMINANT"),] ## only contaminants 

#define data2 as the relative abundance of each contaminant ASV across the mock samples 

data2 = data.frame(ncol=nrow(data1.contaminants), nrow=length(129:135))  ## 129:135 are the mock samples in asvs.table 

for (i in 1: nrow(data1.contaminants))  { ## contaminant ASVs 
  vec = vector()
  a = which(data1.contaminants[i,]$Mock_ASV==colnames(asvs.table))
  
  ## double check it is the right ASV 
  for (j in 129:135) { # mock samples in asvs.table
    
    b = asvs.table[j,a] ## abundance of the same contaminant ASV as i, across the different mock samples as j 
    vec = c(vec, b)
  }
  
  data2 = cbind(data2, vec)
}

data2 <- data2[,-c(1,2)]
rownames(data2) <- rownames(asvs.table)[129:135]
colnames(data2) <- data1.contaminants$Mock_ASV

data2.relative = data2

for (i in 1:nrow(data2)){
  data2.relative[i,] <- data2.relative[i,]/ rowSums(asvs.table[c(129:135),])[i]
}

##  define data3 data.frame of abundances of the true resident (member) ASVs across mock samples 

data1.residents = data1[which(data1$status=="MEMBER"),] ## only true residents 


data3 = data.frame(ncol=nrow(data1.residents), nrow=length(129:135))

for (i in 1: nrow(data1.residents))  { ## true residents 
  vec = vector()
  a = which(data1.residents[i,]$Mock_ASV==colnames(asvs.table))
  
  ## double check it is the right ASV 
  for (j in 129:135) { # mock samples in asvs.table
    
    b = asvs.table[j,a] ## abundance of the same true resident ASV as i, accross the different mock samples as j 
    vec = c(vec, b)
  }
  
  data3 = cbind(data3, vec)
}

data3 <- data3[,-c(1,2)]
rownames(data3) <- rownames(asvs.table)[129:135]
colnames(data3) <- data1.residents$Mock_ASV

data3.relative = data3

for (i in 1:nrow(data3)){
  data3.relative[i,] <- data3.relative[i,]/ rowSums(asvs.table[c(129:135),])[i]
}


# inspect the summary statistics of contaminant ASVs across mock samples 

data2.summary <- data.frame() 

for (i in 1:nrow(data2.relative)) {
  vec = vector()
  
  detected_i = length(which(data2.relative[i,] > 0))
  
  d_i= apply(data2.relative[i, which(data2.relative[i,] > 0)],1,summary)
  sd_i= apply(data2.relative[i, which(data2.relative[i,] > 0)],1,sd)
  vec = c(detected_i, t(d_i), t(sd_i))
  
  d_i= apply(data2.relative[i, which(data2.relative[i,] > 0)],2,summary)
  sd_i= apply(data2.relative[i, which(data2.relative[i,] > 0)],2,sd)
  vec = c(detected_i, t(d_i), t(sd_i))
  
  data2.summary <- rbind(data2.summary, vec)
}

colnames(data2.summary) <- c("detetced_ASV", "Min", "1st_Qu", "Median", "Mean", "3rd_Qu", "Max", "SD")
rownames(data2.summary) <- rownames(asvs.table)[129:135]

write.csv(data2.summary, "contaminants.mocks.summury.relative.abundances.csv")



# inspect the summary statistis of true resident (member) ASVs across mock samples 

#data3.summary <- data.frame(ncol=5, nrow=length(129:135))

data3.summary <- data.frame() 

for (i in 1:nrow(data3.relative)) {
  vec = vector()
  
  detected_i = length(which(data3.relative[i,] > 0))
  d_i= apply(data3.relative[i, which(data3.relative[i,] > 0)],1,summary)
  sd_i= apply(data3.relative[i, which(data3.relative[i,] > 0)],1,sd)
  vec = c(detected_i, t(d_i), t(sd_i))
  
  data3.summary <- rbind(data3.summary, vec)
}

colnames(data3.summary) <- c("detetced_ASV", "Min", "1st_Qu", "Median", "Mean", "3rd_Qu", "Max", "SD")
rownames(data3.summary) <- rownames(asvs.table)[129:135]

write.csv(data3.summary, "true.residents.mocks.summury.relative.abundances.csv")

# define the percentage abundances of contaminants across the mock samples 

vec.per.conta <- vector()

for (i in 1: nrow(asvs.table[129:135,])){
  a = data1[which(data1$status=="CONTAMINANT"),]$Mock_ASV
  g = c(129:135)[i]
  print(g)
  b = sum(asvs.table[g,a]) * 100 / sum(asvs.table[g,])
  vec.per.conta <- c(vec.per.conta,b)
}

contamination.percentages = data.frame(vec.per.conta)
rownames(contamination.percentages) <- rownames(asvs.table)[129:135]
colnames(contamination.percentages) <- "Percentage_Contam"

write.csv(contamination.percentages, "contamination.percentages.mocks.refereces.csv")



# Identification of contaminant OTUs (adjsuted ASVs) using "decontam"
library(phyloseq)
library(decontam)

ps.all.otus <- phyloseq(otu_table(as.matrix(otus.table), taxa_are_rows=FALSE), sample_data(metadata.all.final),  tax_table(as.matrix(taxonomy.otus)))

## identify the contaminant OTUs (adjusted ASVs) using the frequency method, and different threshold paramters 

otus.all.contamdf.freq.0.1 <- isContaminant(ps.all.otus, method="frequency", conc="FA_ng_microl_total", threshold=0.1)
head(otus.all.contamdf.freq.0.1) ## detected 54 contaminants 

otus.all.contamdf.freq.0.5 <- isContaminant(ps.all.otus, method="frequency", conc="FA_ng_microl_total", threshold=0.5)
head(otus.all.contamdf.freq.0.5) ## detected 312 contaminants 

## Assess the accuracy of each tested threshold parameter using the mock dilution samples, defining the fraction of false negatives, and false positive contaminants

# define mock.uc.otus : the vector all detected OTUs (adjusted ASVs) in the mock samples 


mock.uc.otus = vector()

for(i in 1:nrow(otus.table[c(129:135),])){
  
  a = which(otus.table[c(129:135)[i],] > 0)
  b = colnames(otus.table)[a]
  mock.uc.otus = c(mock.uc.otus,b)
}

length(mock.uc.otus)
length(unique(mock.uc.otus))## 124 

mock.uc.otus = unique(mock.uc.otus)


# define the list of identified contaminants for each tested threshold parameter 
otus.all.contaminants.0.1 = rownames(otus.all.contamdf.freq.0.1)[which(otus.all.contamdf.freq.0.1$contaminant==TRUE)]
otus.all.contaminants.0.5 = rownames(otus.all.contamdf.freq.0.5)[which(otus.all.contamdf.freq.0.5$contaminant==TRUE)]

##threshold = 0.1

intersect(mock.uc.otus,otus.all.contaminants.0.1) ## 4 detected in mocks, and defined as contaminants 
setdiff(mock.uc.otus,otus.all.contaminants.0.1)  ## 120 detected in mocks and defined as true residents
setdiff(otus.all.contaminants.0.1,mock.uc.otus)   ## 50 defined as contaminants and not detected in mocks 

# false negatives: defined as contaminant, but is a true resident (based on the reference built based on the Geneious alignement)

intersect(mock.uc.otus,otus.all.contaminants.0.1) ## 4 detected in mocks, and defined as contaminants 

false.negative.0.1 = vector()

for(i in 1:length(intersect(mock.uc.otus,otus.all.contaminants.0.1))){
  
  vec= vector()
  
  a = which(colnames(otus.data) == intersect(mock.uc.otus,otus.all.contaminants.0.1)[i])  ## col name in otus.data for contaminant otu i 
  b = otus.data[,a] ## the ASVs which belong to contaminant otu i 
  b = b[which(b != "none")]
  
  for (j in 1:length(b)){
    g = which(data1$Mock_ASV == b[j])
    x = data1[g,]$status
    vec = c(vec,x)
  }
  
  vec= unique(vec)
  
  false.negative.0.1 = c(false.negative.0.1, vec)  
}

# all identified contaminants are indeed contaminants ==> threshold 0.1  did not yield to false negative (defined a contaminant as a resident) 


# false positives 


setdiff(mock.uc.otus,otus.all.contaminants.0.1)  ## 120 detected in mocks and defined as true residents

false.positive.0.1 = vector()

for (i in 1:length(setdiff(mock.uc.otus,otus.all.contaminants.0.1))){
  vec= vector()
  a = which(colnames(otus.data) == setdiff(mock.uc.otus,otus.all.contaminants.0.1)[i])
  b = otus.data[,a]
  b = b[which(b != "none")]
  
  for (j in 1:length(b)){
    g = which(data1$Mock_ASV == b[j])
    x = data1[g,]$status
    vec = c(vec,x)
  }
  
  vec= unique(vec)
  
  if (length(vec) != 1){
    print(i)
    vec = paste(vec[1], vec[2], sep="_")
    
  }
  
  false.positive.0.1 = c(false.positive.0.1, vec)
  
}

length(false.positive.0.1) == length(setdiff(mock.uc.otus,otus.all.contaminants.0.1))

false.positive.0.1 = data.frame(false.positive.0.1)
rownames(false.positive.0.1) = setdiff(mock.uc.otus,otus.all.contaminants.0.1) 

### 91 contaminants are still NOT identified as as contaminant ==> threshold 0.1 missed 91 contaminants 

### threshold 0.5 

intersect(mock.uc.otus,otus.all.contaminants.0.5) ##  14 detected in mocks, and defined as contaminants 

false.negative.0.5 = vector()


for(i in 1:length(intersect(mock.uc.otus,otus.all.contaminants.0.5))){
  
  vec= vector()
  
  a = which(colnames(otus.data) == intersect(mock.uc.otus,otus.all.contaminants.0.5)[i])
  b = otus.data[,a]
  b = b[which(b != "none")]
  
  for (j in 1:length(b)){
    g = which(data1$Mock_ASV == b[j])
    x = data1[g,]$status
    vec = c(vec,x)
  }
  
  vec= unique(vec)
  
  false.negative.0.5 = c(false.negative.0.5, vec)  
}

## threshold identified OTU_35 as contaminant, while it is a true resident ==> 1 false negative 

# false positive 

setdiff(mock.uc.otus,otus.all.contaminants.0.5)  ##  110 detected in mocks and defined as true residents

false.positive.0.5 = vector()

for (i in 1:length(setdiff(mock.uc.otus,otus.all.contaminants.0.5))){
  vec= vector()
  a = which(colnames(otus.data) == setdiff(mock.uc.otus,otus.all.contaminants.0.5)[i])
  b = otus.data[,a]
  b = b[which(b != "none")]
  
  for (j in 1:length(b)){
    g = which(data1$Mock_ASV == b[j])
    x = data1[g,]$status
    vec = c(vec,x)
  }
  
  vec= unique(vec)
  
  if (length(vec) != 1){
    print(i)
    vec = paste(vec[1], vec[2], sep="_")
    
  }
  
  false.positive.0.5 = c(false.positive.0.5, vec)
  
}

length(false.positive.0.5) == length(setdiff(mock.uc.otus,otus.all.contaminants.0.5))

false.positive.0.5 = data.frame(false.positive.0.5)
rownames(false.positive.0.5) = setdiff(mock.uc.otus,otus.all.contaminants.0.5) 

# there are 82 contaminants which are still unidentified contaminants ==> threshold 0.5 missed 82 true contaminants 

# Test further threshold between 0.1 and 0.5 

otus.all.contamdf.freq.0.2 <- isContaminant(ps.all.otus, method="frequency", conc="FA_ng_microl_total", threshold=0.2)
head(otus.all.contamdf.freq.0.2)## detcted 98 contaminants 


otus.all.contamdf.freq.0.3 <- isContaminant(ps.all.otus, method="frequency", conc="FA_ng_microl_total", threshold=0.3)
head(otus.all.contamdf.freq.0.3)## detcted  144 contaminants 

# Evaluate these added thresholds as described above 

otus.all.contaminants.0.2 = rownames(otus.all.contamdf.freq.0.2)[which(otus.all.contamdf.freq.0.2$contaminant==TRUE)]
otus.all.contaminants.0.3 = rownames(otus.all.contamdf.freq.0.3)[which(otus.all.contamdf.freq.0.3$contaminant==TRUE)]

## threshold 0.2 

intersect(mock.uc.otus,otus.all.contaminants.0.2)##  5 detected in mocks, and defined as contaminants 
setdiff(mock.uc.otus,otus.all.contaminants.0.2) ##   119 detected in mocks and defined as true residents: decontaminated mock otus 
setdiff(otus.all.contaminants.0.2,mock.uc.otus) ##   defined as contaminants and not detected in mocks

# false negatives 

false.negative.0.2 = vector()


for(i in 1:length(intersect(mock.uc.otus,otus.all.contaminants.0.2))){
  
  vec= vector()
  
  a = which(colnames(otus.data) == intersect(mock.uc.otus,otus.all.contaminants.0.2)[i])
  b = otus.data[,a]
  b = b[which(b != "none")]
  
  for (j in 1:length(b)){
    g = which(data1$Mock_ASV == b[j])
    x = data1[g,]$status
    vec = c(vec,x)
  }
  
  vec= unique(vec)
  
  false.negative.0.2 = c(false.negative.0.2, vec)  
}

# all defined contaminants are indeed true contaminants 

# false positive 

setdiff(mock.uc.otus,otus.all.contaminants.0.2)  ##   detected in mocks and defined as true residents

false.positive.0.2 = vector()

for (i in 1:length(setdiff(mock.uc.otus,otus.all.contaminants.0.2))){
  vec= vector()
  a = which(colnames(otus.data) == setdiff(mock.uc.otus,otus.all.contaminants.0.2)[i])
  b = otus.data[,a]
  b = b[which(b != "none")]
  
  for (j in 1:length(b)){
    g = which(data1$Mock_ASV == b[j])
    x = data1[g,]$status
    vec = c(vec,x)
  }
  
  vec= unique(vec)
  
  if (length(vec) != 1){
    print(i)
    vec = paste(vec[1], vec[2], sep="_")
    
  }
  
  false.positive.0.2 = c(false.positive.0.2, vec)
  
}

length(false.positive.0.2) == length(setdiff(mock.uc.otus,otus.all.contaminants.0.2))

false.positive.0.2 = data.frame(false.positive.0.2)
rownames(false.positive.0.2) = setdiff(mock.uc.otus,otus.all.contaminants.0.2) 

# 90 contaminants are still undetected ==> threshold 0.2 missed 90 contaminants 


# threshold 0.3 

intersect(mock.uc.otus,otus.all.contaminants.0.3)##  6 detected in mocks, and defined as contaminants 
setdiff(mock.uc.otus,otus.all.contaminants.0.3) ##   118 detected in mocks and defined as true residents: decontaminated mock  OTUs
setdiff(otus.all.contaminants.0.3,mock.uc.otus) ##   defined as contaminants and not detected in mocks 

# false negatives 

false.negative.0.3 = vector()


for(i in 1:length(intersect(mock.uc.otus,otus.all.contaminants.0.3))){
  
  vec= vector()
  
  a = which(colnames(otus.data) == intersect(mock.uc.otus,otus.all.contaminants.0.3)[i])
  b = otus.data[,a]
  b = b[which(b != "none")]
  
  for (j in 1:length(b)){
    g = which(data1$Mock_ASV == b[j])
    x = data1[g,]$status
    vec = c(vec,x)
  }
  
  vec= unique(vec)
  
  false.negative.0.3 = c(false.negative.0.3, vec)  
}

# all defined contaminants are indeed contaminants 

# false positives 


false.positive.0.3 = vector()

for (i in 1:length(setdiff(mock.uc.otus,otus.all.contaminants.0.3))){
  vec= vector()
  a = which(colnames(otus.data) == setdiff(mock.uc.otus,otus.all.contaminants.0.3)[i])
  b = otus.data[,a]
  b = b[which(b != "none")]
  
  for (j in 1:length(b)){
    g = which(data1$Mock_ASV == b[j])
    x = data1[g,]$status
    vec = c(vec,x)
  }
  
  vec= unique(vec)
  
  if (length(vec) != 1){
    print(i)
    vec = paste(vec[1], vec[2], sep="_")
    
  }
  
  false.positive.0.3 = c(false.positive.0.3, vec)
  
}

length(false.positive.0.3) == length(setdiff(mock.uc.otus,otus.all.contaminants.0.3))
ot
false.positive.0.3 = data.frame(false.positive.0.3)
rownames(false.positive.0.3) = setdiff(mock.uc.otus,otus.all.contaminants.0.3) 

## still 89 contaminants not detected ==> threshold 0.3 missed 89 contaminants 


#### test threshold  0.4

otus.all.contamdf.freq.0.4 <- isContaminant(ps.all.otus, method="frequency", conc="FA_ng_microl_total", threshold=0.4)
head(otus.all.contamdf.freq.0.4) ## detected 206 contaminants 

otus.all.contaminants.0.4 = rownames(otus.all.contamdf.freq.0.4)[which(otus.all.contamdf.freq.0.4$contaminant==TRUE)]

# evaluate threshold 0.4

intersect(mock.uc.otus,otus.all.contaminants.0.4) ##  10  detected in mocks, and defined as contaminants 
setdiff(mock.uc.otus,otus.all.contaminants.0.4) ##    detected in mocks and defined as true residents: decontaminated mock otus 
setdiff(otus.all.contaminants.0.4,mock.uc.otus) ##   defined as contaminants and not detected in mocks 

# false negatives 

false.negative.0.4 = vector()


for(i in 1:length(intersect(mock.uc.otus,otus.all.contaminants.0.4))){
  w
  vec= vector()
  
  a = which(colnames(otus.data) == intersect(mock.uc.otus,otus.all.contaminants.0.4)[i])
  b = otus.data[,a]
  b = b[which(b != "none")]
  
  for (j in 1:length(b)){
    g = which(data1$Mock_ASV == b[j])
    x = data1[g,]$status
    vec = c(vec,x)
  }
  
  vec= unique(vec)
  
  false.negative.0.4 = c(false.negative.0.4, vec)  
}

# all defined contaminants are indeed contaminants, no false negatives: true resident detected as contaminant  

### false positives 


false.positive.0.4 = vector()

for (i in 1:length(setdiff(mock.uc.otus,otus.all.contaminants.0.4))){
  vec= vector()
  a = which(colnames(otus.data) == setdiff(mock.uc.otus,otus.all.contaminants.0.4)[i])
  b = otus.data[,a]
  b = b[which(b != "none")]
  
  for (j in 1:length(b)){
    g = which(data1$Mock_ASV == b[j])
    x = data1[g,]$status
    vec = c(vec,x)
  }
  
  vec= unique(vec)
  
  if (length(vec) != 1){
    print(i)
    vec = paste(vec[1], vec[2], sep="_")
    
  }
  
  false.positive.0.4 = c(false.positive.0.4, vec)
  
}

length(false.positive.0.4) == length(setdiff(mock.uc.otus,otus.all.contaminants.0.4))

false.positive.0.4 = data.frame(false.positive.0.4)
rownames(false.positive.0.4) = setdiff(mock.uc.otus,otus.all.contaminants.0.4) 

## 85 contaminants are still undetected ==> threshold 


#### Conclusion : threshold 0.5 detected the most contaminants, has the least number of false positive, but has one false negative. threshold 0.4 does not result in false negatives, and detected the highest number of true contaminants ==> select threshold 0.4 


### Define the best minimum abundance threshold of OTUs (adjusted ASVs) to exclude spurious ASVs sample-wise, based on the technical replicates (n=16): 4 pure cultures, 11 BALFs, and 1 mock D0 (stock solution) 

# explorative analyses 

# inspect the summary statistcis of the abundances of shared, and unique OTUs (adjusted ASVs) across each all replicates samples (n=16) 
shared.stats.replicates = vector()

for (i in 1:length(all.replicates)) {
  detected_j = vector()
  vec2 = vector()
  a = which(metadata.all.final$Sample_name== all.replicates[i])
  
  for (j in 1:length(a)){
    b = which(otus.table[rownames(metadata.all.final)[a][j],] > 0)
    detected_j= c(detected_j,colnames(otus.table)[b])
  }
  
  shared_i =  names(which(table(detected_j) == length(a)))## shared otus accross all technical replicates 
  #unique_i = names(which(table(detected_j) == 1))
  
  for (f in 1 :length(a)){
    
    d = otus.table[rownames(metadata.all.final)[a][f],shared_i] / sum(otus.table[rownames(metadata.all.final)[a][f],])
    vec2 = c(vec2,t(d)[,1])
  }
  
  shared.stats.replicates = c(shared.stats.replicates, vec2) 
}

summary(shared.stats.replicates)

sd(shared.stats.replicates)



# inspect the summary statistics of the abundances of shared, and unique OTUs (adjusted ASVs) only across the pure culture samples

shared.stats.replicates = vector()

for (i in 1:length(all.replicates[12:15])) {
  detected_j = vector()
  vec2 = vector()
  a = which(metadata.all.final$Sample_name== all.replicates[12:15][i])
  
  for (j in 1:length(a)){
    b = which(otus.table[rownames(metadata.all.final)[a][j],] > 0)
    detected_j= c(detected_j,colnames(otus.table)[b])
  }
  
  shared_i =  names(which(table(detected_j) == length(a)))
  #unique_i = names(which(table(detected_j) == 1))
  
  for (f in 1 :length(a)){
    d = otus.table[rownames(metadata.all.final)[a][f],shared_i] / sum(otus.table[rownames(metadata.all.final)[a][f],])
    vec2 = c(vec2,t(d)[,1] )
  }
  
  shared.stats.replicates = c(shared.stats.replicates,vec2) 
}

summary(shared.stats.replicates)


sd(shared.stats.replicates)

# inspect the summary statistics of the abundances of shared, and unique OTUs (adjusted ASVs) only across the BALFs specimens 

shared.stats.replicates = vector()

for (i in 1:length(all.replicates[1:11])) {
  detected_j = vector()
  vec2 = vector()
  a = which(metadata.all.final$Sample_name== all.replicates[1:11][i])
  
  for (j in 1:length(a)){
    b = which(otus.table[rownames(metadata.all.final)[a][j],] > 0)
    detected_j= c(detected_j,colnames(otus.table)[b])
  }
  
  shared_i =  names(which(table(detected_j) == length(a)))
  #unique_i = names(which(table(detected_j) == 1))
  
  for (f in 1 :length(a)){
    d = otus.table[rownames(metadata.all.final)[a][f],shared_i] / sum(otus.table[rownames(metadata.all.final)[a][f],])
    vec2 = c(vec2, t(d)[,1])
  }
  
  shared.stats.replicates = c(shared.stats.replicates, vec2) 
}

summary(shared.stats.replicates)

sd(shared.stats.replicates)



# inspect the summary statistics of the abundances of shared, and unique OTUs (adjusted ASVs) across BALFs specimens and mock samples 

shared.stats.replicates = vector()

for (i in 1:length(all.replicates[c(1:11,16)])) {
  detected_j = vector()
  vec2 = vector()
  a = which(metadata.all.final$Sample_name == all.replicates[c(1:11,16)][i])
  
  for (j in 1:length(a)){
    b = which(otus.table[rownames(metadata.all.final)[a][j],] > 0)
    detected_j= c(detected_j,colnames(otus.table)[b])
  }
  
  shared_i =  names(which(table(detected_j) == length(a)))
  #unique_i = names(which(table(detected_j) == 1))
  
  for (f in 1 :length(a)){
    d = otus.table[rownames(metadata.all.final)[a][f],shared_i] / sum(otus.table[rownames(metadata.all.final)[a][f],])
    vec2 = c(vec2, t(d)[,1])
  }
  
  shared.stats.replicates = c(shared.stats.replicates, vec2) 
}

summary(shared.stats.replicates)

sd(shared.stats.replicates) * 100



## Test the thresholds defined based on the technical replicates of BALFs, and mock specimens namely: median, mean, 3rd quantile, in addition to the usually recommended thresholds of 1%, and 2% : Sample-wise, if the abundance of a given OTU (adjusted ASV) is below the tested abundance threshold, its abundance is set to 0 within that sample
## and evaluate each tested abundance threshold based on its efficiency to detect contaminants in the mock samples 

# the mean 

round(0.0612843,digits=4) ## 0.0613

data.cuarated = data.frame()

for (i in 1 : nrow(otus.table[129:135,])){
  d2 = vector()
  d1 = otus.table[c(129:135)[i],] / sum(otus.table[c(129:135)[i],])
  
  
  for (j in 1:ncol(d1)){
    
    if(d1[,j] <= 0.0613){
      d1[,j] = 0 } else { 
        
        d1[,j] = d1[,j]
      }
    
    d2 = c(d2,d1[,j])
    
  }
  data.cuarated = rbind(data.cuarated, d2)
  
}

rownames(data.cuarated) = rownames(otus.table)[129:135]
colnames(data.cuarated) = colnames(otus.table)

# evaluation of the efficiency in detecting contaminants based on the reference built in Geneious 

rownames(data.cuarated) == rownames(otus.table)[129:135]

status_0 = vector()
status_1= vector()

for (i in 1:length(c(129:135))){
  
  vec1 = vector()
  vec2 =vector()
  vec3 = vector()
  vec4 = vector()
  
  a = which(otus.table[c(129:135)[i],] > 0)
  a_i = colnames(otus.table[c(129:135)[i],])[a]  ## names of OTUs detected in sample i before any curation 
  b = which(data.cuarated[i,] > 0) 
  b_i = colnames(data.cuarated[i,])[b]  ### names of OTUs detected in sample i after curation based on the selected abundance threshold 
  
  for(j in 1 : length(a_i)){
    x_i= which(colnames(otus.data)==a_i[j])
    x_i = otus.data[,x_i]
    x_i= x_i[which(x_i != "none")]
    x2_i = intersect(x_i, data1$Mock_ASV)
    vec1= c(vec1, x2_i)  ## vector of the detected asvs in sample i
  }
  
  vec1 = unique(vec1) ## asvs of the OTUs detected in sample i before any curation 
  
  for(k in 1:length(b_i)){
    y_i= which(colnames(otus.data)==b_i[k])
    y_i = otus.data[,y_i]
    y_i= y_i[which(y_i != "none")]
    y2_i = intersect(y_i, data1$Mock_ASV)
    vec2 = c(vec2, y2_i) ## vector of the detected asvs in sample i
  }
  
  vec2 = unique(vec2) ## asvs of the OTUs detected in sample i after curation based on the mean 
  
  for (f in 1:length(vec1)){
    #which(data1$Mock_ASV== vec1[f])
    g = data1[which(data1$Mock_ASV== vec1[f]),]$status
    vec3 = c(vec3,g) ## the status contaminant/resident of the asvs of the OTUs detected in sample i before any curation 
  }
  
  for(l in 1:length(vec2)){
    k = data1[which(data1$Mock_ASV== vec2[l]),]$status
    vec4 = c(vec4,k)  ### the status contaminant/resident of the asvs of the OTUs detected in sample i after curation based on the mean 
  } 
  
  status_0 = c(status_0, vec3) # across all mock samples 
  status_1 = c(status_1,vec4) # across all mock samples 
  
}

## 
table(status_0) 
# contaminant   member 
# 263             4655 

table(status_1)
# contaminant   member
#     0           3259

# the median 
round(0.0066327, digits=4)## 0.0066


data.cuarated = data.frame()

for (i in 1 : nrow(otus.table[129:135,])){
  d2 = vector()
  d1 = otus.table[c(129:135)[i],] / sum(otus.table[c(129:135)[i],])
  
  
  
  for (j in 1:ncol(d1)){
    
    if(d1[,j] <= 0.0066){
      d1[,j] = 0 } else { 
        
        d1[,j] = d1[,j]
      }
    
    d2 = c(d2,d1[,j])
    
  }
  data.cuarated = rbind(data.cuarated, d2)
  
}

rownames(data.cuarated) = rownames(otus.table)[129:135]
colnames(data.cuarated) = colnames(otus.table)

#Evaluation

rownames(data.cuarated) == rownames(otus.table)[129:135]

status_0 = vector() ## based on Geneious 
status_1=vector()   ## based on this curation 

for (i in 1:length(c(129:135))){
  vec1 = vector()
  vec2=vector()
  vec3 = vector()
  vec4 = vector()
  
  a = which(otus.table[c(129:135)[i],] > 0)
  a_i = colnames(otus.table[c(129:135)[i],])[a]
  
  b = which(data.cuarated[i,] > 0) 
  b_i = colnames(data.cuarated[i,])[b]
  
  for(j in 1 : length(a_i)){
    x_i= which(colnames(otus.data)==a_i[j])
    x_i = otus.data[,x_i]
    x_i= x_i[which(x_i != "none")]
    x2_i = intersect(x_i, data1$Mock_ASV)
    vec1= c(vec1, x2_i)  ## vector of the detected asvs in sample i
  }
  
  vec1 = unique(vec1)
  
  for(k in 1:length(b_i)){
    y_i= which(colnames(otus.data)==b_i[k])
    y_i = otus.data[,y_i]
    y_i= y_i[which(y_i != "none")]
    y2_i = intersect(y_i, data1$Mock_ASV)
    vec2 = c(vec2, y2_i) ## vector of the detected asvs in sample i
  }
  
  vec2 = unique(vec2)
  
  for (f in 1:length(vec1)){
    #which(data1$Mock_ASV== vec1[f])
    g = data1[which(data1$Mock_ASV== vec1[f]),]$status
    vec3 = c(vec3,g)
  }
  
  for(l in 1:length(vec2)){
    k = data1[which(data1$Mock_ASV== vec2[l]),]$status
    vec4 = c(vec4,k)
  }
  
  status_0 = c(status_0, vec3)
  status_1 = c(status_1,vec4)
  
}


table(status_0) 
# contaminant member
#   263        4655 

table(status_1)
# contaminant member 
#     14       4521


#the 3rd quantile
round(0.0190909, digits=4) ## 0.0191


data.cuarated = data.frame()

for (i in 1 : nrow(otus.table[129:135,])){
  d2 = vector()
  d1 = otus.table[c(129:135)[i],] / sum(otus.table[c(129:135)[i],])
  
  
  
  for (j in 1:ncol(d1)){
    
    if(d1[,j] <= 0.0191){
      d1[,j] = 0 } else { 
        
        d1[,j] = d1[,j]
      }
    
    d2 = c(d2,d1[,j])
    
  }
  data.cuarated = rbind(data.cuarated, d2)
  
}

rownames(data.cuarated) = rownames(otus.table)[129:135]
colnames(data.cuarated) = colnames(otus.table)

# evaluation 

rownames(data.cuarated) == rownames(otus.table)[129:135]

status_0 = vector()
status_1=vector()

for (i in 1:length(c(129:135))){
  vec1 = vector()
  vec2=vector()
  vec3 = vector()
  vec4 = vector()
  
  a = which(otus.table[c(129:135)[i],] > 0)
  a_i = colnames(otus.table[c(129:135)[i],])[a]
  
  b = which(data.cuarated[i,] > 0) 
  b_i = colnames(data.cuarated[i,])[b]
  
  for(j in 1 : length(a_i)){
    x_i= which(colnames(otus.data)==a_i[j])
    x_i = otus.data[,x_i]
    x_i= x_i[which(x_i != "none")]
    x2_i = intersect(x_i, data1$Mock_ASV)
    vec1= c(vec1, x2_i)  ## vector of the detected asvs in sample i
  }
  
  vec1 = unique(vec1)
  
  for(k in 1:length(b_i)){
    y_i= which(colnames(otus.data)==b_i[k])
    y_i = otus.data[,y_i]
    y_i= y_i[which(y_i != "none")]
    y2_i = intersect(y_i, data1$Mock_ASV)
    vec2 = c(vec2, y2_i) ## vector of the detected asvs in sample i
  }
  
  vec2 = unique(vec2)
  
  for (f in 1:length(vec1)){
    #which(data1$Mock_ASV== vec1[f])
    g = data1[which(data1$Mock_ASV== vec1[f]),]$status
    vec3 = c(vec3,g)
  }
  
  for(l in 1:length(vec2)){
    k = data1[which(data1$Mock_ASV== vec2[l]),]$status
    vec4 = c(vec4,k)
  }
  
  status_0 = c(status_0, vec3)
  status_1 = c(status_1,vec4)
  
}


table(status_0)
# contaminant   member 
# 263            4655 


table(status_1)
#Contaminant             Member 
#    2                   4256 


## 2% threshold 
data.cuarated = data.frame()

for (i in 1 : nrow(otus.table[129:135,])){
  d2 = vector()
  d1 = otus.table[c(129:135)[i],] / sum(otus.table[c(129:135)[i],])
  
  
  
  for (j in 1:ncol(d1)){
    
    if(d1[,j] <= 0.0200){
      d1[,j] = 0 } else { 
        
        d1[,j] = d1[,j]
      }
    
    d2 = c(d2,d1[,j])
    
  }
  data.cuarated = rbind(data.cuarated, d2)
  
}

rownames(data.cuarated) = rownames(otus.table)[129:135]
colnames(data.cuarated) = colnames(otus.table)

#Evaluation 

rownames(data.cuarated) == rownames(otus.table)[129:135]

status_0 = vector()
status_1=vector()

for (i in 1:length(c(129:135))){
  vec1 = vector()
  vec2=vector()
  vec3 = vector()
  vec4 = vector()
  
  a = which(otus.table[c(129:135)[i],] > 0)
  a_i = colnames(otus.table[c(129:135)[i],])[a]
  
  b = which(data.cuarated[i,] > 0) 
  b_i = colnames(data.cuarated[i,])[b]
  
  for(j in 1 : length(a_i)){
    x_i= which(colnames(otus.data)==a_i[j])
    x_i = otus.data[,x_i]
    x_i= x_i[which(x_i != "none")]
    x2_i = intersect(x_i, data1$Mock_ASV)
    vec1= c(vec1, x2_i)  ## vector of the detected asvs in sample i
  }
  
  vec1 = unique(vec1)
  
  for(k in 1:length(b_i)){
    y_i= which(colnames(otus.data)==b_i[k])
    y_i = otus.data[,y_i]
    y_i= y_i[which(y_i != "none")]
    y2_i = intersect(y_i, data1$Mock_ASV)
    vec2 = c(vec2, y2_i) ## vector of the detected asvs in sample i
  }
  
  vec2 = unique(vec2)
  
  for (f in 1:length(vec1)){
    #which(data1$Mock_ASV== vec1[f])
    g = data1[which(data1$Mock_ASV== vec1[f]),]$status
    vec3 = c(vec3,g)
  }
  
  for(l in 1:length(vec2)){
    k = data1[which(data1$Mock_ASV== vec2[l]),]$status
    vec4 = c(vec4,k)
  }
  
  status_0 = c(status_0, vec3)
  status_1 = c(status_1,vec4)
  
}


table(status_0)
# contaminant    member 
#  236            4655 

table(status_1)
# contaminant member 
# 2               4371 

# 1% threshold 
data.cuarated = data.frame()

for (i in 1 : nrow(otus.table[129:135,])){
  d2 = vector()
  d1 = otus.table[c(129:135)[i],] / sum(otus.table[c(129:135)[i],])
  
  
  
  for (j in 1:ncol(d1)){
    
    if(d1[,j] <= 0.0100){
      d1[,j] = 0 } else { 
        
        d1[,j] = d1[,j]
      }
    
    d2 = c(d2,d1[,j])
    
  }
  data.cuarated = rbind(data.cuarated, d2)
  
}

rownames(data.cuarated) = rownames(otus.table)[129:135]
colnames(data.cuarated) = colnames(otus.table)

#evaluation 

rownames(data.cuarated) == rownames(otus.table)[129:135]

status_0 = vector()
status_1=vector()

for (i in 1:length(c(129:135))){
  vec1 = vector()
  vec2=vector()
  vec3 = vector()
  vec4 = vector()
  
  a = which(otus.table[c(129:135)[i],] > 0)
  a_i = colnames(otus.table[c(129:135)[i],])[a]
  
  b = which(data.cuarated[i,] > 0) 
  b_i = colnames(data.cuarated[i,])[b]
  
  for(j in 1 : length(a_i)){
    x_i= which(colnames(otus.data)==a_i[j])
    x_i = otus.data[,x_i]
    x_i= x_i[which(x_i != "none")]
    x2_i = intersect(x_i, data1$Mock_ASV)
    vec1= c(vec1, x2_i)  ## vector of the detected asvs in sample i
  }
  
  vec1 = unique(vec1)
  
  for(k in 1:length(b_i)){
    y_i= which(colnames(otus.data)==b_i[k])
    y_i = otus.data[,y_i]
    y_i= y_i[which(y_i != "none")]
    y2_i = intersect(y_i, data1$Mock_ASV)
    vec2 = c(vec2, y2_i) ## vector of the detected asvs in OTU k 
  }
  
  vec2 = unique(vec2)
  
  for (f in 1:length(vec1)){
    #which(data1$Mock_ASV== vec1[f])
    g = data1[which(data1$Mock_ASV== vec1[f]),]$status
    vec3 = c(vec3,g)
  }
  
  for(l in 1:length(vec2)){
    k = data1[which(data1$Mock_ASV== vec2[l]),]$status
    vec4 = c(vec4,k)
  }
  
  status_0 = c(status_0, vec3)
  status_1 = c(status_1,vec4)
  
}


table(status_0)
# contaminant 	member 
#   263         4655 

table(status_1)
# contaminant 	member 
# 3             4486 


## minimum abundance threshold of 2% showed the best ratio of removed true contaminants, and retained true residents ==> 2% threshold was applied sample-wise on the entire dataset as follows:

# set abundance to 0 for all OTUs which minimum abundance is less than 2%, SAMPLE-WISE 

otus.table.curated = otus.table

for(i in 1:nrow(otus.table.curated)){
  
  a = (sum(otus.table.curated[i,])*2)/100 ## identify the 2% threshold of sample i 
  
  for(j in 1:ncol(otus.table.curated[i,])){
    if(otus.table.curated[i,j] < a){
      otus.table.curated[i,j] = 0 } else {
        otus.table.curated[i,j] = otus.table.curated[i,j]
      }
  }
}

## summary sampling depth 
summary(apply(otus.table.curated, 1,sum))


## summary OTUs : detected OTUs 186/2084, no singeltons. No doubletons 


# make a vector of the spurious OTUs 

bad.otus = names(which(apply(otus.table.curated,2,sum)== 0))
taxonomy.otus[bad.otus,] ## 1898 

# make a vector of contaminant, and spurious OTUs 

### go back to the decontam output threshold 0.4 

length(otus.all.contaminants.0.4) # 206 

intersect(bad.otus,otus.all.contaminants.0.4) ## 159 OTUs defined as contaminants by decontam, and as spurious by the 2% abundance threshold 

setdiff(bad.otus,otus.all.contaminants.0.4) ## 1739 OTUs defined as spurious by the 2% threshold, and not as contaminants by decontam 

setdiff(otus.all.contaminants.0.4,bad.otus) ## 47 detected as contaminants by decontam, by not as spurious by the 2% threshold

# define  one vector of the union of contaminant, and spurious OTUs 

bad.otus.all = unique(c(bad.otus, otus.all.contaminants.0.4))

final.otus = setdiff(colnames(otus.table),bad.otus.all)

otus.table.curated[,final.otus]

summary(apply(otus.table.curated[,final.otus],1,sum)) ### sample-wise 

summary(apply(otus.table.curated[,final.otus],2,sum)) ### otu wise


### Evaluate the efficiency of the decontamination steps based on the mock samples 


otus.table.curated[129:135,final.otus] ## only mock samples 

length(mock.uc.otus) ## 124 detected OTUs in the mock samples 

intersect(final.otus,mock.uc.otus) ## 16 OTUs detected in the mocks, and defined as a true resident. This intersection should not contain any left contaminants 


check.intersection = vector()
for (i in 1:length(intersect(final.otus,mock.uc.otus))){
  vec= vector()
  
  x =  which(colnames(otus.data)==intersect(final.otus,mock.uc.otus)[i])
  a = which(otus.data[,x] != "none")
  a = otus.data[,x][a] ## a is the list of ASVs which belong to OTU_i 
  b = intersect(data1$Mock_ASV,a) ## b is the asvs of OTU_i which are detected in the 
  
  for (j in 1:length(b)){
    x2= which(data1$Mock_ASV == b[j])
    g = data1[x2,]$status
    vec = c(vec,g)
  }
  check.intersection = c(check.intersection,vec)
}

table(check.intersection) ## 872 Members, no contaminants in the final selected OTUs present in the mock samples 

setdiff(final.otus, mock.uc.otus) ## 123 defined as true resident, and not detected in the curated/decontaminated mock samples. No further checking required for these OTUs. 

setdiff(mock.uc.otus,final.otus) ## 108 OTUs detected in the mocks, and not retained/defined as true resident in the mock samples. These should not contain any resident. All should be contaminants 

check.diff = vector()
check.abund = vector()

for (i in 1:length(setdiff(mock.uc.otus,final.otus))){
  vec= vector()
  vec2 =vector()
  
  x =  which(colnames(otus.data)==setdiff(mock.uc.otus,final.otus)[i])
  a = which(otus.data[,x] != "none")
  a = otus.data[,x][a] ## a is the list of ASVs which belong to OTU_i 
  b = intersect(data1$Mock_ASV,a) ## b is the asvs of OTU_i which are detected in the 
  
  for (j in 1:length(b)){
    x2= which(data1$Mock_ASV == b[j])
    g = data1[x2,]$status
    
    if (g == "MEMBER"){
      print(i)
      #x3= data1[x2,]$Mock_ASV
    }
    vec = c(vec,g)
    vec2= c(vec2,b[j])
  }
  check.diff = c(check.diff,vec)
  check.abund = c(check.abund,vec2)
}

table(check.diff)
# check.diff
# CONTAMINANT      MEMBER
#185                 17

# These 17 true resident OTUs, were not retained in the final datasets ==> They are present accros the mock with abundances lower than 2%, thus considered spurious, and excluded from the datasets. 

# define the percentage of contamination in mock samples across the contamination 

metadata.all.final[129:135,]$FA_ng_microl_total

otus.final.check = vector()

for(i in 1:length(final.otus)){
  vec = vector()
  
  x =  which(colnames(otus.data)== final.otus[i])
  a = which(otus.data[,x] != "none")
  a = otus.data[,x][a] ## a is the list of ASVs which belong to OTU_i 
  b = intersect(data1$Mock_ASV,a) ## b is the asvs of OTU_i which are detected in the 
  
  for(j in 1:length(b)){
    x2= which(data1$Mock_ASV == b[j])
    g = data1[x2,]$status
    vec = c(vec,g)
  }
  otus.final.check = c(otus.final.check,vec)
  
}

table(otus.final.check)
# Member n= 732 
### the final otus only contain true resident 

summary(metadata.all.final[129:135,]$FA_ng_microl_total)


#  contamination.percentages across the dilutions 

#  Mock Sample                    Percentage_Contam
#Ultra-clean-Mock-D0-DNA_S142         3.0615400
#Ultra-clean-Mock-D0-DNA_S44          3.5204567
#Ultra-clean-Mock-D1-DNA_S45          4.8419598
#Ultra-clean-Mock-D2-DNA_S46          0.8346836
#Ultra-clean-Mock-D3-DNA_S47          5.0113505
#Ultra-clean-Mock-D4-DNA_S48          0.2855404
#Ultra-clean-Mock-D8-DNA_S49          7.7313770 

## Dilution D0 has the lowest concentration, and the highest contamination rate. After all the decontamination steps, only true resident are left in mock sample D8

## inspect the concentration of the BALFs samples, and how thier concentrations compare to the mock samples 
# exclude BALF samples which concentration is less than the concentration of mock sample D8 

summary(metadata.all.final[c(1:128,136),]$FA_ng_microl_total)
which(metadata.all.final[1:128,]$FA_ng_microl_total < 2.00) ## 10 


rownames(metadata.all.final[1:128,])[which(metadata.all.final[1:128,]$FA_ng_microl_total < 2.00)]

#[1] "1-BAL-S-6582_S259"               "BAL-bis-Rnaesy_S216"
#[3] "BAL-S-2075-DNA-Extraction2_S68"  "BAL-S-2077-DNA-Extraction2_S278"
#[5] "BAL-S-6026-DNA-Extraction2_S64"  "BAL-S-6034-DNA-Extraction2_S65"
#[7] "BAL-S-6328-DNA-Extraction2_S66"  "BAL-S-6421-DNA-Extraction2_S67"
#[9] "BAL-S-8000-DNA-Extraction2_S70"  "BAL-S-8383-DNA-Extraction2_S71"


## prepare the final data.frames 

# for OTUs datasset (adjusted ASVs), final otus n=139 

otus.table.final = otus.table.curated[,final.otus]

summary(apply(otus.table.final,1,sum)) ## sample wise 
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#0     11699   15998   16286   19684   43209


## remove too low concentration samples and samples with 0 reads after removal of contaminant and spurious OTUs

which(rownames(otus.table.final)=="1-BAL-S-6582_S259") ## 1 
which(rownames(otus.table.final)=="BAL-bis-Rnaesy_S216") ## 7
which(rownames(otus.table.final)=="BAL-S-2075-DNA-Extraction2_S68") ## 22 
which(rownames(otus.table.final)=="BAL-S-2077-DNA-Extraction2_S278") ## 23 
which(rownames(otus.table.final)=="BAL-S-6026-DNA-Extraction2_S64") ## 49 
which(rownames(otus.table.final)=="BAL-S-6034-DNA-Extraction2_S65") ## 53 
which(rownames(otus.table.final)=="BAL-S-6328-DNA-Extraction2_S66") ## 61 
which(rownames(otus.table.final)=="BAL-S-6421-DNA-Extraction2_S67") ## 62 
which(rownames(otus.table.final)=="BAL-S-8000-DNA-Extraction2_S70") ## 93
which(rownames(otus.table.final)=="BAL-S-8383-DNA-Extraction2_S71") ## 98 
which(rownames(otus.table.final)=="OR-Ultra-clean-DNA_S26") ##112 

otus.table.final= otus.table.final[-c(1,7,22,23,49,53,61,62,93,98,112),]

which(apply(otus.table.final,1,sum)==0)
summary(apply(otus.table.final,2,sum))
length(which(apply(otus.table.final,2,sum)==0))

vec = vector()
for(i in 1:length(names(which(apply(otus.table.final,2,sum)==0)))){
  a = which(colnames(otus.table.final) ==names(which(apply(otus.table.final,2,sum)==0))[i])
  vec = c(vec,a)
}

otus.table.final = data.frame(otus.table.final[,-c(39,53,57,62,65,66,74,101,106,113,116,121,123,131,138,139)])

dim(otus.table.final) ## 123 OTUs and 125 samples 


# for the ASVs

## remove too low concentration samples and samples with 0 reads after removal of contaminant and spurious ASVs

which(rownames(asvs.table)=="1-BAL-S-6582_S259") ##  
which(rownames(asvs.table)=="BAL-bis-Rnaesy_S216") ## 
which(rownames(asvs.table)=="BAL-S-2075-DNA-Extraction2_S68") ##  
which(rownames(asvs.table)=="BAL-S-2077-DNA-Extraction2_S278") ## 
which(rownames(asvs.table)=="BAL-S-6026-DNA-Extraction2_S64") ## 
which(rownames(asvs.table)=="BAL-S-6034-DNA-Extraction2_S65") ## 
which(rownames(asvs.table)=="BAL-S-6328-DNA-Extraction2_S66") ## 
which(rownames(asvs.table)=="BAL-S-6421-DNA-Extraction2_S67") ## 
which(rownames(asvs.table)=="BAL-S-8000-DNA-Extraction2_S70") ## 
which(rownames(asvs.table)=="BAL-S-8383-DNA-Extraction2_S71") ##  

which(rownames(asvs.table)=="OR-Ultra-clean-DNA_S26") ## 

asvs.table.final =  data.frame(asvs.table[-c(1,7,22,23,49,53,61,62,93,98,112),])

dim(asvs.table.final) ## 125 samples, and 12987 asvs 


### Select the final/clean ASVs according to the final OTUs (adjusted ASVs)

for (i in 1:nrow(asvs.table.final)){
  vec = vector()
  vec2 = vector()
  x = which(otus.table.final[i,] > 0) ## select only the otus which are detected in sample i 
  g = colnames(otus.table.final)[x] ## detected otus names 
  
  for(j in 1:length(g)){  ## for each otu, sum up the asvs separately 
    x2 = which(colnames(otus.data) == g[j])
    x3 = otus.data[,x2][which(otus.data[,x2] != "none")] ## asvs which belong to otu j 
    vec = c(vec, x3)
  }
  
  if (sum(otus.table.final[i,g]) != sum(asvs.table.final[i,vec])){
    print(i)
  }
  
  #length(vec)
  vec2 = setdiff(colnames(asvs.table.final[i,]),vec)
  
  ## all asvs which are not in vec, have their  abundances = 0
  for(f in 1:length(vec2)){
    h= which(colnames(asvs.table.final[i,]) == vec2[f])
    asvs.table.final[i,h] = 0
  }
  
  if (sum(otus.table.final[i,]) != sum(asvs.table.final[i,])){
    print(i)
  }
  
}

write.csv(otus.table.final, "otus.table.final.csv")
write.csv(asvs.table.final, "asvs.table.final.csv")

#Prepare the final metadata 

metadata.all.final2 = metadata.all.final

metadata.all.final2 = metadata.all.final2[rownames(asvs.table.final),]

write.csv(metadata.all.final2, "metadata.all.final2.csv")

#Prepare the final taxonomy files 

## ASVs taxonomy 
names(which(apply(asvs.table.final,2,sum) == 0))   ## 8874 ASVs
names(which(apply(asvs.table.final,2,sum) > 0))  


taxonomy.final = taxonomy[names(which(apply(asvs.table.final,2,sum) > 0)),]
dim(taxonomy.final) ## 4113 ASVs 

write.csv(taxonomy.final, "taxonomy.final.csv")


## OTUs (adjusted ASVs) taxonomy 
taxonomy.otus ## 2084 otus 
dim(otus.table.final) ## 123 otus 

taxonomy.otus.final = taxonomy.otus[colnames(otus.table.final),] ## 123 
write.csv(taxonomy.otus.final, "taxonomy.otus.final.csv")


##  Make genus level taxonomy 

# Sum up the abundances of the ASVs that belong to the same genus 

# Substitute "NA" with "Unclassified"

tax.final.genera <- taxonomy.otus.final
tax.final.genera[is.na(tax.final.genera)] <- "Unclassified"

# Assign a unique identifier to the unclassified ASVs at the genus levels 
unclassified.genera <- which(tax.final.genera$Genus=="Unclassified") ##

Genus.new <- vector()

for (i in 1:length(unclassified.genera)){
  
  a <- "Unclassified"
  
  if (tax.final.genera[unclassified.genera[i],5] != "Unclassified")  { #Family 
    a <- tax.final.genera[unclassified.genera[i],5]
  } 
  
  
  if (a == "Unclassified"){
    a <- tax.final.genera[unclassified.genera[i],4]#Order 
  }
  
  
  if (a == "Unclassified"){
    a <- tax.final.genera[unclassified.genera[i],3]# Class
  }
  
  
  if (a == "Unclassified"){ #phylum
    a <- tax.final.genera[unclassified.genera[i],2] #phylum
  }  # after this step the ASVs designed as "Unclassified" are unclassified even at the Phylum level 
  
  
  b <- paste("Unclassified", a, sep="_")## "Unclassified_Unclassfied" are unclassified at the phylum level, can be substitute by "Unclassified_Bacteria"
  
  Genus.new <- c(Genus.new, b)
}

# Substitute "Unclassified_Unclassified" with "Unclassified_Bacteria"
which(Genus.new=="Unclassified_Unclassified") ## none 
length(which(Genus.new=="Unclassified_Unclassified")) == length(which(tax.final.genera[unclassified.genera,2]=="Unclassified"))

Genus.new[which(Genus.new=="Unclassified_Unclassified")] ="Unclassified_Bacteria"

# Insert the new genera names into the taxonomy file 

for (i in 1:length(unclassified.genera)){
  tax.final.genera[unclassified.genera[i],6] <- Genus.new[i]
}

tax.final.genera$Genus


redundant.genera <- names(which(table(tax.final.genera$Genus) > 1)) ## 25 
unique.genera    <- names(which(table(tax.final.genera$Genus) == 1))## 19 

##New asvs.table.final.genera with the abundances of the ASVs that belong to the same genus summed within samples 

data.frame1 <- data.frame(nrow=nrow(otus.table.final), ncol=length(redundant.genera)) #rows are samples, column number now is length of redundant genera 

for (j in 1:length(redundant.genera)){
  vec <- vector()
  g <- which(tax.final.genera$Genus==redundant.genera[j])##select the rows that are classified the that redundant genus j 
  ##since colnames of otus.table,final and tax.final are equal 
  
  #print(rownames(tax.final.genera) [g] == colnames(otus.table.final) [g])
  vec <-apply(otus.table.final[,g],1,sum)
  data.frame1 <- cbind(data.frame1, vec)
} 

data.frame1 <- data.frame1[,-c(1,2)]
colnames(data.frame1) <- redundant.genera

data.frame2 <-data.frame(nrow=nrow(otus.table.final), ncol=length(unique.genera))

for (f in 1:length(unique.genera)){
  gg <- which(tax.final.genera$Genus==unique.genera[f])
  print(rownames(tax.final.genera)[gg] == colnames(otus.table.final)[gg])
  data.frame2 <- cbind(data.frame2,otus.table.final[,gg])
}


data.frame2 <- data.frame2[,-c(1,2)]
colnames(data.frame2) <- unique.genera
rownames(data.frame2) <- rownames(otus.table.final)


tax.final.genera.abundance  <- as.data.frame(cbind(data.frame1,data.frame2))
##In this data.frame all colnames MUST be unique 


str(tax.final.genera.abundance)
dim(tax.final.genera.abundance)

write.csv(tax.final.genera.abundance, "tax.final.genera.abundance.csv")

##  Make phylum level taxonomy 


# Substitute "NA" with "Unclassified"

tax.final.phyla <- taxonomy.otus.final

tax.final.phyla[is.na(tax.final.phyla)] <- "Unclassified"

# Group the ASVs which phylum is unclassified to "Unclassified_Bacteria"
length(unique((tax.final.phyla$Phylum)))##how many detected phyla 7 

which(tax.final.phyla$Phylum=="Unclassified")

tax.final.phyla[which(tax.final.phyla$Phylum=="Unclassified"),]$Phylum <- "Unclassified_Bacteria"

##Sum up the abundances of ASVs that belong to the same Phylum
tax.final.phyla.abundance <- data.frame(nrow=nrow(otus.table.final), ncol=length(unique(tax.final.phyla$Phylum)))

for (i in 1:length(unique(tax.final.phyla$Phylum))){
  
  ggg <- which(tax.final.phyla$Phylum ==unique(tax.final.phyla$Phylum)[i])
  #print(colnames(otus.table.final[,ggg])==rownames(tax.final.phyla)[ggg])
  
  
  if (length(ggg) > 1){
    hhh <- apply(otus.table.final[,ggg],1,sum)
  } else {
    hhh <- otus.table.final[,ggg]
  }
  
  tax.final.phyla.abundance <- cbind(tax.final.phyla.abundance,hhh)
}


tax.final.phyla.abundance <- tax.final.phyla.abundance[,-c(1,2)]
colnames(tax.final.phyla.abundance) <- unique(tax.final.phyla$Phylum)


write.csv(tax.final.phyla.abundance,"tax.final.phyla.abundance.csv")


## Normalization of the sequencing depth : sampling depth = 10 000 reads / sample 

## At the ASVs level 

## exclude samples which total read depth is lower that the sampling threshold 
which(rowSums(asvs.table.final)<10000) ## 19 samples, including 8 BALFs samples, 6 different BALFs samples ( are technical replicates)

asvs.table.final2 = asvs.table.final

which(rowSums(asvs.table.final2)<10000)

which(rowSums(asvs.table.final2[-c(1,2,3,5,12,16,18,19,31,32,33,44,55,95,103,106,110,114,115),])<10000)

asvs.table.final2 = asvs.table.final2[-c(1,2,3,5,12,16,18,19,31,32,33,44,55,95,103,106,110,114,115),]

rarefied.data.once <- data.frame()

for (f in 1:nrow(asvs.table.final2)){  # within every sample 
  
  vector1 <- vector()
  for (i in 1:ncol(asvs.table.final2)){ # for every ASV replicate n abundance times 
    a <- rep(colnames(asvs.table.final2)[i], asvs.table.final2[f,i])
    vector1 <- c(vector1, a)
  }
  print(length(vector1)==sum(asvs.table.final2[f,]))
  
  vector2 <- vector()
  vector2 <- sample(vector1, 10000, replace=F)## 
  
  vector3 <- vector() #put all the ASVs repeated n times together 
  
  for (j in 1:ncol(asvs.table.final2)){
    b <- length(which(vector2==colnames(asvs.table.final2)[j]))
    vector3 <- c(vector3,b)
  }
  rarefied.data.once <- rbind(rarefied.data.once, vector3)
}

colnames(rarefied.data.once) <- colnames(asvs.table.final2)
rownames(rarefied.data.once)<- rownames(asvs.table.final2)
dim(rarefied.data.once)


write.csv(rarefied.data.once, "asvs.table.final.normalized.csv")

### Accordignly adjust the OTUs, genera, and phylum datasets 

## OTUs

otus.table.final.rarefied = otus.table.final 

##remove samples which total read depth is lower that the sampling threshold 
which(rowSums(otus.table.final.rarefied) < 10000)

otus.table.final.rarefied = otus.table.final.rarefied[-c(1,2,3,5,12,16,18,19,31,32,33,44,55,95,103,106,110,114,115),]

for (i in 1:nrow(rarefied.data.once)){
  
  a = which(rarefied.data.once[i,] > 0)
  a = colnames(rarefied.data.once[i,])[a]
  
  for (j in 1:ncol(otus.data)){
    
    b = otus.data[,j][which(otus.data[,j] != "none")] ## b are all the asvs belong to otu j 
    x= intersect(b,a)## x is the asvs which belong to otus j which are detected in sample i 
    g = which(colnames(otus.table.final.rarefied) == colnames(otus.data)[j])
    
    if(length(x)==0){
      otus.table.final.rarefied[i,g] = 0
    } else {
      otus.table.final.rarefied[i,g] = sum(rarefied.data.once[i,x])
    }
  }
  
}

write.csv(otus.table.final.rarefied, "otus.table.final.rarefied.csv")


## Genera 

## Substitute "NA" with "Unclassified"

tax.final.genera <- taxonomy.otus.final
tax.final.genera[is.na(tax.final.genera)] <- "Unclassified"

##2. Assign a unique identifier to the unclassified ASVs at the genus levels 
unclassified.genera <- which(tax.final.genera$Genus=="Unclassified") ## 14 

Genus.new <- vector()

for (i in 1:length(unclassified.genera)){
  
  a <- "Unclassified"
  
  if (tax.final.genera[unclassified.genera[i],5] != "Unclassified")  { #Family 
    a <- tax.final.genera[unclassified.genera[i],5]
  } 
  
  
  if (a == "Unclassified"){
    a <- tax.final.genera[unclassified.genera[i],4]#Order 
  }
  
  
  if (a == "Unclassified"){
    a <- tax.final.genera[unclassified.genera[i],3]# Class
  }
  
  
  if (a == "Unclassified"){ #phylum
    a <- tax.final.genera[unclassified.genera[i],2] #phylum
  }  # after this step the ASVs designed as "Unclassified" are unclassified even at the Phylum level 
  
  
  b <- paste("Unclassified", a, sep="_")## "Unclassified_Unclassfied" are unclassified at the phylum level, can be substitute by "Unclassified_Bacteria"
  
  Genus.new <- c(Genus.new, b)
}


##Substitute "Unclassified_Unclassified" with "Unclassified_Bacteria"
which(Genus.new=="Unclassified_Unclassified") ## none 
length(which(Genus.new=="Unclassified_Unclassified")) == length(which(tax.final.genera[unclassified.genera,2]=="Unclassified"))

Genus.new[which(Genus.new=="Unclassified_Unclassified")] ="Unclassified_Bacteria"


for (i in 1:length(unclassified.genera)){
  tax.final.genera[unclassified.genera[i],6] <- Genus.new[i]
}

tax.final.genera$Genus


redundant.genera <- names(which(table(tax.final.genera$Genus) > 1)) ## 25 
unique.genera    <- names(which(table(tax.final.genera$Genus) == 1))##  19

##New asvs.table.final.genera with the abundances of the ASVs that belong to the same genus summed within samples 

data.frame1 <- data.frame(nrow=nrow(otus.table.final.rarefied), ncol=length(redundant.genera)) #rows are samples, column number now is length of redundant genera 

for (j in 1:length(redundant.genera)){
  vec <- vector()
  g <- which(tax.final.genera$Genus==redundant.genera[j])##select the rows that are classified the that redundant genus j 
  ##since colnames of otus.table,final and tax.final are equal 
  
  #print(rownames(tax.final.genera) [g] == colnames(otus.table.final.rarefied) [g])
  vec <-apply(otus.table.final.rarefied[,g],1,sum)
  data.frame1 <- cbind(data.frame1, vec)
} 

data.frame1 <- data.frame1[,-c(1,2)]
colnames(data.frame1) <- redundant.genera

data.frame2 <-data.frame(nrow=nrow(otus.table.final.rarefied), ncol=length(unique.genera))

for (f in 1:length(unique.genera)){
  gg <- which(tax.final.genera$Genus==unique.genera[f])
  print(rownames(tax.final.genera)[gg] == colnames(otus.table.final.rarefied)[gg])
  data.frame2 <- cbind(data.frame2,otus.table.final.rarefied[,gg])
}


data.frame2 <- data.frame2[,-c(1,2)]
colnames(data.frame2) <- unique.genera
rownames(data.frame2) <- rownames(otus.table.final.rarefied)


tax.final.genera.abundancel.rarefied  <- as.data.frame(cbind(data.frame1,data.frame2))
##In this data.frame all colnames MUST be unique 


str(tax.final.genera.abundancel.rarefied)
dim(tax.final.genera.abundancel.rarefied)

write.csv(tax.final.genera.abundancel.rarefied, "tax.final.genera.abundancel.rarefied.csv")

# Phyla 

#Substitute "NA" with "Unclassified"

tax.final.phyla <- taxonomy.otus.final

tax.final.phyla[is.na(tax.final.phyla)] <- "Unclassified"

##2. Group the ASVs which phylum is unclassified to "Unclassified_Bacteria"
length(unique((tax.final.phyla$Phylum)))##how many detected phyla 7 

which(tax.final.phyla$Phylum=="Unclassified") ## none 

tax.final.phyla[which(tax.final.phyla$Phylum=="Unclassified"),]$Phylum <- "Unclassified_Bacteria"

##Sum up the abundances of ASVs that belong to the same Phylum
tax.final.phyla.abundance.rarefied <- data.frame(nrow=nrow(otus.table.final.rarefied), ncol=length(unique(tax.final.phyla$Phylum)))

for (i in 1:length(unique(tax.final.phyla$Phylum))){
  
  ggg <- which(tax.final.phyla$Phylum ==unique(tax.final.phyla$Phylum)[i])
  #print(colnames(otus.table.final.rarefied[,ggg])==rownames(tax.final.phyla)[ggg])
  
  
  if (length(ggg) > 1){
    hhh <- apply(otus.table.final.rarefied[,ggg],1,sum)
  } else {
    hhh <- otus.table.final.rarefied[,ggg]
  }
  
  tax.final.phyla.abundance.rarefied <- cbind(tax.final.phyla.abundance.rarefied,hhh)
}


tax.final.phyla.abundance.rarefied <- tax.final.phyla.abundance.rarefied[,-c(1,2)]
colnames(tax.final.phyla.abundance.rarefied) <- unique(tax.final.phyla$Phylum)


write.csv(tax.final.phyla.abundance.rarefied,"tax.final.phyla.abundance.rarefied.csv")

## Metadata 

metadata.all.final.rarefied = metadata.all.final2
dim(rarefied.data.once) ## 106 samples 

intersect(rownames(rarefied.data.once), rownames(metadata.all.final.rarefied))### 106 all of the samples 

metadata.all.final.rarefied = metadata.all.final.rarefied[intersect(rownames(rarefied.data.once), rownames(metadata.all.final.rarefied)),]

write.csv(metadata.all.final.rarefied, "metadata.all.final.rarefied.csv")

## Explore further clustering thresholds : 98% and 97% OTUs clustering

# assess the fractions of shared and non shared final OTUs (after all decontamination steps, and reads normalozation)

## in pure cultures shared and unique otus 

dim(otus.table.final) 

final.shared.otus = vector()

for(i in 1:length(pure.cultures.replicated)){
  vec = vector()
  a= rownames(metadata.all.final2)[which(metadata.all.final2$Sample_name==pure.cultures.replicated[i])]
  
  for (j in 1:length(a)){
    
    b = colnames(otus.table.final)[which(otus.table.final[a[j],] > 0)]
    vec = c(vec,b)
  }
  
  final.shared.otus = c(final.shared.otus, length(which(table(vec) > 1)))
  
}

final.shared.otus = data.frame(final.shared.otus)
rownames(final.shared.otus) =  pure.cultures.replicated
colnames(final.shared.otus) = "Final_shared_OTUs"

######## non shared OTUs within each sample accross its technical replicates 

final.non.shared.otus = vector()

for(i in 1:length(pure.cultures.replicated)){
  vec = vector()
  a= rownames(metadata.all.final2)[which(metadata.all.final2$Sample_name==pure.cultures.replicated[i])] ## the technical replicate samples of pure culture i 
  
  for (j in 1:length(a)){
    
    b = colnames(otus.table.final)[which(otus.table.final[a[j],] > 0)] ## detected otus in technical replicate J 
    vec = c(vec,b) ## otus names detected in the technical replicates
  }
  
  final.non.shared.otus = c(final.non.shared.otus, length(which(table(vec) == 1)))
  
}

final.non.shared.otus = data.frame(final.non.shared.otus)
rownames(final.non.shared.otus) =  pure.cultures.replicated
colnames(final.non.shared.otus) = "Non_Shared_OTUs"

# threshold 0.02: all OTUs which lead/OTUs are within a distance 0.02 are collpased to one new 0.02_OTU  

centroid.otus = vector()
for(i in 1:ncol(otus.data)){
  a = otus.data[,i][1]
  centroid.otus = c(centroid.otus,a)
}


dm.asvs.diag.centroids = dm.asvs.diag[centroid.otus,centroid.otus]
dim(dm.asvs.diag.centroids) ## 2084 2084 

### Define the first otu.0.02 

a_1 = which(dm.asvs.diag.centroids[,1] <= 0.0200)
OTU_1 = rownames(dm.asvs.diag.centroids)[a_1]
length(OTU_1) ## 7
b_1 = setdiff(rownames(dm.asvs.diag.centroids),OTU_1)  # b are the ASVs which are NOT yet clustered into an OTU


otus.data.0.02 = data.frame(nrow=nrow(dm.asvs.diag.centroids))
otus.asvs.0.02 = OTU_1 ## vector list of all asvs already clustered into an OTU 
b_i = b_1 

### 
for (i in 2:ncol(dm.asvs.diag.centroids)){
  
  OTU_i = vector()
  
  if (length(intersect(colnames(dm.asvs.diag.centroids)[i],b_i)) == 0){ # if the length of the intersection is 0, means that this centroid otu has already been clustered into an OTU.0.02
    
    OTU_i = NULL
    otus.asvs.0.02 = otus.asvs.0.02
    otus.data.0.02 = otus.data.0.02
    b_i = b_i    # b_i is the vector of the ASVs not yet assigned/clustered to an OTU 
    
  } else { 
    
    #a_i = which(dm.asvs.diag.centroids[,b_i[1]] <= 0.0200)
    a_i = which(dm.asvs.diag.centroids[,i] <= 0.0200)
    OTU_i = rownames(dm.asvs.diag.centroids)[a_i]
    OTU_i= setdiff(OTU_i,otus.asvs.0.02) ### some centroid otus could belong to several OTUs which lead/centroid ASVs are distant more than 0.0200 They are clustered to the most abundant ASV, when several matches occur
    otus.asvs.0.02 = c(otus.asvs.0.02, OTU_i)
    OTU_i = c(OTU_i, rep("none", nrow(dm.asvs.diag.centroids) - length(OTU_i))) 
    otus.data.0.02 = cbind(otus.data.0.02, OTU_i)
    b_i = setdiff(rownames(dm.asvs.diag.centroids), otus.asvs.0.02) # b_i is the vector of the ASVs not yet assigned/clustered to an OTU 
    
  } 
  #print(length(otus.asvs.0.02))
  # when length b_i = 0 all OTUs have been defined 
  if (length(b_i) ==0) {
    print(i)
  }
  
  if(sort(table(otus.asvs.0.02), decreasing = T)[1] != 1){
    print(i)
  }
  
}

length(otus.asvs.0.02) 
length(b_i) # if 0 all centroid otus have been assigned to an OTU cluster 

otus.data.0.02 = otus.data.0.02[,-1]

otus.data.0.02 = data.frame(cbind(c(OTU_1, rep("none", nrow(dm.asvs.diag.centroids) - length(OTU_1))),otus.data.0.02))
ncol(otus.data.0.02) ## there are 1545    0.02 clustered otus 
colnames(otus.data.0.02) = paste(rep("OTU_", ncol(otus.data.0.02)),1:ncol(otus.data.0.02) ,sep="")

colnames(otus.data.0.02) = paste(colnames(otus.data.0.02), "0.02",sep="_")

write.csv(otus.data.0.02,"otus.data.0.02.csv")


## check that all ASVs within an 0.02_OTU have the same taxonomy : take the taxonomy of the centroid otus 

test = data.frame()

for (i in 1:ncol(otus.data.0.02)){
  vec= vector()
  
  a = which(otus.data.0.02[,i] != "none")
  a = otus.data.0.02[,i][a] ## a is the list of ASVs which belong to OTU i 
  
  for (j in 1:length(a)){
    b = taxonomy[a[j],]
    vec = rbind(vec,b)
  }
  
  #if (length(unique(vec$Genus)) > 1){
  if (length(unique(vec$Genus)) > 2){
    print(i)
    #print(unique(vec$Genus))
    
    test = rbind(test, c(unique(vec$Genus), rep("none", 5 - length(unique(vec$Genus))))) ## I assigned 5 as the maximum length of length(unique(vec$Genus))
    
  }
  
}

colnames(test) = c("Genus1","Genus2","Genus3","Genus4","Genus5")

## there are 2 disagreement the same as previousely noticed. Detected in only the mock community samples 

# make 0.02_OTUs taxonomz

taxonomy.otus.0.02 = data.frame()

for(i in 1:ncol(otus.data.0.02)){
  
  a = otus.data.0.02[,i][1]
  b = taxonomy[a,]
  taxonomy.otus.0.02 = rbind(taxonomy.otus.0.02,b)
  
}

dim(taxonomy.otus.0.02)
rownames(taxonomy.otus.0.02) = colnames(otus.data.0.02)

write.csv(taxonomy.otus.0.02,"taxonomy.otus.0.02.csv")


## make the 0.02?OTUs abundance tables  


otus.table.0.02.rarefied = data.frame()

for (i in 1:nrow(rarefied.data.once)){
  print(i)
  
  d_i = data.frame(nrow=ncol(otus.data)) ## dataframe of centroids asvs (corrected OTUs), and the asvs detected in sample i  
  d_f = data.frame(nrow=ncol(otus.data)) ## dataframe of otus.0.02 otus, and the asvs detcted in sample i 
  abund_i = data.frame(nrow=1)
  data_i = data.frame(rep(0, ncol(otus.data.0.02)))
  colnames(data_i) = "Sample_i"
  rownames(data_i) = colnames(otus.data.0.02)
  data_i = data.frame(t(data_i))
  
  ## first step: identify the asvs detected in sample i c
  a = which(rarefied.data.once[i,] > 0)
  a = colnames(rarefied.data.once[i,])[a] ## ASVs detected in sample i 
  
  ## second step: identify to wich corrected-asvs / otu the detected asvs belong to 
  for(j in 1:ncol(otus.data)){
    x = NULL ## asvs of centroid asvs j detected in sample i 
    vec = vector() ## dataframe of x  
    
    b = otus.data[,j][which(otus.data[,j] != "none")]
    x= intersect(b,a) 
    
    if(length(x)!= 0){
      x2 = colnames(otus.data)[j] 
      vec = data.frame(c(x, rep("none", ncol(otus.data)-length(x)))) ## x as a dataframe 
      colnames(vec) = x2
      d_i = cbind(d_i, vec) ## collection of all j´s 
    }
    
  }## j 
  
  
  if(ncol(d_i) == 2){
    x2_new = colnames(d_i)[2]
    d_i = data.frame(d_i[,-1])
    colnames(d_i) = x2_new
  } else {
    d_i = d_i[,-1]
  }
  
  ### step three: identify the lead asv of each detetced corrected.asv / otu 
  centroids_i = vector() ## the centroid asvs of the OTUs in d_i 
  
  for (m in 1:ncol(d_i)){
    x_m = colnames(d_i)[m]
    centroids_i = c(centroids_i,otus.data[1,x_m])
  } # m 
  
  
  ## step four: identify the otus.0.02 to which belong the detected asvs on sample i 
  
  for(f in 1:ncol(otus.data.0.02)){
    vec3 = vector()
    y_k = vector()
    x3 = vector()
    
    g = otus.data.0.02[,f][which(otus.data.0.02[,f] !="none")] ## centroid asvs which belong to otus.0.02 f 
    h = intersect(centroids_i,g) ## h is the asvs centroid which belong to otus.0.002 f, and detected in sample i.    Length indicate to how many otus.0.02 the cenroids belong  ################fix here when h > 1
    
    if(length(h) > 1){
      
      x3 = colnames(otus.data.0.02)[f]
      
      for (k in 1:length(h)){
        y = which(centroids_i==h[k])
        #print(length(d_i[,y][which(d_i[,y] != "none")]))
        y_k = c(y_k, d_i[,y][which(d_i[,y] != "none")])
      } ## k 
      vec3 = data.frame(c(y_k, rep("none", nrow(d_i)-length(y_k))))
      colnames(vec3)= x3
      d_f = cbind(d_f,vec3)
    } 
    
    
    if(length(h) == 1){
      x3 = colnames(otus.data.0.02)[f]
      y = which(centroids_i == h)
      vec3 = data.frame(d_i[,y])
      colnames(vec3)= x3
      d_f = cbind(d_f,vec3)
    }
  }##f
  
  
  if(ncol(d_f) ==2){
    x3_new = colnames(d_f)[2]
    d_f = data.frame(d_f[,-1])
    colnames(d_f) = x3_new
  }else{
    d_f = d_f[,-1]
  }
  
  #### step five: calculate the abudance of of each otu.0.02 withib sample i 
  for(n in 1:ncol(d_f)){
    x4 = d_f[,n][which(d_f[,n] != "none")]
    y4 = sum(rarefied.data.once[i,x4])
    abund_i = cbind(abund_i, y4)
  }##n
  
  
  abund_i = data.frame(abund_i[,-1])
  colnames(abund_i) = colnames(d_f)
  
  
  for (e in 1:ncol(abund_i)){
    x_e = which(colnames(data_i)==colnames(abund_i)[e])
    data_i[,x_e] = abund_i[,e]
  } # e 
  
  if(sum(data_i) != 10000){
    print("ERROR!")
    print(i)
  }
  
  otus.table.0.02.rarefied = rbind(otus.table.0.02.rarefied,data_i)
  
}  ## i 

rownames(otus.table.0.02.rarefied) = rownames(rarefied.data.once)
apply(otus.table.0.02.rarefied,1,sum)

write.csv(otus.table.0.02.rarefied,"otus.table.0.02.rarefied.csv")



# threshold 0.03 : all OTUs which lead/OTUs are within a distance 0.03 are collpased to one new 0.03_OTU  

# Define the first 0.03_OTU 


a_1 = which(dm.asvs.diag.centroids[,1] <= 0.0300)
OTU_1 = rownames(dm.asvs.diag.centroids)[a_1]
length(OTU_1) ## 12
b_1 = setdiff(rownames(dm.asvs.diag.centroids),OTU_1)  # b are the ASVs which are NOT yet clustered into an OTU


otus.data.0.03 = data.frame(nrow=nrow(dm.asvs.diag.centroids))
otus.asvs.0.03 = OTU_1 ## vector list of all asvs already clustered into an OTU 
b_i = b_1 

### 
for (i in 2:ncol(dm.asvs.diag.centroids)){
  
  OTU_i = vector()
  
  if (length(intersect(colnames(dm.asvs.diag.centroids)[i],b_i)) == 0){ # if the length of the intersection is 0, means that this centroid otu has already been clustered into an OTU.0.03
    
    OTU_i = NULL
    otus.asvs.0.03 = otus.asvs.0.03
    otus.data.0.03 = otus.data.0.03
    b_i = b_i    # b_i is the vector of the ASVs not yet assigned/clustered to an OTU 
    
  } else { 
    
    #a_i = which(dm.asvs.diag.centroids[,b_i[1]] <= 0.0300)
    a_i = which(dm.asvs.diag.centroids[,i] <= 0.0300)
    OTU_i = rownames(dm.asvs.diag.centroids)[a_i]
    OTU_i= setdiff(OTU_i,otus.asvs.0.03) ### some centroid otus could belong to several OTUs which lead/centroid ASVs are distant more than 0.0300 They are clustered to the most abundant ASV, when several matches occur
    otus.asvs.0.03 = c(otus.asvs.0.03, OTU_i)
    OTU_i = c(OTU_i, rep("none", nrow(dm.asvs.diag.centroids) - length(OTU_i))) 
    otus.data.0.03 = cbind(otus.data.0.03, OTU_i)
    b_i = setdiff(rownames(dm.asvs.diag.centroids), otus.asvs.0.03) # b_i is the vector of the ASVs not yet assigned/clustered to an OTU 
    
  } 
  #print(length(otus.asvs.0.03))
  # when length b_i = 0 all OTUs have been defined 
  if (length(b_i) ==0) {
    print(i)
  }
  
  if(sort(table(otus.asvs.0.03), decreasing = T)[1] != 1){
    print(i)
  }
  
}

length(otus.asvs.0.03) 
length(b_i) # if 0 all centroid otus have been assigned to an OTU cluster 

otus.data.0.03 = otus.data.0.03[,-1]

otus.data.0.03 = data.frame(cbind(c(OTU_1, rep("none", nrow(dm.asvs.diag.centroids) - length(OTU_1))),otus.data.0.03))
ncol(otus.data.0.03) ## there are 1088    0.03 clustered otus 
colnames(otus.data.0.03) = paste(rep("OTU_", ncol(otus.data.0.03)),1:ncol(otus.data.0.03) ,sep="")

colnames(otus.data.0.03) = paste(colnames(otus.data.0.03), "0.03",sep="_")

write.csv(otus.data.0.03, "otus.data.0.03.csv")




## check that all asvs within an OTU have the same taxonomy 

### take the taxonomy of the centroid otus

test = data.frame()

for (i in 1:ncol(otus.data.0.03)){
  vec= vector()
  
  a = which(otus.data.0.03[,i] != "none")
  a = otus.data.0.03[,i][a] ## a is the list of ASVs which belong to OTU i 
  
  for (j in 1:length(a)){
    b = taxonomy[a[j],]
    vec = rbind(vec,b)
  }
  
  #if (length(unique(vec$Genus)) > 1){
  if (length(unique(vec$Genus)) > 2){
    print(i)
    #print(length(unique(vec$Genus))) ## Max length(unique(vec$Genus)) is 7 
    
    test = rbind(test,  c(unique(vec$Genus), rep("none", 7 - length(unique(vec$Genus))) ))
    
  }
  
}

## 20,91,101,274,382
colnames(test) = c("Genus1","Genus2","Genus3","Genus4","Genus5", "Genus6", "Genus7")
rownames(test) = colnames(otus.data.0.03)[c(20,91,101,274,382)]

## there are 2 disagreements the same as previousely noticed. Detected in the mock community samples 
## there are addiional/novel disagreements OTU_91_0.03, and OTU_382_0.03 
## double check at end of the processing if these otus are present in the final anlyzed BALF samples 

# make the taxonomy 

taxonomy.otus.0.03 = data.frame()

for(i in 1:ncol(otus.data.0.03)){
  
  a = otus.data.0.03[,i][1]
  b = taxonomy[a,]
  taxonomy.otus.0.03 = rbind(taxonomy.otus.0.03,b)
  
}

dim(taxonomy.otus.0.03)
rownames(taxonomy.otus.0.03) = colnames(otus.data.0.03)

write.csv(taxonomy.otus.0.03, "taxonomy.otus.0.03.csv")

## make the otu.0.03 abundance tables 

otus.table.0.03.rarefied = data.frame()

for (i in 1:nrow(rarefied.data.once)){
  
  print(i)
  
  d_i = data.frame(nrow=ncol(otus.data)) ## dataframe of centroids asvs (corrected OTUs), and the asvs detected in sample i  
  d_f = data.frame(nrow=ncol(otus.data)) ## dataframe of otus.0.03 otus, and the asvs detcted in sample i 
  abund_i = data.frame(nrow=1)
  data_i = data.frame(rep(0, ncol(otus.data.0.03)))
  colnames(data_i) = "Sample_i"
  rownames(data_i) = colnames(otus.data.0.03)
  data_i = data.frame(t(data_i))
  
  ## first step: identify the asvs detected in sample i c
  a = which(rarefied.data.once[i,] > 0)
  a = colnames(rarefied.data.once[i,])[a] ## ASVs detected in sample i 
  
  ## second step: identify to wich corrected-asvs / otu the detected asvs belong to 
  for(j in 1:ncol(otus.data)){
    x = NULL ## asvs of centroid asvs j detected in sample i 
    vec = vector() ## dataframe of x  
    
    b = otus.data[,j][which(otus.data[,j] != "none")]
    x= intersect(b,a) 
    
    if(length(x)!= 0){
      x2 = colnames(otus.data)[j] 
      vec = data.frame(c(x, rep("none", ncol(otus.data)-length(x)))) ## x as a dataframe 
      colnames(vec) = x2
      d_i = cbind(d_i, vec) ## collection of all j´s 
    }
    
  }## j 
  
  
  if(ncol(d_i) == 2){
    x2_new = colnames(d_i)[2]
    d_i = data.frame(d_i[,-1])
    colnames(d_i) = x2_new
  } else {
    d_i = d_i[,-1]
  }
  
  ### step three: identify the lead asv of each detetced corrected.asv / otu 
  centroids_i = vector() ## the centroid asvs of the OTUs in d_i 
  
  for (m in 1:ncol(d_i)){
    x_m = colnames(d_i)[m]
    centroids_i = c(centroids_i,otus.data[1,x_m])
  } # m 
  
  
  ## step four: identify the otus.0.03 to which belong the detected asvs on sample i 
  
  for(f in 1:ncol(otus.data.0.03)){
    vec3 = vector()
    y_k = vector()
    x3 = vector()
    
    g = otus.data.0.03[,f][which(otus.data.0.03[,f] !="none")] ## centroid asvs which belong to otus.0.03 f 
    h = intersect(centroids_i,g) ## h is the asvs centroid which belong to otus.0.002 f, and detected in sample i.    Length indicate to how many otus.0.03 the cenroids belong  ################fix here when h > 1
    
    if(length(h) > 1){
      
      x3 = colnames(otus.data.0.03)[f]
      
      for (k in 1:length(h)){
        y = which(centroids_i==h[k])
        #print(length(d_i[,y][which(d_i[,y] != "none")]))
        y_k = c(y_k, d_i[,y][which(d_i[,y] != "none")])
      } ## k 
      vec3 = data.frame(c(y_k, rep("none", nrow(d_i)-length(y_k))))
      colnames(vec3)= x3
      d_f = cbind(d_f,vec3)
    } 
    
    
    if(length(h) == 1){
      x3 = colnames(otus.data.0.03)[f]
      y = which(centroids_i == h)
      vec3 = data.frame(d_i[,y])
      colnames(vec3)= x3
      d_f = cbind(d_f,vec3)
    }
  }##f
  
  
  if(ncol(d_f) ==2){
    x3_new = colnames(d_f)[2]
    d_f = data.frame(d_f[,-1])
    colnames(d_f) = x3_new
  }else{
    d_f = d_f[,-1]
  }
  
  #### step five: calculate the abudance of of each otu.0.03 withib sample i 
  for(n in 1:ncol(d_f)){
    x4 = d_f[,n][which(d_f[,n] != "none")]
    y4 = sum(rarefied.data.once[i,x4])
    abund_i = cbind(abund_i, y4)
  }##n
  
  
  abund_i = data.frame(abund_i[,-1])
  colnames(abund_i) = colnames(d_f)
  
  
  for (e in 1:ncol(abund_i)){
    x_e = which(colnames(data_i)==colnames(abund_i)[e])
    data_i[,x_e] = abund_i[,e]
  } # e 
  
  if(sum(data_i) != 10000){
    print("ERROR!")
    print(i)
  }
  
  otus.table.0.03.rarefied = rbind(otus.table.0.03.rarefied,data_i)
  
}  ## i 

rownames(otus.table.0.03.rarefied) = rownames(rarefied.data.once)
apply(otus.table.0.03.rarefied,1,sum)

write.csv(otus.table.0.03.rarefied, "otus.table.0.03.rarefied.csv")


### double check if the wrongly classified are still in the dataset, and in which samples 

## there are 2 disagreements the same as previousely noticed. Detected in the mock community samples 
## there are addiional/novel disagreements OTU_91_0.03, and OTU_382_0.03 
## Double check at end of the processing if these otus are present in the final anlyzed BALF samples 

which(colnames(otus.table.0.03.rarefied)=="OTU_91_0.03")## 91 

otus.table.0.03.rarefied[,91] ## abundabnce 0 accross all samples 

which(colnames(otus.table.0.03.rarefied)=="OTU_382_0.03")## 382 

otus.table.0.03.rarefied[,382] ## abundabnce 0 accross all samples 


# assess the fractions of shared 



############ for all samples, generate a summary of the detected ASVs	OTUs (adjusted ASVs), OTUs_0.02, and OTUs_0.03


vec_asvs = vector()
vec_asvs.corr = vector() 
vec_0.02 = vector()
vec_0.3 = vector()


for(i in 1:nrow(otus.table.0.03.rarefied)){
  print(i)
  
  a_asvs= length(which(rarefied.data.once[i,] > 0 ))
  vec_asvs = c(vec_asvs,a_asvs)
  
  a_asvs.corr = length(which(otus.table.final.rarefied[i,] > 0))
  vec_asvs.corr = c(vec_asvs.corr,a_asvs.corr )
  
  a_0.02 = length(which(otus.table.0.02.rarefied[i,] >0))
  vec_0.02 = c(vec_0.02, a_0.02)
  
  a_0.03 = length(which(otus.table.0.03.rarefied[i,] >0))
  vec_0.3 = c(vec_0.3,a_0.03)
  
}

contrast_abund = data.frame(cbind(vec_asvs,vec_asvs.corr,vec_0.02,vec_0.3))
rownames(contrast_abund) = rownames(otus.table.0.03.rarefied)
colnames(contrast_abund) = c("ASVs", "Corr_ASVs", "OTUs_0.02", "OTUs_0.03")

write.csv(contrast_abund, "contrast_abund.csv")














