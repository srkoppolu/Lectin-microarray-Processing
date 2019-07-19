
############  Input Section  ##############
###########################################
###########################################

# Set the working directory (place where the input data files are stored).
# Make sure to convert all the "\" to "/" in windows (Mac and Linux most likely 
# give the desired format by itself.)
WD <- "E:/Mahal_Lab/Praveen_Scripts/Lectin_Array_Processing"

# Enter the name of the three input files.
fname.data <- "Zhongyin_Lectin_Array.txt" # Data File
fname.lectins <- "Zhongyin_PrintList.txt" # Printlist of lectins, just the list of lectins used in order.
fname.samples <- "Zhongyin_SamplesList.txt" # Samples file indicating arrays(blocks)

# Desired filename for the output file
# Do not include the ".txt" filetype at the end.
fname.output <- "Zhongyin_output_data.txt"

###########################################

# Set working directory
setwd(WD)

############  Input Preferences  ##############
###############################################
###############################################

# Number of columns in an array. (Default : 15)
n.cols <- 15

# Indicate the number of times a lectin is printed consecutively ( triplicate or replicate, etc.)
n.reps <- 3

# Decide whether to use log ratio values or the background substituted median values.
# If log ratios, enter "logvalues" and if median values, enter "medvalues"
input.type <- "medvalues"

# Select the processing mode that you would like to use ( i.e., single color ("single") or dual color ("dual"))
# If single color, specify "Cy5" or "Cy3"
process.mode <- "dual"
select.channel = "Cy5"

# Select whether median normalization is necessary.
med.norm <- "yes"

# Cut-off value for the Q test scores. (If Zmax or Zmin > z.cutoff, those lectins are removed)
z.cutoff <- 1.15

# Choice of colors for the heatmap of final output data
heatmap.color <- c("blue","gray","yellow")
clustering.row <- TRUE # TRUE if you want the lectins to be clustered in Heatmap
clustering.column <- TRUE # TRUE if you want the samples (celllines) to be clustered in the heatmap

# Display preferences

flag.display <- FALSE     # Display Flag Data 

# Remove nothing (NULL) is default, but include Good (G), Bad(B), Absent(A), 
# Not Found(N), Userdefined (U). Example : remove.flags <- c("B","A","N")
remove.flags <- NULL


##############################################################################################
##############################################################################################



############  Load Packages  ##############

# Define the required packages and install/load them using the custom built PackageLoad2 function
reqd.packages <- c("tm","reshape2","matrixStats","gplots")
source("E:/Dropbox/Git_Files/PackageLoad2.R")
PackageLoad2(reqd.packages)



###########################################




############  Defining Functions  ##############

# Define a function to print a label based on the inputs for particular flag
print.flag <- function(x){
  
  lectin <- data.lectins[as.numeric(x[1]),as.numeric(x[2])]
  sample <- data.samples[as.numeric(x[3]),2]
  label <- x[4]
  cat(paste(sample,lectin,label,sep=" :: "),sep = "\n")
  
}

# Compute the z-scores, both max and min
zmax <- function(x){
  
  y <- as.numeric(x[3:6])
  zscore <- abs((y[1]-y[3])/y[2])
  return(zscore)
  
}

zmin <- function(x){
  
  y <- as.numeric(x[3:6])
  zscore <- abs((y[1]-y[4])/y[2])
  return(zscore)
  
}

# Calculate the EPT max zscores
eptzmax <- function(y){
  
  if (as.numeric(y[7]) > z.cutoff){
    y[5]
  }else 0
}

# Calculate the EPT min zscores
eptzmin <- function(y){
  
  if (y[8] > z.cutoff){
    y[6]
  }else{
    0
  }
}

# Q-test Function
Qtest <- function(data, refdata){
  
  sample <- data[1]
  lectin <- data[2]
  Emax <- refdata$E.max[refdata$Sample == sample & refdata$Lectin == lectin]
  Emin <- refdata$E.min[refdata$Sample == sample & refdata$Lectin == lectin]
  
  if (Emax == data[3] || Emin == data[3]){
    
  }else{as.numeric(data[3])}
}

# Thresholding the Q tested mean values
threshold <- function(x){
  
  if(as.numeric(x[3]) < 1){
    100
  }else{ as.numeric(x[3]) }
  
}


colorgradient <- colorRampPalette(heatmap.color)






################################################



############  Data Processing :: Meta Data  ##############
##########################################################

# Following meta data : 

# Type
# Date and Time
# Settings
# Pixel Size
# Wavelengths
# Image Files
# PMT Gain
# Scan Power
# Laser Power
# Filters


# Read the first few lines of the data file to collect information
# regarding the microarrays as published by the software
con<- file(fname.data)
open(con)

meta.data <- readLines(con, n = 32, ok = FALSE)

close(con)

# Type
meta.type <- sub(".*=","",meta.data[3])
meta.type <- sub("\t.*", "", meta.type)

# Date and Time
meta.date <- sub(".*=","",meta.data[4])
meta.date <- sub("\t.*", "", meta.date)

# Settings
meta.settings <- sub(".*=","",meta.data[5])
meta.settings <- sub("\t.*", "", meta.settings)

# Pixel Size
meta.pixelsize <- sub(".*=","",meta.data[7])
meta.pixelsize <- sub("\t.*", "", meta.pixelsize)

# Wavelengths
meta.twavelengths <- sub(".*=","",meta.data[8])
meta.twavelengths <- sub("\t\t.*", "", meta.twavelengths)
meta.wavelengths <- paste(substr(meta.twavelengths,1,3), ", ", substr(meta.twavelengths,5,7))
rm(meta.twavelengths)

# Image Files
meta.imagefiles <- sub(".*=","",meta.data[9])
meta.imagefiles <- sub("\t.*", "", meta.imagefiles)

# PMT Gain
meta.tpmtgain <- sub(".*=","",meta.data[26])
meta.tpmtgain <- sub("\t\t.*", "", meta.tpmtgain)
meta.pmtgain <- paste(substr(meta.tpmtgain,1,3), ", ", substr(meta.tpmtgain,5,7))
rm(meta.tpmtgain)

# Scan Power
meta.tscanpower <- sub(".*=","",meta.data[27])
meta.tscanpower <- sub("\t\t.*", "", meta.tscanpower)
meta.scanpower <- paste(substr(meta.tscanpower,1,3), ", ", substr(meta.tscanpower,5,7))
rm(meta.tscanpower)

# Laser Power
meta.tlaserpower <- sub(".*=","",meta.data[28])
meta.tlaserpower <- sub("\t\t.*", "", meta.tlaserpower)
meta.laserpower <- paste(substr(meta.tlaserpower,1,4), ", ", substr(meta.tlaserpower,6,9))
rm(meta.tlaserpower)

# Filters
meta.tfilters <- sub(".*=","",meta.data[29])
meta.tfilters <- sub("\t\t.*", "", meta.tfilters)
meta.filters <- meta.tfilters
meta.filters <- sub("\t.*","",meta.filters)
meta.tfilters <- sub(".*\t","",meta.tfilters)
meta.filters <- paste(meta.filters, ", ", meta.tfilters)
rm(meta.tfilters)

# Output filename for the meta data.
# Same as the user given output filename, except _meta added at the end
fname.meta.output <- paste(fname.output,"_meta.txt",sep="")

# Write the meta data to the meta file.
sink(file = fname.meta.output, append = FALSE)
cat("Input Data File = ",fname.data,"\n")
cat("Input Lectins File = ",fname.lectins,"\n")
cat("Input Samples File = ",fname.samples,"\n")
cat("Type = ",meta.type,"\n")
cat("Date and Time = ",meta.date,"\n")
cat("Settings = ",meta.settings,"\n")
cat("Pixel Size = ",meta.pixelsize,"\n")
cat("Wavelengths = ",meta.wavelengths,"\n")
cat("Image File = ",meta.imagefiles,"\n")
cat("PMT Gain = ",meta.pmtgain,"\n")
cat("Scan Power = ",meta.scanpower,"\n")
cat("Laser Power = ", meta.laserpower,"\n")
cat("Filters = ",meta.filters,"\n")
sink()

# Remove any unnecessary objects
rm(con, meta.data, meta.type,meta.date,meta.settings,meta.pixelsize,meta.wavelengths,
   meta.imagefiles,meta.pmtgain,meta.scanpower,meta.laserpower,meta.filters)

##########################################################




############  Data Processing :: Data Initialization  ##############
####################################################################

data <- read.table(fname.data, header = T, sep = "\t", na.strings = "NA", skip = 32,
                   blank.lines.skip = T, quote = "")

#Get the column headers for the data
headers <- colnames(data)

header.test <- substr(headers[1],1,1)

if(header.test == "X"){
  
  # Get the indices for the column headers for the columns of interest (COI)
  blk.ind <- match("X.Block.",headers)
  col.ind <- match("X.Column.",headers)
  row.ind <- match("X.Row.", headers)
  flg.ind <- match("X.Flags.", headers)
  snr.635 <- match("X.SNR.635.", headers)
  snr.532 <- match("X.SNR.532.", headers)
  med.635 <- match("X.F635.Median...B635.", headers)
  med.532 <- match("X.F532.Median...B532.", headers)
  log.rat <- match("X.Log.Ratio..635.532..", headers)
  
}else{
  
  # Get the indices for the column headers for the columns of interest (COI)
  blk.ind <- match("Block",headers)
  col.ind <- match("Column",headers)
  row.ind <- match("Row", headers)
  flg.ind <- match("Flags", headers)
  snr.635 <- match("SNR.635", headers)
  snr.532 <- match("SNR.532", headers)
  med.635 <- match("F635.Median...B635", headers)
  med.532 <- match("F532.Median...B532", headers)
  log.rat <- match("Log.Ratio..635.532.", headers)
  
}

# Form a new shortened data object with the columns of interest (COI)
data.COI <- data.frame(data[,c(blk.ind, col.ind, row.ind, flg.ind, snr.635, snr.532,
                               med.635, med.532)])
colnames(data.COI) <- c("Block","Column","Row", "Flag", "SNR.635","SNR.532","Median.635",
                        "Median.532")


# Remove unnecessary objects
rm(data, blk.ind, col.ind, row.ind, flg.ind, snr.635, snr.532, med.635, med.532, log.rat, headers)


#### Lectins ####
#################

# Read the data file containig the printlist of lectins in order
input_lectins <- read.delim(fname.lectins, header = F, sep = "\t")
lectins.uniq <- input_lectins

# Number of lectins in the printlist
n.lectins <- dim(input_lectins)[1]    

# Calculate the number of rows required
n.rows <- round(n.lectins*n.reps/n.cols)   

# Generate the list of lectins with specified number of times 
lectins <- rep(input_lectins[,1], each = n.reps)

# Generate the pattern of rows and columns
rows <- rep(1:n.rows, each = n.cols, length.out = n.lectins*n.reps)
cols <- rep(1:n.cols, length.out = n.lectins*n.reps)

# Generate the final list of lectins specifying its address in rows and columns
data.lectins <- data.frame(Row = rows, Column = cols, Lectin = lectins)
colnames(data.lectins) <- c("Row","Column","Lectin")

# Re-arrange the final list of lectins in an array format
lec.melt <- melt(data.lectins, id = c("Column","Row"))
data.lectins <- cast(lec.melt, Row~Column)
data.lectins[,1] <- NULL

# Remove unwanted objects
rm(input_lectins, rows, cols, lec.melt)



#### Samples ####
#################

# Read the text file giving the list of samples in order
samples <- read.table(fname.samples, header = F, sep = "\t", na.strings = "NA", skip = 0,
                      blank.lines.skip = T, quote = "")

# The number of samples calculated from the input samples file
n.samples <- dim(samples)[1]

# Generate the block numbers in order, for the samples
blocks <- c(1:n.samples)

# Create the data frame containing the block number and the sample name
data.samples <- data.frame( Block = blocks, Sample = samples)
colnames(data.samples) <- c("Block","Sample")


# Generate a list of total lectins and samples for all samples and arrays.
lectins.total.list <- rep(lectins, times = n.samples)
samples.total.list <- rep(samples[,1], each = length(lectins))


################################################################




############  Data Processing :: Flags  ##############
######################################################

if (flag.display){
  
  # Define Flag values
  good.flag.val <- 100
  bad.flag.val <- -100
  NA.flag.val <- -50
  absent.flag.val <- -75
  unflag.flag.val <- 0
  
  # Isolate the flag values from the data
  data.flag <- data.frame(data.COI$Flag)
  colnames(data.flag) <- c("Flag")
  
  # Change the flag value to flag category based on defined values
  data.flag[data.flag$Flag == good.flag.val,] = "Good"
  data.flag[data.flag$Flag == bad.flag.val,] = "Bad"
  data.flag[data.flag$Flag == NA.flag.val,] = "Not Found"
  data.flag[data.flag$Flag == absent.flag.val,] = "Absent"
  data.flag[data.flag$Flag == unflag.flag.val,] = "Unflagged"
  data.flag[data.flag$Flag != "Good" & data.flag$Flag != "Bad" & data.flag$Flag != "Not Found" & data.flag$Flag != "Absent" & data.flag$Flag != "Unflagged",] = "User-defined"
  
  # Filename for the flag data output file
  fname.flag.output <- paste(fname.output, "_flag.txt", sep="")
  
  if( any(data.flag$Flag == "Good")){
    disp.good <- data.COI[data.flag$Flag == "Good", c("Row","Column","Block")]
    disp.good$label <- 'Good'
  }else{disp.good <- NULL}
  
  if( any(data.flag$Flag == "Bad")){
    disp.bad <- data.COI[data.flag$Flag == "Bad", c("Row","Column","Block")]
    disp.bad$label <- 'Bad'
  }else{disp.bad <- NULL}
  
  if( any(data.flag$Flag == "Not Found")){
    disp.notfound <- data.COI[data.flag$Flag == "Not Found", c("Row","Column","Block")]
    disp.notfound$label <- 'Not Found'
  }else{disp.notfound <- NULL}
  
  if( any(data.flag$Flag == "Absent")){
    disp.absent <- data.COI[data.flag$Flag == "Absent", c("Row","Column","Block")]
    disp.absent$label <- 'Absent'
  }else{disp.absent <- NULL}
  
  if( any(data.flag$Flag == "User-defined")){
    disp.userdefined <- data.COI[data.flag$Flag == "User-defined", c("Row","Column","Block")]
    disp.userdefined$label <- 'User-defined'
  }else{disp.userdefined <- NULL}
  
  # Combine all the flag entries into a single object
  disp.flag <- rbind(disp.good,disp.bad,disp.notfound,disp.absent,disp.userdefined)
  
  
  # Write the flag data to the flag file.
  sink(file = fname.flag.output, append = FALSE)
  apply(disp.flag,1,print.flag)
  sink()
  
  
  
  
}

######################################################



############  Data Processing :: Data Processing  ##############
################################################################

############ Step 1 : Original Values
############ Step 2 : Z-score values
############ Step 3 : Q-test values
############ Step 4 : Q-test passed values

if (process.mode == "dual"){
  
  
  if (input.type == "medvalues"){
    
    # Step 1 is to aggregate the triplicate values for all samples and lectin to 
    # find the mean, standard deviation, maximum and minimum values.
    data.step1 <- data.frame(Sample = samples.total.list, Lectin = lectins.total.list,
                             Signal1 = data.COI$Median.635, Signal2 = data.COI$Median.532)
    colnames(data.step1) <- c("Sample","Lectin","Signal.635","Signal.532")
    
    # Calculate the mean, standard deviation, maximum and minimum values that are used in Z-score calculation
    data.step1.means <- aggregate(data.step1[,3:4], by = list(Lectin = data.step1$Lectin, Sample = data.step1$Sample), mean)
    data.step1.sds <- aggregate(data.step1[,3:4], by = list(Lectin = data.step1$Lectin, Sample = data.step1$Sample), colSds)
    data.step1.maxs <- aggregate(data.step1[,3:4], by = list(Lectin = data.step1$Lectin, Sample = data.step1$Sample), max)
    data.step1.mins <- aggregate(data.step1[,3:4], by = list(Lectin = data.step1$Lectin, Sample = data.step1$Sample), min)
    
    # Step 2 is to calculate the z scores for all the lectins and samples
    data.step2.635 <- data.frame(data.step1.means[,1:3], data.step1.sds[,3], data.step1.maxs[,3], data.step1.mins[,3])
    colnames(data.step2.635) <- c("Lectin","Sample","Mean","SD","Max","Min")
    
    data.step2.532 <- data.frame(data.step1.means[,c(1:2,4)], data.step1.sds[,4], data.step1.maxs[,4], data.step1.mins[,4])
    colnames(data.step2.532) <- c("Lectin","Sample","Mean","SD","Max","Min")
    
    # Computing the Z- scores
    data.step2.635$Z.max <- apply(data.step2.635, 1, zmax)
    data.step2.532$Z.max <- apply(data.step2.532, 1, zmax)
    data.step2.635$Z.min <- apply(data.step2.635, 1, zmin)
    data.step2.532$Z.min <- apply(data.step2.532, 1, zmin)
    
    # Changing all the NAs in Zmax and Zmin to 0 ( source of NAs : SD = 0)
    data.step2.635[data.step2.635$SD == 0, 7:8] <- 0
    data.step2.532[data.step2.532$SD == 0, 7:8] <- 0
    
    # Calulating the EPT max and min scores.
    data.step2.635$E.max <- apply(data.step2.635, 1, eptzmax)
    data.step2.635$E.min <- apply(data.step2.635, 1, eptzmin)
    data.step2.532$E.max <- apply(data.step2.532, 1, eptzmax)
    data.step2.532$E.min <- apply(data.step2.532, 1, eptzmin)
    
    # Calculating the Qtest values
    data.step3.635 <- data.frame(data.step1[,1:3], data.COI$SNR.635)
    colnames(data.step3.635)[3:4] <- c("Signal", "SNR")
    
    data.step3.532 <- data.frame(data.step1[,c(1,2,4)], data.COI$SNR.532)
    colnames(data.step3.532)[3:4] <- c("Signal", "SNR")
    
    data.step3.635$qtest <- apply(data.step3.635, 1, Qtest, refdata = data.step2.635)
    data.step3.532$qtest <- apply(data.step3.532, 1, Qtest, refdata = data.step2.532)
    
    # Remove the Qtest failed lectins
    data.step3.635 <- data.step3.635[!(data.step3.635$qtest == "NULL"),]
    data.step3.532 <- data.step3.532[!(data.step3.532$qtest == "NULL"),]
    
    
    ########
    
    
    #  >>>>> Condition for SNR based cutoff
    
    
    ########
    
    # Find the mean and standard deviation of Qtest passed lectins
    data.step4.635.means <- aggregate(data.step3.635$Signal, by = list(Lectin = data.step3.635$Lectin, Sample = data.step3.635$Sample), mean)
    data.step4.532.means <- aggregate(data.step3.532$Signal, by = list(Lectin = data.step3.532$Lectin, Sample = data.step3.532$Sample), mean)
    
    data.step4.635.sds <- aggregate(data.step3.635$Signal, by = list(Lectin = data.step3.635$Lectin, Sample = data.step3.635$Sample), sd)
    data.step4.532.sds <- aggregate(data.step3.532$Signal, by = list(Lectin = data.step3.532$Lectin, Sample = data.step3.532$Sample), sd)
    
    data.step4.635 <- data.frame(data.step4.635.means, data.step4.635.sds[,3])
    colnames(data.step4.635) <- c("Lectin","Sample","Mean","SD")
    data.step4.532 <- data.frame(data.step4.532.means, data.step4.532.sds[,3])
    colnames(data.step4.532) <- c("Lectin","Sample","Mean","SD")
    
    # Threshold values for both the channels.
    data.step4.635$Thresh <- apply(data.step4.635, 1, threshold)
    data.step4.532$Thresh <- apply(data.step4.532, 1, threshold)
    
    if (med.norm == "yes"){
      
      # Apply log2 to the threshold values of each channel
      data.step4.635$LogThresh <- log2(data.step4.635$Thresh)
      data.step4.532$LogThresh <- log2(data.step4.532$Thresh)
      
      # Calculate the median for the Log Threshold values for each channel
      LogThresh.Median.635 <- aggregate(data.step4.635$LogThresh, by = list(Sample = data.step4.635$Sample), median)
      colnames(LogThresh.Median.635) <- c("Sample","Median")
      LogThresh.Median.532 <- aggregate(data.step4.532$LogThresh, by = list(Sample = data.step4.532$Sample), median)
      colnames(LogThresh.Median.532) <- c("Sample","Median")
      
      LogThresh.Median.635 <- rep(LogThresh.Median.635$Median, each = as.numeric(n.lectins))
      LogThresh.Median.532 <- rep(LogThresh.Median.532$Median, each = as.numeric(n.lectins))
      
      # Perform median centering (Median Normalization) by subtracting the median value from each value in both channels
      data.step4.635$MedCen.LogThresh <- data.step4.635$LogThresh - LogThresh.Median.635
      data.step4.532$MedCen.LogThresh <- data.step4.532$LogThresh - LogThresh.Median.532
      
      # Data in step5 is the final data in a premature data.frame format with log(635/532) channels
      data.step5 <- data.step4.635[,1:2]
      data.step5$LogRatio <- data.step4.635$MedCen.LogThresh - data.step4.532$MedCen.LogThresh
      
      # Finally rearrange the data into an array such that lectins are along the rows and samples across columns.
      data.step5.melt <- melt(data.step5, id = c("Lectin","Sample"))
      data.step5.array <- cast(data.step5.melt, Lectin~Sample)
      
      write.table(data.step5.array, fname.output, row.names = FALSE, col.names = TRUE, sep = "\t")
      
    }
    
    if (med.norm == "no"){
      
      data.step5 <- data.step4.635[,1:2]
      data.step5.ratios <- data.step4.635$Thresh/data.step4.532$Thresh
      data.step5$logratios <- log2(data.step5.ratios)
      
      # Finally rearrange the data into an array such that lectins are along the rows and samples across columns.
      data.step5.melt <- melt(data.step5, id = c("Lectin","Sample"))
      data.step5.array <- cast(data.step5.melt, Lectin~Sample)
      
      write.table(data.step5.array, fname.output, row.names = FALSE, col.names = TRUE, sep = "\t")
      
      
    }
    
  }
  
  if (input.type == "logvalues"){
    
    # Step 1 is to aggregate the triplicate values for all samples and lectin to 
    # find the mean, standard deviation, maximum and minimum values.
    data.step1 <- data.frame(Sample = samples.total.list, Lectin = lectins.total.list,
                             Signal = data.COI$Log.Ratio.635.532)
    colnames(data.step1) <- c("Sample","Lectin","Signal")
    
    # Calculate the mean, standard deviation, maximum and minimum values that are used in Z-score calculation
    data.step1.means <- aggregate(as.numeric(data.step1[,3]), by = list(Lectin = data.step1$Lectin, Sample = data.step1$Sample), mean)
    data.step1.sds <- aggregate(as.numeric(data.step1[,3]), by = list(Lectin = data.step1$Lectin, Sample = data.step1$Sample), sd)
    data.step1.maxs <- aggregate(as.numeric(data.step1[,3]), by = list(Lectin = data.step1$Lectin, Sample = data.step1$Sample), max)
    data.step1.mins <- aggregate(as.numeric(data.step1[,3]), by = list(Lectin = data.step1$Lectin, Sample = data.step1$Sample), min)
    
    # Step 2 is to calculate the z scores for all the lectins and samples
    data.step2 <- data.frame(data.step1.means[,1:3], data.step1.sds[,3], data.step1.maxs[,3], data.step1.mins[,3])
    colnames(data.step2) <- c("Lectin","Sample","Mean","SD","Max","Min")
    
    # Computing the Z- scores
    data.step2$Z.max <- apply(data.step2, 1, zmax)
    data.step2$Z.min <- apply(data.step2, 1, zmin)
    
    # Changing all the NAs in Zmax and Zmin to 0 ( source of NAs : SD = 0)
    data.step2[data.step2$SD == 0, 7:8] <- 0
    
    # Calulating the EPT max and min scores.
    data.step2$E.max <- apply(data.step2, 1, eptzmax)
    data.step2$E.min <- apply(data.step2, 1, eptzmin)
    
    # Calculating the Qtest values
    data.step3 <- data.frame(data.step1[,1:3], data.COI$SNR.635)
    colnames(data.step3)[3:4] <- c("Signal", "SNR")
    
    # Perform the Q test
    data.step3$qtest <- apply(data.step3, 1, Qtest, refdata = data.step2)
    
    # Remove the Qtest failed lectins
    data.step3 <- data.step3[!(data.step3$qtest == "NULL"),]
    
    
    ########
    
    
    #  >>>>> Condition for SNR based cutoff
    
    
    ########
    
    # Find the mean and standard deviation of Qtest passed lectins
    data.step4.means <- aggregate(data.step3$Signal, by = list(Lectin = data.step3$Lectin, Sample = data.step3$Sample), mean)
    data.step4.sds <- aggregate(data.step3$Signal, by = list(Lectin = data.step3$Lectin, Sample = data.step3$Sample), sd)
    
    data.step4 <- data.frame(data.step4.means, data.step4.sds[,3])
    colnames(data.step4) <- c("Lectin","Sample","Mean","SD")
    
    # Data in step5 is the final data in a premature data.frame format with log(635/532) channels
    data.step5 <- data.step4[,c(1:3)]
    
    # Finally rearrange the data into an array such that lectins are along the rows and samples across columns.
    data.step5.melt <- melt(data.step5, id = c("Lectin","Sample"))
    data.step5.array <- cast(data.step5.melt, Lectin~Sample)
    
    write.table(data.step5.array, fname.output, row.names = FALSE, col.names = TRUE, sep = "\t")
    
  }
  
  
}



if (process.mode == "single"){
  
  
  # Step 1 is to aggregate the triplicate values for all samples and lectin to 
  # find the mean, standard deviation, maximum and minimum values.
  
  if (select.channel == "Cy5"){
    
    data.step1 <- data.frame(Sample = samples.total.list, Lectin = lectins.total.list,
                             Signal = data.COI$Median.635)
    colnames(data.step1) <- c("Sample","Lectin","Signal")
    
  }
  
  if (select.channel == "Cy3"){
    
    data.step1 <- data.frame(Sample = samples.total.list, Lectin = lectins.total.list,
                             Signal = data.COI$Median.532)
    colnames(data.step1) <- c("Sample","Lectin","Signal")
    
  }
  
  # Calculate the mean, standard deviation, maximum and minimum values that are used in Z-score calculation
  data.step1.means <- aggregate(as.numeric(data.step1[,3]), by = list(Lectin = data.step1$Lectin, Sample = data.step1$Sample), mean)
  data.step1.sds <- aggregate(as.numeric(data.step1[,3]), by = list(Lectin = data.step1$Lectin, Sample = data.step1$Sample), sd)
  data.step1.maxs <- aggregate(as.numeric(data.step1[,3]), by = list(Lectin = data.step1$Lectin, Sample = data.step1$Sample), max)
  data.step1.mins <- aggregate(as.numeric(data.step1[,3]), by = list(Lectin = data.step1$Lectin, Sample = data.step1$Sample), min)
  
  # Step 2 is to calculate the z scores for all the lectins and samples
  data.step2 <- data.frame(data.step1.means[,1:3], data.step1.sds[,3], data.step1.maxs[,3], data.step1.mins[,3])
  colnames(data.step2) <- c("Lectin","Sample","Mean","SD","Max","Min")
  
  # Computing the Z- scores
  data.step2$Z.max <- apply(data.step2, 1, zmax)
  data.step2$Z.min <- apply(data.step2, 1, zmin)
  
  # Changing all the NAs in Zmax and Zmin to 0 ( source of NAs : SD = 0)
  data.step2[data.step2$SD == 0, 7:8] <- 0
  
  # Calulating the EPT max and min scores.
  data.step2$E.max <- apply(data.step2, 1, eptzmax)
  data.step2$E.min <- apply(data.step2, 1, eptzmin)
  
  # Calculating the Qtest values
  data.step3 <- data.frame(data.step1[,1:3], data.COI$SNR.635)
  colnames(data.step3)[3:4] <- c("Signal", "SNR")
  
  # Perform the Q test
  data.step3$qtest <- apply(data.step3, 1, Qtest, refdata = data.step2)
  
  # Remove the Qtest failed lectins
  data.step3 <- data.step3[!(data.step3$qtest == "NULL"),]
  
  
  ########
  
  
  #  >>>>> Condition for SNR based cutoff
  
  
  ########
  
  # Find the mean and standard deviation of Qtest passed lectins
  data.step4.means <- aggregate(data.step3$Signal, by = list(Lectin = data.step3$Lectin, Sample = data.step3$Sample), mean)
  data.step4.sds <- aggregate(data.step3$Signal, by = list(Lectin = data.step3$Lectin, Sample = data.step3$Sample), sd)
  
  data.step4 <- data.frame(data.step4.means, data.step4.sds[,3])
  colnames(data.step4) <- c("Lectin","Sample","Mean","SD")
  
  # Data in step5 is the final data in a premature data.frame format with log(635/532) channels
  data.step5 <- data.step4[,c(1:3)]
  
  # Finally rearrange the data into an array such that lectins are along the rows and samples across columns.
  data.step5.melt <- melt(data.step5, id = c("Lectin","Sample"))
  data.step5.array <- cast(data.step5.melt, Lectin~Sample)
  
  # Write output data into a text file with the user defined filename in fname.output
  write.table(data.step5.array, fname.output, row.names = FALSE, col.names = TRUE, sep = "\t")
  
  
}

## Generate and save a heatmap of the array in the working directory
png(paste(substring(fname.output,1,nchar(fname.output)-4), "_heatmap.png", sep = ""), width = 4800, height = 3600, res = 300)
heatmap.2(as.matrix(data.step5.array), Colv = clustering.column, Rowv = clustering.row, scale = 'none', col = colorgradient(16), trace = "none", tracecol = "red", cexRow = 0.6, density.info = "density", densadj = 0.5)
dev.off()





########################################################
