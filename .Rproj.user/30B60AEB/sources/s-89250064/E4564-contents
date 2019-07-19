

############  Input Section  ##############

# Set the working directory (place where the input data files are stored).
# Make sure to convert all the "\" to "/" in windows (Mac and Linux most likely 
# give the desired format by itself.)
WD <- "E:/Mahal_Lab/Praveen_Scripts/Lectin_Array_Processing"

# Enter the name of the three input files.
fname.data <- "input_data.txt" # Data File
fname.lectins <- "input_lectins.txt" # Lectins file indicating row/column
fname.samples <- "input_samples.txt" # Samples file indicating arrays(blocks)


###########################################

# Set working directory
setwd(WD)

############  Input Preferences  ##############

# Desired filename for the output file
fname.output <- "Output_data.txt"





###############################################



############  Load Packages  ##############

# install.packages("tm")
# library("tm")




###########################################




############  Defining Functions  ##############







################################################



############  Data Processing :: Meta Data Processing  ##############

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


# Remove any unnecessary objects
rm(con, meta.data)






######################################################################



############  Data Processing :: Block 2  ##############







########################################################

