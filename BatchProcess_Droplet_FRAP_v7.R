#BatchProcess_Droplet_FRAP_v7.R
#David Brown
#2019-07-31
#v2 - Section 1 works well
#v3 - Section 2 works 
#v4 - Loads 'Results.csv' files rather than all '.csv' files to avoid loading analysis files
#v5 - BiExpFitting with for loop handles failure better
#v6 - Attempt to make the script generic
#excludes truncated files
#v7
#need to make it a script executable from command line

#Load "pracma" for "fzero" to calculate the tHalf.
library("pracma", lib.loc="~/R/win-library/3.4")

#[1.0] Set parameters
Bleach_frame<-25
Interval<-0.3 #seconds
#assume 9 column format
col_for <- 9

#Choose a desired column
column <- "FRAP...PB"

#Get the path 
path <- "Y:/DavidB/HP1_manuscript/"
path

rep<-"Rep1"

#Choose Directory
dir<-paste0(path,"/04_Results/",rep,"/")
output<-paste0(path,"/05_Irel/",rep,"/")

################################ Functions ########################

organiseColumns <- function(dir, Bleach_frame=25, col_for=9, column="FRAP...PB"){

  ##[1.1] Get file list
  files<-list.files(dir, pattern = "Results.csv", full.names = TRUE)
  file_count <- length(files)
  
  #[1.2] Read in .csv files as a list of data
  dataFiles <- lapply(Sys.glob(files), read.csv)
  
  #[1.3] Check data is what you expect
  expect_row_num <- dim(as.data.frame(dataFiles))[1]
  
  correct_col_num <- logical(length(dataFiles))
  correct_row_num <- logical(length(dataFiles))
  for (n in (1:length(dataFiles))){
    correct_col_num[n] <- ncol(as.data.frame(dataFiles[n])) == col_for
    correct_row_num[n] <- nrow(as.data.frame(dataFiles[n])) == expect_row_num
  }
  
  #[1.4] Subset dataFiles to those with the expected number of variables
  ValidDataFiles <- dataFiles[correct_col_num & correct_row_num]
  #subset names
  ValidFileNames <- files[correct_col_num & correct_row_num]
  
  #[1.5] Import data for each file
  df<-as.data.frame(ValidDataFiles)
  
  #[1.6] Extract chosen column
  
  #Check the input for column is valid
  column_valid <- is.numeric(column)|is.character(column)
  
  if (is.numeric(column)){
    print(paste("Retrieving column",column))
    #Get every nth column
    columns<-(seq(1,ncol(df), col_for)-1)+column
    columns
    Irel_frame<-df[,columns]
  } else if (is.character(column)) {
    Irel_frame<-df[,grep(pattern = column, names(df))]
  }
  
  #Convert file paths to file/ samples names
  files_names<-basename(ValidFileNames)
  
  #Trim Final Results
  s_names<-sub("_MMStack_Default.ome_Results.csv", "", files_names)
  
  #Name the Irel data
  names(Irel_frame)<-s_names

  return(Irel_frame)
}

#Run the organiseColumns function
Irel_frame<-organiseColumns(dir)

#Save
write.csv(x=Irel_frame, file=paste0(output,paste0("Irel_frame_",rep,".csv")))

##[1.7] Trim the normalised measurements to give a frap recovery curve
trimMeasurement<-function(Measurement, Start=(Bleach_frame+1), Stop=length(Measurement)){
  return(Measurement[Start:Stop])
}

####### [2] Calculate BiExp Fit for each curve ########

#[2.1] 

########## Define the fitting function, with optional start values ###########
BiExpFit<-function(time, frap, A=0.96, B=-0.136, c=-0.0008, d=-0.166){
  
  #Organise data for nls fit
  test <- data.frame(time,frap)
  
  #BiExp Fit
  biexp <- function(t, A, B, c, d, v){A*exp(c*t) + B*exp(d*t)}
  biexp_fit <- nls(frap ~ biexp(time, A, B, c, d), 
                   data=test, 
                   start=list(A=0.96, B=-0.136, c=-0.0008, d=-0.166))
  return(biexp_fit)
}

#Preallocate Fit
Batch_Fit<-function(Irel_frame, Bleach_frame){

  #Work out how long the frap measurment is.
  frap_length<-dim(Irel_frame)[1]-(Bleach_frame)

  BiExpFit_Frame<-zeros(frap_length, dim(Irel_frame)[2])
  
  for (frap_experiment in 1:length(Irel_frame)){
    #print(frap_experiment)
    frap <- trimMeasurement(Irel_frame[,frap_experiment])
    time <- (1:length(frap))*Interval
    biexp_fit<-try(BiExpFit(time, frap), silent = TRUE)
    if (is.character(biexp_fit)){
      biexp_fit<-matrix(NA,length(frap),1)
    } else {
      biexp_fit<-fitted(biexp_fit)
    }
    BiExpFit_Frame[,frap_experiment]<-biexp_fit
  }
  #Name the Irel data
  names(BiExpFit_Frame)<-names(Irel_frame)
  
  return(BiExpFit_Frame)
}

BiExpFit_Frame <- Batch_Fit(Irel_frame, Bleach_frame)

#Save
write.csv(x=BiExpFit_Frame, file=paste0(output, paste0("BiExpFit_frame_",rep,".csv")))

#Check fits
frap_length<-dim(Irel_frame)[1]-(Bleach_frame)
time <- (1:frap_length)*Interval

f=43
plot(time, Irel_frame[26:325,f])
lines(time, BiExpFit_Frame[,f], col='red')

names(BiExpFit_Frame)[f]

thalf_list <- Batch_thalf(Irel_frame)
#Save
write.csv(x=thalf_list, file=paste0(output, paste0("thalf_list_",rep,".csv")))

####### [3] Calculate Thalf for each Fit ##############

n_coefficients<-4

BiExpFit_Coefficients=zeros(dim(Irel_frame)[2], n_coefficients)

getThalf<-function(biexp_fit){
  
  #Get Coefficients
  A<-coef(biexp_fit)[1]
  B<-coef(biexp_fit)[2]
  c<-coef(biexp_fit)[3]
  d<-coef(biexp_fit)[4]
  
  #Only choose values before max recovery
  recovery <- fitted(biexp_fit)[1:which.max( fitted( biexp_fit))]
  recovery_time <- time[1:length(recovery)]
  
  #Get t1/2 of initial recovery
  #first_fit_point
  first_fit_point <- recovery[1]
  max_fit_point <- tail(recovery, 1)
  magnitude <- first_fit_point + max_fit_point
  
  #Define a starting guess as 1/2 your measurement time.
  midway <- tail(recovery_time, 1)/2
  
  #Solve the equation
  thalf <- as.numeric(fzero(function(t){2*(A*exp(c*t) + B*exp(d*t))-(magnitude)}, midway)[1])
  return(thalf)
}

Batch_thalf <- function(Irel_frame){
  thalf_list<-list()
  warning_list<-list()
  for (frap_experiment in 1:ncol(Irel_frame)){
    print(paste("Analysing frap_experiment", frap_experiment))
    frap <- trimMeasurement(Irel_frame[,frap_experiment])
    time <- (1:length(frap))*Interval
    biexp_fit<-try(BiExpFit(time, frap))
    thalf<-try(getThalf(biexp_fit))
    thalf_list[frap_experiment]<-thalf
    warning_list[frap_experiment]
  }
  
  thalf_list<-as.numeric(as.character(thalf_list)) 
  names(thalf_list)<-names(Irel_frame)
  return(thalf_list)
}

####### [4] Generate Boxplot ##########################

#Get the Well
getWell <-function(name){
  #Does not handle vectors
  split_name <- strsplit(name, split = "_")
  Well<-split_name[[1]][1]
  stopifnot(grepl("Well",Well))
  return(split_name[[1]][1])
}

getWells <-function(name_list){
  #Handles Vectors
  Wells<-list(length(name_list))
  for (i in 1:length(name_list)){
    split_name <- strsplit(name_list[i], split = "_")
    Well<-split_name[[1]][1]
    stopifnot(grepl("Well",Well))
    Wells[i]<-Well
  }
  return(unlist(Wells))
}

#Check if HP1 protein is WT
is.WT<-function(name){
  grepl("WT|wt|Wt", x = name)
}

#DNA length
DNAlength<-function(name){
  #Does not work with vector
  if (grepl(pattern = "147bp|150bp", x = name)){
    DNA <- "150bp"
  } else if (grepl(pattern = "3kb", x = name)){
    DNA <- "3kb"
  } else if (grepl(pattern = "8kb", x = name)){
    DNA <- "8kb"
  } else if (grepl(pattern = "50kb|lambda", x = name)){
    DNA <- "50kb"
  }
  return(DNA)
}

#Get Field
getField<-function(name){
  #Does not work with vector
  split_name <- strsplit(name, split = "_")
  #find the last numeric field
  numeric_info<-suppressWarnings(
    split_name[[1]][!is.na(
      as.numeric(
        split_name[[1]]))])
  return(tail(numeric_info[[1]],1))
}

#Vectorized function
getFields<-Vectorize(getField)

#What is labelled?
labelled<-function(name, dye="488"){
  #Does not work with vector
  split_name <- strsplit(name, split = "[+]|[_]")
  #find 488, or other dye
  label<-grep(dye, split_name[[1]], value = TRUE)
  label<-gsub("Atto", "", label)
  return(label)
}

fluorophore<-function(name){
  #Does not work with vector
  if (grepl(pattern = "488", x = name)){
    label <- "Atto488"
  } else if (grepl(pattern = "YOYO", x = name)){
    label <- "YOYO-1"
  } else if (grepl(pattern = "FAM", x = name)){
    label <- "FAM"
  }
  return(label)
}

getConditions<-function(name_list){
  #Get data from file names
  Experimental_Conditions=data.frame()
  for (i in 1:length(name_list)){
    name<-name_list[i]
    Experimental_Conditions[i,1]<-name
    Experimental_Conditions[i,2]<-getWell(name)
    Experimental_Conditions[i,3]<-is.WT(name)
    Experimental_Conditions[i,4]<-DNAlength(name)
    Experimental_Conditions[i,5]<-getFields(name)
    Experimental_Conditions[i,6]<-fluorophore(name)
  }
  names(Experimental_Conditions)<-c("Filename", "Well", "WT", "DNA", "Field of View", "Fluorophore")
  return(Experimental_Conditions)
}

getConditions(names(thalf_list))

####### Intersect FRAP Analysis with Experimental Conditions #########

##Simply add t1/2 column to dataframe

tHalf_data<-getConditions(names(thalf_list))
tHalf_data["thalf"]<-thalf_list
#Supply the alternative to "WT"
tHalf_data["WT"]<-factor(as.logical(unlist(tHalf_data["WT"])),levels = c(FALSE,TRUE), labels=c("Hinge","WT"))

#Title by Folder
boxplot(thalf ~ Well, tHalf_data)
title(main=dir)

#Average by well
mean_thalf<-aggregate(x = tHalf_data$thalf, by = list(tHalf_data$Well), FUN = 'mean', na.rm=T)

###################### Call The Functions ######################

#Run the organiseColumns function
Irel_frame<-organiseColumns(dir)
#BiExpFit_Frame <- Batch_Fit(Irel_frame, Bleach_frame) #optional curve fits
thalf_list <- Batch_thalf(Irel_frame)
tHalf_data<-getConditions(names(thalf_list))
tHalf_data["thalf"]<-thalf_list

#Optional - In Exp0133 all protein is WT.
for (i in 1:length(thalf_list)){
  tHalf_data[i,"New Name"]<-paste0("WT+",
                                   tHalf_data[i,"DNA"],
                                   "_",
                                   tHalf_data[i,"Fluorophore"])
}

ProcessDir <-function(dir){
  Irel_frame<-organiseColumns(dir)
  thalf_list <- Batch_thalf(Irel_frame)
  tHalf_data<-getConditions(names(thalf_list))
  tHalf_data["thalf"]<-thalf_list
  for (i in 1:length(thalf_list)){
    tHalf_data[i,"New Name"]<-paste0("WT+",
                                     tHalf_data[i,"DNA"],
                                     "_",
                                     tHalf_data[i,"Fluorophore"])
  }
  return(tHalf_data)
}

######## Process all Directories ########################

Thalf_data_day2<-ProcessDir("Y:/DavidB/Data/Exp0133_HP1_FRAP/Output Day 2/")
Thalf_data_day3<-ProcessDir("Y:/DavidB/Data/Exp0133_HP1_FRAP/Output Day 3/")

########## Label Data Set Results #######################

Thalf_data_day2$Day <- 2
Thalf_data_day3$Day <- 3

############ Combine Data Sets ##########################

Thalf_data_Combined <- rbind(Thalf_data_day2,Thalf_data_day3)
#Make DNA an ordered factor
Thalf_data_Combined$DNA <- factor(Thalf_data_Combined$DNA, ordered = TRUE, levels = c("150bp", "3kb", "8kb", "50kb"))

mean_thalf<-aggregate(x = Thalf_data_Combined$thalf, by = list(Thalf_data_Combined$`New Name`), FUN = 'mean', na.rm=T)

#488
Atto488<-Thalf_data_Combined[Thalf_data_Combined$Fluorophore=="Atto488",]
boxplot(thalf ~ DNA, Atto488)

#YOYO
YOYO<-Thalf_data_Combined[Thalf_data_Combined$Fluorophore=="YOYO-1",]
boxplot(thalf ~ DNA, YOYO, outline=F)

#FAM
FAM<-Thalf_data_Combined[Thalf_data_Combined$Fluorophore=="FAM",]
boxplot(thalf ~ DNA, FAM)

#output results
#write.csv(x=tHalf_data, file=paste0(dir,"Thalf_data.csv"))
