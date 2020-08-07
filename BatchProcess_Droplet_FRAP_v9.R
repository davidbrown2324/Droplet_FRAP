#BatchProcess_Droplet_FRAP_v8.R
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
#v9 - commented for publication

#Load "pracma" for "fzero" to calculate the tHalf.
library("pracma", lib.loc="~/R/win-library/3.4")

#[1.0] Set parameters
Bleach_frame<-25
Interval<-0.3 #seconds
#assume 9 column format
col_for <- 9
rep<-"Rep1"

#Choose a desired column
column <- "FRAP...PB"

#Get the path
#path <- "Y:/DavidB/HP1_manuscript/"
path <- choose.dir(default = "Y:/DavidB/HP1_manuscript/", caption = "Select project folder")

#Choose Directory
dir<-paste0(path,"/04_Results/",rep,"/")
output<-paste0(path,"/05_Irel/",rep,"/")
plts<-paste0(path,"/06_Plots/")

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

##[1.7] Trim the normalised measurements to give a frap recovery curve
trimMeasurement<-function(Measurement, Start=(Bleach_frame+1), Stop=length(Measurement)){
  return(Measurement[Start:Stop])
}

#Define the fitting function, with optional start values
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

#Preallocate Fit - Works
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
  #Copy names from the Irel data
  BiExpFit_Frame<-as.data.frame(BiExpFit_Frame)
  names(BiExpFit_Frame)<-names(Irel_frame)
  
  return(BiExpFit_Frame)
}

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
  for (frap_experiment in 1:ncol(Irel_frame)){
    print(paste("Analysing frap_experiment", frap_experiment))
    frap <- trimMeasurement(Irel_frame[,frap_experiment])
    time <- (1:length(frap))*Interval
    biexp_fit<-try(BiExpFit(time, frap))
    thalf<-try(getThalf(biexp_fit))
    thalf_list[frap_experiment]<-thalf
  }
  
  thalf_list<-as.numeric(as.character(thalf_list)) 
  names(thalf_list)<-names(Irel_frame)
  return(thalf_list)
}

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

#To add error bars on barcharts
error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

#To add error bars to line graphs
error.line<-function(x, y, y_sd, stride, color){
  error_x0 = x[seq(to=length(x),by=stride)]
  error_y0 = y[seq(to=length(y),by=stride)]
  error_y1 = (y+y_sd)[seq(to=length(y),by=stride)]
  error_y2 = (y-y_sd)[seq(to=length(y),by=stride)]
  arrows(error_x0, 
         error_y0, 
         error_x0, 
         error_y1,
         angle=90,
         length=0.04,
         col = color)
  arrows(error_x0, 
         error_y0, 
         error_x0, 
         error_y2,
         angle=90,
         length=0.04,
         col = color)
}

#Plot average FRAP recoveries with error bars
group_plot<-function(x, y, y_error, stride){
  cl <- rainbow(dim(y)[2])
  for (i in 1:dim(y)[2]){
    if (i == 1){
      plot(x, y[,i], type='l', col=cl[i], xlim=c(0,60), ylim=c(0.7,1), xlab = "", ylab= "")
    } else {
      lines(x, y[,i], col=cl[i])
    }
    error_bars(x, y[,i], y_error[,i], stride, col=cl[i])
  }
}

####### [1] Import data and run the organiseColumns function
Irel_frame<-organiseColumns(dir)

##I want to give these shorter names
#Optional - In Exp0133 all protein is WT.
Experimental_Conditions$WT <- TRUE
for (i in 1:length(Experimental_Conditions$Filename)){
  Experimental_Conditions[i,"New Name"]<-paste0("WT+",
                                  Experimental_Conditions[i,"DNA"],
                                   "_",
                                   Experimental_Conditions[i,"Fluorophore"])
}
Experimental_Conditions$`New Name`<-as.factor(Experimental_Conditions$`New Name`)


#Calculate the average recovery per group
mean_Irel<-as.data.frame(lapply(split(as.list(test),f = colnames(test)), function(x) Reduce(`+`,x) / length(x)))

#Calculate the standard deviation per group
sd_Irel<-as.data.frame(sapply(split.default(test, names(test)), function(x) apply(x, 1, sd, na.rm=TRUE)))

#Group averages by experiment
Atto_means<-mean_Irel[,grep('488', names(mean_Irel))]
Atto_sd<-sd_Irel[,grep('488', names(sd_Irel))]

##Trim prebleach
time = (frame-25)*0.3
frap_time = time[26:325]
Atto_frap_means = Atto_means[26:325,]
Atto_frap_sd = Atto_sd[26:325,]


#Figure 5.b 
#This will plot groups, with error bars at a given spacing
figure<-"Figure_5b.pdf"
plt<-paste0(plts,figure)
pdf(plt, width=7, height=7)
group_plot(frap_time, Atto_frap_means, Atto_frap_sd, 10)
title(main="HP1 FRAP", xlab = "Time (s)", ylab = "Fluorescence Intensity")
dev.off()

#Figure 5.e 
YOYO_means<-mean_Irel[,grep('YOYO', names(mean_Irel))]
YOYO_sd<-sd_Irel[,grep('YOYO', names(sd_Irel))]
YOYO_frap_means = YOYO_means[26:325,]
YOYO_frap_sd = YOYO_sd[26:325,]

figure<-"Figure_5e.pdf"
plt<-paste0(plts,figure)
pdf(plt, width=7, height=7)
group_plot(frap_time, YOYO_frap_means, YOYO_frap_sd, 10)
title(main="HP1 FRAP", xlab = "Time (s)", ylab = "Fluorescence Intensity")
dev.off()

YOYO_index<-grep('YOYO', names(mean_Irel))
FAM_index<-grep('FAM', names(mean_Irel))

#Save
#write.csv(x=Irel_frame, file=paste0(output,paste0("Irel_frame_",rep,".csv")))

####### [2] Calculate BiExp Fit for each curve ########

#[2.1] 
BiExpFit_Frame <- Batch_Fit(Irel_frame, Bleach_frame)
##I want to give these shorter names too
Experimental_Conditions<-getConditions(names(BiExpFit_Frame))

#Save
write.csv(x=BiExpFit_Frame, file=paste0(output, paste0("BiExpFit_frame_",rep,".csv")))

####### [3] Calculate Thalf for each Fit ##############
thalf_list <- Batch_thalf(Irel_frame)

#Save
#write.csv(x=thalf_list, file=paste0(output, paste0("thalf_list_",rep,".csv")))

n_coefficients<-4

BiExpFit_Coefficients=zeros(dim(Irel_frame)[2], n_coefficients)

####### [4] Intersect FRAP Analysis with Experimental Conditions #########

##Add t1/2 column to dataframe

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

######## Process all Directories ########################

rep<-"Rep1"
dir<-paste0(path,"/04_Results/",rep,"/")
output<-paste0(path,"/05_Irel/",rep,"/")
Thalf_data_day2<-ProcessDir(dir)

rep<-"Rep2"
dir<-paste0(path,"/04_Results/",rep,"/")
output<-paste0(path,"/05_Irel/",rep,"/")
Thalf_data_day3<-ProcessDir(dir)

########## Label Data Set Results #######################

Thalf_data_day2$Day <- 2
Thalf_data_day3$Day <- 3

############ Combine Data Sets ##########################

Thalf_data_Combined <- rbind(Thalf_data_day2,Thalf_data_day3)
#Make DNA an ordered factor
Thalf_data_Combined$DNA <- factor(Thalf_data_Combined$DNA, ordered = TRUE, levels = c("150bp", "3kb", "8kb", "50kb"))

mean_thalf<-aggregate(x = Thalf_data_Combined$thalf, by = list(Thalf_data_Combined$`New Name`), FUN = 'mean', na.rm=T)
names(mean_thalf)<-c("New Name","mean_thalf")

sd_thalf<-aggregate(x = Thalf_data_Combined$thalf, by = list(Thalf_data_Combined$`New Name`), FUN = 'sd', na.rm=T)
names(sd_thalf)<-c("New Name","sd_thalf")

count_thalf<- as.numeric(by(Thalf_data_Combined$thalf, list(Thalf_data_Combined$`New Name`), length, simplify = TRUE))
count_thalf

Thalf_data_Final<-merge(mean_thalf,sd_thalf)
Thalf_data_Final['n']<-count_thalf
Thalf_data_Final['se']<-Thalf_data_Final['sd_thalf']/sqrt(Thalf_data_Final['n'])

#Reorder by DNA - This step is not transportable to other cases!!!
Thalf_data_Final<-Thalf_data_Final[c(2,1,4,8,6,3,5,9,7),]

########### Bar Plots ###################################

#Figure 5c. 488
Atto488_agg<-Thalf_data_Final[grep("488", Thalf_data_Final$`New Name`),]
figure<-"Figure_5c.pdf"
plt<-paste0(plts,figure)
pdf(plt, width=7, height=7)
atto_bar_plot<-barplot(Atto488_agg$mean_thalf, names.arg = c("150", "3k", "8k", "50k"), ylim = c(0,25))
error.bar(atto_bar_plot, Atto488_agg$mean_thalf, Atto488_agg$se)
title(main = 'WT HP1-Atto488', xlab ='DNA length (bp)', ylab='t1/2 (s)')
dev.off()

#Figure 5f. YOYO
YOYO_agg<-Thalf_data_Final[grep("YOYO", Thalf_data_Final$`New Name`),]
figure<-"Figure_5f.pdf"
plt<-paste0(plts,figure)
pdf(plt, width=7, height=7)
yoyo_bar_plot<-barplot(YOYO_agg$mean_thalf, names.arg = c("150", "3k", "8k", "50k"), ylim = c(0,25))
error.bar(yoyo_bar_plot, YOYO_agg$mean_thalf, YOYO_agg$se)
title(main = 'WT HP1-Atto488', xlab ='DNA length (bp)', ylab='t1/2 (s)')
dev.off()



