//01_BatchFRAP_Script.ijm

//Define fixed input and output folders
//Select tiff file
//input = getDirectory("Choose a Directory");

//Or
//input = "Y:/Shuqin/Clover/Microscopy/20190821/00analysis/Images Average of 10 frames/"

//Define the script you want to run
//Or choose folder containing all scripts
script = File.openDialog("Choose an .ijm script to run")

//script = "C:/Users/David Brown/Desktop/Shuqin_2019-08-22_MaskNuclei_Script.ijm"

print(script)

//list all files in input folder (should all be .nd2s)
filelist = getFileList(input)

//for loop, iterating through folder 
for (i=0; i<filelist.length; i++){

	run(script+"...", filelist[i])

}