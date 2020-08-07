/* DropletFRAPmacro1
 *  David Brown
 *  3rd of March 2018
 */

//Select tiff file
path = File.openDialog("Choose a .tiff File");
//print(path)

//Choose an output directory
output=getDirectory("Choose output folder");
//print(output)


open(path)
setSlice(10) //Show the Bleach Frame (in this case 10)
file=getInfo("image.filename");
//print(file)

//Remove Micromanager file name tag (last 21 characters in this case)
new_file=substring(file, 0, lengthOf(file)-21);
//print(new_file)
saveAs("tiff", output+new_file)

//Make Mask
//Copy image
run("Duplicate...", "duplicate");
//filter high and low frequencies
run("Bandpass Filter...", "filter_large=210 filter_small=3 suppress=None tolerance=5 autoscale saturate process");
//threshold the image
setAutoThreshold("Huang dark");
run("Convert to Mask", "method=Huang background=Dark");
run("Fill Holes", "stack");
run("Erode", "stack");
run("Erode", "stack");
run("Erode", "stack");
run("Dilate", "stack");
run("Dilate", "stack");
run("Dilate", "stack");

//Save Mask
saveAs("tiff", output+new_file+"_Nice_Mask")

/*Get filename
name=getInfo("image.filename");
ROIname=replace(name, ".tif", "_ROIs.zip");
*/

//Clear ROI manager if autosaving ROIs!!
roiManager("deselect")
roiManager("delete")

//loop through stack to detect droplets in each frame
for (i = 1; i < nSlices+1; i++) {
		setSlice(i);				//Go to next slice
        run("Create Selection");	//Generate droplet area ROI
        Roi.setName("Frame"+i);		//Sensible Name
		roiManager("Add");			//Add to manager
}

//Combine droplet ROIs to get max, min and background
roiManager("select", Array.getSequence(nSlices));	//This is by far the best way to do this
roiManager("AND");
Roi.setName("min_Droplet_set");
roiManager("Add");
min_ds=roiManager("index");
roiManager("OR");
Roi.setName("max_Droplet_set");
roiManager("Add");
max_ds=roiManager("index");
run("Make Inverse");
Roi.setName("Background");
roiManager("Add");
bgd=roiManager("index");

/*Save based on filename
roiManager("Select", 0);
run("Select All");
roiManager("Save", ROIname);
*/

roiManager("Open", "Y:/DavidB/Data/Exp0072-HP1_FRAP/Day3/Combined_and_Unbleached_3x1_RoiSet.zip");
//open("Y:/DavidB/Data/Exp0072-HP1_FRAP/Day3/Combined_and_Unbleached_3x1_RoiSet.zip"); not necessary

//Calculate Final FRAP ROI set (always in a droplet)
roiManager("deselect");
ub=roiManager("Count")-1
raw_frap=roiManager("Count")-2;
bgd=roiManager("Count")-3;
min_ds=roiManager("Count")-5;
roiManager("Select",newArray(raw_frap, min_ds)); // select raw_frap and min_ds
roiManager("AND");
Roi.setName("final_FRAP");
roiManager("Add");

//Calculate Final Control ROI set (always in a droplet but never bleached)
roiManager("deselect")
roiManager("Select",newArray(min_ds, ub)); // select min_ds and unbleached
roiManager("AND");
Roi.setName("final_Control");
roiManager("Add");

//Calculate Final Background (never in a droplet and never bleached)
roiManager("deselect")
roiManager("Select", newArray(bgd, ub)); // select background and unbleached
roiManager("AND");
Roi.setName("final_Background");
roiManager("Add");

//Save RoiSet
roiManager("save", output+new_file+"_Final_Roiset.zip")

//Select original image by closing active Mask
close()

//Measure
roiManager("deselect")
roiManager("Select", newArray(roiManager("Count")-3, roiManager("Count")-2, roiManager("Count")-1)); 
//Set Measurements
run("Set Measurements...", "mean redirect=None decimal=3");
run("Clear Results");
 
//Multimeasure
roiManager("Multi Measure");
selectWindow("Results");
saveAs("Measurements", output+new_file+"_Final_Results.csv");

//Close Results Window
selectWindow("Results"); 
     run("Close");

//Close Image      
close()
