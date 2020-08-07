/*
 * David Brown
 * 2/Feb/2019
 */

name=getInfo("image.filename");
Maskname=replace(name, ".tif", "_Mask.tif");
ROIname=replace(name, ".tif", "_ROIs.zip");
Resultsname=replace(name, ".tif", "_Results.csv");

//Choose the bleach ROI used
FRAPname="E://DAVE_BROWN_RAW_DATA_BACKUP/2019/Input Day 3/FRAPline_3pix.roi"
BleachFrame=25;

//Use this output directory
output="E://DAVE_BROWN_RAW_DATA_BACKUP/2019/Output Day 3/";

run("Duplicate...", "duplicate");
//Maybe Try Smoothing
run("Smooth", "stack");
setAutoThreshold("Huang dark");
//run("Threshold...");
setOption("BlackBackground", true);
run("Convert to Mask", "method=Huang background=Dark calculate black");
run("Fill Holes", "stack");

/*May not need this
run("Erode", "stack");
run("Erode", "stack");
run("Erode", "stack");
run("Dilate", "stack");
run("Dilate", "stack");
run("Dilate", "stack");
*/
 
run("Invert", "stack");

saveAs("tiff", output+Maskname)

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
roiManager("Add"); //"Add" does not "Select"
min_ds=roiManager("Count")-1;
roiManager("OR");
Roi.setName("max_Droplet_set");
roiManager("Add");
max_ds=roiManager("Count")-1;
run("Make Inverse");
Roi.setName("Background");
roiManager("Add");
bgd=roiManager("Count")-1;

//Load FRAP region
roiManager("Open", FRAPname);
raw_frap=roiManager("Count")-1;

//Calculate unbleached region
roiManager("Select", raw_frap);
run("Make Inverse");
Roi.setName("Unbleached");
roiManager("Add");
ub=roiManager("Count")-1;

//Calculate Final FRAP ROI set (always in a droplet)
roiManager("deselect")
roiManager("Select",newArray(raw_frap, min_ds)); // select raw_frap and min_ds
roiManager("AND");
Roi.setName("final_FRAP");
roiManager("Add");
final_FRAP=roiManager("Count")-1;

//Calculate Final Control ROI set (always in a droplet but never bleached)
roiManager("deselect")
roiManager("Select",newArray(min_ds, ub)); // select min_ds and unbleached
roiManager("AND");
Roi.setName("final_Control");
roiManager("Add");
final_control=roiManager("Count")-1;

//Calculate Final Background (never in a droplet and never bleached)
roiManager("deselect")
roiManager("Select", newArray(bgd, ub)); // select background and unbleached
roiManager("AND");
Roi.setName("final_Background");
roiManager("Add");
final_bgd=roiManager("Count")-1;

//Need to split Final FRAP into distinct ROIs.

//Save based on filename
roiManager("Select", 0);
run("Select All");
roiManager("Save", output+ROIname);

//Select original image by closing active Mask
close();

//Measure
//roiManager("deselect");
roiManager("Select", newArray(roiManager("Count")-3, roiManager("Count")-2, roiManager("Count")-1)); 
//Set Measurements
run("Set Measurements...", "mean redirect=None decimal=3");
run("Clear Results");

//Multimeasure
roiManager("Multi Measure");
selectWindow("Results");

//PROCESS MEASURMENTS
//Subtract background 'BGD'
for (i = 0; i < nResults; i++) {
setResult("FRAP - BGD", i, getResult("Mean(final_FRAP)", i)-getResult("Mean(final_Background)", i));
setResult("CONTROL - BGD", i, getResult("Mean(final_Control)", i)-getResult("Mean(final_Background)", i));
}
updateResults();

//Calculate Initial Conditions
FRAPi = newArray(BleachFrame-1);
CONTROLi = newArray(BleachFrame-1);

for(i=0; i<BleachFrame-1; i++) {
	FRAPi[i] = getResult("FRAP - BGD", i);
	CONTROLi[i] = getResult("CONTROL - BGD", i);
}
Array.getStatistics(FRAPi, min, max, mean, std);
FRAP_AVG=mean;
Array.getStatistics(CONTROLi, min, max, mean, std);
CONTROL_AVG=mean;
//print(FRAP_AVG);
//print(CONTROL_AVG);

//Normalise for Initial Conditions
for (i = 0; i < nResults; i++) {
setResult("FRAP / Initial", i, getResult("FRAP - BGD", i)/FRAP_AVG);
setResult("CONTROL / Initial", i, getResult("CONTROL - BGD", i)/CONTROL_AVG);

//Correct for Photobleaching
setResult("FRAP - PB", i, getResult("FRAP / Initial", i)/getResult("CONTROL / Initial", i));
}
updateResults();

//Save Results
saveAs("Measurements", output+Resultsname);
