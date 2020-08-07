/*FRAP_Normalisation1.ijm
David Brown
2019-11-03

This is an imageJ macro to process multimeasure results from FRAP ROIs.
It expects "final_FRAP", "final_Background" and "final_Control" as columns in the imageJ Results window.
*/


//Subtract background 'BGD'
for (i = 0; i < nResults; i++) {
setResult("FRAP - BGD", i, getResult("Mean(final_FRAP)", i)-getResult("Mean(final_Background)", i));
setResult("CONTROL - BGD", i, getResult("Mean(final_Control)", i)-getResult("Mean(final_Background)", i));
}
updateResults();

//Calculate Initial Conditions
FRAPi = newArray(8);
CONTROLi = newArray(8);

for(i=0; i<8; i++) {
	FRAPi[i] = getResult("FRAP - BGD", i);
	CONTROLi[i] = getResult("CONTROL - BGD", i);
}
Array.getStatistics(FRAPi, min, max, mean, std);
FRAP_AVG=mean;
Array.getStatistics(CONTROLi, min, max, mean, std);
CONTROL_AVG=mean;
print(FRAP_AVG);
print(CONTROL_AVG);

//Normalise for Initial Conditions
for (i = 0; i < nResults; i++) {
setResult("FRAP / Initial", i, getResult("FRAP - BGD", i)/FRAP_AVG);
setResult("CONTROL / Initial", i, getResult("CONTROL - BGD", i)/CONTROL_AVG);

//Correct for Photobleaching
setResult("FRAP - PB", i, getResult("FRAP / Initial", i)/getResult("CONTROL / Initial", i));
}
updateResults();
