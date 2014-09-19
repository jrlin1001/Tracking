//script for creating mask and measure for live imaging
//Jerry 2014/09/15
// for Operetta imaging
// optimize segmenation for operetta imaging

rowarray = newArray("000","001","002","003","004","005","006","007","008");
colarray =newArray("000","001","002","003","004","005","006","007","008","009","010","011","012");

for (site=2;site<61;site++)
{

//----------------------------------------------------------------
start = getTime;
while (nImages>0) { 
          selectImage(nImages);
          close(); 
}

start = getTime;
if (isOpen("ROI Manager")) { 
      selectWindow("ROI Manager"); 
      run("Close");
} 
setBatchMode(true); 

//--open image sequence of nuclear mask and duplicate--

run("Image Sequence...", "open=Z:\\sorger\\data\\NIC\\Jerry\\20140916FOXO-EGFIGF\\20140916-FOXO-EGFIGF_w1YFP_s45_t197.TIF number=29161 starting=1 increment=1 scale=100 file=w2RFP_s"+site+"_ sort");

//run("Image Sequence...", "open=Z:\\sorger\\data\\Operetta\\Jerry\\20140826-FOXO-live-EGF-insulin[1371]\\20140826-live-2[5948]\\2014-08-26T154535-0400[7197]\\003010-1-034001001.tif increment=2 file="+rowarray[row]+colarray[col]+"-1- sort");

run("Subtract Background...", "rolling=25 stack");
setOption("BlackBackground", false);
rename("original");
run("Duplicate...", "title=mask duplicate range=1-300");
//---------------------------------------------------
print("Site ",site," Building Stacks",(getTime-start)/1000); 


//--processing to binary mask--
selectWindow("mask");

run("Gaussian Blur...", "sigma=1.5 stack");
run("Enhance Contrast...", "saturated=0.5");
run("Make Binary", "method=Li background=Dark calculate");
run("Fill Holes", "stack");
run("Watershed", "stack");


rename("nuclear-mask");
//------------------------------------
print("Site ",site," Convert to Mask",(getTime-start)/1000); 

//--analysis particle---
selectWindow("nuclear-mask");
run("Set Scale...", "distance=0 known=0 pixel=1 unit=pixel global");
run("Analyze Particles...", "size=150-1500 circularity=0.10-1.00 show=Masks clear include add stack");
rename("mask");
run("Clear Results");
selectWindow("mask");
close();
//selectWindow("nuclear-mask");
//close();
//----------------------------
print("Site ",site," Creating Mask",(getTime-start)/1000); 

//--Measure--
selectWindow("original");
run("Set Measurements...", "area mean standard min centroid perimeter shape median area_fraction stack add redirect=None decimal=3");
roiManager("Measure");
saveAs("Results", "C:\\tracking_20140916\\measure-mask-s"+site+".csv");
//----------------
close();
roiManager("Show None");
run("Clear Results");
print("Site ",site," Measuring Mask",(getTime-start)/1000);
 

//--Measure 2nd YFP-FOXO channel---

run("Image Sequence...","open=Z:\\sorger\\data\\NIC\\Jerry\\20140916FOXO-EGFIGF\\20140916-FOXO-EGFIGF_w1YFP_s45_t197.TIF number=29161 starting=1 increment=1 scale=100 file=w1YFP_s"+site+"_ sort");

//run("Image Sequence...", "open=Z:\\sorger\\data\\Operetta\\Jerry\\20140826-FOXO-live-EGF-insulin[1371]\\20140826-live-2[5948]\\2014-08-26T154535-0400[7197]\\003010-1-034001001.tif starting=2 increment=2 file="+rowarray[row]+colarray[col]+"-1- sort");
run("Subtract Background...", "rolling=50 stack");

roiManager("Show None"); 
counts=roiManager("count");
for(i=0; i<counts; i++) {
    roiManager("Select", i);
    run("Measure");
    run("Enlarge...", "enlarge=3");
    run("Measure");
}
saveAs("Results", "C:\\tracking_20140916\\measure-cytosol-s"+site+".csv");
run("Clear Results");

setBatchMode(false); 

print("Site ",site," Final measure",(getTime-start)/1000); 

//-----------------------------------

}


