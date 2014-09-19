//script for creating mask and measure for live imaging 
//Jerry 2014/09/12
// for Operetta imaging
// optimize segmenation for operetta imaging

rowarray = newArray("000","001","002","003","004","005","006","007","008");
colarray =newArray("000","001","002","003","004","005","006","007","008","009","010","011","012");

for (row=2;row<3;row++)
{
for (col=2;col<3;col++)
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
//run("Image Sequence...", "open=Z:\\sorger\\data\\NIC\\Jerry\\20140822-NIC-foxo\\FOXO-inhibitors_w1RFP_s1_t1.TIF number=29161 starting=1 increment=1 scale=100 file=w1RFP_s"+site+"_ sort");
run("Image Sequence...", "open=G:\\20140704_Live_foxo[1254]\\0704-1[5173]\\2014-07-04T160249-0400[6209]\\002002-1-004001001.tif starting=2 increment=2 file="+rowarray[row]+colarray[col]+"-1- sort");

run("Subtract Background...", "rolling=25 stack");
setOption("BlackBackground", false);
rename("original");
run("Duplicate...", "title=mask duplicate range=1-300");
//---------------------------------------------------
print("Well ",row,"-",col," Building Stacks",(getTime-start)/1000); 


//--processing to binary mask--
selectWindow("mask");

run("Gaussian Blur...", "sigma=1.5 stack");
run("Enhance Contrast...", "saturated=0.5");
run("Make Binary", "method=Li background=Dark calculate");
run("Fill Holes", "stack");
run("Watershed", "stack");

rename("nuclear-mask");
//------------------------------------
print("Well",row,"-",col," Convert to Mask",(getTime-start)/1000); 

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
print("Well",row,"-",col," Creating Mask",(getTime-start)/1000); 

//--Measure--
selectWindow("original");
run("Set Measurements...", "area mean standard min centroid perimeter shape median area_fraction stack add redirect=None decimal=3");
roiManager("Measure");
saveAs("Results", "C:\\tracking_20140820-2\\measure-mask-"+rowarray[row]+colarray[col]+".csv");
//----------------
close();
roiManager("Show None");
run("Clear Results");
print("Well",row,"-",col," Measuring Mask",(getTime-start)/1000);
 

//--Measure 2nd YFP-FOXO channel---
//run("Image Sequence...","open=Z:\\sorger\\data\\NIC\\Jerry\\20140822-NIC-foxo\\FOXO-inhibitors_w2YFP_s1_t1.TIF number=29161 starting=1 increment=1 scale=100 file=w2YFP_s"+site+"_ sort");
run("Image Sequence...", "open=G:\\20140704_Live_foxo[1254]\\0704-1[5173]\\2014-07-04T160249-0400[6209]\\002002-1-004001001.tif starting=1 increment=2 file="+rowarray[row]+colarray[col]+"-1- sort");
run("Subtract Background...", "rolling=50 stack");

roiManager("Show None"); 
counts=roiManager("count");
for(i=0; i<counts; i++) {
    roiManager("Select", i);
    run("Measure");
    run("Enlarge...", "enlarge=3");
    run("Measure");
}
saveAs("Results", "C:\\tracking_20140820-2\\measure-cytosol-"+rowarray[row]+colarray[col]+".csv");
run("Clear Results");

setBatchMode(false); 

print("Well",row,"-",col," Final measure",(getTime-start)/1000); 

//-----------------------------------

}

}
