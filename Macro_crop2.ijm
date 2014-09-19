//Macro for crop tracking images
//Jerry Lin 20140910

pathfile=File.openDialog("Choose the file to Open:");

filestring=File.openAsString(pathfile);

rows=split(filestring, "\n");

//s=newArray(rows.length);
//x=newArray(rows.length);
//y=newArray(rows.length);

for(i=0; i<rows.length; i++){
columns=split(rows[i],"\t");
s=parseInt(columns[0]);
x=parseInt(columns[1]);
y=parseInt(columns[2]);
t=parseFloat(columns[3]);
print(s,x,y,t);

selectWindow("image1");
setSlice(s);
x=x-50;
y=y-50;

makeRectangle(x,y,100,100);
run("Duplicate...","title=frame"+i);
x=4;
y=14;
setForegroundColor(0, 0, 0);
setFont("SansSerif", 14, " antialiased");
setColor("yellow");
drawString("NC="+t,x,y);
x=4;
y=98;
drawString("Frame="+(i+1),x,y);
}

run("Images to Stack", "method=[Copy (center)] name=crop1 title=frame use");
