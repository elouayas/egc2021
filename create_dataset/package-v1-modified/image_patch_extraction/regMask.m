function mask=regMask(xPos,yPos,diskSize)

xRow=yPos;
yCol=xPos;

minRow=min(xRow);
maxRow=max(xRow);
minCol=min(yCol);
maxCol=max(yCol);

rowSize=maxRow-minRow+1;
colSize=maxCol-minCol+1;

imgMarker=zeros(rowSize,colSize);
rowMarker=xRow-minRow+1;
colMarker=yCol-minCol+1;
indexTemp=(colMarker-1)*rowSize+rowMarker;
imgMarker(indexTemp)=1;

m1=imgMarker;
disk1=strel('disk',diskSize);
imgMarker=imdilate(imgMarker,disk1);

m2=imgMarker;
maskTemp=imfill(imgMarker,'holes');

subTemp=maskTemp-imgMarker;
checkVal=sum(sum(subTemp));

if checkVal>1200 || sum(sum(m1))==sum(sum(m2))
    maskTemp=imerode(maskTemp,disk1);
    mask=maskTemp;
else
    diskSize=diskSize+11;
    mask=regMask(xPos,yPos,diskSize);
end
end

    
