function regionBoundary=extract_boundary_points(xmlFile, regionName)
% read the xmlFile and region, then return the boundary position of region named as regionName
% input: xmlFile- XML file which has the region labelled as regionName.
% output: regionBoundary are the boundary position along x and y axis. Each Cell is corrresponding one region named as regionName
% note: There may be more than one region named as regionName, the first column denotes the x value along x axis (Column in image) and the second column denotes the y value along y axis (Row in image).
% Example: regionBoundary=extract_boundary_points('Test1.xml','tumor cells');
% read Test1.xml file and get the boundary points of region labelled with 'tumor cells'
%=======================================================================================================================================
% This function is developped by Faliu Yi at UT Southwestern Medical Center in 2015.
%=======================================================================================================================================

xmlDoc=xmlread(xmlFile);   %read the .xml file.
RegionArray=xmlDoc.getElementsByTagName('Region');

indexNum=0;  % number of region named as regionName

for i=0:RegionArray.getLength-1           % find element node of the tumor region in the XML file
    thisItem=RegionArray.item(i);
    
    for ii=1:14  % find the index of attribute "Text" 
        RegionAttributesName=char(thisItem.getAttributes.item(ii).getNodeName);  % Note: the Attributes are ordered automatically
        index1=ii;
        if strcmp(RegionAttributesName,'Text')
           break;
        end
    end
    
    RegionAttributesValue=char(thisItem.getAttributes.item(index1).getValue);  % find the attribute named as regionName such as "tumor region"
    if strcmp(RegionAttributesValue,regionName)
        indexNum=indexNum+1;
        regionIndex(indexNum)=i;
    end
end


if indexNum~=0
      
regionNum=length(regionIndex);
for iii=1:regionNum
  xPos=0;
  yPos=0;
  index2=regionIndex(iii);
  vertexArray=RegionArray.item(index2).getElementsByTagName('Vertex');
  vertexNum=vertexArray.getLength-1;

  for i=0:vertexNum
      thisItem=vertexArray.item(i);
      xPos(i+1)=str2num(thisItem.getAttributes.item(0).getValue);  % value in x axis should be the column value in image
      yPos(i+1)=str2num(thisItem.getAttributes.item(1).getValue);  % value in y axis should be the row value in image
  end
  
  regionBoundary{iii}=[xPos' yPos'];
  
end

else
  regionBoundary='';
end
