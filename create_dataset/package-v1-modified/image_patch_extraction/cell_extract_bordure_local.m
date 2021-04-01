function cell_extract(slidefolder,codefolder,savefolder,roi)

%slidefolder='~/shared_folder_ubuntu/cellules/40/TCGA-38-4625-01Z-00-DX1';
%codefolder='~/Bureau/package-v1/image_patch_extraction';
%savefolder ='~/Bureau/trash';
%roi='internet';
cd(slidefolder);
addpath(genpath(codefolder));
curdate=savefolder;

archivedir=dir;
archivenames={archivedir.name};
logisdir=[archivedir.isdir]==1;
isnotdir=~logisdir;
archivefilename=archivenames(isnotdir);
selectfile=regexp(archivefilename,'.xml');
archivefilename=archivefilename(:,~cellfun(@isempty,selectfile));

for ifileNum=1:length(archivefilename)   
    filename=char(archivefilename(ifileNum));
    baseName=filename(1:end-4)
    sprintf('The %dth image',ifileNum)
    
    currsvs=[baseName '.svs'];
    rhead=bfGetReader(currsvs);
    omeMeta=rhead.getMetadataStore();
    magnification = omeMeta.getObjectiveNominalMagnification(0,0);
    magnification = num2str(double(magnification));
    magnification = str2num(magnification)
    xmlFile=[baseName '.xml'];

        interestedRegion=roi;

        %============ Step 1: extract the interested region =============
        regionBoundary=extract_boundary_points(xmlFile,interestedRegion);  
        interestedRegion(find(isspace(interestedRegion))) = '_';   
        num=length(regionBoundary);     
        
        %======Step 2: segment the cells within interested region========
        if num~=0
           temp=regionBoundary{1};       
           xPos=round(temp(:,1));
           yPos=round(temp(:,2));
                
           minRowTmp=min(yPos);            
           minColTmp=min(xPos);
           maxRowTmp=max(yPos);            
           maxColTmp=max(xPos);
           finRowTmp=floor((maxRowTmp+minRowTmp)/2);            
           finColTmp=floor((maxColTmp+minColTmp)/2);
                        
                
           xPosTmp=ceil(temp(:,1)/10);     
           yPosTmp=ceil(temp(:,2)/10);    
                
           diskSize=5;
           objMask=regMask(xPosTmp,yPosTmp,diskSize);
             
           if magnification==40
              se=strel('disk',120);
              nRow=5000;
              nCol=5000;
           else
              se=strel('disk',50);
              nRow=3000;
              nCol=3000;
           end
                
           if size(objMask,1)<nRow/10 || size(objMask,2)<nCol/10
  	       sprintf('major2')
              continue;
           end
                
           objMaskCover=zeros(size(objMask,1)+2,size(objMask,2)+2);
           objMaskCover(2:end-1,2:end-1)=objMask;
           objMaskTmp2=imerode(objMaskCover,se);
           objMaskTmp=objMaskTmp2(2:end-1,2:end-1);
                
           if sum(sum(objMaskTmp))<10   
               objMaskTmp=objMask;
           end
              
           objMaskCenter=zeros(size(objMask,1),size(objMask,2));
           objMaskCenter((nRow/20)+1:end-nRow/20,(nCol/20)+1:end-nCol/20)=1; 
           objMaskTmp2=objMaskTmp & objMaskCenter;
              
           % random select 10 points wihtin the region.
           [rowIndex,colIndex]=find(objMaskTmp2==1);
           %randomP=randperm(length(rowIndex));
             
           isubTmp=1;
       %for iTen=1:10                       
           %cRowShk=rowIndex(randomP(iTen));  
           %cColShk=colIndex(randomP(iTen));  
           cRowShk=floor((max(rowIndex)+min(rowIndex))/2);
           cColShk=floor((max(colIndex)+min(rowIndex))/2);
           cRow=finRowTmp;  % amplify 10 times at y axis
           cCol=finColTmp;  % amplify 10 times at x axis
           minRow=cRow-(nRow/2)+1;
           maxRow=cRow+(nRow/2);
           minCol=cCol-(nCol/2)+1;
           maxCol=cCol+(nCol/2);

           cropRegionImg(:,:,1)=bfGetPlane(rhead,1,minCol,minRow,nCol,nRow);
           cropRegionImg(:,:,2)=bfGetPlane(rhead,2,minCol,minRow,nCol,nRow);
           cropRegionImg(:,:,3)=bfGetPlane(rhead,3,minCol,minRow,nCol,nRow);

           %=========Corresponding masks are defined==================
           minRowSub=cRowShk-(nRow/20)+1;
	   if minRowSub<1
		   minRowSub=1;
	   end
           maxRowSub=cRowShk+(nRow/20);
           minColSub=cColShk-(nCol/20)+1;
	   if minColSub<1
		   minColSub=1;
	   end
           maxColSub=cColShk+(nCol/20);
	   if maxColSub==0 || maxRowSub==0
		   sprintf('3 major')
		   continue
	   end
	   try
		   objMaskSubRegion=objMask(minRowSub:maxRowSub,minColSub:maxColSub);
	   catch
		   sprintf('1 major')
		   continue
	   end 
           
           NormRH=mat2gray(cropRegionImg);   

           totalRow=maxRow-minRow+1;
           totalCol=maxCol-minCol+1;

           tileRow=500;       
           tileCol=500;

           rowNum=ceil(totalRow/tileRow);
           colNum=ceil(totalCol/tileCol);

           isub=1;
           index=1;
           excelIndex=1;
           saveInfo{1,1}='Image Name';    %image name
           saveInfo{1,2}='Patch Name';    %patch name
           saveInfo{1,3}='Location in X axis';  %location of patch center in x axis
           saveInfo{1,4}='Location in Y axis';  %location of pathc center in y axis

           for i=1:rowNum
               for j=1:colNum
                  rowStart=(i-1)*tileRow+1;
                  rowEnd=i*tileRow;
                  colStart=(j-1)*tileCol+1;
                  colEnd=j*tileCol;

                  imgTemp=cropRegionImg(rowStart:rowEnd,colStart:colEnd,1:end);
                  imgTempNorm=NormRH(rowStart:rowEnd,colStart:colEnd,1:end);         %normalized block

                  objCtrMask=segNuclei(imgTemp,magnification);     %segment the extracted block

                  rowStartBlock=ceil(rowStart/10);
                  rowEndBlock=round(rowEnd/10);
                  colStartBlock=ceil(colStart/10);
                  colEndBlock=round(colEnd/10);
                  rse = rowStartBlock:rowEndBlock;
                  cse = colStartBlock:colEndBlock;		  


                  objMaskTmp2=objMaskSubRegion(rse,cse);
                  objMaskTmp3=imresize(objMaskTmp2,[tileRow tileCol]);  %amplification  
                  objMaskTmp4=objMaskTmp3>0.10;
                  filterObjCtrMask=objCtrMask & objMaskTmp4;

                  if sum(sum(filterObjCtrMask))<20
                      sprintf('2')	
                      continue;
                  end


                 % ========= Save the marked Regions===================
                  baseName2=baseName;
                  baseName2(find(baseName2==' '))='_';

                  interestedRegionMarked=[curdate '/MarkedImageInfo/' baseName2 '_' num2str(isubTmp) '_' interestedRegion '_Marked_Image'];

                  if ~isdir(interestedRegionMarked)
                     mkdir(interestedRegionMarked);
                  end

                 %============== Extract the image patch =============
                  [cellPatch,locRow,locCol]=patchExtROI160(NormRH,filterObjCtrMask,magnification,rowStart,colStart);  

                  patchNum=length(cellPatch); 
                  if patchNum~=0
                      dirTumorRegion=[curdate '/ImagePatchInfo/' baseName2 '_' num2str(isubTmp) '_' interestedRegion '_Image_Patch'];

                      if ~isdir(dirTumorRegion)
                          mkdir(dirTumorRegion);
                      end

                      for iNum=1:patchNum
                        cellPatchTemp=cellPatch{iNum};
                        eval(['imwrite(cellPatchTemp,' '''' dirTumorRegion '/' interestedRegion '_' baseName2 '_' num2str(isubTmp) '_Patch_' num2str(index) '.tif' '''',');'])

                        saveInfo{excelIndex+1,1}=baseName2;
                        saveInfo{excelIndex+1,2}=[interestedRegion '_' baseName2 '_' num2str(isubTmp) '_Patch_' num2str(index) '.tif'];
                        saveInfo{excelIndex+1,3}=locCol(iNum)+(minCol-1)+(colStart-1);
                        saveInfo{excelIndex+1,4}=locRow(iNum)+(minRow-1)+(rowStart-1);

                        index=index+1;
                        excelIndex=excelIndex+1;
                      end
                  end
                  imgTemp=ctrDraw2(imgTemp,filterObjCtrMask,3);
                  eval(['imwrite(imgTemp,' '''' interestedRegionMarked '/' baseName2 '_Area_' num2str(isubTmp) '_subArea_' num2str(isub) '.tif' '''',');']) 
                  isub=isub+1;
                  clear imgTemp   
               end
           end

           csvfile=[curdate '/cellLocationInfo'];
           if ~isdir(csvfile)
              mkdir(csvfile);
           end
           fileName=[csvfile '/' baseName2 '_' num2str(isubTmp) '.csv'];
           cell2csv(fileName,saveInfo);

           clear cropRegionImg;
           clear saveInfo;
           isubTmp=isubTmp+1;     
       end         
        end   
    disp('Well done!');
end

