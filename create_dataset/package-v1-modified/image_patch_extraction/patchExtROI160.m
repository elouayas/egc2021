function [cell locRow locCol]=patchExtROI160(normImg,ctrPointMask,magnification,rowStartBlock, colStartBlock)
% input: normImg-- normalized image area, ctrPointMask: centroid point of
% each extracted block within image area, (rowStartBlock,colStartBlock)--the start
% point of the image block within the image area normImg.

% output: cell: image patch, (locRow,locCol)-location of each image patch


[cellCtrdRow cellCtrdCol]=find(ctrPointMask==1);
num=length(cellCtrdRow);

[row col]=size(ctrPointMask);
[trow tcol tvol]=size(normImg);


if magnification==40  % extract 160 * 160 and shrink to 80 *80
    if trow<160 || tcol<160
        cell=[];
        locRow=[];
        locCol=[];
    else
   
    for i=1:num  % num: number of image patch at this block
        rowStatImg=cellCtrdRow(i)+rowStartBlock-1-79;
        rowEndImg=cellCtrdRow(i)+rowStartBlock-1+80;
        if rowStatImg<1
            rowEndImg=rowEndImg+abs(rowStatImg)+1;
            rowStatImg=1;
        elseif rowEndImg>trow
            rowStatImg=rowStatImg-(rowEndImg-trow);
            rowEndImg=trow;
        end
        
        colStatImg=cellCtrdCol(i)+colStartBlock-1-79;
        colEndImg=cellCtrdCol(i)+colStartBlock-1+80;
        if colStatImg<1
            colEndImg=colEndImg+abs(colStatImg)+1;
            colStatImg=1;
        elseif colEndImg>tcol
            colStatImg=colStatImg-(colEndImg-tcol);
            colEndImg=tcol;
        end
        patchTemp=normImg(rowStatImg:rowEndImg,colStatImg:colEndImg,1:3);
        patchTemp=imresize(patchTemp,[80 80]);  %shrunk to 80 * 80
        temp{i}=patchTemp;
    end
    
    if num~=0
        cell=temp;
        locRow=cellCtrdRow;
        locCol=cellCtrdCol;
    else
        cell=[];
        locRow=[];
        locCol=[];
    end
    end
    
elseif magnification==20 % extract 80 * 80 directly
    if trow<80 || tcol<80
        cell=[];
        locRow=[];
        locCol=[];
    else
        for i=1:num
           rowStatImg=cellCtrdRow(i)+rowStartBlock-1-39;
           rowEndImg=cellCtrdRow(i)+rowStartBlock-1+40;
           if rowStatImg<1
             rowEndImg=rowEndImg+abs(rowStatImg)+1;
             rowStatImg=1;
           elseif rowEndImg>trow
             rowStatImg=rowStatImg-(rowEndImg-trow);
             rowEndImg=trow;
           end
           
           colStatImg=cellCtrdCol(i)+colStartBlock-1-39;
           colEndImg=cellCtrdCol(i)+colStartBlock-1+40;
           if colStatImg<1
              colEndImg=colEndImg+abs(colStatImg)+1;
              colStatImg=1;
           elseif colEndImg>tcol
              colStatImg=colStatImg-(colEndImg-tcol);
              colEndImg=tcol;
           end
           patchTemp=normImg(rowStatImg:rowEndImg,colStatImg:colEndImg,1:3);
           temp{i}=patchTemp;
        end
        
        if num~=0
          cell=temp;
          locRow=cellCtrdRow;
          locCol=cellCtrdCol;
        else
          cell=[];
          locRow=[];
          locCol=[];
        end
    end
end