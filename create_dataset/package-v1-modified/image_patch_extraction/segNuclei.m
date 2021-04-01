function objCtrMask=segNucleiLSRF(img,magnification)

% input: the image of interested region (img)
% output: the centroid point of cell
% Faliu.Yi@UTSouthwestern.edu

%========================0: read image=====================================
%==========================================================================
sampleRGB=img;
[height width channel]=size(sampleRGB);

%========================1: color deconvolution============================
%==========================================================================
SourceImage=sampleRGB;
H_deinterlace = [0 1 0; 0 2 0; 0 1 0] ./4;

sample_deinterlace = zeros(height, width, channel);
for k=1:channel
    sample_deinterlace(:,:,k) = filter2(H_deinterlace,double(sampleRGB(:,:,k)),'same');
end

%=== Convert RGB intensity to optical density (absorbance)
sampleRGB_OD = -log((sample_deinterlace+1)./256);

%% Construct color deconvolution matrix
H2 = ones(10,10) ./ 100;
sampleRGB_OD_Blur = zeros(height,width,channel);

for k=1:channel
    sampleRGB_OD_Blur(:,:,k) = filter2(H2,sampleRGB_OD(:,:,k),'same');
end

% Standard values from literature
He = [0.550 0.758 0.351]';
Eo = [0.398 0.634 0.600]';
Bg = [0.754 0.077 0.652]';

% Create Deconvolution matrix
M = [He/norm(He) Eo/norm(Eo) Bg/norm(Bg)];
D = inv(M);

% Apply Color Deconvolution
sampleHEB_OD = zeros(height, width, channel);
for i=1:height
    for j=1:width
        RGB = reshape(sampleRGB_OD(i,j,:),channel,1);
        HEB = D * RGB;

      	sampleHEB_OD(i,j,1) = HEB(1);
       	sampleHEB_OD(i,j,2) = HEB(2);
       	sampleHEB_OD(i,j,3) = HEB(3);
    end
end

% Extract tumor cells that are stained by hematoxylin only
hematoxylin = sampleHEB_OD(:,:,1);
hematoxylin = (hematoxylin - 0.05) .* (hematoxylin > 0.05);
hematoxylin = hematoxylin ./ max(max(hematoxylin));
h=fspecial('sobel');
hematoxylin_grad=sqrt(imfilter(hematoxylin,h,'replicate').^2+imfilter(hematoxylin,h','replicate').^2);


%=================2: morphological operstions==============================
%==========================================================================
hImg=hematoxylin;

se1=strel('disk',4);              %opening by reconstruction: remove small object
hImgEro=imerode(hImg,se1);
hImgEroRec=imreconstruct(hImgEro,hImg);

se2=strel('disk',7);
hImgDia=imdilate(hImgEroRec,se2);  %closing by reconstruction: remove small black spot
hImgDiaRec=imreconstruct(imcomplement(hImgDia),imcomplement(hImgEroRec));
hImgDiaRec=imcomplement(hImgDiaRec);
hImgFill=imfill(hImgDiaRec,'holes');


%=======================3:Thresholding=====================================
%==========================================================================
T1=graythresh(hImgFill);
hImgThres=(hImgFill>T1);
hImgThresFill=imfill(hImgThres,'holes');
exMaskAdd=imcomplement(hImgThresFill);  % used as additional external markers
se22=strel('disk',3);
exMaskAdd=imerode(exMaskAdd,se22);
se3=strel('disk',2);
hImgClose=imopen(hImgThresFill,se3);   

hImgAreaOpen=bwareaopen(hImgClose,30); 



distInfo=bwdist(imcomplement(hImgAreaOpen));
inMark=distInfo>2;
inMark=bwareaopen(inMark,5);          
dist2=bwdist(inMark);
inSeg=watershed(dist2);
exMark=inSeg==0;                      


se4=strel('disk',1);                  
inMark=imerode(inMark,se4);
markers=inMark | exMark | exMaskAdd;              

%===========================4: Level Set==================================
%=========================================================================

Img = double(hImgFill);
[row,col] = size(Img);

phi = ones(row,col);
phi(exMark|exMaskAdd)=-1;
u=phi;

sigma = 1;
G = fspecial('gaussian', 11, sigma);   %5

delt = 1;
Iter = 30;     
mu=60;

for n = 1:Iter
    [ux, uy] = gradient(u);
   
    c1 = sum(sum(Img.*(u<0)))/(sum(sum(u<0)));
    c2 = sum(sum(Img.*(u>=0)))/(sum(sum(u>=0)));
    
    spf = Img - (c1 + c2)/2;
    spf = spf/(max(abs(spf(:))));
    
    u = u + delt*(mu*spf.*sqrt(ux.^2 + uy.^2));
    u = (u >= 0) - ( u< 0);
    u = conv2(u, G, 'same');
end

temp=im2bw(u,0.8);
temp=bwareaopen(temp,26);  %remove small region



imgMask=temp;   % from level set
imgMaskDist=bwdist(imcomplement(imgMask));
imgDistMax=imregionalmax(imgMaskDist);

if magnification==40
se5=strel('disk',8);
imgDistMax=imdilate(imgDistMax,se5);

elseif magnification==20
se5=strel('disk',4);
imgDistMax=imdilate(imgDistMax,se5);
end



[label1,num1]=bwlabel(imgDistMax);
prop1=regionprops(label1,'Centroid');
ctrPos1=cat(1,prop1.Centroid);

xp=ctrPos1(:,1);
yp=ctrPos1(:,2);

ctrRow=round(yp);
ctrCol=round(xp);

objCtrMaskTemp=zeros(height,width);
indexTemp=(ctrCol-1)*height+ctrRow;
objCtrMaskTemp(indexTemp)=1;

objCtrMask=objCtrMaskTemp;        




%=========================separate overlapping region ==============================
if magnification==40
   S1=fastradial(hImgFill,7,2,0.01);
   S1=imcomplement(S1);

   regMin1=imhmin(S1,0.01);
   regMin1=imregionalmin(regMin1);
   se=strel('disk',1);
   regMin1=imdilate(regMin1,se);
%    figure
%    imshow(regMin1)
%    title('radial symmetry with size 8')

   S2=fastradial(hImgFill,9,2,0.01);
   S2=imcomplement(S2);

   regMin2=imhmin(S2,0.01);
   regMin2=imregionalmin(regMin2);
   se=strel('disk',1);
   regMin2=imdilate(regMin2,se);
%    figure
%    imshow(regMin2)
%    title('radial symmetry with size 11')

   S3=fastradial(hImgFill,13,2,0.01);
   S3=imcomplement(S3);

   regMin3=imhmin(S3,0.01);
   regMin3=imregionalmin(regMin3);
   se=strel('disk',1);
   regMin3=imdilate(regMin3,se);
%    figure
%    imshow(regMin3)
%    title('radial symmetry with size 15')

   S4=fastradial(hImgFill,17,2,0.01);
   S4=imcomplement(S4);

   regMin4=imhmin(S4,0.01);
   regMin4=imregionalmin(regMin4);
   se=strel('disk',1);
   regMin4=imdilate(regMin4,se);
%    figure
%    imshow(regMin4)
%    title('radial symmetry with size 19')

   S5=fastradial(hImgFill,21,2,0.01);
   S5=imcomplement(S5);

   regMin5=imhmin(S5,0.01);
   regMin5=imregionalmin(regMin5);
   se=strel('disk',1);
   regMin5=imdilate(regMin5,se);
%    figure
%    imshow(regMin5)
%    title('radial symmetry with size 25')

   for iRF=1:4
      eval(['rfMarkerTemp=regMin' num2str(iRF) '& regMin' num2str(iRF+1) ';']);
      eval(['im' num2str(iRF) '=imreconstruct(rfMarkerTemp,regMin' num2str(iRF) ');']);
      eval(['im' num2str(iRF+1) '=imreconstruct(rfMarkerTemp,regMin' num2str(iRF+1) ');']);
      eval(['regMin' num2str(iRF+1) '=rfMarkerTemp+(regMin' num2str(iRF) '-im' num2str(iRF) ')+(regMin' num2str(iRF+1) '-im' num2str(iRF+1) ');'])
      eval(['regMin' num2str(iRF+1) '=logical(regMin' num2str(iRF+1) ');'])
   end
   rFinMarker=regMin5;
%    imgPerm=bwperim(rFinMarker);
%    overlay2=imoverlay(sampleRGB,imgPerm,[.3 1 .3]);
%    figure
%    imshow(overlay2)
%    title('internal markers from RF')

elseif magnification==20
   S1=fastradial(hImgFill,4,2,0.01);
   S1=imcomplement(S1);

   regMin1=imhmin(S1,0.01);
   regMin1=imregionalmin(regMin1);
   se=strel('disk',1);
   regMin1=imdilate(regMin1,se);
%    figure
%    imshow(regMin1)
%    title('radial symmetry with size 4')

   S2=fastradial(hImgFill,6,2,0.01);
   S2=imcomplement(S2);

   regMin2=imhmin(S2,0.01);
   regMin2=imregionalmin(regMin2);
   se=strel('disk',1);
   regMin2=imdilate(regMin2,se);
%    figure
%    imshow(regMin2)
%    title('radial symmetry with size 6')

   S3=fastradial(hImgFill,8,2,0.01);
   S3=imcomplement(S3);

   regMin3=imhmin(S3,0.01);
   regMin3=imregionalmin(regMin3);
   se=strel('disk',1);
   regMin3=imdilate(regMin3,se);
%    figure
%    imshow(regMin3)
%    title('radial symmetry with size 8')

   S4=fastradial(hImgFill,10,2,0.01);
   S4=imcomplement(S4);

   regMin4=imhmin(S4,0.01);
   regMin4=imregionalmin(regMin4);
   se=strel('disk',1);
   regMin4=imdilate(regMin4,se);
%    figure
%    imshow(regMin4)
%    title('radial symmetry with size 10')

   S5=fastradial(hImgFill,11,2,0.01);
   S5=imcomplement(S5);

   regMin5=imhmin(S5,0.01);
   regMin5=imregionalmin(regMin5);
   se=strel('disk',1);
   regMin5=imdilate(regMin5,se);
%    figure
%    imshow(regMin5)
%    title('radial symmetry with size 11')

   for iRF=1:4
      eval(['rfMarkerTemp=regMin' num2str(iRF) '& regMin' num2str(iRF+1) ';']);
      eval(['im' num2str(iRF) '=imreconstruct(rfMarkerTemp,regMin' num2str(iRF) ');']);
      eval(['im' num2str(iRF+1) '=imreconstruct(rfMarkerTemp,regMin' num2str(iRF+1) ');']);
      eval(['regMin' num2str(iRF+1) '=rfMarkerTemp+(regMin' num2str(iRF) '-im' num2str(iRF) ')+(regMin' num2str(iRF+1) '-im' num2str(iRF+1) ');'])
      eval(['regMin' num2str(iRF+1) '=logical(regMin' num2str(iRF+1) ');'])
   end
   rFinMarker=regMin5;
end

if magnification==40
   se1=strel('disk',5);
   imgMaskRF=imdilate(rFinMarker,se1);
elseif magnification==20
    se1=strel('disk',3);
    imgMaskRF=imdilate(rFinMarker,se1);
end

[label1 num1]=bwlabel(imgMaskRF);
prop1=regionprops(label1,'Centroid');
ctrPos1=cat(1,prop1.Centroid);

xp=ctrPos1(:,1);
yp=ctrPos1(:,2);

ctrRow=round(yp);
ctrCol=round(xp);

objCtrMaskTemp=zeros(height,width);
indexTemp=(ctrCol-1)*height+ctrRow;
objCtrMaskTemp(indexTemp)=1;

objCtrMask2=objCtrMaskTemp;              % dot representation from RF

[ypf xpf]=find(objCtrMask2==1);


[label1,num1]=bwlabel(imgMask);
cellProp1=regionprops(label1,'area');
cellArea=[cellProp1.Area];
[aValue aIndex]=sort(cellArea,'descend');

bigAreaMask=logical(zeros(height,width));

if length(aIndex)~=0
   if magnification==40
       bigIndex=aValue>10000;
   elseif magnification==20
       bigIndex=aValue>2500;
   end
   
   bigNum=sum(bigIndex);
   
   if bigNum~=0
%       bigAreaMask=logical(zeros(height,width));
      for ibigNum=1:bigNum
          bigAreaIndex=aIndex(ibigNum);
          bigAreaMaskTemp=label1==bigAreaIndex;
          bigAreaMask=bigAreaMask | bigAreaMaskTemp;
      end
   end   
end

objCtrMask3=objCtrMask2 & bigAreaMask;


bigAreaMaskComp=imcomplement(bigAreaMask);
objCtrMask4=bigAreaMaskComp & objCtrMask;
objCtrMask=objCtrMask3 | objCtrMask4;
if magnification==40
    se=strel('disk',5);
    objCtrMask=imdilate(objCtrMask,se);
elseif magnification==20
    se=strel('disk',3);
    objCtrMask=imdilate(objCtrMask,se);
end

if sum(sum(objCtrMask))~=0

[label1,num1]=bwlabel(objCtrMask);
prop1=regionprops(label1,'Centroid');
ctrPos1=cat(1,prop1.Centroid);

xp=ctrPos1(:,1);
yp=ctrPos1(:,2);

ctrRow=round(yp);
ctrCol=round(xp);

objCtrMaskTemp=zeros(height,width);
indexTemp=(ctrCol-1)*height+ctrRow;
objCtrMaskTemp(indexTemp)=1;

objCtrMask=objCtrMaskTemp;  
end

% figure
% imshow(objCtrMask)
% title('final dot representation')

% [yp xp]=find(objCtrMask==1);
% figure
% colormap(gray(256)), imagesc(sampleRGB);
% axis off
% hold on
% plot(xp,yp,'.r')
% hold off




