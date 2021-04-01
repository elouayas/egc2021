function outImg=ctrDraw2(inImg,mask,type)
R=inImg(:,:,1);
G=inImg(:,:,2);
B=inImg(:,:,3);

cross=zeros(7,7);
cross(:,4)=1;
cross(4,:)=1;

square=strel('square',3);
diamond=strel('diamond',2);
bry=strel('disk',1);

if type==1    
    index=imdilate(mask,cross);
    index=logical(index);
    R(index)=255;
    G(index)=0;
    B(index)=0;
    
elseif type==2
    index=imdilate(mask,square);
    index=logical(index);
    R(index)=255;
    G(index)=255;
    B(index)=0;
    
elseif type==3
    index=imdilate(mask,cross);
    index=logical(index);
    R(index)=0;
    G(index)=255;
    B(index)=0;
    
else
    index=mask;
    index=logical(index);
    R(index)=255;
    G(index)=0;
    B(index)=0;
end

clear inImg

outImg(:,:,1)=R;
outImg(:,:,2)=G;
outImg(:,:,3)=B;