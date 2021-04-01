% FRST Test
clc
clear all
close all

img=imread('g.jpg');
S=FRST(img,3,0.12,2,1);

figure
imshow(img)
title('original image')

figure
imshow(S)
title('After FRST')