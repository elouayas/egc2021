clc
clear all
close all

slidefolder='/home/faliu/Chinese_dataset/slide_image';     % the directory including slide images
codefolder='/home/faliu/take_over/image_patch_extraction'; % the directory including matlab codes
savefolder='/home/faliu/20170203';                         % the folder used to save the extracted cells and location infos
roi='ROI';                                                 % the name of interested region manually labelled

cell_extract(slidefolder,codefolder,savefolder,roi);     
% this function returns three sub-folders under savefolder
% =========================================================================
% ================== 1: ImagePatchInfo=====================================
% === including extracted image pathes for each slide image================
% ===================2: CellLocationInfo===================================
% ====including location info of each extracted image patch================
% ===================3: MarkedImageInfo====================================
% ====including images with detected cell labelled ========================