README
by Tao Wang, Faliu Yi & Shidan Wang, University of Texas Southwestern Medical Center, 2.1.2017
Thanks for any suggestions or bug report to help us improve our package!
*********
Current version only works under unix system. 
Anytime you want to stop this program, just press Ctrl+Z.
*********
==================== DEPENDENCIES ====================

Matlab:
1) Install matlab. (tested environment: R2015a and R2016b)
2) Export matlab path with something like: $ export PATH=$PATH:/Applications/MATLAB_R2015a.app/bin/;

R:
1) Needs to have Rscript utility
2) Install R with the following packages: kimisc, magic, abind, grid, tiff, xlsx;

Data:
Create a file folder, and put your image (.svs file) and region of interest (ROI) definition (.xml file) there. Each pair of svs file and xml file should share the same filename. All of your results will be in that folder.

==================== RUNTIME OPTIONS ====================

Usage: perl main.pl [OPTIONS]

***abbreviation to uniqueness allowed***
-dd:       Setting DataDir, where you put your image files *Required*
-roi:      Setting Region of Interest
-partition:Setting size of Deeplearning bunch of image patches; recommanded 1000
-version:  Report Version infomation
-help:     Display help message
-citation: Display Citation format

==================== OUTPUT ====================

-cell_border:	
	Images of predicted cell type and region border for sampling region 
		*Tumor cells: green
                *Lymphocytes: blue
                *Stroma cells: red

-cell_Info:	
	Cell center locations and predicted types

-deeplearning_results:
	Txt files with predicted probabilities of three cell types

-extracted_features:
	Perimeter and size infomation for each sampling region

-patch_extraction_results:
    -cellLocationInfo: 
	Cell center coordinates, one file for each sampling region
    -ImagePatchInfo: 
	All image patches, one folder for each sampling region
    -MarkedImageInfo:
	Original sample regions with detected cell center marked

====================  CONTACT/CITATION/... ==================== 

Email: tao.wang@utsouthwestern.edu
If you ever used XXX in your publication, please cite: **********
Version: 1.0

