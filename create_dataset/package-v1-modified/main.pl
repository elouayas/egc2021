#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use Cwd 'abs_path';

##########PART 0, setting environment#############
my $DataDir = "*Please define your data dir with -dd*";
my $roi = "ROI";
my $partition_n = 1000;
my $Version = "";
my $Help = "";
my $citation = "";
my $i = 1;
GetOptions ("dd=s" => \$DataDir, "roi=s" => \$roi, "partition=i" => \$partition_n, "version!" => \$Version, "help!" => \$Help, "citation!" => \$citation);

#***check option***s
if ($Version){
	print "This is Pathologicla Image Analysis Package, version: 1.0\n\n";
	print "Copyright University of Texas Southwestern Medical Center\n\n";
	print "You can get updated version at our Home Page: http://.....\n";
	exit;
}
if ($Help){
	print "Usage: perl main.pl [options]\n\n";
	print "Options (abbreviation to uniqueness allowed):\n";
	print "    -dd        Setting DataDir\n";
	print "    -roi       Setting Region of Interest\n";
	print "    -partition Setting size of Deeplearning patch bunch\n";
	print "    -version   Report Version\n";
	print "    -help      Display this help message\n";
	print "    -citation  Display Citation format\n";
	print "\n                [Pathological Image Analysis Package v1.0]\n";
	exit;
}
if ($citation){
	print "If you used our program, please cite us as ********\n";
	print "Contact us at Tao.Wang\@UTSouthwestern.edu\n";
	exit;
}

print "Welcome to use our pathological image processing program!\n";
#****check if Dirs exist****
my $WorkingDir = abs_path($0);
$WorkingDir =~ s/\/\w+.pl//;
$DataDir =~ s/\/$//;
if (!-e $WorkingDir){
	die "WorkingDir $WorkingDir cannot be found!";
}
if (!-e $DataDir){
	die "DataDir $DataDir cannot be found!";
}
#****making required dir.****
my @directory = qw(patch_extraction_results deeplearning_results cell_Info cell_border extracted_features);
for my $i (@directory){
#	unless(-e "$DataDir/$i" or mkdir("$DataDir/$i", 0755)){
#		die "Unable to create $DataDir/$i\n";
#	}
#	print("Making directory: $DataDir/$i...\n");
	if (-e "$DataDir/$i"){
		print("Existed directory: $DataDir/$i...\n");
	}else{
		mkdir("$DataDir/$i", 0755) or die "Unable to create $DataDir/$i\n";
	}
}
#get all original slide names
my $slidefolder = $DataDir;
opendir my $dir, $slidefolder or die "Cannot open directory: $slidefolder!";
my @slides = grep {/\.svs$/} readdir $dir;
if (scalar @slides < 1){
	die "No image slide is detected!";
}
print("The following image files are detected:\n");
foreach my $slide (@slides){
	print("$slide\n");
	#print(substr $slide, 0, -4);
}
closedir $dir;

###########PART 1, image patch extraction###########
print("********************\n");
print("Start 10 sampling regions extraction and cell segmentation\n");
print("********************\n");
my $codefolder = "$WorkingDir/image_patch_extraction";
my $savefolder = "$DataDir/patch_extraction_results";
system("/usr/local/MATLAB/R2018a/bin/matlab -nodisplay -nodesktop -r \"cd('".$codefolder."'); cell_extract_bordure_local('".$slidefolder."', '".$codefolder."', '".$savefolder."', '".$roi."'); quit;\" >> ".$DataDir."/log.txt");

###########PART 2, deep learning patch classification############
foreach my $slide (@slides){
	my $cellinfooutput = "$DataDir/cell_Info";
	my $slidename = substr $slide, 0, -4;
	
	print("********************\n");
	print("Deeplearning on the $i th sampling region of slide $slidename\n");
	print("********************\n");
	my $deepinput = "$savefolder/ImagePatchInfo/".$slidename."_".$i.
"_".$roi."_Image_Patch";
	my $deepoutput = "$DataDir/deeplearning_results/".$slidename."_".$i;
	if (!-e $deepoutput){
        	mkdir("$deepoutput", 0755) or die "Unable to create $DataDir/$i\n";
	}
	system("Rscript $WorkingDir/deeplearning/deeplearning.R \"$deepinput\" \"$deepoutput\" \"$partition_n\"");
	#summarize deeplearning output
	print("Summerizing deeplearning output...\n");
	system("Rscript $WorkingDir/deeplearning_summarize/txtToCSV.R \"".$slidename."_".$i."\" \"$DataDir\"");
}


##########PART 3, features extraction##############
print("********************\n");
print("Starting extracting features\n");
print("********************\n");
system("/usr/local/MATLAB/R2018a/bin/matlab -nodisplay -nodesktop -r \"cd(['".$WorkingDir."', '/image_border/']); main('"."$DataDir"."'); quit;\" >> ".$DataDir."/log.txt");
print("Well done!\n");

