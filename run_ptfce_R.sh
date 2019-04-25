#!/bin/bash
#
# Author: Tamas Spisak
# Requires:
# - FSL (https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/)
# - R (https://www.r-project.org)
# - pTFCE (https://spisakt.github.io/pTFCE/)

# A POSIX variable
OPTIND=1         # Reset in case getopts has been used previously in the shell.

while getopts "hr:d:" opt; do
    case "$opt" in
    h)
        #show_help
        echo "This script runs pTFCE R-package from the command line"
        echo "It requires FSL and R with pTFCE installed."
        echo ""
        echo "Usage: run_ptfce.sh [options] <INPUT_Z_STAT_IMAGE> <INPUT_MASK_IMAGE>" 
        echo ""
        echo "Command line options:"
        echo "-r Use the 4D residual file for smoothness estimation (more accurate)."
        echo "-d Specify degrees of freedom for 4D residual based smoothness estimation (obligatory with -r)."
        echo ""
        echo "Examples:" 
        echo "run_ptfce.sh zstat.nii.gz mask.nii.gz"
        echo "run_ptfce.sh -r res4d.nii.gz -d 32 zstat.nii.gz mask.nii.gz"
        echo ""
        echo "Author: Tamas Spisak"
        echo "https://spisakt.github.io/pTFCE/"
        exit 0
        ;;
    r)  residual=$OPTARG
        ;;
    d)  dof=$OPTARG
        ;;
    esac
done

shift $(($OPTIND - 1))

if [ "$#" -ne 2 ]; then
    echo "Illegal number of parameters."
    echo "
    "
    echo "Usage: run_ptfce.sh [options] <INPUT_Z_STAT_IMAGE> <INPUT_MASK_IMAGE>" 
        echo ""
        echo "Command line options:"
        echo "-r Use the 4D residual file for smoothness estimation (more accurate)."
        echo "-d Specify degrees of freedom for 4D residual based smoothness estimation (obligatory with -r)."
        echo ""
        echo "Examples:" 
        echo "run_ptfce.sh zstat.nii.gz mask.nii.gz"
        echo "run_ptfce.sh -r res4d.nii.gz -d 32 zstat.nii.gz mask.nii.gz"
        echo ""
        echo "Author: Tamas Spisak"
        echo "https://spisakt.github.io/pTFCE/"
    exit 1
fi


IN=$1 # Z-score image
MASK=$2 # maks image

BASE=`basename $IN .nii.gz`

if [ -z "$residual" ]
then
	smoothest -z $IN -m $MASK > smoothness_full_$BASE.txt
else
	if [ -z "$dof" ]
	then
		echo "Please specify degrees-of-freedom with the -d option!"
		echo "Process terminated."
		exit 1
	fi
	smoothest -d $dof -r $residual -m $MASK > smoothness_full_$BASE.txt
fi

# for compatibility with FSL5
head -n3 smoothness_full_$BASE.txt > smoothness_$BASE.txt

#calculate number of resels
Rscript -e "x=read.table(\"smoothness_$BASE.txt\");numRES=x[2,2]/x[3,2];write.table(numRES, \"numRESELs_$BASE.txt\", row.names=F, col.names=F)"
#calculate z threshold			
ptoz 0.05 -g `cat numRESELs_$BASE.txt` > thres_z_$BASE.txt



# run ptfce 
echo "pTFCE running:"

time Rscript -e "library(pTFCE);library(oro.nifti);Z=readNIfTI(\"$IN\", reorient = FALSE);MASK=readNIfTI(\"$MASK\", reorient = FALSE);x=read.table(\"smoothness_$BASE.txt\");V=x[2,2];Rd=x[1,2]*x[2,2];write.table(Rd, \"Rd_$BASE.txt\", row.names=F, col.names=F);PTFCE=ptfce(Z, MASK, V=V, Rd=Rd, resels=x[3,2]);writeNIfTI(PTFCE\$Z, \"pTFCE_Z_$BASE\");write.table(PTFCE\$fwer0.05.Z, \"ptfce_zthr.txt\", row.names = F, col.names = F)"

echo "
"
echo "FWER Z threshold: " `cat thres_z_$BASE.txt`
echo "View results overlaid on the standard temaplate with:"
echo "fsleyes  $FSLDIR/data/standard/MNI152_T1_1mm_brain.nii.gz -cm greyscale pTFCE_Z_$BASE  -dr `cat thres_z_$BASE.txt` `fslstats pTFCE_Z_$BASE -p 100` -cm red-yellow&"


