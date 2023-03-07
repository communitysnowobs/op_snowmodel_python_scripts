#! /bin/bash

#this is a test file, changing dates to reduce latency.

#This script will execute the complete assimilative snowmodel run for the central OR
#domain. See inline comments for further explanation. This version of the code does a 
#complete run from Oct 1 of current water year to today - (3days). Only the last time
#slice is extracted, saved, and uploaded.

#Note, we need to do a number of things that are related to the date of the model run.
#this includes a reset the time axis to start 1 Oct of present water year. Get today's
#date info and figure out year of current water year. Let's do that up front in
#this script. We do this since the model run takes a long time (12+ hours) and may
#conclude TOMORROW and not today. We want today's date.
day=$(date '+%d')
month=$(date '+%b')
monthnum=$(date '+%m')
year=$(date '+%Y')
if [ $((10#$monthnum)) -lt 10 ]
then
	year=$(($year - 1))
elif [ $((10#$monthnum)) -eq 10 ]	
then
if [ $day -lt 3 ] #changed from 4 to 3 since reducing latency by one day. Jan 2023
then
    year=$(($year - 1))
fi
fi

#also, we are ultimately going to extract just the 'last' time slice, so let's figure
#out the time stamp for that. Get date string from two days ago

# changed from 3 days to 2 days - change on 23 Nov 2022
d=$(date --date="2 days ago" '+%d')
m=$(date --date="2 days ago" '+%m')
y=$(date --date="2 days ago" '+%Y')
STAMP="${y}_${m}_${d}"

################################
#run the query for met data
cd /nfs/depot/cce_u1/hill/dfh/op_snowmodel/op_snowmodel_python_scripts
source /nfs/attic/dfh/miniconda/bin/activate ee
#changed name of script on line below from or.py to or2.py
ipython met_data_or2.py
conda deactivate

################################
#run assim snowmodel...(review that script closely for details, there is a lot going on)
cd /nfs/depot/cce_u1/hill/dfh/op_snowmodel/op_snowmodel_python_scripts
source /nfs/attic/dfh/miniconda/bin/activate snowmodelcal
#changed name of script on line below from or.py to or2.py
ipython assim_or2.py
conda deactivate

echo
echo "snowmodel has finished"
echo

################################
#time to clean up a bit...delete ssmt and sspr grads files, both in 
#the wo_assim and the wi_assim folders. We just don't need them. 

smpath="/nfs/depot/cce_u1/hill/dfh/op_snowmodel/or_snowmodel/"
#define output path on scratch
outpath="/scratch/op_snowmodel_outputs/OR2/"

rm "${smpath}outputs/wi_assim/ssmt.gdat"
rm "${smpath}outputs/wi_assim/sspr.gdat"
rm "${smpath}outputs/wo_assim/ssmt.gdat"
rm "${smpath}outputs/wo_assim/sspr.gdat"

################################
#Next, a TON of data repackaging.
#convert the grads output to .nc
/scratch/cdo/bin/cdo -f nc import_binary "${smpath}ctl_files/wo_assim/swed.ctl" "${smpath}ctl_files/wo_assim/swed.nc"
/scratch/cdo/bin/cdo -f nc import_binary "${smpath}ctl_files/wi_assim/swed.ctl" "${smpath}ctl_files/wi_assim/swed.nc"
/scratch/cdo/bin/cdo -f nc import_binary "${smpath}ctl_files/wo_assim/snod.ctl" "${smpath}ctl_files/wo_assim/snod.nc"
/scratch/cdo/bin/cdo -f nc import_binary "${smpath}ctl_files/wi_assim/snod.ctl" "${smpath}ctl_files/wi_assim/snod.nc"

echo
echo "done converting to nc"
echo

#clean up more...we can deled snod and swed gdats (huge files) since we have the .nc now
rm "${smpath}outputs/wi_assim/swed.gdat"
rm "${smpath}outputs/wi_assim/snod.gdat"
rm "${smpath}outputs/wo_assim/swed.gdat"
rm "${smpath}outputs/wo_assim/snod.gdat"

#use cdo to reset the time axis.
/scratch/cdo/bin/cdo settaxis,$year-10-01,00:00:00,1days "${smpath}ctl_files/wo_assim/swed.nc" "${smpath}ctl_files/wo_assim/swed2.nc"
/scratch/cdo/bin/cdo settaxis,$year-10-01,00:00:00,1days "${smpath}ctl_files/wi_assim/swed.nc" "${smpath}ctl_files/wi_assim/swed2.nc"
/scratch/cdo/bin/cdo settaxis,$year-10-01,00:00:00,1days "${smpath}ctl_files/wo_assim/snod.nc" "${smpath}ctl_files/wo_assim/snod2.nc"
/scratch/cdo/bin/cdo settaxis,$year-10-01,00:00:00,1days "${smpath}ctl_files/wi_assim/snod.nc" "${smpath}ctl_files/wi_assim/snod2.nc"

echo
echo " done resetting time axis"
echo

#next, we are going to mess with the .nc metadata. CDO treats files as lon / lat, whereas
#we have projected coords out of SnowModel.
# DO SWED (WO_ASSIM) FIRST
################################
#First, change variable and dimension names
infile1="${smpath}ctl_files/wo_assim/swed2.nc"
outfile1="${smpath}ctl_files/wo_assim/swed3.nc"
ncrename -O -v .lat,projection_y_coordinate -d .lat,projection_y_coordinate -v .lon,projection_x_coordinate -d .lon,projection_x_coordinate "${infile1}" "${outfile1}"

#Next, change attributes
ncatted -O -a long_name,projection_y_coordinate,o,c,y "${outfile1}"
ncatted -O -a standard_name,projection_y_coordinate,o,c,y "${outfile1}"
ncatted -O -a long_name,projection_x_coordinate,o,c,x "${outfile1}"
ncatted -O -a standard_name,projection_x_coordinate,o,c,x "${outfile1}"

#finally, to change units requires that we know something from the .par file. See
#my writeup on code changes regarding the grid options for the ctl files. Currently,
#line 151 of the par file has the flag that controls the grid (index, m, km, etc.)
parfile="${smpath}snowmodel.par"
unitchoice=$(sed -n '151p' "${parfile}")
unitchoice="${unitchoice:0:1}"

#check on unitchoice and change .nc file appropriately
if [ $unitchoice -eq 1 ]
then
ncatted -O -a units,projection_y_coordinate,o,c,index "${outfile1}"
ncatted -O -a units,projection_x_coordinate,o,c,index "${outfile1}"
elif [ $unitchoice -eq 2 ]
then
ncatted -O -a units,projection_y_coordinate,o,c,meters "${outfile1}"
ncatted -O -a units,projection_x_coordinate,o,c,meters "${outfile1}"
elif [ $unitchoice -eq 3 ]
then
ncatted -O -a units,projection_y_coordinate,o,c,kilometers "${outfile1}"
ncatted -O -a units,projection_x_coordinate,o,c,kilometers "${outfile1}"
elif [ $unitchoice -eq 4 ]
then
ncatted -O -a units,projection_y_coordinate,o,c,meters "${outfile1}"
ncatted -O -a units,projection_x_coordinate,o,c,meters "${outfile1}"
elif [ $unitchoice -eq 5 ]
then
ncatted -O -a units,projection_y_coordinate,o,c,kilometers "${outfile1}"
ncatted -O -a units,projection_x_coordinate,o,c,kilometers "${outfile1}"
fi

# DO SWED (WI_ASSIM) next
################################
#First, change variable and dimension names
infile2="${smpath}ctl_files/wi_assim/swed2.nc"
outfile2="${smpath}ctl_files/wi_assim/swed3.nc"
ncrename -O -v .lat,projection_y_coordinate -d .lat,projection_y_coordinate -v .lon,projection_x_coordinate -d .lon,projection_x_coordinate "${infile2}" "${outfile2}"

#Next, change attributes
ncatted -O -a long_name,projection_y_coordinate,o,c,y "${outfile2}"
ncatted -O -a standard_name,projection_y_coordinate,o,c,y "${outfile2}"
ncatted -O -a long_name,projection_x_coordinate,o,c,x "${outfile2}"
ncatted -O -a standard_name,projection_x_coordinate,o,c,x "${outfile2}"

#check on unitchoice and change .nc file appropriately
if [ $unitchoice -eq 1 ]
then
ncatted -O -a units,projection_y_coordinate,o,c,index "${outfile2}"
ncatted -O -a units,projection_x_coordinate,o,c,index "${outfile2}"
elif [ $unitchoice -eq 2 ]
then
ncatted -O -a units,projection_y_coordinate,o,c,meters "${outfile2}"
ncatted -O -a units,projection_x_coordinate,o,c,meters "${outfile2}"
elif [ $unitchoice -eq 3 ]
then
ncatted -O -a units,projection_y_coordinate,o,c,kilometers "${outfile2}"
ncatted -O -a units,projection_x_coordinate,o,c,kilometers "${outfile2}"
elif [ $unitchoice -eq 4 ]
then
ncatted -O -a units,projection_y_coordinate,o,c,meters "${outfile2}"
ncatted -O -a units,projection_x_coordinate,o,c,meters "${outfile2}"
elif [ $unitchoice -eq 5 ]
then
ncatted -O -a units,projection_y_coordinate,o,c,kilometers "${outfile2}"
ncatted -O -a units,projection_x_coordinate,o,c,kilometers "${outfile2}"
fi

# DO SNOD (WO_ASSIM) next
################################
#First, change variable and dimension names
infile3="${smpath}ctl_files/wo_assim/snod2.nc"
outfile3="${smpath}ctl_files/wo_assim/snod3.nc"
ncrename -O -v .lat,projection_y_coordinate -d .lat,projection_y_coordinate -v .lon,projection_x_coordinate -d .lon,projection_x_coordinate "${infile3}" "${outfile3}"

#Next, change attributes
ncatted -O -a long_name,projection_y_coordinate,o,c,y "${outfile3}"
ncatted -O -a standard_name,projection_y_coordinate,o,c,y "${outfile3}"
ncatted -O -a long_name,projection_x_coordinate,o,c,x "${outfile3}"
ncatted -O -a standard_name,projection_x_coordinate,o,c,x "${outfile3}"

#check on unitchoice and change .nc file appropriately
if [ $unitchoice -eq 1 ]
then
ncatted -O -a units,projection_y_coordinate,o,c,index "${outfile3}"
ncatted -O -a units,projection_x_coordinate,o,c,index "${outfile3}"
elif [ $unitchoice -eq 2 ]
then
ncatted -O -a units,projection_y_coordinate,o,c,meters "${outfile3}"
ncatted -O -a units,projection_x_coordinate,o,c,meters "${outfile3}"
elif [ $unitchoice -eq 3 ]
then
ncatted -O -a units,projection_y_coordinate,o,c,kilometers "${outfile3}"
ncatted -O -a units,projection_x_coordinate,o,c,kilometers "${outfile3}"
elif [ $unitchoice -eq 4 ]
then
ncatted -O -a units,projection_y_coordinate,o,c,meters "${outfile3}"
ncatted -O -a units,projection_x_coordinate,o,c,meters "${outfile3}"
elif [ $unitchoice -eq 5 ]
then
ncatted -O -a units,projection_y_coordinate,o,c,kilometers "${outfile3}"
ncatted -O -a units,projection_x_coordinate,o,c,kilometers "${outfile3}"
fi

# DO SNOD (WI_ASSIM) last
################################
#First, change variable and dimension names
infile4="${smpath}ctl_files/wi_assim/snod2.nc"
outfile4="${smpath}ctl_files/wi_assim/snod3.nc"
ncrename -O -v .lat,projection_y_coordinate -d .lat,projection_y_coordinate -v .lon,projection_x_coordinate -d .lon,projection_x_coordinate "${infile4}" "${outfile4}"

#Next, change attributes
ncatted -O -a long_name,projection_y_coordinate,o,c,y "${outfile4}"
ncatted -O -a standard_name,projection_y_coordinate,o,c,y "${outfile4}"
ncatted -O -a long_name,projection_x_coordinate,o,c,x "${outfile4}"
ncatted -O -a standard_name,projection_x_coordinate,o,c,x "${outfile4}"

#check on unitchoice and change .nc file appropriately
if [ $unitchoice -eq 1 ]
then
ncatted -O -a units,projection_y_coordinate,o,c,index "${outfile4}"
ncatted -O -a units,projection_x_coordinate,o,c,index "${outfile4}"
elif [ $unitchoice -eq 2 ]
then
ncatted -O -a units,projection_y_coordinate,o,c,meters "${outfile4}"
ncatted -O -a units,projection_x_coordinate,o,c,meters "${outfile4}"
elif [ $unitchoice -eq 3 ]
then
ncatted -O -a units,projection_y_coordinate,o,c,kilometers "${outfile4}"
ncatted -O -a units,projection_x_coordinate,o,c,kilometers "${outfile4}"
elif [ $unitchoice -eq 4 ]
then
ncatted -O -a units,projection_y_coordinate,o,c,meters "${outfile4}"
ncatted -O -a units,projection_x_coordinate,o,c,meters "${outfile4}"
elif [ $unitchoice -eq 5 ]
then
ncatted -O -a units,projection_y_coordinate,o,c,kilometers "${outfile4}"
ncatted -O -a units,projection_x_coordinate,o,c,kilometers "${outfile4}"
fi

echo
echo "done editing metadata"
echo

################################
#cool, so our entire water year of SWED and SNOD is now converted over properly to .nc. 
#Let us now extract just the final day and rename it appropriately. First, figure out
#number of time steps in the file.
numsteps=$(/scratch/cdo/bin/cdo -ntime "${outfile4}")
echo $numsteps

#do swed wo_assim
singleday="${smpath}ctl_files/wo_assim/${STAMP}_swed_wo_assim.nc"
/scratch/cdo/bin/cdo -seltimestep,$numsteps $outfile1 $singleday

#do swed wi_assim
singleday="${smpath}ctl_files/wi_assim/${STAMP}_swed_wi_assim.nc"
/scratch/cdo/bin/cdo -seltimestep,$numsteps $outfile2 $singleday

#do snod wo_assim
singleday="${smpath}ctl_files/wo_assim/${STAMP}_snod_wo_assim.nc"
/scratch/cdo/bin/cdo -seltimestep,$numsteps $outfile3 $singleday

#do snod wi_assim
singleday="${smpath}ctl_files/wi_assim/${STAMP}_snod_wi_assim.nc"
/scratch/cdo/bin/cdo -seltimestep,$numsteps $outfile4 $singleday

#clean up a bit
rm "${outfile1}"
rm "${infile1}"
rm "${outfile2}"
rm "${infile2}"
rm "${outfile3}"
rm "${infile3}"
rm "${outfile4}"
rm "${infile4}"

################################
#next, we want to convert this .nc to a geotiff. We can do this with gdal. The synatx:
#>>gdal_translate -of GTiff -a_srs EPSG:xxxx file.nc file.tif
#the -of GTiff requests tiff as output format. The -a_srs EPSG:xxxx sets the projection
#the -a_ullr fixes the weird 'shift' issue we have been having!
fin="${smpath}ctl_files/wo_assim/${STAMP}_swed_wo_assim.nc"
fout="${smpath}ctl_files/wo_assim/${STAMP}_swed_wo_assim.tif"
gdal_translate -q -of GTiff -a_srs EPSG:32610 -a_ullr 570350 4955850 652450 4832750 $fin $fout
#clean up
rm "${fin}"

fin="${smpath}ctl_files/wi_assim/${STAMP}_swed_wi_assim.nc"
fout="${smpath}ctl_files/wi_assim/${STAMP}_swed_wi_assim.tif"
gdal_translate -q -of GTiff -a_srs EPSG:32610 -a_ullr 570350 4955850 652450 4832750 $fin $fout
#clean up
rm "${fin}"

fin="${smpath}ctl_files/wo_assim/${STAMP}_snod_wo_assim.nc"
fout="${smpath}ctl_files/wo_assim/${STAMP}_snod_wo_assim.tif"
gdal_translate -q -of GTiff -a_srs EPSG:32610 -a_ullr 570350 4955850 652450 4832750 $fin $fout
#clean up
rm "${fin}"

fin="${smpath}ctl_files/wi_assim/${STAMP}_snod_wi_assim.nc"
fout="${smpath}ctl_files/wi_assim/${STAMP}_snod_wi_assim.tif"
gdal_translate -q -of GTiff -a_srs EPSG:32610 -a_ullr 570350 4955850 652450 4832750 $fin $fout
#clean up
rm "${fin}"

echo
echo "done converting to geotiff"
echo

################################
#next, let's set values of swe less than 0.001 m to be 'nodata' values. In this way, 
#when we plot the tif, those cells will be transparent. I tested this in QGIS and they do show up 
#transparent. Need to deactivate conda and then reactivate it, in order to access
#gdal_calc. Do same for hs less than 0.01 m
source /nfs/attic/dfh/miniconda/bin/activate cso
fin="${smpath}ctl_files/wo_assim/${STAMP}_swed_wo_assim.tif"
fout="${smpath}ctl_files/wo_assim/${STAMP}_mask_swed_wo_assim.tif"
python /nfs/attic/dfh/miniconda/envs/cso/bin/gdal_calc.py -A $fin --outfile=$fout --calc="A*(A>0.001)" --NoDataValue=0 --quiet
#clean up
rm "${fin}"

fin="${smpath}ctl_files/wi_assim/${STAMP}_swed_wi_assim.tif"
fout="${smpath}ctl_files/wi_assim/${STAMP}_mask_swed_wi_assim.tif"
python /nfs/attic/dfh/miniconda/envs/cso/bin/gdal_calc.py -A $fin --outfile=$fout --calc="A*(A>0.001)" --NoDataValue=0 --quiet
#clean up
rm "${fin}"

fin="${smpath}ctl_files/wo_assim/${STAMP}_snod_wo_assim.tif"
fout="${smpath}ctl_files/wo_assim/${STAMP}_mask_snod_wo_assim.tif"
python /nfs/attic/dfh/miniconda/envs/cso/bin/gdal_calc.py -A $fin --outfile=$fout --calc="A*(A>0.01)" --NoDataValue=0 --quiet
rm "${fin}"

fin="${smpath}ctl_files/wi_assim/${STAMP}_snod_wi_assim.tif"
fout="${smpath}ctl_files/wi_assim/${STAMP}_mask_snod_wi_assim.tif"
python /nfs/attic/dfh/miniconda/envs/cso/bin/gdal_calc.py -A $fin --outfile=$fout --calc="A*(A>0.01)" --NoDataValue=0 --quiet
rm "${fin}"
conda deactivate

echo
echo "done masking"
echo

################################
#cloud optimized geotiff...
fin="${smpath}ctl_files/wo_assim/${STAMP}_mask_swed_wo_assim.tif"
fout="${smpath}ctl_files/wo_assim/${STAMP}_mask_cog_swed_wo_assim.tif"
gdal_translate -q $fin $fout -co TILED=YES -co COPY_SRC_OVERVIEWS=YES -co COMPRESS=DEFLATE
rm -f $fin
fin="${smpath}ctl_files/wo_assim/${STAMP}_swed_wo_assim.tif"
mv $fout $fin
rm -f $fout
gsutil cp $fin gs://cso_test_upload/or_domain/swed_wo_assim/
#let's move it to /scratch and get it off of depot
mv $fin "${outpath}${STAMP}_swed_wo_assim.tif"

fin="${smpath}ctl_files/wi_assim/${STAMP}_mask_swed_wi_assim.tif"
fout="${smpath}ctl_files/wi_assim/${STAMP}_mask_cog_swed_wi_assim.tif"
gdal_translate -q $fin $fout -co TILED=YES -co COPY_SRC_OVERVIEWS=YES -co COMPRESS=DEFLATE
rm -f $fin
fin="${smpath}ctl_files/wi_assim/${STAMP}_swed_wi_assim.tif"
mv $fout $fin
rm -f $fout
gsutil cp $fin gs://cso_test_upload/or_domain/swed_wi_assim/
#let's move it to /scratch and get it off of depot
mv $fin "${outpath}${STAMP}_swed_wi_assim.tif"

fin="${smpath}ctl_files/wo_assim/${STAMP}_mask_snod_wo_assim.tif"
fout="${smpath}ctl_files/wo_assim/${STAMP}_mask_cog_snod_wo_assim.tif"
gdal_translate -q $fin $fout -co TILED=YES -co COPY_SRC_OVERVIEWS=YES -co COMPRESS=DEFLATE
rm -f $fin
fin="${smpath}ctl_files/wo_assim/${STAMP}_snod_wo_assim.tif"
mv $fout $fin
rm -f $fout
gsutil cp $fin gs://cso_test_upload/or_domain/snod_wo_assim/
#let's move it to /scratch and get it off of depot
mv $fin "${outpath}${STAMP}_snod_wo_assim.tif"

fin="${smpath}ctl_files/wi_assim/${STAMP}_mask_snod_wi_assim.tif"
fout="${smpath}ctl_files/wi_assim/${STAMP}_mask_cog_snod_wi_assim.tif"
gdal_translate -q $fin $fout -co TILED=YES -co COPY_SRC_OVERVIEWS=YES -co COMPRESS=DEFLATE
rm -f $fin
fin="${smpath}ctl_files/wi_assim/${STAMP}_snod_wi_assim.tif"
mv $fout $fin
rm -f $fout
gsutil cp $fin gs://cso_test_upload/or_domain/snod_wi_assim/
#let's move it to /scratch and get it off of depot
mv $fin "${outpath}${STAMP}_snod_wi_assim.tif"

echo
echo "done uploading to GCS"
echo

