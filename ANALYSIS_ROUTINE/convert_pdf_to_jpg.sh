#!/bin/bash

FFSCT=0.047
folder=/data/daniellin/RAMSES_RUNS/EOS_GAMMA_HEATING_COOLING/TR_FLOOR_10K/ISM_FFSCT_${FFSCT}_10000M
slice_source_folder=$folder/figures
slice_target_folder=$folder/figures_jpg
profile_source_folder=$folder/profiles
profile_target_folder=$folder/profiles_jpg

cd $slice_source_folder

for file in *.pdf
do
  echo $file
  output="${file%.*}"
  pdftoppm -singlefile -jpeg -r 300 $file $slice_target_folder/$output
done

cd $slice_target_folder
ffmpeg -i rho_x_0_%05d.jpg -vcodec libx264 \
 -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" -r 24 \
 -y -an $slice_target_folder/rho_vel.mp4

ffmpeg -i amr_x_0_%05d.jpg -vcodec libx264 \
 -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" -r 24 \
 -y -an $slice_target_folder/amr.mp4

ffmpeg -i P_x_0_%05d.jpg -vcodec libx264 \
 -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" -r 24 \
 -y -an $slice_target_folder/pressure.mp4

ffmpeg -i T_x_0_%05d.jpg -vcodec libx264 \
 -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" -r 24 \
 -y -an $slice_target_folder/temperature.mp4


cd $profile_source_folder

for file in *.pdf
do
  echo $file
  output="${file%.*}"
  pdftoppm -singlefile -jpeg -r 300 $file $profile_target_folder/$output
done

ffmpeg -i $profile_target_folder/profiles_%05d.jpg -vcodec libx264 \
 -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" -r 15 \
 -y -an $profile_target_folder/profile.mp4 

