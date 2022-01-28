#!/bin/bash
# SBATCH -J plot_3D
# Florent Brient
# > est lancé par la commande >crontab table_cron_FB (suivre les crons avec crontab -l)
# manuellement se lance >sqsub run_plot_LES_3D.sh

path=`pwd` #/Users/florentbrient/Dropbox/MESO-NH/Github/objects-LES/
cd $path'/scripts/'
pathdata=$path'/data/'

echo $1
echo $2
echo $3
simus=($1 $2 $3 'V0301')

#python2.7 plot_LES_3D.py FIRE Ls2x0 024 V0301

hours=${simus[2]}
hours=${hours//,/' '}

for hour in $hours 
do
 echo $hour
 simus2=($pathdata ${simus[0]} ${simus[1]} $hour ${simus[3]})
 python plot_LES_3D.py ${simus2[@]}
done

#simus2=($pathdata ${simus[0]} ${simus[1]} $hour ${simus[3]})
