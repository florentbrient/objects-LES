#!/bin/bash
#SBATCH --export=NONE
#SBATCH -N 1            # nodes number (=NBP)
#SBATCH -n 1            # CPUs number (on all nodes) (=NBP*TPN)
#SBATCH -t 04:00:00     # time limit

# Echo des commandes
#ulimit -c 0
#ulimit -s unlimited
# Arrete du job des la premiere erreur
#set -e
# Nom de la machine
#hostname

path=`pwd` 
cd $path'/scripts/'
pathdata=$path'/data/'

# case sens hour
echo $1
echo $2
echo $3
echo $4
case=$1
sens=$2
hours=$3
name='V0301'
svt=($4 $5 $6)

typs=(downdraft)
#tracer=SVT003 #_corr
#tracer=WT
tracer=SVT006
for typ in "${typs[@]}"
do
for file in "${files[@]}"
do
echo '1'
#python stats_objects_tophat_v2.py $typ $tracer $simu $file $cas $name
done
done


typs=(updraft)
#tracer=SVT001_WT
#tracer=WT
tracer=SVT004_WT
for typ in "${typs[@]}"
do
for file in "${files[@]}"
do
echo '2'
#python stats_objects_tophat_v2.py $typ $tracer $simu $file $cas $name
done
done

typs=(downdraft)
tracer=SVT004_WT
#tracer=WT
#tracer=SVT004_WT
for typ in "${typs[@]}"
do
for file in "${files[@]}"
do
echo '2 bis'
#python stats_objects_tophat_v2.py $typ $tracer $simu $file $cas $name
done
done


typs=(updraft2)
tracer=SVT002_WT
#tracer=WT
#tracer=SVT004_WT
for typ in "${typs[@]}"
do
for file in "${files[@]}"
do
echo '4'
#python stats_objects_tophat_v2.py $typ $tracer $simu $file $cas $name
done
done

typs=(downdraft)
tracer=SVT005_WT
#tracer=WT
#tracer=SVT004_WT
for typ in "${typs[@]}"
do
for file in "${files[@]}"
do
echo '5'
#python stats_objects_tophat_v2.py $typ $tracer $simu $file $cas $name
done
done

typs=(updraft,downdraft,downdraft,downdraft)
#tracer=SVT003,SVT001,SVT002_WT
#tracer=WT,WT
#tracer=SVT006_WT,SVT004_WT,SVT004_WT,SVT005_WT
tracer=${svt[0]}_WT,${svt[0]}_WT,${svt[1]}_WT,${svt[2]}_WT
for typ in "${typs[@]}"
do
for hour in "${hours[@]}"
do
echo 'All'
echo $PWD
echo $tracer
python stats_flux.py $typ $tracer $sens $hour $case $name 'or'
done
done

