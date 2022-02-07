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

subcloud=0
cloud=0

typs=(updraft)
tracer=${svt[0]}_WT
for typ in "${typs[@]}"
do
for hour in "${hours[@]}"
do
echo '1'
echo $PWD
echo $tracer
python stats_flux.py $typ $tracer $sens $hour $case $name $subcloud $cloud
done
done


typs=(updraft,downdraft,downdraft,downdraft)
tracer=${svt[0]}_WT,${svt[0]}_WT,${svt[1]}_WT,${svt[2]}_WT
for typ in "${typs[@]}"
do
for hour in "${hours[@]}"
do
echo 'All'
echo $PWD
echo $tracer
#python stats_flux.py $typ $tracer $sens $hour $case $name 'or'
done
done

