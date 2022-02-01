#!/bin/sh
#SBATCH --export=NONE
#SBATCH -N 1            # nodes number (=NBP)
#SBATCH -n 1            # CPUs number (on all nodes) (=NBP*TPN)
#SBATCH -t 04:00:00     # time limit

# Echo des commandes
ulimit -c 0
ulimit -s unlimited
# Arrete du job des la premiere erreur
set -e
# Nom de la machine
hostname


cas='RICO' #'IHOP' #'RICO' #'RICO' #'FIRE'
simu='Ru0x0' #'REF05' #'Ru0x0' # 'Ls2x0' #'L25.6'
files=(002 004 006 008 010 012 014)

name=V0301

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

typs=(downdraft,updraft,downdraft,downdraft)
#tracer=SVT003,SVT001,SVT002_WT
#tracer=WT,WT
tracer=SVT006,SVT004_WT,SVT004_WT,SVT005_WT
for typ in "${typs[@]}"
do
for file in "${files[@]}"
do
echo 'All'
echo $PWD
python stats_objects_tophat_v2.py $typ $tracer $simu $file $cas $name 'or'
done
done

