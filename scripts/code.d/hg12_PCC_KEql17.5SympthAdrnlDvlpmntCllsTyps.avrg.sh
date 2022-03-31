#!/bin/bash -l

#SBATCH -A snic2020-5-457
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 3:00:00
#SBATCH -M rackham
#SBATCH -J hg12_PCC_KEql17.5SympthAdrnlDvlpmntCllsTyps.avrg
#SBATCH -o hg12_PCC_KEql17.5SympthAdrnlDvlpmntCllsTyps.avrg.log

source /crex2/proj/sllstore2017016/bioinformatics.d/parallel_versions.d/.profile 
export MPLCONFIGDIR='/crex/proj/sllstore2017016/bioinformatics.d/code.d/figs.d/NB.d.1/CllTypsExtnd.Ps_7.d/code.d/pheos_case.d/hg12_PCC_KEql17.5SympthAdrnlDvlpmntCllsTyps.avrg.tmp/.config/matplotlib'

##########################
##########################
##########################
#####   Variables   ######
##########################
##########################
##########################
#Input file with gene of interest to plot
gnsIntrstFl='gnsIntrst.tsv'
#Input file with indexed database files
indxDBsFl='databases.index'
#Path to the folder with h5ad files
h5ad_db_path='h5ad.d'
#Path to the folder with pagoda rds files
pagodaApp_path='pagodaApp.d'
#Path to the folder with png icon files
icons_path='icons.d'
#Path to the folder with png illustration files
illstrns='illstrtns.d'
#Path to the folder with pairwise comparison folders
prwsCmprsnFldr='prwsCmprsn.avrg.d'
#Case study name for the target dataset (left plot)
case_study_left='hg12_PCC_KEql17'
#Case study species for the target dataset (left plot)
case_study_spp='human'
#Enrichment type prefix
enrchmnt_prfx='5SympthAdrnlDvlpmntCllsTyps'
#Enrichment folder
enrchmnt_fldr='5SympthAdrnlDvlpmntCllsTyps.avrg.d'
#Output plot folder
outFig_path='figs.d'
#Name of the png illustration file

#Special parameters
pltSgnfcncThrshld='0.05'
slctTopSignfcnt=FALSE

python CllTypsExtnd.Ps_7.py -g=$gnsIntrstFl -d=$indxDBsFl -5=$h5ad_db_path -p=$pagodaApp_path -i=$icons_path -l=$illstrns -w=$prwsCmprsnFldr -c=$case_study_left -s=$case_study_spp -P=$enrchmnt_prfx -F=$enrchmnt_fldr -o=$outFig_path  -T=$pltSgnfcncThrshld -X=$slctTopSignfcnt

rm -r hg12_PCC_KEql17.5SympthAdrnlDvlpmntCllsTyps.avrg.tmp
