#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  CrtxMsHmn.Ps_6.py
#  
#  Copyright 2019 oscar <oscar@oscar-J53kOiSr>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  

import argparse,os,sys
import scanpyOri as sc

from numpy import ceil
from singlecell.plots_py3 import wrpAll


"""
Script to make figure 1.
"""

#################################
#         Parse inputs          #                
#################################
def str2bool(v):
	if isinstance(v, bool):
		return v
	if v.lower() in ('yes', 'true', 't', 'y', '1'):
		return True
	elif v.lower() in ('no', 'false', 'f', 'n', '0'):
		return False
	else:
		raise argparse.ArgumentTypeError('Boolean value expected.')


parser = argparse.ArgumentParser()

######################
#  Input parameters  #
######################
parser.add_argument('-g','--gnsIntrstFl',help='Input file with gene of interest to plot',default=None)
parser.add_argument('-d','--indxDBsFl',help='Input file with indexed database files',default=None)
parser.add_argument('-5','--h5ad_db_path',help='Path to the folder with h5ad files',default=None)
parser.add_argument('-p','--pagodaApp_path',help='Path to the folder with pagoda rds files',default=None)
parser.add_argument('-i','--icons_path',help='Path to the folder with svg icon files',default=None)
parser.add_argument('-l','--illstrns',help='Path to the folder with png illustration files',default=None)
parser.add_argument('-w','--prwsCmprsnFldr',help='Path to the folder with pairwise comparison folders',default=None)
parser.add_argument('-c','--case_study_left',help='Case study name for the target dataset (left plot)',default=None)
parser.add_argument('-s','--case_study_spp',help='Case study species for the target dataset (left plot)',default=None)
parser.add_argument('-P','--enrchmnt_prfx',help='Enrichment type prefix',default=None)
parser.add_argument('-F','--enrchmnt_fldr',help='Enrichment folder',default=None)
parser.add_argument('-o','--outFig_path',help='Output plot folder',default=None)
parser.add_argument('-L','--illstrtn_case_study_left',help='Name of the png illustration file',default=None)
#########################
#  Optional parameters  #
#########################
parser.add_argument('-t','--FDRTrhshld',help='FDRTrhshld is the FDR threshold to plot the venn diagrams',default=1,type=float)
parser.add_argument('-b','--incldBst',help='Number of best results to include',default=4,type=int)
parser.add_argument('-G','--gnrlScl',help='General scale for subplots',default=10000,type=int)
parser.add_argument('-T','--pltSgnfcncThrshld',help='if not None will includ only significant results for the right reference with FDR<=pltSgnfcncThrshld',default=0.01,type=float)
parser.add_argument('-D','--dpi',help='DPI for output figure file',default=100,type=int)
parser.add_argument('-n','--fontsize',help='Fontsize for output figure file',default=170,type=int)
parser.add_argument('-N','--fontsizeEnrchmnt',help='Fontsize for enrichment inset of the output figure file',default=160,type=int)
parser.add_argument('-M','--fontsizeSmall',help='Small fontsize for output figure file',default=160,type=int)
parser.add_argument('-m','--fontsizeSuperSmall',help='Super small fontsize for output figure file',default=160,type=int)
parser.add_argument('-I','--fontsizeLetterInsert',help='Insert fontsize for output figure file',default=300,type=int)
parser.add_argument('-W','--txtWrpWidth',help='Text wrapper width for output figure file',default=20,type=int)
parser.add_argument('-X','--slctTopSignfcnt',help='Plot only significant genes',default=True,type=str2bool)
parser.add_argument('-k','--pickLwstP',help='Plot only significant enrichments (bottom insert)',default=True,type=str2bool)
parser.add_argument('-Z','--incldFll',help='Include all the significant results (bottom insert)',default=False,type=str2bool)

args = parser.parse_args()

######################
######################
######################
#  Input parameters  #
######################
######################
######################
gnsIntrstFl = args.gnsIntrstFl
indxDBsFl = args.indxDBsFl
h5ad_db_path = args.h5ad_db_path
pagodaApp_path = args.pagodaApp_path
icons_path = args.icons_path
illstrns = args.illstrns
case_study_left = args.case_study_left
case_study_spp = args.case_study_spp
enrchmnt_prfx = args.enrchmnt_prfx
enrchmnt_fldr = args.enrchmnt_fldr
outFig_path = args.outFig_path
prwsCmprsnFldr = args.prwsCmprsnFldr
illstrtn_case_study_left = args.illstrtn_case_study_left
#########################
#  Optional parameters  #
#########################
FDRTrhshld = args.FDRTrhshld
incldBst = args.incldBst
gnrlScl = args.gnrlScl
pltSgnfcncThrshld = args.pltSgnfcncThrshld
dpi = args.dpi
fontsize = args.fontsize
fontsizeEnrchmnt = args.fontsizeEnrchmnt
fontsizeSmall = args.fontsizeSmall
fontsizeSuperSmall = args.fontsizeSuperSmall
fontsizeLetterInsert = args.fontsizeLetterInsert
txtWrpWidth = args.txtWrpWidth
slctTopSignfcnt = args.slctTopSignfcnt
pickLwstP = args.pickLwstP
incldFll = args.incldFll
#~ d_leftClstrScl = {0:12.6*20,1:3.2*20,2:1*20*0.00001,3:1*20*0.00001, \
#~ 4:2*20*0.00001,5:1*20,6:3.6*20,7:1*20}
d_leftClstrScl = None
locationOfInsertD=[1,2]#location in (horizontal,vertical) positions in the output plot
width_ratios=[5,5,8]
height_ratios=[5,15,3]
figsize=(300, 300)

assert os.path.exists(gnsIntrstFl)
assert os.path.exists(indxDBsFl)
assert os.path.exists(h5ad_db_path)
assert os.path.exists(pagodaApp_path)
assert os.path.exists(icons_path)
assert os.path.exists(outFig_path)
assert case_study_spp in {'human','mouse'}
if illstrtn_case_study_left is not None:
	assert os.path.exists(illstrtn_case_study_left)

print('Log info:')
print('\tgnsIntrstFl -> ',gnsIntrstFl)
print('\tindxDBsFl -> ',indxDBsFl)
print('\th5ad_db_path -> ',h5ad_db_path)
print('\tpagodaApp_path -> ',pagodaApp_path)
print('\ticons_path -> ',icons_path)
print('\toutFig_path -> ',outFig_path)
print('\tcase_study_spp -> ',case_study_spp)
print('\tenrchmnt_prfx -> ',enrchmnt_prfx)
print('\tenrchmnt_fldr -> ',enrchmnt_fldr)
print('\toutFig_path -> ',outFig_path)
print('\tprwsCmprsnFldr -> ',prwsCmprsnFldr)
print('\tillstrtn_case_study_left -> ',illstrtn_case_study_left)

dOriClstNmsCmpr={'E13_NC_K6':{'SCP':'SCP_blue_1','Chromaffin_cells': \
'Chromaffin_cells','Sympathoblast':'Sympathoblast','Bridge_cells': \
'Bridge_cells','Unmentioned_1':'SCP_yellow_4','Unmentioned_2':'gray_6'}}

#################
#################
#################
###  Execute  ###
#################
#################
#################
#Read gene list of interest
lGnsIntrstAll = [l.splitlines()[0] for l in open(gnsIntrstFl,'r') if \
l.strip()]
#Dicitonary to add orthologues
dSpplSrtdGnsAdd={'human':['TP63'],'mouse':['Trp63']}
#Read index file
dbsInfo = [[v for v in l.split() if v.strip()] for l in open(indxDBsFl).read(). \
splitlines() if l.strip() and l.find('# ~ ')==-1]
hdr = dbsInfo.pop(0)
nDBs=len(dbsInfo)
#Define positions in file from database
db_case_study_pos = hdr.index('#Case_study')
illstrtnInsrtB_pos = hdr.index('illstrtnInsrt')
spp_pos = hdr.index('spp')
###
#Data for left
###
in_h5adDB_left = os.path.join(h5ad_db_path,'%s.h5ad'%case_study_left)
illstrtn_left = os.path.join(icons_path,'%s.icon.png'%case_study_left)
PAGODAApp_left = os.path.join(pagodaApp_path,'%s.pagodaApp.rds'% \
case_study_left)
lFlsEnrchmnt_left = [os.path.join(enrchmnt_fldr,'%s.%s.tsv'%(case_study_left, \
enrchmnt_prfx))]#Change to any db
spp_left = case_study_spp
outFldr = os.path.join(outFig_path,'%s.d'%enrchmnt_prfx,'%s.d'%case_study_left)
if not os.path.exists(os.path.split(outFldr)[0]):
	os.mkdir(os.path.split(outFldr)[0])
if not os.path.exists(outFldr):
	os.mkdir(outFldr)
#
print('\tcase_study_left -> ',case_study_left)
print('\tspp_left -> ',spp_left)
print('\tin_h5adDB_left -> ',in_h5adDB_left)
print('\tillstrtn_left -> ',illstrtn_left)
print('\tPAGODAApp_left -> ',PAGODAApp_left)
###
#Process left data
###
adata_left=sc.read(in_h5adDB_left)
#Numbers (0-based) of the target clusters of interest in order of plotting 
#from top to bottom
########################################################################
########################################################################
########################################################################
# In case of tsv file (TSV_file_with_cluster_number_0-column_and_cluster_name_1) 
# with cluster order and name use: 
#~ 
#~ dClstrClrs_left = dict([(l.split('\t')[0],l.split('\t')[1]) for l in \
#~ open(TSV_file_with_cluster_number_0-column_and_cluster_name_1).read().splitlines()
#~ if l.strip() and l[0]!='#'])
#~ strd_adata_clstrs_left = adata_left.obs['leiden'].cat.categories.tolist()
#~ srtdClstrs_left = [dClstrClrs_left[c] for c in strd_adata_clstrs_left]
#~ # ~ srtdClstrs_left = strd_adata_clstrs_left
#~ dClstrsClrs = dict([(srtdClstrs_left[pos],adata_left.uns['leiden_colors'][pos]) \
#~ for pos in range(len(strd_adata_clstrs_left))])
#~ adata_left.obs['PAGODA_hc'] = [dClstrClrs_left[c] for c in adata_left.obs['leiden'].tolist()]
#~ adata_left.obs['PAGODA_hc'] = adata_left.obs['PAGODA_hc'].astype('category')
#~ adata_left.obs['cluster_color']=[dClstrsClrs[clstr] for clstr in adata_left.obs['PAGODA_hc'].tolist()]
#~ adata_left.obs['cluster_color']=adata_left.obs['cluster_color'].astype('category')
#~ 
#~ and silence the following:
strd_adata_clstrs_left = adata_left.uns['dendrogram_PAGODA_hc'] \
['categories_idx_ordered'].tolist()
srtdClstrs_left = adata_left.obs['PAGODA_hc'].cat.categories \
[adata_left.uns['dendrogram_PAGODA_hc']['categories_idx_ordered']]. \
tolist()
########################################################################
########################################################################
########################################################################
dClstrClrs_left = dict([(c,adata_left.obs[adata_left.obs \
['PAGODA_hc']==c]['cluster_color'][0]) for c in srtdClstrs_left])

#
if case_study_left=='E13_NC_K6':
	dOriNmsAltrvNms = dOriClstNmsCmpr[case_study_left]
	srtdClstrs_left = [dOriNmsAltrvNms[c] for c in srtdClstrs_left]
	dClstrClrs_left = dict([(dOriNmsAltrvNms[c],v) for c,v in \
	dClstrClrs_left.items()])
#
srtdLClstrs_left = [str(c) for c in srtdClstrs_left]
###
#Process right data
###
for db_right_pos in range(nDBs):
	db_right = dbsInfo[db_right_pos]
	case_study_right = db_right[db_case_study_pos]
	spp_right = db_right[spp_pos]
	assert spp_right in {'human','mouse'}
	###
	#Data for right
	###
	in_h5adDB_right = os.path.join(h5ad_db_path,'%s.h5ad'%case_study_right)
	illstrtn_right = os.path.join(icons_path,'%s.icon.png'%case_study_right)
	PAGODAApp_right = os.path.join(pagodaApp_path,'%s.pagodaApp.rds'% \
	case_study_right)
	lFlsEnrchmnt_right=[os.path.join(enrchmnt_fldr,'%s.%s.tsv'% \
	(case_study_right,enrchmnt_prfx))]#Change to any db
	###
	#Pairwise comparison results
	###
	cllOvrlpngFl = os.path.join(prwsCmprsnFldr,'%s.d'%case_study_left, \
	'%s.tsv'%case_study_right)
	###
	#Further data for right
	###
	adata_right=sc.read(in_h5adDB_right)
	#Numbers (0-based) of the _righterence clusters of interest in order of 
	#plotting from top to bottom
	strd_adata_clstrs_right = adata_right.uns['dendrogram_PAGODA_hc'] \
	['categories_idx_ordered'].tolist()
	srtdClstrs_right = adata_right.obs['PAGODA_hc'].cat.categories \
	[adata_right.uns['dendrogram_PAGODA_hc']['categories_idx_ordered']]. \
	tolist()
	dClstrClrs_right = dict([(c,adata_right.obs[adata_right.obs \
	['PAGODA_hc']==c]['cluster_color'][0]) for c in srtdClstrs_right])
	#
	if case_study_right=='E13_NC_K6':
		dOriNmsAltrvNms = dOriClstNmsCmpr[case_study_right]
		srtdClstrs_right = [dOriNmsAltrvNms[c] for c in srtdClstrs_right]
		dClstrClrs_right = dict([(dOriNmsAltrvNms[c],v) for c,v in \
		dClstrClrs_right.items()])
	#
	srtdLClstrs_right = [str(c) for c in srtdClstrs_right]
	strdClstrLst_right = srtdClstrs_right
	outplt=os.path.join(outFldr,'%s.svg'%case_study_right)
	#
	print('\tcase_study_right -> ',case_study_right)
	print('\tspp_right -> ',spp_right)
	print('\toutplt -> ',outplt)
	print('\tin_h5adDB_right -> ',in_h5adDB_right)
	print('\tillstrtn_right -> ',illstrtn_right)
	print('\tPAGODAApp_right -> ',PAGODAApp_right)
	print('\tcllOvrlpngFl -> ',cllOvrlpngFl)
	#
	if not os.path.exists(outplt) or not os.path.getsize(outplt):
		###
		#Define genes of interest
		###
		avlGns = set([v.upper() for v in adata_left.var_names.tolist()]). \
		intersection(set([v.upper() for v in adata_right.var_names.tolist()]))
		lGnsIntrst = [g for g in lGnsIntrstAll if g.upper() in avlGns]
		print('Genes to plot... %s'%len(lGnsIntrst))
		###
		#Run
		###
		wrpAll(outplt=outplt,db_left=case_study_left,in_h5adDB_left= \
		in_h5adDB_left,illstrtn_left=illstrtn_left,PAGODAApp_left=PAGODAApp_left, \
		lFlsGOEnrchmnt_left=lFlsEnrchmnt_left,db_case_study_left=case_study_left, \
		illstrtnInsrtB=illstrtn_case_study_left,adata_left=adata_left, \
		strd_adata_clstrs_left=strd_adata_clstrs_left,srtdClstrs_left= \
		srtdClstrs_left,dClstrClrs_left=dClstrClrs_left,srtdLClstrs_left= \
		srtdLClstrs_left,db_right_pos=db_right_pos,db_right=db_right, \
		in_h5adDB_right=in_h5adDB_right,illstrtn_right=illstrtn_right, \
		PAGODAApp_right=PAGODAApp_right,lFlsGOEnrchmnt_right= \
		lFlsEnrchmnt_right,db_case_study_right=case_study_right,cllOvrlpngFl= \
		cllOvrlpngFl,adata_right=adata_right,strd_adata_clstrs_right= \
		strd_adata_clstrs_right,srtdClstrs_right=srtdClstrs_right, \
		dClstrClrs_right=dClstrClrs_right,srtdLClstrs_right=srtdLClstrs_right, \
		strdClstrLst_right=strdClstrLst_right,lGnsIntrst=lGnsIntrst, \
		width_ratios=width_ratios,height_ratios=height_ratios,figsize= \
		figsize,dpi=dpi,fontsize=fontsize,fontsizeEnrchmnt=fontsizeEnrchmnt, \
		fontsizeSmall=fontsizeSmall,fontsizeSuperSmall=fontsizeSuperSmall, \
		txtItrvls=int(ceil((len(srtdLClstrs_right)+len(srtdLClstrs_left))/2)), \
		fontsizeLetterInsert=fontsizeLetterInsert,txtWrpWidth=txtWrpWidth, \
		spp_left=spp_left,spp_right=spp_right,d_leftClstrScl=d_leftClstrScl, \
		FDRTrhshld=FDRTrhshld,incldBst=incldBst,gnrlScl=gnrlScl, \
		pltSgnfcncThrshld=pltSgnfcncThrshld,locationOfInsertD= \
		locationOfInsertD,dSpplSrtdGnsAdd=dSpplSrtdGnsAdd,slctTopSignfcnt= \
		slctTopSignfcnt,pickLwstP=pickLwstP,incldFll=incldFll)
