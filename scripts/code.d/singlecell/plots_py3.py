#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  plots_py3.py
#  
#  Copyright 2018 oscar <oscar@oscar-J53kOiSr>
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

"""
These methods were build to store various method to plots.
"""
#Assert is running on python 3
import sys
assert sys.version_info.major>2#Assert is running on python 3

import colorsys
import matplotlib
matplotlib.use('Agg')
import matplotlib.cm as cm
import matplotlib.colors as clrs
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.gridspec as gridspec
import mpl_toolkits.axisartist as axisartist
import scanpy as sc
import scanpy.external as sce
import seaborn as sns
import textwrap 
import os
import rpy2.robjects.numpy2ri

from numpy import arange,array,ceil,cumsum,float32,floor,isnan,ma,mean,median,nan,sum,random,divide,linspace,object,sqrt,unique,where,zeros
from scipy.stats import mannwhitneyu
from pandas import DataFrame
from numpy import array,float32,log,inf
from matplotlib_venn import venn2,venn2_circles
from rpy2 import robjects
from rpy2.robjects.packages import importr
from scanpy.plotting import heatmapAdjstd



plt.rcParams.update({
    "font.family": 'sans-serif',
    "font.sans-serif": ['Liberation Sans'],  # use a specific sans-serif font
})


#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################                                          ######################
#################      Plotting methods for NB figures     ######################
#################                                          ######################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################


#################################################################################
##################          To compose the figure           #####################
#################################################################################
def wrpAll(outplt,db_left,in_h5adDB_left,illstrtn_left,PAGODAApp_left, \
	lFlsGOEnrchmnt_left,db_case_study_left,illstrtnInsrtB,adata_left, \
	strd_adata_clstrs_left,srtdClstrs_left,dClstrClrs_left,srtdLClstrs_left, \
	db_right_pos,db_right,in_h5adDB_right,illstrtn_right,PAGODAApp_right, \
	lFlsGOEnrchmnt_right,db_case_study_right,cllOvrlpngFl,adata_right, \
	strd_adata_clstrs_right,srtdClstrs_right,dClstrClrs_right,srtdLClstrs_right, \
	strdClstrLst_right,lGnsIntrst,width_ratios,height_ratios,figsize,dpi, \
	fontsize,fontsizeEnrchmnt,fontsizeSmall,fontsizeSuperSmall,txtItrvls, \
	fontsizeLetterInsert,txtWrpWidth,spp_left,spp_right,d_leftClstrScl=None, \
	pltE10=False,FDRTrhshld=1,incldBst=4,gnrlScl=10000,pltSgnfcncThrshld= \
	0.01,locationOfInsertD=[1,2],dSpplSrtdGnsAdd=None,slctTopSignfcnt=True, \
	pickLwstP=True,incldFll=False):
	"""
	Wrapper to plot the figures
	Input: Standard inputs. dSpplSrtdGnsAdd is a dicitonary with gene 
	names to add in case of an paralogous/orthologous name difference.
	"""
	#################################
	### Initial figure parameters ###
	#################################
	print('\tmatplotlib.get_configdir() -> ',matplotlib.get_configdir())
	fig = plt.figure(figsize=figsize,dpi=dpi)#,constrained_layout=True)
	gs = gridspec.GridSpec(ncols=3,nrows=3,width_ratios=width_ratios, \
	height_ratios=height_ratios,figure=fig)
	matplotlib.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']
	matplotlib.rcParams['text.usetex']=True
	matplotlib.rcParams['text.latex.unicode']=True
	matplotlib.rcParams['font.size'] = '140'
	matplotlib.rcParams['font.sans-serif'] = ['Liberation Sans']
	matplotlib.rcParams['font.family'] = 'sans-serif'
	fontsize = fontsize
	fontsizeEnrchmnt = fontsizeEnrchmnt
	fontsizeSmall =fontsizeSmall
	fontsizeSuperSmall = fontsizeSuperSmall
	txtItrvls = txtItrvls
	fontsizeLetterInsert = fontsizeLetterInsert
	wrapper = textwrap.TextWrapper(width=txtWrpWidth)#length text in insert E
	#####################
	### Plot insert B ###
	#####################
	#Insert the figure
	gs0=gridspec.GridSpecFromSubplotSpec(ncols=1,nrows=1,subplot_spec=gs[0,2])
	ax = fig.add_subplot(gs0[0])
	ax.axis('off')
	if illstrtnInsrtB is not None:
		opndIllstrtnInsrtB = mpimg.imread(illstrtnInsrtB)
		ax.imshow(opndIllstrtnInsrtB,interpolation='nearest')
	for sp in ax.spines.values():
		sp.set_visible(False)
	#####################
	### Plot insert A ###
	#####################
	#Insert the small illustration left
	gs0=gridspec.GridSpecFromSubplotSpec(ncols=2,nrows=1,width_ratios=[1,8], \
	subplot_spec=gs[0,0])
	ax = fig.add_subplot(gs0[0,0])
	ax.axis('off')
	ax.imshow(mpimg.imread(illstrtn_left),interpolation='nearest')
	#Make the plots left
	adata_left = sc.read(in_h5adDB_left)
	sc.pl.tsne(adata_left,color=['PAGODA_hc'],sgmntPlt=gs0[0,1], \
	fig=fig,pltName=False,s=20000,use_raw=False, \
	legend_fontsize=fontsize,legend_loc='on data')
	#Insert the small illustration right
	gs0=gridspec.GridSpecFromSubplotSpec(ncols=2,nrows=1,width_ratios=[1,5], \
	subplot_spec=gs[0,1])
	ax = fig.add_subplot(gs0[0,0])
	ax.axis('off')
	ax.imshow(mpimg.imread(illstrtn_right),interpolation='nearest')
	#Make the plots right
	adata_right = sc.read(in_h5adDB_right)
	sc.pl.tsne(adata_right,color=['PAGODA_hc'],sgmntPlt=gs0[0,1], \
	fig=fig,pltName=False,s=20000,use_raw=False, \
	legend_fontsize=fontsize,legend_loc='on data')
	###########################
	###   Plot inserts C    ###
	###########################
	#Gene list to plot
	# ~ plt.savefig(outplt,dpi=fig.dpi,bbox_inches='tight',quality=100,optimize=True)
	# ~ #Make the plots left
	if spp_left=='human':
		lGnsIntrst_left = [g.upper() for g in lGnsIntrst]
		if dSpplSrtdGnsAdd is not None:
			lGnsIntrst_left.extend([v for v in dSpplSrtdGnsAdd['human'] if v in adata_left.var_names])
	else:
		lGnsIntrst_left = [g for g in lGnsIntrst]
		if dSpplSrtdGnsAdd is not None:
			lGnsIntrst_left.extend([v for v in dSpplSrtdGnsAdd['mouse'] if v in adata_left.var_names])
	heatmapAdjstd(adata_left,lGnsIntrst_left,groupby='PAGODA_hc', \
	swap_axes=True,show=False,srtByDndrgm=True,sgmntPlt=gs[1,0],fig=fig, \
	use_raw=False,labelsize=fontsize,**{'cmap':'seismic','vmin':-7, \
	'vmax':7})
	#Make the plots right
	if spp_right=='human':
		lGnsIntrst_right = [g.upper() for g in lGnsIntrst]
		if dSpplSrtdGnsAdd is not None:
			lGnsIntrst_right.extend([v for v in dSpplSrtdGnsAdd['human'] if v in adata_right.var_names])
	else:
		lGnsIntrst_right = [g for g in lGnsIntrst]
		if dSpplSrtdGnsAdd is not None:
			lGnsIntrst_right.extend([v for v in dSpplSrtdGnsAdd['mouse'] if v in adata_right.var_names])
	if pltE10:#only to plot E10 in a different color
		heatmapAdjstd(adata_right,lGnsIntrst_right,groupby='PAGODA_hc', \
		swap_axes=True,show=False,srtByDndrgm=True,sgmntPlt=gs[1,1],fig=fig, \
		show_gene_labels=False,use_raw=False,labelsize=fontsize,**{'cmap': \
		'autumn_r','vmin':0,'vmax':1})
	else:
		heatmapAdjstd(adata_right,lGnsIntrst_right,groupby='PAGODA_hc', \
		swap_axes=True,show=False,srtByDndrgm=True,sgmntPlt=gs[1,1],fig=fig, \
		show_gene_labels=False,use_raw=False,labelsize=fontsize,**{'cmap': \
		'seismic','vmin':-7,'vmax':7})
	#####################
	### Plot insert D ###
	#####################
	pltVnnDgrmsCmpsd(cllOvrlpngFl,'',srtdClstrs_left,adata_left, \
	strdClstrLst_right=strdClstrLst_right,adata_right=adata_right,illstrtn_left= \
	illstrtn_left,illstrtn_right=illstrtn_right,srtdClstrs_left=srtdClstrs_left, \
	srtdClstrs_right=srtdClstrs_right,incldBst=incldBst,FDRTrhshld=FDRTrhshld, \
	d_leftClstrScl=d_leftClstrScl,gnrlScl=gnrlScl,strd_adata_clstrs_left= \
	strd_adata_clstrs_left,strd_adata_clstrs_right=strd_adata_clstrs_right, \
	pltSgnfcncThrshld=pltSgnfcncThrshld,rtrnPlot=True,sgmntPlt= \
	gs[locationOfInsertD[0],locationOfInsertD[1]],fig=fig,fontsize=fontsize, \
	fontsizeSmall=fontsizeSmall,incldFDR=True,sense='uprgltd',slctTopSignfcnt= \
	slctTopSignfcnt)
	##########################
	### Plot insert E bars ###
	##########################
	#FDR for each GO and cluster:
	dUpRgltd_left = mkDctryClstrGOsgnfcnt(lFlsGOEnrchmnt_left, \
	sense='uprgltd',pickLwstP=pickLwstP,incldFll=incldFll)
	dDownRgltd_left = mkDctryClstrGOsgnfcnt(lFlsGOEnrchmnt_left, \
	sense='dwnrgltd',pickLwstP=pickLwstP,incldFll=incldFll)
	dUpRgltd_right = mkDctryClstrGOsgnfcnt(lFlsGOEnrchmnt_right, \
	sense='uprgltd',pickLwstP=pickLwstP,incldFll=incldFll)
	dDownRgltd_right = mkDctryClstrGOsgnfcnt(lFlsGOEnrchmnt_right, \
	sense='dwnrgltd',pickLwstP=pickLwstP,incldFll=incldFll)
	#Total number of cells for each cell type
	#Make the plots left
	adata_left = sc.read(in_h5adDB_left)
	adata_right = sc.read(in_h5adDB_right)
	dClstrCllCnt_left = mkDctryClstrCllCnt(adata_left)
	dClstrCllCnt_right = mkDctryClstrCllCnt(adata_right)
	#Make the plots left
	gs0=gridspec.GridSpecFromSubplotSpec(ncols=1,nrows=2,subplot_spec=gs[2,0])
	####################################################################
	####################################################################
	#############  -Silence- only use to make plots in fig. 3  ##########
	####################################################################
	####################################################################
	####################################################################
	xlimsLmt,lPosGOnames_up = mkPlotEnrchmntBrs(dict([(c,1) for c in \
	dClstrCllCnt_left if dClstrCllCnt_left[c]]),dUpRgltd_left,srtdLClstrs_left, \
	'',dClstrClrs_left,sgmntPlt=gs0,fig=fig,posIn=0,pltNmbrs=True,fontsize= \
	fontsizeEnrchmnt)
	lPosGOnames = [v for v in lPosGOnames_up if v.strip()]
	dummy,lPosGOnames_down = mkPlotEnrchmntBrs(dict([(c,1) for c in \
	dClstrCllCnt_left if dClstrCllCnt_left[c]]),dDownRgltd_left, \
	srtdLClstrs_left,'',dClstrClrs_left,sgmntPlt=gs0,fig=fig,posIn=1,xlims= \
	xlimsLmt,pltNmbrs=True,crrntCntr=len(lPosGOnames),sense='negative', \
	fontsize=fontsizeEnrchmnt)
	lPosGOnames.extend([v for v in lPosGOnames_down if v.strip()])
	####################################################################
	####################################################################
	####################################################################
	# ~ xlimsLmt,lPosGOnames_up = mkPlotEnrchmntBrs(dClstrCllCnt_left,dUpRgltd_left,srtdLClstrs_left,'',dClstrClrs_left,sgmntPlt=gs0,fig=fig,posIn=0,pltNmbrs=True,fontsize=fontsizeEnrchmnt)
	# ~ lPosGOnames = [v for v in lPosGOnames_up if v.strip()]
	# ~ dummy,lPosGOnames_down = mkPlotEnrchmntBrs(dClstrCllCnt_left,dDownRgltd_left,srtdLClstrs_left,'',dClstrClrs_left,sgmntPlt=gs0,fig=fig,posIn=1,xlims=xlimsLmt,pltNmbrs=True,crrntCntr=len(lPosGOnames),sense='negative',fontsize=fontsizeEnrchmnt)
	# ~ lPosGOnames.extend([v for v in lPosGOnames_down if v.strip()])
	#Make the plots right
	gs0=gridspec.GridSpecFromSubplotSpec(ncols=1,nrows=2,subplot_spec=gs[2,1])
	xlimsLmt,lPosGOnames_up = mkPlotEnrchmntBrs(dClstrCllCnt_right, \
	dUpRgltd_right,srtdLClstrs_right,'',dClstrClrs_right,sgmntPlt=gs0,fig= \
	fig,posIn=0,pltNmbrs=True,crrntCntr=len(lPosGOnames),fontsize= \
	fontsizeEnrchmnt)
	lPosGOnames.extend([v for v in lPosGOnames_up if v.strip()])
	dummy,lPosGOnames_down = mkPlotEnrchmntBrs(dClstrCllCnt_right, \
	dDownRgltd_right,srtdLClstrs_right,'',dClstrClrs_right,sgmntPlt=gs0, \
	fig=fig,posIn=1,xlims=xlimsLmt,pltNmbrs=True,crrntCntr=len(lPosGOnames), \
	sense='negative',fontsize=fontsizeEnrchmnt)
	lPosGOnames.extend([v for v in lPosGOnames_down if v.strip()])
	###########################
	### Plot insert E names ###
	###########################
	assert lPosGOnames_up is not None and lPosGOnames_down is not None
	n_cls = int(ceil(len(lPosGOnames)/txtItrvls))
	gs0 = gridspec.GridSpecFromSubplotSpec(ncols=n_cls,nrows=1, \
	subplot_spec=gs[2,2])
	print('\t\tPloting insert E names')
	for cClmn in range(n_cls):
		ax = fig.add_subplot(gs0[cClmn])
		print('\t\tcClmn,len(lPosGOnames),txtItrvls,n_cls -> ',cClmn, \
		len(lPosGOnames),txtItrvls,n_cls)
		indxAdd = cClmn*txtItrvls
		txt_insrtE = '\n'.join([r'\textbf{[%s]}  %s'%(pos+1+indxAdd, \
		str(wrapper.fill(text=textwrap.dedent(text=value.replace('_',' '). \
		split(' -')[0])))) for pos,value in enumerate(lPosGOnames[indxAdd: \
		indxAdd+txtItrvls])])
		#~ txt_insrtE = '\n'.join(['[%s]  %s'%(pos+1+indxAdd,str(wrapper.fill(text=textwrap.dedent(text=value.replace('_',' '))))) for pos,value in enumerate(lPosGOnames[indxAdd:indxAdd+txtItrvls])])
		ax.text(0, 0.5, txt_insrtE, ha='left',va='center', wrap=True, \
		fontsize=fontsizeSuperSmall)
		ax.axis('off')
	#####################
	#### Save figure ####
	#####################
	plt.savefig(outplt,dpi=fig.dpi,bbox_inches='tight',quality=100,optimize=True)
	plt.close()


#################################################################################
######################### To make venn diagrams in D ############################
#################################################################################
robjects.r('''
cntCllsPAGODARDS = function(PAGODAAppFl)
	{
	#Method to count the number of cells in an application
	#file. Returns a number of cells with colors and the 
	#colors sorted
	#
	library(colorRamps)
	print('Reading cells results in PAGODA RDS file...')
	opndFl = readRDS(PAGODAAppFl)
	clls = opndFl$results$colcol
	l2cols <- matlab.like2(length(unique(c(clls))))
	return(list(clls,l2cols))
	}
''')

cntCllsPAGODARDS = robjects.r['cntCllsPAGODARDS']

def mxClrs(clr1,clr2):
	"""
	Method to mix colors.
	Input: clr1 and clr2 are colors.
	Output: rgb_to_hex(mxdClr) is the hexadecimal from the color mix,
	"""
	def hex_to_rgb(value):
		value = value.lstrip('#')
		lv = len(value)
		return tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, \
		lv // 3))
	#
	def rgb_to_hex(rgb):
		return '#%02x%02x%02x' % rgb
	#
	try:
		r1,g1,b1 = hex_to_rgb(clr1)
	except:
		if clr1[0]=='#':
			r1,g1,b1 = hex_to_rgb(clr1[:7])
		else:
			r1,g1,b1 = hex_to_rgb(clrs.cnames[clr1])
	try:
		r2,g2,b2 = hex_to_rgb(clr2)
	except:
		if clr2[0]=='#':
			r2,g2,b2 = hex_to_rgb(clr2[:7])
		else:
			r2,g2,b2 = hex_to_rgb(clrs.cnames[clr2])
	mxdClr = min(int((r1+r2)/2), 256), min(int((g1+g2)/2), 256), \
	min(int((b1+b2)/2), 256) 
	return rgb_to_hex(mxdClr)

# Take negative and positive data apart and cumulate
def get_cumulated_array(data, **kwargs):
	"""
	Method to get cumulative array
	"""
	cum = data.clip(**kwargs)
	cum = cumsum(cum, axis=0)
	d = zeros(data.shape)
	d[1:] = cum[:-1]
	return d

def rtrnDctryClstrNmbrClrCllCnts(adata,updtNmE13K6=True):
	"""
	Method to return colors and counts of cells by color.
	Input: adata is the input adata_ application.
	Output: dClstrNmbrClr is a dictionary of cluster number and colors.
	dClstrNmbrCllsCnts is a dictionary of cluster number and cell 
	counts.
	"""
	allClstrs=adata.obs['PAGODA_hc'].cat.categories.tolist()
	dClstrNmbrCllsCnts = dict([(c,len(adata.obs[adata.obs['PAGODA_hc']==c])) \
	for c in allClstrs])
	dClstrNmbrClr = dict([(c,adata.obs[adata.obs['PAGODA_hc']==c] \
	['cluster_color'][0]) for c in allClstrs])
	if updtNmE13K6 and 'Unmentioned_1' in allClstrs and 'Unmentioned_2' \
		in allClstrs:
		dOriNmsAltrvNms ={'SCP':'SCP_blue_1','Chromaffin_cells': \
		'Chromaffin_cells','Sympathoblast':'Sympathoblast','Bridge_cells': \
		'Bridge_cells','Unmentioned_1':'SCP_yellow_4','Unmentioned_2': \
		'gray_6'}
		dClstrNmbrCllsCnts = dict([(dOriNmsAltrvNms[c],v) for c,v in \
		dClstrNmbrCllsCnts.items()])
		dClstrNmbrClr = dict([(dOriNmsAltrvNms[c],v) for c,v in \
		dClstrNmbrClr.items()])
	return dClstrNmbrClr,dClstrNmbrCllsCnts

def rtrnDctryClstrsStsts(cllOvrlpngFl,incldBst=1,FDRTrhshld=0.05, \
	strd_adata_clstrs_left=None,sense=None):
	"""
	Method to make a dictionary of clusters and significant overlapping 
	{target cluster:{_righterence clusters:(gene counts target, gene counts 
	_righterence, gene counts intersecting and probabilities)}}.
	Input: cllOvrlpngFl is the fil with statistics (standard format).
	incldBst is the number of top most significant _righterence clusters to
	return. FDRTrhshld is the threshold to include in the results.
	strd_adata_clstrs_left is a list of clusters in order of the stats 
	file but _rightering to the PAGODA colors for the target dataset. 
	Output: A dictionary of clusters and significant "incldBst" 
	overlapping {target cluster:{_righterence clusters:(gene counts target,
	gene counts _righterence, gene counts intersecting and probabilities)}}. 
	"""
	if sense is not None:
		assert sense in {'dwnrgltd','uprgltd','upOdwn'}
	#
	dctryClstrsStstsTmp,dClstrFDRClstr_right = {},{}
	dClstrNmbrPAGODAClr = dict([(str(c+1),n) for c,n in \
	enumerate(strd_adata_clstrs_left)])
	#work with the input file
	hdr=True
	for l in open(cllOvrlpngFl,'r'):
		if l.strip():
			l=l.splitlines()[0]
			l=l.split('\t')
			if hdr:
				posIntrsctn = l.index( \
				'genes cluster of interest in cell type of interest')
				posInClstrNtCllTyp = l.index( \
				'genes cluster of interest NOT in cell type of interest')
				posInCllTypNtClstr = l.index( \
				'genes NOT in cluster of interest in cell type of interest')
				posClstr = l.index('cluster')
				posCllTyp = l.index('cell type')
				posSgnfnc = l.index('significance')
				posFDR = l.index('FDR')
				hdr=False
			elif sense is None or (l[posClstr].split('_')[-1].strip()==sense \
				and l[posCllTyp].split('_')[-1].strip()==sense):
				FDR = float32(l[posFDR])
				if FDR<=FDRTrhshld:
					clstr_left = l[posClstr]
					clstr_right = l[posCllTyp]
					if sense is not None:
						clstr_left=clstr_left.split('_%s'%sense)[0]
						clstr_right=clstr_right.split('_%s'%sense)[0]
					nGns_left = int(l[posInClstrNtCllTyp])
					nGns_right = int(l[posInCllTypNtClstr])
					nGnsIntrsct = int(l[posIntrsctn])
					#Add values
					if clstr_left in dctryClstrsStstsTmp:
						dctryClstrsStstsTmp[clstr_left][clstr_right]=(nGns_left, \
						nGns_right,nGnsIntrsct,FDR)
						dClstrFDRClstr_right[clstr_left].append((FDR, \
						clstr_right))
					else:
						dctryClstrsStstsTmp[clstr_left] = {clstr_right: \
						(nGns_left,nGns_right,nGnsIntrsct,FDR)}
						dClstrFDRClstr_right[clstr_left]=[(FDR,clstr_right)]
	#Make dictionary with incldBst	
	dctryClstrsStsts = {}
	for clstr_left in dClstrFDRClstr_right.keys():
		lFDRs_right = sorted(dClstrFDRClstr_right[clstr_left])
		for clstr_rightTmp in lFDRs_right[:incldBst]:
			clstr_right = clstr_rightTmp[1]
			if clstr_left in dctryClstrsStsts:
				dctryClstrsStsts[clstr_left][clstr_right] = \
				dctryClstrsStstsTmp[clstr_left][clstr_right]
			else:
				dctryClstrsStsts[clstr_left] = {clstr_right: \
				dctryClstrsStstsTmp[clstr_left][clstr_right]}
	#
	return dctryClstrsStsts

def mkDillstrtn(outPlt,dClstrNmbrClr_left,dClstrNmbrCllsCnts_left, \
	dClstrNmbrNm_left,dClstrNmbrClr_right,dClstrNmbrCllsCnts_right, \
	dctryClstrsStsts,dClstrNmNmbr_right,illstrtn_left=None,illstrtn_right=None, \
	srtdClstrs_left=None,srtdClstrs_right=None,d_leftClstrScl=None,gnrlScl=1, \
	fontsize=60,fontsizeSmall=50,incldOnlySignfcntThrsld=None,sclOnGns= \
	True,rtrnPlot=True,sgmntPlt=None,fig=None,incldFDR=True,stckPlt=True, \
	slctTopSignfcnt=False):
	"""
	Core method to compose the figure of insert D.
	Input: outPlt is the output plot that includes all the elements.
	dClstrNmbrClr_left is a dictionary of cluster number and colors for
	the target dataset. dClstrNmbrCllsCnts_left is a dictionary of 
	cluster number and cell counts for the target dataset. 
	dClstrNmbrNm_left is a dictionary of cluster number and name for the
	target dataset. dClstrNmbrClr_right is a dictionary of cluster number 
	and colors for the _righterence dataset. dClstrNmbrCllsCnts_right is a 
	dictionary of cluster number and cell counts for the _righterence 
	dataset. dctryClstrsStsts is dictionary {target cluster:{_righterence 
	clusters:(gene counts target,gene counts _righterence, gene counts 
	intersecting and probabilities)}}. dClstrNmNmbr_right is a dictionary
	of cluster names and number for the _righterence dataset. illstrtn_left 
	is a(n optional) illustration to be included in the taget column, 
	illstrtn_right is a(n optional) illustration to be included in the 
	_righterence column, srtdClstrs_left is a(n optional) sorted list of 
	clusters in the target data. srtdClstrs_right is a(n optional) sorted 
	list of categories/clusters in the _righterence data. d_leftClstrScl is 
	the scale (xTimes) of each cluster. gnrlScl is the scale to multiply
	the axis. incldOnlySignfcntThrsld if not None will includ only 
	significant results for the _righterence with FDR<=
	incldOnlySignfcntThrsld. sclOnGns will scale the _righterence sizes
	based on genes.
	"""
	#Set up initial values
	if srtdClstrs_left is None:
		srtdClstrs_left = sorted(dClstrNmbrNm_left.keys())
	else:
		dClstrNmNmbr_left = dict([(v,k) for k,v in \
		dClstrNmbrNm_left.items()])
	if srtdClstrs_right is None:
		srtdClstrs_right = sorted(dClstrNmNmbr_right.keys())
	#Initiate heights with icons space
	heights = [100*gnrlScl]
	#Estimate heights
	if d_leftClstrScl is not None:
		heights.extend([d_leftClstrScl[c]*gnrlScl for c in \
		srtdClstrs_left])
	else:
		heights.extend([80*gnrlScl for c in srtdClstrs_left])
	#Close heights with space for the rule and calculate number of rows
	heights.append(15*gnrlScl)#1/10 of the space
	nRows = len(heights)
	#Establish widths and calcualte number of columns
	if incldFDR:
		if stckPlt:
			widths = [4*gnrlScl,5*gnrlScl,4*gnrlScl,1*gnrlScl]#Venn diagram width
		else:
			widths = [4*gnrlScl,10*gnrlScl,4*gnrlScl,1*gnrlScl]#Venn diagram width
	else:
		if stckPlt:
			widths = [4*gnrlScl,5*gnrlScl,4*gnrlScl]#Venn diagram width
		else:
			widths = [4*gnrlScl,10*gnrlScl,4*gnrlScl]#Venn diagram width
	nClmns = len(widths)
	#Set up the template
	if sgmntPlt is None:
		fig = plt.figure(figsize=(30, 110),dpi=100)#,constrained_layout=True)
		gs0=gridspec.GridSpec(ncols=nClmns,nrows=nRows,width_ratios=widths, \
		height_ratios=heights,figure=fig)#,left=0.05*gnrlScl,right=0.48*gnrlScl, \
	else:
		gs0=gridspec.GridSpecFromSubplotSpec(ncols=nClmns,nrows=nRows, \
		width_ratios=widths,height_ratios=heights,subplot_spec=sgmntPlt)
		#,left=0.05*gnrlScl,right=0.48*gnrlScl, \
		#wspace=0.05*gnrlScl)#1/20 of the space for margins
	#Include icons in image
	cRow = 0
	ax = fig.add_subplot(gs0[cRow,0])
	ax.axis('off')
	if illstrtn_left is not None:
		img_leftIcon = mpimg.imread(illstrtn_left)
		ax.imshow(img_leftIcon,interpolation='nearest')
	for sp in ax.spines.values():
		sp.set_visible(False)
	ax = fig.add_subplot(gs0[cRow,2])
	ax.axis('off')
	if illstrtn_right is not None:
		img_rightIcon = mpimg.imread(illstrtn_right)
		ax.imshow(img_rightIcon)
	for sp in ax.spines.values():
		sp.set_visible(False)
	#Make plots for bars representing cell number in target
	mxVl_left=int(floor(max(dClstrNmbrCllsCnts_left.values())/100)+1)*100#set maximum value Target
	cClmn = 0 #current column 
	for clstr_left in srtdClstrs_left:
		cRow+=1
		ax = fig.add_subplot(gs0[cRow,cClmn])
		ax.axis('off')
		#~ if dctryClstrsStsts.has_key(clstr_left):#Deprecated in python 3
		if clstr_left in dctryClstrsStsts:
			nClls = dClstrNmbrCllsCnts_left[clstr_left]
			_leftClr = dClstrNmbrClr_left[clstr_left]
			eqvlnt = nClls/float32(mxVl_left)
			scld = (1-eqvlnt)*gnrlScl
			ax.axvspan(scld,1*gnrlScl,color=_leftClr)
			ax.annotate(clstr_left.replace('_',' '), \
			(0.5*gnrlScl,0.5),ha='center',va='center', \
			fontsize=fontsize)
			ax.set_xlim(0,1*gnrlScl)
	#Add bar with counts in target
	cRow+=1
	ax = fig.add_subplot(gs0[cRow,cClmn])
	ax.set_xlim([0,gnrlScl])
	ax.set_xlabel(r'\textbf{Number of cells}',fontsize=fontsize)
	ax.hlines(0.1,xmin=0,xmax=gnrlScl,linewidth=20)
	ax.axes.get_yaxis().set_ticklabels([])
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.spines['bottom'].set_visible(False)
	ax.spines['left'].set_visible(False)
	ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, pos: \
	"{:,}".format(mxVl_left-int(x*mxVl_left/float32(gnrlScl)))))
	#Make plots for bars representing cell number in _righterences
	cRow = 0#reset row
	if dClstrNmbrCllsCnts_right is not None:
		cClmn=2
		mxVl_right=int(floor(max(dClstrNmbrCllsCnts_right.values())/100)+1)*100#set maximum value _righterence
		for clstr_left in srtdClstrs_left:
			cRow+=1
			#~ if dctryClstrsStsts.has_key(clstr_left):#Deprecated in python 3
			if clstr_left in dctryClstrsStsts:
				d_rightClstrs = dctryClstrsStsts[clstr_left]
				#########
				strRghtByFDR = sorted([(d_rightClstrs[c][3],c) for c in \
				set(d_rightClstrs.keys())])
				srtdClstrs_right = [c[1] for c in strRghtByFDR]
				dClstrNmNmbr_right = dict([(c,n) for n,c in \
				enumerate(srtdClstrs_right)])
				########
				if incldOnlySignfcntThrsld is not None:
					intrsctng_rightC = [c for c in srtdClstrs_right if c in \
					set(d_rightClstrs.keys()) if d_rightClstrs[c][3]<= \
					incldOnlySignfcntThrsld]
				else:
					intrsctng_rightC = [c for c in srtdClstrs_right if c in \
					set(d_rightClstrs.keys())]
				if slctTopSignfcnt:
					minFDRvl = min([d_rightClstrs[c][3] for c in \
					srtdClstrs_right if c in set(d_rightClstrs.keys())])
					if incldOnlySignfcntThrsld is not None:
						intrsctng_rightC = [c for c in srtdClstrs_right if c in \
						set(d_rightClstrs.keys()) and d_rightClstrs[c][3]<= \
						incldOnlySignfcntThrsld and d_rightClstrs[c][3]== \
						minFDRvl]
					else:
						intrsctng_rightC = [c for c in srtdClstrs_right if c in \
						set(d_rightClstrs.keys()) and d_rightClstrs[c][3]== \
						minFDRvl]
				n_rightClstrs = len(intrsctng_rightC)
				#Work on cell counts				
				if sclOnGns:
					gnCntTtl = [max(d_rightClstrs[c]) for c in intrsctng_rightC]
					#
					if stckPlt:
						height_ratios = [1 for gnCnt in gnCntTtl]
					else:
						height_ratios = [gnCnt/float32(sum(gnCntTtl)) \
						for gnCnt in gnCntTtl]
					#
					if gnCntTtl:
						axx = gridspec.GridSpecFromSubplotSpec(n_rightClstrs,1, \
						subplot_spec=gs0[cRow,cClmn],wspace=0.001, \
						height_ratios=height_ratios)
					else:
						axx = None
				else:
					if n_rightClstrs:
						axx = gridspec.GridSpecFromSubplotSpec(n_rightClstrs,1, \
						subplot_spec=gs0[cRow,cClmn],wspace=0.001)
					else:
						axx = None
				#
				cnt_rightC = -1
				for _rightCNm in intrsctng_rightC:
					cnt_rightC+=1
					_rightCNbr  = dClstrNmNmbr_right[_rightCNm]
					try:
						nClls_right = dClstrNmbrCllsCnts_right[_rightCNm]
					except:
						import pdb
						pdb.set_trace()
					#
					_rightClr = dClstrNmbrClr_right[_rightCNm]
					eqvlnt = nClls_right/float32(mxVl_right)
					scld = eqvlnt*gnrlScl
					axxx = fig.add_subplot(axx[cnt_rightC])
					axxx.set_xlim(0,gnrlScl)
					axxx.axvspan(0,scld,color=_rightClr)
					axxx.annotate(_rightCNm.replace('_',' '),(0.5*gnrlScl, \
					0.5),ha='center',va='center',fontsize=fontsize)
					axxx.axis('off')
		#Add bar with counts in _righterence
		cRow+=1
		axx = gridspec.GridSpecFromSubplotSpec(1,1,subplot_spec= \
		gs0[cRow,cClmn],wspace=0.001)
		axxx = fig.add_subplot(axx[0])
		axxx.set_xlim([0,gnrlScl])
		axxx.set_xlabel(r'\textbf{Number of cells}',fontsize=fontsize)
		axxx.hlines(0.1,xmin=0,xmax=gnrlScl,linewidth=20)
		axxx.axes.get_yaxis().set_ticklabels([])
		axxx.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, pos: \
		"{:,}".format(int(x*mxVl_right/float32(gnrlScl)))))
		axxx.spines['top'].set_visible(False)
		axxx.spines['right'].set_visible(False)
		axxx.spines['bottom'].set_visible(False)
		axxx.spines['left'].set_visible(False)
	#Get stack limits
	if stckPlt:#plot stack bars
		xLimMin,xLimMax=inf,-inf
		cClmn=1
		cRow = 0#reset row
		for clstr_left,clstr_leftNm in enumerate(srtdClstrs_left):
			cRow+=1
			_leftClr = dClstrNmbrClr_left[clstr_leftNm]
			if clstr_leftNm in dctryClstrsStsts:
				d_rightClstrs = dctryClstrsStsts[clstr_leftNm]
				#########
				strRghtByFDR = sorted([(d_rightClstrs[c][3],c) for c in \
				set(d_rightClstrs.keys())])
				srtdClstrs_right = [c[1] for c in strRghtByFDR]
				dClstrNmNmbr_right = dict([(c,n) for n,c in \
				enumerate(srtdClstrs_right)])
				########
				if incldOnlySignfcntThrsld is not None:
					intrsctng_rightC = [c for c in srtdClstrs_right if c in \
					set(d_rightClstrs.keys()) if d_rightClstrs[c][3]<= \
					incldOnlySignfcntThrsld]
				else:
					intrsctng_rightC = [c for c in srtdClstrs_right if c in \
					set(d_rightClstrs.keys())]
				if slctTopSignfcnt:
					minFDRvl = min([d_rightClstrs[c][3] for c in \
					srtdClstrs_right if c in set(d_rightClstrs.keys())])
					if incldOnlySignfcntThrsld is not None:
						intrsctng_rightC = [c for c in srtdClstrs_right if c in \
						set(d_rightClstrs.keys()) and d_rightClstrs[c][3]<= \
						incldOnlySignfcntThrsld and d_rightClstrs[c][3]== \
						minFDRvl]
					else:
						intrsctng_rightC = [c for c in srtdClstrs_right if c in \
						set(d_rightClstrs.keys()) and d_rightClstrs[c][3]== \
						minFDRvl]
				#
				for _rightCNm in intrsctng_rightC:
					gnCnts_left,gnCnts_right,gnCntsIntrsctn, \
					FDR=d_rightClstrs[_rightCNm]
					#
					if xLimMin>(-gnCnts_left-gnCntsIntrsctn):
						xLimMin=(-gnCnts_left-gnCntsIntrsctn)
					if xLimMax<(gnCnts_right+gnCntsIntrsctn):
						xLimMax=(gnCnts_right+gnCntsIntrsctn)
		#Get the closest bottom
		if xLimMin<inf or xLimMax>-inf:
			xLimMin = int(floor(xLimMin/100)*100)
			xLimMax = int(ceil(xLimMax/100)*100)
	#Make venn diagrams
	cClmn=1
	cRow = 0#reset row
	for clstr_left,clstr_leftNm in enumerate(srtdClstrs_left):
		cRow+=1
		_leftClr = dClstrNmbrClr_left[clstr_leftNm]
		if clstr_leftNm in dctryClstrsStsts:
			d_rightClstrs = dctryClstrsStsts[clstr_leftNm]
			#########
			strRghtByFDR = sorted([(d_rightClstrs[c][3],c) for c in \
			set(d_rightClstrs.keys())])
			srtdClstrs_right = [c[1] for c in strRghtByFDR]
			dClstrNmNmbr_right = dict([(c,n) for n,c in \
			enumerate(srtdClstrs_right)])
			########
			if incldOnlySignfcntThrsld is not None:
				intrsctng_rightC = [c for c in srtdClstrs_right if c in \
				set(d_rightClstrs.keys()) if d_rightClstrs[c][3]<= \
				incldOnlySignfcntThrsld]
			else:
				intrsctng_rightC = [c for c in srtdClstrs_right if c in \
				set(d_rightClstrs.keys())]
			if slctTopSignfcnt:				
				minFDRvl = min([d_rightClstrs[c][3] for c in \
				srtdClstrs_right if c in set(d_rightClstrs.keys())])
				if incldOnlySignfcntThrsld is not None:
					intrsctng_rightC = [c for c in srtdClstrs_right if c in \
					set(d_rightClstrs.keys()) and d_rightClstrs[c][3]<= \
					incldOnlySignfcntThrsld and d_rightClstrs[c][3]== \
					minFDRvl]
				else:
					intrsctng_rightC = [c for c in srtdClstrs_right if c in \
					set(d_rightClstrs.keys()) and d_rightClstrs[c][3]== \
					minFDRvl]
			n_rightClstrs = len(intrsctng_rightC)
			#Work on cell counts
			if sclOnGns:
				gnCntTtl = [max(d_rightClstrs[c]) for c in intrsctng_rightC]
				#
				if stckPlt:
					height_ratios = [1 for gnCnt in gnCntTtl]
				else:
					height_ratios = [gnCnt/float32(sum(gnCntTtl)) \
					for gnCnt in gnCntTtl]
				#
				if gnCntTtl:
					nColsBars=1
					if stckPlt:
						nColsBars=2
					axx = gridspec.GridSpecFromSubplotSpec(n_rightClstrs, \
					nColsBars,subplot_spec=gs0[cRow,cClmn],wspace=0.001, \
					height_ratios=height_ratios)
				else:
					axx = None
			else:
				if n_rightClstrs:
					nColsBars=1
					if stckPlt:
						nColsBars=2
					axx = gridspec.GridSpecFromSubplotSpec(n_rightClstrs, \
					nColsBars,subplot_spec=gs0[cRow,cClmn],wspace=0.001)
				else:
					axx = None
			cnt_rightC = -1
			for _rightCNm in intrsctng_rightC:
				cnt_rightC+=1
				gnCnts_left,gnCnts_right,gnCntsIntrsctn,FDR = \
				d_rightClstrs[_rightCNm]
				#
				if dClstrNmbrClr_right is not None:
					_rightCNbr  = dClstrNmNmbr_right[_rightCNm]
					# ~ _rightClr = dClstrNmbrClr_right[_rightCNbr]
					_rightClr = dClstrNmbrClr_right[_rightCNm]
					mxdClr = mxClrs(_leftClr,_rightClr)
					if stckPlt:#plot stack bars
						halfIntrsctn = gnCntsIntrsctn
						data = array([[-halfIntrsctn],[-gnCnts_left], \
						[halfIntrsctn],[gnCnts_right]])
						data_shape = data.shape						
						cumulated_data = get_cumulated_array(data, min=0)
						cumulated_data_neg = get_cumulated_array(data,max=0)
						row_mask = (data<0)
						cumulated_data[row_mask] = cumulated_data_neg[row_mask]
						data_stack = cumulated_data
						cols=[mxdClr,_leftClr,mxdClr, _rightClr]
						# Re-merge negative and positive data.
						axxx = fig.add_subplot(axx[cnt_rightC,0])
						axxx.set_xlim(xLimMin,0)
						for i in arange(0, 2):
							axxx.barh(arange(data.shape[1]), data[i], \
							left=data_stack[i],color=cols[i],linestyle= \
							'solid',alpha=0.6)
						axxx.annotate(-data[1][0],xy=(data[0][0]+data[1][0],0), \
						ha='left',va='center',fontsize=fontsizeSmall)
						axxx.axis('off')
						#
						axxx = fig.add_subplot(axx[cnt_rightC,1])
						axxx.set_xlim(0,xLimMax)
						for i in arange(2, 4):
							axxx.barh(arange(data.shape[1]), data[i], \
							left=data_stack[i],color=cols[i],linestyle= \
							'solid',alpha=0.6)
						axxx.annotate(data[3][0],xy=(data[3][0],0),ha='right', \
						va='center',fontsize=fontsizeSmall)
						axxx.annotate(data[2][0],xy=(0,0),ha='center', \
						va='bottom',fontsize=fontsizeSmall)
						axxx.axis('off')
					else:#plot venn diagrams
						axxx = fig.add_subplot(axx[cnt_rightC])
						v=venn2(subsets=(gnCnts_left,gnCnts_right, \
						gnCntsIntrsctn),set_colors=[_leftClr,_rightClr], \
						set_labels=['',''])
						venn2_circles(subsets=(gnCnts_left,gnCnts_right, \
						gnCntsIntrsctn),linestyle='solid')
				else:
					if stckPlt:#plot stack bars
						halfIntrsctn = gnCntsIntrsctn
						data = array([[-halfIntrsctn],[-gnCnts_left], \
						[halfIntrsctn],[gnCnts_right]])
						data_shape = data.shape						
						cumulated_data = get_cumulated_array(data, min=0)
						cumulated_data_neg = get_cumulated_array(data, max=0)
						row_mask = (data<0)
						cumulated_data[row_mask] = cumulated_data_neg[row_mask]
						data_stack = cumulated_data
						# Re-merge negative and positive data.
						axxx = fig.add_subplot(axx[cnt_rightC,0])
						axxx.set_xlim(xLimMin,0)
						for i in arange(0, 2):
							if i==1:
								axxx.barh(arange(data.shape[1]),data[i], \
								left=data_stack[i],color=_leftClr, \
								linestyle='solid',alpha=0.6)
							else:
								axxx.barh(arange(data.shape[1]), data[i], \
								left=data_stack[i],linestyle='solid',alpha=0.6)
						axxx.annotate(-data[1][0],xy=(data[0][0]+data[1][0],0), \
						ha='left',va='center',fontsize=fontsizeSmall)
						axxx.axis('off')
						axxx = fig.add_subplot(axx[cnt_rightC,1])
						axxx.set_xlim(0,xLimMax)
						for i in arange(2, 4):
							if i==2:
								axxx.barh(arange(data.shape[1]), data[i], \
								left=data_stack[i], color=_rightClr, \
								linestyle='solid',alpha=0.6)
							else:
								axxx.barh(arange(data.shape[1]), data[i], \
								left=data_stack[i],linestyle='solid',alpha=0.6)
						axxx.annotate(data[3][0],xy=(data[3][0],0),ha='right', \
						va='center',fontsize=fontsizeSmall)
						axxx.annotate(data[2][0],xy=(0,0),ha='center', \
						va='bottom',fontsize=fontsizeSmall)
						axxx.axis('off')
					else:#plot venn diagrams
						axxx = fig.add_subplot(axx[cnt_rightC])
						v=venn2(subsets=(gnCnts_left,gnCnts_right, \
						gnCntsIntrsctn),set_labels=['',''])
						venn2_circles(subsets=(gnCnts_left,gnCnts_right, \
						gnCntsIntrsctn),set_labels=['',''],linestyle='solid')
				if not stckPlt:#plot stack bars
					v.get_label_by_id('11').set_text('')#remove number left
					v.get_label_by_id('10').set_text('')#remove number left
					v.get_label_by_id('01').set_text('')#remove number right
					axxx.annotate(gnCntsIntrsctn,xy=v.get_label_by_id('11'). \
					get_position(),ha='center',va='bottom',fontsize=fontsizeSmall)
					axxx.annotate(gnCnts_left,xy=v.get_label_by_id('10'). \
					get_position(),ha='right',va='center',fontsize=fontsizeSmall)
					axxx.annotate(gnCnts_right,xy=v.get_label_by_id('01'). \
					get_position(),ha='left',va='center',fontsize=fontsizeSmall)
	#Add legend in gene counts
	cRow+=1
	if stckPlt and (xLimMin<inf or xLimMax>-inf):#plot stack bars
		#left
		axx = gridspec.GridSpecFromSubplotSpec(1,2,subplot_spec= \
		gs0[cRow,1],wspace=0.001)
		axxx = fig.add_subplot(axx[0])
		axxx.set_xlim([xLimMin,0])
		axxx.set_xlabel(r'\textbf{Number  }',ha='left',fontsize=fontsize)
		axxx.hlines(0.1,xmin=xLimMin,xmax=0,linewidth=20)
		axxx.axes.get_yaxis().set_ticklabels([])
		axxx.set_xticks([xLimMin,0])
		axxx.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, pos: \
		"{:,}".format(int(-x))))
		axxx.spines['top'].set_visible(False)
		axxx.spines['right'].set_visible(False)
		axxx.spines['bottom'].set_visible(False)
		axxx.spines['left'].set_visible(False)
		#right
		axxx = fig.add_subplot(axx[1])
		axxx.set_xlim([0,xLimMax])
		axxx.set_xlabel(r'\textbf{  of genes}',ha='right',fontsize=fontsize)
		axxx.hlines(0.1,xmin=0,xmax=xLimMax,linewidth=20)
		axxx.axes.get_yaxis().set_ticklabels([])
		axxx.set_xticks([xLimMax])
		axxx.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, pos: \
		"{:,}".format(int(x))))
		axxx.spines['top'].set_visible(False)
		axxx.spines['right'].set_visible(False)
		axxx.spines['bottom'].set_visible(False)
		axxx.spines['left'].set_visible(False)
	elif (xLimMin<inf or xLimMax>-inf):
		ax = fig.add_subplot(gs0[cRow,1])
		ax.set_xlim([0,gnrlScl])
		ax.annotate(r'\textbf{Number of\ngenes}',(0.5*gnrlScl,0.5),ha='center', \
		va='center',fontsize=fontsize)
		ax.axis('off') 
	#Make FDR addition
	if incldFDR and (xLimMin<inf or xLimMax>-inf):
		cClmn=3	
		cRow = 0#reset row	
		for clstr_left,clstr_leftNm in enumerate(srtdClstrs_left):
			cRow+=1
			_leftClr = dClstrNmbrClr_left[clstr_leftNm]
			if clstr_leftNm in dctryClstrsStsts:
				d_rightClstrs = dctryClstrsStsts[clstr_leftNm]
				#########
				strRghtByFDR = sorted([(d_rightClstrs[c][3],c) for c in \
				set(d_rightClstrs.keys())])
				srtdClstrs_right = [c[1] for c in strRghtByFDR]
				dClstrNmNmbr_right = dict([(c,n) for n,c in \
				enumerate(srtdClstrs_right)])
				########
				if incldOnlySignfcntThrsld is not None:
					intrsctng_rightC = [c for c in srtdClstrs_right if c in \
					set(d_rightClstrs.keys()) if d_rightClstrs[c][3]<= \
					incldOnlySignfcntThrsld]
				else:
					intrsctng_rightC = [c for c in srtdClstrs_right if c in \
					set(d_rightClstrs.keys())]
				if slctTopSignfcnt:
					minFDRvl = min([d_rightClstrs[c][3] for c in \
					srtdClstrs_right if c in set(d_rightClstrs.keys())])
					if incldOnlySignfcntThrsld is not None:
						intrsctng_rightC = [c for c in srtdClstrs_right if c in \
						set(d_rightClstrs.keys()) and d_rightClstrs[c][3]<= \
						incldOnlySignfcntThrsld and d_rightClstrs[c][3]== \
						minFDRvl]
					else:
						intrsctng_rightC = [c for c in srtdClstrs_right if c in \
						set(d_rightClstrs.keys()) and d_rightClstrs[c][3]== \
						minFDRvl]

				n_rightClstrs = len(intrsctng_rightC)
				#Work on cell counts
				if sclOnGns:
					gnCntTtl = [max(d_rightClstrs[c]) for c in intrsctng_rightC]
					#
					if stckPlt:
						height_ratios = [1 for gnCnt in gnCntTtl]
					else:
						height_ratios = [gnCnt/float32(sum(gnCntTtl)) \
						for gnCnt in gnCntTtl]
					#
					if gnCntTtl:
						axx = gridspec.GridSpecFromSubplotSpec(n_rightClstrs,1, \
						subplot_spec=gs0[cRow,cClmn],wspace=0.001, \
						height_ratios=height_ratios)
					else:
						axx = None
				else:
					axx = gridspec.GridSpecFromSubplotSpec(n_rightClstrs,1, \
					subplot_spec=gs0[cRow,cClmn],wspace=0.001)
				cnt_rightC = -1
				for _rightCNm in intrsctng_rightC:
					cnt_rightC+=1
					gnCnts_left,gnCnts_right,gnCntsIntrsctn,FDR = \
					d_rightClstrs[_rightCNm]
					axxx = fig.add_subplot(axx[cnt_rightC])
					axxx.annotate('%.g'%FDR,(0.5,0),xycoords='axes fraction', \
					ha='center',fontsize=fontsizeSmall-5)
					axxx.axis('off')
	#Add legend in gene counts
	cRow+=1
	if incldFDR and (xLimMin<inf or xLimMax>-inf):
		ax = fig.add_subplot(gs0[cRow,3])
		ax.set_xlim([0,gnrlScl])
		ax.annotate(r'\textit{FDR}',(0.5*gnrlScl,0.5),ha='center', \
		va='center',fontsize=fontsize)
		ax.axis('off') 
	#Save figure
	if rtrnPlot:
		return plt
	else:
		if outPlt:
			plt.savefig(outPlt,dpi=fig.dpi,bbox_inches='tight',quality=100, \
			optimize=True)
			plt.close()
		return 0

def pltVnnDgrmsCmpsd(cllOvrlpngFl,outPlt,strdClstrLst_left,adata_left, \
	strdClstrLst_right=None,adata_right=None,illstrtn_left=None, \
	illstrtn_right=None,srtdClstrs_left=None,srtdClstrs_right=None, \
	incldBst=1,FDRTrhshld=0.05,d_leftClstrScl=None,gnrlScl=1, \
	strd_adata_clstrs_left=None,strd_adata_clstrs_right=None, \
	pltSgnfcncThrshld=None,rtrnPlot=True,sgmntPlt=None,fig=None, \
	fontsize=60,fontsizeSmall=50,incldFDR=True,sense=None, \
	slctTopSignfcnt=False):
	"""
	Method to make composite figures in D including venn diagrams, bars
	with colors, and two dummy images.
	Input: cllOvrlpngFl is the enrichment/depletion statistics file. outPlt 
	is the output plot that includes all the elements. strdClstrLst_left 
	is a list of sorted names of targets (sorted as in color order). 
	adata_left is the adata_ application for the Target. 
	strdClstrLst_right is a list of sorted names of _righterence. adata_right 
	is (optionally) the adata_ application for the _righterence. 
	illstrtn_left is a(n optional) illustration to be included in the 
	taget column, illstrtn_right is a(n optional) illustration to be 
	included in the _righterence column, srtdClstrs_left is a(n optional) 
	sorted list of clusters in the target data. srtdClstrs_right is a(n 
	optional) sorted list of categories/clusters in the _righterence data. 
	incldBst is the maximum number of _righterence significant groups to 
	plot the comparisons. FDRTrhshld is the FDR threshold to plot the 
	venn diagrams. d_leftClstrScl is the scale (xTimes) of each cluster.
	gnrlScl is the scale to multiply the axis. strd_adata_clstrs_left is
	a list of clusters in order of the stats file but _rightering to the 
	adata_ colors for the target dataset. strd_adata_clstrs_right is a list 
	of clusters in order of the stats file but _rightering to the adata_ 
	colors for the _righterence dataset. incldOnlySignfcntThrsld if not 
	None will includ only significant results for the _righterence with 
	FDR<=pltSgnfcncThrshld.
	Output: outPlt is the output plot that includes all the elements. 
	This is in the standard format of insert D in comparison figures.
	"""
	#Set up values cel numbers
	dClstrNmbrClr_left,dClstrNmbrCllsCnts_left = \
	rtrnDctryClstrNmbrClrCllCnts(adata_left)
	dClstrNmbrNm_left = dict([(c,n) for c,n in \
	enumerate(strdClstrLst_left)])
	if adata_right is not None:
		assert strdClstrLst_right is not None
		dClstrNmbrClr_right,dClstrNmbrCllsCnts_right = \
		rtrnDctryClstrNmbrClrCllCnts(adata_right)
	else:
		dClstrNmbrClr_right,dClstrNmbrCllsCnts_right = None,None
	#Return statistics
	dctryClstrsStsts = rtrnDctryClstrsStsts(cllOvrlpngFl,incldBst=incldBst, \
	FDRTrhshld=FDRTrhshld,strd_adata_clstrs_left=strd_adata_clstrs_left,sense=sense)
	#Until here good
	if strdClstrLst_right is not None:
		dClstrNmNmbr_right = dict([(n,c) for c,n in \
		enumerate(strdClstrLst_right)])
	else:
		dClstrNmNmbr_right = {}
		for clstr_leftNm in dctryClstrsStsts:
			dClstrNmNmbr_right.update(dict([(k,k) for k in \
			dctryClstrsStsts[clstr_leftNm].keys()]))
	#Make illustration
	pltVnnDgrm = mkDillstrtn(outPlt,dClstrNmbrClr_left,dClstrNmbrCllsCnts_left, \
	dClstrNmbrNm_left,dClstrNmbrClr_right,dClstrNmbrCllsCnts_right, \
	dctryClstrsStsts,dClstrNmNmbr_right,illstrtn_left,illstrtn_right, \
	srtdClstrs_left,srtdClstrs_right,d_leftClstrScl,gnrlScl=gnrlScl, \
	incldOnlySignfcntThrsld=pltSgnfcncThrshld,rtrnPlot=rtrnPlot, \
	sgmntPlt=sgmntPlt,fig=fig,fontsize=fontsize,fontsizeSmall=fontsizeSmall, \
	incldFDR=incldFDR,slctTopSignfcnt=slctTopSignfcnt)
	if rtrnPlot:
		return pltVnnDgrm
	else:
		return 0


#################################################################################
############################# To make barplot in E ##############################
#################################################################################
#Method to compute the values of interest
def mkDctryClstrGOsgnfcnt(lFlsGO,sense=None,pickLwstP=True,incldFll=False):
	"""
	Method to compute a dictionary of significant categories.
	"""
	if sense is not None:
		assert sense in {'dwnrgltd','uprgltd','upOdwn'}
	#
	dClstr3SgnfnctTmp = {}
	for inFl in lFlsGO:
		hdr = True
		for l in open(inFl,'r'):
			if l.strip():
				l=l.splitlines()[0]
				l=l.split('\t')
				if hdr:
					hdr = False
					indxTrm = l.index('cell type')
					indxClstr = l.index('cluster')
					indxPvl = l.index('p-value')
					indxFDR = l.index('FDR')
					indxSgnfc = l.index('significance')
				elif sense is None or l[indxClstr].split('_')[-1]. \
				strip()==sense:
					sgnfc = l[indxSgnfc]
					if sgnfc.find('*')>-1:
						# ~ trms = '%s (%s)'%(l[indxTrm],os.path.split(inFl)[1].split('_')[1].split('.')[0])
						trms = '%s (%s)'%(l[indxTrm],os.path.split(inFl)[1]. \
						split('.')[1].replace('sc','d')[:2].upper())
						clstr = l[indxClstr]
						if sense is not None:
							clstr=clstr.split('_%s'%sense)[0]
						pVl = float32(l[indxPvl])
						FDRvl = float32(l[indxFDR])
						if FDRvl: 
							FDR = -log(FDRvl)
						else:
							FDR = 200.00
						if clstr in dClstr3SgnfnctTmp:
							dClstr3SgnfnctTmp[clstr].append([pVl,FDR, \
							trms])
						else:
							dClstr3SgnfnctTmp[clstr]=[[pVl,FDR,trms]]
	#Make the final dictionary
	dClstr3Sgnfnct = {}
	sClstr = set(dClstr3SgnfnctTmp.keys())
	while(sClstr):
		clstr = sClstr.pop()
		vals = dClstr3SgnfnctTmp[clstr]
		vals.sort()
		if pickLwstP:
			maxP = max([v[1] for v in vals])
			vals = [v for v in vals if v[1]==maxP]
			if len(vals)>1:
				vals = [(v[0],v[1]-(n*1e-10),v[2]) for n,v in \
				enumerate(vals)]
			dClstr3Sgnfnct[clstr] = [(v[1],v[2]) for v in vals]
		elif incldFll:
			dClstr3Sgnfnct[clstr] = [(v[1],v[2]) for v in vals]
		elif len(vals)>2:
			dClstr3Sgnfnct[clstr] = [(v[1],v[2]) for v in vals[:1]]
			# ~ print(clstr,vals[:3])
		else:
			dClstr3Sgnfnct[clstr] = [(v[1],v[2]) for v in vals]
	#
	return dClstr3Sgnfnct

#Method to compute the cell number in each cluster
def mkDctryClstrCllCnt(adata,updtNmE13K6=True):
	"""
	Method to return a number of cells for each cluster.
	"""
	allClstrs=adata.obs['PAGODA_hc'].cat.categories.tolist()
	dClstrCllCnt = dict([(c,len(adata.obs[adata.obs['PAGODA_hc']==c])) \
	for c in allClstrs])
	if updtNmE13K6 and 'Unmentioned_1' in allClstrs and 'Unmentioned_2' \
		in allClstrs:
		dOriNmsAltrvNms ={'SCP':'SCP_blue_1','Chromaffin_cells': \
		'Chromaffin_cells','Sympathoblast':'Sympathoblast','Bridge_cells': \
		'Bridge_cells','Unmentioned_1':'SCP_yellow_4','Unmentioned_2': \
		'gray_6'}
		dClstrCllCnt = dict([(dOriNmsAltrvNms[c],v) for c,v in \
		dClstrCllCnt.items()])
	return dClstrCllCnt

#Method to plot the figure
def mkPlotEnrchmntBrs(dClstrCllCnt,dDffRgltd,srtdLClstrs,outPltRgltd, \
	dClstrClrs,sgmntPlt=None,fig=None,posIn=None,sense='positive', \
	xlims=None,pltNmbrs=False,crrntCntr=None,fontsize=None):
	"""
	Method to plot the bar enrichment.
	"""
	lHeight,lGONames,lYpos,lWidth,lClrs = [],[],[],[],[]
	crrntPos = 0
	print('\t\tPloting enrichment bars...')
	for clstr in srtdLClstrs:
		totalWidth = float32(dClstrCllCnt[clstr])
		if clstr in dDffRgltd:
			lInfos = dDffRgltd[clstr]
		else:
			lInfos = [[0,'']]
		lenGOt = float32(len(lInfos))
		GOtPosBFnl = totalWidth/lenGOt
		crrntPos += GOtPosBFnl/float32(2)
		#~ crrntPos += 1
		GOtWidth = GOtPosBFnl-2
		for FDR,trm in lInfos:
			print('\t\tclstr,FDR,trm -> ',clstr,FDR,trm)
			lYpos.append(crrntPos)
			lGONames.append(trm)
			lWidth.append(GOtWidth)
			lHeight.append(FDR)
			lClrs.append(dClstrClrs[clstr])
			crrntPos+=GOtPosBFnl
	#Make array
	lHeight,lGONames,lYpos,lWidth,lClrs = array(lHeight),array(lGONames), \
	array(lYpos),array(lWidth),array(lClrs)
	if sense=='negative':
		lHeight=-lHeight
	#Plot
	if outPltRgltd:
		plt = fig.add_subplot(sgmntPlt[0,0])
		plt.subplots(figsize=(62,16))
		plt.bar(lYpos,lHeight,width=lWidth,color=lClrs)
		plt.xticks(lYpos, lGONames,rotation=90)
		plt.savefig(outPltRgltd)
		plt.close()
		return 0
	else:
		ax = fig.add_subplot(sgmntPlt[posIn])
		ax.spines['right'].set_color('none')
		ax.bar(lYpos,lHeight,width=lWidth,color=lClrs)
		ax.set_ylabel('-log(FDR)')
		ax.set_xticklabels('')
		if crrntCntr is None:
			crrntCntr = 0
		for pos in range(len(lYpos)):
			if pltNmbrs:
				if lGONames[pos]:
					crrntCntr+=1
					if sense=='negative':
						ax.text(lYpos[pos],lHeight[pos],'[%s]'%crrntCntr, \
						va='top',ha='center',fontsize=fontsize,rotation=90)
					else:
						ax.text(lYpos[pos],lHeight[pos],'[%s]'%crrntCntr, \
						va='bottom',ha='center',fontsize=fontsize,rotation=90)
			else:
				ax.text(lYpos[pos],lHeight[pos],lGONames[pos],rotation=10, \
				va='bottom')
		fig.canvas.draw()
		if sense=='negative':
			ax.spines['bottom'].set_color('none')
			allVlsYtcks = [v.get_text().replace('-','').replace('−',''). \
			replace('$','') for v in ax.get_yticklabels()]
			ax.set_yticklabels([str(float(v)) for v in allVlsYtcks])
		else:
			ax.set_yticklabels([str(float(v.get_text().replace('$',''))) \
			for v in ax.get_yticklabels()])
			ax.spines['top'].set_color('none')
		if xlims is None:
			xlims=ax.get_xlim()
		else:
			ax.set_xlim(xlims[0],xlims[1])
		ax.xaxis.set_tick_params(rotation=90)
		return (xlims,lGONames)
		

###########################
###########################
###### General plots ######
###########################
###########################

def pltBoxplot(ar_VrblNms,ar_VrblVals,outPlt,posX=0,yLabl="Counts", \
	prsAsIs=False,szX=None,szY=None):
	"""
	Method for plotting box plots.
	Input: ar_VrblNms is an array of variable names. ar_VrblVals is an 
	array of size len(cellNames) x N where N are all values for each
	variable name. posX is the position in ar_VrblNms for the X value.
	prsAsIs is a switch, if True will pass the variables as is.
	Output: Will plot an violin plot for the values in ar_VrblVals.
	"""
	#Make plot
	if szX is None or szY is None:
		f, ax = plt.subplots(figsize=(len(ar_VrblVals)/5,6))
	else:
		f, ax = plt.subplots(figsize=(szX,szY))
	if prsAsIs:
		sns.boxplot(x=ar_VrblNms,y=ar_VrblVals)
	else:
		dctInput = {}
		for vrblNmPos in range(len(ar_VrblNms)):
			dctInput[ar_VrblNms[vrblNmPos]] = ar_VrblVals[:,vrblNmPos]
		#Set dataframe	
		dataF = DataFrame(dctInput)
		vrblNm = ar_VrblNms[posX]
		#~ clrPltt=sns.color_palette(n_colors=len(ar_VrblVals)-1)
		pltPos=-1
		for vrblNmPos in range(len(ar_VrblNms)):
			if vrblNmPos!=posX:
				pltPos+=1
				vrblNmY = ar_VrblNms[vrblNmPos]
				sns.boxplot(x=dataF[vrblNm],y=float32(dataF[vrblNmY]))#, \
				#label=vrblNmY)#,color=clrPltt[pltPos])
	#Write plot
	locs, labels = plt.xticks()
	plt.setp(labels, rotation=90)
	plt.ylabel(yLabl)
	#~ plt.xlabel("Variable")
	ax.legend(ncol=7, loc="upper right", frameon=True)
	#~ sns.despine(left=True, bottom=True)
	plt.tight_layout()
	plt.savefig(outPlt)
	plt.close()
	return 0
