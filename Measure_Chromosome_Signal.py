def VAK(path,channel):
    #Measure value at kinetochore using histone signal to segment. Save two-channel videos and one-channel bg video in the same folder.

	import numpy as np
	import math
	from skimage.external.tifffile import imread
	from skimage.io import imsave
	from scipy import ndimage
	from skimage.morphology import binary_opening as opening
	from skimage.filters import gaussian
	from skimage.filters import threshold_isodata
	from skimage.morphology import binary_dilation as dilation
	import os
	import xlsxwriter
    
	movielist=os.listdir(path)
	print(movielist)
	nbMovies=(len(movielist)-1)/2
	NEBDArray=(np.asarray(open(path+"ZNEBD.txt").read().split())).astype(int)

	header=np.tile(np.arange(1,nbMovies+1)[None,:,None],[1,11])
    
	table=np.zeros([150,nbMovies,11])
    
	for e in range(nbMovies):
		stack=imread(path+movielist[2*e])
		bg=imread(path+movielist[2*e+1])
		signalchannel=stack[:,channel]
		chromosomechannel=stack[:,np.asarray([0,1])[np.asarray([0,1])!=channel][0]]
        
		for t in range(np.shape(chromosomechannel)[0]):
			chromosomechannel[t]=65635*gaussian(chromosomechannel[t],sigma=1)        
		segchromosomes=chromosomechannel>threshold_isodata(chromosomechannel)
		segchromosomes_expanded=dilation(dilation(dilation(segchromosomes)))
		bg_surround=segchromosomes_expanded^segchromosomes
        
		segchrimage,nb=ndimage.measurements.label(segchromosomes)
		imsave(path+'chromosomes_'+movielist[2*e], segchrimage)
		bg_surroundimage,nb=ndimage.measurements.label(bg_surround)
		imsave(path+'bg_surround_'+movielist[2*e], bg_surroundimage)
        
        
		for t in range(np.shape(chromosomechannel)[0]):
			table[t,e,0]=np.count_nonzero(segchromosomes[t])
			table[t,e,1]=np.mean(bg[t])
			table[t,e,2]=np.sum(signalchannel[t][np.where(segchromosomes[t]>0)])
			table[t,e,3]=np.mean(signalchannel[t][np.where(segchromosomes[t]>0)])
			table[t,e,7]=np.mean(signalchannel[t][np.where(bg_surround[t]>0)])
                
	table[:,:,4]=np.multiply(table[:,:,0],table[:,:,1])
	table[:,:,5]=np.divide((table[:,:,2]-table[:,:,4]),table[:,:,1])
	table[:,:,6]=np.divide((table[:,:,3]-table[:,:,1]),table[:,:,1])
	table[:,:,8]=np.multiply(table[:,:,0],table[:,:,7])
	table[:,:,9]=np.divide((table[:,:,2]-table[:,:,8]),table[:,:,1])
	table[:,:,10]=np.divide((table[:,:,3]-table[:,:,7]),table[:,:,1])    
    
    
	#Concatenate Header and Output Table
	table=np.concatenate((header,table),axis=0)
    
	alignedtable=np.zeros([150,nbMovies,11])
	for e in range(nbMovies):
		alignedtable[:50,e,:]=table[NEBDArray[e]:NEBDArray[e]+50,e,:]

	#Concatenate Header and Aligned Output Table
	alignedtable=np.concatenate((header,alignedtable),axis=0)
    
    
	#Export to Excel Sheet
	workbook = xlsxwriter.Workbook(path+"data.xlsx",options={'nan_inf_to_errors': True})
	titles=('number of pixels','BG Average','Signal Sum','Signal Average','BG sum','Norm Signal Sum, -cBG','Norm Signal Average, -cBG','surrounding BG Average','surrounding BG Sum','Norm Signal Sum, -sBG','Norm Signal Average, -sg BG')
	for s in range(11):
		worksheet = workbook.add_worksheet(titles[s])
		col=0
		for row, data in enumerate(table[:,:,s]):
			worksheet.write_row(row, col, data)
            
	#Export to Excel Sheet
	workbook = xlsxwriter.Workbook(path+"aligned_data.xlsx",options={'nan_inf_to_errors': True})
	titles=('number of pixels','BG Average','Signal Sum','Signal Average','BG sum','Norm Signal Sum, -cBG','Norm Signal Average, -cBG','surrounding BG Average','surrounding BG Sum','Norm Signal Sum, -sBG','Norm Signal Average, -sBG')
	for s in range(11):
		worksheet = workbook.add_worksheet(titles[s])
		col=0
		for row, data in enumerate(alignedtable[:,:,s]):
			worksheet.write_row(row, col, data)
