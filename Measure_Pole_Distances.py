def PoleDistances(path,pole_channel,threshold,exclusion):
     #Measure distance between poles. Save videos and ZNEBD.txt file in the same folder. ZNEBD lists in a column the timing of NEBD, with numbering of frames starting from 0. Specify path, pole signal channel (channel numbering from 0), threshold value for segmentation, and if exclusion of signal in the same channel as the poles needs to be performed.
    
	import numpy as np
	import math
	import os
	from skimage.io import imsave
	from skimage.external.tifffile import imread
	from skimage.measure import regionprops
	from skimage.morphology import binary_closing as close
	from skimage.morphology import binary_dilation as dilation
	from skimage.filters import gaussian
	from scipy import ndimage
	import xlsxwriter
    
	movielist=os.listdir(path)
	print(movielist)
	header=np.tile(np.arange(1,len(movielist))[None,:,None],[1,2])
	table=np.zeros([150,len(movielist)-1,2])
	NEBDArray=(np.asarray(open(path+"ZNEBD.txt").read().split())).astype(int)
    
	for e in range(len(movielist)-1):
		moviepath=path+movielist[e]
		stack=imread(moviepath)   
		poles=stack[:,pole_channel]
		for t in range(np.shape(poles)[0]):
			poles[t]=65635*gaussian(poles[t],sigma=2)
		segpoles=poles>threshold
		poleimage=np.zeros([np.shape(segpoles)[0],np.shape(segpoles)[1],np.shape(segpoles)[2]])
        
		for t in range(np.shape(poles)[0]):
			pole_labels,nbpoles=ndimage.measurements.label(close(segpoles[t]))
			label_list=[]
			for n in range(nbpoles):
				if regionprops(pole_labels)[n].area>20 & regionprops(pole_labels)[n].area<150 :
					label_list.append(n+1)
			polelist=ndimage.measurements.center_of_mass(segpoles[t],pole_labels,label_list)
			polelist=np.asarray(polelist)
			if len(label_list)>1:
				distance=math.sqrt(np.square(polelist[0,0]-polelist[1,0])+np.square(polelist[0,1]-polelist[1,1]))*50/234
			else:
				distance=0
			table[t,e,0]=distance
			pole_labels_exp=np.expand_dims(np.asarray(pole_labels),0)
			poleimage=np.concatenate((poleimage,pole_labels_exp), axis=0)

		imsave(path+'poles_'+movielist[e], poleimage)
		table[:50,e,1]=table[NEBDArray[e]:NEBDArray[e]+50,e,0]
        
	table=np.concatenate((header,table),axis=0)
    
	#Export to Excel Sheet
	workbook = xlsxwriter.Workbook(path+"data.xlsx",options={'nan_inf_to_errors': True})
	titles=('distances','aligned distances')
	for s in range(2):
		worksheet = workbook.add_worksheet(titles[s])
		col=0
		for row, data in enumerate(table[:,:,s]):
			worksheet.write_row(row, col, data)
