import helpFunctions as hf

import numpy as np
from scipy import misc
from PIL import Image
	
dirIn = '/home/mehmet/Desktop/DTU_Spring/Math_Modelling/cn/Exercises/Exercise_1/ex1_campusnet/data/'

multiIm, annotationIm = hf.loadMulti('multispectral_day01.mat' , 'annotation_day01.png', dirIn)
[fatPix, fatR, fatC] = hf.getPix(multiIm, annotationIm[:,:,1]);
[meatPix, meatR, meatC] = hf.getPix(multiIm, annotationIm[:,:,2]);


imageFat = np.zeros((multiIm.shape[0], multiIm.shape[1]))
imageMeat = np.zeros((multiIm.shape[0], multiIm.shape[1]))
middle = (annotationIm[:,:,2]*1 + annotationIm[:,:,1]*1 + annotationIm[:,:,0]*1)


def lda_amount(typePix, meatPix, fatPix, multiIm, x, y):
	from numpy.linalg import inv
	covMatrixFat = np.cov(np.transpose(fatPix)) #cov matrix
	covMatrixMeat = np.cov(np.transpose(meatPix))
	
	numo = ( (fatPix.shape[0]-1)*covMatrixFat+(meatPix.shape[0]-1)*covMatrixMeat)
	deno = fatPix.shape[0] + meatPix.shape[0]-2
	covMatrix = numo/deno

	mean = np.mean(typePix,0)
	prior = typePix.shape[0]/(meatPix.shape[0] + fatPix.shape[0])
	
	first1 = np.transpose(multiIm[x,y,:])
	first2 = inv(covMatrix)
	first = np.dot(np.dot( first1, first2  ),mean)

	second = np.dot(1/2*np.dot(np.transpose(mean), inv(covMatrix)), mean)

	third = np.log(prior)

	return  (first - second +third)





def lda_apply(meatPix, fatPix, multiIm, dirIn):
	multiIm, annotationIm = hf.loadMulti('multispectral_day01.mat' , 'annotation_day01.png', dirIn)
	annoFat = annotationIm[:,:,1]
	annoMeat = annotationIm[:,:,2] 
	for i in range(multiIm.shape[0]):
		for j in range(multiIm.shape[1]):
		
			if not middle[i,j] ==1:
				print(i,j)
				continue
			
			valueFat = lda_amount(fatPix, meatPix, fatPix, multiIm, i, j)
			valueMeat = lda_amount(meatPix, meatPix, fatPix, multiIm, i, j)
			print(valueMeat , valueFat)
			if middle[i,j] ==1:
				if valueFat > valueMeat :
					imageFat[i,j] = 1
				else:
					imageMeat[i,j] = 1

	from scipy.misc import toimage
	img = toimage(imageFat)
	img.save("day1Fat.png")
	img2 = toimage(imageMeat)
	img2.save("day1Meat.png")
	
	errorFatArray = annoFat[(imageFat ==0) & (annoFat==1)] 
	
	errorMeatArray = annoMeat[(imageMeat ==0) & (annoMeat==1)] 
	print(errorFatArray.shape)

lda_apply(meatPix, fatPix, multiIm, dirIn)

