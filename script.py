import helpFunctions as hf
import matplotlib.pyplot as plt
import numpy as np
from scipy import misc
from PIL import Image

dirIn = '/home/mehmet/Desktop/DTU_Spring/Math_Modelling/cn/Exercises/Exercise_1/ex1_campusnet/data/'

multiIm, annotationIm = hf.loadMulti('multispectral_day20.mat' , 'annotation_day20.png', dirIn)

my_array = annotationIm[:,:,0]
im = Image.fromarray(annotationIm[:,:,0].reshape((514,514)).astype('uint8')*255)
im.save("grey_scale.png")

my_array2 = annotationIm[:,:,1]
im = Image.fromarray(my_array2.reshape((514,514)).astype('uint8')*255)
im.save("second.png")

my_array3 = annotationIm[:,:,2]
im = Image.fromarray(my_array3.reshape((514,514)).astype('uint8')*255)
im.save("third.png")


##############################################################################################################################
#1. Calculate the threshold value t for all spectral bands for day 1
#to do: for day 1.
multiIm, annotationIm = hf.loadMulti('multispectral_day01.mat' , 'annotation_day01.png', dirIn)
[fatPix, fatR, fatC] = hf.getPix(multiIm, annotationIm[:,:,1])
[meatPix, meatR, meatC] = hf.getPix(multiIm, annotationIm[:,:,2])
thresholds_day20 = np.array([20, 13, 19, 22, 22, 22, 22, 25, 41, 52, 57, 64, 70, 69, 67, 59, 59, 49, 2.5])
thresholds_day01 = np.array([19, 13, 20, 23, 23, 27, 25, 29, 48, 54, 58, 64, 69, 69, 67, 63, 61, 53, 4 ])
def error_rate_fat(fatPix,layerNo):
        #returns the failure rate of the assumption
        return (fatPix[:,0].shape[0]-fatPix[fatPix[:,layerNo]>=20,0].shape[0])/fatPix[:,0].shape[0]
##############################################################################################################################
#2. Calculate the error rate for each spectral band.
for i in range(19):
        error_rate_fat(fatPix,i)
        h = hf.showHistograms(multiIm, annotationIm[:,:,1:3], i, 1)
##############################################################################################################################
#3. Identify the spectral band, which has the best discriminative properties for meat and fat
[fatPix, fatR, fatC] = hf.getPix(multiIm, annotationIm[:,:,1]);
[meatPix, meatR, meatC] = hf.getPix(multiIm, annotationIm[:,:,2]);
plt.plot(np.mean(meatPix,0),'b')
plt.plot(np.mean(fatPix,0),'r')
plt.show()#spectral distrubution
"""
We guess that the third mean difference is the biggest one. So it is a good candidate for training day.
"""


"""
We can see the biggest difference between the means is seen at day 3. However, we can do another better
determination via both increasing the difference between means and decreasing the variance in the groups.
Because the whole idea behind the LDA is this 2 principles. (mu1 - mu2)²/(s1² + s2²). The largest spectrum
is defined as our training set. 
"""
def best_band(meatPix, fatPix):
	numerator = (np.mean(meatPix,0) - np.mean(fatPix,0))**2
	denomunator = np.std(meatPix,0)**2 + np.std(fatPix,0)**2
	
	values = numerator/denomunator
	return np.argmax(values)

"""
We were right. It is the third day has the best indicators.
"""
##############################################################################################################################
#4. Classify the entire image of the salami for day 1, and visualise it


"""
SIMPLEST MODEL (thresholding)
We implement our simplest model to day1. While doing this we do it on the third picture(450 nm).

"""
middle = (annotationIm[:,:,2]*1 + annotationIm[:,:,1]*1 + annotationIm[:,:,0]*1)

def classify(threshold, multiIm):
	from scipy.misc import toimage	
	meat_anno = np.zeros((514,514))
	fat_anno =  np.zeros((514,514))
	for x in range(514):
		for y in range(514):
			if middle[x][y] ==1:
				if multiIm[:,:,2][x][y] >=20:
					meat_anno[x][y] = 0
					fat_anno[x][y] = 1
				else:
					meat_anno[x][y] = 1
					fat_anno[x][y] = 0


thresholds_day01 = np.array([19, 13, 20, 23, 23, 27, 25, 29, 48, 54, 58, 64, 69, 69, 67, 63, 61, 53, 4 ])
for value in thresholds_day01:
	classify(value, multiIm, )

from scipy.misc import toimage				
img = toimage(fat_anno)
img.save("third.png")

img2 = toimage(meat_anno)
img2.save("meatDay1Layer3.png")
##############################################################################################################################
#1. Calculate the multivariate linear discriminant function as described in (23) for day 1. 

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


##############################################################################################################################
#2. Calculate the error rate (disagreement between the model and the annotations) for the training set.




def error_rate():
	import numpy as np
	from scipy import misc
	from PIL import Image
	import helpFunctions as hf
	import matplotlib.pyplot as plt
	dirIn = '/home/mehmet/Desktop/DTU_Spring/Math_Modelling/cn/Exercises/Exercise_1/ex1_campusnet/data/'

	multiIm, annotationIm = hf.loadMulti('multispectral_day01.mat' , 'annotation_day01.png', dirIn)
	[fatPix, fatR, fatC] = hf.getPix(multiIm, annotationIm[:,:,1]);
	[meatPix, meatR, meatC] = hf.getPix(multiIm, annotationIm[:,:,2]);

	imageFat = np.zeros((multiIm.shape[0], multiIm.shape[1]))
	imageMeat = np.zeros((multiIm.shape[0], multiIm.shape[1]))
	
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

	

	errorFatArray = annoFat[(imageFat ==0) & (annoFat==1)] 

	return errorFatArray.shape















