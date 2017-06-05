#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 11 10:26:36 2017

@author: RozyckiM
"""

import os, sys
import numpy as np
from sklearn.decomposition import PCA
if os.uname()[1] == 'sbia-pc93': #If running program locally
    sys.path.append('/mnt/sbiasfw/external/python/canopy/2.7.9/Canopy_64bit/User/lib/python2.7/site-packages/')
import nibabel as nb

 
 
def run_PCA_Perfusion(inPath,refPath,outDir):
    
    #Read nifti file
    print "Loading input nifti files"
    inNifti = nb.nifti1.load(inPath)
    refNifti = nb.nifti1.load(refPath) #The refernce image that perfusion is registered to.
    refImage = refNifti.get_data()
    #Need to reshape data to have a vector at every timepoint
    I = inNifti.get_data() #Need to flatten, since each voxel counts as a sample
    N = len(np.ravel(I[:,:,:,0]))
    numT = I.shape[3]
    inArray = np.zeros([N,numT]) #The timepoints are the features. We will reduce these to 5 or so.
    
    for t in range(numT):
        inArray[:,t] = np.ravel(I[:,:,:,t])
    
    print "Running PCA"
    p = PCA(n_components=5 )
    p.fit(inArray) #Calculate principal components
    outArray = p.fit_transform(inArray) #Our original data projected onto the new components
    
    #Now write out nifti for first 5 components
    bn=os.path.basename(inPath)
    bn = bn.replace('.nii.gz','')
    
    for i in range(5):
        #Should load header of t1ce image
        outImage = nb.nifti1.Nifti1Image(outArray[:,i].reshape(refImage.shape),refNifti.affine,header=refNifti.header)
        nb.nifti1.save(outImage,outDir + '/' + bn + '_PCA_' + str(i) + '.nii.gz')
    
    
    return 0
    
    
    
    
#Takes Input Array from Nifti File and Scales it Between 0 and 255, exlcuding voxels above 95.5 percentile.
def scaleData(inArray):
        newInd = np.where(inArray != 0 ) #Find all non background voxels.
#        maxThresh = np.percentile(inArray[posInd],99.9) #Exclude voxels above this percentile
        #TODO: Excllusions looked weird. Seems to be excluding normal voxels. Not doing this at the momemnt.
#        newInd = np.where( (inArray != 0) & (inArray < maxThresh ))
        
        #Array to hold scaled values
        newf = np.zeros(inArray.shape,dtype='uint8')
        minF = np.min(inArray[newInd])
        maxF = np.max(inArray[newInd])
        rangeF = maxF - minF        
        
        #Scale data to byte range
        N = len(newInd[0])
        for i in range(N):
            loc = (newInd[0][i],newInd[1][i],newInd[2][i])
            newVal =  float( inArray[loc] - minF ) / rangeF * 254 + 1 #So that minimum signal is 1, and max is 255
            newf[newInd[0][i],newInd[1][i],newInd[2][i]] = newVal
            
        return newf