#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 11 2017

@author: RozyckiM
"""

import sys,os
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
if os.uname()[1] == 'sbia-pc93': #If running program locally
    sys.path.append('/mnt/sbiasfw/external/python/canopy/2.7.9/Canopy_64bit/User/lib/python2.7/site-packages/')
    sys.path.append('/mnt/rozyckim/Scripts/BrainTumorProject')
else:
    sys.path.append('/cbica/home/rozyckim/Scripts/BrainTumorProject')
    
import BT_UsefulFunctions as bt
import nibabel as nb
import getopt

    


def main(sub,outDir):

    #need to calculate vector of means for each feature, as well as covariance matrix between all features
    #This will be used to calculate Mahalanobis distance for each new voxel
    Features = ['t1','t2','t1ce','flair','AX','FA','RAD','TR','PCA1','PCA2','PCA3','PCA4','PCA5']
    
    
    #Load and Scale Each Modality
        #Ignore Perf for now
    for m in range(8):
        mod = Features[m]
        if os.uname()[1] == 'sbia-pc93': #If running program locally
            inNifti = nb.nifti1.load('/mnt/brain_tumor/Brain_Tumor_2015/Projects/Hamed_Recurrence_65/Data/' + sub + '/' + sub + '_PreOp_' + mod + '_pp.nii.gz')
        else:
            inNifti = nb.nifti1.load('/cbica/projects/brain_tumor/Brain_Tumor_2015/Projects/Hamed_Recurrence_65/Data/' + sub + '/' + sub + '_PreOp_' + mod + '_pp.nii.gz')
        f= bt.scaleData(inNifti.get_data())    
        outImg = nb.nifti1.Nifti1Image(f,inNifti.affine)
        nb.nifti1.save(outImg,outDir + '/' + sub + '_' + mod + '_scaled.nii.gz')
        
        
    #Lets load the perfusion PCA Data    
#    for pca in range(5):
#        inNifti =  nb.nifti1.load('/mnt/brain_tumor/Brain_Tumor_2015/Projects/Hamed_Recurrence_65/Protocols/PCA_Perfusion/' + sub + '/' + sub + '_PreOp_perf_pp_PCA_' + str(pca) + '.nii.gz').get_data()
#        f= bt.scaleData(inNifti.getData())
#        outImg = nb.nifti1.Nifti1Image(f,inNifti.affine)
#        nb.nifti1.save(outImg,outDir + '/' + sub + '_PCA' + str(pca) + '_scaled.nii.gz') 
        
    return 0
        
        

        


def parseArguments(arguments):
	print "arguments:",arguments
	try:
		opts, args = getopt.getopt(arguments,"s:o:")
	except getopt.GetoptError as err:
	    	print str(err) # will print something like "option -a not recognized"
       		sys.exit(2)

	for o,a in opts:
         if o == "-s":
             inputSubject = a
         elif o == "-o":
             outputDir = a 				
         else:
             assert False, "unhandled option"

	return inputSubject,outputDir



if __name__ == '__main__':

     S,O = parseArguments(sys.argv[1:])
     print "inputSubject:",S
     print "outputDir:",O



     main(S,O)










    


