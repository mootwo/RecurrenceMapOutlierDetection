# -*- coding: utf-8 -*-
"""
Created on Wed Apr 12 12:02:38 2017

This will examine effectiveness of byte scaling for each subject.

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

    

#Get Spreadsheet of subjects
subs = pd.read_csv('/mnt/rozyckim/Brain_tumor_work/Projects/Recurrence_2017/RecurrenceMapOutlierDetection/retroSubs.csv',index_col=0).index
#For testing, subs is only 5 subjects
subs = subs[0:5]


#need to calculate vector of means for each feature, as well as covariance matrix between all features
#This will be used to calculate Mahalanobis distance for each new voxel

Features = ['t1','t2','t1ce','flair','AX','FA','RAD','TR','PCA1','PCA2','PCA3','PCA4','PCA5']
featN = len(Features) #13
N = len(subs) #31


#For modality, plot a boxplot of the distruction of intensity for each subject.
#Hopefully we will see the normalization equalize things here.
#Ignore Perf for now
for m in range(8):
    mod = Features[m]

    rawDict = {}
    scaleDict = {}
    fig, ax = plt.subplots(nrows=2,ncols=1)
    for s in subs:

        #Get Pure and Inf ROI
        P = nb.nifti1.load('/mnt/brain_tumor/Brain_Tumor_2015/Projects/Hamed_Recurrence_65/Data/' + s + '/' + s + '_PreOp_Pure.nii.gz').get_data()
        P = np.where(P == 1)
#        I = nb.nifti1.load('/mnt/brain_tumor/Brain_Tumor_2015/Projects/Hamed_Recurrence_65/Data/' + s + '/' + s + '_PreOp_Infiltrated.nii.gz').get_data()
#        I = np.where(I == 1)

        f = nb.nifti1.load('/mnt/brain_tumor/Brain_Tumor_2015/Projects/Hamed_Recurrence_65/Data/' + s + '/' + s + '_PreOp_' + mod + '_pp.nii.gz').get_data()
        posInd = np.where(f != 0 )
        fScaled = nb.nifti1.load('/mnt/rozyckim/Brain_tumor_work/Projects/Recurrence_2017/RecurrenceMapOutlierDetection/Data/' + s + '/' + s + '_' + mod + '_scaled.nii.gz').get_data()
        posIndScaled = np.where(fScaled !=0 )
        rawDict[s] =  f[P]
        #rawDict[s] =  np.ravel(f[posInd])
        scaleDict[s] = fScaled[P]
        #scaleDict[s] = np.ravel(fScaled[posIndScaled])

    ax[0].boxplot([rawDict[x] for x in subs])
    ax[0].set_title('Raw ' + mod)
    ax[1].boxplot([scaleDict[x] for x in subs])
    ax[1].set_title('Scale ' + mod)
    plt.show()
   
    
    
        











    


