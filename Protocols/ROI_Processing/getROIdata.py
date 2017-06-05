#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon May 15 16:24:26 2017

@author: RozyckiM

"""

import sys,os
if os.uname()[1] == 'sbia-pc93': #If running program locally
    homeDir='/mnt/rozyckim'
    btDir = '/mnt/brain_tumor'
else:
    homeDir='/cbica/home/rozyckim'
    btDir = '/cbica/projects/brain_tumor'
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
if os.uname()[1] == 'sbia-pc93': #If running program locally
    sys.path.append('/mnt/sbiasfw/external/python/canopy/2.7.9/Canopy_64bit/User/lib/python2.7/site-packages/')
    
import nibabel as nb
from sklearn.decomposition import PCA


#Get Spreadsheet of subjects
subs = pd.read_csv(homeDir + '/Brain_tumor_work/Projects/Recurrence_2017/RecurrenceMapOutlierDetection/Lists/retroSubs.csv',index_col=0).index
subsTest='AAMG'
subs = subs.drop(['AAAR','AAKQ','AABA','AAMG']) #Three of these don't have segmentations

#need to calculate vector of means for each feature, as well as covariance matrix between all features
#This will be used to calculate Mahalanobis distance for each new voxel

Features = ['t1','t2','t1ce','flair','AX','FA','RAD','TR','PCA1','PCA2','PCA3','PCA4','PCA5']
featN = len(Features) #13
N = len(subs) #31
print "Number of Subjects:",N


#These ROIs are pretty small. Can probably load all values into memory for each subject
#Calculate Number of Voxels in each training set (Pure vs Infiltrated)
pCount = 0
iCount = 0
dCount = 0 #For dilated edema array
pDict = {} #These will hold length of roi for each subject
iDict = {}  
dDict = {}
print "Calculating Size of ROIs"
for s in subs:

    P = nb.nifti1.load(btDir + '/Brain_Tumor_2015/Projects/Hamed_Recurrence_65/Data/' + s + '/' + s + '_PreOp_Pure.nii.gz').get_data()
    pCount = pCount + np.sum(P)
    pDict[s] = np.sum(P)
    
    I = nb.nifti1.load(btDir + '/Brain_Tumor_2015/Projects/Hamed_Recurrence_65/Data/' + s + '/' + s + '_PreOp_Infiltrated.nii.gz').get_data()
    iCount = iCount + np.sum(I)
    iDict[s] = np.sum(I)
    
    D = nb.nifti1.load(homeDir + '/Brain_tumor_work/Projects/Recurrence_2017/RecurrenceMapOutlierDetection/Protocols/WM_ROI_Drawing/DataToDrawOn/' + s + '/' + s + '_WM_ROI.nii.gz').get_data()
    dCount = dCount + np.sum(D).astype(np.int64)
    dDict[s] = np.sum(D).astype((np.int64))
    
print "Size of Pure Array"
print pCount
print "Size of Infiltrated Array"
print iCount
print "Size of Dilated Edema Array"
print dCount


print "Loading Basic Modality Data"

pureArray = np.zeros([featN,pCount],dtype='float')
infArray = np.zeros([featN,iCount],dtype='float')
wmArray = np.zeros([featN,dCount],dtype='float')

iLoc = 0
pLoc = 0 #Location where we're inputing data into feature matrix
dLoc = 0
for s in subs:

    #Get Pure and Inf ROI
    P = nb.nifti1.load(btDir + '/Brain_Tumor_2015/Projects/Hamed_Recurrence_65/Data/' + s + '/' + s + '_PreOp_Pure.nii.gz').get_data()
    P = np.where(P == 1)
    I = nb.nifti1.load(btDir + '/Brain_Tumor_2015/Projects/Hamed_Recurrence_65/Data/' + s + '/' + s + '_PreOp_Infiltrated.nii.gz').get_data()
    I = np.where(I == 1)
    D = nb.nifti1.load(homeDir + '/Brain_tumor_work/Projects/Recurrence_2017/RecurrenceMapOutlierDetection/Protocols/WM_ROI_Drawing/DataToDrawOn/' + s + '/' + s + '_WM_ROI.nii.gz').get_data()
    D = np.where(D == 1)

    #Ignore Perf for now
    for m in range(8):
        mod = Features[m]
        f = nb.nifti1.load(homeDir + '/Brain_tumor_work/Projects/Recurrence_2017/RecurrenceMapOutlierDetection/Data/' + s + '/' + s + '_' + mod + '_scaled.nii.gz').get_data()
        pureArray[m,range(pLoc,(pLoc+len(P[0])))] = f[P]
        infArray[m,range(iLoc,(iLoc+len(I[0])))] = f[I]
        wmArray[m,range(dLoc,(dLoc+len(D[0])))] = f[D]
        
    pLoc = pLoc + len(P[0])
    iLoc = iLoc + len(I[0])
    dLoc = dLoc + len(D[0])
    
     
     
#print "Loading Perfusion Data"
#    
##Calculate PCA for Perfusion Data. Should load coefficients just like Hamed did. Going through the perfusion data takes too
##much memory.  
#    
print "Loading Pure ROI"
#    
#Load perfusion data. But load all pure data first. Then all infiltration data. To allow easier splitting later.
perfArray = np.zeros([pCount + iCount,45],dtype='short')
Loc = 0 #Location where we're inputing data into feature matrix
for s in subs:
    #Get Pure and Inf ROI
    P = nb.nifti1.load(btDir + '/Brain_Tumor_2015/Projects/Hamed_Recurrence_65/Data/' + s + '/' + s + '_PreOp_Pure.nii.gz').get_data()
    P = np.where(P == 1)
    #Read Perfusion
    f = nb.nifti1.load(btDir + '/Brain_Tumor_2015/Projects/Hamed_Recurrence_65/Data/' + s + '/' + s + '_PreOp_perf_pp.nii.gz').get_data()
    Np = len(P[0])
    for t in range(45):
        perfArray[range(Loc,(Loc+Np)),t] = f[P[0],P[1],P[2],t]
        
    Loc = Loc + Np
    
print "Loading Infiltrated ROI"

#Now get load infiltrated roi.
for s in subs:
    #Get Pure and Inf ROI
    I = nb.nifti1.load(btDir + '/Brain_Tumor_2015/Projects/Hamed_Recurrence_65/Data/' + s + '/' + s + '_PreOp_Infiltrated.nii.gz').get_data()
    I = np.where(I == 1)
    #Read Perfusion
    f = nb.nifti1.load(btDir + '/Brain_Tumor_2015/Projects/Hamed_Recurrence_65/Data/' + s + '/' + s + '_PreOp_perf_pp.nii.gz').get_data()
    Ni = len(I[0])
    for t in range(45):
        perfArray[range(Loc,(Loc+Ni)),t] = f[I[0],I[1],I[2],t]
        
    Loc = Loc + Ni
    
    
print "Loading Dilated ROI"
Loc = 0
perfArrayDilated = np.zeros([dCount,45],dtype='short')
for s in subs:
    #Get Pure and Inf ROI
    D = nb.nifti1.load(homeDir + '/Brain_tumor_work/Projects/Recurrence_2017/RecurrenceMapOutlierDetection/Protocols/WM_ROI_Drawing/DataToDrawOn/' + s + '/' + s + '_WM_ROI.nii.gz').get_data()
    D = np.where(D == 1)
    #Read Perfusion
    f = nb.nifti1.load(btDir + '/Brain_Tumor_2015/Projects/Hamed_Recurrence_65/Data/' + s + '/' + s + '_PreOp_perf_pp.nii.gz').get_data()
    Nd = len(D[0])
    for t in range(45):
        perfArrayDilated[range(Loc,(Loc+Nd)),t] = f[D[0],D[1],D[2],t]
        
    Loc = Loc + Nd
    
    
print "Applying PCA Transform to Perfusion Data"
p = PCA(n_components=5 )
p.fit(perfArray) #Calculate principal components.
outFile = homeDir + '/Brain_tumor_work/Projects/Recurrence_2017/RecurrenceMapOutlierDetection/Protocols/ROI_Processing/PCA_27subs.pkl'
from sklearn.externals import joblib
if not os.path.isfile(outFile):
    joblib.dump(p,outFile)
else:
    print "Output Model Already Exisits"
    

#This is the pca-transformed array
coefPerf =p.transform(perfArray) #Nvoxelsx5components
print "Shape PCA Transofmred array",coefPerf.shape
pureCoef = coefPerf[0:pCount,:]
infCoef = coefPerf[pCount:,:]
wmCoef = p.transform(perfArrayDilated)

print "Placing PCA Results into Feature Arrays"

#Now place pca output into pure and infiltrated modality arrays
iLoc = 0
pLoc = 0 #Location where we're inputing data into feature matrix
dLoc = 0
for s in subs:
    for m in range(8,13):
        mod = Features[m]
        pureArray[m,range(pLoc,(pLoc+pDict[s]))] = pureCoef[range(pLoc,(pLoc+pDict[s])),m - 8] #minus 8 since pureCoef have only 5 columns
        infArray[m,range(iLoc,(iLoc+iDict[s]))] = infCoef[range(iLoc,(iLoc+iDict[s])), m - 8]
        wmArray[m,range(dLoc,(dLoc+dDict[s]))] = wmCoef[range(dLoc,(dLoc+dDict[s])), m - 8]
        
    pLoc = pLoc + pDict[s]
    iLoc = iLoc + iDict[s]
    dLoc = dLoc + dDict[s]
    

joinedArray = np.concatenate((infArray,pureArray),axis=1)
joinedArray = np.transpose(joinedArray)
Y = np.zeros(joinedArray.shape[0])
Y[range(infArray.shape[1])] = 1

#Save Joined Array
joinedArray.tofile(homeDir + '/Brain_tumor_work/Projects/Recurrence_2017/RecurrenceMapOutlierDetection/Protocols/ROI_Processing/joinedFeatures_27subs.csv',sep=',') #shape: (41027, 13) We have left a subject out here.
Y.tofile(homeDir + '/Brain_tumor_work/Projects/Recurrence_2017/RecurrenceMapOutlierDetection/Protocols/ROI_Processing/joinedFeatures_27subs_labels.csv',sep=',') #shape: (41027,_

wmArray = np.transpose(wmArray)
wmArray.tofile(homeDir + '/Brain_tumor_work/Projects/Recurrence_2017/RecurrenceMapOutlierDetection/Protocols/ROI_Processing/wmRoiFeatures_27subs.csv',sep=',')











    


