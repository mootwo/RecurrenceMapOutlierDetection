#!/bin/bash

for i in `cat ~/Brain_tumor_work/Projects/Recurrence_2017/RecurrenceMapOutlierDetection/Lists/retroSubs.csv | grep -vw "ID" `
do
	outDir=/cbica/home/rozyckim/Brain_tumor_work/Projects/Recurrence_2017/RecurrenceMapOutlierDetection/Protocols/WM_ROI_Drawing/DataToDrawOn/${i}
	if [ ! -d ${outDir} ]
	then
		mkdir ${outDir}
	fi

	echo ${i}

	qsub -l h_vmem=15G -l short -j y -o ${outDir}/\$JOB_NAME-\$JOB_ID.log \
	/cbica/home/rozyckim/Brain_tumor_work/Projects/Recurrence_2017/RecurrenceMapOutlierDetection/Protocols/WM_ROI_Drawing/GetPerfusionMask.py \
	 -i /cbica/projects/brain_tumor/Brain_Tumor_2015/Projects/Hamed_Recurrence_65/Data/${i}/${i}_PreOp_perf_pp.nii.gz \
	 -r /cbica/home/rozyckim/Brain_tumor_work/Projects/Recurrence_2017/RecurrenceMapOutlierDetection/Data/${i}/${i}_t1ce_scaled.nii.gz \
	 -o ${outDir} 

 done
		


