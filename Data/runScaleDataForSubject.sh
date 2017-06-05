#!/bin/bash

for i in `cat ~/Brain_tumor_work/Projects/Recurrence_2017/RecurrenceMapOutlierDetection/prospectSubs.csv | grep -vw "ID" `
do
	outDir=${PWD}/${i}
	if [ ! -d ${outDir} ]
	then
		mkdir ${outDir}
	fi

	qsub -l short -j y -o ${outDir}/\$JOB_NAME-\$JOB_ID.log \
	/cbica/home/rozyckim/Brain_tumor_work/Projects/Recurrence_2017/RecurrenceMapOutlierDetection/Data/ScaleDataForSubject.py \
	 -s ${i} \
	 -o ${outDir} 

 done
		


