#!/bin/bash

#FILES=/icgc/dkfzlsdf/analysis/B260/projects/HipSci/openAccess/endoderm_differentation/data_raw/*/fastq
FILES=/icgc/dkfzlsdf/analysis/B260/projects/HipSci/openAccess/endoderm_differentation/data_raw/incoming/run_25475/fastq

for f in $FILES

do

if [ -z "$(ls -A $f)" ]; then
   echo "Empty folder: $f"
else
   :
fi

done

#FILES=/icgc/dkfzlsdf/analysis/B260/projects/HipSci/openAccess/endoderm_differentation/data_raw/*/featureCounts/unique_counts
FILES=/icgc/dkfzlsdf/analysis/B260/projects/HipSci/openAccess/endoderm_differentation/data_raw/incoming/run_25475/featureCounts/unique_counts

for f in $FILES

do

if [ -z "$(ls -A $f)" ]; then
   echo "Empty Folder: $f"
else
   :
fi

done


#FOLDERS=/icgc/dkfzlsdf/analysis/B260/projects/HipSci/openAccess/endoderm_differentation/data_raw/*/star/*
FOLDERS=/icgc/dkfzlsdf/analysis/B260/projects/HipSci/openAccess/endoderm_differentation/data_raw/incoming/run_25475/star/* 

for f in $FOLDERS

do

if [ -z "$(ls -A $f)" ]; then
   echo "Empty Folder: $f"
else
   :
fi

done
