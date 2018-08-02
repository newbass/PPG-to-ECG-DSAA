#!/bin/bash
# echo $1
# echo $2
path_to_train=$1
path_to_test=$2
path_to_save=$3
# echo $path_to_train
# echo $path_to_test
# echo $path_to_save
matlab  -nosplash  -nodesktop  -nodisplay  -nojvm  -r  "ppg2ecg('$path_to_train', '$path_to_test', '$path_to_save'); quit;"
