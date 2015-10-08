#!/usr/bin/env bash

## This short script copies the prepared training and test data to babbage 
## for use with Weka.

scp ./training/training.arff martin@130.216.54.245:~/rf/
scp ./training/training-sub.arff martin@130.216.54.245:~/rf/
scp ./testing/testing.arff martin@130.216.54.245:~/rf/
scp ./testing/testing-sub.arff martin@130.216.54.245:~/rf/
