#!/bin/bash

#PBS -A dp033
#PBS -l walltime=12:00:00
#PBS -l nodes=1:ppn=1

./exspec

gzip -v packets*.out
