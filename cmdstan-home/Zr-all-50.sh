#!/bin/bash

#Tom LaBone
#September 4, 2021


echo -e '\n submitted a job'
echo 'hostname is'
hostname

DIR=

# Define a timestamp function
timestamp() {
  date +"%T" # current time
}

timestamp 

echo Zirconium 

Zr/zircon-all \
   sample num_warmup=3000 num_samples=10000 \
   data file=Zr/Zr612-50.data.R \
   init=Zr/Zr612-50.init.R \
   output file=Zr/Zr612-50.csv
  

timestamp 
