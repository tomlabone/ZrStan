#!/bin/bash

#Tom LaBone
#September 2, 2021

Zr/zircon-ind \
   sample num_warmup=4000 num_samples=15000 \
   data file=Zr/Zr_10.data.R \
   init=Zr/Zr_10.init.R \
   output file=Zr/Zr_10.csv
   
Zr/zircon-ind \
  optimize \
  data file=Zr/Zr_10.data.R \
  init=Zr/Zr_10.init.R \
  output file=Zr/Zr_10opt.csv

