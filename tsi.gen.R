# R Code for one of the functions provided in the
# Least Angle Regression Companion (LARC) project:
# https://github.com/brfitzpatrick/larc.git 

# Copyright (C) 2016 Benjamin R. Fitzpatrick, David W. Lamb and Kerrie Mengersen.

# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along with this program, in the form a text file titled 'LICENSE.
# If not, see http://www.gnu.org/licenses/.

# The program author, Benjamin R. Fitzpatrick, may be contacted via email at <ben.r.fitzpatrick@gmail.com>

# This function is used in the analysis discussed in:

# `Ultrahigh Dimensional Variable Selection for Interpolation of Point Referenced Spatial Data: A Digital Soil Mapping Case Study' currently in press at the Journal PLOS ONE 10.1371/journal.pone.0162489

tsi.gen <- function(n = 10, all = 60, train = 35){
  require('randtoolbox')  
  tsg = function(){
  qunt = quantile(x = 1:all, probs = seq(from = 0, to = 1, length.out = 11)) # 10 intervals # so pick one validation point from each DivBy interval
  A60 = ceiling(all*SFMT(n = 200, dim = 1, mexp = 19937, usepset = TRUE, withtorus = 0.5, usetime = TRUE))
  B60 = c(A60[150:200], A60[1:50], A60[100:149], A60[51:99])
  U60 = unique(B60)
  #U60 = unique(ceiling(60*WELL(n=200, dim = 1, order = 44497, temper = TRUE, version = "a")))
    if(length(U60) < all){
      while(length(U60)<all){
        U60 = unique(c(U60, ceiling(all*SFMT(n = 200, dim = 1, mexp = 19937, usepset = TRUE, withtorus = 0.5, usetime = TRUE))))}}
    U60 = U60[1:all]
    tdf = data.frame(DivBy = 1:all, DataInd = U60) # jumbling DataInd so that partitioning by DivBy doesn't bias things
    # ceiling(qunt)
    #  0%  10%  20%  30%  40%  50%  60%  70%  80%  90% 100% 
    #  1    7   13   19   25   31   37   43   49   55   60
    n.drop = all-train
    drop.index = numeric()
    while(length(unique(drop.index)) < n.drop){
        drop.index = unique(c(drop.index,ceiling(n.drop* WELL(n = n.drop, dim = 1, order = 44497, temper = TRUE, version = "a")) + ceiling(seq(from = 0, to = train, length.out = n.drop))))}
    drop.index = drop.index[1:n.drop]
    return(tdf[-drop.index, 'DataInd'])}  
  rand.start = proc.time()
  n.val = train
  index = 1:all
  CV.VI <- matrix(0, nrow = n.val, ncol = n) # CV.Train.Indices
  i = 1
  while(i < n+1){ 
    cvi = sort(tsg())
    if(i == 1){
      CV.VI[ ,i] <- cvi
      print(paste(i))
      i <- i+1} else{
        if(i == 2){
          if(!0 == sum(CV.VI[ , 1:(i-1)] - cvi)){ # ensuring uniqueness
            CV.VI[ , i] <- cvi
            time.to.here = round((proc.time()[3] - rand.start[3])/60,2)
            print(c(paste(100*i/n),'%',paste(time.to.here),'min'),sep = ' ', quote = FALSE)
            i <- i + 1}} else{
          if(!0 %in% colSums(CV.VI[ , 1:(i-1)] - cvi)){ # ensuring uniqueness 
            CV.VI[ , i] <- cvi
            time.to.here = round((proc.time()[3] - rand.start[3])/60,2)
            print(c(paste(round(100*i/n,4)),'%',paste(time.to.here),'min'),sep = ' ', quote = FALSE)
            i <- i + 1}}}}
  return(CV.VI)} 
