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

# this function replicates the function 'lar.cv.R' but substitutes LAR variable selection for variable selection with one of: exhaustive search, backward stepwise selection, forwards stepwise selection or sequential replacement (forward/backwards stepwise) variable selection for comparison purposes
# these variable selection algorithms are implemented through the R package `leaps'

lm.vs.cv <- function(Response = Response.df$Carbon,
                   DM = n60.P4I2.DM.F40$DM,
                   CVTI = CVTI_TS35_VS25_n500,
                   vs_type = c('exhaustive','backward', 'forward', 'seqrep')){
    require('leaps')
    n = ncol(CVTI)
    Coef.ls = vector(mode = 'list', length = n)
    ideal.ss.size = vector(mode = 'numeric', length = n)
    LL.VPSS.ls = vector(mode = 'list', length = n)
    vspe.mat = matrix(ncol = n, nrow = (nrow(DM)-nrow(CVTI)))
    start = proc.time()
    DM.TI.C.means.store = matrix(data = NA, nrow = n, ncol = ncol(DM))
    DM.TI.C.norms.store = matrix(data = NA, nrow = n, ncol = ncol(DM))
    Response.TI.means    = vector(mode = 'numeric', length = n)  
    for(i in 1:n){
        # ReCenter & ReScale Training Set Design Matrix
        Response.TI.means[i] = mean(Response[CVTI[ , i]])
        Response.TI.RC = (Response[CVTI[ , i]] - Response.TI.means[i])
        DM.TI = as.matrix(DM[ CVTI[ , i], ])
        DM.TI.C.means.store[i,] = colMeans(DM.TI)
        DM.TI.RR = scale(DM.TI, center = DM.TI.C.means.store[i,], scale = FALSE)
        DM.TI.C.norms.store[i,] = sqrt(colSums(DM.TI.RR^2))        
        DM.TI.RR = scale(DM.TI.RR, center = FALSE, scale = DM.TI.C.norms.store[i, ]) # each column divided by matching entry in normx        
        LL.i = regsubsets(x = DM.TI.RR,
                          y = Response.TI.RC,
                          method = vs_type,
                          nvmax=ncol(DM.TI.RR),
                          intercept = FALSE)
        max.steps = (length(summary(LL.i)$obj) + 1)
        vspe.ss = vector(mode = 'numeric', length = max.steps)        
        vspe.mat.i = matrix(nrow = (nrow(DM)-nrow(CVTI)), ncol = max.steps)
        # ReCenter & ReScale Validation Set Design Matrix with the same transformations as were applied to the associated Training Set Design Matrix
        DM.VI = as.matrix(DM[ -CVTI[ , i], ])
        DM.VI.RR = scale(DM.VI, center = DM.TI.C.means.store[i,], scale = FALSE)
        DM.VI.RR = scale(DM.VI.RR, center = FALSE, scale = DM.TI.C.norms.store[i, ])        
        for(j in 1:max.steps){
            if(j == 1){
                # intercept only model
                LL.i.pred = Response.TI.means[i]} else{
                    # use covariate(s)
                    Coef.i.M = coef(object = LL.i, id = (j-1))
                    DM.VI.RR.i.sc = as.matrix(DM.VI.RR[, names(Coef.i.M)])
                    LL.i.pred = (( as.matrix(DM.VI.RR.i.sc) %*% matrix(data = Coef.i.M, ncol = 1)) + Response.TI.means[i])}
                    #
                    #  MA.Pred[1:60,1] = Pred[1:60, From.CV1 From.CV2 ... From.CVn] %*% | Weight Model 1 |
                    #                                                                   | Weight Model 2 |
                    #                                                                   |        .       |
                    #                                                                   |        .       |
                    #                                                                   |        .       |
                    #                                                                   | Weight Model n |
            VSPE = (Response[-CVTI[ , i]] - (LL.i.pred))
            vspe.ss[j]   = sum((VSPE)^2)
            vspe.mat.i[,j] = VSPE}
        min.vspe.step = min(which( vspe.ss == min(vspe.ss))) # this is j value of min vsep.ss element
        vspe.mat[,i] = vspe.mat.i[,min.vspe.step]
        if(min.vspe.step == 1){
            Coef = 0
            ideal.ss.size[i] = 0}
        if(min.vspe.step > 1){
            Coef =  coef(object = LL.i, id = (min.vspe.step-1))
            ideal.ss.size[i] = length(names(Coef))}
        Coef.ls[[i]] = Coef        
        LL.VPSS.ls[[i]] = vspe.ss
        here = proc.time()
        print(paste(round(100*i/n,1), '% Complete, Predict', round((((here[3] - start[3])/i)*(n-i))/60,2),'min remaining'))}
    if(NA %in% ideal.ss.size){return(print('NA in Ideal Subset Size'))} else{
        if(NA %in% vspe.ss){return(print('NA in Validation Set Prediction Error Matrix'))} else{
            output.ls = list(Design.mat = DM, VSPE.mat = vspe.mat, SSS.num = ideal.ss.size, SC.ls = Coef.ls, Response = Response, CV.TI = CVTI, Response.TI.means = Response.TI.means, DM.TI.C.means.store = DM.TI.C.means.store, DM.TI.C.norms.store = DM.TI.C.norms.store)
            return(output.ls)}}}
