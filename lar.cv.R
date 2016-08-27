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

lar.cv <- function(Response = Response.df$Carbon,
                   DM = n60.P4I2.DM$DM,
                   CVTI = CVTI_TS35_VS25_n500){    
    require('lars')
    n = ncol(CVTI)
    Coef.ls = vector(mode = 'list', length = n)
    ideal.ss.size = vector(mode = 'numeric', length = n)
    # LL.ls = list()
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
        # check:
        # summary(abs(colMeans(DM.TI.RR)) < 1e-10)
        # and check:
        # summary(round((colSums(DM.TI.RR^2)),10) == 1)
        LL.i = lars(x = DM.TI.RR,
                    y = Response.TI.RC,
                    type = 'lasso',
                    trace = FALSE,
                    normalize = FALSE,
                    intercept = FALSE,
                    eps = .Machine$double.eps,
                    max.steps = 2.5e5,
                    use.Gram = TRUE) # one execution requires 4.3 Gb of RAM (in addition to the amount the OS is using to run and R is already using) for a DM.TI.RR of 35 rows by 2205 columns
        # LL.ls[[i]] = LL.i # Keeping this list in RAM for n = 500 consumes more RAM than I have...
        max.steps = nrow(summary(LL.i))
        vspe.ss = vector(mode = 'numeric', length = max.steps) 
        # ReCenter & ReScale Validation Set Design Matrix with the same transformations as were applied to the associated Training Set Design Matrix        
        DM.VI = as.matrix(DM[ -CVTI[ , i], ])
        DM.VI.RR = scale(DM.VI, center = DM.TI.C.means.store[i,], scale = FALSE)
        DM.VI.RR = scale(DM.VI.RR, center = FALSE, scale = DM.TI.C.norms.store[i, ])
        # vspe.ls = vector(mode = 'list', length = max.steps)
        vspe.mat.i = matrix(nrow = nrow(DM.VI), ncol = max.steps)
        for(j in 1:max.steps){ # j = 1 is intercept only model
            # LL.i.pred.0 = (predict(object =  LL.i,
            #                    type = c("fit"),
            #                    newx = DM.VI.RR,
            #                    s = j,
            #                    mode = c("step"))$fit)  + Response.TI.means[i]
            LL.i.pred.coef = (predict(object =  LL.i, type = 'coefficients', s = j, mode = 'step'))
            LL.i.pred.coef.val = LL.i.pred.coef$coefficients[abs(LL.i.pred.coef$coefficients) > 0]
            LL.i.pred = (DM.VI.RR[,names(LL.i.pred.coef.val)] %*% matrix(data = LL.i.pred.coef.val, ncol = 1) ) + Response.TI.means[i]
            # summary(LL.i.pred.0 - LL.i.pred  ) # check
            VSPE = (Response[-CVTI[ , i]] - (LL.i.pred))
            vspe.ss[j]   = sum((VSPE)^2)
            # vspe.ls[[j]] = VSPE
            vspe.mat.i[,j] = VSPE}
        min.vspe.step = min(which( vspe.ss == min(vspe.ss))) # the j value of the min vspe.ss
        # vspe.mat[,i] = vspe.ls[[min.vspe.step]]
        vspe.mat[,i] = vspe.mat.i[,min.vspe.step]
        Coef = predict(object =  LL.i, type = 'coefficients', s = min.vspe.step, mode = 'step')
        Val.Coef.ls = (Coef$coefficients[abs(Coef$coefficients) > 0])
        Coef.ls[[i]] = Val.Coef.ls
        LL.VPSS.ls[[i]] = vspe.ss
        ideal.ss.size[i] = length(names(Val.Coef.ls))
        here = proc.time()
        print(paste(round(100*i/n,1), '% Complete, Predict', round((((here[3] - start[3])/i)*(n-i))/60,2),'min remaining'))}
    if(NA %in% ideal.ss.size){return(print('NA in Ideal Subset Size'))} else{
        if(NA %in% vspe.ss){return(print('NA in Validation Set Prediction Error Matrix'))} else{
            output.ls = list(Design.mat = DM, VSPE.mat = vspe.mat, SSS.num = ideal.ss.size, SC.ls = Coef.ls, Response = Response, CV.TI = CVTI, Response.TI.means = Response.TI.means, DM.TI.C.means.store = DM.TI.C.means.store, DM.TI.C.norms.store = DM.TI.C.norms.store)
            return(output.ls)}}} 
