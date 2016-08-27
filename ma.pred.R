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

# Calculates Model Averaged Predictions from Variable Selection Results from a Cross Validatation Scheme
# will work for both predicting at the locations where the response was observed (and all covariates are available) and at other locations where the all covariates are available but the response may or may not be observed
# note to calculate an R2 you need response observations to accompany the predictions so you'll only be able to calculate an R2 for predictions at the locations at which the response was observed

# now supersedes arma.pred.R 

ma.pred <- function(cotpa.mat = LAR.F40.TS35.n500$Design.mat,
                    vspe.mat = LAR.F40.TS35.n500$VSPE.mat,
                    sc.ls = LAR.F40.TS35.n500$SC.ls,
                    Response = LAR.F40.TS35.n500$Response,
                    CV.TI = LAR.F40.TS35.n500$CV.TI,
                    Response.TI.means = LAR.F40.TS35.n500$Response.TI.means,
                    DM.TI.C.means.store = LAR.F40.TS35.n500$DM.TI.C.means.store,
                    DM.TI.C.norms.store = LAR.F40.TS35.n500$DM.TI.C.norms.store,
                    pred.disp.est = c('variance','interval.width'),
                    interval = 0.95,
                    plot.wgf = TRUE,
                    calc.R.2 = TRUE){    
    # cotpa.mat = matrix of Covariate Observations To Predict At (i.e. at which to calculate the model averaged predictions)
    # vspe.mat = matrix of Validation Set Prediction Errors (columns are each the validation set prediction errors from a LAR run on a different training set
    # sc.ls = list of Selected Covariates from the selected LAR run on each training set from the cross validation scheme
    require(ggplot2)
    n = ncol(vspe.mat)
    vspess = colSums(vspe.mat^2)
    weight = (1/vspess)/sum((1/vspess))
    # sum(weight) == 1 # check
    Weight.df = data.frame(vspess, weight)
    W.p = ggplot(aes(x = vspess, y = weight), data = Weight.df) + geom_point(size = 2, alpha = 0.25) + labs(x = 'Validation Set Prediction Error Sum of Squares', y = 'Model Weighting for Model Averaging', title = paste('Weightings of Each of the', n, 'Selected Models from the Cross Validation Scheme'))
        #
        # (Pred.i = as.matrix(cotpa.mat[, names(sc.ls[[1]])]) %*%  | sc.ls[[1]][1] | ) + intercept
        #                                                          | sc.ls[[1]][2] |
        #                                                          |      ...      |
        #                                                          | sc.ls[[1]][p] |
        #
        # summary(matrix(data = sc.ls[[3]], byrow = FALSE, ncol = 1) == matrix(data = sc.ls[[3]], byrow = TRUE, ncol = 1)) # check: byrow makes no difference when it's a 1 row or 1 column matrix
        Pred = matrix(nrow = nrow(cotpa.mat), ncol = length(sc.ls))
        for(i in 1:length(sc.ls)){
            DM.RR.i = cotpa.mat
            # recenter and rescale design matrix being predicted from with the same transformations that were applied to the training set to which the LAR algorithm was applied that produced the model we are using to predict at this iteration
            # check:
            # summary((abs(colMeans(scale(cotpa.mat[CV.TI[,i],] , center = DM.TI.C.means.store[i,],scale = FALSE))))< 1e-10)
            if(length(sc.ls[[i]]) == 1){
                if(sc.ls[[i]] == 0){
                    Pred[,i] = Response.TI.means[i]} else{
                        if(FALSE %in% (names(sc.ls[[i]])) %in% colnames(cotpa.mat)){return(print(x = 'Error: Covariate Names in sc.ls not found in colnames of cotpa.mat'))} else{
                            DM.RR.i = scale(DM.RR.i, center = DM.TI.C.means.store[i,], scale = FALSE)
                            DM.RR.i = scale(DM.RR.i, center = FALSE, scale = DM.TI.C.norms.store[i, ])
                            DM.RR.i.sc = as.matrix(DM.RR.i[, names(sc.ls[[i]])])
                            Coef.i.M = matrix(data = sc.ls[[i]], byrow = FALSE, ncol = 1)
                            Pred[,i] = (( as.matrix(DM.RR.i.sc) %*% Coef.i.M ) + Response.TI.means[i])}}} else{
                                if(FALSE %in% (names(sc.ls[[i]]) %in% colnames(cotpa.mat))){return(print(x = 'Error: Covariate Names in sc.ls not found in colnames of cotpa.mat'))} else{
                            DM.RR.i = scale(DM.RR.i, center = DM.TI.C.means.store[i,], scale = FALSE)
                            DM.RR.i = scale(DM.RR.i, center = FALSE, scale = DM.TI.C.norms.store[i, ])
                            DM.RR.i.sc = as.matrix(DM.RR.i[, names(sc.ls[[i]])])
                            Coef.i.M = matrix(data = sc.ls[[i]], byrow = FALSE, ncol = 1)
                            Pred[,i] = (( as.matrix(DM.RR.i.sc) %*% Coef.i.M ) + Response.TI.means[i])}}}
        #
        #  MA.Pred[1:60,1] = Pred[1:60, From.CV1 From.CV2 ... From.CVn] %*% | Weight Model 1 |
        #                                                                   | Weight Model 2 |
        #                                                                   |        .       |
        #                                                                   |        .       |
        #                                                                   |        .       |
        #                                                                   | Weight Model n |
        MA.Pred = Pred %*% matrix(data = weight, nrow = n, ncol = 1, byrow = TRUE)
        Results = data.frame(MA.Pred)
        if('interval.width' %in% pred.disp.est){
            if( (interval > 0) & (interval < 1)){
                Lwr = (1 - interval)/2
                Upr = 1 - Lwr
                Results$Int.Width = numeric(length = nrow(Results))
                for(i in 1:nrow(Pred)){
                    Upper.i = quantile(x = Pred[i,], probs = Upr)
                    Lower.i = quantile(x = Pred[i,], probs = Lwr)
                    Results[i,'Int.Width'] = (Upper.i - Lower.i)}} else{ print(c('Please specify the percentage of the data you would like contained in the interval as a decimal'))}}
        if('variance' %in% pred.disp.est){
            Results$Var = numeric(length = nrow(Results))
            for(i in 1:nrow(Pred)){
                Results[i,'Var'] = var(Pred[i,])}}        
        output.ls = list(Results = Results)
        if(calc.R.2 == TRUE){
           if( length(Response) == length(MA.Pred) ){
               MA.Pred.Residual = Response - MA.Pred
               # Coefficient of Determination = 1 - (Residual SS/Total SS)
               # Residual Sum of Squares = sum((Response - Prediction)^2)
               # Total Sum of Squares = sum((Response - mean(Response))^2)
               R2 = 1 - (sum((Response - MA.Pred)^2)/sum((Response-mean(Response))^2))
               output.ls$R2 = R2} else{ R2 = c('You can not calculate the coefficient of determination for preditions of the response at locations other than those at which the response was observed. Note: the model averaged predictions have still been calculated at the locations requested.')}}
        if(plot.wgf == TRUE){output.ls$plot = W.p}
        return(output.ls)}
