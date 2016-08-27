# R Code for one of the functions provided in the
# Least Angle Regression Companion (LARC) project:
# https://github.com/brfitzpatrick/larc.git 

# Copyright (C) 2016 Benjamin R. Fitzpatrick, David W. Lamb and Kerrie Mengersen.

# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along with this program, in the form a text file titled 'LICENSE.
# If not, see http://www.gnu.org/licenses/.

# The program author, Benjamin R. Fitzpatrick, may be contacted via email at ben.r.fitzpatrick@gmail.com

# This function is used in the analysis discussed in:

# `Ultrahigh Dimensional Variable Selection for Interpolation of Point Referenced Spatial Data: A Digital Soil Mapping Case Study' currently in press at the Journal PLOS ONE 10.1371/journal.pone.0162489

sp.dm.build <- function(LME = Response.df[,c('Easting', 'Northing')], mpc = 0.95, mppo = 12){
    # SDMB.IF.PFO = Spatial Design Matrix Builder with Intelligent Filtering 
    # LME = Linear Main Effects design matrix (as dataframe)
    # mpc = Maximum Permissible Correlation between design matrix columns
    # mppo = Maximum Proposed Polynomial Order
    # RCRS = ReCenter each column on 0 and ReScale each column to have magnitude 1 as for lars() function
    # Keep These: tick
    # Center all Covariates on Zero & Scale all Covariates to have Magnitude 1
    LME.x = as.matrix(LME)
    LME.meanx = colMeans(LME.x)    
    LME.x = scale(LME.x, center = LME.meanx, scale = FALSE)
    normx = sqrt(colSums(LME.x^2))
    LME.normx = sqrt(colSums(LME.x^2))    
    names(normx) = NULL
    LME.x = scale(LME.x, center = FALSE, scale = normx)
    if(max(abs(colMeans(LME.x))) < 1e-9 & max(abs(1 - colSums(LME.x^2))) < 1e-9){
        LME = LME.x} else{print('RCRS Error')}    
    # Filtering Linear Main Effects Design Matrix so remainder has no pairs of covariates with correlatation coefficients magnitude greater than mpc
    Data = data.frame(LME)    
    # Polynomial Proposals
    East.Poly.C.means = numeric()
    East.Poly.C.norms = numeric()
    North.Poly.C.means = numeric()
    North.Poly.C.norms = numeric()    
    for(i in 2:mppo){
        Propose = Data[,'Easting']^i
        # Keep these
        East.Poly.C.means.i = mean(Propose)
        Propose = Propose - East.Poly.C.means.i
        East.Poly.C.norms.i = sqrt(sum(Propose^2))
        Propose = Propose/East.Poly.C.norms.i        
        if((TRUE %in% (abs(cor(Propose,Data)) > mpc)) == FALSE){
            East.Poly.C.means <- c(East.Poly.C.means, East.Poly.C.means.i)
            East.Poly.C.norms <- c(East.Poly.C.norms, East.Poly.C.norms.i)            
            Data = data.frame(Data,Propose)
            colnames(Data) = c(colnames(Data)[1:(ncol(Data)-1)], paste('Easting', i, sep = '.'))
            names(East.Poly.C.means)[East.Poly.C.means == East.Poly.C.means.i] <- paste('Easting', i, sep = '.')
            names(East.Poly.C.norms)[East.Poly.C.norms == East.Poly.C.norms.i] <- paste('Easting', i, sep = '.') }
        Propose = Data[,'Northing']^i
        # Keep these
        North.Poly.C.means.i = mean(Propose)        
        Propose = Propose - North.Poly.C.means.i
        North.Poly.C.norms.i = sqrt(sum(Propose^2))
        Propose = Propose/North.Poly.C.norms.i        
        if((TRUE %in% (abs(cor(Propose,Data)) > mpc)) == FALSE){
            North.Poly.C.means <- c(North.Poly.C.means, North.Poly.C.means.i)
            North.Poly.C.norms <- c(North.Poly.C.norms, North.Poly.C.norms.i)            
            Data = data.frame(Data,Propose)
            colnames(Data) = c(colnames(Data)[1:(ncol(Data)-1)], paste('Northing', i, sep = '.'))
            names(North.Poly.C.means)[North.Poly.C.means == North.Poly.C.means.i] <- paste('Northing', i, sep = '.')
            names(North.Poly.C.norms)[North.Poly.C.norms == North.Poly.C.norms.i] <- paste('Northing', i, sep = '.') } }
    # Interaction Proposals:
    # # Note: for interaction of polynomial terms of two variables
    #         all possible 3rd and higher order interactoins can be
    #         expressed as a 2nd order (pairwise)  interaction
    #         of polynomial terms of these two variables
    P.ng = expand.grid(colnames(LME), (2:(mppo/2-1)))
    colnames(P.ng) = c('Var','P')
    PDM = data.frame(matrix(0, nrow = nrow(LME), ncol = nrow(P.ng)))
    P.n = character(nrow(P.ng))
    for(i in 1:nrow(P.ng)){P.n[i] = paste(P.ng[i,'Var'], P.ng[i,'P'], sep = '.')}
    colnames(PDM) = P.n
    for(i in 1:nrow(P.ng)){PDM[,paste(P.ng[i,'Var'], P.ng[i,'P'], sep = '.')] = (LME[,paste(P.ng[i,'Var'])])^(P.ng[i,'P'])}    
    # Center all Covariates on Zero & Scale all Covariates to have Magnitude 1
    PDM.x = as.matrix(PDM)
    PDM.meanx = colMeans(PDM.x)
    PDM.x = scale(PDM.x, center = PDM.meanx, scale = FALSE)
    PDM.normx = sqrt(colSums(PDM.x^2))
    normx = sqrt(colSums(PDM.x^2))
    names(normx) = NULL
    PDM.x = scale(PDM.x, center = FALSE, scale = normx)
    if(max(abs(colMeans(PDM.x))) < 1e-9 & max(abs(1 - colSums(PDM.x^2))) < 1e-9){
        PDM = PDM.x} else{print('RCRS Error')}
    PPDM = data.frame(LME, PDM)    
    P.ng.acron = data.frame(P.ng, Acron = paste(P.ng[,'Var'],P.ng[,'P'],sep='.'))
    P.ng.acron = rbind(data.frame(Var = c('Easting', 'Northing'), P = c(1,1), Acron = c('Easting', 'Northing')), P.ng.acron)    
    AI.2 = combn(x = 1:ncol(PPDM), m = 2)
    Int.C.meanx = vector(mode = 'numeric')
    Int.C.normx = vector(mode = 'numeric')    
    for(i in 1:ncol(AI.2)){
        AC = data.frame(Acron = factor(colnames(PPDM[,AI.2[,i]]), levels = levels(P.ng.acron$Acron)), PO = numeric(length(colnames(PPDM[,AI.2[,i]]))))
        for(k in 1:nrow(AC)){
            AC[k, 'PO'] = P.ng.acron[P.ng.acron$Acron == AC[k, 'Acron'], 'P']}
        if(sum(AC$PO) <= floor(mppo/2)){
            Propose = apply(X = PPDM[,AI.2[,i]], MARGIN = 1, FUN = prod)
            Int.C.meanx.i = mean(Propose)
            Propose = Propose - Int.C.meanx.i
            Int.C.normx.i = sqrt(sum(Propose^2))            
            Propose = Propose/Int.C.normx.i
            if((TRUE %in% (abs(cor(Propose,Data)) > mpc)) == FALSE){
                Data = data.frame(Data,Propose)
                colnames(Data) = c(colnames(Data)[1:(ncol(Data)-1)], paste(colnames(PPDM[,AI.2[,i]]), collapse = '_X_'))                
                Int.C.meanx <- c(Int.C.meanx, Int.C.meanx.i) 
                names(Int.C.meanx)[Int.C.meanx == Int.C.meanx.i] <- paste(colnames(PPDM[,AI.2[,i]]), collapse = '_X_')                
                Int.C.normx <- c(Int.C.normx, Int.C.normx.i)
                names(Int.C.normx)[Int.C.normx == Int.C.normx.i] <- paste(colnames(PPDM[,AI.2[,i]]), collapse = '_X_') }
        }
    }    
    Output = list(DM = Data,
                  LME.meanx = LME.meanx,
                  LME.normx = LME.normx,
                  Poly.C.means = c(East.Poly.C.means, North.Poly.C.means), # use in contructing single term polynomials
                  Poly.C.norms = c(East.Poly.C.norms, North.Poly.C.norms),        
                  PDM.meanx = PDM.meanx, # used in constructing single term polynomials for use in constructing the interaciton terms
                  PDM.normx = PDM.normx,
                  Int.C.meanx = Int.C.meanx,
                  Int.C.normx = Int.C.normx)
    return(Output)}
