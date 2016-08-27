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


######################################################################
#                                                                    #
#   as sp.dm.build but recenter with n60.colmeans and n60.colnorms   #
#                                                                    #
######################################################################

rast.sp.dm.build <- function(LME = B1.CRS.df[,c('Easting', 'Northing')],
                        mppo = 12,
                        n60.DM = SpEr.DM$DM,     
                        n60.LME.meanx = SpEr.DM$LME.meanx,
                        n60.LME.normx = SpEr.DM$LME.normx,
                        n60.Poly.C.means = SpEr.DM$Poly.C.means,
                        n60.Poly.C.norms = SpEr.DM$Poly.C.norms,
                        n60.PDM.meanx = SpEr.DM$PDM.meanx,
                        n60.PDM.normx = SpEr.DM$PDM.normx,
                        n60.Int.C.meanx = SpEr.DM$Int.C.meanx,
                        n60.Int.C.normx = SpEr.DM$Int.C.normx){    
    # SDMB.IF.PFO = Spatial Design Matrix Builder with Intelligent Filtering 
    # LME = Linear Main Effects design matrix (as dataframe)
    # mpc = Maximum Permissible Correlation between design matrix columns
    # mppo = Maximum Proposed Polynomial Order
    # RCRS = ReCenter each column on 0 and ReScale each column to have magnitude 1 as for lars() function
    # Keep These: tick
    # Center all Covariates on Zero & Scale all Covariates to have Magnitude 1
    LME.x = as.matrix(LME)
    LME.x = scale(LME.x, center = n60.LME.meanx, scale = FALSE)
    names(n60.LME.normx) = NULL
    LME.x = scale(LME.x, center = FALSE, scale = n60.LME.normx)
    LME = LME.x
    Data = data.frame(LME)    
    # Polynomial Proposals
    # don't actually do the filtering just produce the full design matrix as the n60 function would for a mpccm of 1 and take the subset of columns from the n60 design matrix
    for(i in 2:mppo){
        East.Prop.Name = paste('Easting', i, sep = '.')
        if( East.Prop.Name %in% colnames(n60.DM) ){
            Propose = Data[,'Easting']^i
            East.Poly.C.means.i = n60.Poly.C.means[East.Prop.Name] 
            Propose = Propose - East.Poly.C.means.i
            East.Poly.C.norms.i = n60.Poly.C.norms[East.Prop.Name]             
            East.Poly.C.norms.i = sqrt(sum(Propose^2))
            Propose = Propose/East.Poly.C.norms.i
            Data = data.frame(Data,Propose)                        
            colnames(Data) = c(colnames(Data)[1:(ncol(Data)-1)], paste('Easting', i, sep = '.'))}
        North.Prop.Name = paste('Northing', i, sep = '.')
        if( North.Prop.Name %in% colnames(n60.DM) ){
            Propose = Data[,'Northing']^i
            North.Poly.C.means.i = n60.Poly.C.means[North.Prop.Name] 
            Propose = Propose - North.Poly.C.means.i
            North.Poly.C.norms.i = n60.Poly.C.norms[North.Prop.Name]             
            North.Poly.C.norms.i = sqrt(sum(Propose^2))
            Propose = Propose/North.Poly.C.norms.i
            Data = data.frame(Data,Propose)                        
            colnames(Data) = c(colnames(Data)[1:(ncol(Data)-1)], paste('Northing', i, sep = '.'))}}
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
    # ReCenter and ReScale Covariate terms with the same transformations that were used to center the n60 covariate terms on Zero & Scale the n60 Covariates to have Magnitude 1
    ##
    PDM.x = as.matrix(PDM)
    PDM.x = scale(PDM.x, center = n60.PDM.meanx, scale = FALSE)
    names(n60.PDM.normx) = NULL
    PDM.x = scale(PDM.x, center = FALSE, scale = n60.PDM.normx)
    PDM = PDM.x
    PPDM = data.frame(LME, PDM)
    P.ng.acron = data.frame(P.ng, Acron = paste(P.ng[,'Var'],P.ng[,'P'],sep='.'))
    P.ng.acron = rbind(data.frame(Var = c('Easting', 'Northing'), P = c(1,1), Acron = c('Easting', 'Northing')), P.ng.acron)    
    AI.2 = combn(x = 1:ncol(PPDM), m = 2)
    for(i in 1:ncol(AI.2)){
        AC = data.frame(Acron = factor(colnames(PPDM[,AI.2[,i]]), levels = levels(P.ng.acron$Acron)), PO = numeric(length(colnames(PPDM[,AI.2[,i]]))))
        for(k in 1:nrow(AC)){
            AC[k, 'PO'] = P.ng.acron[P.ng.acron$Acron == AC[k, 'Acron'], 'P']}
        if(sum(AC$PO) <= floor(mppo/2)){
            Propose.name = paste(colnames(PPDM[,AI.2[,i]]), collapse = '_X_')
            if(Propose.name %in% colnames(n60.DM)){
                Propose = apply(X = PPDM[,AI.2[,i]], MARGIN = 1, FUN = prod)
                Int.C.meanx.i = n60.Int.C.meanx[Propose.name]
                Propose = Propose - Int.C.meanx.i
                Int.C.normx.i = n60.Int.C.normx[Propose.name]
                Propose = Propose/Int.C.normx.i
                Data = data.frame(Data,Propose)
                colnames(Data) = c(colnames(Data)[1:(ncol(Data)-1)], paste(colnames(PPDM[,AI.2[,i]]), collapse = '_X_'))}}}
    Output = list(DM = Data)
    return(Output)}
