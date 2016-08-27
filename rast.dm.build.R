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

# This function will build the full raster design matrix recentering and rescaling it identically to how the full n = 60 design matrix was recentered and rescaled.

# In particular it will:

#    #  1) start with the raw covariate values for all raster pixels and recenter them with the n60.LME.meanx then rescale them with the n60.LME.normx to create LME.rast.RR

#    #  2) create the raster polynomial design matrix (PDM.rast) from the LME.rast.RR design matrix then recenter PDM.rast with n60.PDM.meanx then rescale it with n60.PDM.normx

#    #  2) create the raster, order two interaction design matrix (IO2.rast) from the LME.rast.RR design matrix then recenter IO2.rast with n60.IO2.meanx then rescale it with n60.IO2.normx

# Note: this function supersedes fcp.pred( )

rast.dm.build <- function(LME.rast = B1.CRS.df[,3:ncol(B1.CRS.df)],
                          n60.LME.meanx = n60.P4I2.DM$LME.meanx,
                          n60.LME.normx = n60.P4I2.DM$LME.normx,
                          n60.PDM.meanx = n60.P4I2.DM$PDM.meanx,
                          n60.PDM.normx = n60.P4I2.DM$PDM.normx,
                          n60.IO2.meanx = n60.P4I2.DM$IO2.meanx,
                          n60.IO2.normx = n60.P4I2.DM$IO2.normx){
    # adapt dm.build.R to use provided n60.xxx.meanx and n60.xxx.normx at each step where the design matrix is recentered and rescaled
    # don't to the filtering to enforce a maximum permitted correlation between covarites in this function instead calculate the full design matrix with this function and take subsets of this design matrix as necessary by selecting subsets of covariates by name    
    mpo = 4 # note the facility to alter this value would require editing of polynomial term filtering code below (currently a value of 4 is hard coded in) e.g. |2|3|4        
    # ReCenter all Covariates in Raster Design Matrix with the same transformation that was used to Recenter all covariates the n = 60 design matrix to have mean zero
    LME.rast.RR = as.matrix(LME.rast)
    LME.rast.RR = scale(x = LME.rast.RR, center = n60.LME.meanx, scale = FALSE)
    # ReScale all Covariates in Raster Design Matrix with the same transformation that was used to ReScale all covariates in the n = 60 design matrix to have magnitude 1
    names(n60.LME.normx) = NULL
    LME.rast.RR = scale(x = LME.rast.RR, center = FALSE, scale = n60.LME.normx)
    Data = LME.rast.RR 
    n.LME = ncol(Data)
    P.ng = expand.grid(colnames(Data), as.character(2:mpo))
    colnames(P.ng) = c('Var','P')
    PDM = data.frame(matrix(0, nrow = nrow(Data), ncol = nrow(P.ng))) # Polynomial Terms Dataframe
    P.n = character(nrow(P.ng))
    for(i in 1:nrow(P.ng)){
        P.n[i] = paste(P.ng[i,'Var'], P.ng[i,'P'], sep = '.')}
    colnames(PDM) = P.n
    for(i in 1:nrow(P.ng)){
        PDM[,paste(P.ng[i,'Var'], P.ng[i,'P'], sep = '.')] = (Data[,paste(P.ng[i,'Var'])])^(as.numeric(paste(P.ng[i,'P'])))}
    # ReCenter all Polynomial Terms in Raster Polynomial Term Design Matrix with the same transformation that was used to Recenter all polynomial terms constructed from the n = 60 design matrix to have mean zero    
    PDM.x = as.matrix(PDM)
    PDM.x = scale(PDM.x, center = n60.PDM.meanx, scale = FALSE)
    # ReScale all Polynomial Terms in Raster Polynomial Term Design Matrix with the same transformation that was used to ReScale all polynomial terms constructed from the n = 60 design matrix to have magnitude 1    
    names(n60.PDM.normx) = NULL
    PDM.x = scale(PDM.x, center = FALSE, scale = n60.PDM.normx)    
    PDM = PDM.x
    # Expand Design Matrix to include Columns for Interactions to Order = 1/2 Max Polynomial Order
    # Order 2 Interactions: Products of two distinct linear terms
    IO2.2L.i = t(combn(x = 1:n.LME, m = 2))
    IO2.2L.n = character(nrow(IO2.2L.i))
    IO2.2L.DM = data.frame(matrix(0, nrow = nrow(Data), ncol = length(IO2.2L.n)))
    for(i in 1:nrow(IO2.2L.i)){
        IO2.2L.n[i] = paste(colnames(Data)[IO2.2L.i[i,1]], colnames(Data)[IO2.2L.i[i,2]], sep = '_X_')}
    colnames(IO2.2L.DM) = IO2.2L.n
    for(i in 1:ncol(IO2.2L.DM)){
        cols = strsplit(x = colnames(IO2.2L.DM)[i], split = '_X_')[[1]]
        IO2.2L.DM[,i] = Data[,cols[1]]*Data[,cols[2]]}
    # ReCenter all interaction terms in Raster Interaction term Design Matrix with the same transformation that was used to Recenter all interaction terms constructed from the n = 60 design matrix to have mean zero    
    IO2.2L.DM.x = as.matrix(IO2.2L.DM)    
    IO2.2L.DM.x = scale(IO2.2L.DM.x, center = n60.IO2.meanx, scale = FALSE)
    # ReScale all interaction terms in Raster Interaction term Design Matrix with the same transformation that was used to ReScale all interaction terms constructed from the n = 60 design matrix to have magnitude one
    names(n60.IO2.normx) = NULL    
    IO2.2L.DM.x = scale(IO2.2L.DM.x, center = FALSE, scale = n60.IO2.normx)    
    IO2.2L.DM = IO2.2L.DM.x    
    Rast.P4I2.DM.RR = data.frame(Data, IO2.2L.DM, PDM)    
    output = list(DM = Rast.P4I2.DM.RR)
    return(output)}
