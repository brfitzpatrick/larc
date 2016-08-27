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

DbC.PSA <- function(LME = B1.DM.LME[,3:ncol(B1.DM.LME)], Active.Columns = colnames(B1.LME[,3:ncol(B1.LME)]), mpc = 0.6, Prefer.Set = c('ECA','NIR','RED', 'SR', 'DVI', 'NDVI', 'SAVI', 'NLVI', 'MNLVI', 'MSR', 'TVI', 'RDVI')){
    # Drop by Correlation Prefer Set A
    # Examines all pairs of columns in LME
    # for those with correlation coefficient magnitudes greater than mpc
    # retains any with Prefer value in column name and discards other member of pair
    # Note: No rescaling or recentering in this function 
    CC = data.frame(t(combn(x = Active.Columns, m = 2)))
    colnames(CC) = c('Var1','Var2')
    CC$Cor = numeric(nrow(CC))
    for(i in 1:nrow(CC)){
        CC[i, 'Cor'] = cor(LME[,paste(CC[i,'Var1'])], LME[,paste(CC[i,'Var2'])])}
    CC.vc = CC[abs(CC$Cor) >= mpc,]
    # From Correlated Pairs choose:
    # Keep over any other covariate since no other covariates are derrived from Keep
    Covariates.drop = Active.Columns
    Keep = rep(NA, nrow(CC.vc))
    Drop = rep(NA, nrow(CC.vc))
    for(i in 1:nrow(CC.vc)){
        test.pair = c(paste(CC.vc[i, 'Var1']), paste(CC.vc[i, 'Var2']))
        if((TRUE %in% (Prefer.Set %in% unlist(strsplit(x = test.pair[1], split = '[.]')))) &
           (TRUE %in% (Prefer.Set %in% unlist(strsplit(x = test.pair[2], split = '[.]'))))){
            Drop[i] = paste(CC.vc[i, 'Var2'])
            Keep[i] = paste(CC.vc[i, 'Var1'])} else{
                if(TRUE %in% (Prefer.Set %in% unlist(strsplit(x = test.pair, split = '[.]')))){
                    if(TRUE %in% (Prefer.Set %in% unlist(strsplit(x = test.pair[1], split = '[.]')))){
                        Keep[i] = paste(test.pair[1])
                        Drop[i] = paste(test.pair[2])} else{
                            Keep[i] = paste(test.pair[2])
                            Drop[i] = paste(test.pair[1])}}}}
    CC.Keep.Filter = data.frame(Keep = factor(x = Keep, levels = c(NA,Active.Columns)), Drop = factor(x = Drop, levels = c(NA,Active.Columns)))
    CC.Keep.Filter = na.omit(CC.Keep.Filter)
    levels(CC.Keep.Filter$Keep) = c(levels(CC.Keep.Filter$Keep), levels(CC.Keep.Filter$Drop))
    levels(CC.Keep.Filter$Drop) = c(levels(CC.Keep.Filter$Drop), levels(CC.Keep.Filter$Keep))
    if(nrow(CC.Keep.Filter)>0){
    for(i in 2:nrow(CC.Keep.Filter)){
        if(CC.Keep.Filter[i, 'Keep'] %in% CC.Keep.Filter[1:(i-1), 'Drop']){
            Switch = c(paste(CC.Keep.Filter[i, 'Keep']), paste(CC.Keep.Filter[i, 'Drop']))
            CC.Keep.Filter[i, 'Keep'] = Switch[2]
            CC.Keep.Filter[i, 'Drop'] = Switch[1]}}
    Covariates.drop1 = unique(CC.Keep.Filter$Drop)
    Covariates.keep1 = Active.Columns[!(Active.Columns %in% Covariates.drop1)]} else{
    Covariates.keep1 = Active.Columns
    Covariates.drop1 = NULL}        
    Output = list(Covariates.keep1,Covariates.drop1)
    names(Output) = c('Keep','Drop')
    return(Output)}

DMB.IF <- function(LME = B1.DM.LME[,3:ncol(B1.DM.LME)], mpccm = 0.6, RCRS = TRUE){    
    # DMB.IF = Design Matrix Builder with Intelligent Filtering
    # LME = Linear Main Effects design matrix (as dataframe)
    # mpo = Maximum Polynomial Order
    # (this function will only interactions of linear main effects)
    # mpccm = Maximum Permissible Correlation between design matrix columns
    # RCRS = ReCenter each column on 0 and ReScale each column to have magnitude 1 as for lars() function
    mpo = 4 # note the facility to alter this value would require editing of polynomial term filtering code below (currently a value of 4 is hard coded in) e.g. |2|3|4
    if(RCRS == TRUE){
        # Center all Covariates on Zero & Scale all Covariates to have Magnitude 1
        # Keep these: # tick
        x = as.matrix(LME)
        LME.meanx = colMeans(x)
        x = scale(x, center = LME.meanx, scale = FALSE)
        normx = sqrt(colSums(x^2))
        LME.normx = sqrt(colSums(x^2))
        names(normx) = NULL
        x = scale(x, center = FALSE, scale = normx)          
        if(unique(round(colMeans(x),10)) == 0 & unique(round(colSums(x^2),10)) == 1){
            LME = x} else{print('RCRS Error')}}
    # Filtering Linear Main Effects Design Matrix so remainder has no pairs of covariates with correlatation coefficients magnitude greater than mpccm
    LME1 = DbC.PSA(LME = LME, Active.Columns = colnames(LME), mpc = mpccm, Prefer.Set = c('ECA'))
    LME2 = DbC.PSA(LME = LME, Active.Columns = LME1$Keep, mpc = mpccm, Prefer.Set = c('SR', 'DVI', 'NDVI', 'SAVI', 'NLVI', 'MNLVI', 'MSR', 'TVI', 'RDVI'))
    LME3 = DbC.PSA(LME = LME, Active.Columns = LME2$Keep, mpc = mpccm, Prefer.Set = c('RED', 'NIR'))
    LME4 = DbC.PSA(LME = LME, Active.Columns = LME3$Keep, mpc = mpccm, Prefer.Set = c('FPCI'))
    LME5 = DbC.PSA(LME = LME, Active.Columns = LME4$Keep, mpc = mpccm, Prefer.Set = c('CosAsp', 'CatAr', 'CatHe', 'CatSl', 'LSF', 'PlanC', 'ProfC', 'SVF', 'Slp', 'SPI', 'TRI', 'TPI', 'VTR', 'VS', 'WI'))
    LME6 = DbC.PSA(LME = LME, Active.Columns = LME5$Keep, mpc = mpccm, Prefer.Set = c('Elev'))
    Data = LME[,LME6$Keep]
    # Expand Design Matrix to include Columns for Polynomial Terms to Max Polynomial Order
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
    # Recenter and Rescale:
    PDM.x = as.matrix(PDM)
    # Keep these: tick    
    PDM.meanx = colMeans(PDM)
    PDM.x = scale(PDM.x, center = PDM.meanx, scale = FALSE)
    PDM.normx = sqrt(colSums(PDM.x^2))
    normx = sqrt(colSums(PDM.x^2))
    names(normx) = NULL
    PDM.x = scale(PDM.x, center = FALSE, scale = normx)    
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
    # Recenter and Rescale
    IO2.2L.DM.x = as.matrix(IO2.2L.DM)
    # Keep these: tick    
    IO2.meanx = colMeans(IO2.2L.DM)
    IO2.2L.DM.x = scale(IO2.2L.DM.x, center = IO2.meanx, scale = FALSE)
    normx = sqrt(colSums(IO2.2L.DM.x^2))
    IO2.normx = sqrt(colSums(IO2.2L.DM.x^2))
    names(normx) = NULL
    IO2.2L.DM.x = scale(IO2.2L.DM.x, center = FALSE, scale = normx)    
    IO2.2L.DM = IO2.2L.DM.x
    B1.RR.P4I2.DM = data.frame(Data, IO2.2L.DM, PDM)
    # Choose from correlated pairs of covariate terms in the design matrix based on simplicity of term
    # e.g. choose polynomial term over interaction term from a correlated pair
    #      or lower order polynomial over higher order polynomial#
    # Choose Polynomial term over Interaction term from a Correlated Pair of Covariates
    B1.CC.7 = data.frame(t(combn(x = colnames(B1.RR.P4I2.DM), m = 2)))
    colnames(B1.CC.7) = c('Var1','Var2')
    B1.CC.7$Cor = numeric(nrow(B1.CC.7))
    for(i in 1:nrow(B1.CC.7)){
        B1.CC.7[i, 'Cor'] = cor(B1.RR.P4I2.DM[,paste(B1.CC.7[i,'Var1'])], B1.RR.P4I2.DM[,paste(B1.CC.7[i,'Var2'])])}
    B1.CC.7.vc = B1.CC.7[abs(B1.CC.7$Cor) >= mpccm,]
    B1.CC.7.vc = data.frame(B1.CC.7.vc, Row = 1:nrow(B1.CC.7.vc))
    Covariates.keep6 = colnames(B1.RR.P4I2.DM)
    Keep7 = rep(NA, nrow(B1.CC.7.vc))
    Drop7 = rep(NA, nrow(B1.CC.7.vc))
    for(i in 1:nrow(B1.CC.7.vc)){
        test.pair = c(paste(B1.CC.7.vc[i, 'Var1']), paste(B1.CC.7.vc[i, 'Var2']))
        if(length(grep(pattern = '_X_', x = test.pair)) == 1 &
           length(grep(pattern = '2|3|4', x = test.pair)) == 1){
            Drop7[i] = test.pair[grep(pattern = '_X_', x = test.pair)]
            Keep7[i] = test.pair[grep(pattern = '2|3|4', x = test.pair)]} else{
                if(length(grep(pattern = '_X_', x = test.pair)) == 1 &
                   length(grep(pattern = '2|3|4', x = test.pair)) == 0){
                    Drop7[i] = test.pair[grep(pattern = '_X_', x = test.pair)]
                    Keep7[i] = test.pair[-grep(pattern = '_X_', x = test.pair)]}}}
    CC.Int.Poly.Filter = data.frame(Keep = factor(x = Keep7, levels = c(NA,colnames(B1.RR.P4I2.DM))), Drop = factor(x = Drop7, levels = c(NA,colnames(B1.RR.P4I2.DM))))
    CC.Int.Poly.Filter = na.omit(CC.Int.Poly.Filter)
    for(i in 2:nrow(CC.Int.Poly.Filter)){
        if(CC.Int.Poly.Filter[i, 'Keep'] %in% CC.Int.Poly.Filter[1:(i-1), 'Drop']){
            Switch = c(paste(CC.Int.Poly.Filter[i, 'Keep']), paste(CC.Int.Poly.Filter[i, 'Drop']))
            CC.Int.Poly.Filter[i, 'Keep'] = Switch[2]
            CC.Int.Poly.Filter[i, 'Drop'] = Switch[1]}}
    Covariates.drop7 = unique(CC.Int.Poly.Filter$Drop)
    Covariates.keep6 = colnames(B1.RR.P4I2.DM)
    Covariates.keep7 = Covariates.keep6[!(Covariates.keep6 %in% Covariates.drop7)]
    # Choose Lower Order Polynomial term over Higher Order Polynomial term from a Correlated Pair of Covariates
    B1.CC.8 = data.frame(t(combn(x = Covariates.keep7, m = 2)))
    colnames(B1.CC.8) = c('Var1','Var2')
    B1.CC.8$Cor = numeric(nrow(B1.CC.8))
    for(i in 1:nrow(B1.CC.8)){
        B1.CC.8[i, 'Cor'] = cor(B1.RR.P4I2.DM[,paste(B1.CC.8[i,'Var1'])], B1.RR.P4I2.DM[,paste(B1.CC.8[i,'Var2'])])}
    B1.CC.8.vc = B1.CC.8[abs(B1.CC.8$Cor) >= mpccm,]
    B1.CC.8.vc = data.frame(B1.CC.8.vc, Row = 1:nrow(B1.CC.8.vc))
    Keep8 = rep(NA, nrow(B1.CC.8.vc))
    Drop8 = rep(NA, nrow(B1.CC.8.vc))
    for(i in 1:nrow(B1.CC.8.vc)){
        test.pair = c(paste(B1.CC.8.vc[i, 'Var1']), paste(B1.CC.8.vc[i, 'Var2']))
        if(length(grep(pattern = '4', x = test.pair)) == 1 &
           length(grep(pattern = '3', x = test.pair)) == 1){
            Drop8[i] = test.pair[grep(pattern = '4', x = test.pair)]
            Keep8[i] = test.pair[grep(pattern = '3', x = test.pair)]}
        if(length(grep(pattern = '4', x = test.pair)) == 1 &
           length(grep(pattern = '2', x = test.pair)) == 1){
            Drop8[i] = test.pair[grep(pattern = '4', x = test.pair)]
            Keep8[i] = test.pair[grep(pattern = '2', x = test.pair)]}
        if(length(grep(pattern = '4', x = test.pair)) == 1 &
           length(grep(pattern = '3|2|_X_', x = test.pair)) == 0){
            Drop8[i] = test.pair[grep(pattern = '4', x = test.pair)]
            Keep8[i] = test.pair[-grep(pattern = '4', x = test.pair)]}
        if(length(grep(pattern = '3', x = test.pair)) == 1 &
           length(grep(pattern = '2', x = test.pair)) == 1){
            Drop8[i] = test.pair[grep(pattern = '3', x = test.pair)]
            Keep8[i] = test.pair[grep(pattern = '2', x = test.pair)]}
        if(length(grep(pattern = '3', x = test.pair)) == 1 &
           length(grep(pattern = '4|2|_X_', x = test.pair)) == 0){
            Drop8[i] = test.pair[grep(pattern = '3', x = test.pair)]
            Keep8[i] = test.pair[-grep(pattern = '3', x = test.pair)]}
        if(length(grep(pattern = '2', x = test.pair)) == 1 &
           length(grep(pattern = '4|3|_X_', x = test.pair)) == 0){
            Drop8[i] = test.pair[grep(pattern = '2', x = test.pair)]
            Keep8[i] = test.pair[-grep(pattern = '2', x = test.pair)]}
        if(length(grep(pattern = '2|3|4', x = test.pair)) == 1 &
           length(grep(pattern = '_X_', x = test.pair)) == 1){
            Drop8[i] = test.pair[grep(pattern = '2|3|4', x = test.pair)]
            Keep8[i] = test.pair[grep(pattern = '_X_', x = test.pair)]}
        if(length(grep(pattern = '_X_', x = test.pair)) == 2){
            Drop8[i] = test.pair[2]
            Keep8[i] = test.pair[1]}}
    # Finally check all |correlations| < mpccm and drop members of pairs at random if there are any remaining that have |cor| > mpccm
    CC.Poly.Poly.Filter = data.frame(Keep = factor(x = Keep8, levels = c(NA,colnames(B1.RR.P4I2.DM))), Drop = factor(x = Drop8, levels = c(NA,colnames(B1.RR.P4I2.DM))))
    CC.Poly.Poly.Filter = na.omit(CC.Poly.Poly.Filter)
    for(i in 2:nrow(CC.Poly.Poly.Filter)){
        if(CC.Poly.Poly.Filter[i, 'Keep'] %in% CC.Poly.Poly.Filter[1:(i-1), 'Drop']){
            Switch = c(paste(CC.Poly.Poly.Filter[i, 'Keep']), paste(CC.Poly.Poly.Filter[i, 'Drop']))
            CC.Poly.Poly.Filter[i, 'Keep'] = Switch[2]
            CC.Poly.Poly.Filter[i, 'Drop'] = Switch[1]}}
    Covariates.drop8 = unique(CC.Poly.Poly.Filter$Drop)
    Covariates.keep8 = Covariates.keep7[!(Covariates.keep7 %in% Covariates.drop8)]
    FDM = B1.RR.P4I2.DM[, Covariates.keep8] # Filtered Design Matrix, potential covariates
    FDM.C = colnames(FDM)
    CC.FF = data.frame(t(combn(x = colnames(FDM), m = 2)))
    colnames(CC.FF) = c('Var1','Var2')
    CC.FF$Cor = numeric(nrow(CC.FF))
    for(i in 1:nrow(CC.FF)){
        CC.FF[i, 'Cor'] = cor(FDM[,paste(CC.FF[i,'Var1'])], FDM[,paste(CC.FF[i,'Var2'])])}
    CC.FF.vc = CC.FF[abs(CC.FF$Cor) > mpccm, ]
    CC.FF.Filter = CC.FF.vc[,c(1,2)]
    colnames(CC.FF.Filter) = c('Keep', 'Drop')
    levels(CC.FF.Filter$Keep) = c(levels(CC.FF.Filter$Keep), levels(CC.FF.Filter$Drop))
    levels(CC.FF.Filter$Drop) = c(levels(CC.FF.Filter$Drop), levels(CC.FF.Filter$Keep))
    if(nrow(CC.FF.Filter)>0){
    for(i in 2:nrow(CC.FF.Filter)){
        if(CC.FF.Filter[i, 'Keep'] %in% CC.FF.Filter[1:(i-1), 'Drop']){
            Switch = c(paste(CC.FF.Filter[i, 'Keep']), paste(CC.FF.Filter[i, 'Drop']))
            CC.FF.Filter[i, 'Keep'] = Switch[2]
            CC.FF.Filter[i, 'Drop'] = Switch[1]}}    
    Covariates.dropF = unique(CC.FF.Filter$Drop)
    Covariates.keepF = FDM.C[!(FDM.C %in% Covariates.dropF)]} else{
    Covariates.keepF = FDM.C
    Covariates.dropF = NULL}
    output = list(DM = FDM[,Covariates.keepF],
                  LME.meanx = LME.meanx,
                  LME.normx = LME.normx,

                  PDM.meanx = PDM.meanx,
                  PDM.normx = PDM.normx,

                  IO2.meanx = IO2.meanx,
                  IO2.normx = IO2.normx)    
    return(output)}
 
DMB.UF <- function(LME = B1.DM.LME[,3:ncol(B1.DM.LME)], RCRS = TRUE){    
    # DMB.UF = Design Matrix Builder with UnFiltered
    # LME = Linear Main Effects design matrix (as dataframe)
    # (this function will only create interactions of linear main effects)
    # mpccm = Maximum Permissible Correlation between design matrix columns
    # RCRS = ReCenter each column on 0 and ReScale each column to have magnitude 1 as for lars() function
    mpo = 4 # note the facility to alter this value would require editing of polynomial term filtering code below (currently a value of 4 is hard coded in) e.g. |2|3|4
    if(RCRS == TRUE){
        # Center all Covariates on Zero & Scale all Covariates to have Magnitude 1
        x = as.matrix(LME)
        # Keep These: tick
        LME.meanx = colMeans(x)
        x = scale(x, center = LME.meanx, scale = FALSE)
        LME.normx = sqrt(colSums(x^2))
        normx = sqrt(colSums(x^2))
        names(normx) = NULL
        x = scale(x, center = FALSE, scale = normx)
        if(unique(round(colMeans(x),10)) == 0 & unique(round(colSums(x^2),10)) == 1){
            LME = x} else{print('RCRS Error')}}
    Data = LME 
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
    # Recenter and Rescale:
    # Keep These: tick
    PDM.x = as.matrix(PDM)
    PDM.meanx = colMeans(PDM)
    PDM.x = scale(PDM.x, center = PDM.meanx, scale = FALSE)
    normx = sqrt(colSums(PDM.x^2))
    PDM.normx = sqrt(colSums(PDM.x^2))
    names(normx) = NULL
    PDM.x = scale(PDM.x, center = FALSE, scale = normx)
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
    # Keep these: tick
    # Recenter and Rescale
    IO2.2L.DM.x = as.matrix(IO2.2L.DM)
    IO2.meanx = colMeans(IO2.2L.DM)
    IO2.2L.DM.x = scale(IO2.2L.DM.x, center = IO2.meanx, scale = FALSE)
    normx = sqrt(colSums(IO2.2L.DM.x^2))
    IO2.normx = sqrt(colSums(IO2.2L.DM.x^2))
    names(normx) = NULL
    IO2.2L.DM.x = scale(IO2.2L.DM.x, center = FALSE, scale = normx)
    IO2.2L.DM = IO2.2L.DM.x
    B1.RR.P4I2.DM = data.frame(Data, IO2.2L.DM, PDM)    
    output = list(DM = B1.RR.P4I2.DM,
                  LME.meanx = LME.meanx,
                  LME.normx = LME.normx,
                  PDM.meanx = PDM.meanx,
                  PDM.normx = PDM.normx,
                  IO2.meanx = IO2.meanx,
                  IO2.normx = IO2.normx)
    return(output)}

dm.build <- function(LMEDM = B1.DM.LME[,3:ncol(B1.DM.LME)], RCRS.TF = TRUE, MPCCM = 0.95){
    if(MPCCM == 1){
        DMB.UF(LME = LMEDM, RCRS = RCRS.TF)} else{
            DMB.IF(LME = LMEDM, mpccm = MPCCM, RCRS = RCRS.TF)}}
