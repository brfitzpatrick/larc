# R Code to conduct an example analysis with functions provided in the
# Least Angle Regression Companion (LARC) project:
# https://github.com/brfitzpatrick/larc.git 

# Copyright (C) 2016 Benjamin R. Fitzpatrick, David W. Lamb and Kerrie Mengersen.

# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along with this program, in the form a text file titled 'LICENSE.
# If not, see http://www.gnu.org/licenses/.

# The program author, Benjamin R. Fitzpatrick, may be contacted via email at <ben.r.fitzpatrick@gmail.com>

# This analysis forms the case study discussed in the article:

# `Ultrahigh Dimensional Variable Selection for Interpolation of Point Referenced Spatial Data: A Digital Soil Mapping Case Study' currently in press at the Journal PLOS ONE 10.1371/journal.pone.0162489

###############################
###############################
##                           ##
##                           ##
##   Beginning of Analysis   ##
##                           ##
##                           ##
###############################
###############################

# load the required packages:
library('ggplot2')
library('grid')
library('xtable')

# Set the Working Directory to where you extracted the files

# setwd()

# Load the Data

load('Newholme_B1_Data')

# workspace should include

ls()

B1.CRS.df # the design matrix containing the means of each covariate over each 25m^2 pixel in the full cover raster, the first two columns are spatial coordinates of pixel centroids

B1.DM.LME # the design matrix consisting of the means of each covariate over each 25m^2 area around each of the locations at which the response was observed (soil core locations) with the associated spatial coordinates in the first two columns

CVTI_TS35_VS25_n500 # 500 unique samples of 35 numbers from the sequence 1:60 - these are used as indices to create the 500 unique 35 observation training sets by subsetting the response vector and design matrix 

CVTI_TS45_VS15_n500 # 500 unique samples of 45 numbers from the sequence 1:60 - these are used as indices to create the 500 unique 45 observation training sets by subsetting the response vector and design matrix 

CVTI_TS55_VS5_n500 # 500 unique samples of 55 numbers from the sequence 1:60 - these are used as indices to create the 500 unique 55 observation training sets by subsetting the response vector and design matrix 

Response.df # reponse observations accompanied by the spatial coordiantes in these observations in first two columns

######################################
######################################
##                                  ##
##                                  ##
##         Load the Functions       ##
##                                  ##
##                                  ##
######################################
######################################

source('dm.build.R')
source('sp.dm.build.R')
source('tsi.gen.R')
source('lar.cv.R')
source('lm.vs.cv.R')
source('cd.plot.R')
source('ma.pred.R')
source('rast.dm.build.R')
source('rast.sp.dm.build.R')
source('pred.rast.plot.R')

######################################
######################################
##                                  ##
##                                  ##
##     Build the Design Matrices    ##
##                                  ##
##                                  ##
######################################
######################################

# Build the expanded design matrix from the design matrix containing observations of each covariate as linear main effects at response observation locations
# The resulting design matrix will have each covariate as a linear main effect, as polynomial terms up to order four along with all possible pairwise interactions of linear terms for covariates (i.e. all possible products of pairs of linear terms for the covariates, some 63 choose 2 terms interaction terms)

# a) Full Design Matrix
  
  # i) for observations at the 60 soil core sample locations

system.time(n60.P4I2.DM <- dm.build(LMEDM = B1.DM.LME[,3:ncol(B1.DM.LME)],
                                    RCRS.TF = TRUE,
                                    MPCCM = 1)) # 0.385 sec
  # arguments:
  # LMEDM = Linear Main Effects Design Matrix (without the spatial coordinates)
  # RCRS.TF = Re-Center and Re-Scale Design Matrix TRUE/FALSE
  # MPCCM = Maximum Permitted Correlation Coefficient Magnitude between covariate pairs

  # ii) for the raster full design matrix (i.e. the observations accompanying the pixels of the full cover raster for the study site) expanded,  recentered and rescaled with identical transformations as were applied to the design matrix composed of the observations at the 60 soil core sample locations

system.time(Rast.P4I2.DM <- rast.dm.build(LME.rast = B1.CRS.df[,3:ncol(B1.CRS.df)],
                                          n60.LME.meanx = n60.P4I2.DM$LME.meanx,
                                          n60.LME.normx = n60.P4I2.DM$LME.normx,
                                          n60.PDM.meanx = n60.P4I2.DM$PDM.meanx,
                                          n60.PDM.normx = n60.P4I2.DM$PDM.normx,
                                          n60.IO2.meanx = n60.P4I2.DM$IO2.meanx,
                                          n60.IO2.normx = n60.P4I2.DM$IO2.normx)) # 0.976 

# Re-attach the pixel centroid coordinates:
Rast.P4I2.DM.Coords <- data.frame(B1.CRS.df[,1:2],Rast.P4I2.DM$DM)



# b) Design Matrix Filtered to enforce a maximum permitted correlation coefficient magnitude between remaining covariate terms of 0.95

  # i) for observations at the 60 soil core sample locations

system.time(n60.P4I2.DM.F95 <- dm.build(LMEDM = B1.DM.LME[,3:ncol(B1.DM.LME)],
                                        RCRS.TF = TRUE,
                                        MPCCM = 0.95)) # 24 mins 33 sec

# c) Design Matrix Filtered to enforce a maximum permitted correlation coefficient magnitude between remaining covariate terms of 0.8

  # i) for observations at the 60 soil core sample locations

system.time(n60.P4I2.DM.F80 <- dm.build(LMEDM = B1.DM.LME[,3:ncol(B1.DM.LME)],
                                        RCRS.TF = TRUE,
                                        MPCCM = 0.8)) # 2 min 31 sec

# d) Design Matrix Filtered to enforce a maximum permitted correlation coefficient magnitude between remaining covariate terms of 0.6

  # i) for observations at the 60 soil core sample locations

system.time(n60.P4I2.DM.F60 <- dm.build(LMEDM = B1.DM.LME[,3:ncol(B1.DM.LME)],
                                        RCRS.TF = TRUE,
                                        MPCCM = 0.6)) # 17.297 sec

# e) Design Matrix Filtered to enforce a maximum permitted correlation coefficient magnitude between remaining covariate terms of 0.4

  # i) for observations at the 60 soil core sample locations

system.time(n60.P4I2.DM.F40  <- dm.build(LMEDM = B1.DM.LME[,3:ncol(B1.DM.LME)],
                                         RCRS.TF = TRUE,
                                         MPCCM = 0.4)) # 2.676 sec


######################################
######################################
##                                  ##
##                                  ##
##     Generate (or Load) the       ##
##  Training Set Indices used in    ##
##   the Cross Validation Scheme    ##
##                                  ##
##                                  ##
######################################
######################################

# Generate the Indices for constructing 500 unique Training sets of 35 observations each
  ## note these indices are generated stochastically so if you want to reproduce our results exactly please use the collections of indices: 'CVTI_TS35_VS25_n500', 'CVTI_TS45_VS15_n500' and 'CVTI_TS55_VS5_n500' provided in the workspace

# system.time(CVTI_TS35_VS25_n500_V2 <- tsi.gen(n = 500, all = 60, train = 35)) # 6 min 15 sec 
  # arguments
  # n = number of unique training set indices to generate
  # all = number of observations to divide into training sets and validations sets
  # train = number of observations each training set should contain

########################################
########################################
##                                    ##
##                                    ##
##  Compare the different variable    ##
##  selection algorithms              ##  
##  (each used within identical       ##
##  cross validation and model        ##
##  averaging schemes)                ##
##                                    ##
##                                    ##
########################################

 # a) Least Angle Regression Variable Selection on 500 unique training sets constructed from a design matrix filtered to enforce a maximum permitted correlation coefficient magnitude between remaining covariate pairs of 0.95

  # i) Run the variable selection within the cross validation scheme

system.time(LAR.F95.TS35.n500 <- lar.cv(Response = Response.df$Carbon,
                                        DM = n60.P4I2.DM.F95$DM ,
                                        CVTI = CVTI_TS35_VS25_n500))
  # arguments:
  # Response = the vector of response observations
  # DM = the design matrix
  # CVTI = the set of indices with which to subset (by rows) the full design matrix into training and validation sets
  # time elapsed: predicted 21 mins 7 sec took 21 mins 18 sec

  # ii) calculate the model averaged predictions at the response sample locations

MAP.LAR.F95.TS35.n500 <- ma.pred(cotpa.mat = LAR.F95.TS35.n500$Design.mat,
                                vspe.mat = LAR.F95.TS35.n500$VSPE.mat,
                                sc.ls = LAR.F95.TS35.n500$SC.ls,
                                Response = LAR.F95.TS35.n500$Response,
                                CV.TI = LAR.F95.TS35.n500$CV.TI,
                                Response.TI.means = LAR.F95.TS35.n500$Response.TI.means,
                                DM.TI.C.means.store = LAR.F95.TS35.n500$DM.TI.C.means.store,
                                DM.TI.C.norms.store = LAR.F95.TS35.n500$DM.TI.C.norms.store,
                                plot.wgf = TRUE,
                                calc.R.2 = TRUE)

MAP.LAR.F95.TS35.n500$R2 # 0.596329

summary(LAR.F95.TS35.n500$VSPE.mat) # ingredient in boxplots comparing VSEPE distributions among different variable selection techniques

 # b) Least Angle Regression Variable Selection on 500 unique training sets constructed from a design matrix filtered to enforce a maximum permitted correlation coefficient magnitude between remaining covariate pairs of 0.4

  # i) Run the variable selection within the cross validation scheme

system.time(LAR.F40.TS35.n500 <- lar.cv(Response = Response.df$Carbon,
                                        DM = n60.P4I2.DM.F40$DM, 
                                        CVTI = CVTI_TS35_VS25_n500)) # 52.5 sec

  # ii) calculate the model averaged predictions at the response sample locations

MAP.LAR.F40.TS35.n500 <- ma.pred(cotpa.mat = LAR.F40.TS35.n500$Design.mat,
                                vspe.mat = LAR.F40.TS35.n500$VSPE.mat,
                                sc.ls = LAR.F40.TS35.n500$SC.ls,
                                Response = LAR.F40.TS35.n500$Response,
                                CV.TI = LAR.F40.TS35.n500$CV.TI,
                                Response.TI.means = LAR.F40.TS35.n500$Response.TI.means,
                                DM.TI.C.means.store = LAR.F40.TS35.n500$DM.TI.C.means.store,
                                DM.TI.C.norms.store = LAR.F40.TS35.n500$DM.TI.C.norms.store,
                                plot.wgf = TRUE,
                                calc.R.2 = TRUE)

MAP.LAR.F40.TS35.n500$R2 # 0.3666025

 # c) Branch and Bound Exhaustive Search Variable Selection on 500 unique training sets constructed from a design matrix filtered to enforce a maximum permitted correlation coefficient magnitude between remaining covariate pairs of 0.4

  # i) Run the variable selection within the cross validation scheme

system.time(Exh.F40.TS35.n500 <- lm.vs.cv(Response = Response.df$Carbon,
                                          DM = n60.P4I2.DM.F40$DM,
                                          CVTI = CVTI_TS35_VS25_n500,
                                          vs_type = c('exhaustive'))) # 2 mins 2 sec

  # ii) calculate the model averaged predictions at the response sample locations
        
MAP.Exh.F40.TS35.n500 <- ma.pred(cotpa.mat = Exh.F40.TS35.n500$Design.mat,
                                 vspe.mat = Exh.F40.TS35.n500$VSPE.mat,
                                 sc.ls = Exh.F40.TS35.n500$SC.ls,
                                 Response = Exh.F40.TS35.n500$Response,
                                 CV.TI = Exh.F40.TS35.n500$CV.TI,
                                 Response.TI.means = Exh.F40.TS35.n500$Response.TI.means,
                                 DM.TI.C.means.store = Exh.F40.TS35.n500$DM.TI.C.means.store,
                                 DM.TI.C.norms.store = Exh.F40.TS35.n500$DM.TI.C.norms.store,
                                 plot.wgf = TRUE,
                                 calc.R.2 = TRUE)

MAP.Exh.F40.TS35.n500$R2 #  0.2882439

 # d) Stepwise Sequential Replacement Variable Selection on 500 unique training sets constructed from a design matrix filtered to enforce a maximum permitted correlation coefficient magnitude between remaining covariate pairs of 0.4

  # i) Run the variable selection within the cross validation scheme

system.time(Seq.F40.TS35.n500 <- lm.vs.cv(Response = Response.df$Carbon,
                                          DM = n60.P4I2.DM.F40$DM,
                                          CVTI = CVTI_TS35_VS25_n500,
                                          vs_type = c('seqrep'))) # 43 sec

  # ii) calculate the model averaged predictions at the response sample locations

MAP.Seq.F40.TS35.n500 <- ma.pred(cotpa.mat = Seq.F40.TS35.n500$Design.mat,
                                vspe.mat = Seq.F40.TS35.n500$VSPE.mat,
                                sc.ls = Seq.F40.TS35.n500$SC.ls,
                                Response = Seq.F40.TS35.n500$Response,
                                CV.TI = Seq.F40.TS35.n500$CV.TI,
                                Response.TI.means = Seq.F40.TS35.n500$Response.TI.means,
                                DM.TI.C.means.store = Seq.F40.TS35.n500$DM.TI.C.means.store,
                                DM.TI.C.norms.store = Seq.F40.TS35.n500$DM.TI.C.norms.store,
                                plot.wgf = TRUE,
                                calc.R.2 = TRUE)

MAP.Seq.F40.TS35.n500$R2 # 0.3054765

 # e) Stepwise Forwards Variable Selection on 500 unique training sets constructed from a design matrix filtered to enforce a maximum permitted correlation coefficient magnitude between remaining covariate pairs of 0.4

  # i) Run the variable selection within the cross validation scheme

system.time(Fwd.F40.TS35.n500 <- lm.vs.cv(Response = Response.df$Carbon,
                                          DM = n60.P4I2.DM.F40$DM,
                                          CVTI = CVTI_TS35_VS25_n500,
                                          vs_type = c('forward'))) # 39 sec

  # ii) calculate the model averaged predictions at the response sample locations

MAP.Fwd.F40.TS35.n500 <- ma.pred(cotpa.mat = Fwd.F40.TS35.n500$Design.mat,
                                vspe.mat = Fwd.F40.TS35.n500$VSPE.mat,
                                sc.ls = Fwd.F40.TS35.n500$SC.ls,
                                Response = Fwd.F40.TS35.n500$Response,
                                CV.TI = Fwd.F40.TS35.n500$CV.TI,
                                Response.TI.means = Fwd.F40.TS35.n500$Response.TI.means,
                                DM.TI.C.means.store = Fwd.F40.TS35.n500$DM.TI.C.means.store,
                                DM.TI.C.norms.store = Fwd.F40.TS35.n500$DM.TI.C.norms.store,
                                plot.wgf = TRUE,
                                calc.R.2 = TRUE)

MAP.Fwd.F40.TS35.n500$R2 # 0.3045973

 # f) Stepwise Backwards Variable Selection on 500 unique training sets constructed from a design matrix filtered to enforce a maximum permitted correlation coefficient magnitude between remaining covariate pairs of 0.4

  # i) Run the variable selection within the cross validation scheme

system.time(Bwd.F40.TS35.n500 <- lm.vs.cv(Response = Response.df$Carbon,
                                          DM = n60.P4I2.DM.F40$DM,
                                          CVTI = CVTI_TS35_VS25_n500,
                                          vs_type = c('backward'))) # 40 sec

  # ii) calculate the model averaged predictions at the response sample locations

MAP.Bwd.F40.TS35.n500 <- ma.pred(cotpa.mat = Bwd.F40.TS35.n500$Design.mat,
                                vspe.mat = Bwd.F40.TS35.n500$VSPE.mat,
                                sc.ls = Bwd.F40.TS35.n500$SC.ls,
                                Response = Bwd.F40.TS35.n500$Response,
                                CV.TI = Bwd.F40.TS35.n500$CV.TI,
                                Response.TI.means = Bwd.F40.TS35.n500$Response.TI.means,
                                DM.TI.C.means.store = Bwd.F40.TS35.n500$DM.TI.C.means.store,
                                DM.TI.C.norms.store = Bwd.F40.TS35.n500$DM.TI.C.norms.store,
                                plot.wgf = TRUE,
                                calc.R.2 = TRUE)

MAP.Bwd.F40.TS35.n500$R2  # 0.2382116

###################
###################
##               ##
##               ##
##    Table 2    ##
##               ##
##               ##
###################
###################

Table.2 <- rbind(summary(abs(as.vector(LAR.F95.TS35.n500$VSPE.mat)))[1:6],
                 summary(abs(as.vector(LAR.F40.TS35.n500$VSPE.mat)))[1:6],
                 summary(abs(as.vector(Exh.F40.TS35.n500$VSPE.mat)))[1:6],
                 summary(abs(as.vector(Seq.F40.TS35.n500$VSPE.mat)))[1:6],
                 summary(abs(as.vector(Fwd.F40.TS35.n500$VSPE.mat)))[1:6],
                 summary(abs(as.vector(Bwd.F40.TS35.n500$VSPE.mat)))[1:6])

Table.2 <- data.frame(Method = c('LAR','LAR','Exh','Seq','Fwd','Bwd'),
                      r.max = c(0.95, 0.4, 0.4, 0.4, 0.4, 0.4),
                      Table.2)
                      
Table.2$R2 <- round(x = c(MAP.LAR.F95.TS35.n500$R2,
                          MAP.LAR.F40.TS35.n500$R2,
                          MAP.Exh.F40.TS35.n500$R2,
                          MAP.Seq.F40.TS35.n500$R2,
                          MAP.Fwd.F40.TS35.n500$R2,
                          MAP.Bwd.F40.TS35.n500$R2),
                    digits = 4)

Table.2

Table.2.xtab <- xtable(x = Table.2,
       digits = 4,
       display = c('s','s','g','g','g','g','g','g','g','g'),
       caption = 'Distribution Summary Statistics of Validation Set Element Prediction Errors from the Cross Validation Scheme and Coefficient of Determinations of the Model Averaged Predictions of the Entire Response Vector ',
       label = 'Tab:7')

print(x = Table.2.xtab, include.rownames = FALSE)

####################
####################
##                ##
##                ##
##    Table 3     ##
##                ##
##                ##
####################
####################

summary(factor(names(unlist(LAR.F95.TS35.n500$SC.ls))))

# If you've used the provided version of 'CVTI_TS35_VS25_n500' you should have
#          ECA.Nov.4             LSF.3           DVI.May                WI 
#                219               139               102               100 
#      ECA.Feb_X_Slp     Mag.II_X_FPCI      SVF_X_Mag.IV             Slp.2 
#                 95                95                94                89 
#   ECA.Feb_X_SR.May         LSF_X_SVF ECA.Nov_X_DVI.Nov        Elev_X_SVF 
#                 88                82                78                76 
#  ECA.Feb_X_DVI.Nov         ECA.Nov.3    ECA.Feb_X_Elev            Elev.3 
#                 74                73                72                72 
#               Elev      SR.May_X_Slp     DVI.May_X_Slp      Elev_X_PlanC 
#                 71                70                65                63 


Table.3 <- data.frame(Covariate = names(summary(factor(names(unlist(LAR.F95.TS35.n500$SC.ls)))))[1:15],
                      Freq = as.vector(summary(factor(names(unlist(LAR.F95.TS35.n500$SC.ls)))))[1:15])

Table.3 <- data.frame(Table.3,
                      CC1 = factor(x = character(length = nrow(Table.3)), levels = colnames(n60.P4I2.DM$DM)),
                      CC2 = factor(x = character(length = nrow(Table.3)), levels = colnames(n60.P4I2.DM$DM)),
                      CC3 = factor(x = character(length = nrow(Table.3)), levels = colnames(n60.P4I2.DM$DM)),
                      CC4 = factor(x = character(length = nrow(Table.3)), levels = colnames(n60.P4I2.DM$DM)),
                      CC5 = factor(x = character(length = nrow(Table.3)), levels = colnames(n60.P4I2.DM$DM)),
                      CC6 = factor(x = character(length = nrow(Table.3)), levels = colnames(n60.P4I2.DM$DM)),
                      CC7 = factor(x = character(length = nrow(Table.3)), levels = colnames(n60.P4I2.DM$DM)),
                      CC8 = factor(x = character(length = nrow(Table.3)), levels = colnames(n60.P4I2.DM$DM)),
                      CC9 = factor(x = character(length = nrow(Table.3)), levels = colnames(n60.P4I2.DM$DM)),
                      CC10 = factor(x = character(length = nrow(Table.3)), levels = colnames(n60.P4I2.DM$DM)))

for(i in 1:nrow(Table.3)){
    covar.i = Table.3[i,'Covariate']
    cor.test.mat = cor(n60.P4I2.DM$DM[,paste(covar.i)], n60.P4I2.DM$DM)
    cor.test.df = data.frame(t(cor.test.mat))
    colnames(cor.test.df) = c('cor')
    cor.test.df$covar = rownames(cor.test.df)
    cor.cov.names = cor.test.df[abs(cor.test.df$cor) > 0.95 & round(x = cor.test.df$cor, digits = 6) < 1.000000, 'covar']
    if(length(cor.cov.names) > 0){
        if(length(cor.cov.names) <= (ncol(Table.3) - 2)){
            Table.3[i,3:(2+length(cor.cov.names)) ] <-  cor.cov.names} else{
                n.extra.cols =  length(cor.cov.names) - (ncol(Table.3) - 2)
                for(j in 1:n.extra.cols){
                    Table.3 <- data.frame(Table.3, factor(x = character(length = nrow(Table.3)), levels = colnames(n60.P4I2.DM$DM)))}
                colnames(Table.3) <- c('Covariate','Frequ', paste('CC', 1:n.extra.cols, sep = ''))
                Table.3[i,3:(2+length(cor.cov.names)) ] <-  cor.cov.names}}}

Drop.These <- numeric()

for(i in 1:ncol(Table.3)){
    if( length(na.omit(Table.3[,i])) == 0 ){ Drop.These <- c(Drop.These,i)}}

Table.3

Table.3 <- Table.3[,-Drop.These]

Table.3.xtab <- xtable(x = Table.3,
                       digits = 4,
                       display = c('s','g','g','g','g','g','g','g'),
                       label = 'Tab:8',
                       caption = '15 most frequently selected covariate terms across the cross validation scheme, frequency of selection and the covariates were filtered out in order to retain the covariate in question')
                  

print(x = Table.3.xtab, include.rownames = FALSE)


###################
###################
##               ##
##               ##
##   Figure 1    ##
##               ##
##               ##
###################
###################
#
# Histograms depicting the distributions of subset sizes selected by each of
# the variable selection techniques trialed here on the 27 covariate term
# design matrix constructed by filtering the full design matrix to enforce
# a maximum permitted correlation coefficient magnitude between remaining
# covariate pairs of 0.4

Figure.1.df.nrow = sum(c(length(LAR.F40.TS35.n500$SSS.num),
                             length(Exh.F40.TS35.n500$SSS.num),
                             length(Seq.F40.TS35.n500$SSS.num),
                             length(Fwd.F40.TS35.n500$SSS.num),
                             length(Bwd.F40.TS35.n500$SSS.num)))

Tech.levels <- c('LAR', 'Exh', 'Seq', 'Fwd', 'Bwd')

Figure.1.df <- data.frame(Tech = factor(x = character(length = Figure.1.df.nrow), levels = Tech.levels), SSS = numeric(length =  Figure.1.df.nrow))   

nr.1 = length(LAR.F40.TS35.n500$SSS.num)

Figure.1.df[1:nr.1,'Tech'] <- 'LAR'

nr.2 = nr.1 + length(Exh.F40.TS35.n500$SSS.num)

Figure.1.df[(nr.1+1):(nr.2),'Tech'] <- 'Exh'

nr.3 = nr.2 + length(Seq.F40.TS35.n500$SSS.num)

Figure.1.df[(nr.2+1):(nr.3),'Tech'] <- 'Seq'

nr.4 = nr.3 + length(Fwd.F40.TS35.n500$SSS.num)

Figure.1.df[(nr.3+1):(nr.4),'Tech'] <- 'Fwd'

nr.5 = nr.4 + length(Bwd.F40.TS35.n500$SSS.num)

Figure.1.df[(nr.4+1):(nr.5),'Tech'] <- 'Bwd'

Figure.1.df[Figure.1.df$Tech == 'LAR', 'SSS'] <- LAR.F40.TS35.n500$SSS.num

Figure.1.df[Figure.1.df$Tech == 'Exh', 'SSS'] <- Exh.F40.TS35.n500$SSS.num

Figure.1.df[Figure.1.df$Tech == 'Seq', 'SSS'] <- Seq.F40.TS35.n500$SSS.num

Figure.1.df[Figure.1.df$Tech == 'Fwd', 'SSS'] <- Fwd.F40.TS35.n500$SSS.num

Figure.1.df[Figure.1.df$Tech == 'Bwd', 'SSS'] <- Bwd.F40.TS35.n500$SSS.num

Figure.1.p <- ggplot(aes(x = SSS), data = Figure.1.df)

Figure.1.p <- Figure.1.p +
    geom_histogram(fill = 'black', colour = 'white', origin = -0.5, binwidth = 1) +
    facet_wrap(facets = ~Tech, nrow = 2) +
    labs(x = 'Selected Subset Size', y = 'Frequency of Selection') + # , title = 'Cross Validation Runs on 500 Unique Training Sets of 35 Each Filtered to have a Maximum Covariate Pair Correlation Magnitude of 0.4'
    theme(text = element_text(size = 18, colour = 'black'),
          axis.text.x = element_text(size = 18, colour = 'black', vjust = 0.5, angle = 90),
          axis.text.y = element_text(size = 18, angle = 0, colour = 'black'),
          legend.position = 'right',
          legend.key.height = unit(x = 1.4, units = 'inches'),
          panel.background = element_rect(fill = 'white'),
          panel.grid.major = element_line(colour = 'darkgrey'),
          panel.grid.minor = element_line(colour = 'grey'),
          axis.ticks = element_line(colour = 'darkgrey'))

Figure.1.p

ggsave(plot = Figure.1.p, filename = 'Figure_1.pdf', width = 12, height = 10, units = 'in')

ggsave(plot = Figure.1.p, filename = 'Figure_1.tiff', width = 12, height = 10, units = 'in', dpi = 600)

##############################
##############################
##                          ##
##                          ##
##         Figure 2         ##
##                          ##
##                          ##
##############################
##############################

# Chord Diagram Depicting covariate term selection frequencies
# across the 500 selected models obtained from applying the
# Least Angle Regression variable selection algorithm to the
# 35 observation training sets constructed from design matrices
# produced by filtering the full design matrix to enforce a maximum
# permitted correlation coefficient magnitude between covariate pairs
# of 0.95.
# Poincare segments represent interaction terms between the covariates they connect.
# Covariate Acronyms are expanded in Table 1.

Figure.2 <- cd.plot(P4I2_UnFiltered_DM = n60.P4I2.DM$DM,
                     LAR_Results = LAR.F95.TS35.n500,
                     graph.title = c(''),
                     Covariate.Lab.Size = 1,
                     Point.Size = 1,
                     Title.Size = 4,
                     Legend.Height.Width = unit(2,'mm'),
                     Legend.Position = 'right',
                     Line.Width = 0.1)

Figure.2

###################################
###################################
##                               ##
##                               ##
##        Areal Inference        ## 
##                               ##
##                               ##
###################################
###################################

# for these Figures with multiple plots as panels we will use the multiplot function from:    
# <http://wiki.stdout.org/rcookbook/Graphs/Multiple>

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  if (numPlots==1) {
    print(plots[[1]])
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

################################################################################
#                                                                              #
#    Figure 3                                                                  #
#                                                                              #
################################################################################
#                                                                              #
#    The observed soil organic carbon percentages (%SOC) at the soil core      #
#    locations have been represented by the shade filling the circles located  #
#    at each of the soil core sample locations.                                #
#    The observed %SOC values have been represented with the same grey scale   #
#    as the predicted %SOC values and associated uncertainties in the rasters. #
#                                                                              #
#    Figure 3 (a)                                                              #
#    The sum of the covariate based predictions and the predictions from the   #
#    model for the spatial component of the errors from the covariate based    #
#    model.                                                                    #
#    The more westerly pixel annotated with a vertical cross represents a      #
#    predicted %SOC value of 17.92 and the more easterly pixel annotated       #
#    with a vertical cross represents a predicted %SOC value of 9.54.          #
#                                                                              #
#    Figure 3 (b)                                                              #
#    The uncertainty estimated to accompany the %SOC predictions.              #
#    The three pixels annotated with vertical crosses represent estimates of   #
#    the uncertainty associated with the model-averaged predicted %SOC values  #
#    of 20.57, 21.66 and 43.66 units on the predicted %SOC scale.              #
#    The estimated uncertainty of 43.66 being the most westerly of these       #
#    three pixels and the estimated uncertainty of 20.57 being the most        #
#    northerly of these three pixels.                                          #
#                                                                              #
################################################################################

Figure.3.ls <- list()

system.time(
Areal.MAP.LAR.F95.TS35.n500 <- ma.pred(cotpa.mat = Rast.P4I2.DM$DM[,colnames(LAR.F95.TS35.n500$Design.mat)],
                                       vspe.mat = LAR.F95.TS35.n500$VSPE.mat,
                                       sc.ls = LAR.F95.TS35.n500$SC.ls,
                                       Response = LAR.F95.TS35.n500$Response,
                                       CV.TI = LAR.F95.TS35.n500$CV.TI,
                                       Response.TI.means = LAR.F95.TS35.n500$Response.TI.means,
                                       DM.TI.C.means.store = LAR.F95.TS35.n500$DM.TI.C.means.store,
                                       DM.TI.C.norms.store = LAR.F95.TS35.n500$DM.TI.C.norms.store,
                                       pred.disp.est = c('variance','interval.width'),
                                       interval = 0.95,
                                       plot.wgf = FALSE,
                                       calc.R.2 = FALSE)) # 45 sec

# Fit a model for the spatial component of the residuals from the first round of modelling Soil Carbon with the environmental covariates:

MAP.LAR.F95.TS35.n500.Resid <- Response.df$Carbon - MAP.LAR.F95.TS35.n500$Results[,'MA.Pred']

MAP.LAR.F95.TS35.n500.Resid <- Response.df$Carbon - MAP.LAR.F95.TS35.n500$Results[,'MA.Pred']

SpEr.DM <- sp.dm.build(LME = Response.df[,c('Easting', 'Northing')], mpc = 0.95, mppo = 12)

system.time(SpEr.LAR.F95.TS35.n500 <- lar.cv(Response = MAP.LAR.F95.TS35.n500.Resid,
                                             DM = SpEr.DM$DM,
                                             CVTI = CVTI_TS35_VS25_n500)) # 45 sec

SpEr.MAP.LAR.F95.TS35.n500 <- ma.pred(cotpa.mat = SpEr.LAR.F95.TS35.n500$Design.mat,
                                      vspe.mat = SpEr.LAR.F95.TS35.n500$VSPE.mat,
                                      sc.ls = SpEr.LAR.F95.TS35.n500$SC.ls,
                                      Response = SpEr.LAR.F95.TS35.n500$Response,
                                      CV.TI = SpEr.LAR.F95.TS35.n500$CV.TI,
                                      Response.TI.means = SpEr.LAR.F95.TS35.n500$Response.TI.means,
                                      DM.TI.C.means.store = SpEr.LAR.F95.TS35.n500$DM.TI.C.means.store,
                                      DM.TI.C.norms.store = SpEr.LAR.F95.TS35.n500$DM.TI.C.norms.store,
                                      plot.wgf = TRUE,
                                      calc.R.2 = TRUE)

SpEr.MAP.LAR.F95.TS35.n500$R2 # 0.143

system.time(
Rast.SpEr.DM <- rast.sp.dm.build(LME = B1.CRS.df[,c('Easting', 'Northing')],
                                 mppo = 12,
                                 n60.DM = SpEr.DM$DM,
                                 n60.LME.meanx = SpEr.DM$LME.meanx,
                                 n60.LME.normx = SpEr.DM$LME.normx,
                                 n60.Poly.C.means = SpEr.DM$Poly.C.means,
                                 n60.Poly.C.norms = SpEr.DM$Poly.C.norms,
                                 n60.PDM.meanx = SpEr.DM$PDM.meanx,
                                 n60.PDM.normx = SpEr.DM$PDM.normx,
                                 n60.Int.C.meanx = SpEr.DM$Int.C.meanx,
                                 n60.Int.C.normx = SpEr.DM$Int.C.normx)) # 0.087 sec

Rast.SpEr.MAP.LAR.F95.TS35.n500 <- ma.pred(cotpa.mat = Rast.SpEr.DM$DM,
                                           vspe.mat = SpEr.LAR.F95.TS35.n500$VSPE.mat,
                                           sc.ls = SpEr.LAR.F95.TS35.n500$SC.ls,
                                           Response = SpEr.LAR.F95.TS35.n500$Response,
                                           CV.TI = SpEr.LAR.F95.TS35.n500$CV.TI,
                                           Response.TI.means = SpEr.LAR.F95.TS35.n500$Response.TI.means,
                                           DM.TI.C.means.store = SpEr.LAR.F95.TS35.n500$DM.TI.C.means.store,
                                           DM.TI.C.norms.store = SpEr.LAR.F95.TS35.n500$DM.TI.C.norms.store,
                                           plot.wgf = TRUE,
                                           calc.R.2 = FALSE)

Rast.SpEr.MAP.LAR.F95.TS35.n500.Coords.df <- data.frame(B1.CRS.df[,1:2],
                                                        Rast.SpEr.MAP.LAR.F95.TS35.n500$Result)

SpEr.MAP <- Rast.SpEr.MAP.LAR.F95.TS35.n500.Coords.df[,c('Easting', 'Northing', 'MA.Pred')]

CVLAR.MAP <- data.frame(B1.CRS.df[,c('Easting','Northing')], MA.Pred = Areal.MAP.LAR.F95.TS35.n500$Results[,'MA.Pred'])

MAP.p.SpEr.df <- data.frame(B1.CRS.df[,c('Easting','Northing')],
                            MAP.p.SpEr = (Areal.MAP.LAR.F95.TS35.n500$Results[,'MA.Pred'] + SpEr.MAP$MA.Pred))

MAP.p.SpEr.LAR.F95.TS35.n500.grey.p <- pred.rast.plot(
    rast.df = MAP.p.SpEr.df,
    fill.name = '% SOC', 
    points.df = Response.df[,c('Easting', 'Northing','Carbon')],
    point.fill.column.name = 'Carbon',  
    colour.scale.type = 'grey', 
    auto.col.lims = TRUE,
    colour.scale.limits = c(0,7), 
    out.of.range.colour = 'white',
    overlay.points = TRUE,
    point.boundary.colour = 'white',  
    point.size = 2.5,
    plot.title = '(a)',
    title.text.size = 25,
    axis.text.size = 25,
    axis.tick.step = 100,
    x.axis.first.tick = 368100,
    x.axis.last.tick = 369800,
    y.axis.first.tick = 663200,
    y.axis.last.tick = 6634200) + annotate(geom = 'point', shape = 3, x = c(368842, 368944), y = c(6633649, 6633649), colour = 'black')

MAP.p.SpEr.LAR.F95.TS35.n500.grey.p

Figure.3.ls[[1]] <- MAP.p.SpEr.LAR.F95.TS35.n500.grey.p

##

Rast.Uncert.LAR.F95.TS35.n500.grey.p <- pred.rast.plot(
    rast.df = data.frame(B1.CRS.df[,c('Easting','Northing')], Fill = Areal.MAP.LAR.F95.TS35.n500$Results[,'Int.Width']),
    fill.name = '% SOC', 
    points.df = Response.df, 
    point.fill.column.name = 'Carbon',  
    colour.scale.type = 'grey', 
    auto.col.lims = TRUE,
    colour.scale.limits = c(0,20),
    out.of.range.colour = 'white',
    overlay.points = TRUE,
    point.boundary.colour = 'white',
    point.size = 2.5,
    plot.title = '(b)',
    title.text.size = 25,
    axis.text.size = 25,
    axis.tick.step = 100,
    x.axis.first.tick = 368100,
    x.axis.last.tick = 369800,
    y.axis.first.tick = 663200,
    y.axis.last.tick = 6634200) +
    annotate(geom = 'point', shape = 3, x = c(368842, 368944, (368944 + 25*3) ), y = c(6633649, 6633649, (6633649 + (25*5))), colour = 'black')

Rast.Uncert.LAR.F95.TS35.n500.grey.p

Figure.3.ls[[2]] <- Rast.Uncert.LAR.F95.TS35.n500.grey.p

pdf(file = 'Figure_3.pdf', height = 20, width = 12)
multiplot(plotlist = Figure.3.ls, cols = 1, layout = matrix(data = c(1,2), nrow = 2, byrow = TRUE))
dev.off()

tiff(filename = 'Figure_3.tiff', height = 20, width = 12, units = "in", pointsize = 12, compression = c('lzw'), bg = "white", res = 600)
multiplot(plotlist = Figure.3.ls, cols = 1, layout = matrix(data = c(1,2), nrow = 2, byrow = TRUE))
dev.off()


################################################################################
################################################################################
##                                                                            ## 
##   End of Analysis Presented in Main Article                                ##
##                                                                            ##
##   Further Tables and Figures are Generated below for use in the            ##
##                                                                            ##
##   Supplementary Materials                                                  ##
##                                                                            ##
################################################################################
################################################################################

# Alternative Colour (heat) Version of Figure X for Revised Manuscript

#     colour.scale.type = 'rainbow', 
#     colour.scale.type = 'heat.custom', 

F1fRM.heat.ls <- list()

MAP.p.SpEr.LAR.F95.TS35.n500.heat.p <- pred.rast.plot(
    rast.df = MAP.p.SpEr.df,
    fill.name = '% SOC', 
    points.df = Response.df[,c('Easting', 'Northing','Carbon')],
    point.fill.column.name = 'Carbon',  
    colour.scale.type = 'heat.custom', 
    auto.col.lims = TRUE,
    colour.scale.limits = c(0,7), 
    out.of.range.colour = 'white',
    overlay.points = TRUE,
    point.boundary.colour = 'white',  
    point.size = 2.5,
    plot.title = '(a)',
    title.text.size = 25,
    axis.text.size = 25,
    axis.tick.step = 100,
    x.axis.first.tick = 368100,
    x.axis.last.tick = 369800,
    y.axis.first.tick = 663200,
    y.axis.last.tick = 6634200) +
    annotate(geom = 'point', shape = 3, x = c(368842, 368944), y = c(6633649, 6633649), colour = 'black')

MAP.p.SpEr.LAR.F95.TS35.n500.heat.p

F1fRM.heat.ls[[1]] <- MAP.p.SpEr.LAR.F95.TS35.n500.heat.p

##

Rast.Uncert.LAR.F95.TS35.n500.heat.p <- pred.rast.plot(
    rast.df = data.frame(B1.CRS.df[,c('Easting','Northing')], Fill = Areal.MAP.LAR.F95.TS35.n500$Results[,'Int.Width']),
    fill.name = '% SOC', 
    points.df = Response.df, 
    point.fill.column.name = 'Carbon',  
    colour.scale.type = 'heat.custom', 
    auto.col.lims = TRUE,
    colour.scale.limits = c(0,20),
    out.of.range.colour = 'white',
    overlay.points = TRUE,
    point.boundary.colour = 'white',
    point.size = 2.5,
    plot.title = '(b)',
    title.text.size = 25,
    axis.text.size = 25,
    axis.tick.step = 100,
    x.axis.first.tick = 368100,
    x.axis.last.tick = 369800,
    y.axis.first.tick = 663200,
    y.axis.last.tick = 6634200) +
    annotate(geom = 'point', shape = 3, x = c(368842, 368944, (368944 + 25*3) ), y = c(6633649, 6633649, (6633649 + (25*5))), colour = 'black')

Rast.Uncert.LAR.F95.TS35.n500.heat.p

F1fRM.heat.ls[[2]] <- Rast.Uncert.LAR.F95.TS35.n500.heat.p

pdf(file = 'Figure_3_Heat_Colours.pdf', height = 20, width = 12)
multiplot(plotlist = F1fRM.heat.ls, cols = 1, layout = matrix(data = c(1,2), nrow = 2, byrow = TRUE))
dev.off()

tiff(filename = 'Figure_3_Heat_Colours.tiff', height = 20, width = 12, units = "in", pointsize = 12, compression = c('lzw'), bg = "white", res = 600)
multiplot(plotlist = F1fRM.heat.ls, cols = 1, layout = matrix(data = c(1,2), nrow = 2, byrow = TRUE))
dev.off()

# Alternative Colour (rainbow) Version of Figure X for Revised Manuscript
#     colour.scale.type = 'rainbow', 


F1fRM.rainbow.ls <- list()

MAP.p.SpEr.LAR.F95.TS35.n500.rainbow.p <- pred.rast.plot(
    rast.df = MAP.p.SpEr.df,
    fill.name = '% SOC', 
    points.df = Response.df[,c('Easting', 'Northing','Carbon')],
    point.fill.column.name = 'Carbon',  
    colour.scale.type = 'rainbow', 
    auto.col.lims = TRUE,
    colour.scale.limits = c(0,7), 
    out.of.range.colour = 'white',
    overlay.points = TRUE,
    point.boundary.colour = 'black',  
    point.size = 2.5,
    plot.title = '(a)',
    title.text.size = 25,
    axis.text.size = 25,
    axis.tick.step = 100,
    x.axis.first.tick = 368100,
    x.axis.last.tick = 369800,
    y.axis.first.tick = 663200,
    y.axis.last.tick = 6634200) +
    annotate(geom = 'point', shape = 3, x = c(368842, 368944), y = c(6633649, 6633649), colour = 'black')

MAP.p.SpEr.LAR.F95.TS35.n500.rainbow.p

F1fRM.rainbow.ls[[1]] <- MAP.p.SpEr.LAR.F95.TS35.n500.rainbow.p

##

Rast.Uncert.LAR.F95.TS35.n500.rainbow.p <- pred.rast.plot(
    rast.df = data.frame(B1.CRS.df[,c('Easting','Northing')], Fill = Areal.MAP.LAR.F95.TS35.n500$Results[,'Int.Width']),
    fill.name = '% SOC', 
    points.df = Response.df, 
    point.fill.column.name = 'Carbon',  
    colour.scale.type = 'rainbow', 
    auto.col.lims = TRUE,
    colour.scale.limits = c(0,20),
    out.of.range.colour = 'white',
    overlay.points = TRUE,
    point.boundary.colour = 'black',
    point.size = 2.5,
    plot.title = '(b)',
    title.text.size = 25,
    axis.text.size = 25,
    axis.tick.step = 100,
    x.axis.first.tick = 368100,
    x.axis.last.tick = 369800,
    y.axis.first.tick = 663200,
    y.axis.last.tick = 6634200) +
    annotate(geom = 'point', shape = 3, x = c(368842, 368944, (368944 + 25*3) ), y = c(6633649, 6633649, (6633649 + (25*5))), colour = 'black')

Rast.Uncert.LAR.F95.TS35.n500.rainbow.p

F1fRM.rainbow.ls[[2]] <- Rast.Uncert.LAR.F95.TS35.n500.rainbow.p

pdf(file = 'Figure_3_Rainbow_Colours.pdf', height = 20, width = 12)
multiplot(plotlist = F1fRM.rainbow.ls, cols = 1, layout = matrix(data = c(1,2), nrow = 2, byrow = TRUE))
dev.off()

tiff(filename = 'Figure_3_Rainbow_Colours.tiff', height = 20, width = 12, units = "in", pointsize = 12, compression = c('lzw'), bg = "white", res = 600)
multiplot(plotlist = F1fRM.ls, cols = 1, layout = matrix(data = c(1,2), nrow = 2, byrow = TRUE))
dev.off()


################################
################################
##                            ##
##                            ##
##    Supplementary Figure    ##
##                            ##
##                            ##
################################
################################

## Comparing Validation Set Element Prediction Error Distributions:

VS = c( rep('F0.40 Exh', nrow(Exh.F40.TS35.n500$VSPE.mat)*ncol(Exh.F40.TS35.n500$VSPE.mat)),
        rep('F0.40 Fwd', nrow(Fwd.F40.TS35.n500$VSPE.mat)*ncol(Fwd.F40.TS35.n500$VSPE.mat)),
        rep('F0.40 Bwd', nrow(Bwd.F40.TS35.n500$VSPE.mat)*ncol(Bwd.F40.TS35.n500$VSPE.mat)),
        rep('F0.40 Seq', nrow(Seq.F40.TS35.n500$VSPE.mat)*ncol(Seq.F40.TS35.n500$VSPE.mat)),
        rep('F0.40 LAR', nrow(LAR.F40.TS35.n500$VSPE.mat)*ncol(LAR.F40.TS35.n500$VSPE.mat)),
        rep('F0.95 LAR', nrow(LAR.F95.TS35.n500$VSPE.mat)*ncol(LAR.F95.TS35.n500$VSPE.mat)))

VSPE.comp.df <- data.frame(VS = factor(x = VS, levels = c('F0.40 Fwd','F0.40 Bwd','F0.40 Seq','F0.40 Exh','F0.40 LAR','F0.95 LAR')), VSEPE = numeric(length = length(VS)))

VSPE.comp.df[VSPE.comp.df$VS == 'F0.40 Exh', 'VSEPE'] <- as.vector(Exh.F40.TS35.n500$VSPE.mat)

VSPE.comp.df[VSPE.comp.df$VS == 'F0.40 Fwd', 'VSEPE'] <- as.vector(Fwd.F40.TS35.n500$VSPE.mat)

VSPE.comp.df[VSPE.comp.df$VS == 'F0.40 Bwd', 'VSEPE'] <- as.vector(Bwd.F40.TS35.n500$VSPE.mat)

VSPE.comp.df[VSPE.comp.df$VS == 'F0.40 Seq', 'VSEPE'] <- as.vector(Seq.F40.TS35.n500$VSPE.mat)

VSPE.comp.df[VSPE.comp.df$VS == 'F0.40 LAR', 'VSEPE'] <- as.vector(LAR.F40.TS35.n500$VSPE.mat)

VSPE.comp.df[VSPE.comp.df$VS == 'F0.95 LAR', 'VSEPE'] <- as.vector(LAR.F95.TS35.n500$VSPE.mat)

VSEPE.comp.p <- ggplot(aes(x = VS, y = VSEPE), data = VSPE.comp.df)

VSEPE.comp.p + geom_boxplot() + coord_flip() 

VSEPE.comp.p <- VSEPE.comp.p +
        geom_jitter(position = position_jitter(height = 0), size = 0.5, colour = 'blue') +
        coord_flip() +
        geom_boxplot(fill = NA, outlier.size = 0, size = 1, colour = 'black') +
        labs(x = '', y = 'Validation Set Element Prediction Error') + #, title = 'Cross Validation Runs on 500 Unique Training Set Selections of 35 from 60 Observations'), y = 'Validation Set Element Prediction Error',  +,x = 'Subset Selection Method & Maximum Permitted Covariate Pair Correlation Magnitude'
        scale_y_continuous(breaks = seq(from = -4, to = 4, by = 0.5)) +
        theme(text = element_text(size = 10, colour = 'black'),
              axis.text.x = element_text(size = 10, angle = 0, colour = 'black', vjust = 0.5),
              axis.text.y = element_text(size = 10, angle = 0, colour = 'black'),
              legend.position = 'right',
              panel.background = element_rect(fill = 'white'),
              panel.grid.major = element_line(colour = 'grey'),
              panel.grid.minor = element_line(colour = 'grey'),
              axis.ticks = element_line(colour = 'grey'))

##############################
##############################
##                          ##
##                          ##
##   Supplementary Figure   ##
##                          ##
##                          ##
##############################
##############################

# A panel of chord diagrams depicting the covariate selection frequencies from
# each of the following combinations of variable selection algorithm and
# Maximum Permitted Correlation Coefficient Magnitude in the filtered design matrix:

## (a) Exhaustive Search Variable Selection with MPCCM = 0.4,
## (b) Sequential Replacement Variable Selection with MPCCM = 0.4,
## (c) Forwards Stepwise Variable Selection with MPCCM = 0.4,
## (d) Backwards Stepwise Variable Selection with MPCCM = 0.4,
## (e) Least Angle Regression with MPCCM = 0.4,
## (f) Least Angle Regression with MPCCM = 0.95.

# arraged like so:
#  (a)   (b)
#  (c)   (d)
#  (e)   (f)

Figure.Chord.Panel.a <- cd.plot(P4I2_UnFiltered_DM = n60.P4I2.DM$DM,
                     LAR_Results = Exh.F40.TS35.n500,
                     graph.title = c('a'),
                     Covariate.Lab.Size = 1,
                     Point.Size = 1,
                     Title.Size = 4,
                     Legend.Height.Width = unit(2,'mm'),
                     Legend.Position = 'right',
                     Line.Width = 0.1)

Figure.Chord.Panel.b <- cd.plot(P4I2_UnFiltered_DM = n60.P4I2.DM$DM,
                     LAR_Results = Seq.F40.TS35.n500,
                     graph.title = c('b'),
                     Covariate.Lab.Size = 1,
                     Point.Size = 1,
                     Title.Size = 4,
                     Legend.Height.Width = unit(2,'mm'),
                     Legend.Position = 'right',
                     Line.Width = 0.1)

Figure.Chord.Panel.c <- cd.plot(P4I2_UnFiltered_DM = n60.P4I2.DM$DM,
                     LAR_Results = Fwd.F40.TS35.n500,
                     graph.title = c('c'),
                     Covariate.Lab.Size = 1,
                     Point.Size = 1,
                     Title.Size = 4,
                     Legend.Height.Width = unit(2,'mm'),
                     Legend.Position = 'right',
                     Line.Width = 0.1)

Figure.Chord.Panel.d <- cd.plot(P4I2_UnFiltered_DM = n60.P4I2.DM$DM,
                     LAR_Results = Bwd.F40.TS35.n500,
                     graph.title = c('d'),
                     Covariate.Lab.Size = 1,
                     Point.Size = 1,
                     Title.Size = 4,
                     Legend.Height.Width = unit(2,'mm'),
                     Legend.Position = 'right',
                     Line.Width = 0.1)

Figure.Chord.Panel.e <- cd.plot(P4I2_UnFiltered_DM = n60.P4I2.DM$DM,
                     LAR_Results = LAR.F40.TS35.n500,
                     graph.title = c('e'),
                     Covariate.Lab.Size = 1,
                     Point.Size = 1,
                     Title.Size = 4,
                     Legend.Height.Width = unit(2,'mm'),
                     Legend.Position = 'right',
                     Line.Width = 0.1)

Figure.Chord.Panel.f <- cd.plot(P4I2_UnFiltered_DM = n60.P4I2.DM$DM,
                     LAR_Results = LAR.F95.TS35.n500,
                     graph.title = c('f'),
                     Covariate.Lab.Size = 1,
                     Point.Size = 1,
                     Title.Size = 4,
                     Legend.Height.Width = unit(2,'mm'),
                     Legend.Position = 'right',
                     Line.Width = 0.1)

Figure.Chord.Panel.ls <- vector(mode = 'list', length = 6)

Figure.Chord.Panel.ls[[1]] <- Figure.Chord.Panel.a
Figure.Chord.Panel.ls[[2]] <- Figure.Chord.Panel.b
Figure.Chord.Panel.ls[[3]] <- Figure.Chord.Panel.c
Figure.Chord.Panel.ls[[4]] <- Figure.Chord.Panel.d
Figure.Chord.Panel.ls[[5]] <- Figure.Chord.Panel.e
Figure.Chord.Panel.ls[[6]] <- Figure.Chord.Panel.f


multiplot(plotlist = Figure.Chord.Panel.ls, layout = matrix(data = c(1,2,3,4,5,6),nrow = 3,byrow = TRUE))

tiff(filename = 'Figure_Chord_Panel_16x20.tiff',
     width = 16,
     height = 20,
     units = "cm",
     pointsize = 12,
     compression = c('lzw'),
     bg = "white",
     res = 600)
multiplot(plotlist = Figure.Chord.Panel.ls, layout = matrix(data = c(1,2,3,4,5,6),nrow = 3,byrow = TRUE))
dev.off() 

pdf(file = 'Figure_Chord_Panel_6.3x7.9.pdf',
     width = 16/2.54,
     height = 20/2.54,
     bg = "white")
multiplot(plotlist = Figure.Chord.Panel.ls, layout = matrix(data = c(1,2,3,4,5,6),nrow = 3,byrow = TRUE))
dev.off() 

tiff(filename = 'Figure_Chord_Panel_dim_16x20.tiff',
     width = 16,
     height = 20,
     units = "cm",
     pointsize = 12,
     compression = c('lzw'),
     bg = "white",
     res = 600)
multiplot(plotlist = Figure.Chord.Panel.ls, layout = matrix(data = c(1,2,3,4,5,6),nrow = 3,byrow = TRUE))
dev.off()  

pdf(file = 'Figure_Chord_Panel_dim_6.3x7.9.pdf',
     width = 16/2.54,
     height = 20/2.54,
     bg = "white")
multiplot(plotlist = Figure.Chord.Panel.ls, layout = matrix(data = c(1,2,3,4,5,6),nrow = 3,byrow = TRUE))
dev.off()

tiff(filename = 'Figure_Chord_Panel_dim_7.5x8.tiff',
     width = 7.5,
     height = 8,
     units = "in",
     pointsize = 12,
     compression = c('lzw'),
     bg = "white",
     res = 600)
multiplot(plotlist = Figure.Chord.Panel.ls, layout = matrix(data = c(1,2,3,4,5,6),nrow = 3,byrow = TRUE))
dev.off()

pdf(file = 'Figure_Chord_Panel_dim_7.5x8.pdf',
     width = 7.5,
     height = 8,
     bg = "white")
multiplot(plotlist = Figure.Chord.Panel.ls, layout = matrix(data = c(1,2,3,4,5,6),nrow = 3,byrow = TRUE))
dev.off()


####################################################
####################################################
##                                                ##
##                                                ##
##              Supplementary Figure              ##
##                                                ##
##                                                ##
####################################################
####################################################

####################################################
#                                                  #
#                                                  #
#    3 by 2 panel                                  #
#                                                  #
#                                                  #
#   Distributions of the coefficient estimates     #
#   for the 6 most frequently selected covariate   #
#   terms from the Least Angle Regression          #
#   variable selection on the 500 unique           #
#   training sets of 35 observations constructed   #
#   from the full design matrix filtered to        #
#   enforce a maximum permitted correlation        #
#   coefficient magnitude between covariate        #
#   pairs of 0.95:                                 #
#                                                  #
#   (a) ECA.Nov$^4$                                #
#   (b) LSF$^3$                                    #
#   (c) DVI.May                                    #
#   (d) WI                                         #
#   (e) ECA.Feb:Slp                                #
#   (f) Mag.II:FPCI                                #
#                                                  #
#                                                  #
####################################################

summary(unlist(LAR.F95.TS35.n500$SC.ls))

names(unlist(LAR.F95.TS35.n500$SC.ls))

L.F95.T35.Sel.Cov <- data.frame(Covariate = names(unlist(LAR.F95.TS35.n500$SC.ls)),
                                Value = unlist(LAR.F95.TS35.n500$SC.ls))


head(L.F95.T35.Sel.Cov)

Figure.4.Panel.ls <- vector(mode = 'list', length = 6)

#   Figure 4
#   (a) ECA.Nov$^4$                                #

ECA.Nov.4.df <- L.F95.T35.Sel.Cov[L.F95.T35.Sel.Cov$Covariate == 'ECA.Nov.4',]

ECA.Nov.4.hist <- ggplot(aes(x = Value), data = ECA.Nov.4.df)

ECA.Nov.4.hist <- ECA.Nov.4.hist +
    geom_histogram(colour = 'black', fill = 'grey') +
    ggtitle('(a)') +
    ylab('Frequency') +
    theme(text = element_text(size = 25, colour = 'black'),
          axis.text.x = element_text(size = 25, angle = 0, colour = 'black', vjust = 0.5),
          axis.text.y = element_text(size = 25, angle = 0, colour = 'black'),
          legend.position = 'right',
          panel.background = element_rect(fill = 'white'),
          panel.grid.major = element_line(colour = 'black'),
          panel.grid.minor = element_line(colour = 'black'),
          axis.ticks = element_line(colour = 'black'))

Figure.4.Panel.ls[[1]] <- ECA.Nov.4.hist 

#   Figure 4
#   (b) LSF$^3$                                    #

LSF.3.df <- L.F95.T35.Sel.Cov[L.F95.T35.Sel.Cov$Covariate == 'LSF.3',]

LSF.3.hist <- ggplot(aes(x = Value), data = LSF.3.df)

LSF.3.hist <- LSF.3.hist +
    geom_histogram(colour = 'black', fill = 'grey') +
    ggtitle('(b)') +
    ylab('Frequency') +
    theme(text = element_text(size = 25, colour = 'black'),
          axis.text.x = element_text(size = 25, angle = 0, colour = 'black', vjust = 0.5),
          axis.text.y = element_text(size = 25, angle = 0, colour = 'black'),
          legend.position = 'right',
          panel.background = element_rect(fill = 'white'),
          panel.grid.major = element_line(colour = 'black'),
          panel.grid.minor = element_line(colour = 'black'),
          axis.ticks = element_line(colour = 'black'))

LSF.3.hist 

Figure.4.Panel.ls[[2]] <- LSF.3.hist 

#   Figure 4
#   (c) DVI.May                                    #

DVI.May.df <- L.F95.T35.Sel.Cov[L.F95.T35.Sel.Cov$Covariate == 'DVI.May',]

DVI.May.hist <- ggplot(aes(x = Value), data = DVI.May.df)

DVI.May.hist <- DVI.May.hist +
    geom_histogram(colour = 'black', fill = 'grey') +
    ggtitle('(c)') +
    ylab('Frequency') +
    theme(text = element_text(size = 25, colour = 'black'),
          axis.text.x = element_text(size = 25, angle = 0, colour = 'black', vjust = 0.5),
          axis.text.y = element_text(size = 25, angle = 0, colour = 'black'),
          legend.position = 'right',
          panel.background = element_rect(fill = 'white'),
          panel.grid.major = element_line(colour = 'black'),
          panel.grid.minor = element_line(colour = 'black'),
          axis.ticks = element_line(colour = 'black'))

DVI.May.hist 

Figure.4.Panel.ls[[3]] <- DVI.May.hist 

#   Figure 4
#   (d) WI                                         #

WI.df <- L.F95.T35.Sel.Cov[L.F95.T35.Sel.Cov$Covariate == 'WI',]

WI.hist <- ggplot(aes(x = Value), data = WI.df)

WI.hist <- WI.hist +
    geom_histogram(colour = 'black', fill = 'grey') +
    ggtitle('(d)') +
    ylab('Frequency') +
    theme(text = element_text(size = 25, colour = 'black'),
          axis.text.x = element_text(size = 25, angle = 0, colour = 'black', vjust = 0.5),
          axis.text.y = element_text(size = 25, angle = 0, colour = 'black'),
          legend.position = 'right',
          panel.background = element_rect(fill = 'white'),
          panel.grid.major = element_line(colour = 'black'),
          panel.grid.minor = element_line(colour = 'black'),
          axis.ticks = element_line(colour = 'black'))

WI.hist 

Figure.4.Panel.ls[[4]] <- WI.hist

#   Figure 4
#   (e) ECA.Feb:Slp                                #

ECA.Feb_X_Slp.df <- L.F95.T35.Sel.Cov[L.F95.T35.Sel.Cov$Covariate == 'ECA.Feb_X_Slp',]

ECA.Feb_X_Slp.hist <- ggplot(aes(x = Value), data = ECA.Feb_X_Slp.df)

ECA.Feb_X_Slp.hist <- ECA.Feb_X_Slp.hist +
    geom_histogram(colour = 'black', fill = 'grey') +
    ggtitle('(e)') +
    ylab('Frequency') +
    theme(text = element_text(size = 25, colour = 'black'),
          axis.text.x = element_text(size = 25, angle = 0, colour = 'black', vjust = 0.5),
          axis.text.y = element_text(size = 25, angle = 0, colour = 'black'),
          legend.position = 'right',
          panel.background = element_rect(fill = 'white'),
          panel.grid.major = element_line(colour = 'black'),
          panel.grid.minor = element_line(colour = 'black'),
          axis.ticks = element_line(colour = 'black'))

ECA.Feb_X_Slp.hist 

Figure.4.Panel.ls[[5]] <- ECA.Feb_X_Slp.hist

#   Figure 4
#   (f) Mag.II:FPCI                                #

Mag.II_X_FPCI.df <- L.F95.T35.Sel.Cov[L.F95.T35.Sel.Cov$Covariate == 'Mag.II_X_FPCI',]

Mag.II_X_FPCI.hist <- ggplot(aes(x = Value), data = Mag.II_X_FPCI.df)

Mag.II_X_FPCI.hist <- Mag.II_X_FPCI.hist +
    geom_histogram(colour = 'black', fill = 'grey') +
    ggtitle('(f)') +
    ylab('Frequency') +
    theme(text = element_text(size = 25, colour = 'black'),
          axis.text.x = element_text(size = 25, angle = 0, colour = 'black', vjust = 0.5),
          axis.text.y = element_text(size = 25, angle = 0, colour = 'black'),
          legend.position = 'right',
          panel.background = element_rect(fill = 'white'),
          panel.grid.major = element_line(colour = 'black'),
          panel.grid.minor = element_line(colour = 'black'),
          axis.ticks = element_line(colour = 'black'))

Mag.II_X_FPCI.hist 

Figure.4.Panel.ls[[6]] <- Mag.II_X_FPCI.hist 

multiplot(plotlist = Figure.4.Panel.ls, cols = 2, layout = matrix(data = c(1,2,3,4,5,6), nrow = 3, byrow = TRUE))

## 

pdf(file = 'Figure_4.pdf',
    height = 20,
    width = 19)
multiplot(plotlist = Figure.4.Panel.ls, cols = 2, layout = matrix(data = c(1,2,3,4,5,6), nrow = 3, byrow = TRUE))
dev.off()

tiff(filename = 'Figure_4.tiff',
     height = 20,
     width = 19,
     units = "in",
     pointsize = 12,
     compression = c('lzw'),
     bg = "white",
     res = 600)
multiplot(plotlist = Figure.4.Panel.ls, cols = 2, layout = matrix(data = c(1,2,3,4,5,6), nrow = 3, byrow = TRUE))
dev.off()


####################################################
####################################################
##                                                ##
##                                                ##
##  LAR on design matrices produced with variety  ##
##  of levels of stringency in the filtering out  ##
##  of correlated covariates and training set     ##
##  sizes                                         ##
##                                                ##
##                                                ##
####################################################
####################################################


###########################################################
#                                                         #
#                                                         #
#   Design Matrix Filtered to Enforce Maximum Permitted   #
#   Correlation Coefficient Magnitude between Covariate   #
#   Pairs of 0.40                                         #
#                                                         #
#                                                         #
###########################################################


 # a) Least Angle Regression Variable Selection on 500 unique 35 observation training sets constructed from a design matrix filtered to enforce a maximum permitted correlation coefficient magnitude between remaining covariate pairs of 0.4 

  # i) Run the variable selection within the cross validation scheme

system.time(LAR.F40.TS35.n500 <- lar.cv(Response = Response.df$Carbon,
                                        DM = n60.P4I2.DM.F40$DM, 
                                        CVTI = CVTI_TS35_VS25_n500)) # 46 sec

  # ii) calculate the model averaged predictions at the response sample locations

MAP.LAR.F40.TS35.n500 <- ma.pred(cotpa.mat = LAR.F40.TS35.n500$Design.mat,
                                vspe.mat = LAR.F40.TS35.n500$VSPE.mat,
                                sc.ls = LAR.F40.TS35.n500$SC.ls,
                                Response = LAR.F40.TS35.n500$Response,
                                CV.TI = LAR.F40.TS35.n500$CV.TI,
                                Response.TI.means = LAR.F40.TS35.n500$Response.TI.means,
                                DM.TI.C.means.store = LAR.F40.TS35.n500$DM.TI.C.means.store,
                                DM.TI.C.norms.store = LAR.F40.TS35.n500$DM.TI.C.norms.store,
                                plot.wgf = TRUE,
                                calc.R.2 = TRUE)

MAP.LAR.F40.TS35.n500$R2

 # b) Least Angle Regression Variable Selection on 500 unique 45 training sets constructed from a design matrix filtered to enforce a maximum permitted correlation coefficient magnitude between remaining covariate pairs of 0.4

  # i) Run the variable selection within the cross validation scheme

system.time(LAR.F40.TS45.n500 <- lar.cv(Response = Response.df$Carbon,
                                        DM = n60.P4I2.DM.F40$DM, 
                                        CVTI = CVTI_TS45_VS15_n500)) # 43 sec

  # ii) calculate the model averaged predictions at the response sample locations

MAP.LAR.F40.TS45.n500 <- ma.pred(cotpa.mat = LAR.F40.TS45.n500$Design.mat,
                                 vspe.mat = LAR.F40.TS45.n500$VSPE.mat,
                                 sc.ls = LAR.F40.TS45.n500$SC.ls,
                                 Response = LAR.F40.TS45.n500$Response,
                                 CV.TI = LAR.F40.TS45.n500$CV.TI,
                                 Response.TI.means = LAR.F40.TS45.n500$Response.TI.means,
                                 DM.TI.C.means.store = LAR.F40.TS45.n500$DM.TI.C.means.store,
                                 DM.TI.C.norms.store = LAR.F40.TS45.n500$DM.TI.C.norms.store,
                                 plot.wgf = TRUE,
                                 calc.R.2 = TRUE)

MAP.LAR.F40.TS45.n500$R2 # 0.4506636

 # c) Least Angle Regression Variable Selection on 500 unique training sets constructed from a design matrix filtered to enforce a maximum permitted correlation coefficient magnitude between remaining covariate pairs of 0.4

  # i) Run the variable selection within the cross validation scheme

system.time(LAR.F40.TS55.n500 <- lar.cv(Response = Response.df$Carbon,
                                        DM = n60.P4I2.DM.F40$DM, 
                                        CVTI = CVTI_TS55_VS5_n500)) # 42 sec

  # ii) calculate the model averaged predictions at the response sample locations

MAP.LAR.F40.TS55.n500 <- ma.pred(cotpa.mat = LAR.F40.TS55.n500$Design.mat,
                                 vspe.mat = LAR.F40.TS55.n500$VSPE.mat,
                                 sc.ls = LAR.F40.TS55.n500$SC.ls,
                                 Response = LAR.F40.TS55.n500$Response,
                                 CV.TI = LAR.F40.TS55.n500$CV.TI,
                                 Response.TI.means = LAR.F40.TS55.n500$Response.TI.means,
                                 DM.TI.C.means.store = LAR.F40.TS55.n500$DM.TI.C.means.store,
                                 DM.TI.C.norms.store = LAR.F40.TS55.n500$DM.TI.C.norms.store,
                                 plot.wgf = TRUE,
                                 calc.R.2 = TRUE)

MAP.LAR.F40.TS55.n500$R2 # 0.5593292

###########################################################
#                                                         #
#                                                         #
#   Design Matrix Filtered to Enforce Maximum Permitted   #
#   Correlation Coefficient Magnitude between Covariate   #
#   Pairs of 0.60                                         #
#                                                         #
#                                                         #
###########################################################


 # a) Least Angle Regression Variable Selection on 500 unique 35 observation training sets constructed from a design matrix filtered to enforce a maximum permitted correlation coefficient magnitude between remaining covariate pairs of 0.6 

  # i) Run the variable selection within the cross validation scheme

system.time(LAR.F60.TS35.n500 <- lar.cv(Response = Response.df$Carbon,
                                        DM = n60.P4I2.DM.F60$DM, 
                                        CVTI = CVTI_TS35_VS25_n500)) # 2 min 12 sec

  # ii) calculate the model averaged predictions at the response sample locations

MAP.LAR.F60.TS35.n500 <- ma.pred(cotpa.mat = LAR.F60.TS35.n500$Design.mat,
                                vspe.mat = LAR.F60.TS35.n500$VSPE.mat,
                                sc.ls = LAR.F60.TS35.n500$SC.ls,
                                Response = LAR.F60.TS35.n500$Response,
                                CV.TI = LAR.F60.TS35.n500$CV.TI,
                                Response.TI.means = LAR.F60.TS35.n500$Response.TI.means,
                                DM.TI.C.means.store = LAR.F60.TS35.n500$DM.TI.C.means.store,
                                DM.TI.C.norms.store = LAR.F60.TS35.n500$DM.TI.C.norms.store,
                                plot.wgf = TRUE,
                                calc.R.2 = TRUE)

MAP.LAR.F60.TS35.n500$R2

 # b) Least Angle Regression Variable Selection on 500 unique 45 training sets constructed from a design matrix filtered to enforce a maximum permitted correlation coefficient magnitude between remaining covariate pairs of 0.6

  # i) Run the variable selection within the cross validation scheme

system.time(LAR.F60.TS45.n500 <- lar.cv(Response = Response.df$Carbon,
                                        DM = n60.P4I2.DM.F60$DM, 
                                        CVTI = CVTI_TS45_VS15_n500)) # 2 mins 25 sec

  # ii) calculate the model averaged predictions at the response sample locations

MAP.LAR.F60.TS45.n500 <- ma.pred(cotpa.mat = LAR.F60.TS45.n500$Design.mat,
                                 vspe.mat = LAR.F60.TS45.n500$VSPE.mat,
                                 sc.ls = LAR.F60.TS45.n500$SC.ls,
                                 Response = LAR.F60.TS45.n500$Response,
                                 CV.TI = LAR.F60.TS45.n500$CV.TI,
                                 Response.TI.means = LAR.F60.TS45.n500$Response.TI.means,
                                 DM.TI.C.means.store = LAR.F60.TS45.n500$DM.TI.C.means.store,
                                 DM.TI.C.norms.store = LAR.F60.TS45.n500$DM.TI.C.norms.store,
                                 plot.wgf = TRUE,
                                 calc.R.2 = TRUE)

MAP.LAR.F60.TS45.n500$R2 # 

 # c) Least Angle Regression Variable Selection on 500 unique training sets constructed from a design matrix filtered to enforce a maximum permitted correlation coefficient magnitude between remaining covariate pairs of 0.6

  # i) Run the variable selection within the cross validation scheme

system.time(LAR.F60.TS55.n500 <- lar.cv(Response = Response.df$Carbon,
                                        DM = n60.P4I2.DM.F60$DM, 
                                        CVTI = CVTI_TS55_VS5_n500)) # 2 mins 38 sec

  # ii) calculate the model averaged predictions at the response sample locations

MAP.LAR.F60.TS55.n500 <- ma.pred(cotpa.mat = LAR.F60.TS55.n500$Design.mat,
                                 vspe.mat = LAR.F60.TS55.n500$VSPE.mat,
                                 sc.ls = LAR.F60.TS55.n500$SC.ls,
                                 Response = LAR.F60.TS55.n500$Response,
                                 CV.TI = LAR.F60.TS55.n500$CV.TI,
                                 Response.TI.means = LAR.F60.TS55.n500$Response.TI.means,
                                 DM.TI.C.means.store = LAR.F60.TS55.n500$DM.TI.C.means.store,
                                 DM.TI.C.norms.store = LAR.F60.TS55.n500$DM.TI.C.norms.store,
                                 plot.wgf = TRUE,
                                 calc.R.2 = TRUE)

MAP.LAR.F60.TS55.n500$R2 # 


###########################################################
#                                                         #
#                                                         #
#   Design Matrix Filtered to Enforce Maximum Permitted   #
#   Correlation Coefficient Magnitude between Covariate   #
#   Pairs of 0.80                                         #
#                                                         #
#                                                         #
###########################################################


 # a) Least Angle Regression Variable Selection on 500 unique 35 observation training sets constructed from a design matrix filtered to enforce a maximum permitted correlation coefficient magnitude between remaining covariate pairs of 0.8 

  # i) Run the variable selection within the cross validation scheme

system.time(LAR.F80.TS35.n500 <- lar.cv(Response = Response.df$Carbon,
                                        DM = n60.P4I2.DM.F80$DM, 
                                        CVTI = CVTI_TS35_VS25_n500)) # 7 min 39 sec

  # ii) calculate the model averaged predictions at the response sample locations

MAP.LAR.F80.TS35.n500 <- ma.pred(cotpa.mat = LAR.F80.TS35.n500$Design.mat,
                                vspe.mat = LAR.F80.TS35.n500$VSPE.mat,
                                sc.ls = LAR.F80.TS35.n500$SC.ls,
                                Response = LAR.F80.TS35.n500$Response,
                                CV.TI = LAR.F80.TS35.n500$CV.TI,
                                Response.TI.means = LAR.F80.TS35.n500$Response.TI.means,
                                DM.TI.C.means.store = LAR.F80.TS35.n500$DM.TI.C.means.store,
                                DM.TI.C.norms.store = LAR.F80.TS35.n500$DM.TI.C.norms.store,
                                plot.wgf = TRUE,
                                calc.R.2 = TRUE)

MAP.LAR.F80.TS35.n500$R2

 # b) Least Angle Regression Variable Selection on 500 unique 45 training sets constructed from a design matrix filtered to enforce a maximum permitted correlation coefficient magnitude between remaining covariate pairs of 0.8

  # i) Run the variable selection within the cross validation scheme

system.time(LAR.F80.TS45.n500 <- lar.cv(Response = Response.df$Carbon,
                                        DM = n60.P4I2.DM.F80$DM, 
                                        CVTI = CVTI_TS45_VS15_n500)) # 8 min 8 sec

  # ii) calculate the model averaged predictions at the response sample locations

MAP.LAR.F80.TS45.n500 <- ma.pred(cotpa.mat = LAR.F80.TS45.n500$Design.mat,
                                 vspe.mat = LAR.F80.TS45.n500$VSPE.mat,
                                 sc.ls = LAR.F80.TS45.n500$SC.ls,
                                 Response = LAR.F80.TS45.n500$Response,
                                 CV.TI = LAR.F80.TS45.n500$CV.TI,
                                 Response.TI.means = LAR.F80.TS45.n500$Response.TI.means,
                                 DM.TI.C.means.store = LAR.F80.TS45.n500$DM.TI.C.means.store,
                                 DM.TI.C.norms.store = LAR.F80.TS45.n500$DM.TI.C.norms.store,
                                 plot.wgf = TRUE,
                                 calc.R.2 = TRUE)

MAP.LAR.F80.TS45.n500$R2 # 

 # c) Least Angle Regression Variable Selection on 500 unique training sets constructed from a design matrix filtered to enforce a maximum permitted correlation coefficient magnitude between remaining covariate pairs of 0.8

  # i) Run the variable selection within the cross validation scheme

system.time(LAR.F80.TS55.n500 <- lar.cv(Response = Response.df$Carbon,
                                        DM = n60.P4I2.DM.F80$DM, 
                                        CVTI = CVTI_TS55_VS5_n500)) # 9 min 8 sec

  # ii) calculate the model averaged predictions at the response sample locations

MAP.LAR.F80.TS55.n500 <- ma.pred(cotpa.mat = LAR.F80.TS55.n500$Design.mat,
                                 vspe.mat = LAR.F80.TS55.n500$VSPE.mat,
                                 sc.ls = LAR.F80.TS55.n500$SC.ls,
                                 Response = LAR.F80.TS55.n500$Response,
                                 CV.TI = LAR.F80.TS55.n500$CV.TI,
                                 Response.TI.means = LAR.F80.TS55.n500$Response.TI.means,
                                 DM.TI.C.means.store = LAR.F80.TS55.n500$DM.TI.C.means.store,
                                 DM.TI.C.norms.store = LAR.F80.TS55.n500$DM.TI.C.norms.store,
                                 plot.wgf = TRUE,
                                 calc.R.2 = TRUE)

MAP.LAR.F80.TS55.n500$R2 #


###########################################################
#                                                         #
#                                                         #
#   Design Matrix Filtered to Enforce Maximum Permitted   #
#   Correlation Coefficient Magnitude between Covariate   #
#   Pairs of 0.95                                         #
#                                                         #
#                                                         #
###########################################################


 # a) Least Angle Regression Variable Selection on 500 unique 35 observation training sets constructed from a design matrix filtered to enforce a maximum permitted correlation coefficient magnitude between remaining covariate pairs of 0.95

  # i) Run the variable selection within the cross validation scheme

system.time(LAR.F95.TS35.n500 <- lar.cv(Response = Response.df$Carbon,
                                        DM = n60.P4I2.DM.F95$DM, 
                                        CVTI = CVTI_TS35_VS25_n500)) #  21 min 14 sec

  # ii) calculate the model averaged predictions at the response sample locations

MAP.LAR.F95.TS35.n500 <- ma.pred(cotpa.mat = LAR.F95.TS35.n500$Design.mat,
                                 vspe.mat = LAR.F95.TS35.n500$VSPE.mat,
                                 sc.ls = LAR.F95.TS35.n500$SC.ls,
                                 Response = LAR.F95.TS35.n500$Response,
                                 CV.TI = LAR.F95.TS35.n500$CV.TI,
                                 Response.TI.means = LAR.F95.TS35.n500$Response.TI.means,
                                 DM.TI.C.means.store = LAR.F95.TS35.n500$DM.TI.C.means.store,
                                 DM.TI.C.norms.store = LAR.F95.TS35.n500$DM.TI.C.norms.store,
                                 plot.wgf = TRUE,
                                 calc.R.2 = TRUE)

MAP.LAR.F95.TS35.n500$R2

 # b) Least Angle Regression Variable Selection on 500 unique 45 training sets constructed from a design matrix filtered to enforce a maximum permitted correlation coefficient magnitude between remaining covariate pairs of 0.8

  # i) Run the variable selection within the cross validation scheme

system.time(LAR.F95.TS45.n500 <- lar.cv(Response = Response.df$Carbon,
                                        DM = n60.P4I2.DM.F95$DM, 
                                        CVTI = CVTI_TS45_VS15_n500)) # 23 min 43 sec

  # ii) calculate the model averaged predictions at the response sample locations

MAP.LAR.F95.TS45.n500 <- ma.pred(cotpa.mat = LAR.F95.TS45.n500$Design.mat,
                                 vspe.mat = LAR.F95.TS45.n500$VSPE.mat,
                                 sc.ls = LAR.F95.TS45.n500$SC.ls,
                                 Response = LAR.F95.TS45.n500$Response,
                                 CV.TI = LAR.F95.TS45.n500$CV.TI,
                                 Response.TI.means = LAR.F95.TS45.n500$Response.TI.means,
                                 DM.TI.C.means.store = LAR.F95.TS45.n500$DM.TI.C.means.store,
                                 DM.TI.C.norms.store = LAR.F95.TS45.n500$DM.TI.C.norms.store,
                                 plot.wgf = TRUE,
                                 calc.R.2 = TRUE)

MAP.LAR.F95.TS45.n500$R2 # 

 # c) Least Angle Regression Variable Selection on 500 unique training sets constructed from a design matrix filtered to enforce a maximum permitted correlation coefficient magnitude between remaining covariate pairs of 0.8

  # i) Run the variable selection within the cross validation scheme

system.time(LAR.F95.TS55.n500 <- lar.cv(Response = Response.df$Carbon,
                                        DM = n60.P4I2.DM.F95$DM, 
                                        CVTI = CVTI_TS55_VS5_n500)) # min  1615.021  sec

  # ii) calculate the model averaged predictions at the response sample locations

MAP.LAR.F95.TS55.n500 <- ma.pred(cotpa.mat = LAR.F95.TS55.n500$Design.mat,
                                 vspe.mat = LAR.F95.TS55.n500$VSPE.mat,
                                 sc.ls = LAR.F95.TS55.n500$SC.ls,
                                 Response = LAR.F95.TS55.n500$Response,
                                 CV.TI = LAR.F95.TS55.n500$CV.TI,
                                 Response.TI.means = LAR.F95.TS55.n500$Response.TI.means,
                                 DM.TI.C.means.store = LAR.F95.TS55.n500$DM.TI.C.means.store,
                                 DM.TI.C.norms.store = LAR.F95.TS55.n500$DM.TI.C.norms.store,
                                 plot.wgf = TRUE,
                                 calc.R.2 = TRUE)

MAP.LAR.F95.TS55.n500$R2 # 


#############################
#############################
##                         ##
##                         ##
##   Supplementary Table   ##
##                         ##
##                         ##
#############################
#############################

Table.11 <- rbind(summary(abs(as.vector(LAR.F95.TS35.n500$VSPE.mat)))[1:6],
                  summary(abs(as.vector(LAR.F95.TS45.n500$VSPE.mat)))[1:6],
                  summary(abs(as.vector(LAR.F95.TS55.n500$VSPE.mat)))[1:6],                  

                  summary(abs(as.vector(LAR.F80.TS35.n500$VSPE.mat)))[1:6],
                  summary(abs(as.vector(LAR.F80.TS45.n500$VSPE.mat)))[1:6],
                  summary(abs(as.vector(LAR.F80.TS55.n500$VSPE.mat)))[1:6],                  

                  summary(abs(as.vector(LAR.F60.TS35.n500$VSPE.mat)))[1:6],
                  summary(abs(as.vector(LAR.F60.TS45.n500$VSPE.mat)))[1:6],
                  summary(abs(as.vector(LAR.F60.TS55.n500$VSPE.mat)))[1:6],                  

                  summary(abs(as.vector(LAR.F40.TS35.n500$VSPE.mat)))[1:6],
                  summary(abs(as.vector(LAR.F40.TS45.n500$VSPE.mat)))[1:6],
                  summary(abs(as.vector(LAR.F40.TS55.n500$VSPE.mat)))[1:6])

Table.11 <- data.frame(r.max = c(0.95, 0.95, 0.95, 0.8, 0.8, 0.8, 0.6, 0.6, 0.6, 0.4, 0.4, 0.4),
                       tss = c(35, 45, 55, 35, 45, 55, 35, 45, 55, 35, 45, 55),
                       Table.11)
                      
Table.11$R2 <- round(c(MAP.LAR.F95.TS35.n500$R2,
                      MAP.LAR.F95.TS45.n500$R2,
                      MAP.LAR.F95.TS55.n500$R2,                  

                      MAP.LAR.F80.TS35.n500$R2,
                      MAP.LAR.F80.TS45.n500$R2,
                      MAP.LAR.F80.TS55.n500$R2,                  

                      MAP.LAR.F60.TS35.n500$R2,
                      MAP.LAR.F60.TS45.n500$R2,
                      MAP.LAR.F60.TS55.n500$R2,                  

                      MAP.LAR.F40.TS35.n500$R2,
                      MAP.LAR.F40.TS45.n500$R2,
                      MAP.LAR.F40.TS55.n500$R2),
                      digits = 4)

Table.11

Table.11.xtab <- xtable(x = Table.11,
                        digits = 4,
                        display = c('g','g','g','g','g','g','g','g','g','g'),
       caption = 'Distribution Summary Statistics of Validation Set Element Prediction Errors from the Cross Validation Scheme and Coefficient of Determinations of the Model Averaged Predictions of the Entire Response Vector',
       label = 'Tab:11')

print(x = Table.11.xtab, include.rownames = FALSE)

############################
############################
##                        ##
##                        ##
##  Supplementary Figure  ##
##                        ##
##                        ##
############################
############################
#
################################################################################
#                                                                              #
#    S13 Figure 1 (a)                                                          #
#      - model-averaged predictions for each pixel in the full cover raster    #
#      - predictions from Least Angle Regression variable selection results    #
#        from the cross validation scheme cross validation scheme 500          #
#        unique 35 observation training                                        #
#      - training set design matrices constructed from the design matrix       #
#        filtered to enforce a maximum permitted correlation coefficient       #
#        magnitude between covariate pairs of 0.95                             #
#                                                                              #
################################################################################

S13.Figure.1.Panel.ls <- list()

Rast.LAR.F95.TS35.n500.rainbow.p <- pred.rast.plot(
    rast.df = data.frame(B1.CRS.df[,c('Easting','Northing')], Fill = Areal.MAP.LAR.F95.TS35.n500$Results[,'MA.Pred']),
    fill.name = '% SOC', 
    points.df = Response.df, 
    point.fill.column.name = 'Carbon',  
    colour.scale.type = 'rainbow', 
    auto.col.lims = TRUE,
    colour.scale.limits = c(0.5,6.5),
    out.of.range.colour = 'black',
    overlay.points = TRUE,
    point.boundary.colour = 'white',
    point.size = 2.5,
    plot.title = '(a)',
    title.text.size = 25,
    axis.text.size = 25,
    axis.tick.step = 100,
    x.axis.first.tick = 368100,
    x.axis.last.tick = 369800,
    y.axis.first.tick = 663200,
    y.axis.last.tick = 6634200)

Rast.LAR.F95.TS35.n500.rainbow.p

S13.Figure.1.Panel.ls[[1]] <- Rast.LAR.F95.TS35.n500.rainbow.p

################################################################################
#                                                                              #
#    S13 Figure 1 (b)                                                          #
#      - uncertainty estimated to accompany these predictions                  #
#                                                                              #  
#                                                                              #
################################################################################

# widths of intervals containing the middle 95% of the predicted values for each pixel
# from the cross validation scheme

Rast.Uncert.LAR.F95.TS35.n500.rainbow.p <- pred.rast.plot(
    rast.df = data.frame(B1.CRS.df[,c('Easting','Northing')], Fill = Areal.MAP.LAR.F95.TS35.n500$Results[,'Int.Width']),
    fill.name = '% SOC', 
    points.df = Response.df, 
    point.fill.column.name = 'Carbon',  
    colour.scale.type = 'rainbow', 
    auto.col.lims = TRUE,
    colour.scale.limits = c(0,20),
    out.of.range.colour = 'black',
    overlay.points = TRUE,
    point.boundary.colour = 'white',
    point.size = 2.5,
    plot.title = '(b)',
    title.text.size = 25,
    axis.text.size = 25,
    axis.tick.step = 100,
    x.axis.first.tick = 368100,
    x.axis.last.tick = 369800,
    y.axis.first.tick = 663200,
    y.axis.last.tick = 6634200)

Rast.Uncert.LAR.F95.TS35.n500.rainbow.p

S13.Figure.1.Panel.ls[[2]] <- Rast.Uncert.LAR.F95.TS35.n500.rainbow.p

pdf(file = 'S13_Figure_1.pdf', height = 20, width = 12)
multiplot(plotlist = S13.Figure.1.Panel.ls, cols = 1, layout = matrix(data = c(1,2), nrow = 2, byrow = TRUE))
dev.off()

tiff(filename = 'S13_Figure_1.tiff', height = 20, width = 12, units = "in", pointsize = 12, compression = c('lzw'), bg = "white", res = 600)
multiplot(plotlist = S13.Figure.1.Panel.ls, cols = 1, layout = matrix(data = c(1,2), nrow = 2, byrow = TRUE))
dev.off()

############################
############################
##                        ##
##                        ##
##  Supplementary Figure  ##
##                        ##
##                        ##
############################
############################


################################################################################
#                                                                              #
#    S14 Figure 1 (a)                                                          #
#      - model-averaged predictions for each pixel in the full cover raster    #
#      - predictions from Least Angle Regression variable selection results    #
#        from the cross validation scheme cross validation scheme 500          #
#        unique 35 observation training                                        #
#      - training set design matrices constructed from the design matrix       #
#        filtered to enforce a maximum permitted correlation coefficient       #
#        magnitude between covariate pairs of 0.95                             #
#                                                                              #
################################################################################

S14.Figure.1.Panel.ls <- list()

Rast.LAR.F95.TS35.n500.grey.p <- pred.rast.plot(
    rast.df = data.frame(B1.CRS.df[,c('Easting','Northing')], Fill = Areal.MAP.LAR.F95.TS35.n500$Results[,'MA.Pred']),
    fill.name = '% SOC', 
    points.df = Response.df, 
    point.fill.column.name = 'Carbon',  
    colour.scale.type = 'grey', 
    auto.col.lims = TRUE,
    colour.scale.limits = c(0.5,6.5),
    out.of.range.colour = 'red',
    overlay.points = TRUE,
    point.boundary.colour = 'white',
    point.size = 2.5,
    plot.title = '(a)',
    title.text.size = 25,
    axis.text.size = 25,
    axis.tick.step = 100,
    x.axis.first.tick = 368100,
    x.axis.last.tick = 369800,
    y.axis.first.tick = 663200,
    y.axis.last.tick = 6634200)

Rast.LAR.F95.TS35.n500.grey.p

S14.Figure.1.Panel.ls[[1]] <- Rast.LAR.F95.TS35.n500.grey.p

################################################################################
#                                                                              #
#    S14 Figure 1 (b)                                                          #
#      - uncertainty estimated to accompany these predictions                  #
#                                                                              #  
#                                                                              #
################################################################################

# widths of intervals containing the middle 95% of the predicted values for each pixel
# from the cross validation scheme

Rast.Uncert.LAR.F95.TS35.n500.grey.p <- pred.rast.plot(
    rast.df = data.frame(B1.CRS.df[,c('Easting','Northing')], Fill = Areal.MAP.LAR.F95.TS35.n500$Results[,'Int.Width']),
    fill.name = '% SOC', 
    points.df = Response.df, 
    point.fill.column.name = 'Carbon',  
    colour.scale.type = 'grey', 
    auto.col.lims = TRUE,
    colour.scale.limits = c(0,20),
    out.of.range.colour = 'red',
    overlay.points = TRUE,
    point.boundary.colour = 'white',
    point.size = 2.5,
    plot.title = '(b)',
    title.text.size = 25,
    axis.text.size = 25,
    axis.tick.step = 100,
    x.axis.first.tick = 368100,
    x.axis.last.tick = 369800,
    y.axis.first.tick = 663200,
    y.axis.last.tick = 6634200)

Rast.Uncert.LAR.F95.TS35.n500.grey.p

S14.Figure.1.Panel.ls[[2]] <- Rast.Uncert.LAR.F95.TS35.n500.grey.p

pdf(file = 'S14_Figure_1.pdf', height = 20, width = 12)
multiplot(plotlist = S14.Figure.1.Panel.ls, cols = 1, layout = matrix(data = c(1,2), nrow = 2, byrow = TRUE))
dev.off()

tiff(filename = 'S14_Figure_1.tiff', height = 20, width = 12, units = "in", pointsize = 12, compression = c('lzw'), bg = "white", res = 600)
multiplot(plotlist = S14.Figure.1.Panel.ls, cols = 1, layout = matrix(data = c(1,2), nrow = 2, byrow = TRUE))
dev.off()

############################
############################
##                        ##
##                        ##
##  Supplementary Figure  ##
##                        ##
##                        ##
############################
############################

# The frequencies of spatial positioning covariate term selection across the
# 500 models selected by Least Angle Regression conducted to model to the
# spatial component of the residuals from the covariate based predictions.
# Polynomial terms are indicated by the numbers following the term name
# (e.g. Easting.2 is the quadratic term for Easting). Interaction terms for a
# pair of terms are represented by the two covariate term names separated
# by ‘ X ’ (e.g. the interaction between Easting and Northing would be
# Easting X Northing).

S19.df <- data.frame(Term = names(unlist(SpEr.LAR.F95.TS35.n500$SC.ls)))

S19.df.sort <- data.frame(Term = levels(S19.df$Term), Frequency = numeric(length = length(levels(S19.df$Term))))
for(i in 1:nrow(S19.df.sort)){
    S19.df.sort[i,'Frequency'] <- length(S19.df[S19.df$Term == S19.df.sort[i,'Term'],])}
S19.df.sorted <- S19.df.sort[order(S19.df.sort$Frequency, decreasing = TRUE),]
S19.df.sorted$Rank <- 1:nrow(S19.df.sorted)

###

xbreak <- seq(from = 1, to = nrow(S19.df.sorted), by = 1)
xlab = S19.df.sorted$Term
length(xbreak) == length(xlab)

S19.Figure.1.p <- ggplot(aes(xmin = Rank-0.5, xmax = Rank+0.5, ymin = 0, ymax = Frequency), data = S19.df.sorted) +
    geom_rect(fill = 'black', colour = 'white') +
    scale_x_continuous(breaks = seq(from = 0, to = 60, by = 10)) +
    theme(text = element_text(size = 20, colour = 'black'),
          axis.text.x = element_text(size = 20, colour = 'black', vjust = 0.5, angle = 90),
          axis.text.y = element_text(size = 20, angle = 0, colour = 'black'),
          legend.position = 'right',
          legend.key.height = unit(x = 1.4, units = 'inches'),
          panel.background = element_rect(fill = 'white'),
          panel.grid.major = element_line(colour = 'darkgrey'),
          panel.grid.minor = element_line(colour = 'grey'),
          axis.ticks = element_line(colour = 'darkgrey')) +
   scale_x_continuous(breaks = xbreak, labels = xlab) + labs(y = 'Selection Frequency') #title = 'Spatial Coordinate Term Selection Frequency from 500 CV LAR UF TSS 35 Runs Selecting Spatial Polynomial Terms to\nPredict Residuals from the Model Averaged Predictions of Soil Carbon')              
S19.Figure.1.p

ggsave(plot = S19.Figure.1.p, file = 'S19_Figure_1_Revised.pdf', width = 19, height = 10, units = 'in') 

ggsave(plot = S19.Figure.1.p, file = 'S19_Figure_1/S19_Figure_1_Revised.tiff', width = 19, height = 10, units = 'in', dpi = 600) 

############################
############################
##                        ##
##                        ##
##     Striking Image     ##
##                        ##
##                        ##
############################
############################

Striking.Image <- cd.plot(P4I2_UnFiltered_DM = n60.P4I2.DM$DM,
                          LAR_Results = LAR.F95.TS35.n500,
                          graph.title = NULL,
                          Covariate.Lab.Size = 1,
                          Point.Size = 1,
                          Title.Size = 4,
                          Legend.Height.Width = unit(2,'mm'),
                          Legend.Position = 'right',
                          Line.Width = 0.1)

Striking.Image

ggsave(plot = Striking.Image, file = 'Striking_Image.tiff', dpi = 600, height = 8, width = 8, units = 'cm')

ggsave(plot = Striking.Image, file = 'Striking_Image.pdf', height = 8, width = 8, units = 'cm')

################################################################################
################################################################################
##                                                                            ## 
##   End of Supplementary Tables and Figures                                  ##
##                                                                            ##
################################################################################
################################################################################

