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

cd.plot <- function(P4I2_UnFiltered_DM = B1_P4I2_UnFiltered_DM,
                    LAR_Results = LAR_F95_TS35_n500_V2,
                    graph.title = c('CV LAR Selected Covariates from C95 Filtered P4I2'),
                    Covariate.Lab.Size = 4,
                    Point.Size = 4,
                    Title.Size = 4,
                    Legend.Height.Width = unit(4,'mm'),
                    Legend.Position = 'right',
                    Line.Width = 1){
    require('ggplot2')
    require('lars')
    require('grid')
    # so I need 63 divisions around the circle of the chord diagram
    Coef.ls = LAR_Results$SC.ls
    Coef.Names.char = character()
    for(i in 1:length(Coef.ls)){ Coef.Names.char = c(Coef.Names.char, names(Coef.ls[[i]]))}
    Uniq.Coef.Names.char = unique(Coef.Names.char)
    sort(Uniq.Coef.Names.char)
    Uniq.Int.Splt.Coef.Names.char = unique(unlist(strsplit(unique(Coef.Names.char), "_X_")))
    sort(Uniq.Int.Splt.Coef.Names.char)
    Selected.Covariates = data.frame(Covariate = factor(Coef.Names.char))
    FSC = data.frame(Covariate = levels(Selected.Covariates$Covariate), Frequency = numeric(length(levels(Selected.Covariates$Covariate))))
    for(i in 1:nrow(FSC)){FSC[i,'Frequency'] = length(which(FSC[i,'Covariate'] == Coef.Names.char))}
    FSCS = FSC[order(c(FSC$Frequency), decreasing = TRUE), ]
    FSCS$Rank = 1:nrow(FSCS)
    cn = c('NIR.Feb', 'RED.Feb', 'ECA.Feb', 'NIR.May', 'RED.May', 'ECA.May', 'NIR.Nov', 'RED.Nov', 'ECA.Nov', 'SR.Feb', 'DVI.Feb', 'NDVI.Feb', 'SAVI.Feb', 'NLVI.Feb', 'MNLVI.Feb', 'MSR.Feb', 'TVI.Feb', 'RDVI.Feb', 'SR.May', 'DVI.May', 'NDVI.May', 'SAVI.May', 'NLVI.May', 'MNLVI.May', 'MSR.May', 'TVI.May', 'RDVI.May', 'SR.Nov', 'DVI.Nov', 'NDVI.Nov', 'SAVI.Nov', 'NLVI.Nov', 'MNLVI.Nov', 'MSR.Nov', 'TVI.Nov', 'RDVI.Nov', 'CosAsp', 'CatAr', 'CatHe', 'CatSl', 'Elev', 'LSF', 'PlanC', 'ProfC', 'SVF', 'Slp', 'SPI', 'TRI', 'TPI', 'VTR', 'VS', 'WI', 'Mag.I', 'Mag.II', 'Mag.III', 'Mag.IV', 'Mag.V', 'Mag.VI', 'Radio.K', 'Radio.U', 'Radio.Th', 'FPCII', 'FPCI')
    Acron.G1 = data.frame(Covariate = cn[grep(x = cn, pattern = 'Feb')],
                         Acron = unlist(strsplit(x = cn[grep(x = cn, pattern = 'Feb')], split = '.Feb')),
                         Group = factor(x = rep('Feb', times = length(grep(x = cn, pattern = 'Feb'))), levels = c('Feb','May','Nov','DEM','Mag','Radio','FPC')))
    Acron.G2 = data.frame(Covariate = cn[grep(x = cn, pattern = 'May')],
                          Acron = unlist(strsplit(x = cn[grep(x = cn, pattern = 'May')], split = '.May')),
                         Group = factor(x = rep('May', times = length(grep(x = cn, pattern = 'May'))), levels = c('Feb','May','Nov','DEM','Mag','Radio','FPC')))
    Acron.G3 = data.frame(Covariate = cn[grep(x = cn, pattern = 'Nov')],
                          Acron = unlist(strsplit(x = cn[grep(x = cn, pattern = 'Nov')], split = '.Nov')),
                          Group = factor(x = rep('Nov', times = length(grep(x = cn, pattern = 'Nov'))), levels = c('Feb','May','Nov','DEM','Mag','Radio','FPC')))
    Acron.G4 = data.frame(Covariate = c('CosAsp', 'CatAr', 'CatHe', 'CatSl', 'Elev', 'LSF', 'PlanC', 'ProfC', 'SVF', 'Slp', 'SPI', 'TRI', 'TPI', 'VTR', 'VS', 'WI'),
                          Acron = c('CosAsp', 'CatAr', 'CatHe', 'CatSl', 'Elev', 'LSF', 'PlanC', 'ProfC', 'SVF', 'Slp', 'SPI', 'TRI', 'TPI', 'VTR', 'VS', 'WI'),
                          Group = factor(x = rep('DEM', times = length(c('CosAsp', 'CatAr', 'CatHe', 'CatSl', 'Elev', 'LSF', 'PlanC', 'ProfC', 'SVF', 'Slp', 'SPI', 'TRI', 'TPI', 'VTR', 'VS', 'WI')), levels = c('Feb','May','Nov','DEM','Mag','Radio','FPC'))))
    Acron.G5 = data.frame(Covariate = cn[grep(x = cn, pattern = 'Mag')],
                          Acron = cn[grep(x = cn, pattern = 'Mag')],
                          Group = factor(x = rep('Mag', times = length(grep(x = cn, pattern = 'Mag'))), levels = c('Feb','May','Nov','DEM','Mag','Radio','FPC')))
    Acron.G5$Acron = c('MagI', 'MagII', 'MagIII', 'MagIV', 'MagV', 'MagVI')
    Acron.G6 = data.frame(Covariate = cn[grep(x = cn, pattern = 'Radio')],
                          Acron = c('K','U','Th'),
                          Group = factor(x = rep('Radio', times = 3), levels = c('Feb','May','Nov','DEM','Mag','Radio','FPC')))
    Acron.G7 = data.frame(Covariate = cn[grep(x = cn, pattern = 'FPC')],
                          Acron = cn[grep(x = cn, pattern = 'FPC')],
                          Group = factor(x = rep('FPC', times = 2), levels = c('Feb','May','Nov','DEM','Mag','Radio','FPC')))
    Acron = rbind(Acron.G1, Acron.G2, Acron.G3, Acron.G4, Acron.G5, Acron.G6, Acron.G7)
    colnames(Acron) = c("Covariate", "Acron", "Label.Group" )
    P4.I2.Acron = data.frame(Full = colnames(P4I2_UnFiltered_DM),
                             Type = factor(x = rep(NA, ncol(P4I2_UnFiltered_DM)), levels = c('Linear', 'Polynomial', 'Interaction')),
                             Var1 = factor(x = rep(NA, ncol(P4I2_UnFiltered_DM)), levels = unique(Acron$Acron)),
                             Var2 = factor(x = rep(NA, ncol(P4I2_UnFiltered_DM)), levels = unique(Acron$Acron)),
                             Poly = numeric(length = ncol(P4I2_UnFiltered_DM)))
    P4.I2.Acron[grep(pattern = '_X_', x = P4.I2.Acron$Full), 'Type'] = 'Interaction'
    P4.I2.Acron[grep(pattern = '2|3|4', x = P4.I2.Acron$Full), 'Type'] = 'Polynomial'
    P4.I2.Acron[-grep(pattern = '_X_|2|3|4', x = P4.I2.Acron$Full), 'Type'] = 'Linear'
    for(i in 1:nrow(P4.I2.Acron)){
        if(P4.I2.Acron[i,'Type'] == 'Linear'){
            P4.I2.Acron[i,'Var1'] =  paste(Acron[Acron$Covariate == paste(P4.I2.Acron[i,'Full']), 'Acron'])}
        if(P4.I2.Acron[i,'Type'] == 'Polynomial'){
            V1 = paste(strsplit(x = paste(P4.I2.Acron[i, 'Full']), split = '[.]2|[.]3|[.]4'))
            VFL = unlist(strsplit(x = paste(P4.I2.Acron[i, 'Full']), split = '[.]'))
            P4.I2.Acron[i,'Var1'] = Acron[Acron$Covariate == V1,'Acron']
            P4.I2.Acron[i,'Poly'] = as.numeric(VFL[grep(pattern = '2|3|4', x = VFL)])}
        if(P4.I2.Acron[i,'Type'] == 'Interaction'){
            VFL = unlist(strsplit(x = paste(P4.I2.Acron[i, 'Full']), split = '_X_'))
            P4.I2.Acron[i,'Var1'] = Acron[Acron$Covariate == VFL[1],'Acron']
            P4.I2.Acron[i,'Var2'] = Acron[Acron$Covariate == VFL[2],'Acron']}}
    Acron.b = character(length = nrow(P4.I2.Acron))
    for(i in 1:nrow(P4.I2.Acron)){
        if(P4.I2.Acron[i,'Type'] == 'Linear'){
            Acron.b[i] = paste(P4.I2.Acron[i,'Var1'])}
        if(P4.I2.Acron[i,'Type'] == 'Polynomial'){
            Acron.b[i] = paste(P4.I2.Acron[i,'Var1'], P4.I2.Acron[i,'Poly'], sep = '.')}
        if(P4.I2.Acron[i,'Type'] == 'Interaction'){
            Acron.b[i] = paste(P4.I2.Acron[i,'Var1'], ':', P4.I2.Acron[i,'Var2'], sep = '')}}
    P4.I2.Acron$Acron = Acron.b
    poincare_segment = function(u1, u2, v1, v2, resolution = 100){
        # Courtesy of Vincent Zoonekynd via
        # http://stackoverflow.com/questions/14599150/chord-diagram-in-r
        # Check that the points are sufficiently different
        if( abs(u1-v1) < 1e-6 && abs(u2-v2) < 1e-6 )
            return( list(x=c(u1,v1), y=c(u2,v2)) )
        # Check that we are in the circle
        stopifnot( u1^2 + u2^2 - 1 <= 1e-6 )
        stopifnot( v1^2 + v2^2 - 1 <= 1e-6 )
        # Check it is not a diameter
        if( abs( u1*v2 - u2*v1 ) < 1e-6 )
            return( list(x=c(u1,v1), y=c(u2,v2)) )
        # Equation of the line: x^2 + y^2 + ax + by + 1 = 0 (circles orthogonal to the unit circle)
        a = ( u2 * (v1^2+v2^2) - v2 * (u1^2+u2^2) + u2 - v2 ) / ( u1*v2 - u2*v1 )
        b = ( u1 * (v1^2+v2^2) - v1 * (u1^2+u2^2) + u1 - v1 ) / ( u2*v1 - u1*v2 ) # Swap 1's and 2's
        # Center and radius of the circle
        cx = -a/2
        cy = -b/2
        radius = sqrt( (a^2+b^2)/4 - 1 )
        # Which portion of the circle should we draw?
        theta1 = atan2( u2-cy, u1-cx )
        theta2 = atan2( v2-cy, v1-cx )
        if( theta2 - theta1 > pi )
            theta2 = theta2 - 2 * pi
        else if( theta2 - theta1 < - pi )
            theta2 = theta2 + 2 * pi
        theta = seq( theta1, theta2, length = resolution)
        x = cx + radius * cos( theta )
        y = cy + radius * sin( theta )
        return(data.frame(x=x, y=y))}
    B1.CCP = data.frame(FSCS,
                       Interaction = logical(length = nrow(FSCS)),
                       Covariate1 = factor(x = rep(NA,nrow(FSCS)), levels = cn),
                       Covariate2 = factor(x = rep(NA,nrow(FSCS)), levels = cn),
                       Polynomial = logical(length = nrow(FSCS)),
                       Poly.Covariate = factor(x = rep(NA,nrow(FSCS)), levels = cn),
                       Poly.Num = rep(1, nrow(FSCS)),
                       Angle = numeric(length = nrow(FSCS)))
    for(i in 1:nrow(B1.CCP)){
        if(length(grep(pattern = '_X_', x = B1.CCP[i,'Covariate'])) > 0){
            B1.CCP[i,'Interaction'] = TRUE
            B1.CCP[i,'Covariate1'] = strsplit(x = paste(B1.CCP[i,'Covariate']), split = '_X_')[[1]][1]
            B1.CCP[i,'Covariate2'] = strsplit(x = paste(B1.CCP[i,'Covariate']), split = '_X_')[[1]][2]}
        if(length(grep(pattern = '2', x = paste(B1.CCP[i,'Covariate']))) > 0){
            B1.CCP[i,'Polynomial'] = TRUE
            B1.CCP[i,'Poly.Covariate'] = strsplit(x = paste(B1.CCP[i,'Covariate']), split = '.2')
            B1.CCP[i,'Poly.Num'] = 2}
        if(length(grep(pattern = '3', x = paste(B1.CCP[i,'Covariate']))) > 0){
            B1.CCP[i,'Polynomial'] = TRUE
            B1.CCP[i,'Poly.Covariate'] = strsplit(x = paste(B1.CCP[i,'Covariate']), split = '.3')
            B1.CCP[i,'Poly.Num'] = 3}
        if(length(grep(pattern = '4', x = paste(B1.CCP[i,'Covariate']))) > 0){
            B1.CCP[i,'Polynomial'] = TRUE
            B1.CCP[i,'Poly.Covariate'] = strsplit(x = paste(B1.CCP[i,'Covariate']), split = '.4')
            B1.CCP[i,'Poly.Num'] = 4}
        if(length(grep(pattern = '4|3|2|_X_', x = paste(B1.CCP[i,'Covariate']))) == 0){
            B1.CCP[i,'Poly.Covariate'] = B1.CCP[i,'Covariate']}
        if(B1.CCP[i,'Interaction'] == FALSE){
            B1.CCP[i,'Angle'] = which(B1.CCP[i,'Poly.Covariate'] == Acron$Covariate)}} # position around the circle divided into 63 segments
    Lines = B1.CCP[B1.CCP$Interaction == 'TRUE',]
    Lines = data.frame(Lines, Angle1 = numeric(nrow(Lines)), Angle2 = numeric(nrow(Lines)), Angle.Mid = numeric(nrow(Lines)))
    for(i in 1:nrow(Lines)){
        Lines[i,'Angle1'] = which(Lines[i,'Covariate1'] == Acron$Covariate)
        Lines[i,'Angle2'] = which(Lines[i,'Covariate2'] == Acron$Covariate)
        Lines[i,'Angle.Mid'] = mean(c(Lines[i,'Angle1'], Lines[i,'Angle2']))}
    angle = seq(from = (2*pi/(2*63)), to = ((2*63-1)*2*pi/(2*63)), by = (2*pi/63))
    UC.63p = data.frame(angle, x = cos(pi/2 - angle), y = sin(pi/2 - angle))
    CNP = data.frame(UC.63p, Acron, Colour = 1:63)
    npoints = 100
    r = nrow(Lines)*npoints
    I.Store = data.frame(x = numeric(),
                         y = numeric(),
                         Interaction.Number = numeric(),
                         Frequency = numeric())
    counter = 0
    for(i in 1:nrow(Lines)){
        counter = counter + 1
        c1.x = CNP[which(CNP[,'Covariate'] == Lines[i,'Covariate1']), 'x']
        c1.y = CNP[which(CNP[,'Covariate'] == Lines[i,'Covariate1']), 'y']
        c2.x = CNP[which(CNP[,'Covariate'] == Lines[i,'Covariate2']), 'x']
        c2.y = CNP[which(CNP[,'Covariate'] == Lines[i,'Covariate2']), 'y']
        Frequency.i = Lines[i,'Frequency']
        if(abs((which(CNP[,'Covariate'] == Lines[i,'Covariate1']) - which(CNP[,'Covariate'] == Lines[i,'Covariate2']))) == 63/2){ #
            # i.e. line will be a diameter line from point to point
            if(c1.x < c2.x){
                gradient = (c2.y - c1.y)/(c2.x - c1.x)
                x.path = seq(from = c1.x, to = c2.x, length.out = npoints)
                y.path = gradient*x.path} else{ # all diameters pass through the center of the circle at (0,0)
                    gradient = (c1.y - c2.y)/(c1.x - c2.x)
                    x.path = seq(from = c2.x, to = c1.x, length.out = npoints)
                    y.path = gradient*x.path}
            I1 = data.frame(x = x.path, y = y.path)} else{
                I1 = poincare_segment(u1 = c1.x, u2 = c1.y, v1 = c2.x, v2 = c2.y, resolution = npoints)}
        ID = data.frame(I1,Interaction.Number = i)
        I.Store[(1+(counter-1)*npoints):(counter*npoints),'x'] = ID$x
        I.Store[(1+(counter-1)*npoints):(counter*npoints),'y'] = ID$y
        I.Store[(1+(counter-1)*npoints):(counter*npoints),'Interaction.Number'] = ID$Interaction.Number
        I.Store[(1+(counter-1)*npoints):(counter*npoints),'Frequency'] = Frequency.i
        print(c(paste(round(npoints*i/nrow(Lines),2)),'%', 'Complete'), quote = FALSE, sep = '')}
       ## Points.xy
    Points = B1.CCP[B1.CCP$Interaction == 'FALSE',]
    Points.xy = data.frame(Points, x = numeric(nrow(Points)), y = numeric(nrow(Points)))
    for(i in 1:nrow(Points.xy)){
        Points.xy[i, 'x'] =  CNP[CNP$Covariate == Points.xy[i, 'Poly.Covariate'], 'x']
        Points.xy[i, 'y'] =  CNP[CNP$Covariate == Points.xy[i, 'Poly.Covariate'], 'y']}
    Points.xy = data.frame(Points.xy, Interaction.Number = rep(0,nrow(Points.xy)))
    ## Arcs
    Feb.Angle = seq(from = min(CNP[CNP$Label.Group == 'Feb', 'angle']),
                    to = max(CNP[CNP$Label.Group == 'Feb', 'angle']),
                    length.out = 100)
    Feb.Arc = data.frame(x = cos(pi/2 - Feb.Angle), y = sin(pi/2 - Feb.Angle), Label = rep('Feb',length(Feb.Angle)))
    May.Angle = seq(from = min(CNP[CNP$Label.Group == 'May', 'angle']),
                    to = max(CNP[CNP$Label.Group == 'May', 'angle']),
                    length.out = 100)
    May.Arc = data.frame(x = cos(pi/2 - May.Angle), y = sin(pi/2 - May.Angle), Label = rep('May',length(May.Angle)))
    Nov.Angle = seq(from = min(CNP[CNP$Label.Group == 'Nov', 'angle']),
                    to = max(CNP[CNP$Label.Group == 'Nov', 'angle']),
                    length.out = 100)
    Nov.Arc = data.frame(x = cos(pi/2 - Nov.Angle), y = sin(pi/2 - Nov.Angle), Label = rep('Nov',length(Nov.Angle)))
    DEM.Angle = seq(from = min(CNP[CNP$Label.Group == 'DEM', 'angle']),
                    to = max(CNP[CNP$Label.Group == 'DEM', 'angle']),
                    length.out = 100)
    DEM.Arc = data.frame(x = cos(pi/2 - DEM.Angle), y = sin(pi/2 - DEM.Angle), Label = rep('DEM',length(DEM.Angle)))
    Mag.Angle = seq(from = min(CNP[CNP$Label.Group == 'Mag', 'angle']),
                    to = max(CNP[CNP$Label.Group == 'Mag', 'angle']),
                    length.out = 100)
    Mag.Arc = data.frame(x = cos(pi/2 - Mag.Angle), y = sin(pi/2 - Mag.Angle), Label = rep('Mag',length(Mag.Angle)))
    Radio.Angle = seq(from = min(CNP[CNP$Label.Group == 'Radio', 'angle']),
                      to = max(CNP[CNP$Label.Group == 'Radio', 'angle']),
                      length.out = 100)
    Radio.Arc = data.frame(x = cos(pi/2 - Radio.Angle), y = sin(pi/2 - Radio.Angle), Label = rep('Radio',length(Radio.Angle)))
    FPC.Angle = seq(from = min(CNP[CNP$Label.Group == 'FPC', 'angle']),
                    to = max(CNP[CNP$Label.Group == 'FPC', 'angle']),
                    length.out = 100)
    FPC.Arc = data.frame(x = cos(pi/2 - FPC.Angle), y = sin(pi/2 - FPC.Angle), Label = rep('FPC',length(FPC.Angle)))
    Arcs = rbind(Feb.Arc, May.Arc, Nov.Arc, DEM.Arc, Mag.Arc, Radio.Arc, FPC.Arc)
    Arcs = data.frame(Arcs, Frequency = rep(0,nrow(Arcs)))
    Label.Angle = c(mean(Feb.Angle), mean(May.Angle), mean(Nov.Angle), mean(DEM.Angle), mean(Mag.Angle), mean(Radio.Angle), mean(FPC.Angle))
    Label.Arc = data.frame(x = cos(pi/2 - Label.Angle), y = sin(pi/2 - Label.Angle), Label = c('Feb', 'May', 'Nov', 'DEM', 'Mag', 'Radio', 'FPC'), Frequency = rep(0,length(Label.Angle)))
    ## Ploting the Chord Diagram
    Chord.p = ggplot(aes(x = x , y = y, group = factor(Interaction.Number), alpha = 100*Frequency/500), data = I.Store)
    Chord.p = Chord.p + geom_path(size = Line.Width) + coord_equal()
    Chord.p = Chord.p + geom_point(aes(x = x*(1+(Poly.Num-1)/8), y = y*(1+(Poly.Num-1)/8), shape = factor(Poly.Num), alpha = 100*Frequency/500), size = Point.Size, data = Points.xy)
    Chord.p = Chord.p + scale_alpha(name = expression(paste('% of Selected \n Models')), breaks = c(20,40,60,80,100), limits = c(0,100)) + labs(x = '', y = '', title = graph.title)
    Chord.p = Chord.p + geom_path(aes(x = (15/8)*x, y = (15/8)*y, group = Label), alpha = 1, data = Arcs, size = Line.Width)
    Chord.p = Chord.p + annotate(geom = 'text', x = (16.5/8)*Label.Arc$x, y = (16.5/8)*Label.Arc$y, label = Label.Arc$Label, size = Covariate.Lab.Size) + scale_shape(name = expression(paste('Single Term \n Polynomial \n Order'))) #, size = Covariate.Lab.Size
    Chord.p = Chord.p + theme(panel.background = element_rect(fill = 'white'), panel.grid = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), text = element_text(size = Title.Size), legend.position = Legend.Position, legend.key.size = Legend.Height.Width)
    Chord.p = Chord.p + annotate(geom = 'text', x = (13/8)*CNP$x, y = (13/8)*CNP$y, label = CNP$Acron, size = Covariate.Lab.Size, angle = c(rep(x = 90, times = ceiling(63/2)), rep(x = -90, times = floor(63/2))) - CNP$angle*360/(2*pi))
    return(Chord.p)}
