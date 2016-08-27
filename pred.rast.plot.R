# R Code to conduct an example analysis with functions provided in the
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

# Function to take a data frame containing the coordinates and values (columns) of the pixels (rows) in a raster and plot them with the ggplot2 geom_raster() then overlay some point observations on the plot and use the same colour scale for the raster fill as for the fill of the points
# This functions is just for convenience in quickly producing raster plots with the particular ggplot2 aesthetics and scales I have found useful for this purpose (you could definately just use geom_raster() and your own selection of scales and aesthetics if you'd rather).

pred.rast.plot <- function(
    rast.df = data.frame(B1.CRS.df[,c('Easting','Northing')], Fill = Areal.MAP.LAR.F95.TS35.n500$Results[,'MA.Pred']),
    fill.name = c('% SOC','Int.Width','Var'), # choose one of these
    points.df = Response.df, # three columns, Easting, Northing, then fill
    point.fill.column.name = c('Carbon','Residual'),  # choose one of these
    colour.scale.type = c('rainbow', 'grey', 'heat', 'heat.custom', 'cyan.black.green'), # choose one of these
    auto.col.lims = TRUE,
    colour.scale.limits = c(min(c(min(Areal.MAP.LAR.F95.TS35.n500$Results[,'MA.Pred']), min(Response.df$Carbon))), max(c(max(Areal.MAP.LAR.F95.TS35.n500$Results[,'MA.Pred']), max(Response.df$Carbon)))),
    out.of.range.colour = 'black',
    overlay.points = TRUE,
    point.boundary.colour = c('black','white','cyan'),  # choose one of these
    point.size = 3,
    plot.title = '63C P4I2 LAR F95 TS35 CV500 MAP SOC (raster) & Observed SOC (points)',
    title.text.size = 14,
    axis.text.size = 14,
    axis.tick.step = 100,
    x.axis.first.tick = c(368100, round(x = quantile(x = rast.df$Easting, prob = 0.1), digits = -3)),
    x.axis.last.tick = c(369800,  round(x = quantile(x = rast.df$Easting, prob = 0.9), digits = -3)),
    y.axis.first.tick = c(663200, round(x = quantile(x = rast.df$Northing, prob = 0.1), digits = -3)),
    y.axis.last.tick = c(6634200, round(x = qunatile(x = rast.df$Northing, prob = 0.9), digits = -3))){    
    require('ggplot2')
    require('grid')
    if(!(ncol(points.df) == 3)){return(c('Please supply a point.df with Eastings in the first column, Northings in the second column and the value of the point (e.g. %SOC or Residual) in the third column'))}
    if(!(colnames(points.df)[1] == 'Easting')){return(c('Please supply a point.df with Eastings in the first column, Northings in the second column and the value of the point (e.g. %SOC or Residual) in the third column'))}
    if(!(colnames(points.df)[2] == 'Northing')){return(c('Please supply a point.df with Eastings in the first column, Northings in the second column and the value of the point (e.g. %SOC or Residual) in the third column'))}
    if(!(ncol(rast.df) == 3)){return(c('Please supply a rast.df with Eastings in the first column, Northings in the second column and the value of the point (e.g. %SOC or Residual) in the third column'))}
    if(!(colnames(rast.df)[1] == 'Easting')){return(c('Please supply a rast.df with Eastings in the first column, Northings in the second column and the value of the point (e.g. %SOC or Residual) in the third column'))}
    if(!(colnames(rast.df)[2] == 'Northing')){return(c('Please supply a rast.df with Eastings in the first column, Northings in the second column and the value of the point (e.g. %SOC or Residual) in the third column'))}
    colnames(points.df) = c('Easting', 'Northing', 'Fill')
    colnames(rast.df) = c('Easting', 'Northing', 'Fill')
    Plot = ggplot(aes(x = Easting, y = Northing, fill = Fill), data = rast.df) +
            geom_raster() +
            coord_equal() +
            theme(text = element_text(size = title.text.size, colour = 'black'),
                  axis.text.x = element_text(size = axis.text.size, angle = 90, colour = 'black', vjust = 0.5),
                  axis.text.y = element_text(size = axis.text.size, angle = 0, colour = 'black'),
                  legend.position = 'right',
                  legend.key.height = unit(x = 1.4, units = 'inches'),
                  panel.background = element_rect(fill = 'white'),
                  panel.grid.major = element_line(colour = 'black'),
                  panel.grid.minor = element_line(colour = 'black'),
                  axis.ticks = element_line(colour = 'black')) +
            labs(title = plot.title, fill = fill.name) +
            scale_x_continuous(breaks = seq(from = x.axis.first.tick, to = x.axis.last.tick, by = axis.tick.step)) +
            scale_y_continuous(breaks = seq(from = y.axis.first.tick, to = y.axis.last.tick, by = axis.tick.step))
    if(overlay.points){
        Plot = Plot + geom_point(aes(x = Easting, y = Northing, fill = Fill), shape = 21, colour = point.boundary.colour, size = point.size, data = points.df)}
    if(colour.scale.type == 'rainbow'){
        rbc.v = rainbow(n = 1e4, start = 0, end = 0.7)
        rbc.v = rev(rbc.v)
        Plot = Plot + scale_fill_gradientn(colours =  rbc.v, limits = colour.scale.limits, na.value = out.of.range.colour)}
    if(colour.scale.type == 'grey'){
        Plot = Plot + scale_fill_gradient(low = 'black', high = 'white', limits = colour.scale.limits, na.value = out.of.range.colour)}
    if(colour.scale.type == 'heat'){
        heat.v = heat.colors(n = 1e4)
        Plot = Plot + scale_fill_gradientn(colours =  heat.v, limits = colour.scale.limits, na.value = out.of.range.colour)}
    if(colour.scale.type == 'heat.custom'){
        heat.v = heat.colors(n = 1e4)
        Plot = Plot + scale_fill_gradientn(colours =  c('black','darkred','red','orange','yellow','white'), limits = colour.scale.limits, na.value = out.of.range.colour)}
    if(colour.scale.type == 'cyan.black.green'){
        heat.v = heat.colors(n = 1e4)
        Plot = Plot + scale_fill_gradientn(colours =  c('cyan','black','green'), limits = colour.scale.limits, na.value = out.of.range.colour)}
    return(Plot)}
