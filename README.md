<h1> larc - Least Angle Regression Companion </h1>

This repository contains the data and code necessary to replicate the analysis described in the PLOS ONE article:

['Ultrahigh Dimensional Variable Selection for Interpolation of Point Referenced Spatial Data: A Digital Soil Mapping Case Study'](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0162489)

by [Benjamin R. Fitzpatrick (BRF)](http://orcid.org/0000-0003-1916-0939), David W. Lamb (DWL) and Kerrie Mengersen (KM).

Code and repository authorship was the sole responsibility of BRF.

The code file ```example_analysis.R``` illustrates how the functions included in this repository may be used to replicate the analysis described in the article.
The article discusses the relevant theory and demonstrates the application of these methods to a geostatistical case study.
This repository contains a set of functions written in the [R Language for Statistical Computing](https://cran.r-project.org/).
The analysis this repository enables makes heavy use of the Least Angle Regression (LAR) algorithm for finding Least Absolute Shrinkage Selection Operator (LASSO) regularised solutions to multiple linear regression problems.
An R package for conducting Least Absolute Shrinkage Selection Operator (LASSO) variable selection with the LAR algorithm already exists and is hosted on the Comprehensive R Archive Network under the name '[lars](https://cran.r-project.org/web/packages/lars/index.html)'.
This repository makes heavy use of functions from the 'lars' package.

This repository contains functions that: 
* randomly generate unique divisions of a sequence of numbers into two groups of user specified sizes (the intent being that these two groups of numbers are used as row indices to create training and validation sets from a full dataframe)
* use the LAR algorithm within a cross validation scheme in a manner that permits greater control of the particulars than is provided by the ```cv.lars( )``` function from the 'lars' package
* use chord diagrams to visualise the covariate selection frequencies that result from conducting LAR within a cross validation scheme
* model average the predictions from the models selected for each of the training sets in the cross validation scheme
* interpolate a geostatistical response variable to a full cover predicted raster via such model averaged predictions. 

The functions provided here depend on the R packages:
 * [lars](https://cran.r-project.org/web/packages/lars/index.html)
 * [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html)
 * [randtoolbox](https://cran.r-project.org/web/packages/randtoolbox/index.html)
 * [raster](https://cran.r-project.org/web/packages/raster/index.html)
 * [leaps](https://cran.r-project.org/web/packages/leaps/index.html)

Copyright (C) 2016 Benjamin R. Fitzpatrick, David W. Lamb and Kerrie Mengersen.

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program, in the form a text file titled 'LICENSE.
If not, see http://www.gnu.org/licenses/.
