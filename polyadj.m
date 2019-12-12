## Copyright (C) 2008-2012 Ben Abbott
##
## This file is part of Octave.
##
## Octave is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or (at
## your option) any later version.
##
## Octave is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with Octave; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {Function File} {} polyadj ()
##   Adjust to a polynomium a vector, with a minimum correlation factor of 0.1.
##
## @code{padj = polyadj (@var{x},@var{fx},@var{rmin})} returns polynomium that adjusts to the vector
## @var{fx} with @var{x} values with a correlation factor larger than @var{rmin}.
## @code{[padj, rcor, padjstruc] = polyadj (@var{x},@var{fx},@var{rmin})} returns a polynomium that adjusts to the vector
## @var{fx} with @var{x} values with a correlation factor larger than @var{rmin} in @var{padj}, the correlation factor in @var{rcor} and the
## polynomium structure in @var{padjstruc}.
## @end deftypefn

## Author: Gonzalo Rodríguez Prieto (gonzalo.rprietoATuclm.es)
## Created: Nov 2013
## Modified: Jun 2015



function [padj, padjstruc, t] = polyadj(x,fx)


#########
# Control error part
#########
if (nargin!=2) #If not enough parameters are given.
 error("poladj: 2 parameters are needed.");
endif;

if (isscalar(x) ==1) #Works when a vector is passed, not a escalar or other type.
 error("poladj: Variable x must be a vector.");
endif;

if (isscalar(fx) ==1) #Works when a vector is passed, not a escalar or other type.
 error("poladj: Variable fx must be a vector.");
endif;

if ( length(x) != length(fx) ) #Length of both vectors are not the same.
 error("poladj: Lengths of vectors x and fx must coincide.");
endif;



#########
# Adjusting loop part
#########
n = 1;
t_vector = [];
t_tot = 666;

#Fitting loop :
while ( (n<10) && (t_tot>0) )
	#Fitting:
	[padj, padjstruc] = polyfit(x,fx, n-1);

	#Calculation of the statistical errors of the fit:
		%From http://www.facstaff.bucknell.edu/maneval/help211/fitting.html
	R = padjstruc.R; %The "R", whatever this it...
	d = ( R' * R)\eye(n); %The covariance matrix
	d = diag(d)'; %ROW vector of the diagonal elements.
	MSE = (padjstruc.normr^2)/ padjstruc.df; %Variance of the residuals.
	st_er = sqrt(MSE*d); %Standard errors
	t = padj./st_er; %Observed t-values. The bigger, the better.
	t_tot = sum(t);
	t_vector = [t_vector,t_tot];
	n = n +1;
endwhile;

[tmax,pos] = max(t_vector); %Better adjustment to the data when t is maximum.

[padj, padjstruc] = polyfit(x,fx, pos-1);

#Calculation of the statistical errors of the fit:
	%From http://www.facstaff.bucknell.edu/maneval/help211/fitting.html
R = padjstruc.R; %The "R", whatever this it...
d = ( R' * R)\eye(pos); %The covariance matrix
d = diag(d)'; %ROW vector of the diagonal elements.
MSE = (padjstruc.normr^2)/ padjstruc.df; %Variance of the residuals.
st_er = sqrt(MSE*d); %Standard errors
t = padj./st_er; %Observed t-values. The bigger, the better.

endfunction;

#That's...that's all, folks! 
