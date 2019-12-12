function [rad_pol,erro,number,aprox,struc] = linear(vector) %The introduced vector has 2 rows.

  aprox = [];

  for i=20:10:length(vector)
    [rad_pol,struc] = polyfit(vector(1:i,1),vector(1:i,2).^2,1); %polinomial fit

	R = struc.R; %The "R", whatever this it...
	d = ( R' * R ) \ eye(2); %The covariance matrix
	d = diag(d)'; %ROW vector of the diagonal elements.
	MSE = (struc.normr^2)/struc.df; %Variance of the residuals.
	ster = sqrt(MSE*d); %Standard errors
	t = rad_pol./ster; %Observed t-values. The bigger, the better.
	erro = sum(abs(t));
    %erro = (sqrt (diag (struc.C)/struc.df)*struc.normr)';
    aprox = [aprox;i,erro];
  endfor

  [mini,posi] = max(aprox(:,2)); %Minimum value of error (bigger t)
  number = aprox(posi,1); %Number of points necessary for the minimum value.
  [rad_pol,struc] = polyfit(vector(1:number,1),vector(1:number,2).^2,1);
  erro = (sqrt (diag (struc.C)/struc.df)*struc.normr)';%deviations of the parameters.

endfunction 

#That's...that's all, folks! 
