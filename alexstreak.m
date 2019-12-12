
more off; %To make the lines display when they appear in the script, not at the end of it.

clear;

tic; %Time control

sweep = 50; %µs of sweeping time in the streak images.

dist = 96.15e-6; %m per pixel in distance.

%shot = "ALEX121"; %Shot to read.
 fflush (stdout);
shot = input("Shot?","s");

####
## Opening the picture and choosing some parameters.
####
filename = horzcat(shot,".dat"); %Compatibility with previous code and data treatment.

disp("Opening streak image.");

[file, msg] = fopen(filename, "r");
if (file == -1) 
   error ("Unable to open file name: %s, %s",filename, msg); 
endif; 

#Charge the image as a matrix: (IT MUST BE ALREADY FORMATTED AS SUCH)
Streak = dlmread(file);
fclose(file);

Streak = Streak(:,1:(end-1)); %To remove the extra added final column...(Bug 38711: http://savannah.gnu.org/bugs/?38711)

ti = sweep/rows(Streak); %Time constant in ¡¡¡¡¡µs!!!!

disp("Processing it.");

#Remove noise:
noise = mean(mean(Streak(1:100,:),2)); %First on the rows, then the row of means.

Streak = Streak - noise; %To remove the noise on the image.

light = 0.05 * max(max(Streak)); %Light level to find the borders ( (first number)% of total light maximum)

radius = [];
posi = [];

for i=1:rows(Streak)
vector = Streak(i,:); %Row i, all columns.
	rad = find(vector>light); %Positions of the light borders

	if length(rad)>1 %To avoid error of just one point...
		radius = [radius;i,rad(1),rad(end)]; %Storing them as time - radius values in pixels.
	endif;
endfor;

%Initiating the radius in zero:
radius(:,2) = radius(:,2)-radius(1,2);
radius(:,3) = radius(:,3)-radius(1,3);

radius = [radius(:,1).*ti,(-radius(:,2)+radius(:,3)).*0.5.*dist.*1e3]; %Conversion to µs and milimeters.



###
# Data processing:
###

#Adjusting the self-similar part:
	%I obtain the polynomium, its components error, the number of points necessary for it, and the vector with points vs errors.
[rad_self,rad_self_err,rad_self_number,rad_aprox_vec,rad_self_struc] = linear([radius(:,1),radius(:,2)]);

#Adjusting ot a polynomium the other part:
	%Using polyadj function to adjust the radius from end of self-similarity to maximum expansion to a polynomium:
[rad_other,rad_other_struc] = polyadj(radius(rad_self_number+1:length(radius),1),radius(rad_self_number+1:length(radius),2));
%%rad_other_err = (sqrt (diag (rad_other_struc.C)/rad_other_struc.df)*rad_other_struc.normr)';

%Making the 2 polynomia:
r_self = polyval(rad_self,radius(1:rad_self_number,1)); %Self similar part.
r_other = polyval(rad_other,radius(rad_self_number+1:length(radius),1)); %Other part.

#Calculation of the statistical errors of the fits:
	%From http://www.facstaff.bucknell.edu/maneval/help211/fitting.html
R_self = rad_self_struc.R; %The "R", whatever this it...
d_self = ( R_self' * R_self)\eye(2); %The covariance matrix
d_self = diag(d_self)'; %ROW vector of the diagonal elements.
MSE_self = (rad_self_struc.normr^2)/rad_self_struc.df; %Variance of the residuals.
ster_self = sqrt(MSE_self*d_self); %Standard errors
t_self = rad_self./ster_self; %Observed t-values. The bigger, the better.
width_self = 3.5*ster_self; %Confidence values (99% confidence for 7 points or more...)
	%From http://www.sjsu.edu/faculty/gerstman/StatPrimer/t-table.pdf

R_other = rad_other_struc.R; %The "R", whatever this it...
d_other = ( R_other' * R_other)\eye(length(rad_other)); %The covariance matrix
d_other = diag(d_other)'; %ROW vector of the diagonal elements.
MSE_other = (rad_other_struc.normr^2)/rad_other_struc.df; %Variance of the residuals.
ster_other = sqrt(MSE_other*d_other); %Standard errors
t_other = rad_other./ster_other; %Observed t-values The bigger, the better.
width_other = 3.5*ster_other; %Confidence values (99% confidence for 7 points or more...)
	%From http://www.sjsu.edu/faculty/gerstman/StatPrimer/t-table.pdf


###
# Self-similar energy extraction:
###
 En_self = sqrt(rad_self(1))^4*1.22*0.02; %Joules

disp("Self-similar energy:");
disp(En_self);

###
# Data imaging:
###

%Making the image:
plot(radius(:,1),radius(:,2),"*k"); %Black data
hold; plot(radius(1:rad_self_number,1),r_self.^(0.5),"+-r");%Red self similar part.
plot(radius(rad_self_number+1:length(radius),1),r_other,"*-g");%Green other part.
title_shot = horzcat(shot," Self Energy: ",num2str(En_self));
title(title_shot);
ylabel("radius (mm)");
xlabel("time (microseconds)");
legend({"exp. data","self-sim.","other fit."});
print(horzcat(shot,"-radius.pdf"),"-dpdf")

###
# Saving files:
###
%Radius points.
name = horzcat(shot,"_radio.txt"); %Adding the right sufix to the shot name.
output = fopen(name,"w"); %Opening the file.
fdisp(output,"descriptor t(µs)  r(m)"); %Veusz format
redond = [4 3]; %Saved precision 
rad = [radius(:,1) radius(:,2)]; %Making the vector a column
display_rounded_matrix(rad, redond, output); 
fclose(output); %Closing the file.

%Self-similar points.
name = horzcat(shot,"_self_points.txt"); %Adding the right sufix to the shot name.
output = fopen(name,"w"); %Opening the file.
fdisp(output,"descriptor t_s(µs)  r_self(m)"); %Veusz format
redond = [4 3]; %Saved precision 
rad_po_self = [radius(1:rad_self_number,1) r_self.^0.5]; %Making the vector a column
display_rounded_matrix(rad_po_self, redond, output); 
fclose(output); %Closing the file.

%Other points.
name = horzcat(shot,"_other_points.txt"); %Adding the right sufix to the shot name.
output = fopen(name,"w"); %Opening the file.
fdisp(output,"descriptor t_o(µs)  r_other(m)"); %Veusz format
redond = [4 3]; %Saved precision 
rad_ot_self = [radius(rad_self_number+1:end,1) r_other]; %Making the vector a column
display_rounded_matrix(rad_ot_self, redond, output); 
fclose(output); %Closing the file.

%Fits data
name = horzcat(shot,"_fits.txt"); %Adding the right sufix to the shot name.
output = fopen(name,"w"); %Opening the file.
fdisp(output,shot); fdisp(output," "); %Shot name
fdisp(output,"Self similar solution (mm/µs^0.5) (mm)");
fdisp(output,sqrt(rad_self));
fdisp(output,"	error");
fdisp(output,sqrt(width_self));
fdisp(output,"T values of the approx.:");
fdisp(output,t_self);
fdisp(output," ");
fdisp(output,"Self similar energy(Joules)");
fdisp(output,En_self);
fdisp(output,"	en. error");
fdisp(output,sqrt(width_self(1)));
fdisp(output," "); fdisp(output," ");
fdisp(output,"Other solution (polynomia, first left value higher)");
fdisp(output,rad_other);
fdisp(output,"	error");
fdisp(output,width_other);
fdisp(output,"T values of the other solution:");
fdisp(output,t_other);
fclose(output); %Closing the file.


###
# Total processing time
###
timing = toc;
disp("Script alexstreak execution time:")
disp(timing)
disp(" seconds")

more on; %To make the lines display when they appear in the script, not at the end of it.

#That's...that's all, folks! 

