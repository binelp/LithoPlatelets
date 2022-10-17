%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ETH Zurich, Switzerland
% Separation Processes Laboratory
%
% Project:  Lithoplatelets
% Year:     2022
% MATLAB:   R2019b, Windows 64bit
% Authors:  Anna Jaeggi (AJ)
%
% Purpose:
% Generate the distribution of the physical angles theta and phi from
% unifrom orientations, in the binned format
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [f_bin] = physAngleDistFromUniform(nbins, normalization)
    delta = (pi/2)/nbins;
    for ii=1:nbins
        for jj=1:nbins
            a = (ii-1)*delta;
            b = (ii)*delta;
            c = (jj-1)*delta;
            d = (jj)*delta;

            theta_bin(ii,jj) = (a+b)/2;
            phi_bin(ii,jj) = (c+d)/2;
            intPhi = @(theta) -asin(cos(d)./sin(theta)).*sin(theta) + asin(cos(c)./sin(theta)).*sin(theta);       
            intTheta = integral(intPhi,a,b);
            intTheta = real(intTheta);
            f_bin(ii,jj) = intTheta;
        end
    end
    if normalization == "integral"
        f_bin = f_bin/(nansum(f_bin,'all')*delta^2);
    end
    if normalization == "sum"
        f_bin = f_bin/(nansum(f_bin,'all'));
    end
end