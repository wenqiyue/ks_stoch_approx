function value=adjust_angles(phi)
%This function tries to adjust the input value (as angles) such that they
%would be within the range between -pi and pi by subtracting off (or
%adding) the required numbers of 2pis. It would be useful for visualizing
%results, etc.

n = ceil((phi-pi)/(2*pi));

phi2=phi-n*2*pi;

value=phi2;

end