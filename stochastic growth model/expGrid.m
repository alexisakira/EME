%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% expGrid
% (c) 2019 Alexis Akira Toda
% 
% Purpose: 
%       Construct an exponential grid
%
% Usage:
%       grid = expGrid(a,b,c,N)
%
% Inputs:
% a:    lower endpoint of grid
% b:    upper endpoint of grid
% c:    median grid point
% N:    number of grid points
%
% Output:
% grid: grid
%
% Version 1.0: June 16, 2019
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function grid = expGrid(a,b,c,N)

% some error checking
if (a >= b)||(c <= a)||(c >= (a+b)/2)
    error('it must be a < c < (a+b)/2')
end

s = (c^2-a*b)/(a+b-2*c); % shift parameter
temp = linspace(log(a+s),log(b+s),N+1); % even grid in log scale
grid = exp(temp(2:end))-s;

end

