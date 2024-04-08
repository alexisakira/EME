%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get_cv
% (c) 2022 Alexis Akira Toda
% 
% Purpose: 
%       Update optimal consumption and value functions from Bellman equation
%
% Usage:
%       [Cmat,Vmat] = get_cv(sg,Vmat0)
%
% Inputs:
% sg:       stochastic growth structure
% Vmat0:    matrix of value functions
%
% Output:
% Cmat:     matrix of consumption functions
% Vmat:     matrix of value functions
%
% Version 1.0: June 22, 2022
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function [Cmat,Vmat] = get_cv(sg,Vmat0)
% 

aGrid = sg.aGrid;
N = length(aGrid); % number of grid points
nz = size(sg.P,1); % number of states

Cmat = zeros(nz,N);
Vmat = zeros(nz,N);

for z = 1:nz
    for n = 1:N
        a = aGrid(n);
        if a == 0
            c = 0;
            Cmat(z,n) = c;
            Vmat(z,n) = get_cv_obj(c,a,z,sg,Vmat0);
        else % a > 0
            func = @(c)(-get_cv_obj(c,a,z,sg,Vmat0));
            [c,fval] = fminbnd(func,0,a); % compute optimal consumption
            Cmat(z,n) = c;
            Vmat(z,n) = -fval;
        end
    end
end

end

