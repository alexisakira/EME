%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get_cv
% (c) 2022 Alexis Akira Toda
% 
% Purpose: 
%       Solve stochastic growth model by value function iteration
%
% Usage:
%       sg = solve_sg(sg)
%
% Inputs:
% sg:       stochastic growth structure
%
% Output:
% sg:       stochastic growth structure
%
% Version 1.0: June 22, 2022
% Version 1.1: September 21, 2023
%           Added option for optimistic policy iteration (OPI)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function sg = solve_sg(sg)

nz = size(sg.P,1); % number of states
N = length(sg.aGrid); % number of grid points
MaxIter = sg.MaxIter; % maximum number of iterations
tol = sg.tol; % error tolerance
m = sg.m; % parameter to skip maximization step

VMat = zeros(nz,N,MaxIter); % store value function for plotting later
Vmat0 = sg.Vmat0; % initialize value function
VMat(:,:,1) = Vmat0; % initialization

% optimistic policy iteration
for i = 1:MaxIter
    if rem(i-1,m) == 0 % carry out maximization
        [Cmat,Vmat] = get_cv(sg,Vmat0);
    else % no maximization
        for z = 1:nz
            for n = 1:N
                a = sg.aGrid(n); % current resource
                c = Cmat0(z,n); % current consumption
                Vmat(z,n) = get_cv_obj(c,a,z,sg,Vmat0); % continuation value
            end
        end
    end
    VMat(:,:,i+1) = Vmat; % store value function
    if max(max(abs(Vmat - Vmat0))) < tol % check convergence
        str = 'Converged after %3.0f iterations\n';
        fprintf(str,i) % print message
        sg.Cmat = Cmat; % consumption function
        sg.Vmat = Vmat; % value function
        sg.VMat = VMat(:,:,1:i+1); % remove irrelevant part of VMat
        break
    else % if not converged, update value function
        Cmat0 = Cmat;
        Vmat0 = Vmat;
    end
end

sg.imax = i; % number of iterations required for convergence

end

