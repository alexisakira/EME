%% main file for solving stochastic growth model

clear
close all
clc;

%% figure formatting

set(0,'DefaultTextInterpreter','latex')
set(0,'DefaultAxesTickLabelInterpreter','latex');
set(0,'DefaultLegendInterpreter', 'latex')

set(0,'DefaultTextFontSize', 12)
set(0,'DefaultAxesFontSize', 12)
set(0,'DefaultLineLineWidth',1)

temp = get(gca,'ColorOrder');
c1 = temp(1,:);
c2 = temp(2,:);
c3 = temp(3,:);

set(groot,'defaultAxesToolbarVisible','off')

close all

%% define stochastic growth structure

% preference and technology parameters
beta = 0.95; % discount factor
gamma_h = 1.5; % relative risk aversion (high)
gamma_l = 0.5; % relative risk aversion (low)
A = [1.1 0.9]'; % productivity
nz = length(A); % number of states
alpha = 0.36; % capital share
delta = 0.08; % capital depreciation
epsilon = 0.1; % exogenous consumption
p = 0.2; % transition probability
P = [1-p p; p 1-p]; % transition probability matrix

% define parameters for algorithm
N = 100; % number of grid points
N_plot = 1000; % number of grid points for plotting
MaxIter = 400; % maximum number of iterations
tol = 1e-4; % error tolerance
m = 1; % parameter to skip maximization step

%% model 1

sg1.beta = beta;
gamma = [gamma_h gamma_h]';
sg1.u = @(c,z)((c + epsilon)^(1-gamma(z))/(1-gamma(z))); % utility function
sg1.f = @(k,z)(A(z)*k^alpha + (1-delta)*k); % production function
sg1.P = P; % transition probability matrix

sg1.MaxIter = MaxIter;
sg1.tol = tol;
sg1.m = m;
sg1.Vmat0 = zeros(nz,N); % initialize value function

% grid for computation
kstar = (max(A)/delta)^(1/(1-alpha));
aGrid = expGrid(0,round(2*kstar,-1),kstar/2,N); % construct exponential grid
sg1.aGrid = aGrid;

% grid for plotting
aMin = min(aGrid);
aMax = max(aGrid);
aGrid_plot = linspace(aMin,aMax,N_plot);

tic
sg1 = solve_sg(sg1);
toc

% plot results

% value function
figure
plot(sg1.aGrid,sg1.Vmat(1,:),'-','Color',c1); hold on
plot(sg1.aGrid,sg1.Vmat(2,:),'--','Color',c2);
xlabel('Resource')
ylabel('Lifetime utility')
legend('State 1','State 2','Location','NW')

fig = gcf;
exportgraphics(fig,'fig_sg1_v.pdf')

% consumption function
figure
plot(sg1.aGrid,sg1.Cmat(1,:),'-','Color',c1); hold on
plot(sg1.aGrid,sg1.Cmat(2,:),'--','Color',c2);
xlabel('Resource')
ylabel('Consumption')
legend('State 1','State 2','Location','NW')

fig = gcf;
exportgraphics(fig,'fig_sg1_c.pdf')

% value function along iteration
imax = sg1.imax;

figure
hold on
for i = 1:imax+1
    t = (i-1)/imax;
    v = interp1(sg1.aGrid,sg1.VMat(1,:,i),aGrid_plot,'spline');
    plot(aGrid_plot,v,'-','Color',(1-t)*[0 1 0] + t*[0 0 1]);
end
xlabel('Resource')
ylabel('Lifetime utility')

fig = gcf;
exportgraphics(fig,'fig_sg1_iter.pdf',...
    'Resolution',300,'ContentType','vector','Colorspace','gray')

% plot the error along iterations
maxError = zeros(1,imax+1);
for i = 1:imax
    maxError(i) = max(max(abs(sg1.VMat(:,:,i) - sg1.Vmat)));
end

iter = 0:imax;
figure
semilogy(iter,maxError,'Color',c1); hold on
semilogy(iter,beta.^iter,'k:');
xlabel('Number of iterations')
legend('Maximum error','$\beta^n$')

%% model 2
sg2 = sg1;
gamma = [gamma_l gamma_l]';
sg2.u = @(c,z)((c)^(1-gamma(z))/(1-gamma(z))); % utility function

tic
sg2 = solve_sg(sg2);
toc

% plot results

% value function
figure
plot(sg2.aGrid,sg2.Vmat(1,:),'-','Color',c1); hold on
plot(sg2.aGrid,sg2.Vmat(2,:),'--','Color',c2);
xlabel('Resource')
ylabel('Lifetime utility')
legend('State 1','State 2','Location','NW')

fig = gcf;
exportgraphics(fig,'fig_sg2_v.pdf')

% consumption function
figure
plot(sg2.aGrid,sg2.Cmat(1,:),'-','Color',c1); hold on
plot(sg2.aGrid,sg2.Cmat(2,:),'--','Color',c2);
xlabel('Resource')
ylabel('Consumption')
legend('State 1','State 2','Location','NW')

fig = gcf;
exportgraphics(fig,'fig_sg2_c.pdf')

% value function along iteration
imax = sg2.imax;

figure
hold on
for i = 1:imax+1
    t = (i-1)/imax;
    v = interp1(sg2.aGrid,sg2.VMat(1,:,i),aGrid_plot,'spline');
    plot(aGrid_plot,v,'-','Color',(1-t)*[0 1 0] + t*[0 0 1]);
end
xlabel('Resource')
ylabel('Lifetime utility')

fig = gcf;
exportgraphics(fig,'fig_sg2_iter.pdf',...
    'Resolution',300,'ContentType','vector','Colorspace','gray')

%% model 3
sg3 = sg2;
V0 = 20*(1-sin(8*pi*aGrid/aMax)); % initial guess is sine function
sg3.Vmat0 = repmat(V0,nz,1);

tic
sg3 = solve_sg(sg3);
toc

% plot results

% value function
figure
plot(sg3.aGrid,sg3.Vmat(1,:),'-','Color',c1); hold on
plot(sg3.aGrid,sg3.Vmat(2,:),'--','Color',c2);
xlabel('Resource')
ylabel('Lifetime utility')
legend('State 1','State 2','Location','NW')

fig = gcf;
exportgraphics(fig,'fig_sg3_v.pdf')

% consumption function
figure
plot(sg3.aGrid,sg3.Cmat(1,:),'-','Color',c1); hold on
plot(sg3.aGrid,sg3.Cmat(2,:),'--','Color',c2);
xlabel('Resource')
ylabel('Consumption')
legend('State 1','State 2','Location','NW')

fig = gcf;
exportgraphics(fig,'fig_sg3_c.pdf')

% value function along iteration
imax = sg3.imax;

figure
hold on
for i = 1:imax+1
    t = (i-1)/imax;
    v = interp1(sg3.aGrid,sg3.VMat(1,:,i),aGrid_plot,'spline');
    plot(aGrid_plot,v,'-','Color',(1-t)*[0 1 0] + t*[0 0 1]);
end
xlabel('Resource')
ylabel('Lifetime utility')

fig = gcf;
exportgraphics(fig,'fig_sg3_iter.pdf',...
    'Resolution',300,'ContentType','vector','Colorspace','gray')

% create cover image

figure
hold on
for i = 1:imax+1
    t = (i-1)/imax;
    v = interp1(sg3.aGrid,sg3.VMat(1,:,i),aGrid_plot,'spline');
    plot(aGrid_plot,v,'-','Color',(1-t)*[0 1 0] + t*[0 0 1]);
end
axis off

% color
fig = gcf;
exportgraphics(fig,'fig_cover_color.pdf',...
    'Resolution',300,'ContentType','vector')

% gray
fig = gcf;
exportgraphics(fig,'fig_cover_gray.pdf',...
    'Resolution',300,'ContentType','vector','Colorspace','gray')

%% model 4 (optimistic policy iteration)
sg4 = sg2;
sg4.m = 10;

tic
sg4 = solve_sg(sg4);
toc

% plot results

% value function
figure
plot(sg4.aGrid,sg4.Vmat(1,:),'-','Color',c1); hold on
plot(sg4.aGrid,sg4.Vmat(2,:),'--','Color',c2);
xlabel('Resource')
ylabel('Lifetime utility')
legend('State 1','State 2','Location','NW')

fig = gcf;
exportgraphics(fig,'fig_sg4_v.pdf')

% consumption function
figure
plot(sg4.aGrid,sg4.Cmat(1,:),'-','Color',c1); hold on
plot(sg4.aGrid,sg4.Cmat(2,:),'--','Color',c2);
xlabel('Resource')
ylabel('Consumption')
legend('State 1','State 2','Location','NW')

fig = gcf;
exportgraphics(fig,'fig_sg4_c.pdf')

% value function along iteration
imax = sg4.imax;

figure
hold on
for i = 1:imax+1
    t = (i-1)/imax;
    v = interp1(sg4.aGrid,sg4.VMat(1,:,i),aGrid_plot,'spline');
    plot(aGrid_plot,v,'-','Color',(1-t)*[0 1 0] + t*[0 0 1]);
end
xlabel('Resource')
ylabel('Lifetime utility')

fig = gcf;
exportgraphics(fig,'fig_sg4_iter.pdf',...
    'Resolution',300,'ContentType','vector','Colorspace','gray')

%close all
%save sg.mat