%% simple_energy_system_optimization.m
% Simple techno-economic capacity expansion model
% for a small power system with solar, wind and gas.


clear; clc; close all;

%%  Time resolution and synthetic data
days          = 7;            % model one week
hoursPerDay   = 24;
H             = days * hoursPerDay;
t             = (0:H-1)';     % hour index

% Synthetic load profile [MW]
baseLoad      = 60;
dailyAmp      = 25;
weeklyAmp     = 10;

hourOfDay     = mod(t, 24);
dailyShape    = 0.5 + 0.5*sin(2*pi*(hourOfDay-8)/24);  % peak in the evening
weeklyShape   = 0.8 + 0.2*sin(2*pi*t/(24*7));          % mild weekly variation

loadMW        = baseLoad + dailyAmp.*dailyShape.*weeklyShape ...
                + 4*randn(H,1);                        % small noise
loadMW        = max(loadMW, 20);                       % no negative load

% Normalised resource profiles (0–1)
solarNorm     = max(0, sin((hourOfDay-6)/12*pi));      % day-only solar
windNorm      = 0.6 + 0.3*sin(2*pi*t/(24*3)) ...       % slower wind variation
                + 0.1*randn(H,1);
windNorm      = max(windNorm, 0.05);                   % minimum availability
gasAvail      = ones(H,1);                             % gas always available

% Collect availability into a matrix: rows = hours, cols = technologies
avail         = [solarNorm, windNorm, gasAvail];

%%  Technology data (very simple, illustrative)
techNames     = {'Solar','Wind','Gas'};
nTech         = numel(techNames);

% CAPEX [EUR/kW]
capex_kW      = [600; 1200; 800];       % demo data
lifetime_yr   = [25; 25; 30];
wacc          = 0.07;                   % discount rate

% Annuity factors and annualised CAPEX [EUR/kW-year]
annuity       = zeros(nTech,1);
for i = 1:nTech
    r         = wacc;
    n         = lifetime_yr(i);
    annuity(i)= r*(1+r)^n / ((1+r)^n - 1);
end
capCostAnn_kW = capex_kW .* annuity;    % EUR/kW-year

% Variable generation cost [EUR/MWh]
varCost       = [0; 0; 80];             % fuel-free vs gas

%% Decision variables
% x = [capacities(MW); generation(MW) for each tech and hour]
nCap          = nTech;
nGen          = H * nTech;
nVar          = nCap + nGen;

% Helper: index of generation variable for tech j at hour h
% idx = nCap + (j-1)*H + h

%% Objective function: minimise total annual cost
% Convert CAPEX to EUR/MW-year: multiply by 1000
f             = zeros(nVar,1);
f(1:nCap)     = capCostAnn_kW * 1000;   % capacity part

for j = 1:nTech
    genIdx    = nCap + (j-1)*H + (1:H);
    f(genIdx) = varCost(j);             % variable cost part
end

%%  Constraints
% Generation limited by capacity and resource:
%     gen(h,j) <= cap(j) * avail(h,j)

nConstrCap    = H * nTech;
%  Demand balance:
%     sum_j gen(h,j) >= load(h)  (we implement as -sum_j gen <= -load)

nConstrLoad   = H;
A             = zeros(nConstrCap + nConstrLoad, nVar);
b             = zeros(nConstrCap + nConstrLoad, 1);

row = 0;

% Capacity + availability constraints
for j = 1:nTech
    for h = 1:H
        row       = row + 1;
        genCol    = nCap + (j-1)*H + h;
        A(row,genCol) = 1;                % gen(h,j)
        A(row,j)      = -avail(h,j);      % -cap(j)*avail(h,j)
        b(row)        = 0;
    end
end

% Demand constraints: -sum_j gen(h,j) <= -load(h)
for h = 1:H
    row = row + 1;
    for j = 1:nTech
        genCol          = nCap + (j-1)*H + h;
        A(row,genCol)   = -1;
    end
    b(row) = -loadMW(h);
end

% No equality constraints
Aeq          = [];
beq          = [];

% Lower bounds: no negative capacities or generation
lb           = zeros(nVar,1);
ub           = [];       % no explicit upper bounds

%% Solve the linear program
options      = optimoptions('linprog','Display','none');
[x,fval,exitflag,output] = linprog(f, A, b, Aeq, beq, lb, ub, options);

if exitflag ~= 1
    warning('Solver did not find a proven optimal solution. Exitflag = %d', exitflag);
end

%% Extract results
capMW        = x(1:nCap);

genSolar     = x(nCap + (0*H+1):(0*H+H));
genWind      = x(nCap + (1*H+1):(1*H+H));
genGas       = x(nCap + (2*H+1):(2*H+H));

totalEnergy  = sum(genSolar + genWind + genGas);   % MWh over model period
totalCost    = fval;                               % EUR per year-equivalent
LCOE         = totalCost / totalEnergy;            % EUR/MWh (approx.)

%% Display key numbers
fprintf('=== Optimal capacities ===\n');
for j = 1:nTech
    fprintf('%s: %.1f MW\n', techNames{j}, capMW(j));
end

fprintf('\nTotal annualised system cost: %.1f EUR\n', totalCost);
fprintf('Approximate system LCOE:      %.2f EUR/MWh\n\n', LCOE);

%%  Plot a sample day
sampleDay    = 2;  % choose a day to plot
idx          = (sampleDay-1)*24 + (1:24);

figure;
area(1:24, [genSolar(idx), genWind(idx), genGas(idx)]);
hold on;
plot(1:24, loadMW(idx), 'k--', 'LineWidth', 1.5);
xlabel('Hour of day');
ylabel('Power [MW]');
legend('Solar','Wind','Gas','Load','Location','best');
title('Optimal dispatch for sample day');
grid on;