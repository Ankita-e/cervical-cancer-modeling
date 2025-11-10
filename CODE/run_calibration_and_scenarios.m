%% run_calibration_and_scenarios.m
% Complete cervical cancer model calibration and scenario forecasting
% Date: [Today]
% Requires: cervical_cancer_tidy_numeric.xlsx in working folder

clear; close all; clc;

%% --- Load data (auto detect file)
[fn, pn] = uigetfile('*.xlsx', 'Select cervical_cancer_tidy_numeric.xlsx file');
if isequal(fn,0)
    error('No file selected. Place the dataset in working folder or pick it via dialog.');
end
fullpath = fullfile(pn, fn);
fprintf('Loading data from %s\n', fullpath);
T = readtable(fullpath);

% Expect columns: AgeGroup, YearNum, Incidence, Mortality
if any(strcmpi(T.Properties.VariableNames,'YearNum'))
    years = T.YearNum;
else
    yrs = string(T.Year);
    yrs = replace(yrs, ["–","—"], "-");
    yrs = regexprep(yrs,'\s+','');
    splitYr = split(yrs,'-');
    yearNum = zeros(size(yrs));
    for i=1:numel(yrs)
        p = splitYr(i,:);
        if size(p,2)>=2 && ~missing(p(2))
            a = str2double(p(1)); b = str2double(p(2));
            yearNum(i) = mean([a,b]);
        else
            yearNum(i) = str2double(p(1));
        end
    end
    years = yearNum;
end

% Aggregate across AgeGroup to national-level by YearNum
df = table(years, T.Incidence, T.Mortality, 'VariableNames', {'YearNum','Incidence','Mortality'});
agg = varfun(@mean, df, 'InputVariables', {'Incidence','Mortality'}, 'GroupingVariables','YearNum');
agg = sortrows(agg, 'YearNum');
Year = agg.YearNum;
Inc_obs = agg.mean_Incidence;
Mort_obs = agg.mean_Mortality;
CFR_obs = Mort_obs ./ Inc_obs;

% Clean missing or zero values
Inc_obs = fillmissing(Inc_obs,'linear','EndValues','nearest');
Mort_obs = fillmissing(Mort_obs,'linear','EndValues','nearest');
CFR_obs = Mort_obs ./ Inc_obs;
Inc_obs(Inc_obs<=0) = 1e-6;
CFR_obs(CFR_obs<=0) = 1e-6;

t0 = Year(1);
t = Year - t0;

%% --- Fit baseline param model ---
% I(t)=I0*exp(-r t); C(t)=C_inf+(C0-C_inf)*exp(-k t)
I0 = Inc_obs(1);
C0 = CFR_obs(1);
objfun = @(x) residuals_combined(x, t, Inc_obs, CFR_obs);

x0 = [0.02, 0.03, max(1e-4, min(CFR_obs))]; % initial [r, k, C_inf]
lb = [0, 0, 0];
ub = [1, 1, 1];

useFallback = false;
try
    opts = optimoptions('lsqnonlin', ...
        'Display','iter', ...
        'MaxIterations',500, ...
        'FunctionTolerance',1e-10);
catch
    warning('lsqnonlin not available; switching to fminsearch fallback.');
    useFallback = true;
end

if useFallback
    objfun_scalar = @(x) sum(objfun(x).^2);
    [xhat, resnorm] = fminsearch(objfun_scalar, x0);
    residual = []; exitflag = 1;
else
    [xhat, resnorm, residual, exitflag] = lsqnonlin(objfun, x0, lb, ub, opts);
end

r_hat = xhat(1);
k_hat = xhat(2);
Cinf_hat = xhat(3);

fprintf('\nFitted params:\n r = %.5f\n k = %.5f\n C_inf = %.5f\n', ...
    r_hat, k_hat, Cinf_hat);

%% --- Goodness of fit diagnostics
I_fit = I0 * exp(-r_hat * t);
C_fit = Cinf_hat + (C0 - Cinf_hat) * exp(-k_hat * t);

figure;
subplot(2,2,1);
plot(Year, Inc_obs, 'ro', Year, I_fit, 'b-');
xlabel('Year'); ylabel('Incidence'); legend('Observed','Fit'); title('Incidence Fit');
subplot(2,2,2);
plot(Year, CFR_obs, 'ro', Year, C_fit, 'b-');
xlabel('Year'); ylabel('CFR'); legend('Observed','Fit'); title('CFR Fit');
subplot(2,2,3);
plot(t, (I_fit - Inc_obs)./max(Inc_obs), 'k.-'); xlabel('Time'); ylabel('Residuals'); title('Incidence Residuals');
subplot(2,2,4);
plot(t, (C_fit - CFR_obs)./max(CFR_obs), 'k.-'); xlabel('Time'); ylabel('Residuals'); title('CFR Residuals');
saveas(gcf, 'fit_plots.png');
close(gcf);

fit_tbl = table(Year, Inc_obs, I_fit, CFR_obs, C_fit, 'VariableNames', ...
    {'Year','Incidence_obs','Incidence_fit','CFR_obs','CFR_fit'});
writetable(fit_tbl, 'model_fit_results.csv');
fprintf('Saved model_fit_results.csv and fit_plots.png\n');

%% --- Scenario forecasting
year_end = 2040;
tgrid_years = (t0:1:year_end)'; 
tgrid = tgrid_years - t0;
nT = numel(tgrid);
gamma1_default = 0.05;
gamma2_default = 0.05;

u0 = zeros(nT,1);
u1_enh = 0.5 * ones(nT,1);
u2_enh = zeros(nT,1);
u1_treat = zeros(nT,1);
u2_treat = 0.5 * ones(nT,1);
u1_comb = 0.5 * ones(nT,1);
u2_comb = 0.5 * ones(nT,1);

scenarios = { 'Baseline', u0, u0;
              'Enhanced_Prevention', u1_enh, u2_enh;
              'Enhanced_Treatment', u1_treat, u2_treat;
              'Combined_Controls', u1_comb, u2_comb };

elim_thresh = 4.0;
scenario_rows = [];
scenario_outputs = struct();

for i = 1:size(scenarios,1)
    name = scenarios{i,1};
    u1 = scenarios{i,2};
    u2 = scenarios{i,3};
    [I_sim, C_sim] = simulate_controls(r_hat, k_hat, Cinf_hat, I0, C0, gamma1_default, gamma2_default, u1, u2, nT);
    idx = find(I_sim < elim_thresh, 1, 'first');
    if isempty(idx)
        elim_year = NaN;
    else
        elim_year = tgrid_years(idx);
    end
    scenario_rows = [scenario_rows; {name, elim_year, I_sim(end), C_sim(end)}];
    scenario_outputs.(matlab.lang.makeValidName(name)) = table(tgrid_years, I_sim, C_sim);
end

scenario_summary = cell2table(scenario_rows, 'VariableNames', {'Scenario','EliminationYear','I_2040','C_2040'});
writetable(scenario_summary, 'scenario_summary.csv');
fprintf('Saved scenario_summary.csv\n');

%% --- Scenario Forecasts: INCIDENCE ---
figure; hold on;

figure('Position',[100 100 900 420]); % Taller figure to allow space underneath
hold on;

p1 = plot(scenario_outputs.Baseline.tgrid_years, scenario_outputs.Baseline.I_sim, 'b-', 'LineWidth', 2);
p2 = plot(scenario_outputs.Enhanced_Prevention.tgrid_years, scenario_outputs.Enhanced_Prevention.I_sim, 'r-', 'LineWidth', 2);
p3 = plot(scenario_outputs.Enhanced_Treatment.tgrid_years, scenario_outputs.Enhanced_Treatment.I_sim, 'k--', 'LineWidth', 2);
p4 = plot(scenario_outputs.Combined_Controls.tgrid_years, scenario_outputs.Combined_Controls.I_sim, 'm-', 'LineWidth', 2);

yline(elim_thresh,'k--','Elimination threshold (4 per 100k)');

xlabel('Year');
ylabel('Incidence per 100,000 women-years');
title('Projected Incidence under Intervention Scenarios');

% ---- Move legend below plot, horizontal, centered ----
lgd = legend([p1 p2 p3 p4], ...
    {'Baseline (no intervention)', ...
     'Enhanced Prevention (u_1 increased)', ...
     'Enhanced Treatment (u_2 increased)', ...
     'Combined Controls (u_1 and u_2 increased)'}, ...
     'Location','southoutside', ...
     'Orientation','horizontal', ...
     'FontSize',11, ...        % smaller font so it fits
     'NumColumns',2);          % wrap legend into 2 columns

lgd.Box = 'off';

% ---- Resize the axes upward so the legend has space ----
ax = gca;
ax.Position = [0.10 0.38 0.85 0.55]; 
% [left bottom width height]

% Reduce white padding around axes
set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02));
saveas(gcf,'scenario_incidence.png');
close(gcf);


%% --- Sensitivity on gamma1
gamma1_values = linspace(0,0.2,11);
sens_rows = [];
for gix = 1:numel(gamma1_values)
    g1 = gamma1_values(gix);
    [I_s, ~] = simulate_controls(r_hat, k_hat, Cinf_hat, I0, C0, g1, gamma2_default, u1_comb, u2_comb, nT);
    idx = find(I_s < elim_thresh, 1, 'first');
    if isempty(idx), ey = NaN; else ey = tgrid_years(idx); end
    sens_rows = [sens_rows; {g1, ey, I_s(end)}];
end
sens_tbl = cell2table(sens_rows, 'VariableNames', {'gamma1','EliminationYear','I_2040'});
writetable(sens_tbl, 'gamma1_sensitivity.csv');
fprintf('Saved gamma1_sensitivity.csv\n');

figure;
elim_plot = sens_tbl.EliminationYear;
elim_plot(isnan(elim_plot)) = 2050;
plot(sens_tbl.gamma1, elim_plot, '-o');
xlabel('gamma1'); ylabel('Elimination year (2050=not reached)');
title('Sensitivity: gamma1 vs elimination year');
saveas(gcf,'sensitivity_gamma1.png');
close(gcf);

fprintf('--- DONE ---\n');
disp('Generated files:');
disp('  model_fit_results.csv');
disp('  fit_plots.png');
disp('  scenario_summary.csv');
disp('  scenario_incidence.png');
disp('  scenario_cfr.png');
disp('  gamma1_sensitivity.csv');
disp('  sensitivity_gamma1.png');

%% --- Helper Functions ---
function res = residuals_combined(x, t, Inc_obs, CFR_obs)
    r = x(1); k = x(2); Cinf = x(3);
    I0 = Inc_obs(1); C0 = CFR_obs(1);
    I_pred = I0 * exp(-r*t);
    C_pred = Cinf + (C0 - Cinf)*exp(-k*t);
    rI = (I_pred - Inc_obs)/max(Inc_obs);
    rC = (C_pred - CFR_obs)/max(CFR_obs);
    res = [rI(:); rC(:)];
end

function [I, C] = simulate_controls(r, k, Cinf, I0, C0, gamma1, gamma2, u1, u2, nT)
    I = zeros(nT,1); C = zeros(nT,1);
    I(1) = I0; C(1) = C0;
    dt = 1.0;
    for tt = 1:(nT-1)
        u1_t = u1(tt); u2_t = u2(tt);
        f = @(y,uu1,uu2) [ -(r + gamma1*uu1)*y(1);
                           -(k + gamma2*uu2)*(y(2) - Cinf) ];
        y = [I(tt); C(tt)];
        k1 = f(y,u1_t,u2_t);
        k2 = f(y + 0.5*dt*k1,u1_t,u2_t);
        k3 = f(y + 0.5*dt*k2,u1_t,u2_t);
        k4 = f(y + dt*k3,u1_t,u2_t);
        ynext = y + dt*(k1+2*k2+2*k3+k4)/6;
        I(tt+1) = max(ynext(1),1e-12);
        C(tt+1) = max(ynext(2),1e-12);
    end
end

