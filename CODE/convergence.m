%% -------------------- Step 1: Load dataset --------------------
[filename, pathname] = uigetfile({'*.xlsx','Excel Files (*.xlsx)'}, ...
                                 'Select cervical_cancer_tidy_numeric.xlsx file');
if isequal(filename,0)
    error(' No file selected. Please rerun and choose the Excel file.');
else
    fullpath = fullfile(pathname, filename);
    fprintf(' File selected: %s\n', fullpath);
end
data = readtable(fullpath);
years   = data.YearNum;
% Aggregate over AgeGroups (mean values)
allData = groupsummary(data, "YearNum", "mean", ["Incidence","Mortality"]);
years   = allData.YearNum;
obsInc  = allData.mean_Incidence;
obsMort = allData.mean_Mortality;
obsCFR  = obsMort ./ obsInc;
% Handle NaNs / stability
obsInc(obsInc <= 0 | isnan(obsInc))   = 1e-6;
obsMort(obsMort <= 0 | isnan(obsMort)) = 1e-6;
obsCFR(obsCFR <= 0 | isnan(obsCFR))   = 1e-6;
%% -------------------- Step 2: Model parameters --------------------
t0 = years(1); T = years(end); N = 200; h = (T-t0)/N;
r = 0.02; k = 0.03; Cinf = min(obsCFR);   % epidemiological
gamma1 = 0.5; gamma2 = 0.5;             % control effectiveness
A = 1; B1 = 0.1; B2 = 0.1; P = 10; eps = 4e-5;
I0 = max(obsInc(1),1e-6);
C0 = max(obsCFR(1),1e-6);
u1max = 1; u2max = 1;
tol = 1e-6; max_iter = 500; theta = 0.5;
t = linspace(t0,T,N+1);
u1 = 0.1*ones(1,N+1);
u2 = 0.1*ones(1,N+1);
errControlHist = zeros(max_iter,1);
errStateHist   = zeros(max_iter,1);
I_prev = zeros(1,N+1);
C_prev = zeros(1,N+1);
%% -------------------- Step 3: Forward-Backward Sweep --------------------
for iter = 1:max_iter
    % --- Forward pass ---
    I = zeros(1,N+1); C = zeros(1,N+1);
    I(1) = I0; C(1) = C0;
    for i = 1:N
        y = [I(i); C(i)];
        ui = [u1(i); u2(i)];
        f = @(y,u)[ -(r + gamma1*u(1))*y(1);
                     -(k + gamma2*u(2))*(y(2)-Cinf) ];
        k1 = f(y,ui);
        k2 = f(y+h/2*k1,ui);
        k3 = f(y+h/2*k2,ui);
        k4 = f(y+h*k3,ui);
        y_next = y + h/6*(k1+2*k2+2*k3+k4);
        I(i+1) = max(y_next(1),1e-6);
        C(i+1) = max(y_next(2),1e-6);
    end
    
    % --- Backward pass ---
    M_T = I(end)*C(end); dPhi = (M_T>eps)*2*(M_T-eps);
    lambdaI = zeros(1,N+1); lambdaC = zeros(1,N+1);
    lambdaI(end) = P*dPhi*C(end);
    lambdaC(end) = P*dPhi*I(end);
    
    for i = N:-1:1
        y = [I(i);C(i)];
        lam = [lambdaI(i+1);lambdaC(i+1)];
        ui = [u1(i);u2(i)];
        g = @(y,u,lam)[ -(A*y(2)-lam(1)*(r+gamma1*u(1)));
                         -(A*y(1)-lam(2)*(k+gamma2*u(2))) ];
        k1 = g(y,ui,lam);
        k2 = g(y,ui,lam-h/2*k1);
        k3 = g(y,ui,lam-h/2*k2);
        k4 = g(y,ui,lam-h*k3);
        lam_prev = lam - h/6*(k1+2*k2+2*k3+k4);
        lambdaI(i) = lam_prev(1);
        lambdaC(i) = lam_prev(2);
    end
    
    % --- Control update ---
    u1_new = (gamma1/B1).*lambdaI.*I;
    u2_new = (gamma2/B2).*lambdaC.*(C-Cinf);
    u1_new = min(u1max,max(0,u1_new));
    u2_new = min(u2max,max(0,u2_new));
    
    u1_prev = u1; u2_prev = u2;
    u1 = theta*u1_new + (1-theta)*u1_prev;
    u2 = theta*u2_new + (1-theta)*u2_prev;
    
    % --- Record convergence errors ---
    errControl = max([norm(u1-u1_prev,inf), norm(u2-u2_prev,inf)]);
    errState   = max([norm(I-I_prev,inf), norm(C-C_prev,inf)]);
    errControlHist(iter) = errControl;
    errStateHist(iter) = errState;
    
    I_prev = I; C_prev = C;
    
    fprintf('Iter %d: max control change = %.3e, max state change = %.3e\n', ...
        iter, errControl, errState);
    
    if errControl<tol && errState<tol
        disp('Method converged!');
        errControlHist = errControlHist(1:iter);
        errStateHist   = errStateHist(1:iter);
        break;
    end
end
%% -------------------- Step 4: Plot results --------------------
figure('Name','States and Controls','Position',[100 100 1000 600]);
subplot(2,2,1);
plot(t,I,'b','LineWidth',2); hold on;
scatter(years, obsInc, 50, 'r','filled');
xlabel('Year'); ylabel('Incidence'); legend('Model','Observed');
subplot(2,2,2);
plot(t,C,'r','LineWidth',2); hold on;
scatter(years, obsCFR, 50, 'k','filled');
xlabel('Year'); ylabel('CFR'); legend('Model','Observed');
subplot(2,2,3);
plot(t,u1,'g','LineWidth',2);
xlabel('Year'); ylabel('u1 (Prevention)');
subplot(2,2,4);
plot(t,u2,'m','LineWidth',2);
xlabel('Year'); ylabel('u2 (Treatment)');
%% -------------------- Step 5: Convergence plot --------------------
figure('Name','Convergence','Position',[200 200 700 400]);
semilogy(1:length(errControlHist), errControlHist,'b-o','LineWidth',2); hold on;
semilogy(1:length(errStateHist), errStateHist,'r-s','LineWidth',2);
xlabel('Iteration'); ylabel('Max change');
title('Convergence of Controls and States');
legend('Control change','State change');
grid on;

