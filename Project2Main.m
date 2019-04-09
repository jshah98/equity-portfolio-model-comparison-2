%% MIE377 (Winter 2019) - Project 1
% The purpose of this program is to implement two different stochastic 
% processes to simulate scenarios: 
% (i) a Gaussian process, and 
% (ii) a non-normal stochastic process with higher moments
% 
% and to use these simulations to test the out-of-sample performance of two
% portfolio optimization models: 
% 1. CVaR
% 2. Robust CVaR
%
% Student Name: Joshua Chang, Jinansh Shah, Qingyang Li
% Student ID: 1003083147, 1003062614, 1002741063



clc
clear all
format short
warning('off','all')

% Program Start
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Read input files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



adjClose = readtable('Project1_Data_adjClose.csv');
adjClose.Properties.RowNames = cellstr(datetime(adjClose.Date));
adjClose.Date = [];

% Load the factors weekly returns
factorRet = readtable('Project1_Data_FF_factors.csv');
factorRet.Properties.RowNames = cellstr(datetime(factorRet.Date));
factorRet.Date = [];

riskFree = factorRet(:,4);
factorRet = factorRet(:,1:3);

% Identify the tickers and the dates 
tickers = adjClose.Properties.VariableNames';
dates   = datetime(factorRet.Properties.RowNames);

% Calculate the stocks' weekly EXCESS returns
prices  = table2array(adjClose);
returns = ( prices(2:end,:) - prices(1:end-1,:) ) ./ prices(1:end-1,:);
returns = returns - ( diag( table2array(riskFree) ) * ones( size(returns) ) );
returns = array2table(returns);
returns.Properties.VariableNames = tickers;
returns.Properties.RowNames = cellstr(datetime(factorRet.Properties.RowNames));

riskFree = table2array(riskFree);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Define your initial parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initial budget to invest
initialVal = 100;

% Start of in-sample calibration period 
calStart = datetime('2012-01-01');
calEnd   = calStart + calmonths(12) - days(1);

% Start of out-of-sample test period 
testStart = datetime('2013-01-01');
testEnd   = testStart + calmonths(6) - days(1);

% Number of investment periods (each investment period is 6 months long)
NoPeriods = 6;

% Monte-Carlo simulations 
MCList = {'MC_Gaussian' 'MC_HM'};

MCList = cellfun(@str2func, MCList, 'UniformOutput', false);
NoSimulations = length(MCList);


% Investment strategies
invList = {'CVaR' 'CVaR_robust_box'};

invList = cellfun(@str2func, invList, 'UniformOutput', false);
NoStrats = length(invList);

n = 20; % number of assets
S = 5000; % number of scenarios


% Tags for factor models under different investment strategies
tags = {'CVaR (Gaussian)' 'Robust CVaR (Gaussian)' 'CVaR (Non-normal)' ...
    'Robust CVaR (Non-normal)'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Construct and rebalance your portfolios
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initiate counter for the number of observations per investment period,
% transaction cost variable, and current value of the portfolio.
currentVal = zeros(NoPeriods, NoStrats * NoSimulations); 
toDay = 0;
pers = NaT(NoPeriods, 1);

exp_r = ones (NoPeriods, NoSimulations*NoStrats); 
realized_std_dev = ones (NoPeriods, NoSimulations*NoStrats); 




for t = 1 : NoPeriods
    
    % Subset the returns and factor returns corresponding to the current
    % calibration period.
    periodReturns = table2array( returns( calStart <= dates & dates <= calEnd, :) );
    periodFactRet = table2array( factorRet( calStart <= dates & dates <= calEnd, :) );
    currentPrices = table2array( adjClose( ( calEnd - days(7) ) <= dates ... 
                                                    & dates <= calEnd, :) )';
    
    % Subset the prices corresponding to the current out-of-sample test 
    % period.
    periodPrices = table2array( adjClose( testStart <= dates & dates <= testEnd,:) );
    
    % Set the initial value of the portfolio or update the portfolio value
    if t == 1
        currentVal= ones(1, 6);
        currentVal(t,:) = initialVal;
        
        
    else
        for k = 1 : (NoStrats * NoSimulations)
            temp = transpose(NoShares{k})* currentPrices;
            currentVal(t,k) = temp;
            
        end
    end
    
    % Update counter for the number of observations per investment period
    fromDay = toDay + 1;
    toDay   = toDay + size(periodPrices,1);
    
    [alpha, beta, D, mu, Q{t}] = FF(periodReturns, periodFactRet);
    
    for i = 1:NoSimulations
        rets{i} = MCList{i}(alpha, beta, D, periodFactRet); 
        for j = 1:NoStrats
            k = j + (i - 1) * 2;
            x{k}(:,t) = invList{j}(rets{i}');
            exp_r(t, k) = mu'*x{k}(:,t);
        end
    end
    
    
    % Calculate the optimal number of shares of each stock you should hold
    for k = 1 : NoStrats * NoSimulations
        % Number of shares your portfolio holds per stock
        NoShares{k} = x{k}(:,t) .* currentVal(t,k) ./ currentPrices;
        
        % Weekly portfolio value during the out-of-sample window
        portfValue(fromDay:toDay,k) = periodPrices * NoShares{k};
        
        NoSharesOld{k} = NoShares{k};
        
        
        realized_std_dev(t, k) = std(periodReturns*x{k}(:, t));
    end

    % Update your calibration and out-of-sample test periods
    calStart = calStart + calmonths(6);
    calEnd   = calStart + calmonths(12) - days(1);
    pers(t) = testStart; 
    testStart = testStart + calmonths(6);
    testEnd   = testStart + calmonths(6) - days(1);

end



%--------------------------------------------------------------------------
% 4.1 Analysis of realized returns/variance
%--------------------------------------------------------------------------

% Initializing portfolio returns and variance to 0 - to be calculated later
plotDates = dates(dates >= datetime('2013-01-01'));
portf_rets = zeros (length(plotDates)-1, NoStrats * NoSimulations);
realized_returns = zeros(NoStrats * NoSimulations, 1);
exp_sd = zeros(NoPeriods, NoStrats * NoSimulations);

%Initailizing mean for geometric mean of portfolio return for each period
portfolioPeriodReturns = zeros(6, NoSimulations * NoStrats);
%Initalizing Sharpe Ratio Matrix
sharpeRatio = zeros(6, NoSimulations * NoStrats);

rf = zeros(6, 1);

rf(1) = (geomean(riskFree(53:52+26)+1)^52-1);
rf(2) = (geomean(riskFree(52+26+1:52+52)+1)^52-1);
rf(3) = (geomean(riskFree(52+52+1:52+78)+1)^52-1);
rf(4) = (geomean(riskFree(52+79: 52+104)+1)^52-1);
rf(5) = (geomean(riskFree(52+105:52+130)+1)^52-1);
rf(6) = (geomean(riskFree(52+131:52+156)+1)^52-1);
    


mean_rsd = mean(realized_std_dev) ;

% Calculating annualized return of each portfolio over the 3 years.
for k = 1 : NoSimulations * NoStrats
    
    % Calculating the weekly returns 
    portf_rets(:, k) =  (portfValue(2:(length(plotDates)),k)- portfValue(1:(length(plotDates)-1),k)) ./ portfValue(1:(length(plotDates)-1),k) + ones(156, 1) ; 
    % Converting weekly returns to annualized returns
    realized_returns(k) = (geomean(portf_rets(:, k)))^52 -1;
    
    sharpeRatio(:, k) = ((exp_r(:, k)+1).^52-1 - rf(:, 1))./(sqrt(52)*mean_rsd(k));
    
    
    
    exp_r(1, k) = (geomean(exp_r(:, k) + 1))^52 -1;

    % Calculating expected portfolio variance for each portfolio
    for t = 1:NoPeriods
        exp_sd(t, k) = sqrt(x{k}(:, t)' * Q{t} * x{k}(:, t));
    end
    

    
end


tags_exp_real = {'CVaR (Gaussian) Expected' 'CVaR (Gaussian) Realized'...
    'Robust CVaR (Gaussian) Expected' 'Robust CVaR (Gaussian) Realized'  ...
    'CVaR (Non-normal) Expected' 'CVaR (Non-normal) Realized' ...
    'Robust CVaR (Non-normal) Expected' 'Robust CVaR (Non-normal) Realized'};

fig11 = figure(11);
plot( pers, exp_sd(:,1), 'b--')
hold on
plot( pers, realized_std_dev(:,1), 'b-')
hold on
plot( pers, exp_sd(:,2), 'r--')
hold on
plot( pers, realized_std_dev(:,2), 'r-')
hold on
plot( pers, exp_sd(:,3), 'g--')
hold on
plot( pers, realized_std_dev(:, 3), 'g-')
hold on
plot( pers, exp_sd(:,4), 'k--')
hold on
plot( pers, realized_std_dev(:,4), 'k-')
hold on

legend(tags_exp_real, 'Location', 'eastoutside','FontSize',12);
datetick('x','dd-mmm-yyyy','keepticks','keeplimits');
set(gca,'XTickLabelRotation',30);
title('Expected & Realized Volatility', 'FontSize', 14)
ylabel('Volatility','interpreter','latex','FontSize',12);

% Define the plot size in inches
set(fig11,'Units','Inches', 'Position', [0 0 8, 5]);
pos1 = get(fig11,'Position');
set(fig11,'PaperPositionMode','Auto','PaperUnits','Inches',...
    'PaperSize',[pos1(3), pos1(4)]);

% If you want to save the figure as .pdf for use in LaTeX
% print(fig1,'fileName','-dpdf','-r0');

% If you want to save the figure as .png for use in MS Word
print(fig11,'fileName11','-dpng','-r0');

%--------------------------------------------------------------------------
% 4.2 Plot the portfolio values 
% 
% Note: The code below plots all portfolios onto a single plot. However,
% you may want to split this into multiple plots for clarity, or to
% compare a subset of the portfolios. 
%--------------------------------------------------------------------------
plotDates = dates(dates >= datetime('2013-01-01'));

fig1 = figure(1);

for k = 1 : NoSimulations * NoStrats
    
    plot( plotDates, portfValue(:,k))
    hold on
    
end

legend(tags, 'Location', 'eastoutside','FontSize',12);
datetick('x','dd-mmm-yyyy','keepticks','keeplimits');
set(gca,'XTickLabelRotation',30);
title('Portfolio value', 'FontSize', 14)
ylabel('Value','interpreter','latex','FontSize',12);

% Define the plot size in inches
set(fig1,'Units','Inches', 'Position', [0 0 8, 5]);
pos1 = get(fig1,'Position');
set(fig1,'PaperPositionMode','Auto','PaperUnits','Inches',...
    'PaperSize',[pos1(3), pos1(4)]);

% If you want to save the figure as .pdf for use in LaTeX
% print(fig1,'fileName','-dpdf','-r0');

% If you want to save the figure as .png for use in MS Word
print(fig1,'fileName','-dpng','-r0');

%--------------------------------------------------------------------------
% 4.3 Plot the portfolio weights 
%--------------------------------------------------------------------------

% CVaR (Gaussian) Plot
fig2 = figure(2);
area(x{1}')
legend(tickers, 'Location', 'eastoutside','FontSize',12);
title('CVaR (Gaussian) portfolio weights', 'FontSize', 14)
ylabel('Weights','interpreter','latex','FontSize',12);
xlabel('Rebalance period','interpreter','latex','FontSize',12);

% Define the plot size in inches
set(fig2,'Units','Inches', 'Position', [0 0 8, 5]);
pos1 = get(fig2,'Position');
set(fig2,'PaperPositionMode','Auto','PaperUnits','Inches',...
    'PaperSize',[pos1(3), pos1(4)]);

% If you want to save the figure as .pdf for use in LaTeX
% print(fig2,'fileName2','-dpdf','-r0');

% If you want to save the figure as .png for use in MS Word
print(fig2,'fileName2','-dpng','-r0');



% Robust CVaR (Gaussian) Plot
fig3 = figure(3);
area(x{2}')
legend(tickers, 'Location', 'eastoutside','FontSize',12);
title('Robust CVaR (Gaussian) portfolio weights', 'FontSize', 14)
ylabel('Weights','interpreter','latex','FontSize',12);
xlabel('Rebalance period','interpreter','latex','FontSize',12);

% Define the plot size in inches
set(fig3,'Units','Inches', 'Position', [0 0 8, 5]);
pos1 = get(fig3,'Position');
set(fig3,'PaperPositionMode','Auto','PaperUnits','Inches',...
    'PaperSize',[pos1(3), pos1(4)]);

% If you want to save the figure as .pdf for use in LaTeX
% print(fig2,'fileName2','-dpdf','-r0');

% If you want to save the figure as .png for use in MS Word
print(fig3,'fileName2','-dpng','-r0');   

% Robust CVaR (Non-normal) Plot
fig4 = figure(4);
area(x{3}')
legend(tickers, 'Location', 'eastoutside','FontSize',12);
title('CVaR (Non-normal) portfolio weights', 'FontSize', 14)
ylabel('Weights','interpreter','latex','FontSize',12);
xlabel('Rebalance period','interpreter','latex','FontSize',12);

% Define the plot size in inches
set(fig4,'Units','Inches', 'Position', [0 0 8, 5]);
pos1 = get(fig4,'Position');
set(fig4,'PaperPositionMode','Auto','PaperUnits','Inches',...
    'PaperSize',[pos1(3), pos1(4)]);

% If you want to save the figure as .pdf for use in LaTeX
% print(fig2,'fileName2','-dpdf','-r0');

% If you want to save the figure as .png for use in MS Word
print(fig4,'fileName2','-dpng','-r0'); 



% Robust CVaR (Gaussian) Plot
fig5 = figure(5);
area(x{4}')
legend(tickers, 'Location', 'eastoutside','FontSize',12);
title('Robust CVaR (Non-normal) portfolio weights', 'FontSize', 14)
ylabel('Weights','interpreter','latex','FontSize',12);
xlabel('Rebalance period','interpreter','latex','FontSize',12);

% Define the plot size in inches
set(fig5,'Units','Inches', 'Position', [0 0 8, 5]);
pos1 = get(fig5,'Position');
set(fig5,'PaperPositionMode','Auto','PaperUnits','Inches',...
    'PaperSize',[pos1(3), pos1(4)]);

% If you want to save the figure as .pdf for use in LaTeX
% print(fig2,'fileName2','-dpdf','-r0');

% If you want to save the figure as .png for use in MS Word
print(fig5,'fileName2','-dpng','-r0');   


%--------------------------------------------------------------------------
% 4.4 Plot of loss distribution
%--------------------------------------------------------------------------


loss = -rets{1}' * x{1};

fig6 = figure(6);
histogram(loss, 50);
xlabel('Portfolio losses ($\%$)','interpreter',...
    'latex','FontSize',14);
ylabel('Frequency','interpreter','latex','FontSize',14); 

set(fig6,'Units','Inches', 'Position', [0 0 10, 4]);
pos2 = get(fig6,'Position');
set(fig6,'PaperPositionMode','Auto','PaperUnits','Inches',...
'PaperSize',[pos2(3), pos2(4)])



loss = -rets{1}' * x{2};

fig7 = figure(7);
histogram(loss, 50);
xlabel('Portfolio losses ($\%$)','interpreter',...
    'latex','FontSize',14);
ylabel('Frequency','interpreter','latex','FontSize',14); 

set(fig7,'Units','Inches', 'Position', [0 0 10, 4]);
pos2 = get(fig7,'Position');
set(fig7,'PaperPositionMode','Auto','PaperUnits','Inches',...
'PaperSize',[pos2(3), pos2(4)])


loss = -rets{2}' * x{3};

fig8 = figure(8);
histogram(loss, 50);
xlabel('Portfolio losses ($\%$)','interpreter',...
    'latex','FontSize',14);
ylabel('Frequency','interpreter','latex','FontSize',14); 

set(fig8,'Units','Inches', 'Position', [0 0 10, 4]);
pos2 = get(fig8,'Position');
set(fig8,'PaperPositionMode','Auto','PaperUnits','Inches',...
'PaperSize',[pos2(3), pos2(4)])


loss = -rets{2}' * x{4};

fig9 = figure(9);
histogram(loss, 50);
xlabel('Portfolio losses ($\%$)','interpreter',...
    'latex','FontSize',14);
ylabel('Frequency','interpreter','latex','FontSize',14); 

set(fig9,'Units','Inches', 'Position', [0 0 10, 4]);
pos2 = get(fig9,'Position');
set(fig9,'PaperPositionMode','Auto','PaperUnits','Inches',...
'PaperSize',[pos2(3), pos2(4)])

