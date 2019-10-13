%fit data to get recruitment stretch
% clear all
% clc
%% Young bladder
%Y01
YDp{1} = [0 0 48 65 72];
Ylp{1} = [0 0 0 28 80]; %31
YDs{1} = [0 0.0825 0.1550 0.2090 0.2555];
Yls{1} = [0 0.0800 0.1485 0.1965 0.2840];
%Y02
YDp{2} = [0 70 82 92 100];
Ylp{2} = [0 0 0 0 17];
YDs{2} = [0.0995 0.2215 0.2385 0.2575 0.2805];
Yls{2} = [0.0385 0.0560 0.1800 0.2305 0.3140];
%Y03
YDp{3} = [0 23 64 73 85];
Ylp{3} = [0 0 0 0 20];
YDs{3} = [0.0490 0.1485 0.1895 0.2715 0.2850];
Yls{3} = [0.0445 0.1465 0.1855 0.2340 0.2630];
%Y04
YDp{4} = [0 6 44 60 68];
Ylp{4} = [0 0 0 0 31]; %25
YDs{4} = [0 0.1010 0.1730 0.2400 0.3015]; % 0 0.1000 0.1740 0.2400 0.3010
Yls{4} = [0 0.0970 0.1510 0.2110 0.3080]; % 0 0.0970 0.1545 0.2105 0.3080

%% Y detrusor
for i = 1:1:4
    for j = 1:length(YDs{i})
        YDlambda{i}(j) = sqrt(YDs{i}(j)*2+1);
        if j == 1
        YDdlam{i}(j) =  YDlambda{i}(j)-1;
        else 
         YDdlam{i}(j) =  YDlambda{i}(j)-YDlambda{i}(j-1);   
    end
    end
end
%% Y lamina propria
for i = 1:1:4
    for j = 1:length(Yls{i})
        Yllambda{i}(j) = sqrt(Yls{i}(j)*2+1);
        if j == 1
        Yldlam{i}(j) =  Yllambda{i}(j)-1;
        else 
         Yldlam{i}(j) =  Yllambda{i}(j)-Yllambda{i}(j-1);   
    end
    end
end

%%
for i = 1:1:4
%% fit the functon DT
DOX = YDlambda{i};
DOY = YDp{i}./100;
%Fit the data using fminsearch
cdf3 = @(r,lam) triangular_fit(r,lam);
cost = @(r)norm(DOY - cdf3(r,DOX))
DYR{i} = fminsearch(cost, [YDlambda{i}(1); YDlambda{i}(3); YDlambda{i}(2)])
pd3 = makedist('Triangular','a',DYR{i}(1),'b',DYR{i}(3),'c',DYR{i}(2));
%Plot cdf and pdf
Xplot = linspace(1, 2);
pdf3 = pdf(pd3,Xplot);
figure(i+4)
plot(DOX, DOY, 'or','MarkerSize',20)
hold on
plot(Xplot, cdf3(DYR{i},Xplot),'r','LineWidth',2)
% grid
xlim([1,2])
xlabel('Stretch')
ylabel('Collagen Fiber Recruitment Fraction')
set(gca,'fontsize',15)
figure(i+10)
hold on
plot(Xplot,pdf3,'r','LineWidth',2)
grid
xlim([1,2])
xlabel('Stretch')
ylabel('Probability density')
set(gca,'fontsize',15)


%% fit the functon LP
LOX = Yllambda{i};
LOY = Ylp{i}./100;
%Fit the data using fminsearch
cdf3 = @(r,lam) triangular_fit(r,lam);
cost = @(r)norm(LOY - cdf3(r,LOX))
LYR{i} = fminsearch(cost, [Yllambda{i}(1); Yllambda{i}(3); Yllambda{i}(2)])
pd3 = makedist('Triangular','a',LYR{i}(1),'b',LYR{i}(3),'c',LYR{i}(2));
Xplot = linspace(1, 2);
pdf3 = pdf(pd3,Xplot);
%Plot cdf and pdf
figure(i+4)
plot(LOX, LOY, 'sb','MarkerSize',10)
hold on
plot(Xplot, cdf3(LYR{i},Xplot),'b','LineWidth',2)
grid
xlim([1,2])
xlabel('Collagen recruitment stretch')
ylabel('Cumulative probability density')
set(gca,'fontsize',15)

figure(i+10)
plot(Xplot,pdf3,'b','LineWidth',2)
grid
xlim([1,2])
xlabel('Stretch')
ylabel('Probability density')
set(gca,'fontsize',15)
legend('DSM','LP')
end

% save recruitment_stretch.mat DYR LYR DOR LOR

