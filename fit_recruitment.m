%fit data to get recruitment stretch
% clear all
% clc
%% Obstructed bladder
%Ob01
ObDp{1} = [0 12 61 69 84.6];
Oblp{1} = [0 0 0 2 87];
ObDs{1} = [0.0165 0.2410 0.2580 0.3365 0.3830];
Obls{1} = [0 0.1810 0.2515 0.3080 0.4215];
%Ob02
ObDp{2} = [0 2 21 45 63];
Oblp{2} = [0 0 0 0 52];
ObDs{2} = [0.0620 0.1795 0.2770 0.3460 0.4195];
Obls{2} = [0.0510 0.2645 0.3190 0.3330 0.3935];
%Ob03
ObDp{3} = [0 44 61 71 79];
Oblp{3} = [0 0 0 0 7];
ObDs{3} = [0.1220 0.2595 0.2960 0.3100 0.3185];
Obls{3} = [0.1095 0.1855 0.2455 0.2800 0.3110];
%Ob04
ObDp{4} = [0 5 44 46 74];
Oblp{4} = [0 0 0 8 85];
ObDs{4} = [0.0225 0.1205 0.2285 0.2905 0.3380];
Obls{4} = [0 0.0300 0.1480 0.2200 0.3380]; % 0.2200 0.3245
% %% Young bladder
% %Y01
% YDp{1} = [0 0 48 65 72];
% Ylp{1} = [0 0 0 28];
% YDs{1} = [0.0410 0.1150 0.2575 0.3325 0.3785];
% Yls{1} = [0.1515 0.2140 0.2500 0.3070 0.4170];
% %Y02
% YDp{2} = [35 70 82 92 100];
% Ylp{2} = [0 0 0 0 17];
% YDs{2} = [0.0740 0.2065 0.2985 0.3710 0.4440];
% Yls{2} = [0.0810 0.2895 0.3325 0.3570 0.4185];
% %Y03
% YDp{3} = [0 23 64 73 85];
% Ylp{3} = [0 0 0 0 18];
% YDs{3} = [0.3250 0.3580 0.3945 0.4075 0.4175];
% Yls{3} = [0.2100 0.3155 0.3435 0.3775 0.4085];
% %Y04
% YDp{4} = [0 6 44 60 68];
% Ylp{4} = [0 0 0 0 31];
% YDs{4} = [0.0275 0.1170 0.2235 0.24 0.3355];
% Yls{4} = [0 0.0240 0.1450 0.2140 0.3205];




%% Ob detrusor
for i = 1:1:4
    for j = 1:length(ObDs{i})
        ObDlambda{i}(j) = sqrt(ObDs{i}(j)*2+1);
        if j == 1
        ObDdlam{i}(j) =  ObDlambda{i}(j)-1;
        else 
         ObDdlam{i}(j) =  ObDlambda{i}(j)-ObDlambda{i}(j-1);   
    end
    end
end
%% Ob lamina propria
for i = 1:1:4
    for j = 1:length(Obls{i})
        Obllambda{i}(j) = sqrt(Obls{i}(j)*2+1);
        if j == 1
        Obldlam{i}(j) =  Obllambda{i}(j)-1;
        else 
         Obldlam{i}(j) =  Obllambda{i}(j)-Obllambda{i}(j-1);   
    end
    end
end

%%
for i = 1:1:4
%% fit the functon DT
DOX = ObDlambda{i};
DOY = ObDp{i}./100;
%Fit the data using fminsearch
cdf3 = @(r,lam) triangular_fit(r,lam);
cost = @(r)norm(DOY - cdf3(r,DOX))
DOR{i} = fminsearch(cost, [ObDlambda{i}(1); ObDlambda{i}(3); ObDlambda{i}(2)])
pd3 = makedist('Triangular','a',DOR{i}(1),'b',DOR{i}(3),'c',DOR{i}(2));
%Plot cdf and pdf
Xplot = linspace(1, 2);
pdf3 = pdf(pd3,Xplot);
figure(i+4)
plot(DOX, DOY,'rs','MarkerFaceColor','r','MarkerSize',18)
hold on
plot(Xplot, cdf3(DOR{i},Xplot),'r','LineWidth',4)
xlim([1,2])
xlabel('Stretch')
ylabel('Collagen Fiber Recruitment Fraction')
set(gca,'fontsize',15)
set(gca,'box','off')

figure(i+10)
hold on
plot(Xplot,pdf3,'r','LineWidth',4)
xlim([1,2])
xlabel('Stretch')
ylabel('Probability density')
set(gca,'fontsize',15)


%% fit the functon LP
LOX = Obllambda{i};
LOY = Oblp{i}./100;
%Fit the data using fminsearch
cdf3 = @(r,lam) triangular_fit(r,lam);
cost = @(r)norm(LOY - cdf3(r,LOX))
LOR{i} = fminsearch(cost, [Obllambda{i}(1); Obllambda{i}(3); Obllambda{i}(2)])
pd3 = makedist('Triangular','a',LOR{i}(1),'b',LOR{i}(3),'c',LOR{i}(2));
Xplot = linspace(1, 2);
pdf3 = pdf(pd3,Xplot);
%Plot cdf and pdf
figure(i+4)
plot(LOX, LOY, 'bs','MarkerFaceColor','b','MarkerSize',18)
hold on
plot(Xplot, cdf3(LOR{i},Xplot),'b','LineWidth',4)
xlim([1,2])
xlabel('Stretch')
ylabel('Collagen Fiber Recruitment Fraction')
set(gca,'fontsize',15)
set(gca,'box','off')

figure(i+10)
plot(Xplot,pdf3,'b','LineWidth',4)
xlim([1,2])
xlabel('Stretch')
ylabel('Probability density')
set(gca,'fontsize',20)
legend('DSM','LP')
end

%fit data to get recruitment stretch
% clear all
% clc
%% Young bladder
%Y01
YDp{1} = [0 0 48 65 72];
Ylp{1} = [0 0 0 28 50]; %31
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
figure(i+20)
plot(DOX, DOY, 'rs','MarkerFaceColor','r','MarkerSize',18)
hold on
plot(Xplot, cdf3(DYR{i},Xplot),'r','LineWidth',4)
xlim([1,2])
xlabel('Stretch')
ylabel('Collagen Fiber Recruitment Fraction')
set(gca,'fontsize',15)
set(gca,'box','off')
figure(i+30)
hold on
plot(Xplot,pdf3,'r','LineWidth',4)
xlim([1,2])
xlabel('Stretch')
ylabel('Probability density')
set(gca,'fontsize',15)
legend('DSM','LP')

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
figure(i+20)
plot(LOX, LOY, 'bs','MarkerFaceColor','b','MarkerSize',18)
hold on
plot(Xplot, cdf3(LYR{i},Xplot),'b','LineWidth',4)
xlim([1,2])
xlabel('Stretch')
ylabel('Collagen Fiber Recruitment Fraction')
set(gca,'fontsize',15)
set(gca,'box','off')

figure(i+30)
plot(Xplot,pdf3,'b','LineWidth',4)
xlim([1,2])
xlabel('Stretch')
ylabel('Probability density')
set(gca,'fontsize',20)
legend('DSM','LP')
end

save recruitment_stretch.mat DYR LYR DOR LOR


