%fit for the collagenstiffness parameter
clear all
clc
% load the recruitment stretch: min mod max
load recruitment_stretch.mat
%Plot Ob tissue
for i = 1:1:4
load(sprintf('newOb0%d.mat',i))
%Zero the first point
 X = X-X(1)+1;
 Y = Y-Y(1);
% %Set the upper limit to 60KPa (threshold pressure)
%  for j = 1:length(Y)
%      if Y(j) > 60000
%          n = j;
%          break
%      end
%  end
%  X = X(1:n);
%  Y = Y(1:n);
 
%Fit for the collagen fiber stiffness 
KE = opt_K_Ob;
s = @(K,lam) stress_fit(DOR{i},LOR{i},K,lam,4);
cost = @(K)norm(Y - s(K,X));
opt_K = fminsearchbnd(cost, [1 1],[0 0],[]);
Obk{i} = opt_K;
cost(opt_K)
Xplot = linspace(1, 1.6);
Y_predict = s(opt_K,X)./1000;
%Plot stress curve fitting

figure(i)
plot(X, Y./1000, 'pr','MarkerSize',10)
hold on
plot(X, s(opt_K,X)./1000,'LineWidth',2)
xlabel('Tissue stretch')
ylabel('Stress(KPa)')
legend(['raw data'],['fitted curve'])
xlim([1,1.4])
ylim([-10,50])
set(gca,'fontsize',15)
grid
save(sprintf('finalOb0%d.mat',i))
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Plot Y tissue
for i = 1:1:4
load(sprintf('newY0%d.mat',i))
%Zero the first point
 X = X-X(1)+1;
 Y = Y-Y(1);
% %Set the upper limit to 60KPa (threshold pressure)
%  for j = 1:length(Y)
%      if Y(j) > 60000
%          n = j;
%          break
%      end
%  end
%  X = X(1:n);
%  Y = Y(1:n);
 
%Fit for the collagen fiber stiffness 
KE = opt_K_Y;
s = @(K,lam) stress_fit(DYR{i},LYR{i},K,lam,4);
cost = @(K)norm(Y - s(K,X));
opt_K = fminsearchbnd(cost, [1 1],[0 0],[]);
Yk{i} = opt_K;
cost(opt_K)
Xplot = linspace(1, 1.6);
Y_predict = s(opt_K,X)./1000;
%Plot stress curve fitting

figure(i+4)
plot(X, Y./1000, 'pr','MarkerSize',10)
hold on
plot(X, s(opt_K,X)./1000,'LineWidth',2)
xlabel('Tissue stretch')
ylabel('Stress(KPa)')
legend(['raw data'],['fitted curve'])
xlim([1,1.4])
ylim([-10,50])
set(gca,'fontsize',15)
grid
save(sprintf('finalY0%d.mat',i))
end
