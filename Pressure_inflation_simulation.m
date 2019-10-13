%fit for the collagenstiffness parameter
clear all
clc
%load the recruitment stretch: min mod max
load recruitment_stretch.mat
%load the wall thickness t and radius r.
r = [10.62 10.72 10.63 11.61 5.69 5.91 5.94 5.65];
t = [1.61 1.53 1.98 1.1 0.72 0.71 0.85 0.82];
lambda = 1:0.01:2;
for i = 1:4
    lambda = 1:0.0001:2;
    radius = r(i);
    thickness = t(i);
%Calculate the prefactor.
prefactor = 2.*thickness./radius*1./lambda.^3;
load(sprintf('finalOb0%d.mat',i))

 X = X-X(1)+1;
 Y = Y-Y(1);
%  for j = 1:length(Y)
%      if Y(j) > 60000
%          n = j;
%          break
%      end
%  end
%  X = X(1:n);
%  Y = Y(1:n);
%  KE = opt_K_adult;
% Calculate the pressure of all the constituients.
[sigma_collagen_lp sigma_collagen_dsm sigma_elastin sigma] = pressure_cal(DOR{i},LOR{i},opt_K,lambda,4);
%%---------%%
% plot pressure inflation
figure(i)
hold on
pressure_all = sigma.*prefactor;
pressure_lp = sigma_collagen_lp.*prefactor;
pressure_dsm = sigma_collagen_dsm.*prefactor;
pressure_e = sigma_elastin.*prefactor;
for j= 1:length(pressure_all)
    if pressure_all(j) >= 8924.05 %obstructed bladder
        n = j;
        break
    end
end
pressure_all = pressure_all(1:n);
pressure_lp = pressure_lp(1:n);
pressure_dsm = pressure_dsm(1:n);
pressure_e = pressure_e(1:n);
lambda = lambda(1:n);
plot(lambda,pressure_all,'LineWidth',3)
plot(lambda,pressure_lp,'LineWidth',3)
plot(lambda,pressure_dsm,'LineWidth',3)
plot(lambda,pressure_e,'LineWidth',3)
legend('Overall pressue', 'Lamina propria pressure', 'Detrusor pressure', 'Elastin pressure')
set(gca,'fontsize',15)
grid

end



for i = 1:4
    lambda = 1:0.0001:2;
    radius = r(i+4);
    thickness = t(i+4);
    prefactor = 2.*thickness./radius*1./lambda.^3;
load(sprintf('finalY0%d.mat',i))
% KE = opt_K_aged;

% Calculate the pressure of all the constituients.
[sigma_collagen_lp sigma_collagen_dsm sigma_elastin sigma] = pressure_cal(DYR{i},LYR{i},opt_K,lambda,4);
%%---------%%
% plot pressure inflation
figure(i+4)
hold on
pressure_all = sigma.*prefactor;
pressure_lp = sigma_collagen_lp.*prefactor;
pressure_dsm = sigma_collagen_dsm.*prefactor;
pressure_e = sigma_elastin.*prefactor;
for j= 1:length(pressure_all)
    if pressure_all(j) >= 8924.05 %Young control pressure 3922.66
        m = j;
        break
    end
end
pressure_all = pressure_all(1:m);
pressure_lp = pressure_lp(1:m);
pressure_dsm = pressure_dsm(1:m);
pressure_e = pressure_e(1:m);
lambda = lambda(1:m);
plot(lambda,pressure_all,'LineWidth',3)
plot(lambda,pressure_lp,'LineWidth',3)
plot(lambda,pressure_dsm,'LineWidth',3)
plot(lambda,pressure_e,'LineWidth',3)
ylim([0 10000])
xlabel('Stretch')
ylabel('Pressure (KPa)')
legend('Overall pressue', 'Lamina propria pressure', 'Detrusor pressure', 'Elastin pressure')
set(gca,'fontsize',15)
grid

end
