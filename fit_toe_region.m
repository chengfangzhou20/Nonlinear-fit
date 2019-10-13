%aged 02 fit
clear all
clc
%Toe region end point of Ob
Ob_toe = [0.2 0.2 0.17 0.2];
%Toe region end point of Y
Y_toe = [0.2 0.2 0.2 0.2];
%Convert strain to stretch
Ob_toe = sqrt(Ob_toe.*2+1);
Y_toe = sqrt(Y_toe.*2+1);
load recruitment_stretch.mat
for i = 1:4
%load Ob data
load(sprintf('Ob0%d.mat',i))
X = sqrt(Ex{1}.*2+1);
Y = (Sx{1}+Sy{1})./2;
%Clean the noise at the beginning
if i == 4
    X = X(30:end);
    Y = Y(30:end);
end
if i == 1
    X = X(2:end);
    Y = Y(2:end);
end
%Plot raw data
figure(i)
hold on
%  X = X-X(1)+1;
%  Y = Y-Y(1);
 delta = mean(Y(1:3));
 Y = Y-delta;
 
plot(X,Y,'p');

% Ob_toe(i) = DOR{i}(1);
 %Load toe region: Aged02 to Aged05
 for j = 1:length(Y)
     if X(j) > Ob_toe(i)
         n = j;
         break
     end
 end
 XX = X(1:n);
 YY = Y(1:n);
 
f = fit(X',Y','exp1');

s = @(K,lam) neo_fit(K,lam);
cost = @(K)norm(YY - s(K,XX));
opt_K_Ob = fminsearchbnd(cost,[1],[0],[]);
K_elastin_Ob{i} = opt_K_Ob;
stress = neo_fit(opt_K_Ob,XX);
plot(XX,stress,'-');
hold on
plot(X,f(X));
save(sprintf('newOb0%d.mat',i))
end


for i = 1:4
%load young bladder data
load(sprintf('Y0%d.mat',i))
X = sqrt(Ex{1}.*2+1);
Y = (Sx{1}+Sy{1})./2;
%Clean the noise at the beginning
if i == 3
    X = X(6:end);
    Y = Y(6:end);
end
if i == 1
    X = X(3:end);
    Y = Y(3:end);
end
% %Plot raw data
figure(i+4)
hold on
%  X = X-X(1)+1;
%  Y = Y-Y(1);
 delta = mean(Y(1:3));
 Y = Y-delta;

plot(X,Y,'p');

% Y_toe(i) = DYR{i}(1);
 %Load toe region: Aged02 to Aged05
 for j = 1:length(Y)
     if X(j) > Y_toe(i)
         n = j;
         break
     end
 end
 XX = X(1:n);
 YY = Y(1:n);

f = fit(X',Y','exp1');

s = @(K,lam) neo_fit(K,lam);
cost = @(K)norm(YY - s(K,XX));
opt_K_Y = fminsearchbnd(cost,[1],[0],[]);
K_elastin_Y{i} = opt_K_Y;
stress = neo_fit(opt_K_Y,XX);
plot(XX,stress,'-');
hold on
plot(X,f(X));

save(sprintf('newY0%d.mat',i))
end
