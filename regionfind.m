function [toe, f3] = regionfind(X,Y)
 
%     optimal_x = [u,r];
% Tfitn(i) = exponmodel(optimal_x,Stretchn(i));
f3 = fit(X',Y','exp1');
Stretchn = 1:0.0001:2;


for i = 1:length(Stretchn)
    
Tfitn(i) = f3(Stretchn(i));

if Tfitn(i) > 100000
    n = i;
    break
end
end
x2 = Stretchn(1:n);
y2 = Tfitn(1:n);
for i = 3:length(y2)
        P = polyfit(x2(1:i),y2(1:i),1);
        ypredict = P(1)*x2(1:i)+P(2);
        Rsq = 1 - sum((y2(1:i) - ypredict).^2)/sum((y2(1:i) - mean(y2(1:i))).^2);
        if Rsq <0.90
             disp(['toe ends at ',num2str(x2(i))])
           toe = x2(i);
            break
        end
    end
    
% ypredict = P(1)*x2+P(2);
% Toeslope = P(1);
%   for i = 1:length(x2)
%       if ypredict(i) < y2(i)-20;
%          disp(['toe ends at ',num2str(x2(i))])
%           toe = x2(i);
%          break
%       end
% %      
%   end

% figure(2)
% plot(x2,ypredict,x2,y2)
% hold on
%Find the transition and high stress region.

% for i = 2:length(y2)
%     n = length(y2);
%     P = polyfit(x2(n-i:n),y2(n-i:n),1);
%     ypredict = P(1)*x2(n-i:n)+P(2);
%     Rsq = 1 - sum((y2(n-i:n) - ypredict).^2)/sum((y2(n-i:n) - mean(y2(n-i:n))).^2);
%         if Rsq <0.98
%             disp(['high stress starts at ',num2str(x2(n-i))])
%             tran = x2(n-i);
%             break
%         end
% end
% 
% ypredict = P(1)*x2+P(2);  
% Highslope = P(1);

end