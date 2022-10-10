function plotksdiff(diff12,diff13,diff23,gamma);

%
%
%

KK=size(diff12,2),
figure(16);
clf;
%  for ii=1:KK,
%      subplot(2,2,ii), hold;
%  end
for ii=1:KK,
    subplot(1,3,ii), plot([0 gamma],[0 0],'k'); hold;
end

for ii=1:KK,
    [f12,x12]=ksdensity(diff12(:,ii));
    [f13,x13]=ksdensity(diff13(:,ii));
    [f23,x23]=ksdensity(diff23(:,ii));
    subplot(1,3,ii), plot(f12,x12,'r');
    subplot(1,3,ii), plot(f13,x13,'k--');
    subplot(1,3,ii), plot(f23,x23,'b.');
   
end

subplot(1,3,1), title('EIBV');
subplot(1,3,2), title('Integrated mis-classification prob');
subplot(1,3,3), title('Negative log score');

