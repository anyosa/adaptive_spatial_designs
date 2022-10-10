function y=getdata(xx,nodes);

%
% Extraction of data from true model
% 
% In a real-world setting this must just read the 
% actual sampling nodes and the registered data.
%

eta=xx(nodes,1); 
y=zeros(size(nodes,1),1);
for i=1:size(nodes,1),
    y(i,1)=(rand<(1/(1+exp(-eta(i,1)))));
end
