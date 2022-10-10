%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pick elements of cm matrix from nodes and jj index
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tm=pickcovs(cm,nodes1,nodes2,jj);

tm=zeros(size(nodes1));
[n1,n2]=size(cm);
jj2=ceil(jj/n1); % Index of this jj location, second direction
jj1=jj-(jj2-1)*n1; % Index of this jj location, first direction
for k1=1:size(nodes1,1),
  tm(k1,1)=cm(1+mod2(nodes1(k1,1),jj1,n1),1+mod2(nodes2(k1,1),jj2,n2));
end;

