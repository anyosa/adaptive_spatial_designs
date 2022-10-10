function [nodes,nodes1,nodes2]=findnodes(opstart,optime,vel,domain,ng);

%
% Extract grid nodes of design from 
% the input opstart and optime. 
% 
% Inputs also include currect velocity (vel)
% and domain of grid
% (min, max in east and north or some other coordinate system)
% as well as the grid size ng=[n1;n2].
%

opend=opstart+optime*vel;

% Find design nodes from operation design and domain and grid size
domNorm=domain(:,2)-domain(:,1);
node_start=ceil(ng.*((opstart-domain(:,1))./domNorm));
node_end=ceil(ng.*((opend-domain(:,1))./domNorm));
dd=(node_end-node_start);
Ddist=sqrt(sum(dd.^2));
rD=dd/Ddist;
if (dd(1,1)>0)
    nodes1INI=ceil(node_start(1,1)+[0:rD(1,1):dd(1,1)])';
else
    nodes1INI=ceil(node_start(1,1)+zeros(dd(2,1)+1,1));
end
nodes2INI=ceil(node_start(2,1)+[0:rD(2,1):dd(2,1)])';
% Set total nodes (summarized)
nodes=zeros(size(nodes1INI));

%count east, row by row in NS
nodesINI=(nodes2INI-1)*ng(1,1)+nodes1INI;
[nodes,i0]=unique(nodesINI);
nodes1=nodes1INI(i0,1);
nodes2=nodes2INI(i0,1);
