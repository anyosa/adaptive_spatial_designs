function [delta,min_east,max_east,min_north,max_north,domain,n1,n2,n,nv,vx,vy,vel]=setparLIZARD;


%
% Specification on 
% grid spacing and domain (N-S and E-W, or other coordinates)
% [ this defines grid size]
% 
% Specification also includes velocity of current field (assumed fixed 
% and constant across the grid)
%

delta=1; % m resolution
% north coordinates
% lat / long boundary box
min_north=8374282.6;
min_east=332780;
max_north=8374389.6;
max_east=332938;
domain=[min_east max_east;min_north max_north];

n1=ceil(max_east-min_east)/delta;
n2=ceil(max_north-min_north)/delta;
nv=[n1;n2];
%%
% Total size
n=n1*n2;

%% Set current
%vx=-0.73; vy=-0.68; %tot 1 m/s
vx=0.86; vy=-0.43;
vel=[vx;vy];

