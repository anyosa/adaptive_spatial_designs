%% Readme : description of files.
%
%% MAIN FILE TO RUN:
%
%% droploc_seq_designsSYNT.m
%
% RUNS SYNTHETIC EXAMPLE, calculation of optimal sequential design
% for drop locations. 
%
%
%% droploc_seq_designsLIZARD.m
%
% LIZARD ISLAND example, calculation of optimal sequential design
% for drop locations. 
%
% (Parameter estimation is done on Eastern end, then design on the rest
% end.)
%
%
%
%
%% sethyperpar.m 
%
% specify hyperparameters of the Gaussian process
% model and the regression mean and covariance model
% with covariates.
% For Lizard Island, they are as specified on Eastern flank. 
% For Synthetic case, they are set. 
%
%% findcovariates.m
%
% Specification of covariates. 
% Simulated in synthetic case, based on bathymetry
% inputs relevant for Lizard Island case.
%
%
%% findnodes.m 
%
% Extract grid nodes of design from 
% the input opstart and optime. 
% 
% Input also include currect velocity (vel)
% and domain of grid
% (min, max in east and north or some other coordinate system)
%
%
%% corals_comp_desNEW.m
% 
% 1st stage computation of design evaluation and
% setting the start value
%

%
%% corals_comp_desSEQnew.m
%
% file comparing a few static (non-sequential) designs
%
%  
%% EIBVcalc.m
%
% Calculation of the expected integrated Bernoulli variance
%
%
% 
%% BevalDesSEQnew.m (BevalDes.m)
%
% Approximate calculation of variance reductions in SEQUENTIAL 
% design of the linearized likelihood model.
% Using FFT to solve matrix systems effectively.
% 
%
%% Beval.m
%
% Update Gaussian model with data.
% Iterative linearization fitting Gaussian approximative posterior
% and variance reduction. 
% Using FFT to solve matrix systems effectively.
% 
%
%% approxGauss.m
%
% find a Gaussian approximation for the posterior (sequential) model.
% called in Beval.m

%
%
%% --- LIZARD ISLAND FILES--------
%
%% parameter_specification.m 
%
% Estimate we instead read locations 
% and the sample data at nodes.
%

%% --- SYNTHETIC FILES--------
%
%% getdata.m 
%
% Read the true data sampled along the path
% determined from nodes. 
% In real-world Lizard Island data we instead read locations 
% and the sample data at nodes.
%
%
%
%% compareSYNT.m
%
% compares replicate runs of three different strategies:
% - EIBV adaptive
% - spatial balance design (pre-scripted)
% - P nearest 0.5 adapative strategy.
% 
% Files are currently set up with a seed (rng(7)). 
% By modifications in the data generation, it can be read from file too,
% and similarly for input parameters in the replicates
% ('sethyperparSYNT.m')
% 
%
%
%% Replicate%d.csv ReplicateREF%d.csv
% 
%
% replicate parameter inputs and design results
% ordered by
% [sig Sigmabeta(1,1) nu Sigmabeta(2,1) scorevecADAPT(bb,1) scorevecADAPT(bb,2) scorevecADAPT(bb,3) scorevecSPAT(bb,1) scorevecSPAT(bb,2) scorevecSPAT(bb,3) scorevecP(bb,1) scorevecP(bb,2) scorevecP(bb,3)]);
% results are in files Replic%d.csv and ReplicREF%d.csv
% Spatial presence-absence datasets for the same replicates are in ysimL%d.txt.
