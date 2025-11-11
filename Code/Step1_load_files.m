clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The following variates are needed to be loaded before performing each analysis step (step2-step6).

% set working paths
Pathnew=['D:\Documents\OneDrive - bnu.edu.cn\My_publications\0-8yearspaper\Forsubmission2025\'];% basic dir
statdir=[Pathnew,'Forgithub\Code'];
savedir=[Pathnew,'Forgithub'];
Pathnewplot=[Pathnew,'Forgithub'];

% load brain surface for brain plotings
BNVsurfacetemplate=[Pathnew,'Forgithub\Data','\BrainMesh_ICBM152_smoothed.nv'];

% load coordinates and name for brain nodes
Pathnewxsl=[Pathnew,'Forgithub\Data\Atlas\500atlas.xlsx'];

% load individual network and covariates, using subset1 as an example.
load([Pathnew,'Forgithub\Data\Network\Subset1\Matrixallpar500FN_subset1.mat']);
load([Pathnew,'Forgithub\Data\Cov\Subset1\Age_gender_Brainsize_subset1.mat']) 
load([Pathnew,'Forgithub\Data\Atlas\Her3name.mat']) 