clc
close all
clear all
addpath('..')

%% Options
do_codegen = 1; % print files
pick_dir = 1;
dirs = {'X','Y'};
fname_RM1 = '../ORMS/GoldenBPMResp_I04.mat';
fname_RM2 = '../ORMS/ORM_11.9.2022/GoldenBPMResp_I04.mat';

load(fname_RM1);
RM1 = Rmat(1).Data(:,:);%  * 1e6; % RM(1) eq to RM(1,1)
load(fname_RM2);
RM2 = Rmat(1).Data(:,:);%  * 1e6; % RM(4) eq to RM(2,2)







