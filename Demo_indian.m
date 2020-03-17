%% =================================================================
% This script run LRTF-DFR-based method.
% More detail can be found in [1]
% [1] Yu-Bang Zheng, Ting-Zhu Huang*, Xi-Le Zhao*, Yong Chen, and Wei He.
%     Double Factor-Regularized Low-Rank Tensor Factorization for Mixed Noise Removal in Hyperspectral Image
% Please make sure your data is in range [0, 1].
% Created by Yu-Bang Zheng £¨zhengyubang@163.com£©
% March 14, 2020

%% =================================================================
clc;
clear;
close all;
addpath(genpath('lib'));
addpath(genpath('data'));
methodname  = { 'Noise', 'LRTF-DFR'};
Mnum = length(methodname);

%% Load initial data
% Please run 'generate_noisyHSI.m' in the file 'data' if there have no nosiy data
case_num=5;
load(strcat('indian_case',num2str(case_num),'.mat'))

%%
Nway = size(Ohsi);

%% evaluation indexes
Re_hsi  =  cell(Mnum,1);
psnr    =  zeros(Mnum,1);
ssim    =  zeros(Mnum,1);
sam     =  zeros(Mnum,1);
time    =  zeros(Mnum,1);
%%  corrupted image

i  = 1;
Re_hsi{i} = Nhsi;
[psnr(i), ssim(i), sam(i)] = HSIQA(Ohsi * 255, Re_hsi{i} * 255);
enList = 1;

%% Performing LRTF-FR
i = i+1;
%%%%%
opts=[];
opts.R       = 12;
opts.rho     = 0.1;
opts.tau     = 0.2;
opts.lambda  = 0.01;
opts.beta    = 15000;
opts.mu      = 0.04;
opts.max_it  = 50;
opts.Bmax_it = 10;
opts.tol     = 1e-4;
%.Xtrue   = Ohsi;
%%%%%
fprintf('\n');
disp(['performing ',methodname{i}, ' ... ']);
t0= tic;
[Re_hsi{i},A,B,S,Out] = Mixed_LRTF_DFR(Nhsi, opts);
%[Re_hsi{i},A,B,S,Out] = Mixed_LRTF_DFR_GPU(Nhsi, opts);
%Re_hsi{i}= gather(Re_hsi{i});
time(i) = toc(t0);
[psnr(i), ssim(i), sam(i)] = HSIQA(Ohsi * 255, Re_hsi{i} * 255);
fprintf('LRTF-DFR: PSNR = %5.4f  time = %5.2f\n',  psnr(i),time(i));
enList = [enList,i];


%% Show result
fprintf('\n');
fprintf('================== Result ==================\n');
fprintf(' %8.8s    %5.4s      %5.4s    %5.4s    \n', 'method','PSNR', 'SSIM', 'SAM');
for i = 1:length(enList)
    fprintf(' %8.8s    %5.4f    %5.4f    %5.4f    \n',...
        methodname{enList(i)},psnr(enList(i)), ssim(enList(i)), sam(enList(i)));
end
fprintf('================== Result ==================\n');

%%
figure,
showHSIResult(Re_hsi,Ohsi,0,1,methodname,enList,1,Nway(3))

