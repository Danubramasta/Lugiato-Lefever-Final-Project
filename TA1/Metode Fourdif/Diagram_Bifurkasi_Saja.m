clc; clear all; close all;

load('eta=1_a=1_gamma=1_L=50_N=200.mat')

no_fig = 1;
ls = 2;
stability(Ei,Norm,lamb,no_fig,ls)
xlim([0 3])
xlabel('Ei')
ylabel('|E|')
title(['Diagram Bifurkasi Menggunakan Pseudo-Spectral'])