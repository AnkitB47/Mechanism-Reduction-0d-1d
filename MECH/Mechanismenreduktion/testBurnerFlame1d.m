clear classes
clear all
clc
clf
tic
fl = flame1d('/home/pruefert/data/inputs/gri30.cti','gri30_mix');

fl.compLeft = 'CH4:1,O2:2,N2:3.5';
 
fl.space = [0.0 0.02 0.04 0.06 0.08 0.1 ...
		0.15 0.2 0.4 0.49 0.5];
fl.initialTemperature = 373;
fl.initialPressure =  0.05*oneatm;
fl.mdotLeft =0.06;

fl.tInLeft = 373;

fl.flameType = 'burnerFlame';
   
fl = fl.solve;
 toc
%%
    hold on
    figure(1)
    fl.plotMassFractions([],'--')
    drawnow
    figure(2)
    hold on
    fl.plotTemp;
