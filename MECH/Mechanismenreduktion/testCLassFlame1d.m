clear classes
clear all
clc
clf


fl = flame1d('/home/uwe/Arbeit/mechanismen/gri30.cti','gas_1');


fl.flameType = 'counterFlowFlame';
fl.compLeft = 'CH4:1';
fl.compRight = 'O2:1';
fl.space = linspace(0,.2,50);
fl.initialTemperature = 390;
fl.initialPressure =  1*oneatm;
    fl.mdotLeft = 0.3;
    fl.mdotRight = 0.6;
    fl.tInLeft = 300;
    fl.tInRight = 300; 
    
    
fl.pureMixing;
% fl.plotMoleFractions([],':')

 
fl.solve;
 
%%
    hold on
    figure(1)

    fl.plotMoleFractions([],'--')
    drawnow
    figure(2)
    hold on
    fl.plotTemp;
