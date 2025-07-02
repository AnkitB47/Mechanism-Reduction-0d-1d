function plotDens(obj)
% plots the temperarture over time
switch nargin
    case 1
        plot(obj.times,obj.dens);
    case 2
        try
            plot(obj.times,obj.dens,optional);
        catch ME
            fprintf(['The argument of plotTemp may be not correct,',...
                ' check the result...\n The excemtions was:\n ',ME.message,'\n'])            
            plot(obj.times,obj.dens);
        end
    otherwise
end
xlabel('time [s]','FontSize',14);
ylabel('density [g/cm^3]','FontSize',14);
end