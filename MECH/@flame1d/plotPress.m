function plotPress(obj)
% plots the temperarture over time
switch nargin
    case 1
        plot(obj.space,obj.pres);
    case 2
        try
            plot(obj.space,obj.pres,optional);
        catch ME
            fprintf(['The argument of plotTemp may be not correct,',...
                ' check the result...\n The excemtions was:\n ',ME.message,'\n'])            
            plot(obj.space,obj.pres);
        end
    otherwise
end
xlabel('t [s]','FontSize',20);
ylabel('p [Atm]','FontSize',20);
set(gca,'FontSize',20);
end