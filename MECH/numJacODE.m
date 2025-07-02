function J = numJacODE(fun,t,y,~)
    % Shortcut for numjac
   
    tr = 1e-6;
    thresh = tr(ones(size(y)));         

    [J,~] = numjac(fun,t,y,fun(t,y),thresh,[],0);  
end