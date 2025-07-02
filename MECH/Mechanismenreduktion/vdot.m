function dv = vdot(t,vol,gas)
dv = 1e6*(pressure(gas) - 101325.0);
end