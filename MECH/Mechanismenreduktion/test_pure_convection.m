gr = grid1D([0 1],0.001);
fem = lagrange11D; 

gr.makeBoundaryMatrix( gr.dirichletBC('1','0'),gr.dirichletBC('1','0'));

[K,M,F] = fem.assema(gr,0.001,0,1);
[Q,G,H,R] = fem.assemb(gr); 
B = fem.convection(gr,1);

y = (K+M+B+1e6*(H'*H)+Q)\(F+(1e6*H'*R)+G);

figure(1)
clf
plot(gr.p,y,'-');
 