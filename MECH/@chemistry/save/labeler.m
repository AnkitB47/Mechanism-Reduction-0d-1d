function labeler(z_label);
% macht nichts weiter, als
% die z-achse mit zlabel  und
% die x und y achse mit x_1 bzw x_2 zu beschriften

grid on
axis on
xlabel('x_{1}','FontSize',24);
ylabel('x_{2}','FontSize',24)
zlabel(z_label,'FontSize',24);
