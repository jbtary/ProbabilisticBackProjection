% Plot an ellipse on top of the current figure. Implicitely uses same units
% as those provided in the inputs (no re-scaling)
% 
% In:
%     x0, y0: x and y coordinates of the ellipse center
%     a: radius in x direction
%     b: radius in y direction
% 
% Out:
%     None, plot the ellipse in current figure (black solid line)

function plot_ell(x0,y0,a,b)

t=-pi:0.01:pi;
x=x0+a*cos(t);
y=y0+b*sin(t);
hold on
plot(x,y,'k','LineWidth',1)
