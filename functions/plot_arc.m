function plot_arc(alpha,beta,color,line_width)

u  = [cos(alpha/180*pi);sin(alpha/180*pi)];
v  = [cos(beta/180*pi);sin(beta/180*pi)];
x0 = -(u(2)-v(2))/(u(1)*v(2)-u(2)*v(1));
y0 =  (u(1)-v(1))/(u(1)*v(2)-u(2)*v(1));
r  = sqrt(x0^2 + y0^2 - 1);
thetaLim(1) = atan2(u(2)-y0,u(1)-x0);
thetaLim(2) = atan2(v(2)-y0,v(1)-x0);

if u(1) >= 0 && v(1) >= 0 
  % ensure the arc is within the unit disk
  theta = [linspace(max(thetaLim),pi,50),...
           linspace(-pi,min(thetaLim),50)].';
else
  theta = linspace(thetaLim(1),thetaLim(2)).';
end
line(r*cos(theta)+x0,...
r*sin(theta)+y0,...
'Color',color,...
'LineWidth', line_width,...
 'PickableParts','none');
end