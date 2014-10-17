function test_matlab2()
clear all; close all;

%
% Figure de base :-)
%

figure;
[x, y, z] = egg(2,1,0.05,0);
surf(x,y,z); axis('off'); axis('equal');

%
% Pour faire un peu plus joli...
%

figure('Color',[1 1 1]); 
[x, y, z] = egg(2,1,0.05); 

h = surf(x,y,z,'FaceLighting','phong', 'LineStyle','none','FaceColor',[0.1 0.9 0.1]); 
rotate(h,[0 1 0],0,[0 0 0]); hold on; 

h = surf(1.5*x+4,1.5*y+3,1.5*z+0.5,'FaceLighting','phong', 'LineStyle','none','FaceColor',[0.1 0.1 0.9]); 

h = surf(2*x+3,2*y+5,2*z,'FaceLighting','phong', 'LineStyle','none','FaceColor',[0.9 0.1 0.1]); 
rotate(h,[0 1 0],-60,[0 0 0]); light('Position',[ 0.0 -0.75 0.5]);

light('Position',[-0.5 -0.75 0.5]); 
axis('off'); axis('equal');view([0 0]);
end







