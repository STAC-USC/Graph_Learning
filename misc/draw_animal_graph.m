
function draw_animal_graph(Laplacian,names)

nAnimals = size(Laplacian,1);
%coordinates
r=1; theta = linspace(0,2*pi,nAnimals+1);
coordinates = r*[cos(theta') sin(theta')];

sizeMat = size(Laplacian);
offDiag = ones(sizeMat) - eye(sizeMat); offDiag = offDiag == 1;
A = zeros(sizeMat); A(offDiag) = -Laplacian(offDiag);
% Threshold A
A(A<0) = 0;
%A(A<alpha*10) = 0;
load('ColorMatrix.mat');
figure;
colors = colormap(Fire/255);
wgPlot_animals(A,coordinates,'edgeColorMap',colors,'edgeWidth',2,'vertexColorMap',255*[1 1 1]);
%coordinates_text = coordinates; coordinates_text(:,1) = coordinates_text(:,1)-0.05; coordinates_text(:,2) = coordinates_text(:,2)+0.05;
coordinates_text = (r+0.015)*[cos(theta') sin(theta')];
fontsize = 12;
fontname = 'Helvetica';%'Times New Roman';
for t=1:length(names)
    x1=coordinates_text(t,1); y1=coordinates_text(t,2);
    if isequal(names{t},'Butterfly')
        x1 = x1-0.1;  y1 = y1-0.01;
    end
    if isequal(names{t},'Mouse')
        x1 = x1-0.1;
    end
    if isequal(names{t},'Squirrel')
        x1 = x1+0.04;
    end
    if isequal(names{t},'Alligator')
        x1 = x1+0.04;
    end
    if theta(t)>=0 && theta(t)<=pi/2
        text(x1,y1,names{t},'HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',fontsize,'FontName',fontname);
    elseif theta(t)>pi/2 && theta(t)<=pi
        text(x1,y1,names{t},'HorizontalAlignment','right','VerticalAlignment','bottom','FontSize',fontsize,'FontName',fontname);
    elseif theta(t)>pi && theta(t)<=3*pi/2
        text(x1,y1,names{t},'HorizontalAlignment','right','VerticalAlignment','top','FontSize',fontsize,'FontName',fontname);
    elseif theta(t)>3*pi/2 && theta(t)<=2*pi
        text(x1,y1,names{t},'HorizontalAlignment','left','VerticalAlignment','top','FontSize',fontsize,'FontName',fontname);
    end
end
set(gcf,'color','w');
axis square

end
