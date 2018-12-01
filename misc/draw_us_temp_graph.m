
function draw_us_temp_graph(Laplacian,center_vector)

sizeMat = size(Laplacian); offDiag = ones(sizeMat) - eye(sizeMat); offDiag = offDiag == 1;
A = zeros(sizeMat); A(offDiag) = -Laplacian(offDiag);
A(A<sqrt(eps)) = 0;
load('ColorMatrix.mat');
figure; colors = colormap(Fire/255);
[lat,long] = borders('continental us'); %borders('countries','nomap'); %load coast;
hold on;
for i=1:length(lat)
    plot(long{i},lat{i},'k');
end
wgPlot_us_temp(A,center_vector,'edgeColorMap',colors,'edgeWidth',2,'vertexColorMap',255*[1 1 1]);
set(gcf,'color','w');
axis tight

end
