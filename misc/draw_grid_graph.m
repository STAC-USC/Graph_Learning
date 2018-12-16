%% Draw a grid graph with link weights having different colors 
% Inputs: A: Adjacency matrix
%         BS: block size
function draw_grid_graph(A,BS)
load ColorMatrix;
for i=1:BS
    for j=1:BS 
        pixel = i+(j-1)*BS;
        graph_points(pixel,2)= -i; % height
        graph_points(pixel,1)= j; % width
    end
end

colors = colormap(Fire/255);
wgPlot_grid(A,graph_points,'edgeColorMap',colors,'edgeWidth',2,'vertexColorMap',255*[1 1 1]);
set(gcf,'color','w');
% wgPlot(A,graph_points,'edgeColorMap',jet);
% set(gcf,'color','w');
