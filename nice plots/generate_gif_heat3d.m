function generate_gif_heat3d

% Load parameters
% Number of frames
numFrames=16;
% Time step between 2 frames
step = 0.3;
animated(1,1,1,numFrames) = 0;
ndim = 28;

% Main loop
for l=0:numFrames
file = strcat('temp_',num2str(l),'.dat');
fid=fopen(file,'r');        
raw = fread(fid,(ndim)^3,'float');
fclose(fid);
d = reshape(raw,[ndim, ndim, ndim]);
xslice = [ndim/2,ndim/2,ndim/2]; 
yslice = [ndim/2,ndim/2,ndim/2]; 
zslice = [ndim/2,ndim/2,ndim/2];

hFig = figure(1);
set(hFig, 'Position', [400 400 750 600]);
slice(d,xslice,yslice,zslice);
shading faceted;

view([-42,22]);
hc=colorbar;
set(hc,'position',[0.932 0.3 0.02 0.6]);
caxis([0 40])
xlabel('x domain');
ylabel('y domain');
zlabel('z domain');

pause(step);    
frame = getframe(figure(1));    
if l == 0
  [animated, cmap] = rgb2ind(frame.cdata, 256, 'nodither');
else
  disp(l);
  animated(:,:,1,l) = rgb2ind(frame.cdata, cmap, 'nodither');
end  
end
% Write final animated gif
imwrite(animated,cmap,'Heat_3D.gif','DelayTime',step,'LoopCount',inf); 

end