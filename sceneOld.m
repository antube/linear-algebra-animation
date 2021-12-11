
clf  %This clears the figure, so remove this line if you want to preserve a plot you have already made
% This creates the 'background' axes
ha = axes('units','normalized', 'position',[0 0 1 1]);
% Move the background axes to the bottom
uistack(ha,'bottom');
% Load in a background image and display it using the correct colors
% The image used below, is just a Roadrunner scene I downloaded.
I=imread('NinjaHome.jpg');
hi = imagesc(I);
colormap gray;
% Turn the handlevisibility off so that we don't inadvertently plot into the axes again
% Also, make the axes invisible
set(ha,'handlevisibility','off', 'visible','off')
% Now we can use the figure, as required.
% For example, we can put a plot in an axes
%axes('position',[0.3,0.35,0.4,0.4])

filename = 'NinjaSword1.jpg';
ninjaColor =[1.0, 1, 1];
thresh = 215;
ninjasword1 = imread(filename);
ns1mtx = fJpeg2pointsConverter(ninjasword1, thresh);
[m,n]=size(ns1mtx);
fprintf("%s size (thresh=%i) , [%i,%i]",filename,thresh,m,n);
disp(m);  disp(n); 
ns1mtx = 0.075*ns1mtx;
ns1mtx = ns1mtx - [0;20]*ones(1,n);
ns1mtx_orig = ns1mtx;

%{ 
(sbh)
This "hb" thing is used for clipping the ninja to a specified axis 
rectangle smaller than the background.
If you don't wan't clipping, just make this box line up with the 
background.   To see the box, change axesVisible='on';  
%}
axesVisible = 'on'; 
axesXpos = 0.1;
axesYpos = 0.4;
axesXdim = 0.625;
axesYdim = 0.625; 


for ii=1:0.5:3
    hb = axes('units','normalized', 'position',[0.1 0.4 0.625 0.625]);
    %hb = axes('position',[axesXpos axesYpos axesXdim axesYdim]);
    h_rr = plot(hb,ns1mtx(1,:), ns1mtx(2,:),   '.', 'color', ninjaColor, 'MarkerSize', 1); 
    axis([0 70 0 70]) %This let me set the scale I wanted in the inserted axes
    set(gca,'color','none','handlevisibility',axesVisible,'visible',axesVisible)
    
    
    Shift = [2.5*ii;ii]*ones(1,n);
    ns1mtx = ns1mtx + Shift;
    
    
    pause(0.25)
    set(h_rr,'Visible','off')  % This line erases the image of the Road Runner and Wile E. Coyote
    axis([0 70 0 70]) % This let me set the scale I wanted in the inserted axes
    set( gca, 'color','none','handlevisibility','off','visible','off')
end

hb = axes('units','normalized', 'position',[0.1 0.4 0.625 0.625]);
h_rr = plot(hb, [1] , [1] ,   'yo',  'MarkerSize', 5); 
axis([0 70 0 70]) %This let me set the scale I wanted in the inserted axes
set(gca,'color','none','handlevisibility',axesVisible,'visible',axesVisible)
Shift = [2.5*ii;ii]*ones(1,n);
ns1mtx = ns1mtx + Shift;
pause()


