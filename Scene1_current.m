%% Create background image
clf %This clears the figure, so remove this line if you want to preserve a plot you have already made
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
ninjaColor =[1, 1, 1];
thresh = 215;
ninjasword1 = imread(filename);
ns1mtx = fJpeg2pointsConverter(ninjasword1, thresh);
[m,n]=size(ns1mtx);
fprintf("%s size (thresh=%i) , [%i,%i]",filename,thresh,m,n);
disp(m);  disp(n); 
ns1mtx = [ns1mtx;ones(1,15851)];
S = [0.025 0 0; 0 0.025 0; 0 0 1];  %This is my rescaling matrix to shrink the character to fit the background
ns1mtx = S*ns1mtx;
%ns1mtx = ns1mtx - [0;20]*ones(1,n);
ns1mtx_orig = ns1mtx;
%% Play background music throughout all scenes.
% Time the animation and match the length of animation with length of audio
[y,Fs] = audioread('ninja_music.wav');
player = audioplayer(y,Fs);
play(player)   % Start the player
stop(player)   % Stop whenever you like...
%% Notes from Stephen
%{ 
(sbh)
This "hb" thing is used for clipping the ninja to a specified axis 
rectangle smaller than the background.
If you don't wan't clipping, just make this box line up with the 
background.   To see the box, change boxVisible='on';  
%}
axesVisible = 'off'; 
axesXpos = 0;
axesYpos = 0;
axesXdim = 1;
axesYdim = 1; 

%% Run towards the center of the building (using shear)
ns1mtx = ShearHScene(ns1mtx,0.5);
for i=1:0.25:4.5
    hb = axes('units','normalized', 'position',[-0.2 .0625 1 1]);
    %hb = axes('position',[axesXpos axesYpos axesXdim axesYdim]);
    h_rr = plot(hb,ns1mtx(1,:), ns1mtx(2,:),   '.', 'color', ninjaColor, 'MarkerSize', 1); 
    axis([0 70 0 70]) %This let me set the scale I wanted in the inserted axes
    set(gca,'color','none','handlevisibility',axesVisible,'visible',axesVisible)
    
    
    Shift = [1 0 i; 0 1 0; 0 0 1];
    ns1mtx = Shift*ns1mtx;
    
    
    pause(0.4)
    set(h_rr,'Visible','off')  % This line erases the image of the Road Runner and Wile E. Coyote
    axis([0 70 0 70]) % This let me set the scale I wanted in the inserted axes
    set( gca, 'color','none','handlevisibility','off','visible','off')
end

%% Parkour
ns1mtx = ShearHScene(ns1mtx,-0.5);
ns1mtx = ReflHScene(ns1mtx);
for i=1:0.5:10
    hb = axes('units','normalized', 'position',[-0.2 .0625 1 1]);
    %hb = axes('position',[axesXpos axesYpos axesXdim axesYdim]);
    h_rr = plot(hb,ns1mtx(1,:), ns1mtx(2,:),   '.', 'color', ninjaColor, 'MarkerSize', 1); 
    axis([0 70 0 70]) %This let me set the scale I wanted in the inserted axes
    set(gca,'color','none','handlevisibility',axesVisible,'visible',axesVisible)
    
    
    Shift = [1 0 i; 0 1 0; 0 0 1];
    ns1mtx = Shift*ns1mtx;
    
    
    pause(0.25)
    set(h_rr,'Visible','off')  % This line erases the image of the Road Runner and Wile E. Coyote
    axis([0 70 0 70]) % This let me set the scale I wanted in the inserted axes
    set( gca, 'color','none','handlevisibility','off','visible','off')
end

%% Scale the building
for i=1:0.5:10
    
end

%% Functions



function fpiv = feetPivot(PP)
    % Get a pivot point at the feet of the character.
    uX = max(PP(1,:));
    lX = min(PP(1,:));
    %uY = max(PP(2,:));
    lY = min(PP(2,:));
    fpiv = [ mean([uX,lX])  ; lY ; 0];
end

function cent = centerPivot(PP)
    % Assume these points are moved into a scene frame.
    uX = max(PP(1,:));
    lX = min(PP(1,:));
    uY = max(PP(2,:));
    lY = min(PP(2,:));
    cent = [ mean([uX,lX])  ; mean([uY,lY]) ; 0];
end

function PPshh = ShearHScene(PP,k)
    [Mrows Ncols] = size(PP);
    if Mrows == 2,
        SH = [1 k ; 0 1];
    else ,
        SH = [1 k 0; 0 1 0; 0 0 1];
    end
    center = feetPivot(PP);
    PPz = ShiftScene(PP, -1.0*center(1,1), -1.0*center(2,1));
    if Mrows == 3,
        PPshh = (SH*PPz) + center;
    else ,
        PPshh = (SH*PPz) + center(1:2 , :); 
    end
end


function PPshsc = ShiftScene(PP,xD,yD)
    Shift = [1 0 xD; 0 1 yD; 0 0 1];
    [Mrows Ncols] = size(PP);
    if Mrows == 2,
        N1 = [PP(1,:) ; PP(2,:)  ;  ones(1,Ncols)];
    else ,
        N1 = PP;
    end
    shN1 = Shift*N1;
    if Mrows == 2,
        PPshsc = [shN1(1,:)  ; shN1(2,:)];
    else ,
        PPshsc = shN1;
    end
end


function PPout = fJpeg2pointsConverter(BB,THRESHOLD)
%This function will take in an N x M x 3 matrix that 
% has been imported into the workspace using the 
% imread('filename.jpg') command and stored in a matrix
% - it is called BB inside this converter.  
%USAGE: BBout = Jpeg2pointsConverter(BB,THRESHOLD)

    BB1=BB(:,:,1);
    [M, N]= size(BB1);
    BB1=double(BB1);
    BB2 = 255-BB1; %Invert so white is 0 instead of 255
    %Any point with high value is replaced by 1, and 
    %any point with a low value is replaced by 0
    BB3 = (BB2 > THRESHOLD);                     
    PP=zeros(2,M*N);
    cnt=0;
    for ii=1:M,
        for jj=1:N, 
            if (BB3(ii,jj)>0.5), 
                PP(:,cnt+1)=[jj;N-ii];
                cnt=cnt+1;
            end,
        end,
    end

    PPout = PP(:,1:cnt);
end

function PPrefl = ReflHScene(PP)
    [Mrows Ncols] = size(PP);
    if Mrows == 2,
        RE = [-1 0 ; 0 1];
    else ,
        RE = [-1 0 0; 0 1 0; 0 0 1];
    end
    center = feetPivot(PP);
    PPz = ShiftScene(PP, -1.0*center(1,1), -1.0*center(2,1));
    if Mrows == 3,
        PPrefl = (RE*PPz) + center;
    else ,
        PPrefl = (RE*PPz) + center(1:2 , :); 
    end
end

    

