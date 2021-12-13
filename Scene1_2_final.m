%{
SCENE 1 - Jake Kaplan
%}

%function [out_flag, ns1mtx, characterCenter1] = Scene1_2_final(~)
%% Play background music throughout all scenes.
[y,Fs] = audioread('ninja_music.wav');
player = audioplayer(y,Fs);
play(player)   % Start the music

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
ninjaColor =[0, 0, 1];
thresh = 219;
ninjasword1 = imread(filename);
ns1mtx = fJpeg2pointsConverter(ninjasword1, thresh);
[m,n]=size(ns1mtx);
fprintf("%s size (thresh=%i) , [%i,%i]",filename,thresh,m,n);
disp(m);  disp(n); 
ns1mtx = [ns1mtx;ones(1,n)]; %Make the matrix 3x3 by adding a row of 1s
S = [0.02 0 0; 0 0.02 0; 0 0 1];  %This is my rescaling matrix to shrink the character to fit the background
ns1mtx = S*ns1mtx;
ns1mtx_orig = ns1mtx;
%gif('Scene1_2_final.m.%gif')

ninjaStarColor =[1, 1, 1];
% import the throwing star sprite
throwingStar = fJpeg2pointsConverter(imread("throwing-star.jpg"), thresh);
% get the size and convert the matrix to a set of homogenous coordinates
[m,n]=size(throwingStar);
throwingStar = [throwingStar;ones(1,n)];
% rescale the throwing star to the character
throwingStar = S*throwingStar;



axesVisible = 'off'; 
axesXpos = 0;
axesYpos = 0;
axesXdim = 1.2;
axesYdim = 1; 

%% Run towards the edge of the building (using shear)
ns1mtx = ShearHScene(ns1mtx,0.5);
hb = axes('units','normalized', 'position',[-0.2 .0625 axesXdim 1]);
r = 1/5;
numItr = 17.5;
for i=1:0.5:numItr
    %hb = axes('position',[axesXpos axesYpos axesXdim axesYdim]);
    h_rr = plot(hb,ns1mtx(1,:), ns1mtx(2,:),   '.', 'color', ninjaColor, 'MarkerSize', 1); 
    axis([0 70 0 70]) %This let me set the scale I wanted in the inserted axes
    set(gca,'color','none','handlevisibility',axesVisible,'visible',axesVisible)
    
    Shift = [1 0 1; 0 1 0; 0 0 1];
    ns1mtx = Shift*ns1mtx;
    ns1mtx = RotationScene(ns1mtx,r);
    r = -1*r;
    
    
    pause(0.1)
    set(h_rr,'Visible','off')  % This line erases the image of the Road Runner and Wile E. Coyote
    axis([0 70 0 70]) % This let me set the scale I wanted in the inserted axes
    set( gca, 'color','none','handlevisibility','off','visible','off')
end
ns1mtx = RotationScene(ns1mtx,r);

%% Reflect character and jump to left
ns1mtx = ShearHScene(ns1mtx,-0.5);
ns1mtx = ReflHScene(ns1mtx);
hb = axes('units','normalized', 'position',[-0.2 .0625 axesXdim 1]);
numItr = 12;
for i=1:numItr
    %hb = axes('position',[axesXpos axesYpos axesXdim axesYdim]);
    h_rr = plot(hb,ns1mtx(1,:), ns1mtx(2,:),   '.', 'color', ninjaColor, 'MarkerSize', 1); 
    axis([0 70 0 70]) %This let me set the scale I wanted in the inserted axes
    set(gca,'color','none','handlevisibility',axesVisible,'visible',axesVisible)
    Shift = [1 0 -(6/numItr); 0 1 (6/numItr); 0 0 1];
    ns1mtx = Shift*ns1mtx;
    
    
    pause(0.001)
    set(h_rr,'Visible','off')  % This line erases the image of the Road Runner and Wile E. Coyote
    axis([0 70 0 70]) % This let me set the scale I wanted in the inserted axes
    set( gca, 'color','none','handlevisibility','off','visible','off')
end

%% Character scales the building
hb = axes('units','normalized', 'position',[-0.2 .0625 axesXdim 1]);
r = 1/9;
for i=1:9
    %hb = axes('position',[axesXpos axesYpos axesXdim axesYdim]);
    h_rr = plot(hb,ns1mtx(1,:), ns1mtx(2,:),   '.', 'color', ninjaColor, 'MarkerSize', 1); 
    axis([0 70 0 70]) %This let me set the scale I wanted in the inserted axes
    set(gca,'color','none','handlevisibility',axesVisible,'visible',axesVisible) 
    
    Shift = [1 0 0; 0 1 1; 0 0 1];
    ns1mtx = Shift*ns1mtx;
    ns1mtx = RotationScene(ns1mtx,r);
    r = -1*r;
    
    
    pause(0.2)
    set(h_rr,'Visible','off')  % This line erases the image of the Road Runner and Wile E. Coyote
    axis([0 70 0 70]) % This let me set the scale I wanted in the inserted axes
    set( gca, 'color','none','handlevisibility','off','visible','off')
end 
ns1mtx = RotationScene(ns1mtx,r);

%% Reflect character and jump to right (to reach roof)
ns1mtx = ReflHScene(ns1mtx);
hb = axes('units','normalized', 'position',[-0.2 .0625 axesXdim 1]);
for i=1:numItr
    %hb = axes('position',[axesXpos axesYpos axesXdim axesYdim]);
    h_rr = plot(hb,ns1mtx(1,:), ns1mtx(2,:),   '.', 'color', ninjaColor, 'MarkerSize', 1); 
    axis([0 70 0 70]) %This let me set the scale I wanted in the inserted axes
    set(gca,'color','none','handlevisibility',axesVisible,'visible',axesVisible)
    
    Shift = [1 0 (5/numItr); 0 1 (5/numItr); 0 0 1];
    ns1mtx = Shift*ns1mtx;
    
    
    pause(0.001);
    set(h_rr,'Visible','off');  % This line erases the image of the Road Runner and Wile E. Coyote
    axis([0 70 0 70]) ;% This let me set the scale I wanted in the inserted axes
    set( gca, 'color','none','handlevisibility','off','visible','off');;
end  

characterCenter1 = centerPivot(ns1mtx);

x_final = characterCenter1(1,1);
y_final = characterCenter1(2,1);
fprintf("x_final = %f", x_final);
fprintf("y_final = %f", y_final);

%%
%{
SCENE 2 - Stephen Horn
%}

%ns1mtx = teleportTo(ns1mtx,35,25);

%% Lands on to roof
for i=1:5
    hb = axes('units','normalized', 'position',[-0.2 .0625 1.2 1]);
    h_rr = plot(hb,ns1mtx(1,:), ns1mtx(2,:),   '.', 'color', ninjaColor, 'MarkerSize', 1); 
    axis([0 70 0 70]) ;
    set(gca,'color','none','handlevisibility',axesVisible,'visible',axesVisible)
    
    nS = [1 0 0.5 ; 0 1 -0.1; 0 0 1 ];
    ns1mtx = nS*ns1mtx;
    
    
    pause(0.05);
    set(h_rr,'Visible','off');  
    axis([0 70 0 70]) ;
    set( gca, 'color','none','handlevisibility','off','visible','off');
end

%% sneaks... 
ns1mtx = squatScene(ns1mtx,1.8,0.6);
r=-1;
for i=1:28
    hb = axes('units','normalized', 'position',[-0.2 .0625 1.2 1]);
    h_rr = plot(hb,ns1mtx(1,:), ns1mtx(2,:),   '.', 'color', ninjaColor, 'MarkerSize', 1); 
    axis([0 70 0 70]) ;
    set(gca,'color','none','handlevisibility',axesVisible,'visible',axesVisible)
    
    nS = [1 0 0.5 ; 0 1 0; 0 0 1 ];
    ns1mtx = nS*ns1mtx;
    ns1mtx = squatScene(ns1mtx, 1.0 + (0.2*r) , 1.0);
    r=-1*r;
    
    
    pause(0.05);
    set(h_rr,'Visible','off');  
    axis([0 70 0 70]) ;
    set( gca, 'color','none','handlevisibility','off','visible','off');
end

%% Character stands up from sneak position
algn = alignWith(ns1mtx, ns1mtx_orig); 
ns1mtx = algn;

for i=1:4
    hb = axes('units','normalized', 'position',[-0.2 .0625 1.2 1]);
    h_rr = plot(hb,ns1mtx(1,:), ns1mtx(2,:),   '.', 'color', ninjaColor, 'MarkerSize', 1); 
    axis([0 70 0 70]) ;
    set(gca,'color','none','handlevisibility',axesVisible,'visible',axesVisible)
    
    % sv + c
    nS = [1 0 0.5 ; 0 1 0; 0 0 1 ];
    ns1mtx = nS*ns1mtx;
    
    
    pause(0.05);
    set(h_rr,'Visible','off');  
    axis([0 70 0 70]) ;
    set( gca, 'color','none','handlevisibility','off','visible','off');
end

nt4mtx = loadNinjaTool4('NinjaTool4.jpg');
Z = (-1)*centerPivot(nt4mtx);
nt4mtx = ShiftScene(nt4mtx, Z(1),Z(2));
nt4mtx = [-1 0 0; 0 -1 0; 0 0 1]*nt4mtx;
algn = alignWith(ns1mtx , nt4mtx);
nt4mtx = algn;

%% Frontflip
v=1;
for i=1:19
    hb = axes('units','normalized', 'position',[-0.2 .0625 1.2 1]);
    h_rr = plot(hb,nt4mtx(1,:), nt4mtx(2,:),   '.', 'color', ninjaColor, 'MarkerSize', 1); 
    axis([0 70 0 70]) ;
    set(gca,'color','none','handlevisibility',axesVisible,'visible',axesVisible)
    
    % sv + c
    nS = [1 0 0.4 ; 0 1 (-0.28)*v+3; 0 0 1 ];
    nt4mtx = nS*nt4mtx;
    nt4mtx = RotationScene(nt4mtx, -0.66 );
    v=v+1;
    
    
    pause(0.05);
    set(h_rr,'Visible','off');  
    axis([0 70 0 70]) ;
    set( gca, 'color','none','handlevisibility','off','visible','off');
end

%% Lands and walks on roof
algn = alignWith(nt4mtx , ns1mtx);
ns1mtx = algn;
for i=1:6
    hb = axes('units','normalized', 'position',[-0.2 .0625 1.2 1]);
    h_rr = plot(hb,ns1mtx(1,:), ns1mtx(2,:),   '.', 'color', ninjaColor, 'MarkerSize', 1); 
    axis([0 70 0 70]) ;
    set(gca,'color','none','handlevisibility',axesVisible,'visible',axesVisible)
    
    nS = [1 0 0.5; 0 1 0; 0 0 1];
    ns1mtx = nS*ns1mtx;
    
    
    pause(0.05);
    set(h_rr,'Visible','off');  
    axis([0 70 0 70]) ;
    set( gca, 'color','none','handlevisibility','off','visible','off');
end

%% Jumps off roof to the edge of the screen
for i=1:5
    hb = axes('units','normalized', 'position',[-0.2 .0625 1.2 1]);
    h_rr = plot(hb,ns1mtx(1,:), ns1mtx(2,:),   '.', 'color', ninjaColor, 'MarkerSize', 1); 
    axis([0 70 0 70]) ;
    set(gca,'color','none','handlevisibility',axesVisible,'visible',axesVisible)
    
    nS = [1 0 1.5 ; 0 1 1; 0 0 1 ];
    ns1mtx = nS*ns1mtx;
    
    
    pause(0.05);
    set(h_rr,'Visible','off');  
    axis([0 70 0 70]) ;
    set( gca, 'color','none','handlevisibility','off','visible','off');
end

characterCenter2 = centerPivot(ns1mtx);

x_final = characterCenter2(1,1);
y_final = characterCenter2(2,1);
fprintf("x_final = %f", x_final);
fprintf("y_final = %f", y_final);





%{
SCENE 3 - Andrew Brown   
%}

% Call scene three function
failureFlag = false;
[failureFlag, ns1mtx, characterCenter, throwingStar1, throwingStar2] = third_scene(ns1mtx, [x_final, y_final], throwingStar, throwingStar, ninjaColor, ninjaStarColor, axesVisible);
x_final = characterCenter(1,:);
y_final = characterCenter(2,:);





%%
%{
SCENE 4 - Giovanni 
%}
% =======
CA=imread('NinjaSword1.jpg');
CAout=fJpeg2pointsConverter(CA,219);
CD=imread('NinjaSword3.jpg');
CDout=fJpeg2pointsConverter(CD,219);
CF=imread('NinjaTool2.jpg');
CFout=fJpeg2pointsConverter(CF,219);
CG=imread('NinjaTool3.jpg');
CGout=fJpeg2pointsConverter(CG,219);
CI=imread('SmokeBomb.jpg');
CIout=fJpeg2pointsConverter(CI,219);
CB=imread('ninjalogo1.jpg');
CBout=fJpeg2pointsConverter(CB,219);
A=CAout;
[m,n1]=size(CAout);
disp(m);  disp(n1); 
CAout = [CAout;ones(1,n1)]; 
S = [0.02 0 0; 0 0.02 0; 0 0 1]; 
CAout = S*CAout;
B=CBout;
[m,n2]=size(CBout);
disp(m);  disp(n2); 
CBout = [CBout;ones(1,n2)]; 
S = [0.02 0 0; 0 0.02 0; 0 0 1]; 
CBout = S*CBout;
D=CDout;
[m,n3]=size(CDout);
disp(m);  disp(n3); 
CDout = [CDout;ones(1,n3)]; 
S = [0.02 0 0; 0 0.02 0; 0 0 1]; 
CDout = S*CDout;
F=CFout;
[m,n4]=size(CFout);
disp(m);  disp(n4); 
CFout = [CFout;ones(1,n4)]; 
S = [0.02 0 0; 0 0.02 0; 0 0 1]; 
CFout = S*CFout;
G=CGout;
[m,n5]=size(CGout);
disp(m);  disp(n5); 
CGout = [CGout;ones(1,n5)]; 
S = [0.02 0 0; 0 0.02 0; 0 0 1]; 
CGout = S*CGout;
I=CIout;
[m,n6]=size(CIout);
disp(m);  disp(n6); 
CIout = [CIout;ones(1,n6)]; 
S = [0.02 0 0; 0 0.02 0; 0 0 1];
CAout(3,12878);
Z=zeros(3,279);
CAout_New= [CAout,Z];
CFout(3,12172);
Z1=zeros(3,985);
CFout_New=[CFout,Z1];
CGout(3,12078);
Z2=zeros(3,94);
CGout_New=[CGout,Z2];
Z3=zeros(3,1831);
CGout_New2=[CGout,Z3];
CIout(3,13909);
Z4=zeros(3,12245);
CIout_New=[CIout,Z4];
disp(S)
CAout_New=S*CAout_New;
CDout=S*CDout;
CFout_New=S*CFout_New;
CFout=S*CFout;
CGout_New=S*CGout_New;
CGout_New2=S*CGout_New2;
CIout=S*CIout;
CIout_New=S*CIout_New;
CBout=S*CBout;
for k=0:1/8:1
	B = (1-k)*CAout_New + k*CDout;
    hb = axes('units','normalized', 'position',[0.2 0 0.5 0.5]);
    h_rr = plot(hb,B(1,:),B(2,:),'.');
    set(gca,'color','none','handlevisibility',axesVisible,'visible',axesVisible)
    pause(0.25)
    set(h_rr,'Visible','off')   
end
for k=0:1/8:1
	B = (1-k)*CDout + k*CFout_New;
	hb = axes('units','normalized', 'position',[0.2 0 0.5 0.5]);
    h_rr = plot(hb,B(1,:),B(2,:),'.');
    set(gca,'color','none','handlevisibility',axesVisible,'visible',axesVisible)
    pause(0.25)
    set(h_rr,'Visible','off')   
end
for k=0:1/8:1
	B = (1-k)*CFout + k*CGout_New;
	hb = axes('units','normalized', 'position',[0.2 0 0.5 0.5]);
    h_rr = plot(hb,B(1,:),B(2,:),'.');
    set(gca,'color','none','handlevisibility',axesVisible,'visible',axesVisible)
    pause(0.25)
    set(h_rr,'Visible','off')   
end
for k=0:1/8:1
	B = (1-k)*CGout_New2 + k*CIout;
	hb = axes('units','normalized', 'position',[0.2 0 0.5 0.5]);
    h_rr = plot(hb,B(1,:),B(2,:),'.');
    set(gca,'color','none','handlevisibility',axesVisible,'visible',axesVisible)
    pause(0.25)
    set(h_rr,'Visible','off')    
end
for k=0:1/8:1
    B = (1-k)*CIout_New + k*CBout;
	hb = axes('units','normalized', 'position',[0.2 0 0.5 0.5]);
    h_rr = plot(hb,B(1,:),B(2,:),'.');
    set(gca,'color','none','handlevisibility',axesVisible,'visible',axesVisible)
    pause(0.25)
    set(h_rr,'Visible','off')   
end
stop(player)   % Stop the music after the animation is complete.
disp('script completed');
%{
----------------------------------------------------------
Functions below
%}

function PPt = teleportTo(PP,tx,ty)
    nc = centerPivot(PP); 
    nP = [1 0 -1*nc(1) ; 0 1 -1*nc(2); 0 0 1 ];
    zPP = nP*PP;
    nS = [1 0 tx ; 0 1 ty; 0 0 1 ];
    PPt = nS*zPP;
end


function PPal = alignWith(PPprevmtx , newmtx )
    [Mrows Ncols] = size(PPprevmtx);
    center = feetPivot(newmtx);
    newzzero = ShiftScene(newmtx, -1.0*center(1,1), -1.0*center(2,1));
    prevc = feetPivot(PPprevmtx);
    if Mrows == 3,PPal = newzzero + prevc;
    else, PPal = newzzero + prevc(1:2 , :); 
    end
end


function PPq = squatScene(PP, xq, yq )
    [Mrows Ncols] = size(PP);
    if Mrows == 2, SH = [xq 0 ; 0 yq];
    else , SH = [xq 0 0; 0 yq 0; 0 0 1];
    end
    center = feetPivot(PP);
    PPz = ShiftScene(PP, -1.0*center(1,1), -1.0*center(2,1));
    if Mrows == 3,PPq = (SH*PPz) + center;
    else, PPq = (SH*PPz) + center(1:2 , :); 
    end
end


function nt4mtx = loadNinjaTool4(filename)
    thresh = 219;
    ninjatool4 = imread(filename);
    nt4mtx = fJpeg2pointsConverter(ninjatool4, thresh);
    [m,n]=size(nt4mtx);
    fprintf("%s size (thresh=%i) , [%i,%i]",filename,thresh,m,n);
    disp(m);  disp(n); 
    nt4mtx = [nt4mtx;ones(1,n)];
    %This is my rescaling matrix to shrink the character to fit the background
    S = [0.025 0 0; 0 0.025 0; 0 0 1];  
    nt4mtx = S*nt4mtx;
end


function fpiv = feetPivot(PP)
    % Get a pivot point at the feet of the character.
    uX = max(PP(1,:));
    lX = min(PP(1,:));
    %uY = max(PP(2,:));
    lY = min(PP(2,:));
    fpiv = [ mean([uX,lX])  ; lY ; 0];
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


function cent = centerPivot(PP)
    % Assume these points are moved into a scene frame.
    uX = max(PP(1,:));
    lX = min(PP(1,:));
    uY = max(PP(2,:));
    lY = min(PP(2,:));
    cent = [ mean([uX,lX])  ; mean([uY,lY]) ; 0];
end


function PPrs = RotationScene(PP,radAngle)
    th=radAngle;
    [Mrows Ncols] = size(PP);
    if Mrows == 2 ,
        R = [cos(th) -sin(th); sin(th) cos(th)];
    else ,
        R = [cos(th) -sin(th) 0; sin(th) cos(th) 0 ; 0 0 1];
    end
    center = centerPivot(PP);
    PPz = ShiftScene(PP, -1.0*center(1,1), -1.0*center(2,1));
    Prot = R*PPz;
    PPrs = Prot + center;
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


function PPout = fJpeg2pointsConverter(BB,THRESHOLD)
    BB1=BB(:,:,1);
    [M, N]= size(BB1);
    BB1=double(BB1);
    BB2 = 255-BB1; 
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


function  [failureFlag, character, characterCenter, throwingStar1, throwingStar2] = third_scene(character, characterCenter, throwingStar1, throwingStar2, ninjaColor, ninjaStarColor, axesVisible)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Setup the nessecary matrices
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % have the character fall into the scene
    fallTransformation = [1 0 0.2; 0 1 -1; 0 0 1];

    % landing matrices.
    compressionTransformation = [1 0 0; 0 0.90 0; 0 0 1];
    decompressionTransformation = inv(compressionTransformation);


    % This transformation matrix  is used to move the ninja stars across the scene
    throwingTransformation1 = [1 0 1.05; 0 1 -0.25; 0 0 1];
    throwingTransformation2 = [1 0 1; 0 1  0  ; 0 0 1];

    throwingRotationTransformation = [1 0 0; 0 1 0; 0 0 1];


    % running transformation matrix
    runningTransformation = [1 0 1; 0 1 0; 0 0 1];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Perform the Animation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    

    character = teleportTo(character, 10, 23);


    % character falls into scene from the previous building jump
    % character lands
    % character throws ninja star at target
    % runs past target to the middle of the scene
    % throws ninja star at second target

    
    % have the character fall into scene
    for j = 1:20
        % setup the plot for the animation frame
        hb = axes('units','normalized', 'position',[-0.2 0.0625 1.2 1]);
        h_rr = plot(hb,character(1,:), character(2,:),   '.', 'color', ninjaColor, 'MarkerSize', 1);
        axis([0 70 0 70]) %This let me set the scale I wanted in the inserted axes
        set(gca,'color','none','handlevisibility',axesVisible,'visible',axesVisible)
    
        % perform the transformation
        character = fallTransformation * character;
        pause(0.01);
    
        %gif

        % perform final setup for the animation
        set(h_rr,'Visible','off')  % This line erases the image of the Road Runner and Wile E. Coyote
        axis([0 70 0 70]) % This let me set the scale I wanted in the inserted axes
        set( gca, 'color','none','handlevisibility','off','visible','off')
        
    end

    % upon landing compress the character slightly to mimick a energy capture after landing
    for j = 1:3
        % setup the plot for the animation frame
        hb = axes('units','normalized', 'position',[-0.2 0.0625 1.2 1]);
        h_rr = plot(hb,character(1,:), character(2,:),   '.', 'color', ninjaColor, 'MarkerSize', 1);
        axis([0 70 0 70]) %This let me set the scale I wanted in the inserted axes
        set(gca,'color','none','handlevisibility',axesVisible,'visible',axesVisible)
    
        % perform the transformation
        character = compressionTransformation * character;
        pause(0.01);

        %gif
    
        % perform final setup for the animation
        set(h_rr,'Visible','off')  % This line erases the image of the Road Runner and Wile E. Coyote
        axis([0 70 0 70]) % This let me set the scale I wanted in the inserted axes
        set( gca, 'color','none','handlevisibility','off','visible','off')
    end

    % decompress to stand back up 
    for j = 1:3
        % setup the plot for the animation frame
        hb = axes('units','normalized', 'position',[-0.2 0.0625 1.2 1]);
        h_rr = plot(hb,character(1,:), character(2,:),   '.', 'color', ninjaColor, 'MarkerSize', 1);
        axis([0 70 0 70]) %This let me set the scale I wanted in the inserted axes
        set(gca,'color','none','handlevisibility',axesVisible,'visible',axesVisible)
    
        % perform the transformation
        character = decompressionTransformation * character;
        pause(0.01);
    
        %gif

        % perform final setup for the animation
        set(h_rr,'Visible','off')  % This line erases the image of the Road Runner and Wile E. Coyote
        axis([0 70 0 70]) % This let me set the scale I wanted in the inserted axes
        set( gca, 'color','none','handlevisibility','off','visible','off')
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % the main emphasise is on the throwing star but I need the character to be visible.
    % I plot the character here and later set his visibility to off
    hb = axes('units','normalized', 'position',[-0.2 0.0625 1.2 1]);
    h_rrCharacterBackground = plot(hb, character(1,:), character(2,:), '.', 'color', ninjaColor, 'MarkerSize', 1);
    axis([0 70 0 70]) %This let me set the scale I wanted in the inserted axes
    set(gca,'color','none','handlevisibility',axesVisible,'visible',axesVisible)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    throwingStar1 = moveToCharacterHand(character, throwingStar1);
    
    
    % throw ninja star at first target
    for j = 1:8
        % setup the plot for the animation frame
        hb = axes('units','normalized', 'position',[-0.2 0.0625 1.2 1]);
        h_rr = plot(hb, throwingStar1(1,:), throwingStar1(2,:), '.', 'color', ninjaStarColor, 'MarkerSize', 1);
        axis([0 70 0 70]) %This let me set the scale I wanted in the inserted axes
        set(gca,'color','none','handlevisibility',axesVisible,'visible',axesVisible)
    
        % perform the transformation
        throwingStar1 = throwingTransformation1 * throwingStar1;
        throwingStar1 = RotationScene(throwingStar1, -0.8);
        %throwingStar1 = RotationScene();
        pause(0.01);
    
        %gif

        % perform final setup for the animation
        set(h_rr,'Visible','off')  % This line erases the image of the Road Runner and Wile E. Coyote
        axis([0 70 0 70]) % This let me set the scale I wanted in the inserted axes
        set( gca, 'color','none','handlevisibility','off','visible','off')
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % draw throwingStar1 to the screen as background image
    % I will never touch this again as the throwing star will stay right where it has landed
    hb = axes('units','normalized', 'position',[-0.2 0.0625 1.2 1]);
    h_rrThrowingStar1Background = plot(hb, throwingStar1(1,:), throwingStar1(2,:), '.', 'color', ninjaStarColor, 'MarkerSize', 1);
    axis([0 70 0 70]) %This let me set the scale I wanted in the inserted axes
    set(gca,'color','none','handlevisibility',axesVisible,'visible',axesVisible)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    throwingStar2 = moveToCharacterHand(character, throwingStar2);

    % throw ninja star at second target
    for j = 1:46
         % setup the plot for the animation frame
         hb = axes('units','normalized', 'position',[-0.2 0.0625 1.2 1]);
         h_rr = plot(hb, throwingStar2(1,:), throwingStar2(2,:),   '.', 'color', ninjaStarColor, 'MarkerSize', 1);
         axis([0 70 0 70]) %This let me set the scale I wanted in the inserted axes
         set(gca,'color','none','handlevisibility',axesVisible,'visible',axesVisible)
     
         % perform the transformation
         throwingStar2 = throwingTransformation2 * throwingStar2;
         throwingStar2 = RotationScene(throwingStar2, -0.8);
         pause(0.05);
     
        %gif

         % perform final setup for the animation
         set(h_rr,'Visible','off')  % This line erases the image of the Road Runner and Wile E. Coyote
         axis([0 70 0 70]) % This let me set the scale I wanted in the inserted axes
         set( gca, 'color','none','handlevisibility','off','visible','off')
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % destroy the background character as he is now the main focus

    % perform final setup for the animation
    set(h_rrCharacterBackground,'Visible','off')  % This line erases the image of the Road Runner and Wile E. Coyote
    axis([0 70 0 70]) % This let me set the scale I wanted in the inserted axes
    set( gca, 'color','none','handlevisibility','off','visible','off')

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % draw throwingStar2 to the screen as background image
    % I will never touch this again as the throwing star will stay right where it has landed
    hb = axes('units','normalized', 'position',[-0.2 0.0625 1.2 1]);
    h_rrThrowingStar2Background = plot(hb, throwingStar2(1,:), throwingStar2(2,:), '.', 'color', ninjaStarColor, 'MarkerSize', 1);
    axis([0 70 0 70]) %This let me set the scale I wanted in the inserted axes
    set(gca,'color','none','handlevisibility',axesVisible,'visible',axesVisible)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


     % run to middle of the scene
     for j = 1:33
        % setup the plot for the animation frame
        hb = axes('units','normalized', 'position',[-0.2 0.0625 1.2 1]);
        h_rr = plot(hb, character(1,:), character(2,:),   '.', 'color', ninjaColor, 'MarkerSize', 1);
        axis([0 70 0 70]) %This let me set the scale I wanted in the inserted axes
        set(gca,'color','none','handlevisibility',axesVisible,'visible',axesVisible)
    
        % perform the transformation
        character = runningTransformation * character;
        pause(0.02);
    
        %gif

        % perform final setup for the animation
        set(h_rr,'Visible','off')  % This line erases the image of the Road Runner and Wile E. Coyote
        axis([0 70 0 70]) % This let me set the scale I wanted in the inserted axes
        set( gca, 'color','none','handlevisibility','off','visible','off')
    end


    

    % perform final setup for the animation
    set(h_rrCharacterBackground,'Visible','off')  % This line erases the image of the Road Runner and Wile E. Coyote
    axis([0 70 0 70]) % This let me set the scale I wanted in the inserted axes
    set( gca, 'color','none','handlevisibility','off','visible','off')


    characterCenter = centerPivot(character);
    failureFlag = false;
end


% This function takes a character and morphs into a different shape specified by the caller
% NOTE: this function only performs one step of the morph it must be called regularly untill the desired image is created.
%       input:
%               originalImage : the matrix containing the original Image to be morphed
%               templateImage : the image to transform the original image into
%
%       output:
%               the result of step the of the morph

function outputImage = morph(originalImage, templateImage, mixingProportion)
    outputImage = (1-mixingProportion) * originalImage + mixingProportion * templateImage;
end



% This funtion moves whatever sprite is passed in to the characters hand
%       input:
%               character : the character matrix
%               sprite    : the sprite to move to the characters hand
%                               NOTE: sprite must be in a homogenous coordinate system.
%
%       output:
%               outputSprite = The resulting sprite in the proper position

function sprite = moveToCharacterHand(character, sprite)
    % get the center of the character
    characterCenter = centerPivot(character);

    % move to the center of the character
    sprite = teleportTo(sprite, characterCenter(1,:), characterCenter(2,:));

    %translation matrix for moving from the center of the character to the hand
    translateToHand = [1 0 2; 0 1 1; 0 0 1];

    % translate to the characters hand
    sprite = translateToHand * sprite;
end