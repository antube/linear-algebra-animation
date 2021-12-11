%{
 sceneTwo.m

%}


%% Create background image
clf %
ha = axes('units','normalized', 'position',[0 0 1 1]);
uistack(ha,'bottom');
I=imread('NinjaHome.jpg');
hi = imagesc(I);
colormap gray;
set(ha,'handlevisibility','off', 'visible','off')
filename = 'NinjaSword1.jpg';
ninjaColor =[1, 1, 1];
thresh = 219;
ninjasword1 = imread(filename);
ns1mtx = fJpeg2pointsConverter(ninjasword1, thresh);
[m,n]=size(ns1mtx);
fprintf("%s size (thresh=%i) , [%i,%i]",filename,thresh,m,n);
disp(m);  disp(n); 
ns1mtx = [ns1mtx;ones(1,n)];
S = [0.025 0 0; 0 0.025 0; 0 0 1];  %This is my rescaling matrix to shrink the character to fit the background
ns1mtx = S*ns1mtx;
ns1mtx_orig = ns1mtx;
%% Play background music throughout all scenes.
% Time the animation and match the length of animation with length of audio
% [y,Fs] = audioread('ninja_music.wav');
% player = audioplayer(y,Fs);
% play(player)   % Start the player
% stop(player)   % Stop whenever you like...

axesVisible = 'on'; 
axesXpos = 0;
axesYpos = 0;
axesXdim = 1;
axesYdim = 1; 



ns1mtx = teleportTo(ns1mtx,35,25);

%% Lands on to roof
for i=1:5
    hb = axes('units','normalized', 'position',[-0.2 .0625 1.2 1]);
    h_rr = plot(hb,ns1mtx(1,:), ns1mtx(2,:),   '.', 'color', ninjaColor, 'MarkerSize', 1); 
    axis([0 70 0 70]) ;
    set(gca,'color','none','handlevisibility',axesVisible,'visible',axesVisible)
    
    nS = [1 0 0.5 ; 0 1 -0.4; 0 0 1 ];
    ns1mtx = nS*ns1mtx;
    
    pause(0.1);
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
    
    pause(0.1);
    set(h_rr,'Visible','off');  
    axis([0 70 0 70]) ;
    set( gca, 'color','none','handlevisibility','off','visible','off');
end

algn = alignWith(ns1mtx, ns1mtx_orig);
ns1mtx = algn;
for i=1:4
    hb = axes('units','normalized', 'position',[-0.2 .0625 1.2 1]);
    h_rr = plot(hb,ns1mtx(1,:), ns1mtx(2,:),   '.', 'color', ninjaColor, 'MarkerSize', 1); 
    axis([0 70 0 70]) ;
    set(gca,'color','none','handlevisibility',axesVisible,'visible',axesVisible)
    
    nS = [1 0 0.5 ; 0 1 0; 0 0 1 ];
    ns1mtx = nS*ns1mtx;
    
    pause(0.1);
    set(h_rr,'Visible','off');  
    axis([0 70 0 70]) ;
    set( gca, 'color','none','handlevisibility','off','visible','off');
end

nt4mtx = loadNinjaTool4('NinjaTool4.jpg');
Z = (-1)*centerPivot(nt4mtx);
nt4mtx = ShiftScene(nt4mtx, Z(1),Z(2));
nt4mtx = [-1 ,0 0; 0 -1 0; 0 0 1]*nt4mtx;
algn = alignWith(ns1mtx , nt4mtx);
nt4mtx = algn;


for i=1:40
    hb = axes('units','normalized', 'position',[-0.2 .0625 1.2 1]);
    h_rr = plot(hb,nt4mtx(1,:), nt4mtx(2,:),   '.', 'color', ninjaColor, 'MarkerSize', 1); 
    axis([0 70 0 70]) ;
    set(gca,'color','none','handlevisibility',axesVisible,'visible',axesVisible)
    
    nS = [1 0 0.5 ; 0 1 0; 0 0 1 ];
    nt4mtx = nS*nt4mtx;
    
    pause(0.1);
    set(h_rr,'Visible','off');  
    axis([0 70 0 70]) ;
    set( gca, 'color','none','handlevisibility','off','visible','off');
end

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


