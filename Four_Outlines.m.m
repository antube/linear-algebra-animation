%{
First Scene.
%}

%Input:
    % Character position
    % Background music (Naruto theme song?)
        % Read audio file.
        [y,Fs] = audioread(filename);
        % Play it.
        sound = (y,Fs);
    
    %Output:
    % Character position
    % Flag to represent the code that functioned properly
    
for jj=1:0.5:10
    %Description of event:   
    % The ninja appears on the street.
    % The ninja runs across the street from left to right to
    % the base of a building.
    % Homogeneous coordinates will be used to translate the ninja horizontally
    % from left to right.
    % The ninja arrives at the base of the building.
    % The ninja will scale the building to the top.
    % Homogeneous coordinates will be use to vertically translate the ninja up the
    % building.
    % Also, the ninja will be slightly rotated back and forth as he travels
    % upward to create an allusion that he is climbing.
end

%{ 
Second Scene 

(Stephen Horn)
Ninja jumps from building to building doing flips. 
Warning: This script does not run. This is just boilerplate.
%}

% Inputs 
%   PP matrix of points for standing ninja, in homogeneous coordinates.
%   BKGRNDFILE   image file name of rooftop baackground scenery.
%   PPsz   dimensions of the ninja points, PP.
%   BGROUNDsz  dimensions of the image file.

% ?? BKGRNDFILE = "QQQQQ.jpg";

clf;  
ha = axes('units','normalized', 'position',[0 0 1 1]);
uistack(ha,'bottom');
I=imread(BKGRNDFILE);
hi = imagesc(I);


% Run with a tilt (shear)
PPrun = ShearHScene(PP,k);
for ii=1:0.5:10
    hb = axes('position',[0 0.1 0.625 0.625]);
    h_rr = plot(hb,PPrun(1,:), PPrun(2,:),'r.');
    axis([0 70 0 70]) %This let me set the scale I wanted in the inserted axes
    set(gca,'color','none','handlevisibility','on','visible','on')
    PPt = PPrun;
    PPrun = ShiftScene(PPt,0.1,0);
    PPt = PP;
    PP = ShiftScene(PPt,0.1,0);  % also move the original unsheared.
    pause(0.25)
    set(h_rr,'Visible','off')  
    axis([0 70 0 70])
    set(gca,'color','none','handlevisibility','off','visible','off') 
end

% Leap to next building with front flip.
for ii=1:0.5:10
    hb = axes('position',[0 0.1 0.625 0.625]);
    h_rr = plot(hb,PP(1,:), PP(2,:),'r.');
    axis([0 70 0 70]) %This let me set the scale I wanted in the inserted axes
    set(gca,'color','none','handlevisibility','on','visible','on')
    PPt = PP;
    PP = RotationScene(PPt,0.1,0); 
    PPt = PP;
    PP = ShiftScene(PPt,0.1);
    pause(0.25)
    set(h_rr,'Visible','off') 
    axis([0 70 0 70]) 
    set(gca,'color','none','handlevisibility','off','visible','off') 
end

% Jump and fall out of the frame below.
for ii=1:0.5:10
    hb = axes('position',[0 0.1 0.625 0.625]);
    h_rr = plot(hb,PP(1,:), PP(2,:),'r.');
    axis([0 70 0 70]) %This let me set the scale I wanted in the inserted axes
    set(gca,'color','none','handlevisibility','on','visible','on')
    %  Maybe morph some here? 
    PPt = PP;
    PP = ShiftScene(PPt,0.1);
    pause(0.25)
    set(h_rr,'Visible','off') 
    axis([0 70 0 70])
    set(gca,'color','none','handlevisibility','off','visible','off') % Turns axis off at end of shift
end





%{
Third Scene.
%}
%{
(Sbh)   As inputs, scene 3 will  need extra ninja images for 
    star throwing poses. 
%}

%throwing ninja stars
for j = 0:1
    
    % input
    %   character position
    
    % output
    %   character position
    %   flag to represent that the code actually functioned properly
    
    
    % character falls into scene from the previous building jump
    % character lands
    % character throws ninja stars at target
    % charcter then runs off scene
end

%{
Fourth Scene. 

    here for the final scean the ninja morphs into smake bomb mode
    and dissapers signaling the end of his journey.
%} 

CA=imread(%the before image of the ninja)
CAout=fJpeg2pointsConverter(CA,150);%converts before image of ninja to data points
CB=imread(%the after image of the ninja);
CBout=fJpeg2pointsConverter(CB,150);%converts after image of ninja to data points
CC=imread(%an after image of a blank page since the ninja will be morphing into nothing);
CCout=fJpeg2pointsConverter(CC,150);%converts after image of blank page to data points
Z=zeros(2,%insert the difference of units to insure a matrix multiplication)%makes shure the dimminsions of the multiplication matrices are the same
for i=0:1/4:1%speed of the morph
	B = (1-i)*%before image * i*after image;
    %morphs the ninja into smoke bomb mode
	plot(B(1,:),B(2,:),'.')%plot the new morphed ninja
	pause(0.2);%so we see the morph
end
for i=0:1/4:1%speed of the morph
	B = (1-i)*%the smoke bomb mode is now the before image * i*the blank page after image;
    %the ninja should basickly just dissapear now
	plot(B(1,:),B(2,:),'.')
	pause(0.2);
end


function PPout = fJpeg2pointsConverter(BB,THRESHOLD)
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

function cent = centerPivot(PP)
    % Assume these points are moved into a scene frame.
    uX = max(PP(1,:));
    lX = min(PP(1,:));
    uY = max(PP(2,:));
    lY = min(PP(2,:));
    cent = [ mean([uX,lX])  ; mean([uY,lY]) ; 0];
end

function fpiv = feetPivot(PP)
    % Get a pivot point at the feet of the character.
    uX = max(PP(1,:));
    lX = min(PP(1,:));
    %uY = max(PP(2,:));
    lY = min(PP(2,:));
    fpiv = [ mean([uX,lX])  ; lY ; 0];
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


function PPshh = ShearHScene(PP,k)
    [Mrows Ncols] = size(PP);
    if Mrows == 2,
        SH = [1 k ; 0 1];
    else ,
        SH = [1 k 0; 0 1 0; 0 0 1];
    end
    center = centerPivot(PP);
    PPz = ShiftScene(PP, -1.0*center(1,1), -1.0*center(2,1));
    PPshh = (SH*PPz) + center;
end 
    
    