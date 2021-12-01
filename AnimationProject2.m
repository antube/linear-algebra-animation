CA=imread('NinjaSword1.jpg');%first ninja image
CAout=fJpeg2pointsConverter(CA,150);%converts first ninja to data points
CD=imread('NinjaSword3.jpg');%second ninja image
CDout=fJpeg2pointsConverter(CD,150);%converts second ninja to data points
CE=imread('NinjaTool1.jpg');%third ninja image
CEout=fJpeg2pointsConverter(CE,150);%converts third ninja to data points
CF=imread('NinjaTool2.jpg');%fourth ninja image
CFout=fJpeg2pointsConverter(CF,150);%converts fourth ninja to data points
CG=imread('NinjaTool3.jpg');%fifth ninja image
CGout=fJpeg2pointsConverter(CG,150);%converts fifth ninja to data points
CH=imread('NinjaTool4.jpg');%sixth ninja image
CHout=fJpeg2pointsConverter(CH,150);%converts sixth ninja to data points
CI=imread('SmokeBomb.jpg');%seventh ninja image
CIout=fJpeg2pointsConverter(CI,150);%converts seventh ninja to data points
A=CAout;
D=CDout;
E=CEout;
F=CFout;
G=CGout;
H=CHout;
I=CIout;
Z=zeros(2,297);%makes shure the diminsions of the ninja is the same for matrix multiplication
CAout_New= [A,Z];%makes shure the diminsions of the ninja is the same for matrix multiplication
Z2=zeros(2,948);%makes shure the diminsions of the ninja is the same for matrix multiplication
CEout_New= [E,Z2];%makes shure the diminsions of the ninja is the same for matrix multiplication
Z3=zeros(2,2612);%makes shure the diminsions of the ninja is the same for matrix multiplication
CFout_New=[F,Z3];%makes shure the diminsions of the ninja is the same for matrix multiplication
Z4=zeros(2,1539);%makes shure the diminsions of the ninja is the same for matrix multiplication
CFout_New2=[F,Z4];%makes shure the diminsions of the ninja is the same for matrix multiplication
Z5=zeros(2,1832);%makes shure the diminsions of the ninja is the same for matrix multiplication
CGout_New=[G,Z5];%makes shure the diminsions of the ninja is the same for matrix multiplication
Z6=zeros(2,1832);%makes shure the diminsions of the ninja is the same for matrix multiplication
CGout_New2=[G,Z6];%makes shure the diminsions of the ninja is the same for matrix multiplication
Z7=zeros(2,490);%makes shure the diminsions of the ninja is the same for matrix multiplication
CIout_New=[I,Z7];%makes shure the diminsions of the ninja is the same for matrix multiplication
for k=0:1/4:1%speed of the morph
	B = (1-k)*CAout_New + k*D;%morphs the ninja
	plot(B(1,:),B(2,:),'.')%plot the new morphed ninja
	pause(0.2);
    %the ninja has morphed into sword mode
end
for J=0:1/4:1%speed of the morph
	B2 = (1-J)*D+J*CEout_New;%morphs the ninja
	plot(B2(1,:),B2(2,:),'.')%plot the new morphed ninja
	pause(0.2);
    %the ninja has morphed into shurekin mode
end
for L=0:1/4:1%speed of the morph
	B3 = (1-L)*CEout_New+L*CFout_New;%morphs the ninja
	plot(B3(1,:),B3(2,:),'.')%plot the new morphed ninja
	pause(0.2);
    %the ninja has morphed into kunai mode
end
for M=0:1/4:1%speed of the morph
	B4 = (1-M)*CFout_New2+M*CGout;%morphs the ninja
	plot(B4(1,:),B4(2,:),'.')%plot the new morphed ninja
	pause(0.2);
    %the ninja has morphed into nunchuck mode
end
for I=0:1/4:1%speed of the morph
	B5 = (1-I)*CGout_New2+I*CHout;%morphs the ninja
	plot(B5(1,:),B5(2,:),'.')%plot the new morphed ninja
	pause(0.2);
    %the ninja has morphed into staff mode
end
for Q=0:1/4:1%speed of the morph
	B6 = (1-Q)*CHout+Q*CIout_New;;%morphs the ninja
	plot(B6(1,:),B6(2,:),'.')%plot the new morphed ninja
	pause(0.2);
    %the ninja has morphed into smoke bomb mode
end
function PPout = fJpeg2pointsConverter(BB,THRESHOLD)
%% This function will take in an N x M x 3 matrix that 
% has been imported into the workspace using the 
% imread('filename.jpg') command and stored in a matrix
% - it is called BB inside this converter.  
%USAGE: BBout = fJpeg2pointsConvert(BB,THRESHOLD)
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