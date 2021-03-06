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
S = [0.025 0 0; 0 0.025 0; 0 0 1];  
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

