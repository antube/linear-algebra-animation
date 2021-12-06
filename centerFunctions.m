function PPout = fJpeg2pointsConvert(BB,THRESHOLD)
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

function eraseChar(CB)
    %% Erase a character by filling a white box over it.
    uX = max(CB(1,:))+3;
    lX = min(CB(1,:))-3;
    uY = max(CB(2,:))+3;
    lY = min(CB(2,:))-3;
    erase = [lX lX uX uX; lY uY uY lY];
    fill(erase(1,:),erase(2,:),  'w' , 'EdgeColor', 'w' );
end



function cent = centerPivot(PP)
    %% Find the center of a character. Use this later for rotations.
    % Assume these points are moved into a scene frame.
    uX = max(PP(1,:));
    lX = min(PP(1,:));
    uY = max(PP(2,:));
    lY = min(PP(2,:));
    cent = [ mean([uX,lX])  ; mean([uY,lY]) ; 0];
end

function feet = feetPivot(PP)
    %% Find a center between a character's feet. 
    % Use this later for sway, walking, etc
    % Assume these points are moved into a scene frame.
    uX = max(PP(1,:));
    lX = min(PP(1,:));
    lY = min(PP(2,:));
    cent = [ mean([uX,lX])  ; lY ; 0];
end


