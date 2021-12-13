
%%%%%%%%%%%%%%%%%%%%%%%%%%%   New Function %%%%%%%%%%%%%%%%%%%%%%%%%%%
% third_scene
% Parameters
%       characterCenter: A vector representing the center of the character
%
% Outputs
%       failureFlag: A flag which when set will indicate the running of the method failed
%       characterCenter: A vector representing the center of the charcter after the third scene is complete

function  [failureFlag, character, characterCenter, throwingStar1, throwingStar2] = third_scene(character, characterCenter, throwingStar1, throwingStar2, ninjaColor, axesVisible)
    %gif('third_scene.m.gif');
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Setup the nessecary matrices
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % have the character fall into the scene
    fallTransformation = [1 0 0.2; 0 1 -1; 0 0 1];

    % landing matrices.
    compressionTransformation = [1 0 0; 0 0.90 0; 0 0 1];
    decompressionTransformation = inv(compressionTransformation);


    % This transformation matrix  is used to move the ninja stars across the scene
    throwingTransformation = [1 0 1; 0 1 0; 0 0 1];


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
    for j = 1:11
        % setup the plot for the animation frame
        hb = axes('units','normalized', 'position',[-0.2 0.0625 1.2 1]);
        h_rr = plot(hb, throwingStar1(1,:), throwingStar1(2,:), '.', 'color', ninjaColor, 'MarkerSize', 1);
        axis([0 70 0 70]) %This let me set the scale I wanted in the inserted axes
        set(gca,'color','none','handlevisibility',axesVisible,'visible',axesVisible)
    
        % perform the transformation
        throwingStar1 = throwingTransformation * throwingStar1;
        pause(0.01);
    
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
    % draw throwingStar1 to the screen as background image
    % I will never touch this again as the throwing star will stay right where it has landed
    hb = axes('units','normalized', 'position',[-0.2 0.0625 1.2 1]);
    h_rrThrowingStar1Background = plot(hb, throwingStar1(1,:), throwingStar1(2,:), '.', 'color', ninjaColor, 'MarkerSize', 1);
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


    % setup the plot for the animation frame
    hb = axes('units','normalized', 'position',[-0.2 0.0625 1.2 1]);
    h_rrCharacterBackground = plot(hb, character(1,:), character(2,:), '.', 'color', ninjaColor, 'MarkerSize', 1);
    axis([0 70 0 70]) %This let me set the scale I wanted in the inserted axes
    set(gca,'color','none','handlevisibility',axesVisible,'visible',axesVisible)


    %{
    throwingStar2 = moveToCharacterHand(character, throwingStar2);

    % throw ninja star at second target
    for j = 1:4
         % setup the plot for the animation frame
         hb = axes('units','normalized', 'position',[-0.2 0.0625 1.2 1]);
         h_rr = plot(hb,character(1,:), character(2,:),   '.', 'color', ninjaColor, 'MarkerSize', 1);
         axis([0 70 0 70]) %This let me set the scale I wanted in the inserted axes
         set(gca,'color','none','handlevisibility',axesVisible,'visible',axesVisible)
     
         % perform the transformation
         throwingStar2 = throwingTransformation * throwingStar2;
         pause(0.05);
     
        gif

         % perform final setup for the animation
         set(h_rr,'Visible','off')  % This line erases the image of the Road Runner and Wile E. Coyote
         axis([0 70 0 70]) % This let me set the scale I wanted in the inserted axes
         set( gca, 'color','none','handlevisibility','off','visible','off')
    end
    %}

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

function outputSprite = moveToCharacterHand(character, sprite)
    % get the center of the character
    characterCenter = centerPivot(character);

    % move to the center of the character
    outputSprite = teleportTo(sprite, characterCenter(1,:), characterCenter(2,:));

    %translation matrix for moving from the center of the character to the hand
    translateToHand = [1 0 7; 0 1 5; 0 0 1];

    % translate to the characters hand
    outputsprite = translateToHand * outputSprite;
end





function cent = centerPivot(PP)
    % Assume these points are moved into a scene frame.
    uX = max(PP(1,:));
    lX = min(PP(1,:));
    uY = max(PP(2,:));
    lY = min(PP(2,:));
    cent = [ mean([uX,lX])  ; mean([uY,lY]) ; 0];
end


function PPt = teleportTo(PP,tx,ty)
    nc = centerPivot(PP); 
    nP = [1 0 -1*nc(1) ; 0 1 -1*nc(2); 0 0 1 ];
    zPP = nP*PP;
    nS = [1 0 tx ; 0 1 ty; 0 0 1 ];
    PPt = nS*zPP;
end