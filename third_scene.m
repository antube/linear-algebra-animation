
%%%%%%%%%%%%%%%%%%%%%%%%%%%   New Function %%%%%%%%%%%%%%%%%%%%%%%%%%%
% third_scene
% Parameters
%       characterCenter: A vector representing the center of the character
%
% Outputs
%       failureFlag: A flag which when set will indicate the running of the method failed
%       characterCenter: A vector representing the center of the charcter after the third scene is complete

function  [failureFlag, character, throwingStar1, throwingStar2] = third_scene(character, throwingStar1, throwingStar2, ninjaColor, axesVisible)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Setup the nessecary matrices
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % I need a centerpivot to call the transform and animate function
    starTemplate = throwingStar1


    % have the character fall into the scene
    fallTransformation = [1 0 0.2; 0 1 -1; 0 0 1];

    % landing matrices.
    compressionTransformation = [1 0 0; 0 0.95 0; 0 0 1];
    decompressionTransformation = inv(compressionTransformation);


    % This transformation matrix  is used to move the ninja stars across the scene
    throwingTransformation = [1 0 1; 0 1 0; 0 0 0];


    % running transformation matrix
    runningTransformation = [1 0 1; 0 1 0; 0 0 0];

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
        %h_rr = plot(hb, 35, 35,   '*', 'color', ninjaColor, 'MarkerSize', 35); 
        axis([0 70 0 70]) %This let me set the scale I wanted in the inserted axes
        set(gca,'color','none','handlevisibility',axesVisible,'visible',axesVisible)
    
        % perform the transformation
        character = fallTransformation * character;
        pause(0.05);
    
        % perform final setup for the animation
        set(h_rr,'Visible','off')  % This line erases the image of the Road Runner and Wile E. Coyote
        axis([0 70 0 70]) % This let me set the scale I wanted in the inserted axes
        set( gca, 'color','none','handlevisibility','off','visible','off')
        
    end

    % upon landing compress the character slightly to mimick a energy capture after landing
    for j = 1:4
        % setup the plot for the animation frame
        hb = axes('units','normalized', 'position',[-0.2 0.0625 1.2 1]);
        h_rr = plot(hb,character(1,:), character(2,:),   '.', 'color', ninjaColor, 'MarkerSize', 1);
        %h_rr = plot(hb, 35, 35,   '*', 'color', ninjaColor, 'MarkerSize', 35); 
        axis([0 70 0 70]) %This let me set the scale I wanted in the inserted axes
        set(gca,'color','none','handlevisibility',axesVisible,'visible',axesVisible)
    
        % perform the transformation
        character = compressionTransformation * character;
        pause(0.05);
    
        % perform final setup for the animation
        set(h_rr,'Visible','off')  % This line erases the image of the Road Runner and Wile E. Coyote
        axis([0 70 0 70]) % This let me set the scale I wanted in the inserted axes
        set( gca, 'color','none','handlevisibility','off','visible','off')
    end

    % decompress to stand back up 
    for j = 1:4
        % setup the plot for the animation frame
        hb = axes('units','normalized', 'position',[-0.2 0.0625 1.2 1]);
        h_rr = plot(hb,character(1,:), character(2,:),   '.', 'color', ninjaColor, 'MarkerSize', 1);
        %h_rr = plot(hb, 35, 35,   '*', 'color', ninjaColor, 'MarkerSize', 35); 
        axis([0 70 0 70]) %This let me set the scale I wanted in the inserted axes
        set(gca,'color','none','handlevisibility',axesVisible,'visible',axesVisible)
    
        % perform the transformation
        character = decompressionTransformation * character;
        pause(0.05);
    
        % perform final setup for the animation
        set(h_rr,'Visible','off')  % This line erases the image of the Road Runner and Wile E. Coyote
        axis([0 70 0 70]) % This let me set the scale I wanted in the inserted axes
        set( gca, 'color','none','handlevisibility','off','visible','off')
    end





    %morph first ninja star from point to star
    throwingStar1 = morph(throwingStar1, starTemplate, 1);
    
    % throw ninja star at first target
    for j = 1:4
        % setup the plot for the animation frame
        hb = axes('units','normalized', 'position',[-0.2 0.0625 1.2 1]);
        h_rr = plot(hb,character(1,:), character(2,:),   '.', 'color', ninjaColor, 'MarkerSize', 1);
        %h_rr = plot(hb, 35, 35,   '*', 'color', ninjaColor, 'MarkerSize', 35); 
        axis([0 70 0 70]) %This let me set the scale I wanted in the inserted axes
        set(gca,'color','none','handlevisibility',axesVisible,'visible',axesVisible)
    
        % perform the transformation
        throwingStar1 = throwingTransformation * throwingStar1;
        pause(0.05);
    
        % perform final setup for the animation
        set(h_rr,'Visible','off')  % This line erases the image of the Road Runner and Wile E. Coyote
        axis([0 70 0 70]) % This let me set the scale I wanted in the inserted axes
        set( gca, 'color','none','handlevisibility','off','visible','off')
    end






    % run to middle of scene
    for j = 1:4
        % setup the plot for the animation frame
        hb = axes('units','normalized', 'position',[-0.2 0.0625 1.2 1]);
        h_rr = plot(hb,character(1,:), character(2,:),   '.', 'color', ninjaColor, 'MarkerSize', 1);
        %h_rr = plot(hb, 35, 35,   '*', 'color', ninjaColor, 'MarkerSize', 35); 
        axis([0 70 0 70]) %This let me set the scale I wanted in the inserted axes
        set(gca,'color','none','handlevisibility',axesVisible,'visible',axesVisible)
    
        % perform the transformation
        character = runningTransformation * character;
        pause(0.05);
    
        % perform final setup for the animation
        set(h_rr,'Visible','off')  % This line erases the image of the Road Runner and Wile E. Coyote
        axis([0 70 0 70]) % This let me set the scale I wanted in the inserted axes
        set( gca, 'color','none','handlevisibility','off','visible','off')
    end



    

    % morph second ninja star from point to star
    
    throwingStar2 = morph(throwingStar2, starTemplate, 1);

    % throw ninja star at second target
    for j = 1:4
         % setup the plot for the animation frame
         hb = axes('units','normalized', 'position',[-0.2 0.0625 1.2 1]);
         h_rr = plot(hb,character(1,:), character(2,:),   '.', 'color', ninjaColor, 'MarkerSize', 1);
         %h_rr = plot(hb, 35, 35,   '*', 'color', ninjaColor, 'MarkerSize', 35); 
         axis([0 70 0 70]) %This let me set the scale I wanted in the inserted axes
         set(gca,'color','none','handlevisibility',axesVisible,'visible',axesVisible)
     
         % perform the transformation
         throwingStar2 = throwingTransformation * throwingStar2;
         pause(0.05);
     
         % perform final setup for the animation
         set(h_rr,'Visible','off')  % This line erases the image of the Road Runner and Wile E. Coyote
         axis([0 70 0 70]) % This let me set the scale I wanted in the inserted axes
         set( gca, 'color','none','handlevisibility','off','visible','off')
    end
    


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