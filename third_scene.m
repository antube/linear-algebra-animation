
%%%%%%%%%%%%%%%%%%%%%%%%%%%   New Function %%%%%%%%%%%%%%%%%%%%%%%%%%%
% third_scene
% Parameters
%       characterCenter: A vector representing the center of the character
%
% Outputs
%       failureFlag: A flag which when set will indicate the running of the method failed
%       characterCenter: A vector representing the center of the charcter after the third scene is complete

function  [failureFlag, character, characterCenter] third_scene(character, characterCenter)
    % no one else needs these throwing stars so I do not care if they are destroyed after third scene completion
    throwingStar1 = [];
    throwingStar2 = [];


    % have the character fall into the scene
    fallTransformation = [];


    % landing matrices.
    compressionTransformation = [];
    decompressionTransformation = inv(compressionTransformation);


    % This transformation matrix  is used to move the ninja stars across the scene
    throwingTransformation = [];


      % running transformation matrix
      runningTransformation = [];


    % character falls into scene from the previous building jump
    % character lands
    % character throws ninja star at target
    % runs past target
    % throws ninja star at second target

    

    for j = 1:4
        [character, characterCenter] = transformAndAnimate(character, characterCenter, fallTransformation);
        pause(0.25);
    end

    

    % upon landing compress the character slightly to mimick a energy capture after landing
    for j = 1:4
        [character, characterCenter] = transformAndAnimate(character, characterCenter, compressionTransformation);
        pause(0.25);
    end

    % decompress to stand back up 
    for j = 1:4
        [character, characterCenter] = transformAndAnimate(charcter, characterCenter, decompressionTransformation);
        pause(0.25);
    end


    


    %morph first ninja star from point to star
    throwingStar1 = morph(throwingStar1, starTemplate);
    

    % throw ninja star at first target
    for j = 1:4
        [character, characterCenter] = transformAndAnimate(character, characterCenter, throwingTransformation);
        pause(0.25);
    end


    % run to middle of scene
    for j = 1:4
        [character, characterCenter] = transformAndAnimate(character, characterCenter, runningTransformation);
        pause(0.25);
    end


    % morph second ninja star from point to star
    throwingStar2 = morph(throwingStar2, starTemplate);

    % throw ninja star at second target
    for j = 1:4
        [character, characterCenter] = transformAndAnimate(character, characterCenter, throwingTransformation)
        pause(0.25);
    end


    % run to the middle of the image
    for j = 1:4
        [character, characterCenter] = transformAndAnimate(character, characterCenter, runningTransformation);
        pause(0.25);
    end




    characterCenter = centerPivot(character);



% this function takes care of the drawing and the transforming
%       input:
%               character        : the matrix of the character vectors
%               characterCenter  : the center of the character
%               transformation   : the transformation matrix to use for the transformation
%
%       output:
%               character       : the character matrix after the transformation
%               characterCenter : the center of the character after the transformation

function [character, characterCenter] = transformAndAnimate(character, characterCenter, transformation)
    % setup the plot for the animation frame
    hb = axes('units','normalized', 'position',[-0.2 .0625 1 1]);
    %hb = axes('position',[axesXpos axesYpos axesXdim axesYdim]);
    h_rr = plot(hb,ns1mtx(1,:), ns1mtx(2,:),   '.', 'color', ninjaColor, 'MarkerSize', 1); 
    axis([0 70 0 70]) %This let me set the scale I wanted in the inserted axes
    set(gca,'color','none','handlevisibility',axesVisible,'visible',axesVisible)
    
    % perform the transformation
    character = transformation * character;
    characterCenter = centerPivot(character);
    

    % perform final setup for the animation
    set(h_rr,'Visible','off')  % This line erases the image of the Road Runner and Wile E. Coyote
    axis([0 70 0 70]) % This let me set the scale I wanted in the inserted axes
    set( gca, 'color','none','handlevisibility','off','visible','off')
    


% This function takes a character and morphs into a different shape specified by the caller
% NOTE: this function only performs one step of the morph it must be called regularly untill the desired image is created.
%       input:
%               originalImage : the matrix containing the original Image to be morphed
%               templateImage : the image to transform the original image into
%
%       output:
%               the result of step the of the morph

function outputImage = morph(originalImage, templateImage)
    B2 = (1-j) * originalImage + templateImage;