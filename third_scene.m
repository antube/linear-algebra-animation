
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



    % character falls into scene from the previous building jump
    % character lands
    % character throws ninja star at target
    % runs past target
    % throws ninja star at second target

    % have the character fall into the scene
    fallTransformation = [];

    for j = 1:4
        [character, characterCenter] = transformAndAnimate(character, characterCenter, fallTransformation);

        pause(0.2);
    end

    % landing matrices.
    compressionTransformation = [];
    decompressionTransformation = inv(compressionTransformation);

    % upon landing compress the character slightly to mimick a energy capture after landing
    for j = 1:4
        [character, characterCenter] = transformAndAnimate(character, characterCenter, compressionTransformation);
    end

    % decompress to stand back up 
    for j = 1:4

    end


    % This transformation matrix  is used to move the ninja stars across the scene
    ninjaStarTransformation = [];


    %morph first ninja star from point to star
    %throwingStar1 = morph(throwingStar1, starTemplate);
    

    % throw ninja star at first target
    for j = 1:4

    end


    % running transformation matrix
    runningTransformation = [];


    % run to middle of scene
    for j = 1:4

    end


    % morph second ninja star from point to star
    %throwingStar2 = morph(throwingStar2, starTemplate);

    % throw ninja star at second target
    for j = 1:4

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
    character = transformation * character;

    characterCenter = centerPivot(character);

    % animation magic here.
    
