%% Online/Offline Score Informed Music Signal Decomposition: Application to Minus One
%
% This is the version of the code that reproduces the result in the article
% "Online/Offline Score Informed Music Signal Decomposition: Application to
% Minus One" submitted to EURASIP Journal on Audio, Speech, and Music
% Processing.
%
% Authors: Antonio Jesus Munoz-Montoro, Julio Jose Carabias-Orti, Pedro
%          Vera-Candeas, Francisco Jesus Canadas-Quesada and Nicolas Ruiz-Reyes
clear; clc;
%% Load user parameters
userParameters;

%% Source Separation algorithm
cd mono_sss;
performSS;
