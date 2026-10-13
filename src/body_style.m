function [col, msize] = body_style(name)
%BODY_STYLE  Plot colour and marker size for a Solar System body.
switch name
    case 'Sun',     col = [1.00 0.78 0.10]; msize = 26;
    case 'Mercury', col = [0.60 0.60 0.60]; msize = 8;
    case 'Venus',   col = [0.90 0.70 0.35]; msize = 10;
    case 'Earth',   col = [0.20 0.45 0.90]; msize = 11;
    case 'Moon',    col = [0.45 0.45 0.45]; msize = 5;
    case 'Mars',    col = [0.85 0.30 0.15]; msize = 10;
    case 'Jupiter', col = [0.80 0.55 0.35]; msize = 18;
    case 'Saturn',  col = [0.85 0.75 0.45]; msize = 16;
    case 'Uranus',  col = [0.40 0.80 0.85]; msize = 14;
    case 'Neptune', col = [0.25 0.35 0.85]; msize = 14;
    otherwise,      col = [0.10 0.10 0.10]; msize = 14;
end
end
