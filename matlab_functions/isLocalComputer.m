function islocal = isLocalComputer()

[~, computer] = system('hostname');
screenInfo.computer = computer(1:end-1);

if isequal(screenInfo.computer,'Ariels-MacBook-Pro.local')
    islocal = 1;
else
    islocal = 0;
end