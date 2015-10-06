function meanGMD = wrapper(allparams)
EEG = pop_loadset('filename','violaine.set','filepath','/home/frans/documents/Microstates/');
OUTEEG = pop_getmicrostates(EEG,2,2,20,2, allparams(1), allparams(2)); 
meanGMD = OUTEEG.meanGMD;
end