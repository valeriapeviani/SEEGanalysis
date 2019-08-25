# SEEGanalysis

The code has been written for a block-design experiment with 4 blocks (corresponding to 2 conditions) in which participants performed a line length judgment judgment task.

The available MATLAB scripts allow to: 
1) preprocess the data (PARAMS + Preprocessing)
2) extract power (high-gamma and beta) via filter-hilbert (PARAMS + TFA_Hilbert)
3) run permutation statistics to select 'responsive channels' (channels in which a change in power has been detected after the stimulus onset) (PARAMS + Statistics)
4) run permutation statistics to compare conditions (PARAMS + Statistics)

the PARAMS script collects all the paramters you need to perform the preprocessing and analysis. 
To be noted that the function Check_triggers, Populatetrialinfo and Trialfun_content need to be customized (and they are quite hardcoded, sorry about that!)
