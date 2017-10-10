% lipid bilayer experiment with a population of 250 RyR's

time = 0;
time_end = 3000;
time_step = 0.01;

% Declare model parameters
 global Nclefts_FRU Nstates_FRU Nstates_FRUdep Nstates_LType Nstates_RyR NRyRs_per_cleft Nindepstates_LType 
 global NFRU NFRU_scale N
 Nclefts_FRU = 4;
 Nstates_FRU = (1+Nclefts_FRU);
 Nstates_FRUdep = 3;

 Nstates_LType = 12;
 Nstates_RyR  = 6;
 NRyRs_per_cleft  = 5;
 Nindepstates_LType = 2;
 NFRU = 250;
 NFRU_scale = 250;
 
 N = 37;
% Read in the initial conditions for the RyR
 input_dir = 'ic/vclamp';
 ic_states_file = strcat(input_dir,'/','ic_states_NVC.txt');
 ic_FRU_file = strcat(input_dir, '/','ic_FRU_NVC.txt');
 ic_LCh_file = strcat(input_dir,'/','ic_LCh_NVC.txt');
 ic_RyR_file = strcat(input_dir,'/','ic_RyR_NVC.txt');
 ic_Ito2_file = strcat(input_dir,'/','ic_Ito2_NVC.txt');
 
 [state, FRU_states, LType_state, RyR_state, Ito2_state] = initialize(ic_states_file, ic_FRU_file, ic_LCh_file,ic_RyR_file,ic_Ito2_file);
 
 CaSSclamp_duration = 200.0;
 CaSSclamp_shift =  10.0;
 CaSSclamp_set = 0.001047733420000;
 CaSSclamp_hold  =  1.10064e-04;
 period = 1000.0;
 index = 1;

for time=0:time_step:time_end

    Time_vec(index) = time;
    if(mod(time,period) >= CaSSclamp_shift && mod(time, period) <= (CaSSclamp_duration + CaSSclamp_shift))
        CaSS(index) = CaSSclamp_set;
    else
        CaSS(index) = CaSSclamp_hold;
    end
    index = index + 1;
end
 
figure
plot(Time_vec, CaSS);