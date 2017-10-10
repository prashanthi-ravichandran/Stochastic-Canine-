% Notes: Ramp has been removed 
%% Indices for state variables
 global index_V index_mNa index_jNa index_hNa index_Nai index_Ki index_Cai index_CaNSR index_xKs
 global index_LTRPNCa index_HTRPNCa index_C0Kv43 index_C1Kv43 index_C2Kv43 index_C3Kv43 index_OKv43 
 global index_CI0Kv43 index_CI1Kv43 index_CI2Kv43 index_CI3Kv43 index_OIKv43 index_C0Kv14 index_C1Kv14
 global index_C2Kv14 index_C3Kv14 index_OKv14 index_CI0Kv14 index_CI1Kv14 index_CI2Kv14 index_CI3Kv14 index_OIKv14 
 global index_CaTOT index_C1Herg index_C2Herg index_C3Herg index_OHerg index_IHerg
 
index_V         = 1;
index_mNa       = 2;
index_hNa       = 3;
index_jNa       = 4;
index_Nai       = 5;
index_Ki        = 6;
index_Cai       = 7;
index_CaNSR     = 8;
index_xKs       = 9;
index_LTRPNCa   = 10;
index_HTRPNCa   = 11;
index_C0Kv43    = 12;
index_C1Kv43    = 13;
index_C2Kv43    = 14;
index_C3Kv43    = 15;
index_OKv43     = 16;
index_CI0Kv43   = 17;
index_CI1Kv43   = 18;
index_CI2Kv43   = 19;
index_CI3Kv43   = 20;
index_OIKv43    = 21;
index_C0Kv14    = 22;
index_C1Kv14    = 23;
index_C2Kv14    = 24;
index_C3Kv14    = 25;
index_OKv14     = 26;
index_CI0Kv14   = 27;
index_CI1Kv14   = 28;
index_CI2Kv14   = 29;
index_CI3Kv14   = 30;
index_OIKv14    = 31;
index_CaTOT     = 32;
index_C1Herg    = 33;
index_C2Herg    = 34;
index_C3Herg    = 35;
index_OHerg     = 36;
index_IHerg     = 37;
% Model definitions
 global Nclefts_FRU Nstates_FRU Nstates_FRUdep Nstates_LType Nstates_RyR NRyRs_per_cleft Nindepstates_LType 
 global NFRU_sim NFRU_scale
 
Nclefts_FRU = 4;
Nstates_FRU = (1+Nclefts_FRU);
Nstates_FRUdep = 3;

Nstates_LType = 12;
Nstates_RyR  = 6;
NRyRs_per_cleft  = 5;
Nindepstates_LType = 2;
NFRU_sim_low  = 250;
NFRU_scale_low = 50.0; % ratio of 12500/NFRU_sim_low
NFRU_sim_high = 250;
NFRU_scale_high = 50.0; % ratio of 12500/NFRU_sim_high

NFRU_sim_max = 250;

NFRU_sim =NFRU_sim_high;
NFRU_scale = NFRU_scale_high;
% Specify simulation parameters
global step_min step_max tolrk
num_beats  = 4;
freq       = 1;
ISI        = 1000/freq; %ms
time_start = 0.0;
time_end = num_beats*ISI; %ms
%time_end = 500; %ms
time_step = 1.0; 
step_min = 1e-5; %ms
step_max = 0.1;
tolrk = 1e-6;
% AP Clamp
global pulse_duration pulse_amplitude  period shift APClampflag 
APClampflag = 0;
pulse_duration = 0.5;
pulse_amplitude = -100.0;
period = 1000.0;
shift = 5.0;

% Define parameters for voltage clamp.
 global vclamp_flag vclamp_duration vclamp_set vclamp_shift vclamp_hold vclamp_period
 vclamp_flag = 1;
 vclamp_duration = 200.0;
 vclamp_set  =   0.0;
 vclamp_shift =  10.0;
 vclamp_hold  =  -100.0;
 vclamp_period = 1000.0;
 
% Specify input and output files
% Define input files
 ic_clamp = 1;
 input_dir = 'ic/vclamp';
 ic_states_file = strcat(input_dir,'/','ic_states_NVC.txt');
 ic_FRU_file = strcat(input_dir, '/','ic_FRU_NVC.txt');
 ic_LCh_file = strcat(input_dir,'/','ic_LCh_NVC.txt');
 ic_RyR_file = strcat(input_dir,'/','ic_RyR_NVC.txt');
 ic_Ito2_file = strcat(input_dir,'/','ic_Ito2_NVC.txt');
 
 % Define output files
 output_dir           = 'vclamp(1Hz)';
 mkdir(output_dir);
 filenumber           = 1;
 output_states_file   = strcat(output_dir,'/','states', num2str(filenumber),'.txt');
 output_currents_file = strcat(output_dir,'/','currents', num2str(filenumber),'.txt');
 output_otherstates_file = strcat(output_dir,'/','otherstates', num2str(filenumber),'.txt');
 RyR_occupancy_file = strcat(output_dir,'/','RyR_occupancies', num2str(filenumber),'.txt');strcat(output_dir,'/','otherstates', num2str(filenumber),'.txt');
 fileID = fopen(output_states_file, 'w');
 fprintf(fileID, '%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s \n', ...
		'Time','V','mNa','hNa','jNa','Nai','Ki','Cai','CaNSR','xKs', ...
		'LTRPNCa','HTRPNCa','C0Kv43','C1Kv43','C2Kv43','C3Kv43','OKv43','CI0Kv43','CI1Kv43','CI2Kv43', ...
		'CI3Kv43','OIKv43','C0Kv14','C1Kv14','C2Kv14','C3Kv14','OKv14','CI0Kv14','CI1Kv14','CI2Kv14', ...
		'CI3Kv14','OIKv14','CaToT','C1Herg','C2Herg','C3Herg','OHerg','IHerg');
 fclose(fileID);   
 
 fileID = fopen(output_currents_file, 'w');
 fprintf(fileID, '%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s \n', ...
		'Time','INa','IKr','IKs','Ito1','IK1','IKp','INaCa','INaK','IpCa', ...
		'ICab','INab','ICa','JDHPR','Jup','Jtrpn','Jtr','Jxfer','IKv43','IKv14', ...
		'IKv14_K','IKv14_Na','Ito2','Istim','Itot');
 fclose(fileID);
 
 fileID = fopen(output_otherstates_file, 'w');
 
 fprintf(fileID, '%s %s %s %s %s %s %s %s %s %s %s %s %s \n', ...
			'Time','CaSSavg','CaJSRavg','JRyRtot','PRyR_open','PRyR_ready','PNorm_mode','PnotVinact','PLType_open','CaToT2', ...
		    'PIto2_open','CaJSRtot','CaSStot');
 fclose(fileID);
 
 fileID = fopen(RyR_occupancy_file, 'w');
  fprintf(fileID, '%s %s %s %s %s %s %s %s %s\n', ...
			'Time','1','2','3','4','5','6','7','8');
 fclose(fileID);
%Initialize state variables and gates
if(ic_clamp)
    [state, FRU_states, LType_state, RyR_state, Ito2_state] = initialize(ic_states_file, ic_FRU_file, ic_LCh_file,ic_RyR_file,ic_Ito2_file);
else
    state = [ -100.0100   1.21087e-4    0.999484    0.999480      10.0000     133.24    1.11074e-4   0.72873    0.104829e-3    0.1003...
                 0.9780      0.968277     0.0133601   0.691875e-4   0.159092e-6 0.0       0.0153235    0.00271424   0.243515e-3 ...
                 0.115007e-4    0.163239e-6 0.824239    0.0522865   0.00124396  0.131359e-4 0.522383e-7 0.118010    0.003334011 0.684631e-3...
                 0.136717e-3 0.451249e-4 0.6810488214E+01   0.990   0.008   0.002   0.0     0.0];
    LType_state = zeros(NFRU_sim, Nclefts_FRU,Nindepstates_LType);
    FRU_states = horzcat( 0.728921.*ones(NFRU_sim,1), 0.111074e-3.*ones(NFRU_sim, Nclefts_FRU));
    RyR_state = ones(NFRU_sim, Nclefts_FRU, NRyRs_per_cleft);
    Ito2_state = ones(NFRU_sim, Nclefts_FRU);
    LType_state(:,:,1) = ones(NFRU_sim, Nclefts_FRU);
    LType_state(:,:,2) = 2.*ones(NFRU_sim, Nclefts_FRU);  
end

% Set up arrays to store simulation data
global N Ncur Nother
N = 37;
Ncur = 24;
Nother = 12;
current     = zeros(Ncur,1);
otherstates = zeros(Nother,1);
Fstate = zeros(N,1);
ii = 1;
saveinterval = 1;
time_now = time_start;
oldstepsize = step_max;
fprintf('Time \t Voltage \t Cai\n');

% Define model geometry and physical constants
global Ko Nao Cao
global Faraday Temp Rgas RT_over_F 
global Acap VNSR VJSR VSS Vmyo
% Standard ionic concentrations
    Ko=    4.0;   % extracellular K+   concentration (mM)
    Nao= 138.0;   % extracellular Na+  concentration (mM)
    Cao=   2.0;   % extracellular Ca++ concentration (mM)

% Physical constants
    Faraday=  96.5;     % Faraday's constant (C/mmol)
    Temp=    310.0;     % absolute temperature (K)
    Rgas=      8.314;   % ideal gas constant (J/(mol*K))
    RT_over_F= (Rgas*Temp/Faraday); %  Rgas*Temp/Faraday (mV)

% Cell geometry

  Acap= 1.534e-4; % capacitive membrane area (cm^2)
  VNSR =(0.53*2.10e-6); % junctional SR volume (uL)
  Vmyo= 25.84e-6; % myoplasmic volume (uL)
  VJSR= ((double(Nclefts_FRU))*0.53*0.5*0.2*1.05e-10); % network SR volume (uL)
  VSS= (0.5*0.2*2.03e-12); % subspace volume (uL)
 

%% Main loop

while(time_now <= time_end)
   fprintf('%g \t %g \t %g \n',time_now, state(1), state(7));
   fileID = fopen(output_states_file, 'a');
   fprintf(fileID, '%d %s', time_now,' ');
   dlmwrite(output_states_file, state,'-append','delimiter',' ')
   fclose(fileID); 
   
    % Calculate the RyR state occupancies
    RyR_state_occupancies = zeros(1,8);
    for state_poss  = 1:8
        for iFRU = 1: NFRU_sim
           for icleft = 1: Nclefts_FRU
               for iRyR = 1:NRyRs_per_cleft
                       if(RyR_state(iFRU,icleft,iRyR) == state_poss)
                           RyR_state_occupancies(state_poss) = RyR_state_occupancies(state_poss) + 1;
                       end
               end
           end
        end
    end
    RyR_state_occupancies = RyR_state_occupancies./(NFRU_sim*Nclefts_FRU*NRyRs_per_cleft);
    fileID = fopen(RyR_occupancy_file, 'a');
    fprintf(fileID, '%d %s', time_now,' ');
    dlmwrite(RyR_occupancy_file ,RyR_state_occupancies,'-append','delimiter',' ')
    fclose(fileID);
    
   tstep = time_now + time_step;
   [oldstepsize, state, current, FRU_states, LType_state, RyR_state, Ito2_state] = rk54pd(time_now, tstep, oldstepsize, state,FRU_states, LType_state,Ito2_state, RyR_state, current);
   
   % Populate currents file and otherstates file
    [otherstates] = simulation_data( FRU_states, state, LType_state, RyR_state, Ito2_state);

    fileID = fopen(output_currents_file, 'a');
    fprintf(fileID,'%d %s',time_now,' ');
    dlmwrite(output_currents_file, current','-append','delimiter',' ')
    fclose(fileID);
    
    fileID = fopen(output_otherstates_file,'a');
    fprintf(fileID, '%d %s', time_now,' ');
    dlmwrite(output_otherstates_file,otherstates ,'-append','delimiter',' ')
    fclose(fileID);    
       
    time_now = tstep;
end
