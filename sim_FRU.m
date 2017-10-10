function [FRU_states_FRU, LType_state_FRU, RyR_state_FRU, Ito2_state_FRU] = sim_FRU(st_time,end_time, FRUdep_states0, FRUdep_statesf, FRU_states_FRU, LType_state_FRU, RyR_state_FRU, Ito2_state_FRU)

global Nclefts_FRU NRyRs_per_cleft 

NRVseqs_per_cleft = 2+ NRyRs_per_cleft+ 1;
NRVseqs_per_FRU = Nclefts_FRU *NRVseqs_per_cleft;
 O1_LType = 6;
 O2_LType = 12;
 Oy_LType = 2;
 O1_RyR   = 3;
 O2_RyR   = 4;
 O3_RyR   = 7;
 
 index_LCC_states = 1;
 index_LCC_Vinact = 2;
 
% Perform the monte carlo simulation

time_FRU = st_time;
timef = end_time;
FRUdep_states_interval = FRUdep_statesf - FRUdep_states0;
FRUdep_states = FRUdep_states0;
time_interval =  end_time - st_time;
LType_open = zeros(1, Nclefts_FRU);
NRyR_open  = zeros(1, Nclefts_FRU);

signal1 =zeros(1, Nclefts_FRU);

for icleft = 1:Nclefts_FRU
% set the indicator variables for channels that are open
 if ( ((LType_state_FRU(icleft,index_LCC_states)== O1_LType)||(LType_state_FRU(icleft,index_LCC_states)==O2_LType))... 
			&& (LType_state_FRU(icleft,index_LCC_Vinact)==Oy_LType))  
			LType_open(icleft) = 1.0;
 else 
			LType_open(icleft) = 0.0;
 end
 for i=1:NRyRs_per_cleft
     if(RyR_state_FRU(icleft,i)== O1_RyR ||RyR_state_FRU(icleft,i)== O2_RyR || RyR_state_FRU(icleft,i)== O3_RyR)
         NRyR_open(icleft) = NRyR_open(icleft) + 1;
     end
 end
 signal1(icleft) = 0;
end

% CaSSclamp_hold = 1.099940000000000e-04;
% CaSSclamp_set = 0.003129320163158;
% CaSSclamp_shift = 10;
% CaSSclamp_duration = 200;
% CaSSclamp_period = 1000;
while(time_FRU<end_time)
%     if(mod(time_FRU, CaSSclamp_period) >= CaSSclamp_shift && mod(time_FRU, CaSSclamp_period) <= (CaSSclamp_duration + CaSSclamp_shift))
%         FRU_states_FRU(2) = CaSSclamp_set;
%         FRU_states_FRU(3) = CaSSclamp_set;
%         FRU_states_FRU(4) = CaSSclamp_set;
%         FRU_states_FRU(5) = CaSSclamp_set;
%     else       
%         FRU_states_FRU(2) = CaSSclamp_hold;
%         FRU_states_FRU(3) = CaSSclamp_hold;
%         FRU_states_FRU(4) = CaSSclamp_hold;
%         FRU_states_FRU(5) = CaSSclamp_hold;
%     end
    [max_exit_rate,LType_trans_rates,LType_index,LType_length, LType_Vinact_exitrate,...
     RyR_trans_rates,RyR_index,RyR_length,Ito2_exitrate ] = fru_rates_local(LType_state_FRU, RyR_state_FRU, Ito2_state_FRU, FRUdep_states, FRU_states_FRU);
    
    time_stepFRU = 0.1/max_exit_rate;
    time_stepFRU = min(time_stepFRU, 10.001e-3);
    CaSStemp = max(FRU_states_FRU(2:end));
    if(CaSStemp > 1.5e-3)
        time_stepFRU = min(time_stepFRU, 5.001e-3);
    else
        time_stepFRU = min(time_stepFRU, timef-time_FRU);
    end
    
   FRU_states1 = FRU_states_FRU;  
   [dFRU_states] = fcn_fru(time_FRU,FRU_states1,FRUdep_states,LType_open,NRyR_open);
   k1 = time_stepFRU.*dFRU_states;
   y_1 = FRU_states1 + k1;
   time_FRU = time_FRU + time_stepFRU;
    % Extrapolate the FRU dependent variables
 
   FRUdep_states = FRUdep_states0 + ((time_FRU-st_time)/time_interval).*FRUdep_states_interval;
   [dFRU_states] = fcn_fru(time_FRU,y_1,FRUdep_states,LType_open,NRyR_open);    
   FRU_states_FRU = FRU_states_FRU + (k1 + time_stepFRU.*dFRU_states)./2.0;
    
    
    LType_trans_probs = zeros(Nclefts_FRU,4);
    LType_Vinact_exitprob = zeros(Nclefts_FRU,1);
    Ito2_exitprob = zeros(Nclefts_FRU,1);
    RyR_trans_probs = zeros(Nclefts_FRU, NRyRs_per_cleft,4);
    
    for icleft = 1:Nclefts_FRU
        LType_trans_probs(icleft,1) = 1.0 - time_stepFRU*LType_trans_rates(icleft,1);
        LType_trans_probs(icleft,2) = time_stepFRU*LType_trans_rates(icleft,2);
        for i = 3:LType_length(icleft)
            LType_trans_probs(icleft,i) = time_stepFRU*LType_trans_rates(icleft,i);
        end
        LType_Vinact_exitprob(icleft) = 1.0 - time_stepFRU*LType_Vinact_exitrate(icleft);
        Ito2_exitprob(icleft) = 1.0 - time_stepFRU*Ito2_exitrate(icleft);
        for j = 1:NRyRs_per_cleft
            RyR_trans_probs(icleft,j,1) = 1.0 - time_stepFRU*RyR_trans_rates(icleft,j,1);
			RyR_trans_probs(icleft,j,2) = time_stepFRU*RyR_trans_rates(icleft,j,2);
                for i= 3: RyR_length(icleft,j)
                        RyR_trans_probs(icleft,j,i) = time_stepFRU*RyR_trans_rates(icleft,j,i);
                end
        end
    end
    
    % Determine if transitions took place if any
    unidev = rand(NRVseqs_per_FRU,1);
    for icleft = 1:Nclefts_FRU
		base_ind= (icleft-1)*NRVseqs_per_cleft + 1;
		Accum_Prob = LType_trans_probs(icleft,1);
            if (unidev(base_ind) > Accum_Prob) 
                    count = 2; % We know that state has been changed
                    Accum_Prob = Accum_Prob + LType_trans_probs(icleft,2);
                        while (unidev(base_ind)>Accum_Prob) 
                            count = count + 1;
                            Accum_Prob = Accum_Prob + LType_trans_probs(icleft,count);
                        end
                    LType_state_FRU(icleft,index_LCC_states) = LType_index(icleft,count);
            end
            
            if (unidev(1+base_ind)>=LType_Vinact_exitprob(icleft)) 
                    LType_state_FRU(icleft,index_LCC_Vinact) = 3 - LType_state_FRU(icleft,index_LCC_Vinact);
            end
            
            for i=1:NRyRs_per_cleft
                Accum_Prob = RyR_trans_probs(icleft,i,1);
                    if (unidev(i+1+base_ind)>Accum_Prob) 
                            count = 2; % we know that state has been changed
                            Accum_Prob = Accum_Prob + RyR_trans_probs(icleft,i,2);
                                while (unidev(i+1+base_ind)>Accum_Prob) 
                                    count = count + 1;
                                    Accum_Prob = Accum_Prob + RyR_trans_probs(icleft,i,count);
                                end
                            RyR_state_FRU(icleft,i) = RyR_index(icleft,i,count);
                            signal1(icleft) = 1;
                    end
            end 
            if (unidev((icleft)*NRVseqs_per_cleft)>=Ito2_exitprob(icleft)) 
                    Ito2_state_FRU(icleft) = 3 - Ito2_state_FRU(icleft);
            end
    end 
    % Reset indicators for which channels are open
    LType_open = zeros(1, Nclefts_FRU);
    NRyR_open = zeros(1, Nclefts_FRU);
    for icleft = 1:Nclefts_FRU 
     if ( ((LType_state_FRU(icleft,index_LCC_states)== O1_LType)||(LType_state_FRU(icleft,index_LCC_states)==O2_LType))... 
			&& (LType_state_FRU(icleft,index_LCC_Vinact)==Oy_LType))  
                LType_open(icleft) = 1.0;
     else   
                LType_open(icleft) = 0.0;
     end

     for i=1:NRyRs_per_cleft
         if(RyR_state_FRU(icleft,i)== O1_RyR ||RyR_state_FRU(icleft,i)== O2_RyR || RyR_state_FRU(icleft,i)== O3_RyR)
             NRyR_open(icleft) = NRyR_open(icleft) + 1;
         end
     end
    end
end

end

