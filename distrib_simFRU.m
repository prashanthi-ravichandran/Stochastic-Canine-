function [Jxfer,Jtr,ICa,Ito2,FRU_states, LType_state, RyR_state, Ito2_state] = distrib_simFRU(st_time,end_time, FRUdep_states0, FRUdep_statesf,FRU_states, LType_state, RyR_state, Ito2_state)

global NFRU_sim Nclefts_FRU NRyRs_per_cleft Nindepstates_LType
LType_state_FRU = zeros(Nclefts_FRU, Nindepstates_LType);
RyR_state_FRU = zeros(Nclefts_FRU, NRyRs_per_cleft);
Ito2_state_FRU = zeros(Nclefts_FRU,1);
FRU_states_FRU = zeros((Nclefts_FRU + 1),1);
% Perform the stochastic simulation
for i = 1:NFRU_sim
    for icleft = 1:Nclefts_FRU
        LType_state_FRU(icleft,1)     = LType_state(i,icleft,1);
        LType_state_FRU(icleft,2)     = LType_state(i,icleft,2);
        Ito2_state_FRU(icleft,1)      = Ito2_state(i,icleft);
        for iRyR = 1:NRyRs_per_cleft
            RyR_state_FRU(icleft, iRyR)      = RyR_state(i,icleft,iRyR);
        end
    end
    for icleft = 1:(Nclefts_FRU + 1)
        FRU_states_FRU(icleft) = FRU_states(i, icleft);
    end

[FRU_states_FRU, LType_state_FRU, RyR_state_FRU, Ito2_state_FRU] = sim_FRU(st_time,end_time, FRUdep_states0, FRUdep_statesf,...
                                                                               FRU_states_FRU, LType_state_FRU, RyR_state_FRU, Ito2_state_FRU);
LType_state(i,:,:) = LType_state_FRU;
RyR_state(i,:,:)   = RyR_state_FRU;
Ito2_state(i,:,:)  = Ito2_state_FRU;  
FRU_states(i,:) = FRU_states_FRU;
end

% Calculate the macroscopic fluxes

[Jxfer,Jtr,ICa,Ito2] = send_calc_fru_flux(FRUdep_statesf, LType_state,Ito2_state,FRU_states);

end

