function [oldstepsize, state, current, FRU_states, LType_state, RyR_state, Ito2_state] = rk54pd(time_now, tstep, oldstepsize, state,FRU_states, LType_state,Ito2_state, RyR_state, current)

global index_frudep_V index_frudep_Cai index_frudep_CaNSR
global step_min step_max tolrk N
global period shift APClampflag 
global index_V index_mNa index_jNa index_hNa index_Nai index_Ki index_Cai index_CaNSR index_xKs
global index_LTRPNCa index_HTRPNCa index_C0Kv43 index_C1Kv43 index_C2Kv43 index_C3Kv43 index_OKv43 
global index_CI0Kv43 index_CI1Kv43 index_CI2Kv43 index_CI3Kv43 index_OIKv43 index_C0Kv14 index_C1Kv14
global index_C2Kv14 index_C3Kv14 index_OKv14 index_CI0Kv14 index_CI1Kv14 index_CI2Kv14 index_CI3Kv14 index_OIKv14 
global index_CaTOT index_C1Herg index_C2Herg index_C3Herg index_OHerg index_IHerg
%%
index_frudep_V = 1;
index_frudep_Cai = 2;
index_frudep_CaNSR = 3;

errweight(index_V) = 0; %1.e-2; % not independent
errweight(index_mNa) = 1.0;
errweight(index_hNa) = 1.0;
errweight(index_jNa) = 1.0;
errweight(index_Nai) = 0.1;
errweight(index_Ki) = 1/140.0;
errweight(index_Cai) = 1.e3;    
errweight(index_CaNSR) = 2.5; 
errweight(index_xKs) = 1.0;    
errweight(index_LTRPNCa) = 1.0;
errweight(index_HTRPNCa) = 1.0;
errweight(index_C0Kv43) = 1.0;
errweight(index_C1Kv43) = 1.0;
errweight(index_C2Kv43) = 1.0;
errweight(index_C3Kv43) = 1.0;
errweight(index_OKv43) = 1.0;
errweight(index_CI0Kv43) = 1.0;
errweight(index_CI1Kv43) = 1.0;
errweight(index_CI2Kv43) = 1.0;
errweight(index_CI3Kv43) = 1.0;
errweight(index_OIKv43) = 0.0; % 1.0, not independent
errweight(index_C0Kv14) = 1.0;
errweight(index_C1Kv14) = 1.0;
errweight(index_C2Kv14) = 1.0;
errweight(index_C3Kv14) = 1.0;
errweight(index_OKv14) = 1.0;
errweight(index_CI0Kv14) = 1.0;
errweight(index_CI1Kv14) = 1.0;
errweight(index_CI2Kv14) = 1.0;
errweight(index_CI3Kv14) = 1.0;
errweight(index_OIKv14) = 0.0; % 1.0, not independent
errweight(index_CaTOT) = 0.1;
errweight(index_C1Herg)= 1.0;
errweight(index_C2Herg)= 1.0;
errweight(index_C3Herg)= 1.0;
errweight(index_OHerg)= 1.0;
errweight(index_IHerg)= 0.0; % 1.0, not independent

step_size = min(oldstepsize, tstep - time_now);
step_size = min(step_max, step_size);
oldstepsize = step_size;
start_time = time_now;
notdone = 1;
success = 1;
keepc = 1;
 
DepFlag = 1;
stepsno = 0;
notaccepted = 0;
forcedaccept = 0;
huge_number = 1.e99;
FRUdep_states0 = zeros(3,1);
FRUdep_statesf = zeros(3,1);
Ncur = 24;
current = zeros(Ncur,1);
%%
while(notdone)
    stepsno = stepsno + 1;
    FRUdep_states0(index_frudep_V) = state(index_V);
    FRUdep_states0(index_frudep_Cai) = state(index_Cai);
    FRUdep_states0(index_frudep_CaNSR) = state(index_CaNSR);
    if(success)
        [Jxfer,Jtr,ICa,Ito2] = send_calc_fru_flux(FRUdep_states0, LType_state,Ito2_state,FRU_states);
        [F1,state,current] = fcn(start_time, state, current, keepc, Jxfer, Jtr, ICa, Ito2);
        keepc = 0;
        [FRU_states_hold,LType_state_hold,Ito2_state_hold,RyR_state_hold] = save_state(FRU_states, LType_state, Ito2_state, RyR_state);
    end       
    k1  = step_size.*F1;
    y_1 = state + (k1./5.0);
    FRUdep_statesf(index_frudep_V) = y_1(index_V);
    FRUdep_statesf(index_frudep_Cai) = y_1(index_Cai);
    FRUdep_statesf(index_frudep_CaNSR) = y_1(index_CaNSR);
    end_time = start_time + step_size / 5.0;
    
    %Perform the stochastic simulation and collect the macroscopic fluxes
    %at t + step_size/5.0
    [Jxfer,Jtr,ICa,Ito2,FRU_states, LType_state, RyR_state, Ito2_state] = distrib_simFRU(start_time,end_time, FRUdep_states0, FRUdep_statesf,FRU_states, LType_state, RyR_state, Ito2_state);
    
    FRUdep_states0(index_frudep_V) = state(index_V);
    FRUdep_states0(index_frudep_Cai) = state(index_Cai);
    FRUdep_states0(index_frudep_CaNSR) = state(index_CaNSR);
    
    time = start_time + step_size/5.0;
    [F,y_1,current] = fcn(time, y_1, current, keepc, Jxfer, Jtr, ICa, Ito2);
    k2 = step_size.*F;
    y_1 = state + ( (3.0/40.0).*k1 + (9.0/40.0).*k2);
    
    FRUdep_statesf(index_frudep_V) = y_1(index_V);
    FRUdep_statesf(index_frudep_Cai) = y_1(index_Cai);
    FRUdep_statesf(index_frudep_CaNSR) = y_1(index_CaNSR);
    
    st_time = start_time + step_size / 5.0;
    end_time = start_time + step_size*(3.0/10.0);
    
    %Perform the stochastic simulation and collect the macroscopic fluxes
    %at t + step_size*(3.0/10.0)
    [Jxfer,Jtr,ICa,Ito2,FRU_states, LType_state, RyR_state, Ito2_state] = distrib_simFRU(st_time,end_time, FRUdep_states0, FRUdep_statesf,FRU_states, LType_state, RyR_state, Ito2_state);
    
    FRUdep_states0(index_frudep_V) = state(index_V);
    FRUdep_states0(index_frudep_Cai) = state(index_Cai);
    FRUdep_states0(index_frudep_CaNSR) = state(index_CaNSR);
    
    time = start_time + step_size*(3.0/10.0);
    [F,y_1,current] = fcn(time, y_1, current, keepc, Jxfer, Jtr, ICa, Ito2);
    k3 = step_size.*F;
    y_1 = state + ( (3.0/10.0).*k1 - (9.0/10.0).*k2 + (6.0/5.0)*k3);
    
    FRUdep_statesf(index_frudep_V) = y_1(index_V);
    FRUdep_statesf(index_frudep_Cai) = y_1(index_Cai);
    FRUdep_statesf(index_frudep_CaNSR) = y_1(index_CaNSR);
    
    st_time = start_time + step_size*(3.0/10.0);
    end_time = start_time + step_size*(3.0/5.0);
    
    %Perform the stochastic simulation and collect the macroscopic fluxes
    %at t + step_size*(3.0/10.0)
    [Jxfer,Jtr,ICa,Ito2,FRU_states, LType_state, RyR_state, Ito2_state] = distrib_simFRU(st_time,end_time, FRUdep_states0, FRUdep_statesf,FRU_states, LType_state, RyR_state, Ito2_state);
    
    FRUdep_states0(index_frudep_V) = state(index_V);
    FRUdep_states0(index_frudep_Cai) = state(index_Cai);
    FRUdep_states0(index_frudep_CaNSR) = state(index_CaNSR);
    
    time = start_time + step_size*(3.0/5.0);
    [F,y_1,current] = fcn(time, y_1, current, keepc, Jxfer, Jtr, ICa, Ito2);
    k4 = step_size.*F;
    y_1 = state + ( (226.0/729.0).*k1 - (25.0/27.0).*k2 + (880.0/729.0).*k3 + (55.0/729.0).*k4);
    
    FRUdep_statesf(index_frudep_V) = y_1(index_V);
    FRUdep_statesf(index_frudep_Cai) = y_1(index_Cai);
    FRUdep_statesf(index_frudep_CaNSR) = y_1(index_CaNSR);
    
    st_time = start_time + step_size*(3.0/5.0);
    end_time = start_time + step_size*(2.0/3.0);
    
    %Perform the stochastic simulation and collect the macroscopic fluxes
    %at t + step_size*(3.0/10.0)
    [Jxfer,Jtr,ICa,Ito2,FRU_states, LType_state, RyR_state, Ito2_state] = distrib_simFRU(st_time,end_time, FRUdep_states0, FRUdep_statesf,FRU_states, LType_state, RyR_state, Ito2_state);
    
    FRUdep_states0(index_frudep_V) = state(index_V);
    FRUdep_states0(index_frudep_Cai) = state(index_Cai);
    FRUdep_states0(index_frudep_CaNSR) = state(index_CaNSR);

    time = start_time + step_size*(2.0/3.0);
    [F,y_1,current] = fcn(time, y_1, current, keepc, Jxfer, Jtr, ICa, Ito2);
    k5 = step_size.*F;
    y_1 = state -(181.0/270.0).*k1 + (5.0/2.0).*k2 - (266.0/297.0).*k3 - (91.0/27.0).*k4 + (189/55.0).*k5;
    
    FRUdep_statesf(index_frudep_V) = y_1(index_V);
    FRUdep_statesf(index_frudep_Cai) = y_1(index_Cai);
    FRUdep_statesf(index_frudep_CaNSR) = y_1(index_CaNSR);
    
    st_time = start_time + step_size*(2.0/3.0);
    end_time = start_time + step_size;
    
    %Perform the stochastic simulation and collect the macroscopic fluxes
    %at t + step_size*(3.0/10.0)
    [Jxfer,Jtr,ICa,Ito2,FRU_states, LType_state, RyR_state, Ito2_state] = distrib_simFRU(st_time,end_time, FRUdep_states0, FRUdep_statesf,FRU_states, LType_state, RyR_state, Ito2_state);
    
    FRUdep_states0(index_frudep_V) = state(index_V);
    FRUdep_states0(index_frudep_Cai) = state(index_Cai);
    FRUdep_states0(index_frudep_CaNSR) = state(index_CaNSR);

    time = start_time + step_size;
    [F,y_1,current] = fcn(time, y_1, current, keepc, Jxfer, Jtr, ICa, Ito2);

    k6 = step_size.*F;

    tr_error = 0.0;
    min_tolvalue = 1e-10;
    y_1 = state + (31.0 / 540.0).* k1 + (190.0 / 297.0).* k3 - (145.0 / 108.0).* k4 + (351.0/220.0).* k5 + (1.0 / 20.0).* k6;
    ym = state + (19.0/216.0).*k1 + (1000.0/2079.0).*k3 - (125.0/216.0).*k4 + (81.0/88.0).*k5 + (5.0/56.0).* k6;
    errtmp = abs((ym - y_1).* 0.2.*errweight);
    tr_error = min(max(errtmp),huge_number);
    if(tr_error < tolrk)
        for i = 1:37
            if(abs(y_1(i)) < min_tolvalue)
                y_1(i) = 0;
            end
        end
        state = y_1;
        start_time = start_time + step_size;
        success = 1;
        if(start_time >= tstep)
            notdone = 0;
        else
            step_size = min(0.85*step_size*(tolrk/tr_error)^0.2, 4.0 * step_size);
            notdone = 1;
            oldstepsize = step_size;
        end
        time_to_stim = shift - mod(start_time, period);
        if (time_to_stim < 0.0) 
            time_to_stim = time_to_stim + period;
        end
        if (time_to_stim < step_size) 
          step_size = max(step_min, time_to_stim - 0.5 * step_min);
        end
    else
        if(step_size <= step_min)
            for i = 1:N
                if(abs(y_1(i)) < min_tolvalue)
                    y_1(i) = 0;
                end
            end
            state = y_1;
            start_time = start_time + stap_size;
            forcedaccept = forcedaccept + 1;
            success = 1;
            if(start_time >= tstep)
                 notdone = 0;
            else
                step_size = step_min;
                notdone = 1;           
            end
            oldstepsize = step_size;
        else
            step_size = min(0.85 * step_size * (tolrk/tr_error)^0.2, 4.0 * step_size);
            if(APClampflag)
                time_to_stim = shift - mod(start_time, period);
                if (time_to_stim < 0.0) 
                    time_to_stim = time_to_stim + period;
                end
                if (time_to_stim < step_size) 
                  step_size = max(step_min, time_to_stim - 0.5 * step_min);
                end
            end
            success = 0;
            notdone = 1;
            [FRU_states, LType_state, Ito2_state, RyR_state] = resume_state(FRU_states_hold,LType_state_hold,Ito2_state_hold,RyR_state_hold);
        end
end

end

