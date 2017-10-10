function [Jxfer,Jtr,ICa,Ito2] = send_calc_fru_flux(FRUdep_states, LType_state,Ito2_state,FRU_states)

global index_frudep_V index_frudep_Cai index_frudep_CaNSR
global NFRU_sim Nclefts_FRU NFRU_scale
global O1_LType O2_LType O_Ito2 Oy_LType

O1_LType = 6;
O2_LType = 12;
O_Ito2 = 2;
Oy_LType = 2;

Cao=   2.0; %extracellular Ca++ concentration (mM)
Clo= 150.0;   % extracellular Cl-  concentration (mM)
Cli=  20.0;   % intracellular Cl-  concentration (mM)
  
PCa= (1.5/2.8*0.2*0.9*0.9468e-11); %(cm/s) *uF  
Acap= 1.534e-4; % capacitive membrane area (cm^2)
PCl= 2.65e-15; %(cm/s) *uF 
  
Faraday=  96.5;     % Faraday's constant (C/mmol)
Temp=    310.0;     % absolute temperature (K)
Rgas=      8.314;   % ideal gas constant (J/[mol*K])
RT_over_F= (Rgas*Temp/Faraday); %  Rgas*Temp/Faraday (mV)


V = FRUdep_states(index_frudep_V);
Cai = FRUdep_states(index_frudep_Cai);
CaNSR = FRUdep_states(index_frudep_CaNSR);

VF_over_RT=V/RT_over_F;
VFsq_over_RT=(1000.0*Faraday)*VF_over_RT;
exp_VFRT = exp(VF_over_RT);
exp_2VFRT = exp_VFRT*exp_VFRT;

ICa_numerator = 0.0;
OCa_numerator = 0;
NIto2_Open = 0;

index_LCC_states = 1;
index_LCC_Vinact = 2;
for iFRU = 1:NFRU_sim
    for icleft = 1:Nclefts_FRU
        if((LType_state(iFRU, icleft,index_LCC_states)== O1_LType ||LType_state(iFRU, icleft,index_LCC_states)== O2_LType) && ...
                (LType_state(iFRU, icleft,index_LCC_Vinact)== Oy_LType))        
            ICa_numerator = ICa_numerator + FRU_states(iFRU,icleft+1);
            OCa_numerator = OCa_numerator + 1;
        end
        if(Ito2_state(iFRU,icleft)==O_Ito2)
            NIto2_Open = NIto2_Open + 1;
        end
    end
end

sum_CaSS = sum(sum(FRU_states(:,2:end)));
sum_CaJSR = sum(FRU_states(:,1));

if (abs(V)<1.e-6)  % First order Taylor expansion
    ICa_fac    =ICa_numerator- OCa_numerator*Cao*0.341; 
    ICa_local  = PCa*2.0*1000.0*Faraday*ICa_fac;
    ICa_local = ICa_local/Acap;		% divide by uF(Acap) to get current normalized to surface area
    Ito2_local = (double(NIto2_Open))*PCl*1000.0*Faraday*(Clo-Cli);
    Ito2_local = Ito2_local/Acap;	% divide by uF(Acap) to get current normalized to surface area
else 
    ICa_fac   =ICa_numerator*exp_2VFRT- OCa_numerator*Cao*0.341; 
    ICa_local = PCa*4.0*VFsq_over_RT*ICa_fac/(exp_2VFRT-1.0); 
    ICa_local = ICa_local/Acap;		% divide by uF(Acap) to get current normalized to surface area
    Ito2_local = (double(NIto2_Open))*PCl*VFsq_over_RT*(Cli-Clo*exp_VFRT)/(1.0 - exp_VFRT);
    Ito2_local = Ito2_local/Acap;	% divide by uF(Acap) to get current normalized to surface area
end
    % NSR to JSR (1/ms)
  tautr= 3.0;
	% subspace to cytosol (1/ms)
  tauxfer= 0.005;
  Jtr_local = ((double(NFRU_sim))*CaNSR - sum_CaJSR)/tautr;
  Jxfer_local = (sum_CaSS - (NFRU_sim*Nclefts_FRU)*Cai)/tauxfer;

  Jxfer=NFRU_scale*Jxfer_local;
  Jtr=NFRU_scale*Jtr_local;
  ICa=NFRU_scale*ICa_local;
  Ito2=NFRU_scale*Ito2_local;

end

