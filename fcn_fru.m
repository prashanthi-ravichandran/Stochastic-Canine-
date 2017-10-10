function [ dFRU_states1 ] = fcn_fru(time_FRU, FRU_states1, FRUdep_states,LType_open,NRyR_open)

global index_frudep_V index_frudep_Cai index_frudep_CaNSR 
global Cao
global Faraday  RT_over_F 
global VJSR VSS 

V = FRUdep_states(index_frudep_V);
Cai = FRUdep_states(index_frudep_Cai);
CaNSR = FRUdep_states(index_frudep_CaNSR);
PCa= (1.5/2.8*0.2*0.9*0.9468e-11); %(cm/s) *uF  

JDHconstant=1.0/(2.0*VSS*Faraday*1000.0);
% JSR to subspace through a single RyR (1/ms)
JRyRmax= (0.6*1.5*1.1*3.96);
% NSR to JSR (1/ms)
tautr= 3.0;
% subspace to cytosol (1/ms)
tauxfer= 0.005;
% intersubspace transfer rate (1/ms)
tauss2ss= (10.0*0.005);

BSLtot=   1.124;  % (mM)
CSQNtot= 13.5; 
BSRtot=   0.047; 
KBSL=     0.0087; 
KmCSQN=   0.63; 
KBSR=     0.00087; 

CaJSR = FRU_states1(1);
CaSS_1 = FRU_states1(2);
CaSS_2 = FRU_states1(3);
CaSS_3 = FRU_states1(4);
CaSS_4 = FRU_states1(5);

Jtr = (CaNSR - CaJSR)/tautr;

JRyR_1 = JRyRmax*(NRyR_open(1))*(CaJSR-CaSS_1);
JRyR_2 = JRyRmax*(NRyR_open(2))*(CaJSR-CaSS_2);
JRyR_3 = JRyRmax*(NRyR_open(2))*(CaJSR-CaSS_3);
JRyR_4 = JRyRmax*(NRyR_open(4))*(CaJSR-CaSS_4);
JRyRtot = JRyR_1+JRyR_2+JRyR_3+JRyR_4;

Jxfer_1 = (CaSS_1-Cai)/tauxfer;
Jxfer_2 = (CaSS_2-Cai)/tauxfer;
Jxfer_3 = (CaSS_3-Cai)/tauxfer;
Jxfer_4 = (CaSS_4-Cai)/tauxfer;

Jss2ss_1 = (CaSS_1 - CaSS_2)/tauss2ss;
Jss2ss_2 = (CaSS_2 - CaSS_3)/tauss2ss;
Jss2ss_3 = (CaSS_3 - CaSS_4)/tauss2ss;
Jss2ss_4 = (CaSS_4 - CaSS_1)/tauss2ss;

VF_over_RT=V/RT_over_F;
VFsq_over_RT=(1000.0*Faraday)*VF_over_RT;
exp_VFRT = exp(2.0*VF_over_RT);

if (abs(V)<1.e-6) 
  JDHPR_1 = -PCa*2.0*1000.0*Faraday*(CaSS_1-Cao*0.341)*LType_open(1)*JDHconstant;
  JDHPR_2 = -PCa*2.0*1000.0*Faraday*(CaSS_2-Cao*0.341)*LType_open(2)*JDHconstant;
  JDHPR_3 = -PCa*2.0*1000.0*Faraday*(CaSS_3-Cao*0.341)*LType_open(2)*JDHconstant;
  JDHPR_4 = -PCa*2.0*1000.0*Faraday*(CaSS_4-Cao*0.341)*LType_open(4)*JDHconstant;
else 
  a2 = PCa*4.0*VFsq_over_RT/(exp_VFRT-1.0)*JDHconstant;
  JDHPR_1 = -(CaSS_1*exp_VFRT-Cao*0.341)*a2*LType_open(1);
  JDHPR_2 = -(CaSS_2*exp_VFRT-Cao*0.341)*a2*LType_open(2);
  JDHPR_3 = -(CaSS_3*exp_VFRT-Cao*0.341)*a2*LType_open(2);
  JDHPR_4 = -(CaSS_4*exp_VFRT-Cao*0.341)*a2*LType_open(4);
end


a1 = BSRtot*KBSR/((CaSS_1+KBSR)*(CaSS_1+KBSR));
a2 = BSLtot*KBSL/((CaSS_1+KBSL)*(CaSS_1+KBSL));
beta_SS_1 = 1.0/(1.0+a1+a2) ;
a1 = BSRtot*KBSR/((CaSS_2+KBSR)*(CaSS_2+KBSR));
a2 = BSLtot*KBSL/((CaSS_2+KBSL)*(CaSS_2+KBSL));
beta_SS_2 = 1.0/(1.0+a1+a2) ;
a1 = BSRtot*KBSR/((CaSS_3+KBSR)*(CaSS_3+KBSR));
a2 = BSLtot*KBSL/((CaSS_3+KBSL)*(CaSS_3+KBSL));
beta_SS_3 = 1.0/(1.0+a1+a2) ;
a1 = BSRtot*KBSR/((CaSS_4+KBSR)*(CaSS_4+KBSR));
a2 = BSLtot*KBSL/((CaSS_4+KBSL)*(CaSS_4+KBSL));
beta_SS_4 = 1.0/(1.0+a1+a2) ;

a1 = CSQNtot*KmCSQN/((CaJSR+KmCSQN)*(CaJSR+KmCSQN));
beta_JSR = 1.0/(1.0+a1);

dFRU_states1 = zeros(5,1);
dFRU_states1(1) = beta_JSR*(Jtr - VSS/VJSR*JRyRtot); % dCaJSR
dFRU_states1(2) = beta_SS_1*(JDHPR_1 + JRyR_1 - Jxfer_1 - Jss2ss_1 + Jss2ss_4); % dCaSS1
dFRU_states1(3) = beta_SS_2*(JDHPR_2 + JRyR_2 - Jxfer_2 - Jss2ss_2 + Jss2ss_1); % dCaSS2
dFRU_states1(4) = beta_SS_3*(JDHPR_3 + JRyR_3 - Jxfer_3 - Jss2ss_3 + Jss2ss_2); % dCaSS3
dFRU_states1(5) = beta_SS_4*(JDHPR_4 + JRyR_4 - Jxfer_4 - Jss2ss_4 + Jss2ss_3); % dCaSS4
end

