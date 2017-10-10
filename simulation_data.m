function [ otherstates ] = simulation_data( FRU_state, state, Ltype_state, RyR_state, Ito2_state)

 global NFRU_sim Nclefts_FRU NRyRs_per_cleft NFRU_scale

  VNSR =(0.53*2.10e-6); % junctional SR volume (uL)
  Vmyo= 25.84e-6; % myoplasmic volume (uL)
  VJSR= ((double(Nclefts_FRU))*0.53*0.5*0.2*1.05e-10); % network SR volume (uL)
  VSS= (0.5*0.2*2.03e-12); % subspace volume (uL)
  
    LTRPNtot= 70.0e-3;
	% total troponin high affinity site conc. (mM)
  HTRPNtot= 140.0e-3;
	% total myoplasmic calmodulin concentration (mM)
  CMDNtot= 50.0e-3;
	% total myoplasmic EGTA concentration (mM)
  EGTAtot= 0.0;
	% Ca++ half sat. constant for calmodulin (mM)
  KmCMDN= 2.38e-3;
	% Ca++ half sat. constant for EGTA (mM)
  KmEGTA= 1.5e-4;
  
 index_CaSSavg  = 1;
 index_CaJSRavg = 2;
 index_JRyRtot  = 3;
 index_PRyR_Open = 4;
 index_PRyR_ready = 5;
 index_PNorm_Mode = 6;
 index_PnotVinact = 7;
 index_PLtype_Open = 8;
 index_CaTOT2 = 9;
 index_PIto2_Open = 10;
 index_CaJSRtot = 11;
 index_CaSStot = 12;
 
 index_Cai       = 7;
 index_CaNSR     = 8;
 index_LTRPNCa   = 10;
 index_HTRPNCa   = 11;
 
 O1_Ltype = 6;
 O2_Ltype = 12;
 Oy_Ltype = 2;
 O1_RyR   = 3;
 O2_RyR   = 4;
 O3_RyR   = 7;
 O_Ito2   = 2;
   
  CaSS = 0.0;
  CaJSR = 0.0;
  CaTOT_SS = 0.0;
  CaTOT_JSR = 0.0;
  JRyR = 0.0;
  PNorm_Mode = 0.0;
  PRyR_Open = 0.0;
  PnotVinact = 0.0;
  PLtype_Open = 0.0;
  PIto2_Open = 0.0;
  NRyR_ready=0;
  
  BSLtot=   1.124;  % (mM)
  CSQNtot= 13.5; 
  BSRtot=   0.047; 
  KBSL=     0.0087; 
  KmCSQN=   0.63; 
  KBSR=     0.00087; 
  JRyRmax= (0.6*1.5*1.1*3.96);
   for iFRU = 1: NFRU_sim
    CaJSR = CaJSR + FRU_state(iFRU,1);
    CaTOT_JSR = CaTOT_JSR + FRU_state(iFRU,1)*(1.0+CSQNtot/(KmCSQN+FRU_state(iFRU,1)));
    for icleft = 1: Nclefts_FRU
      CaSS  = CaSS + FRU_state(iFRU, icleft+1);
      CaTOT_SS = CaTOT_SS + FRU_state(iFRU,icleft+1)*(1.0+BSRtot/(KBSR+FRU_state(iFRU,icleft+1))+BSLtot/(KBSL+FRU_state(iFRU,icleft+1)));
      NRyR_open=0;
      for i= 1:NRyRs_per_cleft
        if ((RyR_state(iFRU,icleft,i)==O1_RyR)||(RyR_state(iFRU,icleft,i)==O2_RyR) ||(RyR_state(iFRU,icleft,i)==O3_RyR))     	  
            NRyR_open = NRyR_open + 1;
        end
        if ((RyR_state(iFRU,icleft,i)==1)||(RyR_state(iFRU,icleft,i)==2)) 
	  	    NRyR_ready = NRyR_ready + 1;
        end
      end
      JRyR = JRyR + JRyRmax*NRyR_open*(FRU_state(iFRU,1) -FRU_state(iFRU,icleft+1));
      PRyR_Open = PRyR_Open + NRyR_open;
      if (Ltype_state(iFRU,icleft,1) <=6) 
		  PNorm_Mode = PNorm_Mode + 1.0;
      end
      if (Ltype_state(iFRU, icleft,2) == Oy_Ltype) 
        PnotVinact = PnotVinact + 1.0;
      end
      if (((Ltype_state(iFRU, icleft,1)==O1_Ltype)||(Ltype_state(iFRU,icleft,1) ==O2_Ltype))&&(Ltype_state(iFRU,icleft,2)==Oy_Ltype))
    	  PLtype_Open = PLtype_Open + 1.0;
      end
      if (Ito2_state(iFRU,icleft)==O_Ito2)
        PIto2_Open = PIto2_Open + 1.0;
      end
    end
   end
 
  otherstates(index_CaSSavg)  = CaSS/(NFRU_sim*Nclefts_FRU);
  otherstates(index_CaJSRavg) = CaJSR/NFRU_sim;
  otherstates(index_JRyRtot) = JRyR*NFRU_scale*VSS/Vmyo;
  otherstates(index_PRyR_Open)  = PRyR_Open/(NRyRs_per_cleft*Nclefts_FRU*NFRU_sim);
  otherstates(index_PRyR_ready)  = (NRyR_ready)/(NRyRs_per_cleft*Nclefts_FRU*NFRU_sim);
  otherstates(index_CaTOT2) = 1.e6*((CaTOT_SS*VSS + CaTOT_JSR*VJSR)*NFRU_scale  ...
				  + state(index_CaNSR)*VNSR  ...
				  + state(index_Cai)*Vmyo* (1.0 + CMDNtot/(KmCMDN+state(index_Cai)) + EGTAtot/(KmEGTA+state(index_Cai))) ...
				  + (state(index_LTRPNCa)*LTRPNtot + state(index_HTRPNCa)*HTRPNtot)*Vmyo );  %femtomoles 
  otherstates(index_PNorm_Mode) = PNorm_Mode/(NFRU_sim*Nclefts_FRU);
  otherstates(index_PnotVinact) = PnotVinact/(NFRU_sim*Nclefts_FRU);
  otherstates(index_PLtype_Open) = PLtype_Open/(NFRU_sim*Nclefts_FRU);
  otherstates(index_PIto2_Open) = PIto2_Open/(NFRU_sim*Nclefts_FRU);
  otherstates(index_CaJSRtot) = 1.e6*CaTOT_SS*VSS*NFRU_scale;
  otherstates(index_CaSStot) = 1.e6*CaTOT_JSR*VJSR*NFRU_scale;

end

