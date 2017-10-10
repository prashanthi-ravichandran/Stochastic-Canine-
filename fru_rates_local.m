function [max_rate, LType_rates, LType_index, LType_length, LType_Vdep_exitrate,...
          RyR_rates, RyR_index,RyR_length,Ito2_exitrate] = fru_rates_local(LType_state_FRU, RyR_state_FRU, Ito2_state_FRU, FRUdep_states, FRU_states_FRU)

 global index_frudep_V NRyRs_per_cleft Nclefts_FRU 
 index_LCC_states = 1;
 index_LCC_Vinact = 2;

 fL=0.85; % transition	rate into open state (1/ms)
 gL=2.0; %	transition rate	out	of open	state (1/ms)
 fLprime=0.005;	% transition rate into	Ca mode	open state (1/ms)
 gLprime=7.0; % transition	rate out of	Ca mode	open state (1/ms)
 bL=1.9356;	% mode	transition parameter
 bL2=bL*bL;
 bL3=bL*bL*bL;
 bL4=bL*bL*bL*bL;
 aL=2.0; %	mode transition	parameter
 aL2=aL*aL;
 aL3=aL*aL*aL;
 aL4=aL*aL*aL*aL;
 omega=0.83*2.0*1.3*0.01;  % mode transition parameter	(1/ms)

 alphacf=4.0*1.2*0.416;
 betacf=4.0*0.45*0.049;
 gammacf=0.83*1.9*1.3*0.31*7.5*0.09233;	% (ms-1 mM-1)
 
 CCa0_to_C0	= omega;		% = omega
 CCa1_to_C1	= omega/bL;	    % = omega/bL
 CCa2_to_C2	= omega/bL2;	% = omega/bL^2
 CCa3_to_C3	= omega/bL3;	% = omega/bL^3
 CCa4_to_C4	= omega/bL4;	% = omega/bL^4
 
 KdIto2=0.1502;
 kbIto2=2.0; 
 kfIto2=kbIto2/KdIto2;
 
 k21=250.0;
 k32=9.6;
 k43=0.07/0.06667*13.0;
 k45=0.07;
 k52=0.001235;
 k65=30.0; %	at 10 CaNSR	rises
 k12cf=3000.0;
 k23cf=10.0*30000.0;
 k34cf=0.6*3000.0; %	(ms-1 mM-2)
 k54cf=0.6*0.198;
 k25cf=10.0*300.0;
 k56cf=2.0*4.0*3000.0;
 
 k_rate=0.00127215;	 
 k_rate2=3.4188;
 threshCa34to7=0.0368369379834969;
 threshCa56to8=0.00011447933531005;
 threshMAXCa = 0.0504410547074504; 
 threshMAX = 2.0;
 V = FRUdep_states(index_frudep_V);
 CaSS = FRU_states_FRU(2:end);
 
 alpha =	alphacf	* exp(0.012*(V-35.0)); 
 beta = betacf *	exp(-0.05*(V-35.0));
 alpha_prime	= aL*alpha;
 beta_prime = beta/bL;

 LType_rates = zeros(Nclefts_FRU,4);
 LType_index = zeros(Nclefts_FRU,4);
 LType_length = zeros(Nclefts_FRU,1);
 LType_Vdep_exitrate = zeros(Nclefts_FRU,1);
 RyR_rates = zeros(Nclefts_FRU, NRyRs_per_cleft,4);
 RyR_index = zeros(Nclefts_FRU, NRyRs_per_cleft,4);
 RyR_length = zeros(Nclefts_FRU, NRyRs_per_cleft);
 Ito2_exitrate = zeros(Nclefts_FRU,1); 
  
 for icleft = 1:Nclefts_FRU
    gamma_rate = gammacf.*CaSS(icleft);
    switch(LType_state_FRU(icleft,index_LCC_states))
        case 1
            C0_to_C1 = 4.0*alpha;
            C0_to_CCa0 = gamma_rate;
            LType_rates(icleft,1) = (C0_to_C1+C0_to_CCa0);	% sum of rates	leaving	current	state
            LType_rates(icleft,2) = C0_to_C1;		% rate	from 1 to 2
            LType_rates(icleft,3) = C0_to_CCa0;		% rate	from 1 to 7
            LType_length(icleft) = 3;			% length=1+#states	connected to current state
            LType_index(icleft,1) = 1;
            LType_index(icleft,2) = 2;			% index(2)	state into which rates(2) takes	you	= 2
            LType_index(icleft,3) = 7;			% index(3)	state into which rates(3) takes	you	= 7
        case 2
            C1_to_C2 = 3.0*alpha;
            C1_to_C0 =		beta;
            C1_to_CCa1 = aL*gamma_rate;	% = gamma_rate*aL
            LType_rates(icleft,1) = (C1_to_C0+C1_to_C2+C1_to_CCa1);
            LType_rates(icleft,2) = C1_to_C0;
            LType_rates(icleft,3) = C1_to_C2;
            LType_rates(icleft,4) = C1_to_CCa1;
            LType_length(icleft) = 4;
            LType_index(icleft,1) = 2;
            LType_index(icleft,2) = 1;
            LType_index(icleft,3) = 3;
            LType_index(icleft,4) = 8;
        case 3
            C2_to_C1 = 2.0*beta;
            C2_to_C3 = 2.0*alpha;
            C2_to_CCa2 = aL2*gamma_rate;
            LType_rates(icleft,1) = (C2_to_C1+C2_to_C3+C2_to_CCa2);
            LType_rates(icleft,2) = C2_to_C1;
            LType_rates(icleft,3) = C2_to_C3;
            LType_rates(icleft,4) = C2_to_CCa2;
            LType_length(icleft) = 4;
            LType_index(icleft,1) = 3;
            LType_index(icleft,2) = 2;
            LType_index(icleft,3) = 4;
            LType_index(icleft,4) = 9;
        case 4
            C3_to_C4 =		alpha;
            C3_to_C2 = 3.0*beta;
            C3_to_CCa3 = aL3*gamma_rate;
            LType_rates(icleft,1) = (C3_to_C2+C3_to_C4+C3_to_CCa3);
            LType_rates(icleft,2) = C3_to_C2;
            LType_rates(icleft,3) = C3_to_C4;
            LType_rates(icleft,4) = C3_to_CCa3;
            LType_length(icleft) = 4;
            LType_index(icleft,1) = 4;
            LType_index(icleft,2) = 3;
            LType_index(icleft,3) = 5;
            LType_index(icleft,4) = 10;
        case 5
            C4_to_C3 = 4.0*beta;
            C4_to_CCa4 = aL4*gamma_rate;
            LType_rates(icleft,1) = (C4_to_C3+fL+C4_to_CCa4);
            LType_rates(icleft,2) = C4_to_C3;
            LType_rates(icleft,3) = fL;
            LType_rates(icleft,4) = C4_to_CCa4;
            LType_length(icleft) = 4;
            LType_index(icleft,1) = 5;
            LType_index(icleft,2) = 4;
            LType_index(icleft,3) = 6;
            LType_index(icleft,4) = 11;
        case 6
            LType_rates(icleft,1) = gL;
            LType_rates(icleft,2) = gL;
            LType_length(icleft) = 2;
            LType_index(icleft,1) = 6;
            LType_index(icleft,2) = 5;
        case 7
            CCa0_to_CCa1 = 4.0*alpha_prime;
            LType_rates(icleft,1) = (CCa0_to_CCa1+CCa0_to_C0);
            LType_rates(icleft,2) = CCa0_to_C0;
            LType_rates(icleft,3) = CCa0_to_CCa1;
            LType_length(icleft) = 3;
            LType_index(icleft,1) = 7;
            LType_index(icleft,2) = 1;
            LType_index(icleft,3) = 8;
        case 8
            CCa1_to_CCa2 = 3.0*alpha_prime;
            CCa1_to_CCa0 =		beta_prime;
            LType_rates(icleft,1) = (CCa1_to_CCa0+CCa1_to_CCa2+CCa1_to_C1);
            LType_rates(icleft,2) = CCa1_to_CCa0;
            LType_rates(icleft,3) = CCa1_to_C1;
            LType_rates(icleft,4) = CCa1_to_CCa2;
            LType_length(icleft) = 4;
            LType_index(icleft,1) = 8;
            LType_index(icleft,2) = 7;
            LType_index(icleft,3) = 2;
            LType_index(icleft,4) = 9;
        case 9
            CCa2_to_CCa3 = 2.0*alpha_prime;
            CCa2_to_CCa1 = 2.0*beta_prime;
            LType_rates(icleft,1) = (CCa2_to_CCa1+CCa2_to_CCa3+CCa2_to_C2);
            LType_rates(icleft,2) = CCa2_to_CCa1;
            LType_rates(icleft,3) = CCa2_to_C2;
            LType_rates(icleft,4) = CCa2_to_CCa3;
            LType_length(icleft) = 4;
            LType_index(icleft,1) = 9;
            LType_index(icleft,2) = 8;
            LType_index(icleft,3) = 3;
            LType_index(icleft,4) = 10;
            
        case 10
            CCa3_to_CCa4 =	alpha_prime;
            CCa3_to_CCa2 = 3.0*beta_prime;
            LType_rates(icleft,1) = (CCa3_to_CCa2+CCa3_to_CCa4+CCa3_to_C3);
            LType_rates(icleft,2) = CCa3_to_CCa2;
            LType_rates(icleft,3) = CCa3_to_C3;
            LType_rates(icleft,4) = CCa3_to_CCa4;
            LType_length(icleft) = 4;
            LType_index(icleft,1) = 10;
            LType_index(icleft,2) = 9;
            LType_index(icleft,3) = 4;
            LType_index(icleft,4) = 11;
        case 11
            CCa4_to_CCa3 = 4.0*beta_prime;
            LType_rates(icleft,1) = (CCa4_to_CCa3+CCa4_to_C4+fLprime);
            LType_rates(icleft,2) = CCa4_to_CCa3;
            LType_rates(icleft,3) = CCa4_to_C4;
            LType_rates(icleft,4) = fLprime;
            LType_length(icleft) = 4;
            LType_index(icleft,1) = 11;
            LType_index(icleft,2) = 10;
            LType_index(icleft,3) = 5;
            LType_index(icleft,4) = 12;
        case 12     
            LType_rates(icleft,1) = gLprime;
            LType_rates(icleft,2) = gLprime;
            LType_length(icleft) = 2;
            LType_index(icleft,1) = 12;
            LType_index(icleft,2) = 11;
    end
 end

yCa_frac=0.4;
yCa_inf	= yCa_frac/(1.0+exp((V + 12.5)/5.0)) + (1.0-yCa_frac);
tau_yCa	= 60.0 + 340.0/(1.0	+ exp((V+30.0)/12.0));
Oy_LType = 2;
Cy_LType = 1;
for icleft = 1:Nclefts_FRU
    switch (LType_state_FRU(icleft,index_LCC_Vinact))
      case Oy_LType
	    LType_Vdep_exitrate(icleft)	= (1.0-yCa_inf)/tau_yCa;
	  case Cy_LType
	    LType_Vdep_exitrate(icleft)	= yCa_inf/tau_yCa;
    end
end

for icleft = 1:Nclefts_FRU
    	Sat_term = min(threshMAX,(CaSS(icleft)*CaSS(icleft))/k_rate);
		Sat_term2 =	min(threshMAX,(CaSS(icleft)*CaSS(icleft))/k_rate2);
		k12	= k12cf	* Sat_term2;
		k23	= k23cf	* Sat_term;
		k34	= k34cf	* Sat_term;
		k54	= k54cf	* Sat_term;
		k25	= k25cf	* Sat_term;
		k56	= k56cf	* Sat_term;
        if (min(CaSS(icleft),threshMAXCa)<threshCa56to8)
            for i=1: NRyRs_per_cleft
                switch (RyR_state_FRU(icleft,i))
                  case 1
                      RyR_rates(icleft,i,1)	= k12;
                      RyR_rates(icleft,i,2)	= k12;
                      RyR_length(icleft,i) = 2;
                      RyR_index(icleft,i,1)	= 1;
                      RyR_index(icleft,i,2)	= 2;      		      
                 case 2
                      RyR_rates(icleft,i,1)	= (k21+k23+k25);
                      RyR_rates(icleft,i,2)	= k21;
                      RyR_rates(icleft,i,3)	= k23;
                      RyR_rates(icleft,i,4)	= k25;
                      RyR_length(icleft,i) = 4;
                      RyR_index(icleft,i,1)	= 2;
                      RyR_index(icleft,i,2)	= 1;
                      RyR_index(icleft,i,3)	= 3;
                      RyR_index(icleft,i,4)	= 5;     
                case 3
                      RyR_rates(icleft,i,1)	= (k32+k34);
                      RyR_rates(icleft,i,2)	= k32;
                      RyR_rates(icleft,i,3)	= k34;
                      RyR_length(icleft,i) = 3;
                      RyR_index(icleft,i,1)	= 3;
                      RyR_index(icleft,i,2)	= 2;
                      RyR_index(icleft,i,3)	= 4;				
                case 4
                      RyR_rates(icleft,i,1)	= (k43+k45);
                      RyR_rates(icleft,i,2)	= k43;
                      RyR_rates(icleft,i,3)	= k45;
                      RyR_length(icleft,i) = 3;
                      RyR_index(icleft,i,1)	= 4;
                      RyR_index(icleft,i,2)	= 3;
                      RyR_index(icleft,i,3)	= 5;
                case 5
                      RyR_rates(icleft,i,1)	= (k52+k54+k56);
                      RyR_rates(icleft,i,2)	= k52;
                      RyR_rates(icleft,i,3)	= k54;
                      RyR_rates(icleft,i,4)	= k56;
                      RyR_length(icleft,i) = 4;
                      RyR_index(icleft,i,1)	= 5;
                      RyR_index(icleft,i,2)	= 2;
                      RyR_index(icleft,i,3)	= 4;
                      RyR_index(icleft,i,4)	= 6;   
                case 6
                      RyR_rates(icleft,i,1)	= k65;
                      RyR_rates(icleft,i,2)	= k65;
                      RyR_length(icleft,i) = 2;
                      RyR_index(icleft,i,1)	= 6;
                      RyR_index(icleft,i,2)	= 5; 
                 case 7	% Not a real state, just transition to	state 4	or 3
                      dnum= rand();
                      if (dnum<(k34/(k34+k43))) 
                                RyR_state_FRU(icleft,i) = 4;
                                RyR_rates(icleft,i,1)	= (k43+k45);
                                RyR_rates(icleft,i,2)	= k43;
                                RyR_rates(icleft,i,3)	= k45;
                                RyR_length(icleft,i) = 3;
                                RyR_index(icleft,i,1)	= 4;
                                RyR_index(icleft,i,2)	= 3;
                                RyR_index(icleft,i,3)	= 5;
                      else 
                                RyR_state_FRU(icleft,i) = 3;
                                RyR_rates(icleft,i,1)	= (k32+k34);
                                RyR_rates(icleft,i,2)	= k32;
                                RyR_rates(icleft,i,3)	= k34;
                                RyR_length(icleft,i) = 3;
                                RyR_index(icleft,i,1)	= 3;
                                RyR_index(icleft,i,2)	= 2;
                                RyR_index(icleft,i,3)	= 4;
                      end 
                    case 8
                          dnum = rand();
                          if (dnum< k56/(k56+k65)) 
                                RyR_state_FRU(icleft,i) = 6;
                                RyR_rates(icleft,i,1)	= k65;
                                RyR_rates(icleft,i,2)	= k65;
                                RyR_length(icleft,i) = 2;
                                RyR_index(icleft,i,1)	= 6;
                                RyR_index(icleft,i,2)	= 5;
                         else 
                                RyR_state_FRU(icleft,i) = 5;
                                RyR_rates(icleft,i,1)	= (k52+k54+k56);
                                RyR_rates(icleft,i,2)	= k52;
                                RyR_rates(icleft,i,3)	= k54;
                                RyR_rates(icleft,i,4)	= k56;
                                RyR_length(icleft,i) = 4;
                                RyR_index(icleft,i,1)	= 5;
                                RyR_index(icleft,i,2)	= 2;
                                RyR_index(icleft,i,3)	= 4;
                                RyR_index(icleft,i,4)	= 6;
                          end
                end
            end
        else
            if (min(CaSS(icleft),threshMAXCa)<threshCa34to7)
                for i=1:NRyRs_per_cleft
                     switch (RyR_state_FRU(icleft,i))
                         case 1
                             	RyR_rates(icleft,i,1)	= k12;
                                RyR_rates(icleft,i,2)	= k12;
                                RyR_length(icleft,i) = 2;
                                RyR_index(icleft,i,1)	= 1;
                                RyR_index(icleft,i,2)	= 2;
                         case 2
                             	RyR_rates(icleft,i,1)	= (k21+k23+k25);
                                RyR_rates(icleft,i,2)	= k21;
                                RyR_rates(icleft,i,3)	= k23;
                                RyR_rates(icleft,i,4)	= k25;
                                RyR_length(icleft,i) = 4;
                                RyR_index(icleft,i,1)	= 2;
                                RyR_index(icleft,i,2)	= 1;
                                RyR_index(icleft,i,3)	= 3;
                                RyR_index(icleft,i,4)	= 8;
                         case 3
                                RyR_rates(icleft,i,1)	= (k32+k34);
                                RyR_rates(icleft,i,2)	= k32;
                                RyR_rates(icleft,i,3)	= k34;
                                RyR_length(icleft,i) = 3;
                                RyR_index(icleft,i,1)	= 3;
                                RyR_index(icleft,i,2)	= 2;
                                RyR_index(icleft,i,3)	= 4;
                         case 4 
                                RyR_rates(icleft,i,1)	= (k43+k45);
                                RyR_rates(icleft,i,2)	= k43;
                                RyR_rates(icleft,i,3)	= k45;
                                RyR_length(icleft,i) = 3;
                                RyR_index(icleft,i,1)	= 4;
                                RyR_index(icleft,i,2)	= 3;
                                RyR_index(icleft,i,3)	= 8;
                         case 5
                            RyR_state_FRU(icleft,i) = 8;		
                            a1 = k65/(k56+k65);
                            a2 = a1*k52;
                            a3 = a1*k54;
                            RyR_rates(icleft,i,1)	= (a2+a3);
                            RyR_rates(icleft,i,2)	= a2;
                            RyR_rates(icleft,i,3)	= a3;
                            RyR_length(icleft,i) = 3;
                            RyR_index(icleft,i,1)	= 8;
                            RyR_index(icleft,i,2)	= 2;
                            RyR_index(icleft,i,3)	= 4;
                         case 6
                            RyR_state_FRU(icleft,i) = 8;
                            a1 = k65/(k56+k65);
                            a2 = a1*k52;
                            a3 = a1*k54;
                            RyR_rates(icleft,i,1)	= (a2+a3);
                            RyR_rates(icleft,i,2)	= a2;
                            RyR_rates(icleft,i,3)	= a3;
                            RyR_length(icleft,i) = 3;
                            RyR_index(icleft,i,1)	= 8;
                            RyR_index(icleft,i,2)	= 2;
                            RyR_index(icleft,i,3)	= 4;
                         case 7
                             dnum = rand();
                             if (dnum< k34/(k34+k43)) 
                              RyR_state_FRU(icleft,i) = 4;
                              RyR_rates(icleft,i,1)	= (k43+k45);
                              RyR_rates(icleft,i,2)	= k43;
                              RyR_rates(icleft,i,3)	= k45;
                              RyR_length(icleft,i) = 3;
                              RyR_index(icleft,i,2)	= 4;
                              RyR_index(icleft,i,2)	= 3;
                              RyR_index(icleft,i,3)	= 8;
                            else 
                              RyR_state_FRU(icleft,i) = 3;
                              RyR_rates(icleft,i,1)	= (k32+k34);
                              RyR_rates(icleft,i,2)	= k32;
                              RyR_rates(icleft,i,3)	= k34;
                              RyR_length(icleft,i) = 3;
                              RyR_index(icleft,i,2)	= 3;
                              RyR_index(icleft,i,2)	= 2;
                              RyR_index(icleft,i,3)	= 4;
                            end
                         case 8
                                a1 = k65/(k56+k65);
                                a2 = a1*k52;
                                a3 = a1*k54;
                                RyR_rates(icleft,i,1)	= (a2+a3);
                                RyR_rates(icleft,i,2)	= a2;
                                RyR_rates(icleft,i,3)	= a3;
                                RyR_length(icleft,i) = 3;
                                RyR_index(icleft,i,1)	= 8;
                                RyR_index(icleft,i,2)	= 2;
                                RyR_index(icleft,i,3)	= 4;
                     end
                end                
            else
                for i=1:NRyRs_per_cleft
                     switch (RyR_state_FRU(icleft,i))
                         case 1
                            RyR_rates(icleft,i,1)	= k12;
                            RyR_rates(icleft,i,2)	= k12;
                            RyR_length(icleft,i) = 2;
                            RyR_index(icleft,i,1)	= 1;
                            RyR_index(icleft,i,2)	= 2;
                         case 2
                            RyR_rates(icleft,i,1)	= (k21+k23+k25);
                            RyR_rates(icleft,i,2)	= k21;
                            RyR_rates(icleft,i,3)	= k23;
                            RyR_rates(icleft,i,4)	= k25;
                            RyR_length(icleft,i) = 4;
                            RyR_index(icleft,i,1)	= 2;
                            RyR_index(icleft,i,2)	= 1;
                            RyR_index(icleft,i,3)	= 7;
                            RyR_index(icleft,i,4)	= 8;
                         case 3
                                RyR_state_FRU(icleft,i) = 7;
                                a1 = k34/(k43+k34);
                                a2 = (1.0-a1)*k32;
                                a3 = a1*k45;
                                RyR_rates(icleft,i,1)	= (a2+a3);
                                RyR_rates(icleft,i,2)	= a2;
                                RyR_rates(icleft,i,3)	= a3;
                                RyR_length(icleft,i) = 3;
                                RyR_index(icleft,i,1)	= 7;
                                RyR_index(icleft,i,2)	= 2;
                                RyR_index(icleft,i,3)	= 8;		
                         case 4 
                             RyR_state_FRU(icleft,i) = 7;
                             a1 = k34/(k43+k34);
                             a2 = (1.0-a1)*k32;
                             a3 = a1*k45;
                             RyR_rates(icleft,i,1)	= (a2+a3);
                             RyR_rates(icleft,i,2)	= a2;
                             RyR_rates(icleft,i,3)	= a3;
                             RyR_length(icleft,i) = 3;
                             RyR_index(icleft,i,1)	= 7;
                             RyR_index(icleft,i,2)	= 2;
                             RyR_index(icleft,i,3)	= 8;
                         case 5
                             	RyR_state_FRU(icleft,i) = 8;		
                                a1 = k65/(k56+k65);
                                a2 = a1*k52;
                                a3 = a1*k54;
                                RyR_rates(icleft,i,1)	= (a2+a3);
                                RyR_rates(icleft,i,2)	= a2;
                                RyR_rates(icleft,i,3)	= a3;
                                RyR_length(icleft,i) = 3;
                                RyR_index(icleft,i,1)	= 8;
                                RyR_index(icleft,i,2)	= 2;
                                RyR_index(icleft,i,3)	= 7;
                         case 6
                             RyR_state_FRU(icleft,i) = 8;
                             a1 = k65/(k56+k65);
                             a2 = a1*k52;
                             a3 = a1*k54;
                             RyR_rates(icleft,i,1)	= (a2+a3);
                             RyR_rates(icleft,i,2)	= a2;
                             RyR_rates(icleft,i,3)	= a3;
                             RyR_length(icleft,i) = 3;
                             RyR_index(icleft,i,1)	= 8;
                             RyR_index(icleft,i,2)	= 2;
                             RyR_index(icleft,i,3)	= 7;
                             
                         case 7
                             	a1 = k34/(k43+k34);
                                a2 = (1.0-a1)*k32;
                                a3 = a1*k45;
                                RyR_rates(icleft,i,1)	= (a2+a3);
                                RyR_rates(icleft,i,2)	= a2;
                                RyR_rates(icleft,i,3)	= a3;
                                RyR_length(icleft,i) = 3;
                                RyR_index(icleft,i,1)	= 7;
                                RyR_index(icleft,i,2)	= 2;
                                RyR_index(icleft,i,3)	= 8;
			
                         case 8
                                a1 = k65/(k56+k65);
                                a2 = a1*k52;
                                a3 = a1*k54;
                                RyR_rates(icleft,i,1)	= (a2+a3);
                                RyR_rates(icleft,i,2)	= a2;
                                RyR_rates(icleft,i,3)	= a3;
                                RyR_length(icleft,i) = 3;
                                RyR_index(icleft,i,1)	= 8;
                                RyR_index(icleft,i,2)	= 2;
                                RyR_index(icleft,i,3)	= 7;
                     end
                end                 
            end
        end
end

C_Ito2 = 1;
O_Ito2 = 2;
for icleft=1:Nclefts_FRU
    switch (Ito2_state_FRU(icleft)) 
        case C_Ito2 % closed
            Ito2_exitrate(icleft) =	kfIto2*CaSS(icleft);

        case O_Ito2 % open
            Ito2_exitrate(icleft) =	kbIto2;
    end
end
max_rate = 0.0;
        for icleft = 1:Nclefts_FRU
            max_rate = max(max_rate,LType_rates(icleft,1));	% Calc. max of	sum	of exit	rates
            max_rate = max(max_rate,LType_Vdep_exitrate(icleft));
            max_rate = max(max_rate,Ito2_exitrate(icleft));
            % The same	unrolled
            max_rate = max(max_rate,RyR_rates(icleft,1,1));
            max_rate = max(max_rate,RyR_rates(icleft,2,1));
            max_rate = max(max_rate,RyR_rates(icleft,3,1));
            max_rate = max(max_rate,RyR_rates(icleft,4,1));
            max_rate = max(max_rate,RyR_rates(icleft,5,1));
        end
end




