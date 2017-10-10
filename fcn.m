function [Fstate,state,current] = fcn(time, state, current, keepc, Jxfer, Jtr, ICa, Ito2)

global Nclefts_FRU
global pulse_duration pulse_amplitude period shift
global Ko Nao Cao
global Faraday Temp Rgas RT_over_F 
global Acap VNSR VJSR VSS Vmyo
global vclamp_flag vclamp_duration vclamp_set vclamp_shift vclamp_hold vclamp_period
global index_V index_mNa index_jNa index_hNa index_Nai index_Ki index_Cai index_CaNSR index_xKs
global index_LTRPNCa index_HTRPNCa index_C0Kv43 index_C1Kv43 index_C2Kv43 index_C3Kv43 index_OKv43 
global index_CI0Kv43 index_CI1Kv43 index_CI2Kv43 index_CI3Kv43 index_OIKv43 index_C0Kv14 index_C1Kv14
global index_C2Kv14 index_C3Kv14 index_OKv14 index_CI0Kv14 index_CI1Kv14 index_CI2Kv14 index_CI3Kv14 index_OIKv14 
global index_CaTOT index_C1Herg index_C2Herg index_C3Herg index_OHerg index_IHerg
 
index_INa = 1;
index_IKr = 2;
index_IKs = 3;
index_Ito1 = 4;
index_IK1  = 5;
index_IKp = 6;
index_INaCa = 7;
index_INaK = 8;
index_IpCa = 9;
index_ICab = 10;
index_INab = 11;
index_ICa = 12;
index_JDHPR = 13;
index_Jup = 14;
index_Jtrpn = 15;
index_Jtr = 16;
index_Jxfer = 17;
index_IKv43 = 18;
index_IKv14 = 19;
index_IKv14_K = 20;
index_IKv14_Na = 21;
index_Ito2 = 22;
index_Istim = 23;
index_Itot = 24;

% Buffering parameters
	% total troponin low affinity site conc. (mM)
  LTRPNtot= 70.0e-3;
	% total troponin high affinity site conc. (mM)
  HTRPNtot= 140.0e-3;
	% Ca++ on rate for troponin high affinity sites (1/(mM*ms))
  khtrpn_plus= 20.0;
	% Ca++ off rate for troponin high affinity sites (1/ms)
  khtrpn_minus= 0.066e-3;
	% Ca++ on rate for troponin low affinity sites (1/(mM*ms))
  kltrpn_plus= 40.0;
	% Ca++ off rate for troponin low affinity sites (1/ms)
  kltrpn_minus= 40.0e-3;
	% total myoplasmic calmodulin concentration (mM)
  CMDNtot= 50.0e-3;
	% total myoplasmic EGTA concentration (mM)
  EGTAtot= 0.0;
	% Ca++ half sat. constant for calmodulin (mM)
  KmCMDN= 2.38e-3;
	% Ca++ half sat. constant for EGTA (mM)
  KmEGTA= 1.5e-4;

%	peak IKr  conductance (mS/uF) 
 GKr=0.029; % by Reza

%	peak IKs  conductance (mS/uF) 
 GKs=0.035; % 0.035	G_Ks in Canine C++ code

%	peak IK1  conductance (mS/uF) 
 GK1=3.0; % 2.73 by Reza

%  peak IKp  conductance (mS/uF) 
 GKp=1.2*0.002216; % 0.001 in stable C++

%  peak INa  conductance (mS/uF) 
 GNa=12.8;

%  scaling factor of Na+/Ca++ exchange (uA/uF) 
 kNaCa=0.9*0.30; 

 KmNa=87.5;% Na+  half sat. constant for Na+/Ca++ exch. (mM)
 KmCa=1.38;% Ca++ half sat. constant for Na+/Ca++ exch. (mM)
 KmK1=13.0;% Ca++ half sat. constant for IK1 (mM)
 ksat=0.2;% Na+/Ca++ exch. sat. factor at negative potentials 
 eta=0.35;% controls voltage dependence of Na+/Ca++ exch.
 INaKmax=1.3*0.693; % maximum Na+/K+ pump current (uA/uF)

 KmNai=10.0; % Na+  half sat. constant for Na+/K+ pump (mM)
 KmKo=1.5; % K+   half sat. constant for Na+/K+ pump (mM)
 IpCamax=0.6*0.05; % maximum sarcolemmal Ca++ pump current (uA/uF)
 KmpCa=0.0005; % half sat. constant for sarcolemmal Ca++ pump (mM)

% max. background Ca++ current conductance (mS/uF)
 GCab=3.3*7.684e-5; % 3e-5 in stable C++
%  GCab=20e-5; % 3e-5 in stable C++

% max. background Na+  current conductance (mS/uF)
% GNab=0.85*0.0031; % 0.000395 in stable C++
 GNab=0.001; % 0.000395 in stable C++

%-------SR PARAMETERS---------------------------------------------------------

 Kfb=0.26e-3; % foward half sat. constant for Ca++ ATPase (mM)
 Krb=1.8; % backward half sat. constant for Ca++ ATPase (mM)
 KSR=1.0; % scaling factor for Ca++ ATPase
 Nfb=0.75; % foward cooperativity constant for Ca++ ATPase
 Nrb=0.75; % backward cooperativity constant for Ca++ ATPase
 vmaxf=1.53*137.0e-6; % Ca++ ATPase forward rate parameter (mM/ms)
 vmaxr=1.53*137.0e-6; % Ca++ ATPase backward rate parameter (mM/ms)

%-------Kv43 AND Kv14 CURRENT PARAMETERS FOR Ito1-----------------------------

 KvScale=1.13*1.03*1.55; % Scale factor for Kv4.3 and Kv1.4 currents
 Kv43Frac=0.77; % Fraction of Ito1 which is Kv4.3 current
 GKv43=Kv43Frac*KvScale*0.1;  % Maximum conductance of Kv4.3 channel (mS/uF)
 PKv14=(1.0-Kv43Frac)*KvScale*4.792933e-7; % Permeability of Kv1.4 channel (cm/s)

%	NOTE: For step from -90mV to 40mV  
%	GKv43=0.1 and PKv14=4.792933e-7 ---> same peak current

% alphaa0Kv43,aaKv43,betaa0Kv43,baKv43 % Voltage dependent 
% alphai0Kv43,aiKv43,betai0Kv43,biKv43 % rate parameters

 alphaa0Kv43=0.543708;
 aaKv43=0.028983;
 betaa0Kv43=0.080185;
 baKv43=0.0468437;
 alphai0Kv43=0.0498424;
 aiKv43=0.000373016;
 betai0Kv43=0.000819482;
 biKv43=0.00000005374;
 alphaa0Kv14=6.0*0.31551;
 aaKv14=0.00695;
 betaa0Kv14=6.0*0.001966;
 baKv14=0.08527;
 alphai0Kv14=6.0*0.0004938;
 betai0Kv14=6.0*0.0000176;

% f1Kv43, f2Kv43, f3Kv43, f4Kv43  % Rate scaling factors
 f1Kv43=1.8936;
 f2Kv43=14.224647456;
 f3Kv43=158.574378389;
 f4Kv43=142.936645351;
 b1Kv43=6.77348;
 b2Kv43=15.6212705152;
 b3Kv43=28.7532603313;
 b4Kv43=524.576206679;
 f1Kv14=0.20005;
 f2Kv14=0.320280;
 f3Kv14=13.50909223;
 f4Kv14=1151.7651385;
 b1Kv14=2.2300;
 b2Kv14=12.000299;
 b3Kv14=5.3701338025;
 b4Kv14=5.2396395511;
 
 %-------RM HERG(+MiRP1) current parameters for IKr (10/00), H designates HERG
 	 T_Const_HERG=5.320000001;   %Temp constant from 23 to 37C,with
	 A0_HERG=0.017147641733086;  
	 B0_HERG=0.03304608038835;  
	 A1_HERG=0.03969328381141;  
	 B1_HERG=-0.04306054163980; 
	 A2_HERG=0.02057448605977;
	 B2_HERG=0.02617412715118;
	 A3_HERG=0.00134366604423;
	 B3_HERG=-0.02691385498399;
	 A4_HERG=0.10666316491288;
	 B4_HERG=0.00568908859717;
	 A5_HERG=0.00646393910049;
	 B5_HERG=-0.04536642959543;
	 A6_HERG=0.00008039374403;
	 B6_HERG=0.00000069808924;
	 C2H_to_C3H=T_Const_HERG*0.02608362043337;
	 C3H_to_C2H=T_Const_HERG*0.14832978132145;
     
% Read in state variables
    V  = state(index_V);            %    1  V
    mNa=state(index_mNa);			%	 2	m
	hNa=state(index_hNa);			%	 3	h
	jNa=state(index_jNa);			%	 4	j
	Nai=state(index_Nai);			%	 5	Nai
	Ki=state(index_Ki);		        %	 6	Ki
	Cai=state(index_Cai);			%	 7	Cai
	CaNSR=state(index_CaNSR);		%	 8	CaNSR
	xKs=state(index_xKs);			%	10	xKs
	LTRPNCa=state(index_LTRPNCa);		% 	11	LTRPNCa
	HTRPNCa=state(index_HTRPNCa);		% 	12	HTRPNCa

	C0Kv43=state(index_C0Kv43);		% 	13	Kv4.3 state C1
	C1Kv43=state(index_C1Kv43);		% 	14	Kv4.3 state C2	
	C2Kv43=state(index_C2Kv43);		% 	15	Kv4.3 state C3
	C3Kv43=state(index_C3Kv43);		% 	16	Kv4.3 state C4
	OKv43=state(index_OKv43);		% 	17	Kv4.3 state O (open);
	CI0Kv43=state(index_CI0Kv43);		% 	18	Kv4.3 state I1
	CI1Kv43=state(index_CI1Kv43);		% 	19	Kv4.3 state I2
	CI2Kv43=state(index_CI2Kv43);		% 	20	Kv4.3 state I3
	CI3Kv43=state(index_CI3Kv43);		% 	21	Kv4.3 state I4
	OIKv43 =1.0-CI0Kv43-CI1Kv43-CI2Kv43-CI3Kv43-OKv43-C0Kv43-C1Kv43-C2Kv43-C3Kv43;	% state(index_OIKv43)		% 	32	Kv1.4 state I5

	C0Kv14=state(index_C0Kv14);		% 	23	Kv1.4 state C1
	C1Kv14=state(index_C1Kv14);		% 	24	Kv1.4 state C2
	C2Kv14=state(index_C2Kv14);		% 	25	Kv1.4 state C3
	C3Kv14=state(index_C3Kv14);		% 	26	Kv1.4 state C4
	OKv14=state(index_OKv14);		% 	27	Kv1.4 state O (open);
	CI0Kv14=state(index_CI0Kv14);		% 	28	Kv1.4 state I1
	CI1Kv14=state(index_CI1Kv14);		% 	29	Kv1.4 state I2
	CI2Kv14=state(index_CI2Kv14);		% 	30	Kv1.4 state I3
	CI3Kv14=state(index_CI3Kv14);		% 	31	Kv1.4 state I4
	OIKv14=1.0-CI0Kv14-CI1Kv14-CI2Kv14-CI3Kv14-OKv14-C0Kv14-C1Kv14-C2Kv14-C3Kv14;

	CaTOT=state(index_CaTOT);		%	33	Total Cell Ca

	C1Herg=state(index_C1Herg);		%	34	HERG state C1
	C2Herg=state(index_C2Herg);		%	35	HERG state C2
	C3Herg=state(index_C3Herg);		%	36	HERG state C3
	OHerg=state(index_OHerg);		%	37	HERG state O1
	IHerg=1.0-OHerg-C1Herg-C2Herg-C3Herg;
    
    if(vclamp_flag == 0)
        time_on_Is1 = floor((time-shift)/period)*period;
        time_off_Is1 = time_on_Is1+pulse_duration;
        if ((time-shift)>=time_on_Is1&&(time-shift)<=time_off_Is1)
            Istim = pulse_amplitude;
        else
            Istim = 0.0;
        end
    else
      Istim = 0.0;
    end    
    
    if(vclamp_flag == 1)
      time_vclamp_on = floor((time-vclamp_shift)/vclamp_period)*vclamp_period;
	  time_vclamp_off = time_vclamp_on + vclamp_duration;
          if (((time-vclamp_shift) >= time_vclamp_on)  && (time_vclamp_on>=0.0) && ((time-vclamp_shift) < time_vclamp_off)) 
             ramp = (((time-vclamp_shift)-time_vclamp_on)/0.2)*(vclamp_set-vclamp_hold) + vclamp_hold;
                  if (vclamp_hold<=vclamp_set) 
                    V = min(vclamp_set,ramp); % depol.  steps
                  else 
                    V = max(vclamp_set,ramp); % hyperpol. steps
                  end
             %  V = vclamp_set;
          else 
              if (((time-vclamp_shift)<time_vclamp_on)||(time_vclamp_on<0.0)) 
                  V = vclamp_hold;
              else 
                   ramp = vclamp_set +((time_vclamp_on + vclamp_duration-(time-vclamp_shift))/0.2)*(vclamp_set-vclamp_hold);
                   if (vclamp_hold<=vclamp_set) 
                      V = max(vclamp_hold,ramp); % depol. step
                   else 
                      V = min(vclamp_hold,ramp); % hyper. step
                   end
              end
               % V = vclamp_hold;
          end
    end
%-------COMPUTE REVERSAL POTENTIALS-------------------------------------------

    	ENa = RT_over_F*log(Nao/Nai);
        EK =  RT_over_F*log(Ko/Ki);

        a1 = Ko+0.01833*Nao;
        a2 = Ki+0.01833*Nai;
        EKs = RT_over_F*log(a1/a2);
        if (Cai<=0.0)
            Cai = 1.e-15;
        end
        ECa = 0.5*RT_over_F*log(Cao/Cai);
%-------COMPUTE INa, IKr, IKs, Ito1, IK1, INab, IKp---------------------------
    if (mNa<=0.0) 
        mNa=1.e-15;
    end

    INa = GNa*(mNa*mNa*mNa)*hNa*jNa*(V-ENa);

	fKo = (Ko/4.0)^0.5;
	IKr = GKr*fKo*OHerg*(V-EK);
 
	IKs = GKs*(xKs*xKs)*(V-EKs);

	IKv43 = GKv43*OKv43*(V-EK);

	VF_over_RT=V/RT_over_F;
	VFsq_over_RT=(1000.0*Faraday)*VF_over_RT;

	exp_VFRT=exp(VF_over_RT);

    a1 =  Ki*exp_VFRT-Ko;
	a2 = exp_VFRT-1.0;

	if (abs(V)<1.e-6) 
		IKv14_K = PKv14*OKv14*1000.0*Faraday*(Ki-Ko);
	else 
		IKv14_K = PKv14*OKv14*VFsq_over_RT*(a1/a2);
	end

	a1 =  Nai*exp_VFRT-Nao;

	if (abs(V)<1.e-6) 
		IKv14_Na = 0.02*PKv14*OKv14*1000.0*Faraday*(Nai-Nao);
	else 
		IKv14_Na = 0.02*PKv14*OKv14*VFsq_over_RT*(a1/a2);
	end

	IKv14 = IKv14_K + IKv14_Na;

	Ito1 = IKv43 + IKv14;

	% New formulation endorsed by R. Mazhari. However, we are using GK1=3.0, not 2.73
	K1_inf=1.0/(0.94+exp(1.76/RT_over_F*(V-EK)));

	IK1 = GK1*Ko/(Ko+KmK1)*K1_inf*(V-EK);

	INab = GNab*(V-ENa); % orig

	KpV = 1.0/(1.0+exp((7.488-V)/5.98));
	IKp = GKp*KpV*(V-EK);
%-------COMPUTE INaK, INaCa, ICab, IpCa --------------------------------------
    VF_over_RT=V/RT_over_F;

	sigma = (exp(Nao/67.3)-1.0)/7.0;
	a1 = 1.0+0.1245*exp(-0.1*VF_over_RT);
	a2 = 0.0365*sigma/exp_VFRT;
    fNaK = 1.0/(a1+a2);
	a1 = Ko/(Ko+KmKo);
	a2 = 1.0+ (KmNai/Nai)^1.5;
	INaK = INaKmax*fNaK*(a1/a2);

	exp_etaVFRT=exp(eta*VF_over_RT);
    a1 = exp_etaVFRT*(Nai*Nai*Nai)*Cao;
	a2 = exp_etaVFRT/exp_VFRT*(Nao*Nao*Nao)*Cai;
	a3 = 1.0+ksat*exp_etaVFRT/exp_VFRT;
	a4 = KmCa+Cao;
	a5 = 5000.0/(KmNa*KmNa*KmNa+Nao*Nao*Nao);
	INaCa = kNaCa*a5*(a1-a2)/(a4*a3);

	ICab = GCab*(V-ECa);

	IpCa = IpCamax*Cai/(KmpCa+Cai);

	%-------COMPUTE GATING VARIABLE DERIVATIVES-----------------------------------

        if(abs(V+47.13)<= 1.e-4)
        alpha_m = 0.32/(0.1 - 0.005*(V+47.13));
    else
        a1 = 0.32*(V+47.13);
		alpha_m = a1/(1.0-exp(-0.1*(V+47.13)));
    end
    beta_m = 0.08*exp(-V/11.0);
    
    if (V>=-90.0) 
        dmNa = alpha_m*(1.0-mNa)-beta_m*mNa;	
    else
        dmNa = 0.0;
		mNa = alpha_m/(alpha_m + beta_m);
		state(index_mNa) = mNa;
    end
    
    alpha_h = 0.135*exp((80+V)/-6.8);

	if (V<-38.73809636838782) 
	  % curves do not cross at -40mV, but close
	  beta_h =3.56*exp(0.079*V)+310000*exp(0.35*V);
    else 
	  beta_h = 1/(0.13*(1+exp((V+10.66)/-11.1)));
    end

	if (V<-37.78) 
	  a1 = -127140*exp(0.2444*V);
	  a2 = 3.474E-5*exp(-0.04391*V);
	  a3 = 1.0+exp(0.311*(V+79.23));
	  alpha_j = (a1-a2)*(V+37.78)/a3;
    else 
	  alpha_j=0;
    end

	if (V<-39.82600037702883) 
	  % curves do not cross at -40mV, but close
	  beta_j = 0.1212*exp(-0.01052*V)/(1.0+exp(-0.1378*(V+40.14)));
	else 
	  beta_j = 0.3*exp(-2.535E-7*V)/(1+exp(-0.1*(V+32)));
    end
    
    dhNa = alpha_h*(1.0-hNa)-beta_h*hNa;	
	djNa = alpha_j*(1.0-jNa)-beta_j*jNa;
    
    % Compute derivatives of Kv4.3 channel
      alpha_act43		= alphaa0Kv43*exp(aaKv43*V);
   beta_act43			= betaa0Kv43*exp(-baKv43*V);
   alpha_inact43		= alphai0Kv43*exp(-aiKv43*V);
   beta_inact43			= betai0Kv43*exp(biKv43*V);

   C0Kv43_to_C1Kv43		= 4.0*alpha_act43;
   C1Kv43_to_C2Kv43		= 3.0*alpha_act43;
   C2Kv43_to_C3Kv43		= 2.0*alpha_act43;
   C3Kv43_to_OKv43		= alpha_act43;

   CI0Kv43_to_CI1Kv43	= 4.0*b1Kv43*alpha_act43;
   CI1Kv43_to_CI2Kv43	= 3.0*b2Kv43*alpha_act43/b1Kv43;
   CI2Kv43_to_CI3Kv43	= 2.0*b3Kv43*alpha_act43/b2Kv43;
   CI3Kv43_to_OIKv43	= b4Kv43*alpha_act43/b3Kv43;

   C1Kv43_to_C0Kv43		= beta_act43;
   C2Kv43_to_C1Kv43		= 2.0*beta_act43;
   C3Kv43_to_C2Kv43		= 3.0*beta_act43;
   OKv43_to_C3Kv43		= 4.0*beta_act43;

   CI1Kv43_to_CI0Kv43	= beta_act43/f1Kv43;
   CI2Kv43_to_CI1Kv43	= 2.0*f1Kv43*beta_act43/f2Kv43;
   CI3Kv43_to_CI2Kv43	= 3.0*f2Kv43*beta_act43/f3Kv43;
   OIKv43_to_CI3Kv43	= 4.0*f3Kv43*beta_act43/f4Kv43;

   C0Kv43_to_CI0Kv43	= beta_inact43;
   C1Kv43_to_CI1Kv43	= f1Kv43*beta_inact43;
   C2Kv43_to_CI2Kv43	= f2Kv43*beta_inact43;
   C3Kv43_to_CI3Kv43	= f3Kv43*beta_inact43;
   OKv43_to_OIKv43		= f4Kv43*beta_inact43;

   CI0Kv43_to_C0Kv43	= alpha_inact43;
   CI1Kv43_to_C1Kv43	= alpha_inact43/b1Kv43;
   CI2Kv43_to_C2Kv43	= alpha_inact43/b2Kv43;
   CI3Kv43_to_C3Kv43	= alpha_inact43/b3Kv43;
   OIKv43_to_OKv43		= alpha_inact43/b4Kv43;


   a1					= (C0Kv43_to_C1Kv43+C0Kv43_to_CI0Kv43)*C0Kv43;
   a2					= C1Kv43_to_C0Kv43*C1Kv43 + CI0Kv43_to_C0Kv43*CI0Kv43;
   dC0Kv43				= a2 - a1;
   a1					= (C1Kv43_to_C2Kv43+C1Kv43_to_C0Kv43+C1Kv43_to_CI1Kv43)*C1Kv43;
   a2					= C2Kv43_to_C1Kv43*C2Kv43 + CI1Kv43_to_C1Kv43*CI1Kv43 + C0Kv43_to_C1Kv43*C0Kv43;
   dC1Kv43				= a2 - a1;
   a1					= (C2Kv43_to_C3Kv43+C2Kv43_to_C1Kv43+C2Kv43_to_CI2Kv43)*C2Kv43;
   a2					= C3Kv43_to_C2Kv43*C3Kv43 + CI2Kv43_to_C2Kv43*CI2Kv43 + C1Kv43_to_C2Kv43*C1Kv43;
   dC2Kv43				= a2 - a1;
   a1					= (C3Kv43_to_OKv43+C3Kv43_to_C2Kv43+C3Kv43_to_CI3Kv43)*C3Kv43;
   a2					= OKv43_to_C3Kv43*OKv43 + CI3Kv43_to_C3Kv43*CI3Kv43 + C2Kv43_to_C3Kv43*C2Kv43;
   dC3Kv43				= a2 - a1;
   a1					= (OKv43_to_C3Kv43+OKv43_to_OIKv43)*OKv43;
   a2					= C3Kv43_to_OKv43*C3Kv43 + OIKv43_to_OKv43*OIKv43;
   dOKv43				= a2 - a1;
   a1					= (CI0Kv43_to_C0Kv43+CI0Kv43_to_CI1Kv43)*CI0Kv43;
   a2					= C0Kv43_to_CI0Kv43*C0Kv43 + CI1Kv43_to_CI0Kv43*CI1Kv43;
   dCI0Kv43				= a2 - a1;
   a1					= (CI1Kv43_to_CI2Kv43+CI1Kv43_to_C1Kv43+CI1Kv43_to_CI0Kv43)*CI1Kv43;
   a2					= CI2Kv43_to_CI1Kv43*CI2Kv43 + C1Kv43_to_CI1Kv43*C1Kv43 + CI0Kv43_to_CI1Kv43*CI0Kv43;
   dCI1Kv43				= a2 - a1;
   a1					= (CI2Kv43_to_CI3Kv43+CI2Kv43_to_C2Kv43+CI2Kv43_to_CI1Kv43)*CI2Kv43;
   a2					= CI3Kv43_to_CI2Kv43*CI3Kv43 + C2Kv43_to_CI2Kv43*C2Kv43 + CI1Kv43_to_CI2Kv43*CI1Kv43;
   dCI2Kv43				= a2 - a1;
   a1					= (CI3Kv43_to_OIKv43+CI3Kv43_to_C3Kv43+CI3Kv43_to_CI2Kv43)*CI3Kv43;
   a2					= OIKv43_to_CI3Kv43*OIKv43 + C3Kv43_to_CI3Kv43*C3Kv43 + CI2Kv43_to_CI3Kv43*CI2Kv43;
   dCI3Kv43				= a2 - a1;
   a1					= (OIKv43_to_OKv43+OIKv43_to_CI3Kv43)*OIKv43;
   a2					= OKv43_to_OIKv43*OKv43 + CI3Kv43_to_OIKv43*CI3Kv43;
   dOIKv43				= a2 - a1;
 %Kv1.4 channel transitions
		   alpha_act14			= alphaa0Kv14*exp(aaKv14*V);
		   beta_act14			= betaa0Kv14*exp(-baKv14*V);
		   alpha_inact14		= alphai0Kv14;
		   beta_inact14			= betai0Kv14;
		   
		   C0Kv14_to_C1Kv14		= 4.0*alpha_act14;
		   C1Kv14_to_C2Kv14		= 3.0*alpha_act14;
		   C2Kv14_to_C3Kv14		= 2.0*alpha_act14;
		   C3Kv14_to_OKv14		= alpha_act14;
		   
		   CI0Kv14_to_CI1Kv14	= 4.0*b1Kv14*alpha_act14;
		   CI1Kv14_to_CI2Kv14	= 3.0*b2Kv14*alpha_act14/b1Kv14;
		   CI2Kv14_to_CI3Kv14	= 2.0*b3Kv14*alpha_act14/b2Kv14;
		   CI3Kv14_to_OIKv14	= b4Kv14*alpha_act14/b3Kv14;
		   
		   C1Kv14_to_C0Kv14		= beta_act14;
		   C2Kv14_to_C1Kv14		= 2.0*beta_act14;
		   C3Kv14_to_C2Kv14		= 3.0*beta_act14;
		   OKv14_to_C3Kv14		= 4.0*beta_act14;
		   
		   CI1Kv14_to_CI0Kv14	= beta_act14/f1Kv14;
		   CI2Kv14_to_CI1Kv14	= 2.0*f1Kv14*beta_act14/f2Kv14;
		   CI3Kv14_to_CI2Kv14	= 3.0*f2Kv14*beta_act14/f3Kv14;
		   OIKv14_to_CI3Kv14	= 4.0*f3Kv14*beta_act14/f4Kv14;
		   
		   C0Kv14_to_CI0Kv14	= beta_inact14;
		   C1Kv14_to_CI1Kv14	= f1Kv14*beta_inact14;
		   C2Kv14_to_CI2Kv14	= f2Kv14*beta_inact14;
		   C3Kv14_to_CI3Kv14	= f3Kv14*beta_inact14;
		   OKv14_to_OIKv14		= f4Kv14*beta_inact14;
		   
		   CI0Kv14_to_C0Kv14	= alpha_inact14;
		   CI1Kv14_to_C1Kv14	= alpha_inact14/b1Kv14;
		   CI2Kv14_to_C2Kv14	= alpha_inact14/b2Kv14;
		   CI3Kv14_to_C3Kv14	= alpha_inact14/b3Kv14;
		   OIKv14_to_OKv14		= alpha_inact14/b4Kv14;
		   
		   a1					= (C0Kv14_to_C1Kv14+C0Kv14_to_CI0Kv14)*C0Kv14;
		   a2					= C1Kv14_to_C0Kv14*C1Kv14 + CI0Kv14_to_C0Kv14*CI0Kv14;
		   dC0Kv14				= a2 - a1;
		   a1					= (C1Kv14_to_C2Kv14+C1Kv14_to_C0Kv14+C1Kv14_to_CI1Kv14)*C1Kv14;
		   a2					= C2Kv14_to_C1Kv14*C2Kv14 + CI1Kv14_to_C1Kv14*CI1Kv14 + C0Kv14_to_C1Kv14*C0Kv14;
		   dC1Kv14				= a2 - a1;
		   a1					= (C2Kv14_to_C3Kv14+C2Kv14_to_C1Kv14+C2Kv14_to_CI2Kv14)*C2Kv14;
		   a2					= C3Kv14_to_C2Kv14*C3Kv14 + CI2Kv14_to_C2Kv14*CI2Kv14 + C1Kv14_to_C2Kv14*C1Kv14;
		   dC2Kv14				= a2 - a1;
		   a1					= (C3Kv14_to_OKv14+C3Kv14_to_C2Kv14+C3Kv14_to_CI3Kv14)*C3Kv14;
		   a2					= OKv14_to_C3Kv14*OKv14 + CI3Kv14_to_C3Kv14*CI3Kv14 + C2Kv14_to_C3Kv14*C2Kv14;
		   dC3Kv14				= a2 - a1;
		   a1					= (OKv14_to_C3Kv14+OKv14_to_OIKv14)*OKv14;
		   a2					= C3Kv14_to_OKv14*C3Kv14 + OIKv14_to_OKv14*OIKv14;
		   dOKv14				= a2 - a1;
		   a1					= (CI0Kv14_to_C0Kv14+CI0Kv14_to_CI1Kv14)*CI0Kv14;
		   a2					= C0Kv14_to_CI0Kv14*C0Kv14 + CI1Kv14_to_CI0Kv14*CI1Kv14;
		   dCI0Kv14				= a2 - a1;
		   a1					= (CI1Kv14_to_CI2Kv14+CI1Kv14_to_C1Kv14+CI1Kv14_to_CI0Kv14)*CI1Kv14;
		   a2					= CI2Kv14_to_CI1Kv14*CI2Kv14 + C1Kv14_to_CI1Kv14*C1Kv14 + CI0Kv14_to_CI1Kv14*CI0Kv14;
		   dCI1Kv14				= a2 - a1;
		   a1					= (CI2Kv14_to_CI3Kv14+CI2Kv14_to_C2Kv14+CI2Kv14_to_CI1Kv14)*CI2Kv14;
		   a2					= CI3Kv14_to_CI2Kv14*CI3Kv14 + C2Kv14_to_CI2Kv14*C2Kv14 + CI1Kv14_to_CI2Kv14*CI1Kv14;
		   dCI2Kv14				= a2 - a1;
		   a1					= (CI3Kv14_to_OIKv14+CI3Kv14_to_C3Kv14+CI3Kv14_to_CI2Kv14)*CI3Kv14;
		   a2					= OIKv14_to_CI3Kv14*OIKv14 + C3Kv14_to_CI3Kv14*C3Kv14 + CI2Kv14_to_CI3Kv14*CI2Kv14;
		   dCI3Kv14				= a2 - a1;
		   a1					= (OIKv14_to_OKv14+OIKv14_to_CI3Kv14)*OIKv14;
		   a2					= OKv14_to_OIKv14*OKv14 + CI3Kv14_to_OIKv14*CI3Kv14;
		   dOIKv14				= a2 - a1;
    
    % Compute the derivatives of the HERG channel
    C1H_to_C2H = T_Const_HERG*A0_HERG*exp(B0_HERG*V);
	C2H_to_C1H = T_Const_HERG*A1_HERG*exp(B1_HERG*V);
	C3H_to_OH =  T_Const_HERG*A2_HERG*exp(B2_HERG*V);
	OH_to_C3H =  T_Const_HERG*A3_HERG*exp(B3_HERG*V);
	OH_to_IH =   T_Const_HERG*A4_HERG*exp(B4_HERG*V);
	IH_to_OH =   T_Const_HERG*A5_HERG*exp(B5_HERG*V);
	C3H_to_IH =  T_Const_HERG*A6_HERG*exp(B6_HERG*V);
	IH_to_C3H =  (OH_to_C3H*IH_to_OH*C3H_to_IH)/(C3H_to_OH*OH_to_IH);
    
    dC1Herg = C2H_to_C1H * C2Herg - C1H_to_C2H * C1Herg;
	a1 = C1H_to_C2H * C1Herg + C3H_to_C2H * C3Herg;
	a2 = (C2H_to_C1H + C2H_to_C3H) * C2Herg;
	dC2Herg = a1-a2;
	a1 = C2H_to_C3H*C2Herg + OH_to_C3H*OHerg + IH_to_C3H*IHerg;
	a2 = (C3H_to_IH + C3H_to_OH + C3H_to_C2H) * C3Herg; 
	dC3Herg = a1-a2;			
	a1 = C3H_to_OH * C3Herg + IH_to_OH * IHerg;
	a2 = (OH_to_C3H + OH_to_IH) * OHerg;
	dOHerg = a1-a2;	
    
    xKs_inf = 1.0/(1.0+exp(-(V-24.7)/13.6));

	if(abs(V-10)<1.e-6)  
		a1 = (7.19e-5)/0.148;
		a2 = (1.31e-4)/0.0687;
	else 
		a1 = 7.19e-5*(V-10.0)/(1.0-exp(-0.148*(V-10.0)));
		a2 = 1.31e-4*(V-10.0)/(exp(0.0687*(V-10.0))-1.0);
    end
	tau_xKs = 1.0/(a1+a2);
	dxKs = (xKs_inf-xKs)/tau_xKs;
    
    fb = (Cai/Kfb)^Nfb;
	rb = (CaNSR/Krb)^Nrb;
	Jup = KSR*(vmaxf*fb - vmaxr*rb)/(1.0 + fb + rb);
        
    a1 = kltrpn_minus * LTRPNCa;
 	dLTRPNCa = kltrpn_plus*Cai*(1.0 - LTRPNCa) - a1;


	a1 = khtrpn_minus * HTRPNCa;
 	dHTRPNCa = khtrpn_plus*Cai*(1.0 - HTRPNCa) - a1;


	Jtrpn = LTRPNtot*dLTRPNCa+HTRPNtot*dHTRPNCa;

	a1 = CMDNtot*KmCMDN/((Cai+KmCMDN)^2);
	a2 = EGTAtot*KmEGTA/((Cai+KmEGTA)^2.0);
	beta_i = 1.0/(1.0+a1+a2);

    a1 = Acap/(Vmyo*Faraday*1000.0); 
    dNai = -( INa+INab+3.0*(INaCa+INaK)+IKv14_Na )*a1;

    
    a3 = IKr+IKs+IK1+IKp;
	dKi = -( a3-2.0*INaK+(Ito1-IKv14_Na) +Istim)*a1;

    a3 = ICab -2.0*INaCa + IpCa;
 	dCai = beta_i*(Jxfer*(VSS/Vmyo) -Jup -Jtrpn - a3*0.5*a1);    
    dCaNSR = Jup*(Vmyo/VNSR) - Jtr*(VJSR/VNSR);


	dCaTOT =  -(a3+ICa)*0.5*a1*Vmyo*1.0e6;
    
    a1 = INa+ICa+IKr+IKs;
	a2 = IK1+IKp+INaCa+INaK+Ito1+Ito2;
	a3 = IpCa+ICab+INab;
    
    if(vclamp_flag)
        dV = 0.0;
        state(index_V) = V;
        Itot = a1 + a2 + a3;
    else
        Itot = a1 + a2 + a3 + Istim;
        dV = -Itot;
    end
    
    Fstate(index_V)   = dV;
	Fstate(index_mNa) = dmNa;
	Fstate(index_hNa) = dhNa;
	Fstate(index_jNa) = djNa;
	Fstate(index_Nai) = dNai;
	Fstate(index_Ki) = dKi;
	Fstate(index_Cai) = dCai;
	Fstate(index_CaNSR) = dCaNSR;
	Fstate(index_xKs) = dxKs;
	Fstate(index_LTRPNCa) = dLTRPNCa;
	Fstate(index_HTRPNCa) = dHTRPNCa;
	Fstate(index_C0Kv43) = dC0Kv43;
	Fstate(index_C1Kv43) = dC1Kv43;
	Fstate(index_C2Kv43) = dC2Kv43;	
	Fstate(index_C3Kv43) = dC3Kv43;	
	Fstate(index_OKv43) = dOKv43;	
	Fstate(index_CI0Kv43) = dCI0Kv43;
	Fstate(index_CI1Kv43) = dCI1Kv43;	
	Fstate(index_CI2Kv43) = dCI2Kv43;	
	Fstate(index_CI3Kv43) = dCI3Kv43;	
	Fstate(index_OIKv43) =  dOIKv43;
	Fstate(index_C0Kv14) = dC0Kv14;
	Fstate(index_C1Kv14) = dC1Kv14;	
	Fstate(index_C2Kv14) = dC2Kv14;	
	Fstate(index_C3Kv14) = dC3Kv14;	
	Fstate(index_OKv14) = dOKv14;	
	Fstate(index_CI0Kv14) = dCI0Kv14;
	Fstate(index_CI1Kv14) = dCI1Kv14;	
	Fstate(index_CI2Kv14) = dCI2Kv14;	
	Fstate(index_CI3Kv14) = dCI3Kv14;	
	Fstate(index_OIKv14) = dOIKv14;
	Fstate(index_CaTOT) = dCaTOT;
	Fstate(index_C1Herg) = dC1Herg;
	Fstate(index_C2Herg) = dC2Herg;
	Fstate(index_C3Herg) = dC3Herg;
	Fstate(index_OHerg) = dOHerg;
	Fstate(index_IHerg) = -dC1Herg-dC2Herg-dC3Herg-dOHerg;
if(keepc == 1)    
    current(index_INa) = INa;
    current(index_IKr) = IKr;
    current(index_IKs) = IKs;
    current(index_Ito1) = Ito1;
    current(index_IK1) = IK1;
    current(index_IKp) = IKp;
    current(index_INaCa) = INaCa;
    current(index_INaK) = INaK;
    current(index_IpCa) = IpCa;
    current(index_ICab) = ICab;
    current(index_INab) = INab;
    current(index_ICa) = ICa;
    current(index_JDHPR) = ICa*Acap/(-2.0*Vmyo*Faraday*1000.0);
    current(index_Jup) = Jup;
    current(index_Jtrpn) = Jtrpn ;
    current(index_Jtr) = Jtr*VJSR/Vmyo;
    current(index_Jxfer) = Jxfer*VSS/Vmyo;
    current(index_IKv43) = IKv43;
    current(index_IKv14) = IKv14 ;
    current(index_IKv14_K) = IKv14_K ;
    current(index_IKv14_Na) = IKv14_Na;
    current(index_Ito2) = Ito2;
    current(index_Istim) = Istim;
    current(index_Itot) = Itot;
end

        
end

