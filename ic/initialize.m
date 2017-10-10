function [ state, FRU_state, Ltype_state, RyR_state, Ito2_state] = initialize(ic_states_file, ic_FRU_file, ic_LCh_file,ic_RyR_file,ic_Ito2_file)

% Read in initial conditions from specified text files 
 global Nclefts_FRU Nstates_FRU  NRyRs_per_cleft Nindepstates_LType
 global NFRU_sim N
 
% Define the indices
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

state       = zeros(N,1);
FRU_state   = zeros(NFRU_sim, Nstates_FRU);
Ltype_state = zeros(NFRU_sim, Nclefts_FRU,Nindepstates_LType);
RyR_state   = zeros(NFRU_sim, Nclefts_FRU,NRyRs_per_cleft);
Ito2_state  = zeros(NFRU_sim,Nclefts_FRU);

% Initialize the states
fileID = fopen(ic_states_file,'r');
C_data  = textscan(fileID,'%s');
state(index_V)        = str2double(C_data{1}(index_V));
state(index_mNa)      = str2double(C_data{1}(index_mNa));
state(index_hNa)      = str2double(C_data{1}(index_hNa));
state(index_jNa)      = str2double(C_data{1}(index_jNa));
state(index_Nai)      = str2double(C_data{1}(index_Nai));
state(index_Ki)       = str2double(C_data{1}(index_Ki));
state(index_Cai)      = str2double(C_data{1}(index_Cai));
state(index_CaNSR)    = str2double(C_data{1}(index_CaNSR));
state(index_LTRPNCa)  = str2double(C_data{1}(index_LTRPNCa));
state(index_HTRPNCa)  = str2double(C_data{1}(index_HTRPNCa));
state(index_xKs)      = str2double(C_data{1}(index_xKs));
state(index_C1Herg)   = str2double(C_data{1}(index_C1Herg));
state(index_C2Herg)   = str2double(C_data{1}(index_C2Herg));
state(index_C3Herg)   = str2double(C_data{1}(index_C3Herg));
state(index_OHerg)    = str2double(C_data{1}(index_OHerg));
state(index_IHerg)    = str2double(C_data{1}(index_IHerg));
state(index_C0Kv43)   = str2double(C_data{1}(index_C0Kv43));
state(index_C1Kv43)   = str2double(C_data{1}(index_C1Kv43));
state(index_C2Kv43)   = str2double(C_data{1}(index_C2Kv43));
state(index_C3Kv43)   = str2double(C_data{1}(index_C3Kv43));
state(index_OKv43)    = str2double(C_data{1}(index_OKv43));
state(index_CI0Kv43)  = str2double(C_data{1}(index_CI0Kv43));
state(index_CI1Kv43)  = str2double(C_data{1}(index_CI1Kv43));
state(index_CI2Kv43)  = str2double(C_data{1}(index_CI2Kv43));
state(index_CI3Kv43)  = str2double(C_data{1}(index_CI3Kv43));
state(index_OIKv43)  = 1.0- sum(state(index_C0Kv43:index_CI3Kv43));
state(index_C0Kv14)  = str2double(C_data{1}(index_C0Kv14));
state(index_C1Kv14)  = str2double(C_data{1}(index_C1Kv14));
state(index_C2Kv14)  = str2double(C_data{1}(index_C2Kv14));
state(index_C3Kv14)   = str2double(C_data{1}(index_C3Kv14));
state(index_OKv14)   = str2double(C_data{1}(index_OKv14));
state(index_CI0Kv14)  = str2double(C_data{1}(index_CI0Kv14));
state(index_CI1Kv14)  = str2double(C_data{1}(index_CI1Kv14));
state(index_CI2Kv14)  = str2double(C_data{1}(index_CI2Kv14));
state(index_CI3Kv14) = str2double(C_data{1}(index_CI3Kv14));
state(index_OIKv14)  =  1.0- sum(state(index_C0Kv14:index_CI3Kv14));
state(index_CaTOT)  =   str2double(C_data{1}(index_CaTOT));

state = state';

% Initialize the FRU concentrations

fileID = fopen(ic_FRU_file, 'r');
C_data1  = textscan(fileID,'%s %s %s %s %s');
FRU_state(:,1) = str2double(C_data1{1}(2:end));
FRU_state(:,2) = str2double(C_data1{2}(2:end));
FRU_state(:,3) = str2double(C_data1{3}(2:end));
FRU_state(:,4) = str2double(C_data1{4}(2:end));
FRU_state(:,5) = str2double(C_data1{5}(2:end));

fileID = fopen(ic_LCh_file);
C_data2  = textscan(fileID,'%s %s');
Ltype_state(1:end, 1, 1) = str2double(C_data2{1}(2:251));
Ltype_state(1:end, 2, 1) = str2double(C_data2{1}(252:501));
Ltype_state(1:end, 3, 1) = str2double(C_data2{1}(502:751));
Ltype_state(1:end, 4, 1) = str2double(C_data2{1}(752:1001));

Ltype_state(1:end, 1, 2) = str2double(C_data2{2}(2:251));
Ltype_state(1:end, 2, 2) = str2double(C_data2{2}(252:501));
Ltype_state(1:end, 3, 2) = str2double(C_data2{2}(502:751));
Ltype_state(1:end, 4, 2) = str2double(C_data2{2}(752:1001));

fileID = fopen(ic_RyR_file);
C_data3  = textscan(fileID,'%s %s %s %s %s');
RyR_state(1:end,1,1) = str2double(C_data3{1}(2:251));
RyR_state(1:end,1,2) = str2double(C_data3{2}(2:251));
RyR_state(1:end,1,3) = str2double(C_data3{3}(2:251));
RyR_state(1:end,1,4) = str2double(C_data3{4}(2:251));
RyR_state(1:end,1,5) = str2double(C_data3{5}(2:251));

RyR_state(1:end,2,1) = str2double(C_data3{1}(252:501));
RyR_state(1:end,2,2) = str2double(C_data3{2}(252:501));
RyR_state(1:end,2,3) = str2double(C_data3{3}(252:501));
RyR_state(1:end,2,4) = str2double(C_data3{4}(252:501));
RyR_state(1:end,2,5) = str2double(C_data3{5}(252:501));

RyR_state(1:end,3,1) = str2double(C_data3{1}(502:751));
RyR_state(1:end,3,2) = str2double(C_data3{2}(502:751));
RyR_state(1:end,3,3) = str2double(C_data3{3}(502:751));
RyR_state(1:end,3,4) = str2double(C_data3{4}(502:751));
RyR_state(1:end,3,5) = str2double(C_data3{5}(502:751));

RyR_state(1:end,4,1) = str2double(C_data3{1}(752:1001));
RyR_state(1:end,4,2) = str2double(C_data3{2}(752:1001));
RyR_state(1:end,4,3) = str2double(C_data3{3}(752:1001));
RyR_state(1:end,4,4) = str2double(C_data3{4}(752:1001));
RyR_state(1:end,4,5) = str2double(C_data3{5}(752:1001));

fileID = fopen(ic_Ito2_file );
C_data4  = textscan(fileID,'%s');
Ito2_state(1:end, 1, 1) = str2double(C_data4{1}(2:251));
Ito2_state(1:end, 2, 1) = str2double(C_data4{1}(252:501));
Ito2_state(1:end, 3, 1) = str2double(C_data4{1}(502:751));
Ito2_state(1:end, 4, 1) = str2double(C_data4{1}(752:1001));

end

