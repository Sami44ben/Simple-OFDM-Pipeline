clear; clc; close all;

offlineMode=true;

FFTsize=1024;
df=30e3;
fs=df*FFTsize;
OSR=2;
bits=2;
validBits=[1 2 4 6 8 10];
modNames={'BPSK','QPSK','16QAM','64QAM','256QAM','1024QAM'};   % BPSK (not pi/2) for bits=1
assert(ismember(bits,validBits),'bits must be 1,2,4,6,8,or 10');
modType=modNames{bits==validBits};

NumActiveSubcarriers=FFTsize*2/4;
CPsize=FFTsize/8;
N_msg=10000;
k0=FFTsize/2+1;
Kneg=NumActiveSubcarriers/2;
Kpos=Kneg;

ratio=2;
pilotindx=2:ratio:NumActiveSubcarriers;
dataindx=setdiff(1:NumActiveSubcarriers,pilotindx);
NdPerSym=numel(dataindx);
pilotVal=(0+1j);

targetRate=1/2;
rv=0;
max_iter=25;

SDR_A='sn: 1044730a1997001610002a0036067c6324';
SDR_B='sn: 10447318ac0f001601000e00998717fdbc';
plutoIDs.Rx_ID=SDR_A;
plutoIDs.Tx_ID=SDR_B;

hw_stp.Tx_G=-10;
hw_stp.Rx_G=50;
hw_stp.Fc_tx=2.3e9;
hw_stp.Fc_rx=2.3e9;
hw_stp.R=OSR*fs;
hw_stp.N=400e3;
hw_stp.Tx_ID=plutoIDs.Tx_ID;
hw_stp.Rx_ID=plutoIDs.Rx_ID;

sim.SNRdB=10;            % set to a realistic value to see the coding gain (Inf/1000 = noiseless)
sim.CFOHz=2000;
sim.delay=20*OSR;
sim.h=[1;zeros(3*OSR-1,1);0.25*exp(1j*0.8);zeros(5*OSR-1,1);0.08*exp(-1j*0.4)];

% ---------------- TX ----------------
rng(12345);
MessageBits=randi([0 1],N_msg,1);
cbsInfo=nrDLSCHInfo(N_msg,targetRate);
bgn=cbsInfo.BGN;
tbCRC=nrCRCEncode(MessageBits,cbsInfo.CRC);
cbs=nrCodeBlockSegmentLDPC(tbCRC,bgn);
N_coded=ceil((N_msg/targetRate)/(NdPerSym*bits))*NdPerSym*bits;
NumSymbols=N_coded/(NdPerSym*bits);
codedBits=nrRateMatchLDPC(nrLDPCEncode(cbs,bgn),N_coded,rv,modType,1);

MessageMap=qammod(codedBits,2^(bits),"InputType","bit","UnitAveragePower",true);

activeGrid=zeros(NumActiveSubcarriers,NumSymbols);

activeGrid(dataindx,:)=reshape(MessageMap,NdPerSym,NumSymbols);

activeGrid(pilotindx,:)=pilotVal;

FDgrid_c=zeros(FFTsize,NumSymbols);

FDgrid_c(k0-Kneg:k0-1,:)=activeGrid(1:Kneg,:);

FDgrid_c(k0+1:k0+Kpos,:)=activeGrid(Kneg+1:end,:);

timeSignal=ifft(ifftshift(FDgrid_c,1),FFTsize)*sqrt(FFTsize);

timeSignal_cp=[timeSignal(end-CPsize+1:end,:);timeSignal];

timeSignal_t=timeSignal_cp(:)/rms(timeSignal_cp(:));

N_zc=FFTsize;
qList = primes(N_zc-1);
q_zc = qList(find(gcd(qList,N_zc)==1,1));
n_zc=(0:N_zc-1).';
prmbl=exp(-1j*pi*q_zc*n_zc.*(n_zc+1)/N_zc);
prmbl=prmbl/rms(prmbl);
prmbl_tx=[prmbl;prmbl];
CPsize_pre=CPsize;

% FIX 3b: preamble UNIT RMS (drop the *sqrt(2) that made it dominate the payload)
prmbl_block=[prmbl_tx(end-CPsize_pre+1:end);prmbl_tx];
gaurd=round(0.1*length(timeSignal_t));
timeSignal1=[zeros(gaurd,1);prmbl_block;timeSignal_t;zeros(gaurd,1)];
%
Txsig=repelem(timeSignal1,OSR);
% FIX 3c: peak-normalize the WHOLE frame once (DAC headroom), preserving preamble:payload ratio
Txsig=Txsig/max(abs(Txsig))*0.8;
L_hw=N_zc*OSR;
hw_stp.N = max(hw_stp.N, ceil(3*numel(Txsig)+sim.delay+numel(sim.h)+L_hw));
hw_stp.N = hw_stp.N + mod(hw_stp.N,2);
if ~offlineMode
    Tx=sdrtx('Pluto','Gain',hw_stp.Tx_G,'CenterFrequency',hw_stp.Fc_tx,'BasebandSampleRate',hw_stp.R,'RadioID',hw_stp.Tx_ID);
    Rx=sdrrx('Pluto','CenterFrequency',hw_stp.Fc_rx,'BasebandSampleRate',hw_stp.R,'SamplesPerFrame',hw_stp.N,'GainSource','AGC Fast Attack','OutputDataType','double','RadioID',hw_stp.Rx_ID);
    transmitRepeat(Tx,Txsig);
    disp('Starting real-time processing loop. Press Ctrl+C to stop.');
else
    simRep=ceil((5*hw_stp.N+sim.delay+numel(sim.h))/numel(Txsig));
    simStream=[zeros(sim.delay,1);filter(sim.h,1,repmat(Txsig,simRep,1))];
    n=(0:numel(simStream)-1).';
    simStream=simStream.*exp(1j*2*pi*sim.CFOHz*n/hw_stp.R);
    if isfinite(sim.SNRdB)
        actP=mean(abs(simStream(abs(simStream)>1e-6)).^2);
        simStream=simStream+sqrt(actP/10^(sim.SNRdB/10)/2)*(randn(size(simStream))+1j*randn(size(simStream)));
    end
    simPtr=1;
    disp('Starting continuous simulated channel loop. Press Ctrl+C to stop.');
end

% ---------------- RX loop ----------------
while true
    if offlineMode
        idx=mod((simPtr:simPtr+hw_stp.N-1)-1,numel(simStream))+1;
        r=simStream(idx);
        simPtr=simPtr+hw_stp.N;
    else
        r=Rx();
    end
    r=r/rms(r);

    % ---- FIX 1: bounded Schmidl-Cox metric + energy gate + CFO at peak ----
    prod_seq=conj(r(1:end-L_hw)).*r(L_hw+1:end);
    P  = filter(ones(L_hw,1),1,prod_seq);
    P  =P(L_hw:end);
    Ra = filter(ones(L_hw,1),1,abs(r(1:end-L_hw)).^2); 
    Ra =Ra(L_hw:end);
    Rb = filter(ones(L_hw,1),1,abs(r(L_hw+1:end)).^2); 
    Rb =Rb(L_hw:end);
    M  = abs(P).^2./(Ra.*Rb+eps);          % bounded in [0,1] by Cauchy-Schwarz
    epow=Ra+Rb;
    Mg = M.*(epow>0.5*max(epow));          % gate out low-energy (guard/zero) regions
    peak_detect = 0.1;    
    if max(Mg)<peak_detect
        fprintf('No frame detected in buffer. Listening...\n');
        continue
    end
    [~,dpk]=max(Mg);                       % peak of gated metric = reliable point
    cfo_coarse=angle(P(dpk))/(2*pi*L_hw/hw_stp.R);
    smpls_rx=r.*exp(-1j*2*pi*cfo_coarse*(0:numel(r)-1)'/hw_stp.R);

    % ---- multi-phase symbol-timing ----
    lk=zeros(OSR,1);
    mf=conv(smpls_rx,repelem(prmbl_tx,OSR));
    for ii=1:OSR
        rk=downsample(mf(ii:end),OSR);
        lk(ii)=sum(abs(rk).^2);
    end
    [~,sampling_time]=max(lk);
    symbs=downsample(smpls_rx(sampling_time:end),OSR);

    % ---- fine timing via cross-correlation (reference is unit-RMS prmbl_block) ----
    [c,lags]=xcorr(symbs,prmbl_block);
    [~,ii]=max(abs(c));
    start_idx=lags(ii)+CPsize_pre;
    if start_idx<0 || start_idx+2*N_zc+numel(timeSignal_t)>numel(symbs)
        continue
    end

    my_frame=symbs(start_idx+1:start_idx+2*N_zc+numel(timeSignal_t));
    my_preamble=my_frame(1:2*N_zc);
    phaseDiff=my_preamble(N_zc+1:2*N_zc).*conj(my_preamble(1:N_zc));
    % ---- FIX 2: angle(sum(.)) not mean(angle(.)) ----
    cfo_fine=angle(sum(phaseDiff))/(2*pi*N_zc/(hw_stp.R/OSR));
    finalCorrectedSignal=my_frame.*exp(-1j*2*pi*cfo_fine*(0:numel(my_frame)-1)'/(hw_stp.R/OSR));

    % ---- OFDM demod ----
    data_sym_comp=finalCorrectedSignal(2*N_zc+1:end);
    sym1=reshape(data_sym_comp,FFTsize+CPsize,NumSymbols);
    sym3=fft(sym1(CPsize+1:end,:))/sqrt(FFTsize);
    Ypay_c=fftshift(sym3,1);
    YpayAct=[Ypay_c(k0-Kneg:k0-1,:);Ypay_c(k0+1:k0+Kpos,:)];

    % ---- channel est + equalize ----
    pilots=pilotVal*ones(numel(pilotindx),1);
    Hpil=mean(YpayAct(pilotindx,:)./pilots,2);
    Z=YpayAct./interp1(pilotindx(:),Hpil,(1:NumActiveSubcarriers).','linear','extrap');

    % ---- per-symbol CPE + SFO tracking (pilot WLS) ----
    pilot_k=pilotindx(:);
    k_all=(1:NumActiveSubcarriers).';
    for s=1:NumSymbols
        Zp=Z(pilot_k,s);
        w=max(abs(Zp),1e-3).^2;
        A=[pilot_k,ones(numel(pilot_k),1)];
        p=(A.*sqrt(w))\(unwrap(angle(Zp.*conj(pilots))).*sqrt(w));
        Z(:,s)=Z(:,s).*exp(-1j*(p(1)*k_all+p(2)));
    end

    rxSyms=Z(dataindx,:);
    rxSyms=rxSyms(:);
    pilotErr=reshape(Z(pilotindx,:)-pilots,[],1);
    noiseVarEst=max(mean(abs(pilotErr).^2),1e-12);
    SNRdB_est=10*log10(max(mean(abs(rxSyms).^2)-noiseVarEst,eps)/noiseVarEst);
    llrs = qamdemod(rxSyms,2^bits,"OutputType","approxllr","UnitAveragePower",true, "NoiseVariance",noiseVarEst);
    
    % ---- LDPC decode + CRC ----
    deRm=nrRateRecoverLDPC(llrs(1:N_coded),N_msg,targetRate,rv,modType,1);
    decCbs=nrLDPCDecode(deRm,bgn,max_iter);
    [blk,blkErr]=nrCodeBlockDesegmentLDPC(decCbs,bgn,N_msg+cbsInfo.L);
    [MessageBits_rx,tbErr]=nrCRCDecode(blk,cbsInfo.CRC);

    nCompare=min(N_msg,numel(MessageBits_rx));
    bitErr=sum(double(MessageBits_rx(1:nCompare))~=MessageBits(1:nCompare));
    fprintf('BER=%.4f (%d/%d) | CFO=%.1f Hz | SNRest=%.1f dB | TBCRC=%d | CBCRC=%d\n',...
        bitErr/nCompare,bitErr,nCompare,cfo_coarse+cfo_fine,SNRdB_est,tbErr,any(blkErr));
end
