
% grphn_equivSchem_v03
ptn=1001;
fmx=2.5e12; fmn=0.15e12;
freqn=fmn:(fmx-fmn)/(ptn-1):fmx;

spl=299792458;

%Layer order:
%Vacuum // layer 1 (metal) // later 2 (dielectric) // layer 3 (metal) // layer 4 (silicon) // vacuum
% nn1=1 //      nL1        //        nL2           //         nL3     //      nL4       //  nn2 = 1 A=1

nn1=1; nn2=1;
%Layer 1 graphen
tauL1=7e-14; thcL1=51*1e-9; RshL1=9e2; sig00L1=0*9e9*1/(RshL1*thcL1); n01=1;

%Layer 2 silicon
n02=3.52; tauL2=9.7e-18; thcL2=5.051*1e-4; RshL2=50e12; sig00L2=7.5*9e9;%*1/(RshL2*thcL2);

%Layer 3 gold
tauL3=1e-17; thcL3=1e-9; RshL3=3.70e2; sig00L3=0*9e9*5e7; n03=1; %RshL3=9e9/(sig00L3*thcL3);

%Layer 4 silicon
tauL4=1e-15; thcL4=1*1e-9; RshL4=50e1; sig00L4=0*9e9*1/(RshL4*thcL4); n04=1; %1/(RshL4*thcL4);


for ii=1:ptn

  omega=2*pi*freqn(ii); kkk=omega/spl;

  conL1=sig00L1/(1-1i*omega*tauL1); epsL1=n01^2+4*pi*1i*conL1/omega; nL1=sqrt(epsL1);
  conL2=sig00L2/(1-1i*omega*tauL2); epsL2=n02^2+4*pi*1i*conL2/omega; nL2=sqrt(epsL2);
  conL3=sig00L3/(1-1i*omega*tauL3); epsL3=n03^2+4*pi*1i*conL3/omega; nL3=sqrt(epsL3);
  conL4=sig00L4/(1-1i*omega*tauL4); epsL4=n04^2+4*pi*1i*conL4/omega; nL4=sqrt(epsL4);
  %nn2=nL3; nL4=nL3;


  %1st boudary xxx=0; to the right ampl is 1; reflector
  xxx=0; nnl=nL4; nnr=nn2; aar=1; bbr=0;

  aal=(aar*(1+nnr/nnl)+bbr*(1-nnr/nnl))/2;
  bbl=(aar*(1-nnr/nnl)+bbr*(1+nnr/nnl))/2;


  %2nd boudary xxx=-thcL4;
  xxx=-thcL4; nnl=nL3; nnr=nL4; aar=aal; bbr=bbl;

  aar=exp(-1i*kkk*nnr*thcL4)*aal;
  bbr=exp(1i*kkk*nnr*thcL4)*bbl;
  aal=(aar*(1+nnr/nnl)+bbr*(1-nnr/nnl))/2;
  bbl=(aar*(1-nnr/nnl)+bbr*(1+nnr/nnl))/2;


  %3rd boudary xxx=-thcL3-thcL4;
  xxx=-thcL4-thcL3; nnl=nL2; nnr=nL3; aar=aal; bbr=bbl;

  aar=exp(-1i*kkk*nnr*thcL3)*aal;
  bbr=exp(1i*kkk*nnr*thcL3)*bbl;
  aal=(aar*(1+nnr/nnl)+bbr*(1-nnr/nnl))/2;
  bbl=(aar*(1-nnr/nnl)+bbr*(1+nnr/nnl))/2;

  %4th boudary xxx=-thcL4-thcL3-thcL2;
  xxx=-thcL4-thcL3-thcL2; nnl=nL1; nnr=nL2; aar=aal; bbr=bbl;

  aar=exp(-1i*kkk*nnr*thcL2)*aal;
  bbr=exp(1i*kkk*nnr*thcL2)*bbl;
  aal=(aar*(1+nnr/nnl)+bbr*(1-nnr/nnl))/2;
  bbl=(aar*(1-nnr/nnl)+bbr*(1+nnr/nnl))/2;

 %5th boudary xxx=-thcL4-thcL3-thcL2-thcL1;
  xxx=-thcL4-thcL3-thcL2-thcL1; nnl=nn1; nnr=nL1; aar=aal; bbr=bbl;

  aar=exp(-1i*kkk*nnr*thcL1)*aal;
  bbr=exp(1i*kkk*nnr*thcL1)*bbl;
  aal=(aar*(1+nnr/nnl)+bbr*(1-nnr/nnl))/2;
  bbl=(aar*(1-nnr/nnl)+bbr*(1+nnr/nnl))/2;

  trnsm(ii)=(abs(1/aal))^2;
  reflt(ii)=(abs(bbl/aal))^(2);
  absrp(ii)=1-reflt(ii)-trnsm(ii);

end

%figure
%plot(freqn,trnsm);

%figure
%plot(freqn,reflt);

figure
plot(freqn,trnsm, 'b');%, freqn,reflt, 'g', freqn,absrp, 'r');

samplMt = dlmread('.\res-siwaf.txt')
frqcc=samplMt(:,1)*1e12;
sigcc=samplMt(:,2).^1;



figure
plot(frqcc,sigcc, 'r', freqn, trnsm, "b");
axis(inf, 0,1);


nsw(:,1)=freqn;
nsw(:,2)=reflt;

%dlmwrite('D:\UEF-Science\2022-prac\TDS-lab\Churchill_Lab Section\PPFi_Day4\ppfi1.txt',nsw,'delimiter','\t');



