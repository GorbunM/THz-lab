function reflt = optima(p, thcEXP)  
    ptn = 1;
    for jj = 1:length(thcEXP)
        freqn = 2e12;
        spl=299792458;
        
        nn1=1; nn2=1;
        %Layer 1 Air
        n01=1;          thcL1=1; 
        
        %Layer 2 Lorentz
        thcL2 = thcEXP(jj);

        %Layer 3 Drude
        n03=1;          thcL3=2*1e-6;    tauL3=1e-19;   sig00L3=1*9e9*5e7; 
        
        %Layer 4 Air
        n04=1;          thcL4=1; 
        
        epsDi = zeros(1, ptn);
        nDi = zeros(1, ptn);
        deltaAu = zeros(1, ptn);
        
        for ii=1:ptn
        
          omega=2*pi*freqn(ii); kkk=omega/spl;
          nL1 = n01;
        
          epsL2= p(1)+1j*p(2);
          nL2=sqrt(epsL2);
          epsDi(ii) = epsL2;
          nDi(ii) = nL2;
        
          conL3=sig00L3/(1-1i*omega*tauL3); epsL3 = n03^2 + 4*pi*1i * conL3 ./ omega; nL3=sqrt(epsL3);
          
          nL4=n04;
        
          %1st boudary
          nnl=nL4; nnr=nn2; aar=1; bbr=0;
        
          aal=(aar*(1+nnr/nnl)+bbr*(1-nnr/nnl))/2;
          bbl=(aar*(1-nnr/nnl)+bbr*(1+nnr/nnl))/2;
        
        
          %2nd boudary 
          nnl=nL3; nnr=nL4; 
        
          aar=exp(-1i*kkk*nnr*thcL4)*aal;
          bbr=exp(1i*kkk*nnr*thcL4)*bbl;
          aal=(aar*(1+nnr/nnl)+bbr*(1-nnr/nnl))/2;
          bbl=(aar*(1-nnr/nnl)+bbr*(1+nnr/nnl))/2;
        
        
          %3rd boudary 
          nnl=nL2; nnr=nL3; 
        
          aar=exp(-1i*kkk*nnr*thcL3)*aal;
          bbr=exp(1i*kkk*nnr*thcL3)*bbl;
        
          deltaAu(ii) = kkk*nnr*thcL3;
        
          aal=(aar*(1+nnr/nnl)+bbr*(1-nnr/nnl))/2;
          bbl=(aar*(1-nnr/nnl)+bbr*(1+nnr/nnl))/2;
        
          %4th boudary 
          nnl=nL1; nnr=nL2;
        
          aar=exp(-1i*kkk*nnr*thcL2)*aal;
          bbr=exp(1i*kkk*nnr*thcL2)*bbl;
          aal=(aar*(1+nnr/nnl)+bbr*(1-nnr/nnl))/2;
          bbl=(aar*(1-nnr/nnl)+bbr*(1+nnr/nnl))/2;
        
          %5th boudary 
          nnl=nn1; nnr=nL1;
        
          aar=exp(-1i*kkk*nnr*thcL1)*aal;
          bbr=exp(1i*kkk*nnr*thcL1)*bbl;
          aal=(aar*(1+nnr/nnl)+bbr*(1-nnr/nnl))/2;
          bbl=(aar*(1-nnr/nnl)+bbr*(1+nnr/nnl))/2;
        
          % trnsm(ii)=(abs(1/aal))^2;     
        end
        reflt(jj)=(abs(bbl/aal))^(2);   
        % disp(reflt)
    end
    % figure
    % plot(thcEXP, reflt, 'b');
end
