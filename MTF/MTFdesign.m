%function [M, F0, Fzw, Ff, Ffshf, nh1, N1, N2, lf,ktory]=MTFdesign(kkk, Tu)
% Uproszczona wersja procedury desMTFilters(kkk, Tu, wDel, dfMV, typ_wag, synt_zwykla)
function [M, Fzw, Ff, F0, Ffshf,ktory, nh1, N1, N2, lf]=MTFdesign(kkk, Tu)
%function [M, F0, Fzw, Ff, Ffshf, Bdel, ATu, N1, N2, lf, ktory, nh1]=desMTFilters(kkk, Tu, wDel, dfMV, typ_wag, synt_zwykla)

Fzw=[]; Fz=[]; F0=[]; F00=[]; Fc=[]; Ff=[]; FA=[]; W=[]; WT=[]; U=[]; Fs=[]; F=[];
wDel=-2/3; dfMV=1; typ_wag=0; synt_zwykla=0;
typZ=['z0'; 'z1'; 'z1'; 'z2'; 'z3'; 'm2'; 'm3'; 'f '; 's3'; 's4'; 'zs']; Nh1=[1 1 1.38 2 2*1.112 1.38 2 2 1.38 1.78 1.38]; %1.35
synt_zwykla=0;  % jest uzywana tylko synt_anal, bo jest sprawdz.
synt_analit=~synt_zwykla;
nh1=Nh1(kkk);
ktory=typZ(kkk,:); %'z0';  nh1=1;
M=round(Tu*nh1);
if(M<3) M=3; end
Bdel=zeros(1,4);
ATuk =[0.4180    0.0002    0.0001    0.0000;...
    0.9035    0.3236    0.0791    0.0345;...
    0.7214    0.0501    0.0326    0.0147;...
    0.7734    0.1036    0.0207    0.0088;...
    0.9597    0.1777    0.0581    0.0275;...
    0.7002    0.0684    0.0316    0.0146;...
    0.7752    0.0968    0.0214    0.0094;...
    0.7417    0.1143    0.0242    0.0104;...
    %    0.6486    0.0918    0.0369    0.0163;...
    0.6486    0.09245    0.0369    0.0163;...
    0.7141    0.1215    0.0274    0.0120;...
    0.7734    0.1036    0.0207    0.0088];...
    BdelDane=[  -0.4964   -0.4724    0.1097    0.0857; ...
    -0.1784   -0.4783   -0.6300   -0.3302; ...
    -0.1977   -0.5755   -0.8444   -0.4666; ...
    -0.0942   -0.3126   -0.4977   -0.2793; ...
    -0.0415   -0.0559   -0.0397   -0.0252; ...
    -0.2645   -0.6437   -0.6725   -0.2933; ...
    -0.1577   -0.5613   -0.8383   -0.4347; ...
    -0.1269   -0.4430   -0.6804   -0.3643; ...
    -0.3061   -0.6093   -0.3955   -0.0922; ...
    -0.2066   -0.6277   -0.7551   -0.3340; ...
    -0.3061   -0.6093   -0.3955   -0.0922];
lfiltr=1;
switch(ktory)
    case 'z0', istot=[1 0 0 0 0 0 0];  ATu=0; Bdel=[-0.5044   -0.4799    0.1110    0.0865]; % z1 -srednia ruchoma
    case 'z1', istot=[1 1 0 0 0 0 0];  ATu=0; Bdel=[-0.2031   -0.5360   -0.7381   -0.4053]; % z1 -zwykly trend liniowy
    case 'z2', istot=[1 1 1 0 0 0 0];  ATu=0; Bdel=[-0.0861   -0.1760   -0.2337   -0.1438]; % z2 -zwykly trend 2go rzedu
    case 'z3', istot=[1 1 1 1 0 0 0];  ATu=0; Bdel=[0.0209    0.2522    0.4928    0.2615];% z3 -zwykly trend 3go rzedu
    case 'm2', istot=[1 0 1 0 0 0 0];  ATu=0; Bdel=[-0.2649   -0.6053   -0.5617   -0.2213];% m2 maksimum w chwili n (dla t=0) rzad 2       ; d1(0)=0.;   d2(0)=2*a2; d3(0)=0;    d4(0)=0;
    case 'm3', istot=[1 0 1 1 0 0 0];  ATu=0; Bdel=[-0.1599   -0.4904   -0.6366   -0.3062];% m3 maksimum w chwili n (dla t=0) rzad 3       ; d1(0)=0.;   d2(0)=2*a2; d3(0)=6*a3; d4(0)=0;
    case 'f ', istot=[1 1 0 0 1 0 0];  ATu=0; Bdel=[-0.1349   -0.3596   -0.4615   -0.2368];% f - rzad 4 punkt przegiêcia w chwili n (t=0)  ; d1(0)=a1;   d2(0)=0.0;  d3(0)=0;    d4(0)=24*a4;
    case 's3', istot=[1 0 0 1 0 0 0];  ATu=[0.6486    0.0918    0.0369    0.0163]; Bdel=[-0.3113   -0.5799   -0.2741   -0.0054]; % s - punkt siodlowy w chwili n (dla t=0) rzad 3; d1(0)=0.;d2(0)=0.0;  d3(0)=6*a3; d4=0
        %                                                                                -0.3113   -0.5799   -0.2741   -0.0054
    case 's4', istot=[1 0 0 1 1 0 0];  ATu=0; Bdel=[-0.2094   -0.5921   -0.6509   -0.2682];% s - punkt siodlowy w chwili n (dla t=0) rzad 4; d1(0)=0.;   d2(0)=0.0;  d3(0)=6*a3; d4(0)=24*a4
    case 'zs', istot=[1 1 0 0 0 0 0]; istotf=[1 0 0 1 0 0 0]; lfiltr=2; ATu=[0.6486    0.0918    0.0369    0.0163]; Bdel=[-0.3113   -0.5799   -0.2741   -0.0054]; % wariant mieszany z1 i s3
        %otherwise, istot=[1 1 0 0 0 0 0]; istotf=[1 0 0 0 0 0 0]; lfiltr=2; ATu=[0.6486    0.0918    0.0369    0.0163]; Bdel=[-0.3113   -0.5799   -0.2741   -0.0054]; % wariant mieszany m0 i s3
end
ATu=ATuk(kkk,2:4);
if(wDel<0) Bdel=BdelDane(kkk,:); end% Bdel=1.5*Bdel;
for(nfiltr=1:lfiltr)
    if(nfiltr==2)
        U=[]; W=[]; WT=[]; Ww=[]; UTW=[]; KA=[]; FA=[]; Wsig=[];
        istot=istotf;
    end
    maxR=1; for(m=2:length(istot)) if(istot(m)>0) maxR=m; end, end
    U=[ones(M,1)];
    for(k=2:maxR)
        if(istot(k)>0) U=[U (([-M+1:0]/M).^(k-1))']; end % ([-M+1:0].^2)' ([-M+1:0].^3)'];
    end
    if(typ_wag>0)
        Ww=0; sW=0; W=zeros(M,M);
        for(i=1:M)
            Ww=Ww+1/(M-i+1)/M; w=sqrt(Ww); W(i,i)=w; sW=sW+w; Wd(i)=w;
        end
        UTW=U'*W/sW; sW,
    else UTW=U';
    end
    KA=[];
    KA=inv(UTW*U);
    FA=KA*UTW; Wsig=U*KA; W=Wsig*U'; %W=U*inv(U'*U)*U';
    WT=W';
    if(synt_analit)
        % Projektowanie analityczne -----
        % ---------------------- Filtr koncowy --------------------
        % M=4;
        if(nargout>2)
            Ffa=zeros(2*M-1,M-1);
            for(p=1: M-1)
                np=M-p;
                for(k=1:M+p-1)
                    nk=2*M-2-k+2; %
                    ip=max(1,1+k-M); ikp=min(p,k);
                    if(M<5) fprintf(1,'\np=%2d k=%2d; ip=%2d ik=%2d',p,k,ip,ikp); end
                    for(i=ip:ikp)
                        rW=M-p+i; kW=M-k+i;
                        mm=k-i+1;
                        if(mm>0 & mm<=M)
                            if(M<5) fprintf(1,'; i=%2d', i); end,
                            Ffa(nk,np)=Ffa(nk,np)+W(rW,kW);
                        else if(M<5) fprintf(1,' !! i=%d m=%d !!', i,mm); end
                        end
                    end;
                    Ffa(nk,np)=Ffa(nk,np)/p;
                end
            end
        else Ffa=[]; Ff=[]; F0=[];
        end
    end
    if(synt_analit & nfiltr==1)
        % filtr centralny
        ppm=1;
        Fza=zeros(2*M-1,ppm);
        for(pp=1:ppm)
            p=pp+M-1;
            np=ppm-pp+1;
            for(k=1:2*M-1)
                nk=2*M-1-k+1;
                ip=max(1,1+k-M); ikp=min(M,k);
                if(M<5) fprintf(1,'\npp=%2d k=%2d; ip=%2d ik=%2d',pp,k,ip,ikp); end
                for(i=ip:ikp)
                    rW=M-(M-i+1)+1; kW=M-(k-i+1)+1;
                    rW=i; kW=M-k+i;
                    mm=k-i+1;
                    if(mm>0 & mm<=M)
                        if(M<5) fprintf(1,'; i=%2d', i); end, %fprintf(1,'; i=%2d; m=%2d', i,mm); end
                        Fza(nk,np)=Fza(nk,np)+W(rW,kW);
                    else if(M<5) fprintf(1,' !! i=%d m=%d !!', i,mm); end
                    end
                end;
                Fza(nk,np)=Fza(nk,np)/M;
            end
        end
        if(0)
            % Filtr startowy ---------------
            F0a=zeros(2*M-1,M-1);
            for(p=1: M-1)
                np=p; np=M-p;
                for(k=1:2*M-2)
                    nk=2*M-2-k+1;
                    if(M<5) fprintf(1,'\np=%2d k=%2d',p,k); end
                    ip=p; ikp=min(M-1,k);
                    for(i=p:ikp)
                        rW=i-p+1; kW=M-k+i;
                        mm=k-i+1;
                        if(mm>0 & mm<=M)
                            if(M<5) fprintf(1,'; i=%2d rW=%d kW=%d', i,rW,kW); end, %fprintf(1,'; i=%2d; m=%2d', i,mm); end
                            F0a(nk,np)=F0a(nk,np)+W(rW,kW);
                        else if(M<5) fprintf(1,' !! i=%d m=%d !!', i,mm); end
                        end
                    end;
                    F0a(nk,np)=F0a(nk,np)/(M-p);
                end
            end
            if(M<5) F0a, end
        end
    end
    %[size(Ffa);  size(Fza); size(F0a)], %max(abs(Fza(:,1)-Fza(:,2))),
    if(M<5) return, end
    if(nfiltr==1)
        %F0=F0a;
        Fzw=Fza; lf=length(Fzw); N2=M-1; N1=N2; Ff=Ffa;
        if(nargout>3) F0=Ff(lf:-1:1,N2:-1:1); end
        if(synt_zwykla==0) F0a=[]; Fza=[]; Ffa=[]; end,
        %[lf,N1]=size(F0); [lf,N2]=size(Ff); N1f=N1+1;
        Lf=N2;
    else  % usredniamy filtr koncowy
        K=round(0.4*M); if(K<1) K=1; else if(K>Lf-1) K=Lf-1; end, end
        for(k=1:K) Ff(:,k)=Ff(:,k)+(Ffa(:,k)-Ff(:,k))*(k-1)/K; end
        Ff(:,K+1:Lf)=Ffa(:,K+1:Lf);
    end
end
symFfsk=0; Fapr=0; figFfsk=0; kzn=0; kznf=0; if(nargout>4) FiltrFfsh; end
