% ================== Pobranie filtru wygładzającego ============
% ================= Synteza filtru wygładazającego MTF ============
 podzial = 1; %30
% nTu = 20 % dzielnik długości filtra
LwAm=round(length(X)/podzial); 
%lT=length(x); %Tu=Tud/5; lfA=lT/Tud; % Tu najwyższa częstość f pasma wynikowego
%Tu=LwAm'/nTu; % Zadane na zewnątrz
% nfx=0; %Tu=LwAm'/40; nfx=1; % PrzedziaĹ‚ wygĹ‚adzania widma
%for(nfw=1:1)
% ================ Synteza filtru MTF ==========================
% ................ #Wybor typu filtru .......................
if(isempty(MTF(nfw).Tu))
    typMTF=[5 5 5 3]; ntyp=1; %nfw=1; %[3 3 5 3]; % typ filtru: 3 - Z1 - filtr koncowy zwykly; 5 - z3 - trend 3.go rzÄ™du
    ntypZ=typMTF(1);
    [M, Fzw]=MTFdesign(ntypZ, Tu); %[M, Fzw, Ff, F0]=MTFdesign(ntypZ, Tu);
    MTF(nfw).Tu=Tu; MTF(nfw).M=M; MTF(nfw).Fzw=Fzw; MTF(nfw).F0=[]; MTF(nfw).Ff=[];
else Tu=MTF(nfw).Tu; M=MTF(nfw).M; Fzw=MTF(nfw).Fzw;
end
% ======================================================================
% X - wynik fft
Ayf = abs(X(1:LwAm'));
Nakt=length(Ayf); N1=MTF(1).M-1; N2=N1; Fzw=MTF(1).Fzw; lf=length(Fzw);
%if(jestAyo2) Yo2=[AYor(N1:-1:2).^2 0 0 AYor(2:end).^2]; else Yo2=[AYor(N1:-1:2) 0 0 AYor(2:end)]; end,
% .............. Uzupełniamy widmo dla ujemnych f(-N1:1), aby dobrze wygładzic początek .......  
Y=[Ayf(N1:-1:1)' Ayf']; Ll=1; N=length(Y); % najpierw wygładzamy widmo środkowoprzepustowe
% .............................................................................................
%Af(lwAm-N2)=0; % Y2=Y.^2; 
NfC=Nakt-N2; % koniec segmentu centralnego
nfw=1; %for(nfw=1:1)
    M=MTF(nfw).M; N1=M-1; N2=N1; Fzw=MTF(nfw).Fzw; lf=length(Fzw);
    Lf=M-1; NfCx=Nakt-N2;
    % ---------- Filtracja w segmencie startowym jest pominieta --------------------------------------
    n0Fzw=0;  %PnpC=N1+1; %Lf+1; %=M poczatek segmentu centralnego
    Af=[]; 
    for(n=1:LwAm-N2) %Nakt-N2)
        n1Fzw=n0Fzw+1;
        %if(nfw==1) 
        Af(n)=Y(n1Fzw:n0Fzw+lf)*Fzw; %Af2(n)=Y2(n1Fzw:n0Fzw+lf)*Fzw; end
        %AYo(n)=Yo2(n1Fzw:n0Fzw+lf)*Fzw;
        n0Fzw=n1Fzw;
    end;
    Ldf=n-1;
    % Af wygładzone widmo Ayf; 
    % --------- Sekcja koncowa jest zbędna -----------------
    ldFzw=lf;
    if (nxf)
        figure(201+nfx); plot([0:LwAm-1],Ayf,'c',0:Ldf,Af,'k'); axis('tight'); xlabel(sprintf('Widmo wygładzone w oknie fw/fmax=1/%d',LwAm/Tu));
    end %if(nfw==1) nx=find(Af2<0); Af2(nx)=0; Afm=sqrt(Af2); end,
    %nx=find(AYo<0); AYo(nx)=0;
    %if(jestAyo2) AYo=sqrt(AYo); end, %if(nfw==1) AYo(1)=AYor(1); else AYo(1)=ATr(1); end