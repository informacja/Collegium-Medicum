%function [pvg,pve,pvM,chi2g,chi2e,chi2M,bins,Nemp]=mhistMGE(y,lbins)
function [pvg,pve,pvM,chi2g,chi2e,chi2M,bins,Nemp]=mhistMGE(y,lbins)
% mhist - zlicza licznosc
% pvg - poziom empiryczny istotnosci
% lbins - liczba skrzynek
if(nargin<2) lbins=30; end
ldanych=length(y); %zwracamy dlugosc wektora y
ymin=min(y); %zwrocenie najmniejszej wartosci z wektora y 
ymax=max(y); %zwrocenie najwiekszej wartosci z wektora y 
ysr=mean(y); %wyliczenie wartosci sredniej wektora y
sigy=std(y); %dyspersja - odchylenie standardowe
%Przyblizenie danego rozkladu do rozkladu normalnego i rownomiernego
CN=1/(sqrt(2*pi)*sigy); %pierwszy czlon wzoru na funkcje o rozkladzie normalnym (wz str 11)  
% CN - stala rozkladu Gaussa
Dy=(ymax-ymin)/lbins; %(szerokosc pudelka) krok co ile ma rysowac slupki histogramu
[Nemp,bins]=hist(y,lbins); %Nemp - podaje empiryczne ilosci wystapien, 
                           %bins - zakresy przedzialow poszczegolnych pudelek
bar(bins,Nemp); %rysuje histogram
hold on;
ldx=500; %ilosc punktow w ktorych bedziemy liczyc wartosci funkcji
dx=(ymax-ymin)/ldx; % teoretyczny rozklad dx - krok o ktory bedziemy zwiekszac wartosc argumentu funkcji
x=ymin+dx/2; %zaczynamy od wartosci najmniejszej wektora y
% a = ysr/(2*sqrt(2/pi)); % from mean
% a = sigy*sqrt(pi/(3*pi-8)); % from var
nm=find(Nemp==max(Nemp));
a = bins(nm(1))/sqrt(2);%nbm=i; % from mode (best, seted on max( bins(max(Nemp)) )
amin=a/2; amax=1.5*a; La=2000; da=(amax-amin)/(La-1); xa=amin;
Jmin=1.e40;
for(i=1:La)
    sx=0;
    for(m=1:lbins)
        z=bins(m); fM=pdfMBn(z,xa); NfM=fM*ldanych*Dy; sx=sx+((Nemp(m)-NfM)^2)/NfM;  % NfM - liczba teoretyczna w pude³ku
    end,
    if(~isempty(sx)) if(isinf(sx)||isnan(sx)) sx=0; end; else sx=0; end, J(i)=sx; if(sx<Jmin) Jmin=sx; aOpt=xa; end
        xa=xa+da;
end, % nbm=i;ysr
a=aOpt;
for(i=1:ldx)
    fg=CN*exp(-((x-ysr)/(sqrt(2)*sigy))^2); %rozklad normalny   
    fe=0.5/sigy*exp(-abs(x-ysr)/sigy); 
    fM= pdfMBn(x,a); %rozklad Maxwella-Boltzmanna
 
    prawdxi=fg*Dy; %prawdopodobienstwo dla gaussowskiego
    Nteorg(i)=ldanych*prawdxi; %licznosc teoretyczna wg Gaussa
    xt(i)=x; %zapisanie argumentu dla ktorego wygenerowano wartosc funkcji
    x=x+dx;
    Ntexp(i)=fe*Dy*ldanych; %licznosc teoretyczna wg expon
    NtMBn(i)=fM*Dy*ldanych; %licznosc teoretyczna wg Boltzm-Maxw
end
plot(xt,Nteorg,'r', xt, Ntexp,'g',xt,NtMBn,'b'); %wyrysowanie wynikow
%Sprawdzanie zgodnosci histogramu z zalozonym rozkladem prawdopodobienstwa za pomoca statystyki chi2
% Liczymy wartoœci statystyki chi2
chi2e=0; chi2g=0; chi2M=0; pvg=0; %pvg - poziom empiryczny istotnosci
x=ymin+Dy/2; 
for(i=1:lbins) %dla kazdego slupka histogramu
    fg=CN*exp(-((x-ysr)/(sqrt(2)*sigy))^2); %rozklad normalny 
    fe=0.5/sigy*exp(-abs(x-ysr)/sigy); % rozklad wykladniczy
    fM= pdfMBn(x,a); %rozklad Maxwella-Boltzmanna
    prawdxi=fg*Dy; %prawdopodobienstwo
    Nteorg(i)=ldanych*prawdxi; %licznosc teoretyczna
    chi2g=chi2g+((Nteorg(i)-Nemp(i))^2)/Nteorg(i); 
    %chi2 dla porownania rozkladu z rozkladem normalnym
    xt(i)=x; %zapisanie argumentu dla ktorego wygenerowano wartosc funkcji
    x=x+Dy;
    prawdxe=fe*Dy; Ntexp(i)=prawdxe*ldanych;
    chi2e=chi2e+((Ntexp(i)-Nemp(i))^2)/Ntexp(i);
    %chi2 dla porownania rozkladu z rozkladem rownomiernym
    sign=sqrt(Nteorg(i)*(1-prawdxi));
    hold on; plot(xt(i),Nteorg(i)+sign,'r.',xt(i),Nteorg(i)-sign,'r.'); 
    sige=sqrt(Ntexp(i)*(1-prawdxe));
    plot(xt(i),Ntexp(i)+sige,'g.',xt(i),Ntexp(i)-sige,'g.'); hold off; 
    prawdxM=fM*Dy; NtMBn(i)=prawdxM*ldanych;
    chi2M=chi2M+((NtMBn(i)-Nemp(i))^2)/NtMBn(i);
    %chi2 dla porownania rozkladu z rozkladem rownomiernym
    sigM=sqrt(NtMBn(i)*(1-prawdxM));
    hold on; plot(xt(i),NtMBn(i)+sigM,'b.',xt(i),NtMBn(i)-sigM,'b.');
end
% Liczymy prawdopodobienstwo, ze mozliwe ch2 > chi2e, czyli pv=1-F(chi2e)
% gdzie F(chi2e) jest wart. dystryb.rozkladu chikwadrat dla chi2e; 
% Obliczone jak wyzej pv jest prawdopod. blednego odrzucenia hipotezy,
% ze nasza suma chi2e jest zmienna o rozkladzie chikwadrat, a wiec, ze
% histogram jest ZGODNY z zalozonym rozkladem
% (cdf - cumulated distribution function)
% pvg=1-chi2distr(chi2e,lbins-1); %liczymy p rawdfopodobienstwwo ze mozliwe chi2 jest wieksze od chi2e czyli 1-f(chi2e)
pvg=1-cdf('chi2',chi2g,lbins-3); %liczymy prawdfopodobienstwwo ze mozliwe chi2 jest wieksze od chi2e czyli 1-f(chi2e)
% pvg=1-f gdzie f(chi2e) jest wartoscia dystrybuanty  rozkladu chi2
%pvr jest prawodopodobienstwem odrzucenia hipotezy blednego ze nasza suma chi2e  jest zmienna o rozkladzie chi2 a wiec ze histogram jest zgodny z zalozonym  rozkladem  
% pvr=1-chi2distr(chi2r,lbins-1); %liczymy prawdfopodobienstwwo ze mozliwe chi2 jest wieksze od chi2e czyli 1-f(chi2e)
pve=1-cdf('chi2',chi2e,lbins-2);
pvM=1-cdf('chi2',chi2M,lbins-2);