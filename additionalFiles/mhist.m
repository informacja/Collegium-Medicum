function [chi2e,pe_ist,chi2r,pe_istr]=mhist(y,lbins)
% mhist - zlicza licznosc
% pe_ist - poziom empiryczny istotnosci
% lbins - liczba skrzynek
ldanych=length(y); %zwracamy dlugosc wektora y
ymin=min(y); %zwrocenie najmniejszej wartosci z wektora y 
ymax=max(y); %zwrocenie najwiekszej wartosci z wektora y 
ysr=mean(y); %wyliczenie wartosci sredniej wektora y
sigy=std(y); %dyspersja - odchylenie standardowe
%Przyblizenie danego rozkladu do rozkladu normalnego i rownomiernego
CN=1/(sqrt(2*pi)*sigy); %pierwszy czlon wzoru na funkcje o rozkladzie normalnym (wz str 11)  
% CN - stala rozkladu Gaussa
Dy=(ymax-ymin)/lbins; %(szerokosc pudelka) krok co ile ma rysowac slupki histogramu
[Nemp,bins]=hist(y,lbins); %Nemp - podaje empiryczne ilosci wystapin, 
                           %bins - zakresy przedzialow poszczegolnych pudelek
bar(bins,Nemp); %rysuje histogram
hold on;
ldx=500; %ilosc punktow w ktorych bedziemy liczyc wartosci funkcji
dx=(ymax-ymin)/ldx; % teoretyczny rozklad dx - krok o ktory bedziemy zwiekszac wartosc argumentu funkcji
x=ymin; %zaczynamy od wartosci najmniejszej wektora y
for(i=1:ldx)
    fg=CN*exp(-((x-ysr)/(sqrt(2)*sigy))^2); %rozklad normalny 
    fr= 1/(ymax-ymin);%rozklad rownomierny
    prawdxi=fg*Dy; %prawdopodobienstwo
    Nteorg(i)=ldanych*prawdxi; %licznosc teoretyczna
    xt(i)=x; %zapisanie argumentu dla ktorego wygenerowano wartosc funkcji
    x=x+dx;
    Ntrow(i)=fr*Dy*ldanych;
end
plot(xt,Nteorg,'r', [ymin xt ymax], [0,Ntrow,0],'g'); %wyrysowanie wynikow
%Sprawdzanie zgodnosci histogramu z zalozonym rozkladem prawdopodobienstwa za pomoca statystyki chi2
chi2e=0; chi2r=0; pe_ist=0; %pe_ist - poziom empiryczny istotnosci
x=ymin+Dy/2; 
for(i=1:lbins) %dla kazdego slupka histogramu
    fg=CN*exp(-((x-ysr)/(sqrt(2)*sigy))^2); %rozklad normalny 
    fr= 1/(ymax-ymin);%rozklad rownomierny
    prawdxi=fg*Dy; %prawdopodobienstwo
    Nteorg(i)=ldanych*prawdxi; %licznosc teoretyczna
    chi2e=chi2e+((Nteorg(i)-Nemp(i))^2)/Nteorg(i); 
    %chi2 dla porownania rozkladu z rozkladem normalnym
    xt(i)=x; %zapisanie argumentu dla ktorego wygenerowano wartosc funkcji
    x=x+Dy;
    Ntrow(i)=fr*Dy*ldanych;
    chi2r=chi2r+((Ntrow(i)-Nemp(i))^2)/Ntrow(i);
    %chi2 dla porownania rozkladu z rozkladem rownomiernym
end

% Liczymy prawdopodobienstwo, ze mozliwe ch2 > chi2e, czyli p_ist=1-F(chi2e)
% gdzie F(chi2e) jest wart. dystryb.rozkladu chikwadrat dla chi2e; 
% Obliczone jak wyzej p_ist jest prawdopod. blednego odrzucenia hipotezy,
% ze nasza suma chi2e jest zmienna o rozkladzie chikwadrat, a wiec, ze
% histogram jest ZGODNY z zalozonym rozkladem
% (cdf - cumulated distribution function)
pe_ist=1-cdf('chi2',chi2e,lbins-1); %liczymy prawdfopodobienstwwo ze mozliwe chi2 jest wieksze od chi2e czyli 1-f(chi2e)
%pe_ist=1-f gdzie f(chi2e) jest wartoscia dystrybuanty  rozkladu chi2
%pe_istr jest prawodopodobienstwem odrzucenia hipotezy blednego ze nasza suma chi2e  jest zmienna o rozkladzie chi2 a wiec ze histogram jest zgodny z zalozonym  rozkladem  
pe_istr=1-cdf('chi2',chi2r,lbins-1);