function [xsr,sig,sumG,sumR,chi2kr,pufE]=mhistogr(y,lbins,puf) 
%funkcja przyjmuje wektor i liczbe pudelek oraz przedzial ufnosci
%teraz okreslamy jego minimum, maksimum, odchylenie standardowe,srednia oraz ilosc probek
xmin=min(y);                            %zwraca œredni¹ arytmetyczn¹ wektora liczb
xmax=max(y);                            %odchylenie standartowe od sredniej
sig=std(y);                             %zwraca najmniejsz¹ wartoœæ w wektorze liczb
xsr=mean(y);                            %zwraca najwiêksz¹ wartoœæ w wektorze liczb
ldanych=length(y);                      %d³ugosc y
dy=(xmax-xmin)/(lbins-1);               %liczymy przedzial - po ile bedzie w pudelku
[Nemp,bins]=hist(y,lbins);              %%Ndemp-lczba danych empiryczna,wektor licznosni probek w pudelkach;
%figure;
bar(bins,Nemp);                 %rysujemy histogram
                                %ponizej wyznaczamy i rysujemy histogram teoretyczny oraz rozklad rownomierny
C=1/sqrt(2*pi)/sig;              %stala Gaussa
lxteor=500;                     %liczba x-ow rozkladu
dx=(xmax-xmin)/(lxteor-1);
sumG=0;sumR=0;
for (i=1:lbins)
    x=bins(i);
    fg=C*exp(-(((x-xsr)/sig)^2)/2);             %funckaj gestosci rozkladu normalnego
    fR=1/(xmax-xmin);                           %funkcja gestosci rozkladu rownomiernego
    iprawd=fg*dy;                               %prawdopodobienstwo wystapienia danej probki w pudelku
    indteor=iprawd*ldanych;                     %teoretyczny rozklad normlany (prawdopodobienstwo)
    sumG=sumG+((indteor-Nemp(i))^2)/indteor;
    indteor=fR*dy*ldanych;                      %to samo dla rozkladu rownomiernego
    sumR=sumR+((indteor-Nemp(i))^2)/indteor;    
end
chi2kr=chi2inv(puf,lbins-1);                    %krytyczna wartosc statyskyki chi^2
pufE=cdf('chi2',sumG,lbins-1);                  %empiryczny przedzial ufnosci
i=1;
for(x=xmin:dx:xmax)
    xt(i)=x;                                    %funckaj gestosci rozkladu normalnego
    fn(i)=C*exp(-(((x-xsr)/sig)^2)/2);          %funckaj gestosci rozkladu normalnego
    fr(i)=1/(xmax-xmin);                        %funkcja gestosci rozkladu rownomiernego
    prawd(i)=fn(i)*dy;
    ndteor(i)=prawd(i)*ldanych;                 %teoretyczny rozklad normalny - wyliczony z prawdopodobienstwa
    ntrown(i)=fr(i)*dy*ldanych;                 %teoretyczny rozklad rownomierny
    i=i+1;
end
    hold on;
%     plot(xt,ndteor,'r',[xmin xt xmax],[0 ntrown 0],'g');
