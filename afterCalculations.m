% for Ludwin
figure(1), close 1
fpom % sampl. freq
xf % freq. vector

v;
Widma;
wyglWidma;
CentrWidm;

dists_chebyM;

CC;
CCE;

% eg
figure(1),
plot(v(1).dataB), hold on, plot(v(1).dataR), hold off; legend
nexttile, 
plot(CentrWidm(1,1).AfM,"DisplayName","CentrWidm"); hold on,
plot(wyglWidma(1,1).Af, "DisplayName","wyglWidma"); hold off; legend; axis tight