function plotDistance(nrF, bf, n, j, zakres, k, Psyg, dEM, dCM, dists_chebyM)

    kol = [ "k."; "r."; "b."; "g."];
    figure(nrF); subplot(2,4,1+bf); hold on; plot(n, Psyg(j, zakres), kol(k));
    subplot(2,4,2+bf); hold on; plot(n, dEM(j, zakres), kol(k));
    subplot(2,4,3+bf); hold on; plot(n, dCM(j, zakres), kol(k));
    subplot(2,4,4+bf); hold on; plot(n, dists_chebyM(j, zakres), kol(k)); 
end
