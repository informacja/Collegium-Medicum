%cdfMBn(x,a)
function F=cdfMBn(x,a)
t=x/a; 
F=erf(t/sqrt(2))-sqrt(2/pi).*t.*exp(-(t.^2)/2);
end