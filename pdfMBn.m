%pdfMBn(x,a)
function f=pdfMBn(x,a)
t=x/a; 
f=sqrt(2/pi)/a*(t).^2*exp(-(t.^2)/2);
end