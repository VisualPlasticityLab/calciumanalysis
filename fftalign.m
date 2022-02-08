
function [u,v] = fftalign(A,B)

lines = min(size(A));
lines = lines-mod(lines,2);
if lines>600
    N = lines-200;    % leave out margin, changed by YJS 1/1/2017 to accommodate bidirection
elseif lines>120
    N = lines-40;
else
    N = lines;
end

yidx = round(size(A,1)/2)-N/2 + 1 : round(size(A,1)/2)+ N/2;
xidx = round(size(A,2)/2)-N/2 + 1 : round(size(A,2)/2)+ N/2;

A = A(yidx,xidx);
B = B(yidx,xidx);

C = fftshift(real(ifft2(fft2(A).*fft2(rot90(B,2)))));
[~,i] = max(C(:));
[ii jj] = ind2sub(size(C),i);

u = N/2-ii;
v = N/2-jj;