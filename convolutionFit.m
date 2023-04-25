function t = convolutionFit(x,z,a,b)
    t = conv(a*exp(b*x),z);
    tempLength = length(x);
    t = t(end-tempLength+1:end);
    maxT = max(t);
    %t = t./maxT;
end