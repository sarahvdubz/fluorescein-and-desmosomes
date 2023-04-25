function t = convolutionFit(x,z,a,b,c,d)
    t = conv((a*exp(b*x)+c*exp(d*x)),z);
    tempLength = length(x);
    [maxT, maxTIndex] = max(t);
    t = t(end-tempLength+1:end);
end