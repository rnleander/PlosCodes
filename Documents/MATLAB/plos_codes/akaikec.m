function [AICc rel]=akaikec(ML,k,n)
%vector of maximum likelihoods
%k=vecotor of number of parameters for each model
%n=number of data points
AICc=2*k-2*ML+2*k.*(k+1)./(n-k-1);
rel=exp((min(AICc)-AICc)/2);
end

