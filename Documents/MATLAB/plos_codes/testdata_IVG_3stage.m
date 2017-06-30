function Y=testdata_IVG_3stage(num,m1,s1,m2,s2,m3,s3)
nu1=1/m1;
lambda1=1/s1^2;
nu2=1/m2;
lambda2=1/s2^2;
nu3=1/m3;
lambda3=1/s3^2;
pd1=makedist('InverseGaussian', 'mu',nu1,'lambda',lambda1);
pd2=makedist('InverseGaussian', 'mu',nu2,'lambda',lambda2);
pd3=makedist('InverseGaussian', 'mu',nu3,'lambda',lambda3);
Y=random(pd1,num,1)+random(pd2,num,1)+random(pd3,num,1);
end