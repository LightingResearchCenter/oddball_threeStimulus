function sim=alpha_stable_rand(N,alpha,gamma)
%
% Random numbers with a alpha-stable distribution
%
% sim=alpha_stable_rand(N,alpha,gamma);
%
% INPUT PARAMETERS-----------------------------------------------------
%
% N:            Number of samples in output sim
% alpha:        The scaling parameter of the alpha-stable distribution
% gamma:        The spread of the alpha-stable distribution  
% 
% OUTPUT VARIABLES-----------------------------------------------------
%
% sim:          N alpha-stable distributed random numbers
%
% ---------------------------------------------------------------
% Written by Espen A. F. Ihlen (espen.ihlen@ntnu.no), 2010

Numb1=rand(N,1);
Numb2=rand(N,1);
phi=pi.*(Numb1-1/2);
sim=gamma.*(((-log10(Numb2).*cos(phi))./cos((1-alpha).*phi)).^(1-(1/alpha))).*(sin(alpha.*phi)./cos(phi));
%Gauss_contr=2.*scale.*sqrt(-log10(Numb2)).*sin(phi);
%Cauchy_contr=scale.*tan(phi);
%Levy_contr=-scale.*(tan(phi)./(2.*log10(Numb2).*cos(phi)));