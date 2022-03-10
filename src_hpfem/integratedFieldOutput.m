function [Field] = integratedFieldOutput(Matrices,Unknown,Options,X,dX,dX2)

%-------------------------------------------------------------------------
% Extract the post processing system matrices
%-------------------------------------------------------------------------

Paa      = Matrices.Paa;
Puu      = Matrices.Puu;
Pau      = Matrices.Pau;
kEnergy  = Matrices.kEnergy;
nSubdoms = Matrices.nSubdoms;
lenEM    = Unknown.EM.nunkt;
lenM     = Unknown.Mech.nunkt;

% if Options.nonLinear == 1
%     X(lenEM+(1:lenM))= 1/2*X(lenEM+(1:lenM));
%     dX(lenEM+(1:lenM))= 1/2*dX(lenEM+(1:lenM));
%     dX2(lenEM+(1:lenM))= 1/2*dX2(lenEM+(1:lenM));
% end

%-------------------------------------------------------------------------
% Compute the integrated fields from the solutions
%-------------------------------------------------------------------------
% Extract the individual physical fields
dA  = dX(1:lenEM+Unknown.EM.npec);
dU  = dX(lenEM+Unknown.EM.npec+(1:lenM+Unknown.Mech.npec));



% Compute the output power
dissPower = zeros(nSubdoms,1);

for i = 1:nSubdoms
    dissPower(i) = -dA'*Paa{i}*dA-dU'*Puu{i}*dU-dA'*Pau{i}*dU-dU'*((Pau{i}).')*dA;
end

% Compute the kinetic energy
kinEnergy = zeros(nSubdoms,1);
for i = 1:nSubdoms
%     kinEnergy(i,1) = abs(dU(1:2:end)'*kEnergy{i}(1:2:end,1:2:end)*dU(1:2:end));
%     kinEnergy(i,2) = abs(dU(2:2:end)'*kEnergy{i}(2:2:end,2:2:end)*dU(2:2:end));
    kinEnergy(i) = abs(dU'*kEnergy{i}*dU);
end

% Adjust the time harmonic form - where the time series is averaged
% if Options.timeDomain == 0
%     dissPower = 1/2*dissPower;
%     kinEnergy = 1/2*kinEnergy;
%     
% end


%-------------------------------------------------------------------------
% Store the computed field variables in the structure
%-------------------------------------------------------------------------
Field.dissPower = dissPower;
Field.kE        = kinEnergy;

end