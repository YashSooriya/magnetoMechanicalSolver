function [Sol]= LagrangeInterpolation(freqRange,U,S,V,freq, Options)

% Function for interpolation of POD results using Lagrange polynomials
% Polynomial of order 3 (4 points needed)

%======================================================================
% Inputs
%======================================================================
% freqRange: vector conatining the values of frequency where we computed the
%snapshots
% U: matrix containing spatial modes
% V: matrix containing parametric modes
% S: eigenvalues
% freq: frequency or frequencies where we want the solution

%==========================================================================
% Outputs
%==========================================================================
% Sol: Solution at desired frequencies (Matrix)
%==========================================================================

freqL=length(freq);
[rows ~]=size(U);
Sol=zeros(rows,freqL);
VintWhole=zeros(size(V));
for k=1:freqL
    if any(freq(k)==freqRange)
        j=find(freqRange==freq(k));
        VintWhole(k, :)=V(j,:);
        Sol(:,k)=U*S*V(j,:)';
    else
        freqDiff=abs(freqRange'-freq(k)*ones(length(freqRange),1));
        i=find(freqDiff==min(freqDiff));
        if length(i)>1
            i=i(1);
        end
        if freqRange(i)<freq(k)
           w=freq(k);
           if i>1 && i<length(freqRange)-1
           w0=freqRange(i-1);
           w1=freqRange(i);
           w2=freqRange(i+1);
           w3=freqRange(i+2);
           l0=((w-w1)/(w0-w1))*((w-w2)/(w0-w2))*((w-w3)/(w0-w3));
           l1=((w-w0)/(w1-w0))*((w-w2)/(w1-w2))*((w-w3)/(w1-w3));
           l2=((w-w0)/(w2-w0))*((w-w1)/(w2-w1))*((w-w3)/(w2-w3));
           l3=((w-w0)/(w3-w0))*((w-w1)/(w3-w1))*((w-w2)/(w3-w2));
           Vint=V(i-1,:)*l0+V(i,:)*l1+V(i+1,:)*l2+V(i+2,:)*l3;
           elseif i==1
           w0=freqRange(i);
           w1=freqRange(i+1);
           w2=freqRange(i+2);
           w3=freqRange(i+3);
           l0=((w-w1)/(w0-w1))*((w-w2)/(w0-w2))*((w-w3)/(w0-w3));
           l1=((w-w0)/(w1-w0))*((w-w2)/(w1-w2))*((w-w3)/(w1-w3));
           l2=((w-w0)/(w2-w0))*((w-w1)/(w2-w1))*((w-w3)/(w2-w3));
           l3=((w-w0)/(w3-w0))*((w-w1)/(w3-w1))*((w-w2)/(w3-w2));
           Vint=V(i,:)*l0+V(i+1,:)*l1+V(i+2,:)*l2+V(i+3,:)*l3;
           elseif i==length(freqRange)-1
           w0=freqRange(i-2);
           w1=freqRange(i-1);
           w2=freqRange(i);
           w3=freqRange(i+1);
           l0=((w-w1)/(w0-w1))*((w-w2)/(w0-w2))*((w-w3)/(w0-w3));
           l1=((w-w0)/(w1-w0))*((w-w2)/(w1-w2))*((w-w3)/(w1-w3));
           l2=((w-w0)/(w2-w0))*((w-w1)/(w2-w1))*((w-w3)/(w2-w3));
           l3=((w-w0)/(w3-w0))*((w-w1)/(w3-w1))*((w-w2)/(w3-w2));
           Vint=V(i-2,:)*l0+V(i-1,:)*l1+V(i,:)*l2+V(i+1,:)*l3;
           end
           VintWhole(k, :)=Vint;
           Sol(:,k)=U*S*Vint';       
        else
           w=freq(k);
           if i>2 && i<length(freqRange)
           w0=freqRange(i-2);
           w1=freqRange(i-1);
           w2=freqRange(i);
           w3=freqRange(i+1);
           l0=((w-w1)/(w0-w1))*((w-w2)/(w0-w2))*((w-w3)/(w0-w3));
           l1=((w-w0)/(w1-w0))*((w-w2)/(w1-w2))*((w-w3)/(w1-w3));
           l2=((w-w0)/(w2-w0))*((w-w1)/(w2-w1))*((w-w3)/(w2-w3));
           l3=((w-w0)/(w3-w0))*((w-w1)/(w3-w1))*((w-w2)/(w3-w2));
           Vint=V(i-2,:)*l0+V(i-1,:)*l1+V(i,:)*l2+V(i+1,:)*l3;
           elseif i==2
           w0=freqRange(i-1);
           w1=freqRange(i);
           w2=freqRange(i+1);
           w3=freqRange(i+2);
           l0=((w-w1)/(w0-w1))*((w-w2)/(w0-w2))*((w-w3)/(w0-w3));
           l1=((w-w0)/(w1-w0))*((w-w2)/(w1-w2))*((w-w3)/(w1-w3));
           l2=((w-w0)/(w2-w0))*((w-w1)/(w2-w1))*((w-w3)/(w2-w3));
           l3=((w-w0)/(w3-w0))*((w-w1)/(w3-w1))*((w-w2)/(w3-w2));
           Vint=V(i-1,:)*l0+V(i,:)*l1+V(i+1,:)*l2+V(i+2,:)*l3;
           elseif i==length(freqRange)
           w0=freqRange(i-3);
           w1=freqRange(i-2);
           w2=freqRange(i-1);
           w3=freqRange(i);
           l0=((w-w1)/(w0-w1))*((w-w2)/(w0-w2))*((w-w3)/(w0-w3));
           l1=((w-w0)/(w1-w0))*((w-w2)/(w1-w2))*((w-w3)/(w1-w3));
           l2=((w-w0)/(w2-w0))*((w-w1)/(w2-w1))*((w-w3)/(w2-w3));
           l3=((w-w0)/(w3-w0))*((w-w1)/(w3-w1))*((w-w2)/(w3-w2));
           Vint=V(i-3,:)*l0+V(i-2,:)*l1+V(i-1,:)*l2+V(i,:)*l3;
           end

           VintWhole(k, :)=Vint;
           Sol(:,k)=U*S*Vint'; 
        end
    end
end
saveLocationG = sprintf('data/GCurve_lagrange_Ns%d_m%d', Options.noSnapshots, Options.nModes);
save(saveLocationG, 'VintWhole');
         


