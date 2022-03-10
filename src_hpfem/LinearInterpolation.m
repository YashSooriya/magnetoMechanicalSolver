function [Sol]= LinearInterpolation(freqRange,U,S,V,freq)

% Function for linearly interpolate the results of POD

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
for k=1:freqL
    if any(freq(k)==freqRange)
        j=find(freqRange==freq(k));
        Sol(:,k)=U*S*V(j,:)';
    else
        freqDiff=abs(freqRange'-freq(k)*ones(length(freqRange),1));
        i=find(freqDiff==min(freqDiff));
        if length(i)>1
            i=i(1);
        end
        if freqRange(i)<freq(k)
            Vint=V(i,:)+(freq(k)-freqRange(i))*((V(i+1,:)-V(i,:))/(freqRange(i+1)-freqRange(i)));
            Sol(:,k)=U*S*Vint';
        else
            Vint=V(i-1,:)+(freq(k)-freqRange(i-1))*((V(i,:)-V(i-1,:))/(freqRange(i)-freqRange(i-1)));
            Sol(:,k)=U*S*Vint';
        end
    end
end
            
        
        
    



