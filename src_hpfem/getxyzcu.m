function [x,y,z]= getxyzcu(xy,xi,eta,zeta,le,lf,flag,gorder,eltype)

gesizet=(gorder+1+1)*(gorder+1+2)*(gorder+1+3)/6;
if gesizet>=500
    error(message('increase array size in getzyzcu'));
end

 ph=basish1(gesizet,xi,eta,zeta,gorder+1,eltype);

% zero array
mycoord = zeros(gesizet,3);

% transfer coefficents to locations vertices
mycoord(1:4,1:3) = xy(1:4,1:3);

if flag==1
    % edge functions
    for i=1:6
        for p=1:gorder
            for j=1:3
                mycoord(4+i+6*(p-1),j)=le(i,p,j);
            end
        end
    end
    
    % face functions
    for i=1:4
        for ii=1:(gorder-1)*gorder/2
            for j=1:3
                mycoord(4+6*gorder+(i-1)*gorder*(gorder-1)/2+ii,j)=lf(i,ii,j);
            end
        end
    end
else
    gesizet=4;
end

x = ph(1:gesizet)'*mycoord(1:gesizet,1);
y = ph(1:gesizet)'*mycoord(1:gesizet,2);
z = ph(1:gesizet)'*mycoord(1:gesizet,3);