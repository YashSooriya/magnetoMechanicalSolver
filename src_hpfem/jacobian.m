function [axi,aeta,azeta,asxi,aseta,aszeta,det]=jacobian(xy,xi,eta,zeta,...
                                                    eltype,le,lf,flag,gorder,mycoord)


gesizet=(gorder+1+1)*(gorder+1+2)*(gorder+1+3)/6;
 if flag==0
     gesizet=4;
 end
if gesizet>=500
    error(message('increase array size in jacobian'));
end
axi(1)=1;
axi(2)=0;
axi(3)=0;
aeta(1)=0;
aeta(2)=1;
aeta(3)=0;
azeta(1)=0;
azeta(2)=0;
azeta(3)=1;
      
ph=gbasish1(xi,eta,zeta,axi,aeta,azeta,gorder+1,eltype,500);
      
while 1% 100      
%mycoord = zeros(gesizet,3);

%transfer coefficents to locations vertices
%mycoord(1:4,1:3) = xy(1:4,1:3);

% if flag==1
%     % edge functions
%     for i=1:6
%         for p=1:gorder
%             for j=1:3
%                 mycoord(4+i+6*(p-1),j)=le(i,p,j);
%             end
%         end
%     end
%     
% % face functions
%     for i=1:4
%         for ii=1:(gorder-1)*gorder/2
%             for j=1:3
%                 mycoord(4+6*gorder+(i-1)*gorder*(gorder-1)/2+ii,j)= lf(i,ii,j);
%             end
%         end
%     end
%     
% else
%     gesizet=4;
% end

jac(1:3,1:3)=ph(1:gesizet,1:3)'*mycoord(1:gesizet,1:3);

% jac = zeros(3,3);
% for p=1:gesizet
%     jac(1:3,1:3)=jac(1:3,1:3)+(ph(p,1:3)'*mycoord(p,1:3));
% end

det=(jac(1,1)*((jac(2,2)*jac(3,3))-(jac(2,3)*jac(3,2))))-...
    (jac(1,2)*((jac(2,1)*jac(3,3))-(jac(2,3)*jac(3,1))))+...
    (jac(1,3)*((jac(2,1)*jac(3,2))-(jac(2,2)*jac(3,1))));

if det<0
    display(message('Mapping invalid!!! det < 0!!!!'));
    if flag==1
        flag=0;
        display('correcting')
        continue
    else
        error(message('Mapping invalid!!! det < 0!!!!'));
    end
end
break
end % 100
% compute asxi,aseta,aszeta      
asxi(1:3)=jac(1,1:3);
aseta(1:3)=jac(2,1:3);
aszeta(1:3)=jac(3,1:3);

v=scvectpd(asxi,aseta,aszeta);

if v==0
    error(message('error deviding by zero in jacobian.m'));
end
           
% compute axi,aeta and azeta

axi(1)=((aseta(2)*aszeta(3))-(aseta(3)*aszeta(2)))/v;
axi(2)=-1.d0*((aseta(1)*aszeta(3))-(aseta(3)*aszeta(1)))/v;
axi(3)=((aseta(1)*aszeta(2))-(aseta(2)*aszeta(1)))/v;

aeta(1)=((aszeta(2)*asxi(3))-(aszeta(3)*asxi(2)))/v;
aeta(2)=-1.d0*((aszeta(1)*asxi(3))-(aszeta(3)*asxi(1)))/v;
aeta(3)=((aszeta(1)*asxi(2))-(aszeta(2)*asxi(1)))/v;

azeta(1)=((asxi(2)*aseta(3))-(asxi(3)*aseta(2)))/v;
azeta(2)=-1.d0*((asxi(1)*aseta(3))-(asxi(3)*aseta(1)))/v;
azeta(3)=((asxi(1)*aseta(2))-(asxi(2)*aseta(1)))/v;

%calculate det g for curl evaluation	
g(1,1) = asxi*asxi';
g(1,2) = asxi*aseta';
g(1,3) = asxi*aszeta';
g(2,1) = g(1,2);
g(2,2) = aseta*aseta';
g(2,3) = aseta*aszeta';
g(3,1) = g(1,3);
g(3,2) = g(2,3);
g(3,3) = aszeta*aszeta';
	
detg=(g(1,1)*((g(2,2)*g(3,3))-(g(2,3)*g(3,2))))-...
    (g(1,2)*((g(2,1)*g(3,3))-(g(2,3)*g(3,1))))+...
    (g(1,3)*((g(2,1)*g(3,2))-(g(2,2)*g(1,3))));
     
for i=1:3
    asxi(i)=asxi(i)/sqrt(detg);
    aseta(i)=aseta(i)/sqrt(detg);
    aszeta(i)=aszeta(i)/sqrt(detg);
end