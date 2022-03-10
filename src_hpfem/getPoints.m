function [p,nmf]=getPoints(j)

    p=zeros(3,3);
    v=zeros(4,3);
    nmf=zeros(4,3);
    
    nmf(1,1)=-0.5;
    nmf(1,2)=-(sqrt(3)/6);
    nmf(1,3)=-(sqrt(6)/12);
    
    nmf(2,1)=0.5;
    nmf(2,2)=-(sqrt(3)/6);
    nmf(2,3)=-(sqrt(6)/12);
    
    nmf(3,1)=0;
    nmf(3,2)=(sqrt(3)/3);
    nmf(3,3)=-(sqrt(6)/12);
    
    nmf(4,1)=0;
    nmf(4,2)=0;
    nmf(4,3)=(sqrt(3)/(2*sqrt(2)));
    
    v(1,1)=-1;
    v(1,2)=0;
    v(1,3)=0;
    
    v(2,1)=1;
    v(2,2)=0;
    v(2,3)=0;
    
    v(3,1)=0;
    v(3,2)=sqrt(3);
    v(3,3)=0;
    
    v(4,1)=0;
    v(4,2)=sqrt(3)/3;
    v(4,3)=2*(sqrt(2)/sqrt(3));

            if j==1
                % set up integration point locations
                p(1,1:3)=v(2,1:3);
                p(2,1:3)=v(3,1:3);
                p(3,1:3)=v(4,1:3);
                
            elseif j==2
                p(1,1:3)=v(3,1:3);
                p(2,1:3)=v(1,1:3);
                p(3,1:3)=v(4,1:3);
                
            elseif j==3
                p(1,1:3)=v(1,1:3);
                p(2,1:3)=v(2,1:3);
                p(3,1:3)=v(4,1:3);
                
            else
                p(1,1:3)=v(1,1:3);
                p(2,1:3)=v(3,1:3);
                p(3,1:3)=v(2,1:3);
                
            end


end

