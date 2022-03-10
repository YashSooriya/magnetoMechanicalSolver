function lfedge=get_lfedge(eltype)
if eltype==1
    lfedge(1,1)=1;
    lfedge(1,2)=5;
    lfedge(1,3)=6;
else
    lfedge(1,1)=1;
    lfedge(1,2)=6;
    lfedge(1,3)=5;
end
lfedge(2,1)=2;
lfedge(2,2)=4;
lfedge(2,3)=6;

lfedge(3,1)=3;
lfedge(3,2)=4;
lfedge(3,3)=5;

if eltype==1
    lfedge(4,1)=3;
    lfedge(4,2)=2;
    lfedge(4,3)=1;
else
    lfedge(4,1)=2;
    lfedge(4,2)=3;
    lfedge(4,3)=1;
end