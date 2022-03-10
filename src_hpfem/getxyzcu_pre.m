function [x,y,z]= getxyzcu_pre(ph,mycoord,gesizet)


x = ph(1:gesizet)'*mycoord(1:gesizet,1);
y = ph(1:gesizet)'*mycoord(1:gesizet,2);
z = ph(1:gesizet)'*mycoord(1:gesizet,3);