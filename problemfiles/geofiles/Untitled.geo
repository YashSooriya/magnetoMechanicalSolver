algebraic3d

# Geometry for Toy problem

# Main coil 1

solid main1out = ellipticcylinder (0, 0, 0; 0,594.7e-3,0;900e-3,0,0);
solid main1in = cylinder (0, 0, 0;0,0,194.8e-3; 250e-3);

solid cube=orthobrick(-10,-10,194.7e-3;10,10,300e-3);
solid cube2=orthobrick(-10,-10,130e-3;10,10,194.8e-3);
solid c1=main1in and cube2;
solid c2=main1out and cube;
solid MainCoil1= c2 and not c1 -bc=2;


#tlo c2 -col=[1,0,0];
tlo main1out -col=[0,1,0];
tlo c1 -col=[0,0,1];








