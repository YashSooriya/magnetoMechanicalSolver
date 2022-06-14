algebraic3d

# Geometry for Toy problem

# Main coil 1

solid main1out = cylinder (0, 0, 135e-3; 0,0,194.7e-3;288.4e-3);
solid main1in = cylinder (0, 0, 135e-3;0,0,194.8e-3; 250e-3);

solid cube=orthobrick(-10,-10,135e-3;10,10,194.7e-3);
solid cube2=orthobrick(-10,-10,130e-3;10,10,194.8e-3);
solid c1=main1in and cube2;
solid c2=main1out and cube;
solid MainCoil1= c2 and not c1 -bc=2;


tlo MainCoil1 -col=[1,0,0];

# Main coil 2

solid main2out = cylinder (0, 0, -194.7e-3; 0,0,-135e-3;288.4e-3);
solid main2in = cylinder (0, 0, -194.8e-3;0,0,-135e-3; 250e-3);

solid cube3=orthobrick(-10,-10,-194.7e-3;10,10,-135e-3);
solid cube4=orthobrick(-10,-10,-194.8e-3;10,10,-130e-3);
solid c3=main2in and cube4;
solid c4=main2out and cube3;
solid MainCoil2= c4 and not c3 -bc=2;

tlo MainCoil2 -col=[1,0,0];





# 4K Shield
solid shield4outLeft = cylinder (0, 0, -250e-3; 0,0,250e-3;243e-3);
solid shield4inLeft = cylinder (0, 0, -250.1e-3;0,0,250.1e-3; 240e-3);
solid cube9=orthobrick(-10,-10,-250e-3;10,10,250e-3);
solid cube10=orthobrick(-10,-10,-250.1e-3;10,10,250.1e-3);
solid c9=shield4inLeft and cube10;
solid c10=shield4outLeft and cube9;
solid Shield4KLeft= c10 and not c9 -bc=3;


#4K Shield Middle
solid shield4outM = cylinder (0, 0, -5e-3; 0,0,5e-3;243e-3);
solid shield4inM = cylinder (0, 0, -5.1e-3;0,0,5e-3; 240e-3);
solid cube9M=orthobrick(-10,-10,-5e-3;10,10,5e-3);
solid cube10M=orthobrick(-10,-10,-5e-3;10,10,5e-3);
solid c9M=shield4inM and cube10M;
solid c10M=shield4outM and cube9M;
solid Shield4KM= c10M and not c9M -bc=4;

solid Shield4K= Shield4KLeft and not Shield4KM;


tlo Shield4KM -col=[1,1,1];
tlo Shield4K -col=[1,1,0];





# Air

solid aircylinderhelp = cylinder(0,0,-1200e-3;0,0,1200e-3;900e-3);
solid cube15 = orthobrick(-1,-1,-1200e-3;1,1,1200.1e-3);
solid aircylinder=aircylinderhelp and cube15;
solid aircylinder2=((((aircylinder and not 
MainCoil1)and not MainCoil2) and not Shield4K) and not Shield4KM);

tlo aircylinder2 -col=[0,0,1] -transparent;

