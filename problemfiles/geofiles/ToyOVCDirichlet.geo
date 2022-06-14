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


# Z Gradient coil 1

solid grad1out = cylinder (0, 0, 130.3e-3; 0,0,180.7e-3;194.8e-3) -bco=1;
solid grad1in = cylinder (0, 0, 130.3e-3;0,0,180.8e-3; 189.4e-3) -bco=2;

solid cube5=orthobrick(-10,-10,130.3e-3;10,10,180.7e-3);
solid cube6=orthobrick(-10,-10,128e-3;10,10,180.8e-3);
solid c5=grad1in and cube6;
solid c6=grad1out and cube5;
solid GradientCoil1= c6 and not c5 -bc=2;


tlo GradientCoil1 -col=[0,0,1];

# Z Gradient coil 2

solid grad2out = cylinder (0, 0, -180.7e-3; 0,0,-130.3e-3;194.8e-3) -bco=1;
solid grad2in = cylinder (0, 0, -180.8e-3;0,0,-130.3e-3; 189.4e-3) -bco=2;

solid cube7=orthobrick(-10,-10,-180.7e-3;10,10,-130.3e-3);
solid cube8=orthobrick(-10,-10,-180.8e-3;10,10,-128e-3);
solid c7=grad2in and cube8;
solid c8=grad2out and cube7;
solid GradientCoil2= c8 and not c7 -bc=2;


tlo GradientCoil2 -col=[0,0,1];


# OVC

solid OVCout = cylinder (0, 0, -250e-3; 0,0,250e-3;211.5e-3);
solid OVCin = cylinder (0, 0, -250.1e-3;0,0,250e-3; 209.5e-3);

solid cube13=orthobrick(-10,-10,-250e-3;10,10,250e-3);
solid cube14=orthobrick(-10,-10,-250.1e-3;10,10,250.1e-3);
solid c13=OVCin and cube14;
solid c14=OVCout and cube13;
solid OVCLeft= c14 and not c13 -bc=3;

#OVC Shield Middle
solid shieldOVCoutM = cylinder (0, 0, -5e-3; 0,0,5e-3;211.5e-3);
solid shieldOVCinM = cylinder (0, 0, -5.1e-3;0,0,5e-3; 209.5e-3);
solid cube13M=orthobrick(-10,-10,-5e-3;10,10,5e-3);
solid cube14M=orthobrick(-10,-10,-5.1e-3;10,10,5.1e-3);
solid c13M=shieldOVCinM and cube14M;
solid c14M=shieldOVCoutM and cube13M;
solid ShieldOVCM= c14M and not c13M -bc=4;

#OVC Shield Middle thicker
solid shieldOVCoutMT = cylinder (0, 0, -5e-3; 0,0,5e-3;211.9e-3);
solid shieldOVCinMT = cylinder (0, 0, -5.1e-3;0,0,5e-3; 209.1e-3);
solid cube13MT=orthobrick(-10,-10,-5e-3;10,10,5e-3);
solid cube14MT=orthobrick(-10,-10,-5e-3;10,10,5e-3);
solid c13MT=shieldOVCinMT and cube14MT;
solid c14MT=shieldOVCoutMT and cube13MT;
solid ShieldOVCMT= c14MT and not c13MT;

solid OVC= OVCLeft and not ShieldOVCMT;



tlo ShieldOVCM -col=[1,1,1];
tlo OVC -col=[1,1,0];

# OVC2

solid OVCout2 = cylinder (0, 0, -250e-3; 0,0,250e-3;214.5e-3);
solid OVCin2 = cylinder (0, 0, -250.1e-3;0,0,250e-3; 211.5e-3);

solid cube132=orthobrick(-10,-10,-250e-3;10,10,250e-3);
solid cube142=orthobrick(-10,-10,-250.1e-3;10,10,250.1e-3);
solid c132=OVCin2 and cube142;
solid c142=OVCout2 and cube132;
solid OVCLeft2= c142 and not c132 -bc=3;

#OVC Shield Middle
solid shieldOVCoutM2 = cylinder (0, 0, -5e-3; 0,0,5e-3;214.5e-3);
solid shieldOVCinM2 = cylinder (0, 0, -5.1e-3;0,0,5e-3; 211.5e-3);
solid cube13M2=orthobrick(-10,-10,-5e-3;10,10,5e-3);
solid cube14M2=orthobrick(-10,-10,-5.1e-3;10,10,5.1e-3);
solid c13M2=shieldOVCinM2 and cube14M2;
solid c14M2=shieldOVCoutM2 and cube13M2;
solid ShieldOVCM2= c14M2 and not c13M2 -bc=4;

#OVC Shield Middle thicker
solid shieldOVCoutMT2 = cylinder (0, 0, -5e-3; 0,0,5e-3;214.9e-3);
solid shieldOVCinMT2 = cylinder (0, 0, -5.1e-3;0,0,5e-3; 211.1e-3);
solid cube13MT2=orthobrick(-10,-10,-5e-3;10,10,5e-3);
solid cube14MT2=orthobrick(-10,-10,-5e-3;10,10,5e-3);
solid c13MT2=shieldOVCinMT2 and cube14MT2;
solid c14MT2=shieldOVCoutMT2 and cube13MT2;
solid ShieldOVCMT2= c14MT2 and not c13MT2;

solid OVC2= OVCLeft2 and not ShieldOVCMT2;



tlo ShieldOVCM2 -col=[1,1,1];
tlo OVC2 -col=[1,1,0];



# Air

solid aircylinderhelp = cylinder(0,0,-1200e-3;0,0,1200e-3;900e-3);
solid cube15 = orthobrick(-1,-1,-1200e-3;1,1,1200.1e-3);
solid aircylinder=aircylinderhelp and cube15;
solid aircylinder2=((((((((aircylinder and not 
MainCoil1)and not MainCoil2)and not GradientCoil1)and not GradientCoil2) and not ShieldOVCM)and not OVC)and not ShieldOVCM2)and not OVC2);

tlo aircylinder2 -col=[0,0,1] -transparent;

