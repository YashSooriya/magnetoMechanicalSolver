algebraic3d

# Geometry for Toy problem with one shield and using symmetry

# Build cubes for symmetry cuts

solid SymCube1 = orthobrick(0,-10,-10;10,10,10);
solid SymCube2 = orthobrick(-10,0,-10;10,10,10);

# Main coil 1

solid main1out = cylinder (0, 0, 135e-3; 0,0,194.7e-3;288.4e-3);
solid main1in = cylinder (0, 0, 135e-3;0,0,194.8e-3; 250e-3);

solid cube=orthobrick(-10,-10,135e-3;10,10,194.7e-3);
solid cube2=orthobrick(-10,-10,130e-3;10,10,194.8e-3);
solid c1=main1in and cube2;
solid c2=main1out and cube;
solid MainCoil1= c2 and not c1 -bc=2;
solid MainCoil1Cut1 = MainCoil1 and not SymCube1;
solid MainCoil1Cut2 = MainCoil1Cut1 and not SymCube2;

tlo MainCoil1Cut2 -col=[1,0,0];

# Main coil 2

solid main2out = cylinder (0, 0, -194.7e-3; 0,0,-135e-3;288.4e-3);
solid main2in = cylinder (0, 0, -194.8e-3;0,0,-135e-3; 250e-3);

solid cube3=orthobrick(-10,-10,-194.7e-3;10,10,-135e-3);
solid cube4=orthobrick(-10,-10,-194.8e-3;10,10,-130e-3);
solid c3=main2in and cube4;
solid c4=main2out and cube3;
solid MainCoil2= c4 and not c3 -bc=2;
solid MainCoil2Cut1 = MainCoil2 and not SymCube1;
solid MainCoil2Cut2 = MainCoil2Cut1 and not SymCube2;

tlo MainCoil2Cut2 -col=[1,0,0];


# Z Gradient coil 1

solid grad1out = cylinder (0, 0, 130.3e-3; 0,0,180.7e-3;194.8e-3) -bco=1;
solid grad1in = cylinder (0, 0, 130.3e-3;0,0,180.8e-3; 189.4e-3) -bco=2;

solid cube5=orthobrick(-10,-10,130.3e-3;10,10,180.7e-3);
solid cube6=orthobrick(-10,-10,128e-3;10,10,180.8e-3);
solid c5=grad1in and cube6;
solid c6=grad1out and cube5;
solid GradientCoil1= c6 and not c5 -bc=2;
solid GradientCoil1Cut1 = GradientCoil1 and not SymCube1;
solid GradientCoil1Cut2 = GradientCoil1Cut1 and not SymCube2;

tlo GradientCoil1Cut2 -col=[0,0,1];

# Z Gradient coil 2

solid grad2out = cylinder (0, 0, -180.7e-3; 0,0,-130.3e-3;194.8e-3) -bco=1;
solid grad2in = cylinder (0, 0, -180.8e-3;0,0,-130.3e-3; 189.4e-3) -bco=2;

solid cube7=orthobrick(-10,-10,-180.7e-3;10,10,-130.3e-3);
solid cube8=orthobrick(-10,-10,-180.8e-3;10,10,-128e-3);
solid c7=grad2in and cube8;
solid c8=grad2out and cube7;
solid GradientCoil2= c8 and not c7 -bc=2;
solid GradientCoil2Cut1 = GradientCoil2 and not SymCube1;
solid GradientCoil2Cut2 = GradientCoil2Cut1 and not SymCube2;

tlo GradientCoil2Cut2 -col=[0,0,1];

# 4K

solid Shield4Kout = cylinder (0, 0, -250e-3; 0,0,250e-3;243e-3);
solid Shield4Kin = cylinder (0, 0, -250.1e-3;0,0,250e-3; 240e-3);

solid cube13=orthobrick(-10,-10,-250e-3;10,10,250e-3);
solid cube14=orthobrick(-10,-10,-250.1e-3;10,10,250.1e-3);
solid c13=Shield4Kin and cube14;
solid c14=Shield4Kout and cube13;
solid Shield4KLeft= c14 and not c13 -bc=3;

#OVC Shield Middle
solid shield4KoutM = cylinder (0, 0, -5e-3; 0,0,5e-3;243e-3);
solid shield4KinM = cylinder (0, 0, -5.1e-3;0,0,5e-3; 240e-3);
solid cube13M=orthobrick(-10,-10,-5e-3;10,10,5e-3);
solid cube14M=orthobrick(-10,-10,-5.1e-3;10,10,5.1e-3);
solid c13M=shield4KinM and cube14M;
solid c14M=shield4KoutM and cube13M;
solid Shield4KM= c14M and not c13M -bc=4;
solid Shield4Kcut1 = Shield4KM and not SymCube1;
solid Shield4Kcut2= Shield4Kcut1 and not SymCube2;

#OVC Shield Middle thicker
solid shield4KoutMT = cylinder (0, 0, -5e-3; 0,0,5e-3;243.1e-3);
solid shield4KinMT = cylinder (0, 0, -5.1e-3;0,0,5e-3; 239.9e-3);
solid cube13MT=orthobrick(-10,-10,-5e-3;10,10,5e-3);
solid cube14MT=orthobrick(-10,-10,-5e-3;10,10,5e-3);
solid c13MT=shield4KinMT and cube14MT;
solid c14MT=shield4KoutMT and cube13MT;
solid Shield4KMT= c14MT and not c13MT;

solid Shield4K= Shield4KLeft and not Shield4KMT;

solid Shield4K2cut1= Shield4K and not SymCube1;
solid Shield4K2cut2= Shield4K2cut1 and not SymCube2;



tlo Shield4Kcut2 -col=[1,1,1];
tlo Shield4K2cut2 -col=[1,1,0];

# 77K

solid Shield77Kout = cylinder (0, 0, -250e-3; 0,0,250e-3;230.5e-3);
solid Shield77Kin = cylinder (0, 0, -250.1e-3;0,0,250e-3; 225.5e-3);

solid c15=Shield77Kin and cube14;
solid c16=Shield77Kout and cube13;
solid Shield77KLeft= c16 and not c15 -bc=3;

#77K Shield Middle
solid shield77KoutM = cylinder (0, 0, -5e-3; 0,0,5e-3;230.5e-3);
solid shield77KinM = cylinder (0, 0, -5.1e-3;0,0,5e-3; 225.5e-3);
solid cube15M=orthobrick(-10,-10,-5e-3;10,10,5e-3);
solid cube16M=orthobrick(-10,-10,-5.1e-3;10,10,5.1e-3);
solid c15M=shield77KinM and cube16M;
solid c16M=shield77KoutM and cube15M;
solid Shield77KM= c16M and not c15M -bc=4;
solid Shield77Kcut1 = Shield77KM and not SymCube1;
solid Shield77Kcut2= Shield77Kcut1 and not SymCube2;

#77K Shield Middle thicker
solid shield77KoutMT = cylinder (0, 0, -5e-3; 0,0,5e-3;230.6e-3);
solid shield77KinMT = cylinder (0, 0, -5.1e-3;0,0,5e-3; 225.4e-3);
solid cube15MT=orthobrick(-10,-10,-5e-3;10,10,5e-3);
solid cube16MT=orthobrick(-10,-10,-5e-3;10,10,5e-3);
solid c15MT=shield77KinMT and cube16MT;
solid c16MT=shield77KoutMT and cube15MT;
solid Shield77KMT= c16MT and not c15MT;

solid Shield77K= Shield77KLeft and not Shield77KMT;

solid Shield77K2cut1= Shield77K and not SymCube1;
solid Shield77K2cut2= Shield77K2cut1 and not SymCube2;



tlo Shield77Kcut2 -col=[1,1,1];
tlo Shield77K2cut2 -col=[1,1,0];

# OVC

solid OVCout = cylinder (0, 0, -250e-3; 0,0,250e-3;214.5e-3);
solid OVCin = cylinder (0, 0, -250.1e-3;0,0,250e-3; 209.5e-3);

solid cube17=orthobrick(-10,-10,-250e-3;10,10,250e-3);
solid cube18=orthobrick(-10,-10,-250.1e-3;10,10,250.1e-3);
solid c17=OVCin and cube18;
solid c18=OVCout and cube17;
solid OVCLeft= c18 and not c17 -bc=3;

#OVC Shield Middle
solid shieldOVCoutM = cylinder (0, 0, -5e-3; 0,0,5e-3;214.5e-3);
solid shieldOVCinM = cylinder (0, 0, -5.1e-3;0,0,5e-3; 209.5e-3);
solid cube17M=orthobrick(-10,-10,-5e-3;10,10,5e-3);
solid cube18M=orthobrick(-10,-10,-5.1e-3;10,10,5.1e-3);
solid c17M=shieldOVCinM and cube18M;
solid c18M=shieldOVCoutM and cube17M;
solid ShieldOVCM= c18M and not c17M -bc=4;
solid OVCcut1 = ShieldOVCM and not SymCube1;
solid OVCcut2= OVCcut1 and not SymCube2;

#OVC Shield Middle thicker
solid shieldOVCoutMT = cylinder (0, 0, -5e-3; 0,0,5e-3;214.9e-3);
solid shieldOVCinMT = cylinder (0, 0, -5.1e-3;0,0,5e-3; 209.1e-3);
solid cube17MT=orthobrick(-10,-10,-5e-3;10,10,5e-3);
solid cube18MT=orthobrick(-10,-10,-5e-3;10,10,5e-3);
solid c17MT=shieldOVCinMT and cube14MT;
solid c18MT=shieldOVCoutMT and cube13MT;
solid ShieldOVCMT= c18MT and not c17MT;

solid OVC= OVCLeft and not ShieldOVCMT;

solid OVC2cut1= OVC and not SymCube1;
solid OVC2cut2= OVC2cut1 and not SymCube2;



tlo OVCcut2 -col=[1,1,1];
tlo OVC2cut2 -col=[1,1,0];





# Air

solid aircylinderhelp = cylinder(0,0,-1200e-3;0,0,1200e-3;900e-3);
solid cube15 = orthobrick(-1,-1,-1200e-3;1,1,1200.1e-3);
solid aircylinder=aircylinderhelp and cube15;
solid aircylinder2=((((((((((aircylinder and not 
MainCoil1Cut2)and not MainCoil2Cut2)and not GradientCoil1Cut2)and not GradientCoil2Cut2) and not OVCcut2)and not OVC2cut2)and not Shield4Kcut2) and not Shield4K2cut2) and not Shield77Kcut2) and not Shield77K2cut2);

solid exteriorCut1= aircylinder2 and not SymCube1;
solid exteriorCut2= exteriorCut1 and not SymCube2;

tlo exteriorCut2 -col=[0,0,1] -transparent;

