algebraic3d

# Geometry for Toy problem with one shield and using symmetry



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


#OVC Shield Middle thicker
solid shield4KoutMT = cylinder (0, 0, -5e-3; 0,0,5e-3;243.1e-3);
solid shield4KinMT = cylinder (0, 0, -5.1e-3;0,0,5e-3; 239.9e-3);
solid cube13MT=orthobrick(-10,-10,-5e-3;10,10,5e-3);
solid cube14MT=orthobrick(-10,-10,-5e-3;10,10,5e-3);
solid c13MT=shield4KinMT and cube14MT;
solid c14MT=shield4KoutMT and cube13MT;
solid Shield4KMT= c14MT and not c13MT;

solid Shield4K= Shield4KLeft and not Shield4KMT;





tlo Shield4KM -col=[1,1,1];
tlo Shield4K -col=[1,1,0];

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


#77K Shield Middle thicker
solid shield77KoutMT = cylinder (0, 0, -5e-3; 0,0,5e-3;230.6e-3);
solid shield77KinMT = cylinder (0, 0, -5.1e-3;0,0,5e-3; 225.4e-3);
solid cube15MT=orthobrick(-10,-10,-5e-3;10,10,5e-3);
solid cube16MT=orthobrick(-10,-10,-5e-3;10,10,5e-3);
solid c15MT=shield77KinMT and cube16MT;
solid c16MT=shield77KoutMT and cube15MT;
solid Shield77KMT= c16MT and not c15MT;

solid Shield77K= Shield77KLeft and not Shield77KMT;




tlo Shield77KM -col=[1,1,1];
tlo Shield77K -col=[1,1,0];

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


#OVC Shield Middle thicker
solid shieldOVCoutMT = cylinder (0, 0, -5e-3; 0,0,5e-3;214.9e-3);
solid shieldOVCinMT = cylinder (0, 0, -5.1e-3;0,0,5e-3; 209.1e-3);
solid cube17MT=orthobrick(-10,-10,-5e-3;10,10,5e-3);
solid cube18MT=orthobrick(-10,-10,-5e-3;10,10,5e-3);
solid c17MT=shieldOVCinMT and cube14MT;
solid c18MT=shieldOVCoutMT and cube13MT;
solid ShieldOVCMT= c18MT and not c17MT;

solid OVC= OVCLeft and not ShieldOVCMT;




tlo ShieldOVCM -col=[1,1,1];
tlo OVC -col=[1,1,0];

# Y Gradient coil 1

solid grad1out_xy = cylinder (0, 0, -240e-3; 0,0,-40e-3;178.6e-3) -bco=1;
solid grad1in_xy = cylinder (0, 0, -240.1e-3;0,0,-39.9e-3; 173.2e-3) -bco=2;

solid cube5g=orthobrick(-10,-10,-240.1e-3;10,10,-40e-3);
solid cube6g=orthobrick(-10,-10,-240.2e-3;10,10,-39.9e-3);
solid c5g=grad1in_xy and cube6g;
solid c6g=grad1out_xy and cube5g;
solid GradientCoil1_xy= c6g and not c5g -bc=2;
solid plane1=plane(-0.1018,0.1401,150e-3;0.809017,0.58778,0);
solid GradientCoil2_xy=GradientCoil1_xy and not plane1;
solid plane2=plane(0.1018,0.1401,150e-3;-0.809017,0.58778,0);
solid GradientCoil3_xy=GradientCoil2_xy and not plane2;

solid grad2out_xy = cylinder (0, 0, -189.6e-3; 0,0,-90.4e-3;178.7e-3) -bco=1;
solid grad2in_xy = cylinder (0, 0, -189.7e-3;0,0,-90.3e-3; 173.1e-3) -bco=2;

solid cube7g=orthobrick(-10,-10,-189.7e-3;10,10,-90.4e-3);
solid cube8g=orthobrick(-10,-10,-189.8e-3;10,10,-90.3e-3);
solid c7g=grad2in_xy and cube8g;
solid c8g=grad2out_xy and cube7g;
solid GradientCoil4_xy= c8g and not c7g -bc=2;
solid plane3=plane(-0.05924,0.16275,150e-3;0.9397,0.342,0);
solid GradientCoil5_xy=GradientCoil4_xy and not plane3;
solid plane4=plane(0.05924,0.16275,150e-3;-0.9397,0.342,0);
solid GradientCoil6_xy=GradientCoil5_xy and not plane4;
solid GradientCoil7= GradientCoil3_xy and not GradientCoil6_xy;

tlo GradientCoil7 -col=[1,0,0];

# Y Gradient coil 2

solid grad1out_xy2 = cylinder (0, 0, 40e-3; 0,0,240e-3;178.6e-3) -bco=1;
solid grad1in_xy2 = cylinder (0, 0, 39.9e-3;0,0,240.1e-3; 173.2e-3) -bco=2;

solid cube5g2=orthobrick(-10,-10,39.9e-3;10,10,240e-3);
solid cube6g2=orthobrick(-10,-10,39.8e-3;10,10,240.1e-3);
solid c5g2=grad1in_xy2 and cube6g2;
solid c6g2=grad1out_xy2 and cube5g2;
solid GradientCoil1_xy2= c6g2 and not c5g2 -bc=2;
solid plane1=plane(-0.1018,0.1401,150e-3;0.809017,0.58778,0);
solid GradientCoil2_xy2=GradientCoil1_xy2 and not plane1;
solid plane2=plane(0.1018,0.1401,150e-3;-0.809017,0.58778,0);
solid GradientCoil3_xy2=GradientCoil2_xy2 and not plane2;

solid grad2out_xy2 = cylinder (0, 0, 90.4e-3; 0,0,189.6e-3;178.7e-3) -bco=1;
solid grad2in_xy2 = cylinder (0, 0, 90.3e-3;0,0,189.7e-3; 173.1e-3) -bco=2;

solid cube7g2=orthobrick(-10,-10,90.3e-3;10,10,189.6e-3);
solid cube8g2=orthobrick(-10,-10,90.2e-3;10,10,189.7e-3);
solid c7g2=grad2in_xy2 and cube8g2;
solid c8g2=grad2out_xy2 and cube7g2;
solid GradientCoil4_xy2= c8g2 and not c7g2 -bc=2;
solid plane3=plane(-0.05924,0.16275,150e-3;0.9397,0.342,0);
solid GradientCoil5_xy2=GradientCoil4_xy2 and not plane3;
solid plane4=plane(0.05924,0.16275,150e-3;-0.9397,0.342,0);
solid GradientCoil6_xy2=GradientCoil5_xy2 and not plane4;
solid GradientCoil8= GradientCoil3_xy2 and not GradientCoil6_xy2;

tlo GradientCoil8 -col=[1,0,0];

# Y Gradient coil 3

solid grad1out_xy3 = cylinder (0, 0, -240e-3; 0,0,-40e-3;178.6e-3) -bco=1;
solid grad1in_xy3 = cylinder (0, 0, -240.1e-3;0,0,-39.9e-3; 173.2e-3) -bco=2;

solid cube5g3=orthobrick(-10,-10,-240.1e-3;10,10,-40e-3);
solid cube6g3=orthobrick(-10,-10,-240.2e-3;10,10,-39.9e-3);
solid c5g3=grad1in_xy3 and cube6g3;
solid c6g3=grad1out_xy3 and cube5g3;
solid GradientCoil1_xy3= c6g3 and not c5g3 -bc=2;
solid plane5=plane(-0.1018,-0.1401,150e-3;0.809017,-0.58778,0);
solid GradientCoil2_xy3=GradientCoil1_xy3 and not plane5;
solid plane6=plane(0.1018,-0.1401,150e-3;-0.809017,-0.58778,0);
solid GradientCoil3_xy3=GradientCoil2_xy3 and not plane6;

solid grad2out_xy3 = cylinder (0, 0, -189.6e-3; 0,0,-90.4e-3;178.7e-3) -bco=1;
solid grad2in_xy3 = cylinder (0, 0, -189.7e-3;0,0,-90.3e-3; 173.1e-3) -bco=2;

solid cube7g3=orthobrick(-10,-10,-189.7e-3;10,10,-90.4e-3);
solid cube8g3=orthobrick(-10,-10,-189.8e-3;10,10,-90.3e-3);
solid c7g3=grad2in_xy3 and cube8g3;
solid c8g3=grad2out_xy3 and cube7g3;
solid GradientCoil4_xy3= c8g3 and not c7g3 -bc=2;
solid plane7=plane(-0.05924,-0.16275,150e-3;0.9397,-0.342,0);
solid GradientCoil5_xy3=GradientCoil4_xy3 and not plane7;
solid plane8=plane(0.05924,-0.16275,150e-3;-0.9397,-0.342,0);
solid GradientCoil6_xy3=GradientCoil5_xy3 and not plane8;
solid GradientCoil9= GradientCoil3_xy3 and not GradientCoil6_xy3;

tlo GradientCoil9 -col=[1,0,0];

# Y Gradient coil 4

solid grad1out_xy4 = cylinder (0, 0, 40e-3; 0,0,240e-3;178.6e-3) -bco=1;
solid grad1in_xy4 = cylinder (0, 0, 39.9e-3;0,0,240.1e-3; 173.2e-3) -bco=2;

solid cube5g4=orthobrick(-10,-10,39.9e-3;10,10,240e-3);
solid cube6g4=orthobrick(-10,-10,39.8e-3;10,10,240.1e-3);
solid c5g4=grad1in_xy4 and cube6g4;
solid c6g4=grad1out_xy4 and cube5g4;
solid GradientCoil1_xy4= c6g4 and not c5g4 -bc=2;
solid plane5=plane(-0.1018,-0.1401,150e-3;0.809017,-0.58778,0);
solid GradientCoil2_xy4=GradientCoil1_xy4 and not plane5;
solid plane6=plane(0.1018,-0.1401,150e-3;-0.809017,-0.58778,0);
solid GradientCoil3_xy4=GradientCoil2_xy4 and not plane6;

solid grad2out_xy4 = cylinder (0, 0, 90.4e-3; 0,0,189.6e-3;178.7e-3) -bco=1;
solid grad2in_xy4= cylinder (0, 0, 90.3e-3;0,0,189.7e-3; 173.1e-3) -bco=2;

solid cube7g4=orthobrick(-10,-10,90.3e-3;10,10,189.6e-3);
solid cube8g4=orthobrick(-10,-10,90.2e-3;10,10,189.7e-3);
solid c7g4=grad2in_xy4 and cube8g4;
solid c8g4=grad2out_xy4 and cube7g4;
solid GradientCoil4_xy4= c8g4 and not c7g4 -bc=2;
solid plane7=plane(-0.05924,-0.16275,150e-3;0.9397,-0.342,0);
solid GradientCoil5_xy4=GradientCoil4_xy4 and not plane7;
solid plane8=plane(0.05924,-0.16275,150e-3;-0.9397,-0.342,0);
solid GradientCoil6_xy4=GradientCoil5_xy4 and not plane8;
solid GradientCoil10= GradientCoil3_xy4 and not GradientCoil6_xy4;

tlo GradientCoil10 -col=[1,0,0];

# X Gradient coil 1

solid grad1out_x = cylinder (0, 0, -240e-3; 0,0,-40e-3;178.6e-3) -bco=1;
solid grad1in_x = cylinder (0, 0, -240.1e-3;0,0,-39.9e-3; 173.2e-3) -bco=2;

solid cube5gx=orthobrick(-10,-10,-240.1e-3;10,10,-40e-3);
solid cube6gx=orthobrick(-10,-10,-240.2e-3;10,10,-39.9e-3);
solid c5gx=grad1in_x and cube6gx;
solid c6gx=grad1out_x and cube5gx;
solid GradientCoil1_x= c6gx and not c5gx -bc=2;
solid plane1x=plane(0.1401,-0.1018,150e-3;0.58778,0.809017,0);
solid GradientCoil2_x=GradientCoil1_x and not plane1x;
solid plane2x=plane(0.1401,0.1018,150e-3;0.58778,-0.809017,0);
solid GradientCoil3_x=GradientCoil2_x and not plane2x;

solid grad2out_x = cylinder (0, 0, -189.6e-3; 0,0,-90.4e-3;178.7e-3) -bco=1;
solid grad2in_x = cylinder (0, 0, -189.7e-3;0,0,-90.3e-3; 173.1e-3) -bco=2;

solid cube7gx=orthobrick(-10,-10,-189.7e-3;10,10,-90.4e-3);
solid cube8gx=orthobrick(-10,-10,-189.8e-3;10,10,-90.3e-3);
solid c7gx=grad2in_x and cube8gx;
solid c8gx=grad2out_x and cube7gx;
solid GradientCoil4_x= c8gx and not c7gx -bc=2;
solid plane3x=plane(0.16275,-0.05924,150e-3;0.342,0.9397,0);
solid GradientCoil5_x=GradientCoil4_x and not plane3x;
solid plane4x=plane(0.16275,0.05924,150e-3;0.342,-0.9397,0);
solid GradientCoil6_x=GradientCoil5_x and not plane4x;
solid GradientCoil7x= GradientCoil3_x and not GradientCoil6_x;

tlo GradientCoil7x -col=[1,0,0];

# X Gradient coil 2

solid grad1out_x2 = cylinder (0, 0, 40e-3; 0,0,240e-3;178.6e-3) -bco=1;
solid grad1in_x2 = cylinder (0, 0, 39.9e-3;0,0,240.1e-3; 173.2e-3) -bco=2;

solid cube5g2x=orthobrick(-10,-10,39.9e-3;10,10,240e-3);
solid cube6g2x=orthobrick(-10,-10,39.8e-3;10,10,240.1e-3);
solid c5g2x=grad1in_x2 and cube6g2x;
solid c6g2x=grad1out_x2 and cube5g2x;
solid GradientCoil1_x2= c6g2x and not c5g2x -bc=2;
solid plane1x=plane(0.1401,-0.1018,150e-3;0.58778,0.809017,0);
solid GradientCoil2_x2=GradientCoil1_x2 and not plane1x;
solid plane2x=plane(0.1401,0.1018,150e-3;0.58778,-0.809017,0);
solid GradientCoil3_x2=GradientCoil2_x2 and not plane2x;

solid grad2out_x2 = cylinder (0, 0, 90.4e-3; 0,0,189.6e-3;178.7e-3) -bco=1;
solid grad2in_x2 = cylinder (0, 0, 90.3e-3;0,0,189.7e-3; 173.1e-3) -bco=2;

solid cube7g2x=orthobrick(-10,-10,90.3e-3;10,10,189.6e-3);
solid cube8g2x=orthobrick(-10,-10,90.2e-3;10,10,189.7e-3);
solid c7g2x=grad2in_x2 and cube8g2x;
solid c8g2x=grad2out_x2 and cube7g2x;
solid GradientCoil4_x2= c8g2x and not c7g2x -bc=2;
solid plane3x=plane(0.16275,-0.05924,150e-3;0.342,0.9397,0);
solid GradientCoil5_x2=GradientCoil4_x2 and not plane3x;
solid plane4x=plane(0.16275,0.05924,150e-3;0.342,-0.9397,0);
solid GradientCoil6_x2=GradientCoil5_x2 and not plane4x;
solid GradientCoil8x= GradientCoil3_x2 and not GradientCoil6_x2;

tlo GradientCoil8x -col=[1,0,0];

# X Gradient coil 3

solid grad1out_x3 = cylinder (0, 0, -240e-3; 0,0,-40e-3;178.6e-3) -bco=1;
solid grad1in_x3 = cylinder (0, 0, -240.1e-3;0,0,-39.9e-3; 173.2e-3) -bco=2;

solid cube5gx3=orthobrick(-10,-10,-240.1e-3;10,10,-40e-3);
solid cube6gx3=orthobrick(-10,-10,-240.2e-3;10,10,-39.9e-3);
solid c5gx3=grad1in_x3 and cube6gx3;
solid c6gx3=grad1out_x3 and cube5gx3;
solid GradientCoil1_x3= c6gx3 and not c5gx3 -bc=2;
solid plane1x3=plane(-0.1401,-0.1018,150e-3;-0.58778,0.809017,0);
solid GradientCoil2_x3=GradientCoil1_x and not plane1x3;
solid plane2x3=plane(-0.1401,0.1018,150e-3;-0.58778,-0.809017,0);
solid GradientCoil3_x3=GradientCoil2_x3 and not plane2x3;

solid grad2out_x3 = cylinder (0, 0, -189.6e-3; 0,0,-90.4e-3;178.7e-3) -bco=1;
solid grad2in_x3 = cylinder (0, 0, -189.7e-3;0,0,-90.3e-3; 173.1e-3) -bco=2;

solid cube7gx3=orthobrick(-10,-10,-189.7e-3;10,10,-90.4e-3);
solid cube8gx3=orthobrick(-10,-10,-189.8e-3;10,10,-90.3e-3);
solid c7gx3=grad2in_x3 and cube8gx3;
solid c8gx3=grad2out_x3 and cube7gx3;
solid GradientCoil4_x3= c8gx3 and not c7gx3 -bc=2;
solid plane3x3=plane(-0.16275,-0.05924,150e-3;-0.342,0.9397,0);
solid GradientCoil5_x3=GradientCoil4_x and not plane3x3;
solid plane4x3=plane(-0.16275,0.05924,150e-3;-0.342,-0.9397,0);
solid GradientCoil6_x3=GradientCoil5_x3 and not plane4x3;
solid GradientCoil7x3= GradientCoil3_x3 and not GradientCoil6_x3;

tlo GradientCoil7x3 -col=[1,0,0];

# X Gradient coil 4

solid grad1out_x4 = cylinder (0, 0, 40e-3; 0,0,240e-3;178.6e-3) -bco=1;
solid grad1in_x4 = cylinder (0, 0, 39.9e-3;0,0,240.1e-3; 173.2e-3) -bco=2;

solid cube5g2x4=orthobrick(-10,-10,39.9e-3;10,10,240e-3);
solid cube6g2x4=orthobrick(-10,-10,39.8e-3;10,10,240.1e-3);
solid c5g2x4=grad1in_x4 and cube6g2x4;
solid c6g2x4=grad1out_x4 and cube5g2x4;
solid GradientCoil1_x4= c6g2x4 and not c5g2x4 -bc=2;
solid plane1x4=plane(-0.1401,-0.1018,150e-3;-0.58778,0.809017,0);
solid GradientCoil2_x4=GradientCoil1_x4 and not plane1x4;
solid plane2x4=plane(-0.1401,0.1018,150e-3;-0.58778,-0.809017,0);
solid GradientCoil3_x4=GradientCoil2_x4 and not plane2x4;

solid grad2out_x4 = cylinder (0, 0, 90.4e-3; 0,0,189.6e-3;178.7e-3) -bco=1;
solid grad2in_x4 = cylinder (0, 0, 90.3e-3;0,0,189.7e-3; 173.1e-3) -bco=2;

solid cube7g2x4=orthobrick(-10,-10,90.3e-3;10,10,189.6e-3);
solid cube8g2x4=orthobrick(-10,-10,90.2e-3;10,10,189.7e-3);
solid c7g2x4=grad2in_x4 and cube8g2x4;
solid c8g2x4=grad2out_x4 and cube7g2x4;
solid GradientCoil4_x4= c8g2x4 and not c7g2x4 -bc=2;
solid plane3x4=plane(-0.16275,-0.05924,150e-3;-0.342,0.9397,0);
solid GradientCoil5_x4=GradientCoil4_x4 and not plane3x4;
solid plane4x4=plane(-0.16275,0.05924,150e-3;-0.342,-0.9397,0);
solid GradientCoil6_x4=GradientCoil5_x4 and not plane4x4;
solid GradientCoil8x4= GradientCoil3_x4 and not GradientCoil6_x4;

tlo GradientCoil8x4 -col=[1,0,0];
# Extra shield

solid Extraout = cylinder (0, 0, -250e-3; 0,0,250e-3;205e-3);
solid Extrain = cylinder (0, 0, -250.1e-3;0,0,250e-3; 199e-3);


solid c19=Extrain and cube18;
solid c20=Extraout and cube17;
solid ExtraLeft= c20 and not c19 -bc=3;





tlo ExtraLeft -col=[1,1,0];

# Air

solid aircylinderhelp = cylinder(0,0,-1200e-3;0,0,1200e-3;900e-3);
solid cube15 = orthobrick(-1,-1,-1200e-3;1,1,1200.1e-3);
solid aircylinder=aircylinderhelp and cube15;
solid aircylinder2=(((((((((((((((((((aircylinder and not 
MainCoil1)and not MainCoil2)and not GradientCoil1)and not GradientCoil2) and not ShieldOVCM)and not OVC)and not Shield4KM) and not Shield4K) and not Shield77KM) and not Shield77K)and not GradientCoil7)and not GradientCoil8)and not GradientCoil9)and not GradientCoil10)and not GradientCoil7x)and not GradientCoil8x)and not GradientCoil7x3)and not GradientCoil8x4)and not ExtraLeft);



tlo aircylinder2 -col=[0,0,1] -transparent;

