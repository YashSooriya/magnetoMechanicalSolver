algebraic3d

# Geometry for Toy problem with one shield and using symmetry





# Y Gradient coil 1

solid grad1out_xy = cylinder (0, 0, -240e-3; 0,0,-40e-3;178.6e-3) -bco=1;
solid grad1in_xy = cylinder (0, 0, -240.1e-3;0,0,-39.9e-3; 173.2e-3) -bco=2;

solid cube5g=orthobrick(-10,-10,-240.1e-3;10,10,-40e-3);
solid cube6g=orthobrick(-10,-10,-240.2e-3;10,10,-39.9e-3);
solid c5g=grad1in_xy and cube6g;
solid c6g=grad1out_xy and cube5g;
solid GradientCoil1_xy= c6g and not c5g -bc=2;
solid plane1=plane(-0.1018,0.1401,150e-3;0.86602,0.5,0);
solid GradientCoil2_xy=GradientCoil1_xy and not plane1;
solid plane2=plane(0.1018,0.1401,150e-3;-0.86602,0.5,0);
solid GradientCoil3_xy=GradientCoil2_xy and not plane2;

solid grad2out_xy = cylinder (0, 0, -189.6e-3; 0,0,-90.4e-3;178.7e-3) -bco=1;
solid grad2in_xy = cylinder (0, 0, -189.7e-3;0,0,-90.3e-3; 173.1e-3) -bco=2;

solid cube7g=orthobrick(-10,-10,-189.7e-3;10,10,-90.4e-3);
solid cube8g=orthobrick(-10,-10,-189.8e-3;10,10,-90.3e-3);
solid c7g=grad2in_xy and cube8g;
solid c8g=grad2out_xy and cube7g;
solid GradientCoil4_xy= c8g and not c7g -bc=2;
solid plane3=plane(-0.05924,0.16275,150e-3;0.86602,0.5,0);
solid GradientCoil5_xy=GradientCoil4_xy and not plane3;
solid plane4=plane(0.05924,0.16275,150e-3;-0.86602,0.5,0);
solid GradientCoil6_xy=GradientCoil5_xy and not plane4;
solid GradientCoil7= GradientCoil3_xy and not GradientCoil6_xy;

tlo GradientCoil7 -col=[1,0,0];





