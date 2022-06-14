algebraic3d

# Z Gradient coil 1

solid grad1out = cylinder (0, 0, -135.3e-3; 0,0,185.7e-3;194.8e-3) -bco=1;
solid grad1in = cylinder (0, 0, -135.4e-3;0,0,185.8e-3; 189.4e-3) -bco=2;

solid cube5=orthobrick(-10,-10,-135.4e-3;10,10,185.7e-3);
solid cube6=orthobrick(-10,-10,-135.5e-3;10,10,185.8e-3);
solid c5=grad1in and cube6;
solid c6=grad1out and cube5;
solid GradientCoil1= c6 and not c5 -bc=2;
solid plane1=plane(-0.09,0.15588457,150e-3;0.86602,0.5,0);
solid GradientCoil2=GradientCoil1 and not plane1;
solid plane2=plane(0.09,0.15588457,150e-3;-0.86602,0.5,0);
solid GradientCoil3=GradientCoil2 and not plane2;

solid grad2out = cylinder (0, 0, -90.3e-3; 0,0,140.7e-3;194.9e-3) -bco=1;
solid grad2in = cylinder (0, 0, -90.4e-3;0,0,140.8e-3; 189.3e-3) -bco=2;

solid cube7=orthobrick(-10,-10,-90.4e-3;10,10,140.7e-3);
solid cube8=orthobrick(-10,-10,-90.5e-3;10,10,140.8e-3);
solid c7=grad2in and cube8;
solid c8=grad2out and cube7;
solid GradientCoil4= c8 and not c7 -bc=2;
solid plane3=plane(-0.046587,0.173867,150e-3;0.86602,0.5,0);
solid GradientCoil5=GradientCoil4 and not plane3;
solid plane4=plane(0.046587,0.173867,150e-3;-0.86602,0.5,0);
solid GradientCoil6=GradientCoil5 and not plane4;
solid GradientCoil7= GradientCoil3 and not GradientCoil6;


#======================================================================================
solid grad1out_x = cylinder (0, 0, -135.4e-3; 0,0,185.8e-3;194.9e-3) -bco=1;
solid grad1in_x = cylinder (0, 0, -135.5e-3;0,0,185.9e-3; 189.3e-3) -bco=2;

solid cube5_x=orthobrick(-10,-10,-135.4e-3;10,10,185.8e-3);
solid cube6_x=orthobrick(-10,-10,-135.5e-3;10,10,185.9e-3);
solid c5_x=grad1in_x and cube6_x;
solid c6_x=grad1out_x and cube5_x;
solid GradientCoil1_x= c6_x and not c5_x -bc=2;
solid plane2_x=plane(0.10511,0.14612,150e-3;-0.86602,0.5,0);
solid GradientCoil2_x=GradientCoil1_x and not plane1;
solid GradientCoil3_x=GradientCoil2_x and not plane2_x;

solid grad2out_x = cylinder (0, 0, -90.3e-3; 0,0,140.7e-3;195e-3) -bco=1;
solid grad2in_x = cylinder (0, 0, -90.4e-3;0,0,140.8e-3; 189.2e-3) -bco=2;

solid c7_x=grad2in_x and cube8;
solid c8_x=grad2out_x and cube7;
solid GradientCoil4_x= c8_x and not c7_x -bc=2;
solid GradientCoil5_x=GradientCoil4_x and not plane3;
solid plane4_2=plane(0.046587,0.173867,150e-3;0.86602,-0.5,0);
solid GradientCoil6_x=GradientCoil5_x and not plane4;
solid GradientCoil7_x= GradientCoil3_x and not GradientCoil6_x;


solid GradientCoil8=GradientCoil7_x and not plane4_2;
solid GradientCoil8_2=GradientCoil7 and not plane4_2;
solid plane5_2=plane(0.046588,0.173867,140.7e-3;0,0,1);
solid GradientCoil9=GradientCoil8 and not plane5_2;
solid GradientCoil9_2=GradientCoil8_2 and not plane5_2;
solid sphere1=cylinder(0.04902,-0.1829,140.7e-3;0.050418,0.18816,140.7e-3;0.0475);
solid AuxPiece=GradientCoil9_2 and sphere1;
solid GradientCoil10=GradientCoil7 and not GradientCoil9;
solid GradientCoil11=GradientCoil10 or AuxPiece;

#tlo GradientCoil3 -col=[1,0,0];
#tlo GradientCoil1 -col=[0,1,0];
tlo GradientCoil11 -col=[1,0,0];





