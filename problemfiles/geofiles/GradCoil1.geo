algebraic3d

# Z Gradient coil 1

solid grad1out = cylinder (0, 0, -130.3e-3; 0,0,180.7e-3;194.8e-3) -bco=1;
solid grad1in = cylinder (0, 0, -130.4e-3;0,0,180.8e-3; 189.4e-3) -bco=2;

solid cube5=orthobrick(-10,-10,-130.4e-3;10,10,180.7e-3);
solid cube6=orthobrick(-10,-10,-130.5e-3;10,10,180.8e-3);
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

tlo GradientCoil7 -col=[1,0,0];
