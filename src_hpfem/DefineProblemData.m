function [probdata]=DefineProblemData(problem,pm,option)

if problem==0
    % Simple plane wave propagation (for testing only)
    profun=@problem0;
    probdata = profun(pm);
elseif problem==1
    % unit conducing magnetic sphere, scaled by delta
    profun=@problem1;
    probdata = profun(pm,option);
elseif problem==2
    % ANSYS Testing unit 1 sphere, scaled by delta
    profun=@problem1ansys;
    probdata = profun(pm);
elseif problem==3
    % TEAM 7 problem
    profun=@problemTEAM7;
    probdata=profun(pm);
elseif problem==4
    % Cylinder
    profun=@problem_cylinder;
    probdata=profun(pm);
elseif problem==5
     % TEAM 7 problem test
    profun=@problemTEAM7b;
    probdata=profun(pm);
elseif problem==6
    % Shells
    profun=@problem_shells;
    probdata=profun(pm);
elseif problem==7
    % Torus
    profun=@problemTorus;
    probdata=profun(pm);
elseif problem==8
    % Sphere with approximate BC
    profun=@problem8;
    probdata=profun(pm);
elseif problem==9
    %Solenoid coil and conducting plate
    profun=@problemSolenoid;
    probdata=profun(pm);
elseif problem==10
    profun=@problemCube;
    probdata=profun(pm,option);
elseif problem==11
    profun=@problemSphere;
    probdata=profun(pm,option);
elseif problem==12
    profun=@problemBeam;
    probdata=profun(pm,option);
elseif problem==13
    profun=@problemCylinder;
    probdata=profun(pm,option);
elseif problem==14
    profun=@ProblemToy;
    probdata=profun(pm,option);
elseif problem==15
    profun=@ProblemToy1Shield;
    probdata=profun(pm,option);
elseif problem==16
    profun=@ProblemToyCube;
    probdata=profun(pm,option);
elseif problem==17
    profun=@ProblemToyCubeCoils;
    probdata=profun(pm,option);
end