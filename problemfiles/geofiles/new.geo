
## An elliptic cylinder

algebraic3d


solid fincylz = cylinder ( 0, 0, 10; 0, 0, -10; 4 )
	and plane (0, 0, -5; 0, 0, -1)
	and plane (0, 0, 5; 0, 0, 1);	


# cut cylinder by planes:

solid fincylx = cylinder ( 10, 0, 0; -10, 0, 0; 5 )
	and plane (-5, 0, 0; -1, 0, 0)
	and plane (5, 0, 0; 1, 0, 0);

solid cutcone2=fincylz and  fincylx;

tlo cutcone2 -col=[1,0,0];
#tlo fincyl -col=[0,1,0];
