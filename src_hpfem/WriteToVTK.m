function WriteToVTK(coordn,intman,phen,pen,filename)



% Open the file.
fid = fopen(filename, 'w');

if fid == -1
    error('Cannot open file for writing.');
end
M1 = size(coordn,1);
[M2,N2] = size(intman);

% coordn = [coordn,zeros(M1,1)];
intman = [N2*ones(M2,1),intman - ones(M2,N2)];
% velo = [velo, zeros(M1,1)];

% New line.
nl = sprintf('\n'); % Stupid matlab doesn't interpret \n normally.

% Write the file header 
fwrite(fid, ['# vtk DataFile Version 1.0' nl '3D Unstructured grid' nl 'ASCII' nl nl]);

% Write the points
fwrite(fid, ['DATASET UNSTRUCTURED_GRID' nl 'POINTS ' num2str(M1) ' float' nl]);
for i = 1:M1
    pcoord = coordn(i,:);
    fprintf(fid,'%3.2f %3.2f %3.2f\n',pcoord);
end

% Newline.
fwrite(fid, nl);
    
% Write the cells (conectivity)
fwrite(fid, ['CELLS ' num2str(M2) ' ' num2str(M2*(N2+1)) nl]);
for i = 1:M2
    pintma = [intman(i,:)];
    fprintf(fid,'%d %d %d %d %d\n',pintma);
end

% Newline.
fwrite(fid, nl);

% Write the cells type (5 - triangle 10-tetrahedra)
fwrite(fid, ['CELL_TYPES ' num2str(M2) nl]);
ctype = 10;
for i = 1:M2
    fwrite(fid,num2str(ctype));
    fwrite(fid, nl);
end

% Newline.
%fwrite(fid, nl);

% Write the Cell data (ie which subdomain?)
% fwrite(fid, ['CELL_DATA ' num2str(M2) nl]);
% fwrite(fid, ['SCALARS cell_scalars int 1' nl ...
%             'LOOKUP_TABLE default' nl]);
% ctype=0;
% for i = 1:M2
%     fwrite(fid,num2str(ctype));
%     fwrite(fid, nl);
% end


% Write the pressure
fwrite(fid, ['POINT_DATA ' num2str(M1) nl 'SCALARS pressure float' nl ...
            'LOOKUP_TABLE default' nl]);
for i = 1:M1
%     ppres = press(i,:);
    ppres = 1;
    fprintf(fid,'%5.4e\n',ppres);
     fwrite(fid,num2str(ppres));
     fwrite(fid, nl);
end

% Newline.
fwrite(fid, nl);

% Write Re(H)
fwrite(fid, ['VECTORS remagnetic float' nl]);
for i = 1:M1
    pheno= real(phen(i,:));
    fprintf(fid,'%5.4e %5.4e %5.4e\n',pheno);
end

% Newline.
fwrite(fid, nl);

% Write Im(H)
fwrite(fid, ['VECTORS immagnetic float' nl]);
for i = 1:M1
    pheno= imag(phen(i,:));
    fprintf(fid,'%5.4e %5.4e %5.4e\n',pheno);
end
% Newline.
fwrite(fid, nl);

% Write Re(J^o)  (J^o = Ohmic, eddy currents)
fwrite(fid, ['VECTORS reeddy float' nl]);
for i = 1:M1
    peno= real(pen(i,:));
    fprintf(fid,'%5.4e %5.4e %5.4e\n',peno);
end

% Newline.
fwrite(fid, nl);

% Write Re(J^o)  (J^o = Ohmic, eddy currents)
fwrite(fid, ['VECTORS imeddy float' nl]);
for i = 1:M1
    peno= imag(pen(i,:));
    fprintf(fid,'%5.4e %5.4e %5.4e\n',peno);
end



fclose(fid);