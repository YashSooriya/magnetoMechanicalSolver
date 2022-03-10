function WriteToVTK(coordn,intman,phen,npoin,nelem,nip,filename)



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
%     fwrite(fid,num2str(pcoord),'double');
%     fwrite(fid, nl);
end

% Newline.
fwrite(fid, nl);
    
% Write the cells (conectivity)
fwrite(fid, ['CELLS ' num2str(M2+nip*nelem) ' ' num2str(M2*(N2+1)+2*nip*nelem) nl]);
for i = 1:M2
    pintma = [intman(i,:)];
    fprintf(fid,'%d %d %d %d %d\n',pintma);
%     fwrite(fid,num2str(pintma));
%     fwrite(fid, nl);
end
for i = 1:M2
    for j=1:nip
    pintma = [1,npoin+(i-1)*nip+j-1];
    fprintf(fid,'%d %d\n',pintma);
%     fwrite(fid,num2str(pintma));
%     fwrite(fid, nl);
    end
end


% Newline.
fwrite(fid, nl);

% Write the cells type (5 - triangle 10-tetrahedra)
fwrite(fid, ['CELL_TYPES ' num2str(M2+M2*nip) nl]);
ctype = 10;
for i = 1:M2
    fwrite(fid,num2str(ctype));
    fwrite(fid, nl);
end
ctype = 1;
for i = 1:M2
    for j=1:nip   
    fwrite(fid,num2str(ctype));
    fwrite(fid, nl);
    end
end


% Newline.
fwrite(fid, nl);

% Write the pressure
fwrite(fid, ['POINT_DATA ' num2str(M1) nl 'SCALARS pressure float' nl ...
            'LOOKUP_TABLE default' nl]);
for i = 1:M1
%     ppres = press(i,:);
    ppres = norm(phen(i,:));
    fprintf(fid,'%5.4e\n',ppres);
%     fwrite(fid,num2str(ppres));
%     fwrite(fid, nl);
end

% Newline.
fwrite(fid, nl);

% Write the vector velocity
% fwrite(fid, ['POINT_DATA ' num2str(M1) nl]);
fwrite(fid, ['VECTORS remagnetic float' nl]);
for i = 1:M1
    pheno= real(phen(i,:));
    fprintf(fid,'%5.4e %5.4e %5.4e\n',pheno);
%     fwrite(fid,num2str(pvelo));
%     fwrite(fid, nl);
end

% Newline.
fwrite(fid, nl);
fwrite(fid, ['VECTORS immagnetic float' nl]);
for i = 1:M1
    pheno= imag(phen(i,:));
    fprintf(fid,'%5.4e %5.4e %5.4e\n',pheno);
%     fwrite(fid,num2str(pvelo));
%     fwrite(fid, nl);
end

fclose(fid);