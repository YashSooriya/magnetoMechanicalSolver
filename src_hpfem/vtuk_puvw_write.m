function vtuk_puvw_write ( title, node_num, element_num, ...
   xyz, element_node,pen, phen,displ)

% Open the file.
output_unit = fopen(title, 'w');
if output_unit == -1
    error('Cannot open file for writing.');
end

element_order=4;
%*****************************************************************************80
%
%% VTU_PUVW_WRITE curl theta_1,curl theta_2, curl theta_3 and sigma theta_1, sigma theta_2, sigma theta_3 to a VTU file.
%
%  Discussion:
%
%    The data is assumed to have been computed by a finite element program
%    for a 3D geometry which has been meshed using tetrahedral elements 
%    of 4 or 10 nodes.
%
%    The solution data includes the pressure and velocity vector at each node.
%
%    The VTU format used here is the modern XML-based format used by the
%    Visual Toolkit for unstructured grid data.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    20 December 2010
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer OUTPUT_UNIT, the output unit.
%
%    Input, string TITLE, a title for the data.
%
%    Input, integer NODE_NUM, the number of nodes.
%
%    Input, integer ELEMENT_NUM, the number of elements.
%
%    Input, integer ELEMENT_ORDER, the order of the elements.
%
%    Input, real XYZ(3,NODE_NUM), the node coordinates.
%
%    Input, integer ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM), the
%    nodes that make up each element.  The node indices are zero based.
%
%    Input, real P(1,NODE_NUM), the pressure at each node.
%
%    Input, real UVW(3,NODE_NUM), the velocity at each node.
%
%   if ( element_order == 10 )
%     fprintf ( 1, '\n' );
%     fprintf ( 1, 'VTU_PUVW_WRITE - Note:\n' );
%     fprintf ( 1, '  As a temporary measure, we are handling quadratic tets\n' );
%     fprintf ( 1, '  as though they were linear.  That is, the element information\n' );
%     fprintf ( 1, '  will only use the vertices, not the midside nodes.\n' );
%     element_order = 4;
%   end

  fprintf ( output_unit, '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="BigEndian">\n' );
  fprintf ( output_unit, '  <UnstructuredGrid>\n' );
  fprintf ( output_unit, '    <Piece NumberOfPoints="%d" NumberOfCells="%d">\n', ...
    node_num, element_num );
 fprintf ( output_unit, '      <PointData Scalars="scalars">\n' );
  fprintf ( output_unit, '        <DataArray type="Float32" NumberOfComponents="3" Name="Re(J)Eddy P1" Format="ascii">\n' );
  for node = 1 : node_num
    fprintf ( output_unit, '        %20.16f  %20.16f  %20.16f\n', real(pen(node,1:3)) );
  end
  fprintf ( output_unit, '        </DataArray>\n' );
  fprintf ( output_unit, '        <DataArray type="Float32" NumberOfComponents="3" Name="Im(J)Eddy P1" Format="ascii">\n' );
  for node = 1 : node_num
    fprintf ( output_unit, '        %20.16f  %20.16f  %20.16f\n', imag(pen(node,1:3)) );
  end  
   fprintf ( output_unit, '        </DataArray>\n' );
  fprintf ( output_unit, '        <DataArray type="Float32" NumberOfComponents="1" Name="J Eddy P1" Format="ascii">\n' );
  for node = 1 : node_num
    fprintf ( output_unit, '        %20.16f  %20.16f  %20.16f\n',  pen(node,1:3) );
  end  
%      fprintf ( output_unit, '        </DataArray>\n' );
%   fprintf ( output_unit, '        <DataArray type="Float32" NumberOfComponents="1" Name="Norm(J_y) Eddy P1" Format="ascii">\n' );
%   for node = 1 : node_num
%     fprintf ( output_unit, '        %20.16f  %20.16f  %20.16f\n',  norm(pen(node,2))*(real(pen(node,2))/norm(real(pen(node,2))))) ;
%   end  
   fprintf ( output_unit, '        </DataArray>\n' );
  fprintf ( output_unit, '        <DataArray type="Float32" NumberOfComponents="3" Name="Re(H)Magnetic P1" Format="ascii">\n' );
  for node = 1 : node_num
    fprintf ( output_unit, '        %20.16f  %20.16f  %20.16f\n', real(phen(node,1:3)) );
  end
  fprintf ( output_unit, '        </DataArray>\n' );
  fprintf ( output_unit, '        <DataArray type="Float32" NumberOfComponents="3" Name="Im(H)Magnetic P1" Format="ascii">\n' );
  for node = 1 : node_num
    fprintf ( output_unit, '        %20.16f  %20.16f  %20.16f\n', imag(phen(node,1:3)) );
  end    
     fprintf ( output_unit, '        </DataArray>\n' );
  fprintf ( output_unit, '        <DataArray type="Float32" NumberOfComponents="3" Name="Re(B)Magnetic P1" Format="ascii">\n' );
  for node = 1 : node_num
    fprintf ( output_unit, '        %20.16f  %20.16f  %20.16f\n', 1.256637061435917e-06*real(phen(node,1:3)) );
  end
  fprintf ( output_unit, '        </DataArray>\n' );
  fprintf ( output_unit, '        <DataArray type="Float32" NumberOfComponents="3" Name="Im(B)Magnetic P1" Format="ascii">\n' );
  for node = 1 : node_num
    fprintf ( output_unit, '        %20.16f  %20.16f  %20.16f\n', 1.256637061435917e-06*imag(phen(node,1:3)) );
  end
    %=================================================================================================================================
    % Mechanical solution
    %=================================================================================================================================
    
  fprintf ( output_unit, '        </DataArray>\n' );
  fprintf ( output_unit, '        <DataArray type="Float32" NumberOfComponents="3" Name="Displacement" Format="ascii">\n' );
  for node = 1 : node_num
    fprintf ( output_unit, '        %20.16f  %20.16f  %20.16f\n', real(displ(node,1:3)));
  end



%      fprintf ( output_unit, '        </DataArray>\n' );
%   fprintf ( output_unit, '        <DataArray type="Float32" NumberOfComponents="1" Name="Norm(J_y) Eddy P1" Format="ascii">\n' );
%   for node = 1 : node_num
%     fprintf ( output_unit, '        %20.16f  %20.16f  %20.16f\n',  norm(pen(node,2))*(real(pen(node,2))/norm(real(pen(node,2))))) ;
%   end  

  
%   fprintf ( output_unit, '        </DataArray>\n' );
%   fprintf ( output_unit, '        <DataArray type="Float32" NumberOfComponents="1" Name="Norm(B_y) Magnetic P1" Format="ascii">\n' );
%   for node = 1 : node_num
%     fprintf ( output_unit, '        %20.16f  %20.16f  %20.16f\n', 1.256637061435917e-06*norm(phen(node,2))*(real(phen(node,2))/norm(real(phen(node,2))) ));
%   end    
 fprintf ( output_unit, '        </DataArray>\n' );
 fprintf ( output_unit, '      </PointData>\n' );
 fprintf ( output_unit, '      <CellData Scalars="scalars">\n' );
 fprintf ( output_unit, '        <DataArray type="Int32" Name="Suddomain" Format="ascii">\n' );
 for node = 1 : element_num
    fprintf ( output_unit, '       %d ', 1 );
    fprintf ( output_unit, '\n' );
  end
 fprintf ( output_unit, '        </DataArray>\n' );
 fprintf ( output_unit, '      </CellData>\n' );
 fprintf ( output_unit, '      <Points>\n' );
  fprintf ( output_unit, '        <DataArray type="Float32" NumberOfComponents="3" Format="ascii">\n' );
  for node = 1 : node_num
    fprintf ( output_unit, '          %20.16f  %20.16f  %20.16f\n', xyz(node,1:3) );
  end
  fprintf ( output_unit, '        </DataArray>\n' );
  fprintf ( output_unit, '      </Points>\n' );
  fprintf ( output_unit, '      <Cells>\n' );

  fprintf ( output_unit, '        <DataArray type="Int32" Name="connectivity" Format="ascii">\n' );
  for element = 1 : element_num
    for order = 1 : element_order
        index=element_node(element,order)-1;
      fprintf ( output_unit, '      %d', index );
    end
    fprintf ( output_unit, '\n' );
  end
  fprintf ( output_unit, '        </DataArray>\n' );
  fprintf ( output_unit, '        <DataArray type="Int32" Name="offsets" Format="ascii">\n' );
  offset = 0;
  for element = 1 : element_num
    offset = offset + element_order;
   %offset=element_node(element,element_order)-1;
    fprintf ( output_unit, '        %d ', offset );
  end
  fprintf ( output_unit, '\n' );
 fprintf ( output_unit, '        </DataArray>\n' );
  fprintf ( output_unit, '        <DataArray type="UInt8" Name="types" Format="ascii">\n' );
  if ( element_order == 4 )
    for element = 1 : element_num
      fprintf ( output_unit, '      10 ' );
    end
    fprintf ( output_unit, '\n' );
    
  elseif ( element_order == 10 )
    for element = 1 : element_num
      fprintf ( output_unit, '24\n' );
    end
  end
  fprintf ( output_unit, '        </DataArray>\n' );

  fprintf ( output_unit, '      </Cells>\n' );
  fprintf ( output_unit, '    </Piece>\n' );
  fprintf ( output_unit, '  </UnstructuredGrid>\n' );
  fprintf ( output_unit, '</VTKFile>\n' );
  fclose(output_unit);

  return
end