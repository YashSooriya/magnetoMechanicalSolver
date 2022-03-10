function [Mv,Cv,Cvpre,CvReg,CvpreReg,Kv,Rv,I,J,Idir,Jdir,nmst,nmstDir,Kdv,Cdv,CdvReg,Mdv]=matrixAssembly(K,Ccond,Ccondpre,CReg,CpreReg,M,Res,Xdir,Mv,Cv,Cvpre,CvReg,CvpreReg,Kv,Rv,I,J,Idir,Jdir,bhelp,probstatic,esizet,esizeH1,nmst,nmstDir,omega,Kdv,Cdv,CdvReg,Mdv,freqSweep,NmechBodies)
%=========================================================================
% This function stores the local elemental system matrices into the global
% system vectors Mv,Cv,Kv,Rv with the position numbering as I and J.
% Inputs:
% K         : Local elemental stiffness matrix
% C         : Local elemental damping matrix
% M         : Local elemental mass matrix
% R         : Local elemental residual vector
% Xa        : Local solution array containing the solution and its time
%             derivatives Xa = [X dXdt dXdt2]
% unk       : Local to global elemental unknown numbering
% dir       : Local to global elemental direction flags
% nDOF      : Number of degrees of freedom per element
% Mv        : Global mass matrix vector
% Cv        : Global damping matrix vector
% Kv        : Global stiffness matrix vector
% Rv        : Global residual vector
% I         : Global matrix positioning vector (rows)
% J         : Global matrix positioning vector (columns)
% nz        : Counter for number of non zeros encountered
% nonLinear : Flag to determine if problem is non linear
% probstatic: Flag to determine if problem is static or transient
% Outputs:
% Mv        : Global mass matrix vector
% Cv        : Global damping matrix vector
% Kv        : Global stiffness matrix vector
% Rv        : Global residual vector
% I         : Global matrix positioning vector (rows)
% J         : Global matrix positioning vector (columns)
% nz        : Counter for number of non zeros encountered
%=========================================================================


if probstatic==1
    omega2=0;
elseif probstatic==0
    omega2=omega;
end
%-------------------------------------------------------------------------
% Local to global matrix assembly
%-------------------------------------------------------------------------
% Assemble linear system in vector format
for j = 1:esizet+3*esizeH1  %econt
    row = bhelp(j);
    if row>0
        for k = 1:esizet+3*esizeH1  %econt
            col = bhelp(k);
            if col>0
                nmst=nmst+1;
                I(nmst) = row;
                J(nmst) = col;
                Kv(nmst) = K(j,k);
                for nBody=1:NmechBodies
                    Cv{nBody}(nmst) = Ccond{nBody}(j,k);
                    Cvpre{nBody}(nmst)=  Ccondpre{nBody}(j,k);
                end
                CvReg(nmst) = CReg(j,k);
                CvpreReg(nmst)= CpreReg(j,k);
                Mv(nmst) = M(j,k);
                
            elseif col<0
                if freqSweep==1 && probstatic==0
                    
                    dirDOF=abs(col);
                    
                    if dirDOF>0
                        
                        nmstDir=nmstDir+1;
                        Idir(nmstDir)=row;
                        Jdir(nmstDir)=dirDOF;
                         Mdv(nmstDir)=M(j,k);
%                         Cdv(nmstDir)=C(j,k);
%                         CdvReg(nmstDir)=CReg(j,k);
%                         Kdv(nmstDir)=Kdv(nmstDir)+K(j,k);
%                         Mdv(nmstDir)=M(j,k);
                         for nBody=1:NmechBodies
                             Cdv{nBody}(nmstDir) = Ccond{nBody}(j,k);
% %                             Cdvpre{nBody}(nmstDir)=  Ccondpre{nBody}(j,k);
                         end
                         CdvReg(nmstDir)=CReg(j,k);
% %                         CdvpreReg(nmst)= CpreReg(j,k);
                        Kdv(nmstDir)=K(j,k);

                        
                    end
                else
                    
                    % move Dirihlet columns to right hand side
                    if probstatic==0
                        
                        
                        
                        
                        C=Ccond{1};
                        for nBody=2:NmechBodies
                            C=C+Ccond{nBody};
                        end
                        Rv(row)=Rv(row)+(complex(0,1)*omega*C(j,k)+omega*CReg(j,k)+K(j,k)-omega2^2*M(j,k))*Xdir(k,1);
                    else
                        C=Ccond{1};
                        for nBody=2:NmechBodies
                            C=C+Ccond{nBody};
                        end
                        Rv(row)=Rv(row)+(C(j,k)+CReg(j,k)+K(j,k))*Xdir(k,1);
                    end
                end
                
            end
        end
        % Add right hand side source terms
        Rv(row)=Rv(row)+Res(j);
    end
end