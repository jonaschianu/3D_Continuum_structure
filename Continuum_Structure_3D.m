%%% a 2D continuum structure %%%
clear;
close all;
clc;
% Plot parameters
nms = 3; % node marker size for displacements plot
% Material properties
E = 1; % Young modulus in [N/m^2]
v = 0.3;
% Size of the structure
Lx = 1; % length of the structure in x-direction [m]
Ly = 1;
Lz=1;
% Mesh
nx = 10; % n. elements in x
ny = 10;
nz = 10;
m = nx*ny*nz; % total n. of elements
Dx = Lx/nx; % size of the element along x
Dy = Ly/ny;
Dz = Lz/nz;
% Nodes
N = (nx+1)*(ny+1)*(nz+1); % total n. of nodes
% Node grid - we use here a grid numeration
ix = [1:nx+1];
iy = [1:ny+1];
iz = [1:nz+1];
% Node position
x = (ix-1)*Dx;
y = (iy-1)*Dy;
z = (iz-1)*Dz;
[xn,yn,zn] = meshgrid(x,y,z);
% Stiffness matrix
K = zeros(3*N,3*N);
Em = [(1-v)   v    v   0   0         0;
         v (1-v)   v  0 0     0;
         v    v   (1-v)   0   0     0;
         0    0   0   (1-2*v)/2   0   0;
         0    0   0   0   (1-2*v)/2   0 ;
         0    0   0   0   0   (1-2*v)/2;] * E/((1+v)*(1-2*v)); % plane strain
T = [1,0,0,0,0,0,0,0,0;
     0,0,0,0,1,0,0,0,0;
     0,0,0,0,0,0,0,0,1;
     0,1,0,1,0,0,0,0,0;
     0,0,1,0,0,0,1,0,0;
     0,0,0,0,0,1,0,1,0];
% Gauss points
g(1,:) = 1/sqrt(3)*[-1,-1,-1]; % the matrix g contains the coordinates of the Gaussian points
g(2,:) = 1/sqrt(3)*[-1,-1,1];
g(3,:) = 1/sqrt(3)*[-1,1,-1];
g(4,:) = 1/sqrt(3)*[-1,1,1];
g(5,:) = 1/sqrt(3)*[1,-1,-1];
g(6,:) = 1/sqrt(3)*[1,-1,1];
g(7,:) = 1/sqrt(3)*[1,1,-1];
g(8,:) = 1/sqrt(3)*[1,1,1];
for ez = 1:nz
    for ey = 1:ny
        for ex = 1:nx
            i = (ez-1)*(ny+1)*(nx+1)+(ey-1)*(nx+1)+ex;
            j = i+nx+1;
            k = i+(nx+1)*(ny+1);
            l= k+nx+1;
            x1 = xn(ey,ex,ez); x2 = xn(ey,ex+1,ez); x3 = x1; x4 = x2; x5 = x1; x6=x2; x7=x1; x8=x2; 
            y1 = yn(ey,ex,ez); y2 = y1; y3 = yn(ey+1,ex,ez); y4 = y3; y5=y1;y6=y1;y7=y3;y8=y3;
            z1 = zn(ey,ex,ez); z2 = z1; z3=z1; z4=z1; z5 = zn(ey,ex,ez+1); z6=z5;z7=z5;z8=z5;
            Kq = 0;
            for gg = 1:8 % the eight Gauss points
                xx = g(gg,1);
                yy = g(gg,2);
                zz = g(gg,3);
                Nxx = [-(1/8)*(yy-1)*(zz-1), (1/8)*(yy-1)*(zz-1), (1/8)*(yy+1)*(zz-1),-(1/8)*(yy+1)*(zz-1),(1/8)*(yy-1)*(zz+1),-(1/8)*(yy-1)*(zz+1),-(1/8)*(yy+1)*(zz+1),(1/8)*(yy+1)*(zz+1)];
                Nyy = [-(1/8)*(xx-1)*(zz-1), (1/8)*(xx+1)*(zz-1), (1/8)*(xx-1)*(zz-1),-(1/8)*(xx+1)*(zz-1),(1/8)*(xx-1)*(zz+1),-(1/8)*(xx+1)*(zz+1),-(1/8)*(xx-1)*(zz+1),(1/8)*(xx+1)*(zz+1)];
                Nzz = [-(1/8)*(xx-1)*(yy-1), (1/8)*(xx+1)*(yy-1), (1/8)*(xx-1)*(yy+1),-(1/8)*(xx+1)*(yy+1),(1/8)*(xx-1)*(yy-1),-(1/8)*(xx+1)*(yy-1),-(1/8)*(xx-1)*(yy+1),(1/8)*(xx+1)*(yy+1)];
                J = [Nxx*[x1;x2;x3;x4;x5;x6;x7;x8], Nxx*[y1;y2;y3;y4;y5;y6;y7;y8], Nxx*[z1;z2;z3;z4;z5;z6;z7;z8]; 
                    Nyy*[x1;x2;x3;x4;x5;x6;x7;x8], Nyy*[y1;y2;y3;y4;y5;y6;y7;y8], Nyy*[z1;z2;z3;z4;z5;z6;z7;z8];
                    Nzz*[x1;x2;x3;x4;x5;x6;x7;x8], Nzz*[y1;y2;y3;y4;y5;y6;y7;y8], Nzz*[z1;z2;z3;z4;z5;z6;z7;z8]];
                dJ = det(J);
                iJ = inv(J);
                A = T*[iJ,zeros(3,3),zeros(3,3);
                    zeros(3,3),iJ,zeros(3,3);
                    zeros(3,3),zeros(3,3),iJ];
                G = [Nxx(1)      0  0 Nxx(2)  0     0 Nxx(3)  0    0 Nxx(4)  0    0  Nxx(5)      0  0 Nxx(6)  0     0 Nxx(7)  0    0 Nxx(8)  0    0;
                     Nyy(1)      0  0 Nyy(2)   0   0 Nyy(3)   0   0 Nyy(4)   0   0  Nyy(5)      0  0 Nyy(6)   0   0 Nyy(7)   0   0 Nyy(8)   0   0;
                     Nzz(1)      0  0 Nzz(2)   0   0 Nzz(3)   0   0 Nzz(4)   0   0  Nzz(5)      0  0 Nzz(6)   0   0 Nzz(7)   0   0 Nzz(8)   0   0;
                          0 Nxx(1)   0   0  Nxx(2)   0   0 Nxx(3)  0    0 Nxx(4) 0 0 Nxx(5)   0   0  Nxx(6)   0   0 Nxx(7)  0    0 Nxx(8) 0;
                          0 Nyy(1)   0   0  Nyy(2)   0   0 Nyy(3)  0    0 Nyy(4) 0 0 Nyy(5)   0   0  Nyy(6)   0   0 Nyy(7)  0    0 Nyy(8) 0;
                          0 Nzz(1)   0   0  Nzz(2)   0   0 Nzz(3)  0    0 Nzz(4) 0 0 Nzz(5)   0   0  Nzz(6)   0   0 Nzz(7)  0    0 Nzz(8) 0;
                          0  0 Nxx(1)   0   0  Nxx(2)   0   0 Nxx(3)  0    0 Nxx(4) 0  0 Nxx(5)   0   0  Nxx(6)   0   0 Nxx(7)  0    0 Nxx(8);
                          0  0 Nyy(1)   0   0  Nyy(2)   0   0 Nyy(3)  0    0 Nyy(4) 0  0 Nyy(5)   0   0  Nyy(6)   0   0 Nyy(7)  0    0 Nyy(8);
                          0  0 Nzz(1)   0   0  Nzz(2)   0   0 Nzz(3)  0    0 Nzz(4) 0  0 Nzz(5)   0   0  Nzz(6)   0   0 Nzz(7)  0    0 Nzz(8)];
                B = A*G;
                Kq = Kq + B'*Em*B*dJ;
            end
            
            edofs = [3*i-2,3*i-1,3*i,3*(i+1)-2,3*(i+1)-1,3*(i+1),3*j-2,3*j-1,3*j,3*(j+1)-2,3*(j+1)-1,3*(j+1),3*k-2,3*k-1,3*k,3*(k+1)-2,3*(k+1)-1,3*(k+1),3*l-2,3*l-1,3*l,3*(l+1)-2,3*(l+1)-1,3*(l+1)];

            K(edofs,edofs) = K(edofs,edofs)+Kq;
        end
    end
end

% Boundary conditions
ixfix = [1:nx+1]; % position of the fixed nodes in the node grid
iyfix = 1;
izfix = [1:nz+1];
% (izfix-1)*(ny+1)*(nx+1)+
ifix=[];
for i = 1:size(izfix,2)
    A= (iyfix-1)*(nx+1)+ixfix+(nx+1)*(ny+1)*(izfix(i)-1);
    ifix =[ifix,A];
end

fix_dofs = [3*(ifix)-2 3*(ifix)-1 3*(ifix)]; % both u and v are constrained here
% Loads
f = zeros(3*N,1);
q = -0.1; % this is the distributed load at the top surface in [N/m]
ixload = [1:nx+1]; % position of the loaded nodes in the node grid
iyload = (ny+1);
izload = [1:nz+1];
iload=[];
for i = 1:size(izload,2)
    A= (iyload-1)*(nx+1)+ixload+(nx+1)*(ny+1)*(izload(i)-1);
    iload =[iload,A];
end
iload;
load_dofs = 3*(iload)-1; % we are loading in the y direction
% xload = [xn(iyload,1,izload) xn(iyload,ixload,izload) xn(iyload,end,izload)];
W=[]; %weight
for k=1:(nz+1)
    for i=1:(nx+1)
        if k==1 || k==(nz+1)
            if i==1 || i==(nx+1)
                W=[W,Dx+Dy];
            else
                W=[W,2*Dx+Dy];
            end     
        else
            if i==1 || i==(nx+1)
                W=[W,Dx+2*Dy];
            else
                W=[W,2*Dx+2*Dy];
            end    
            
        end
    end
end
W;
% f(load_dofs) = q*0.5*(xload(3:end,3:end)-xload(1:end-2,1:end-2));
f(load_dofs) = q*0.05*W;
% f(20)=-0.05;
% f(23)=-0.05;
% f(26)=-0.05;
% 
% f(47)=-0.05;
% f(50)=-0.05;
% f(53)=-0.05;
% 
% f(74)=-0.05;
% f(77)=-0.05;
% f(80)=-0.05;

% Solution
free_dofs = setdiff([1:3*N],fix_dofs);
K_free=K(free_dofs,free_dofs);
u = zeros(3*N,1); % the function 'zeros' gives a matrix by defauls, if one of the sizes is 1, that gives a tensor
u(free_dofs) = K(free_dofs,free_dofs)\f(free_dofs);
% Plot results
figure(1)
hold on
grid on
grid minor
view(45,45)
title('Nodes','fontsize',15)

xn_fix=xn(iyfix,ixfix,izfix);
yn_fix=yn(iyfix,ixfix,izfix);
zn_fix=zn(iyfix,ixfix,izfix);

xn_load=xn(iyload,ixload,izload);
yn_load=yn(iyload,ixload,izload);
zn_load=zn(iyload,ixload,izload);

xn_fix_=reshape(xn_fix,1,numel(xn_fix));
yn_fix_=reshape(yn_fix,1,numel(yn_fix));
zn_fix_=reshape(zn_fix,1,numel(zn_fix));

xn_load_=reshape(xn_load,1,numel(xn_load));
yn_load_=reshape(yn_load,1,numel(yn_load));
zn_load_=reshape(zn_load,1,numel(zn_load));

xn_=reshape(xn,1,numel(xn));
yn_=reshape(yn,1,numel(yn));
zn_=reshape(zn,1,numel(zn));
% plot(xn_fix_,zn_fix_)
plot3(xn_fix_,zn_fix_,yn_fix_,'rs','Markersize',6)
plot3(xn_load_,zn_load_,yn_load_,'gv','MarkerEdgeColor',[0,0.7,0],'Markersize',6)
plot3(xn_,zn_,yn_,'ko','Markersize',3)
legend('fixed nodes','loaded nodes','all nodes','Location','NorthEastOutside')
axis equal;
%%%
figure(2)
title('Displacements','fontsize',15)
hold on; grid on; grid minor;
view(45,45)
grey = 0.6*[1 1 1]; % light grey as the color for the undeformed structure
for ez = 1:nz
    for ey = 1:ny
        for ex = 1:nx
%             x1 = xn(ey,ex,ez); x2 = xn(ey,ex+1,ez); x3 = x1; x4 = x2; x5 = x1; x6=x2; x7=x1; x8=x2; 
%             y1 = yn(ey,ex,ez); y2 = y1; y3 = yn(ey+1,ex,ez); y4 = y3; y5=y1;y6=y1;y7=y3;y8=y3;
%             z1 = zn(ey,ex,ez); z2 = z1; z3=z1; z4=z1; z5 = zn(ey,ex,ez+1); z6=z5;z7=z5;z8=z5;
            
            x1 = xn(ey,ex,ez); x2 = xn(ey,ex+1,ez); x3 = x1; x4 = x2; x5 = x1; x6=x2; x7=x1; x8=x2; 
            y1 = yn(ey,ex,ez); y2 = y1; y3 = yn(ey+1,ex,ez); y4 = y3; y5=y1;y6=y1;y7=y3;y8=y3;
            z1 = zn(ey,ex,ez); z2 = z1; z3=z1; z4=z1; z5 = zn(ey,ex,ez+1); z6=z5;z7=z5;z8=z5;
            
            i = (ez-1)*(ny+1)*(nx+1)+(ey-1)*(nx+1)+ex;
            j = i+nx+1;
            k = i+(nx+1)*(ny+1);
            l= k+nx+1;
%[3*i-2,3*i-1,3*i,3*(i+1)-2,3*(i+1)-1,3*(i+1),3*j-2,3*j-1,3*j,3*(j+1)-2,3*(j+1)-1,3*(j+1),3*k-2,3*k-1,3*k,3*(k+1)-2,3*(k+1)-1,3*(k+1),3*l-2,3*l-1,3*l,3*(l+1)-2,3*(l+1)-1,3*(l+1)];
            u1 = u(3*i-2,1); v1 = u(3*i-1,1); w1 = u(3*i,1);
            u2 = u(3*(i+1)-2,1); v2 = u(3*(i+1)-1,1); w2 = u(3*(i+1),1);
            u3 = u(3*j-2,1); v3 = u(3*j-1,1); w3 = u(3*j,1);
            u4 = u(3*(j+1)-2,1); v4 = u(3*(j+1)-1,1); w4 = u(3*(j+1),1);
            u5 = u(3*k-2,1); v5 = u(3*k-1,1); w5 = u(3*k,1);
            u6 = u(3*(k+1)-2,1); v6 = u(3*(k+1)-1,1); w6 = u(3*(k+1),1);
            u7 = u(3*l-2,1); v7 = u(3*l-1,1); w7 = u(3*l,1);
            u8 = u(3*(l+1)-2,1); v8 = u(3*(l+1)-1,1); w8 = u(3*(l+1),1);
            plot3([x1,x2],[z1,z2],[y1,y2],'k-o','Color',grey,'MarkerSize',nms,'markerfacecolor',grey)
            plot3([x1,x3],[z1,z3],[y1,y3],'k-o','Color',grey,'MarkerSize',nms,'markerfacecolor',grey)
            plot3([x2,x4],[z2,z4],[y2,y4],'k-o','Color',grey,'MarkerSize',nms,'markerfacecolor',grey)
            plot3([x3,x4],[z3,z4],[y3,y4],'k-o','Color',grey,'MarkerSize',nms,'markerfacecolor',grey)
            plot3([x5,x6],[z5,z6],[y5,y6],'k-o','Color',grey,'MarkerSize',nms,'markerfacecolor',grey)
            plot3([x5,x7],[z5,z7],[y5,y7],'k-o','Color',grey,'MarkerSize',nms,'markerfacecolor',grey)
            plot3([x6,x8],[z6,z8],[y6,y8],'k-o','Color',grey,'MarkerSize',nms,'markerfacecolor',grey)
            plot3([x7,x8],[z7,z8],[y7,y8],'k-o','Color',grey,'MarkerSize',nms,'markerfacecolor',grey)
            plot3([x2,x6],[z2,z6],[y2,y6],'k-o','Color',grey,'MarkerSize',nms,'markerfacecolor',grey)
            plot3([x4,x8],[z4,z8],[y4,y8],'k-o','Color',grey,'MarkerSize',nms,'markerfacecolor',grey)
            plot3([x1,x5],[z1,z5],[y1,y5],'k-o','Color',grey,'MarkerSize',nms,'markerfacecolor',grey)
            plot3([x3,x7],[z3,z7],[y3,y7],'k-o','Color',grey,'MarkerSize',nms,'markerfacecolor',grey)
            
            plot3([x1+u1,x2+u2],[z1+w1,z2+w2],[y1+v1,y2+v2],'b-o','MarkerSize',nms,'markerfacecolor','b')
            plot3([x1+u1,x3+u3],[z1+w1,z3+w3],[y1+v1,y3+v3],'b-o','MarkerSize',nms,'markerfacecolor','b')
            plot3([x2+u2,x4+u4],[z2+w2,z4+w4],[y2+v2,y4+v4],'b-o','MarkerSize',nms,'markerfacecolor','b')
            plot3([x3+u3,x4+u4],[z3+w3,z4+w4],[y3+v3,y4+v4],'b-o','MarkerSize',nms,'markerfacecolor','b')
            plot3([x5+u5,x6+u6],[z5+w5,z6+w6],[y5+v5,y6+v6],'b-o','MarkerSize',nms,'markerfacecolor','b')
            plot3([x5+u5,x7+u7],[z5+w5,z7+w7],[y5+v5,y7+v7],'b-o','MarkerSize',nms,'markerfacecolor','b')
            plot3([x6+u6,x8+u8],[z6+w6,z8+w8],[y6+v6,y8+v8],'b-o','MarkerSize',nms,'markerfacecolor','b')
            plot3([x7+u7,x8+u8],[z7+w7,z8+w8],[y7+v7,y8+v8],'b-o','MarkerSize',nms,'markerfacecolor','b')
            plot3([x2+u2,x6+u6],[z2+w2,z6+w6],[y2+v2,y6+v6],'b-o','MarkerSize',nms,'markerfacecolor','b')
            plot3([x4+u4,x8+u8],[z4+w4,z8+w8],[y4+v4,y8+v8],'b-o','MarkerSize',nms,'markerfacecolor','b')
            plot3([x1+u1,x5+u5],[z1+w1,z5+w5],[y1+v1,y5+v5],'b-o','MarkerSize',nms,'markerfacecolor','b')
            plot3([x3+u3,x7+u7],[z3+w3,z7+w7],[y3+v3,y7+v7],'b-o','MarkerSize',nms,'markerfacecolor','b')

            % average strains and stresses in the element
            ee11 = 0.5*((u2-u1)/(x2-x1)+(u4-u3)/(x4-x3)+(u6-u5)/(x6-x5)+(u8-u7)/(x8-x7));
            ee22 = 0.5*((v3-v1)/(y3-y1)+(v4-v2)/(y4-y2)+(v7-v5)/(y7-y5)+(v8-v6)/(y8-y6));
            ee33 = 0.5*((w5-w1)/(z5-z1)+(w6-w2)/(z6-z2)+(w7-w3)/(z7-z3)+(w8-w4)/(z8-z4));
            ee12 = 0.5*((u3-u1)/(y3-y1)+(u4-u2)/(y4-y2)+(u7-u5)/(y7-y5)+(u8-u6)/(y8-y6)+(v2-v1)/(x2-x1)+(v4-v3)/(x4-x3)+(v6-v5)/(x6-x5)+(v8-v7)/(x8-x7));
            ee13 = 0.5*((u5-u1)/(z5-z1)+(u6-u2)/(z6-z2)+(u7-u3)/(z7-z3)+(u8-u4)/(z8-z4)+(w2-w1)/(x2-x1)+(w4-w3)/(x4-x3)+(w6-w5)/(x6-x5)+(w8-w7)/(x8-x7));
            ee23 = 0.5*((v5-v1)/(z5-z1)+(v6-v2)/(z6-z2)+(v7-v3)/(z7-z3)+(v8-v4)/(z8-z4)+(w3-w1)/(y3-y1)+(w4-w2)/(y4-y2)+(w7-w5)/(y7-y5)+(w8-w6)/(y8-y6));
            eev(ey,ex,ez,:) = [ee11;ee22;ee33;ee12;ee13;ee23];
            sv(ey,ex,ez,:) = Em*[ee11;ee22;ee33;ee12;ee13;ee23];
            xc(ey,ex,ez) = (1/8)*(x1+x2+x3+x4+x5+x6+x7+x8);
            yc(ey,ex,ez) = (1/8)*(y1+y2+y3+y4+y5+y6+y7+y8);
            zc(ey,ex,ez) = (1/8)*(z1+z2+z3+z4+z5+z6+z7+z8);
        end
    end
end

axis equal
% for splt = 1:3
%     figure(2+splt)
% %     contour3(xc(:,:,1),yc(:,:,1),zc(:,:,1))
%     isosurface(xc, yc, zc, sv(:,:,:,splt))
%     colorbar
%     axis equal
% end