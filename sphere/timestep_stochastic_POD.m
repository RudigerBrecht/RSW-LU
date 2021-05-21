
function  [Vs,Hs,k,dd,ee]=timestep_stochastic_POD(U_svd,S_svd,c_xx_yy_zz,c_xy_yz,c_xz,Dij,Hinit,Vtk,Div_hex,hex2tri,tt,Vk,Vk_old,Hk,Hk_0,Bk,delt,dissipation,tol,Div_tri,Ac2i,Dnor,VAR_Curl_vxmis,VAR_Curl_vxpls,Djim,Dijm,Dijp,Djip,WeightIp,WeightIm,WeightJp,WeightJm,...
    edge_length_hex,VAR_SWE_Ekin,g,f_pls,f_mis)

global nu tang_x tang_y nor_x nor_y nor_z  Av2c Ai2v rec_u rec_v rec_w nedges edge_length_tri Dtan

Vs     = Vk;
Hs     = Hk;

alphaV=1;
alphaH=1;

u = Ac2i*(rec_u*Vk);
v = Ac2i*(rec_v*Vk);
w = Ac2i*(rec_w*Vk);
%% create noise for water depth and velocity field
[sdBdt] = createNoise(S_svd,U_svd);

c_xx_yy_zz=c_xx_yy_zz*delt;
c_xy_yz   =c_xy_yz*delt;
c_xz      =c_xz*delt;


[u_noise,u_diffu]=applynoise(u,sdBdt,c_xx_yy_zz,c_xy_yz,c_xz);
[v_noise,v_diffu]=applynoise(v,sdBdt,c_xx_yy_zz,c_xy_yz,c_xz);
[w_noise,w_diffu]=applynoise(w,sdBdt,c_xx_yy_zz,c_xy_yz,c_xz);

Vn_noise=u_noise.*nor_x(:)+v_noise.*nor_y(:)+w_noise.*nor_z(:);
Vn_diffu=u_diffu.*nor_x(:)+v_diffu.*nor_y(:)+w_diffu.*nor_z(:);

[h_noise,h_diffu]=applynoise(Ac2i*Hk,sdBdt,c_xx_yy_zz,c_xy_yz,c_xz);

% add noise to the height field 
Hs = Hs + delt*(Av2c*Ai2v)*(alphaH^2*h_diffu-alphaH*h_noise);        

%% %%% CALCULATION OF VALUES AT TIME k  %%%%%%%%%%%%

%--- DIVERGENCE TERM/CONSTANT PART
DvHk = - Div_tri*((Ac2i*Hk).*Vk);
Ck   =   Hk + 0.5*delt*1*DvHk;

%— GRADIENT TERM/CONSTANT PART
GrHk   = - g*(Dnor*(Hk + Bk));
GrHs   =   GrHk;

%— NONLINEAR ADVECTION & CORIOLIS TERM/CONSTANT PARTS

VAR_NONLIN = 1;

q_mis   = (VAR_NONLIN*VAR_Curl_vxmis*Vk + f_mis);
q_pls   = (VAR_NONLIN*VAR_Curl_vxpls*Vk + f_pls);


nonterm = (q_mis.*((Djim*Hk).*(WeightIm*Vk) + (Dijm*Hk).*(WeightJm*Vk)) ...
    - q_pls.*((Djip*Hk).*(WeightIp*Vk)     + (Dijp*Hk).*(WeightJp*Vk)))./edge_length_hex;

Advk     = nonterm./(Dij*Hk) - VAR_NONLIN*(     (VAR_SWE_Ekin*(Vk.*Vk) )     ) + dissipation;
Advs = Advk;

%— CONSTANT FLUX TERM

Fk = Vk + 0.5*delt*(1*Advk + 0*GrHk+alphaV^2*Vn_diffu-alphaV*Vn_noise);

%%%%%%   CALCULATION OF VALUES V*,H* FOR k+1 %%%%%%

%— FIRST TIME: NOT CONVERGED YET
k = 1;

currentNorm=1;

%— FIXED POINT LOOP; FIXED POINT ITERATION
while currentNorm>tol
    
    %— KEEP OLD Vs IN MEMORY
    oldVs=Vs; oldHs=Hs;
    
    %— UPDATE OF VELOCITY V* AT t*:
    Vs = Fk + 0.5*delt*(1*Advs + 2*GrHs+alphaV^2*Vn_diffu-alphaV*Vn_noise)-nu.*hyplap(Vs);

    
    %— CONTINUITY EQUATION UPDATE AT t*
    DvHs = - Div_tri*((Ac2i*Hs).*Vs);
    
    Hs     = Ck + 0.5*delt*1*DvHs;
    
    %— GRADIENT TERM UPDATE AT t*
    GrHs = - g*Dnor*(Hs + Bk);
    

    q_mis   = (VAR_NONLIN*VAR_Curl_vxmis*Vs + f_mis);
    q_pls   = (VAR_NONLIN*VAR_Curl_vxpls*Vs + f_pls);

    
    %— NONLINEAR ADVECTION/CORIOLIS UPDATE AT t*
    nonterm = (q_mis.*((Djim*Hs).*(WeightIm*Vs) + (Dijm*Hs).*(WeightJm*Vs)) ...
        - q_pls.*((Djip*Hs).*(WeightIp*Vs)   + (Dijp*Hs).*(WeightJp*Vs)))./edge_length_hex;
    Advs    = nonterm./(Dij*Hs) - VAR_NONLIN*(     (VAR_SWE_Ekin*( Vs.*Vs) )     ) + dissipation;
    
    
    
    k = k + 1; dd= norm(oldVs-Vs,inf); ee= norm(oldHs-Hs,inf);
    currentNorm=norm(oldVs-Vs)+norm(oldHs-Hs,inf);
    
end

end

function [hn,hd]=applynoise(h,sdBdt,c_xx_yy_zz,c_xy_yz,c_xz)

    global Ai2v Av2c Dnor Dtan nor_x nor_y nor_z tang_x tang_y tang_z nedges

    hex_u = Ai2v*h;
    tri_u = Av2c*hex_u;
    
    hx = (Dnor*tri_u(:)).*nor_x(:)+(Dtan*hex_u(:)).*tang_x(:);
    hy = (Dnor*tri_u(:)).*nor_y(:)+(Dtan*hex_u(:)).*tang_y(:);
    hz = (Dnor*tri_u(:)).*nor_z(:)+(Dtan*hex_u(:)).*tang_z(:);
    

    hn = sdBdt(:,1).*hx+sdBdt(:,2).*hy+sdBdt(:,3).*hz; 
    hd = 0.5*(gradxx(c_xx_yy_zz(1:nedges).*h)+gradyy(c_xx_yy_zz(nedges+1:2*nedges).*h)+gradzz(c_xx_yy_zz(2*nedges+1:end).*h))...
         +gradxy(c_xy_yz(1:nedges).*h)+gradyz(c_xy_yz(nedges+1:end).*h)+gradxz(c_xz.*h);
        
end


function ul = hyplap(u)


ul= lap1(lap1(u));

end

function ul = lap1(u)

global Dnor Div_tri Dtan Curl_hex

    ul = Dnor*Div_tri*u-Dtan*Curl_hex*u;

end

function uxx = gradxx(u)

    global Ai2v Av2c Dnor Dtan nor_x tang_x 

    hex_u = Ai2v*u;
    tri_u = Av2c*hex_u;
    
    ux = (Dnor*tri_u(:)).*nor_x(:)+(Dtan*hex_u(:)).*tang_x(:);
    
    hex_u = Ai2v*ux;
    tri_u = Av2c*hex_u;
    
    uxx = (Dnor*tri_u(:)).*nor_x(:)+(Dtan*hex_u(:)).*tang_x(:);
end

function uyy = gradyy(u)

    global Ai2v Av2c Dnor Dtan nor_y tang_y

    hex_u = Ai2v*u;
    tri_u = Av2c*hex_u;
    
    uy = (Dnor*tri_u(:)).*nor_y(:)+(Dtan*hex_u(:)).*tang_y(:);
    
    hex_u = Ai2v*uy;
    tri_u = Av2c*hex_u;
    
    uyy = (Dnor*tri_u(:)).*nor_y(:)+(Dtan*hex_u(:)).*tang_y(:);
end

function uzz = gradzz(u)

    global Ai2v Av2c Dnor Dtan nor_z tang_z


    hex_u = Ai2v*u;
    tri_u = Av2c*hex_u;
    
    uz = (Dnor*tri_u(:)).*nor_z(:)+(Dtan*hex_u(:)).*tang_z(:);
    
    hex_u = Ai2v*uz;
    tri_u = Av2c*hex_u;
    
    uzz = (Dnor*tri_u(:)).*nor_z(:)+(Dtan*hex_u(:)).*tang_z(:);
end

function uxy = gradxy(u)

    global Ai2v Av2c Dnor Dtan nor_x tang_x nor_y tang_y


    hex_u = Ai2v*u;
    tri_u = Av2c*hex_u;
    
    uz = (Dnor*tri_u(:)).*nor_x(:)+(Dtan*hex_u(:)).*tang_x(:);
    
    hex_u = Ai2v*uz;
    tri_u = Av2c*hex_u;
    
    uxy = (Dnor*tri_u(:)).*nor_y(:)+(Dtan*hex_u(:)).*tang_y(:);
end

function uxz = gradxz(u)

    global Ai2v Av2c Dnor Dtan nor_x tang_x nor_z tang_z


    hex_u = Ai2v*u;
    tri_u = Av2c*hex_u;
    
    uz = (Dnor*tri_u(:)).*nor_x(:)+(Dtan*hex_u(:)).*tang_x(:);
    
    hex_u = Ai2v*uz;
    tri_u = Av2c*hex_u;
    
    uxz = (Dnor*tri_u(:)).*nor_z(:)+(Dtan*hex_u(:)).*tang_z(:);
end

function uyz = gradyz(u)

    global Ai2v Av2c Dnor Dtan nor_y tang_y nor_z tang_z


    hex_u = Ai2v*u;
    tri_u = Av2c*hex_u;
    
    uz = (Dnor*tri_u(:)).*nor_y(:)+(Dtan*hex_u(:)).*tang_y(:);
    
    hex_u = Ai2v*uz;
    tri_u = Av2c*hex_u;
    
    uyz = (Dnor*tri_u(:)).*nor_z(:)+(Dtan*hex_u(:)).*tang_z(:);
end

function [sdBdt] = createNoise(S,U)

%%
% multiply eiginvalues and iid Brownian motions
global intersect_x intersect_y intersect_z nedges R

rng('shuffle');

gamma = bsxfun(@times, randn(numel(S),1), S); 

% Karhunen Loeve theorem

sdBdt = reshape(U*gamma, [nedges 3]);


for edge = 1:nedges
    x=intersect_x(edge);
    y=intersect_y(edge);
    z=intersect_z(edge);
    
    P=1/R^2*[R^2-x^2 , -y*x   , -x*z  ;...
               -x*y  ,R^2-y^2 ,-y*z   ;...
               -x*z  ,-z*y    , R^2-z^2 ];
    
    sdBdt(edge,:)=P*sdBdt(edge,:)';
    
    
end





end



