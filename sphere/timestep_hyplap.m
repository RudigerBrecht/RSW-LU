
function  [Vs,Hs,k,dd,ee]=timestep_hyplap(Dij,Hinit,Vtk,Div_hex,hex2tri,tt,Vk,Vk_old,Hk,Hk_0,Bk,delt,dissipation,tol,Div_tri,Ac2i,Dnor,VAR_Curl_vxmis,VAR_Curl_vxpls,Djim,Dijm,Dijp,Djip,WeightIp,WeightIm,WeightJp,WeightJm,...
    edge_length_hex,VAR_SWE_Ekin,g,f_pls,f_mis)

global rec_u rec_v rec_w nor_x nor_y tang_x tang_y Dtan Ai2v Av2c edge_length_tri nu

Vs     = Vk;
Hs     = Hk;

%%%%% CALCULATION OF VALUES AT TIME k  %%%%%%%%%%%%

%--- DIVERGENCE TERM/CONSTANT PART

 DvHk = - Div_tri*((Ac2i*Hk).*Vk);
Ck     =   Hk + 0.5*delt*1*DvHk;

%--- GRADIENT TERM/CONSTANT PART
GrHk   = - g*(Dnor*(Hk + Bk));
GrHs   =   GrHk;

%--- NONLINEAR ADVECTION & CORIOLIS TERM/CONSTANT PARTS

VAR_NONLIN = 1; 

q_mis   = (VAR_NONLIN*VAR_Curl_vxmis*Vk + f_mis);
q_pls   = (VAR_NONLIN*VAR_Curl_vxpls*Vk + f_pls);


    nonterm = (q_mis.*((Djim*Hk).*(WeightIm*Vk) + (Dijm*Hk).*(WeightJm*Vk)) ...
        - q_pls.*((Djip*Hk).*(WeightIp*Vk)     + (Dijp*Hk).*(WeightJp*Vk)))./edge_length_hex;
    Advk     = nonterm./(Dij*Hk) - VAR_NONLIN*(     (VAR_SWE_Ekin*(Vk.*Vk) )     ) + dissipation;
Advs = Advk;

%--- CONSTANT FLUX TERM

Fk = Vk + 0.5*delt*(1*Advk + 0*GrHk);

%%%%%%   CALCULATION OF VALUES V*,H* FOR k+1 %%%%%%

%--- FIRST TIME: NOT CONVERGED YET
 k = 1;

currentNorm=1;

%--- FIXED POINT LOOP; FIXED POINT ITERATION
while currentNorm>tol
    
    %--- KEEP OLD Vs IN MEMORY
    oldVs=Vs; oldHs=Hs; 
    
    %--- UPDATE OF VELOCITY V* AT t*:
    Vs = Fk + 0.5*delt*(1*Advs + 2*GrHs-nu.*hyplap(Vs));
    
    
    %--- CONTINUITY EQUATION UPDATE AT t*    
    DvHs = - Div_tri*((Ac2i*Hs).*Vs);

    Hs     = Ck + 0.5*delt*1*DvHs;
    
    %--- GRADIENT TERM UPDATE AT t*
    GrHs = - g*Dnor*(Hs + Bk);
    
        q_mis   = (VAR_NONLIN*VAR_Curl_vxmis*Vs + f_mis);
        q_pls   = (VAR_NONLIN*VAR_Curl_vxpls*Vs + f_pls);

        
    %--- NONLINEAR ADVECTION/CORIOLIS UPDATE AT t*

        nonterm = (q_mis.*((Djim*Hs).*(WeightIm*Vs) + (Dijm*Hs).*(WeightJm*Vs)) ...
            - q_pls.*((Djip*Hs).*(WeightIp*Vs)   + (Dijp*Hs).*(WeightJp*Vs)))./edge_length_hex;
        Advs    = nonterm./(Dij*Hs) - VAR_NONLIN*(     (VAR_SWE_Ekin*( Vs.*Vs) )     ) + dissipation;
  
    
    k = k + 1; dd= norm(oldVs-Vs,inf); ee= norm(oldHs-Hs,inf);
    currentNorm=norm(oldVs-Vs)+norm(oldHs-Hs,inf);
    
end

end


function ul = hyplap(u)


ul= lap1(lap1(u));

end

function ul = lap1(u)

global Dnor Div_tri Dtan Curl_hex

    ul = Dnor*Div_tri*u-Dtan*Curl_hex*u;

end

