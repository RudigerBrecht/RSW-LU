function [Vs,Hs,k,dd,ee,sdBdt] = ops_st_diffu_homN(Vk,Vk_old,Hk,Bk,delt,Mass_kin,Mass_flat,dissipation,tol,N)

global Div Ac2i Av2i Dnor Curl edge_d_length  Ai2v Av2c tri_A
global g f
global VAR_Curl_vxmis VAR_Curl_vxpls 
global Var_Ac2i_iimis Var_Ac2i_iipls Var_Ac2i_jjmis Var_Ac2i_jjpls
global VAR_Rec_iimis  VAR_Rec_iipls  VAR_Rec_jjmis  VAR_Rec_jjpls
global VAR_Rec_vxmis  VAR_Rec_vxpls
global Rec_t2n_ham
global kinE_loc VAR_SWE_Dnor VAR_SWE_Ekin
global Idneig Rec_t2n_thub tangential_tx normal_tx tangential_ty normal_ty intersect_x intersect_y nablasq_tri nedges 
global Dtan kinE_loc_irreg Avrg_c2v


alphaV=1;
alphaH=1;

Vs     = Vk;      
Hs     = Hk;                         

[sdBdt]=createNoise(N);

 vn    = Rec_t2n_thub*Vk; 
 u     = normal_tx.*vn + tangential_tx.*Vk;
 v     = normal_ty.*vn + tangential_ty.*Vk;

vn_noise = ((sdBdt(:,1)).*gradx(u)+(sdBdt(:,2)).*grady(u)).*normal_tx + ...
           ((sdBdt(:,1)).*gradx(v)+(sdBdt(:,2)).*grady(v)).*normal_ty;
       
hn_noise = (sdBdt(:,1)).*gradx(Ac2i*Hk)+(sdBdt(:,2)).*grady(Ac2i*Hk); 


diffusionV = diffuVnew(u,v);
diffusionH = diffuHnew(Ac2i*Hs);


%%%%% CALCULATION OF VALUES AT TIME k  %%%%%%%%%%%%

%--- DIVERGENCE TERM/CONSTANT PART
DvHk   = - Div*((Ac2i*Hk).*Vk); 
Ck     =   Hk + 0.5*delt*1*DvHk;

%--- GRADIENT TERM/CONSTANT PART
GrHk   = - g*(Dnor*(Hk + Bk));                        
GrHs   =   GrHk;                          

%--- NONLINEAR ADVECTION & CORIOLIS TERM/CONSTANT PARTS
        
VAR_NONLIN = 1; 
q_mis   = (VAR_NONLIN*VAR_Curl_vxmis*Vk + f); 
q_pls   = (VAR_NONLIN*VAR_Curl_vxpls*Vk + f); 
nonterm = (q_mis.*((Var_Ac2i_iimis*Hk).*(VAR_Rec_iimis*Vk) + (Var_Ac2i_jjmis*Hk).*(VAR_Rec_jjmis*Vk)) ...
                  - q_pls.*((Var_Ac2i_iipls*Hk).*(VAR_Rec_iipls*Vk) + (Var_Ac2i_jjpls*Hk).*(VAR_Rec_jjpls*Vk)))./edge_d_length;


Advk     = nonterm./(Ac2i*Hk) - VAR_NONLIN*(     (VAR_SWE_Ekin*((Mass_kin*Vk).*Vk) )     ) + dissipation;
Advs = Advk;

%--- CONSTANT FLUX TERM                          

        Fk_RHS = Mass_flat*Vk + 0.5*delt*(1*Advk + 0*GrHk);
        Fk     = Mass_flat\Fk_RHS;

% add noise to height field
Hs = Hs + delt*(Av2c*Ai2v)*(alphaH^2*diffusionH-alphaH*hn_noise);        
     

     
%%%%%%   CALCULATION OF VALUES V*,H* FOR k+1 %%%%%%

%--- FIRST TIME: NOT CONVERGED YET
oldVs=Vs; oldHs=Hs; k = 1;

currentNorm=1;
%--- FIXED POINT LOOP; FIXED POINT ITERATION
 while currentNorm > tol

    %--- KEEP OLD Vs IN MEMORY
    oldVs=Vs; oldHs=Hs; 

    %--- UPDATE OF VELOCITY V* AT t*:

            Vs_RHS = Mass_flat*Fk + 0.5*delt*(1*Advs + 2*GrHs);
            Vs     = Mass_flat\Vs_RHS+0.5*delt*(alphaV^2*diffusionV-alphaV*vn_noise);

    %--- CONTINUITY EQUATION UPDATE AT t*                          
    DvHs   = - Div*((Ac2i*Hs).*Vs);
    Hs     = Ck + 0.5*delt*1*DvHs;

    %--- GRADIENT TERM UPDATE AT t*
    GrHs = - g*Dnor*(Hs + Bk); 

    %--- NONLINEAR ADVECTION/CORIOLIS UPDATE AT t*
    
            VAR_NONLIN = 1; 
            q_mis   = (VAR_NONLIN*VAR_Curl_vxmis*Vs + f);  q_pls   = (VAR_NONLIN*VAR_Curl_vxpls*Vs + f); 
            nonterm = (q_mis.*((Var_Ac2i_iimis*Hs).*(VAR_Rec_iimis*Vs) + (Var_Ac2i_jjmis*Hs).*(VAR_Rec_jjmis*Vs)) ...
                     - q_pls.*((Var_Ac2i_iipls*Hs).*(VAR_Rec_iipls*Vs) + (Var_Ac2i_jjpls*Hs).*(VAR_Rec_jjpls*Vs)))./edge_d_length;
      
            Advs    = nonterm./(Ac2i*Hs) - VAR_NONLIN*(     (VAR_SWE_Ekin*((Mass_kin*Vs).*Vs) )     ) + dissipation;


    k = k + 1; dd= norm(oldVs-Vs,inf); ee= norm(oldHs-Hs,inf);
    currentNorm=norm(oldVs-Vs)+norm(oldHs-Hs,inf);
 end     
 
end
   


function ux = gradx(u)

    global Ai2v Av2c tangential_tx normal_tx  Dnor Dtan
    
    hex_u = Ai2v*u;
    tri_u = Av2c*hex_u;
    
    ux = (Dnor*tri_u(:)).*normal_tx(:)+(Dtan*hex_u(:)).*tangential_tx(:);

end

function uy = grady(u)

    global Ai2v Av2c tangential_ty normal_ty Dnor Dtan
    
    hex_u = Ai2v*u;
    tri_u = Av2c*hex_u;
    
    uy = (Dnor*tri_u(:)).*normal_ty(:)+(Dtan*hex_u(:)).*tangential_ty(:);

end

function uxx = gradxx(u)

    global Ai2v Av2c tangential_tx normal_tx  Dnor Dtan
    
    hex_u = Ai2v*u;
    tri_u = Av2c*hex_u;
    
    ux = (Dnor*tri_u(:)).*normal_tx(:)+(Dtan*hex_u(:)).*tangential_tx(:);

    hex_u = Ai2v*ux;
    tri_u = Av2c*hex_u;
    
    uxx = (Dnor*tri_u(:)).*normal_tx(:)+(Dtan*hex_u(:)).*tangential_tx(:);
    
end

function uyy = gradyy(u)

    global Ai2v Av2c tangential_ty normal_ty Dnor Dtan
    
    hex_u = Ai2v*u;
    tri_u = Av2c*hex_u;
    
    uy = (Dnor*tri_u(:)).*normal_ty(:)+(Dtan*hex_u(:)).*tangential_ty(:);

    hex_u = Ai2v*uy;
    tri_u = Av2c*hex_u;
    
    uyy = (Dnor*tri_u(:)).*normal_ty(:)+(Dtan*hex_u(:)).*tangential_ty(:);
    
end

function diffusion = diffuVnew(u,v)

global normal_tx normal_ty  a0

 diffusion =0.5*a0*(gradxx(u)+gradyy(u)).*normal_tx+0.5*a0*(gradxx(v)+gradyy(v)).*normal_ty;



end

function diffusion = diffuHnew(h)

global a0

 diffusion = 0.5*a0*(gradxx(h)+gradyy(h));

end

function [sdBdt] = createNoise(N)

global sigma_sqdt  nvertices Av2i nedges

% Fourier transform of white noise
dB_sqdt = fft2(randn([N,N]));

% Multiplication by the Fourier transform of the kernel
ft_sigma_dBdt = bsxfun(@times, sigma_sqdt, dB_sqdt);

% Homogeneous velocity field
sigma_dBdt = real(ifft2(ft_sigma_dBdt));

sdBdt=zeros(nedges,2);

sdBdt(:,1)=Av2i*reshape(sigma_dBdt(:,:,1),nvertices,1);
sdBdt(:,2)=Av2i*reshape(sigma_dBdt(:,:,2),nvertices,1);

end
 
     
