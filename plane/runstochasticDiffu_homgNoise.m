function [] = runstochasticDiffu_homgNoise(runNr,DT)

global Div Ac2i Av2i Dnor Curl edge_d_length  Ai2v Av2c tri_A
global g f
global VAR_Curl_vxmis VAR_Curl_vxpls
global Var_Ac2i_iimis Var_Ac2i_iipls Var_Ac2i_jjmis Var_Ac2i_jjpls
global VAR_Rec_iimis  VAR_Rec_iipls  VAR_Rec_jjmis  VAR_Rec_jjpls
global VAR_Rec_vxmis  VAR_Rec_vxpls
global Rec_t2n_ham
global kinE_loc VAR_SWE_Dnor VAR_SWE_Ekin
global Idneig Rec_t2n_thub tangential_tx normal_tx tangential_ty normal_ty intersect_x intersect_y nablasq_tri nedges
global Dtan kinE_loc_irreg Avrg_c2v sigma_sqdt  nvertices Av2i nedges a0
%%
myfolder = ['homNoise/',num2str(DT),'/run-',num2str(runNr)];
N=128;
load(['init-',num2str(N),'.mat'])
delt=delt/4;
dissipation=0;
tol=1e-6;
%%
delt=1/DT*delt;
ptime  = 1;%0.1;
iptime = int32(ptime/delt) ;
days   = 2;
nsteps = days/delt;
%%



load(['homNoise-',num2str(N),'.mat'])

scale=500;
sigma_sqdt = sigma_sqdt*sqrt(scale);
a0=a0_dt*delt*scale;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% predefine diagnostic arrays:

mass                = zeros(int32(nsteps+1),1);
energy_pot          = zeros(int32(nsteps+1),1);
energy_kin_u        = zeros(int32(nsteps+1),1);
energy_tot          = zeros(int32(nsteps+1),1);

v_prog=vn;
h_prog=h_cell;
v_prog_old  = v_prog;


for i = 1:nsteps
    
    
    [Vs,Hs,k,dd,ee,sdBdt] =  ops_st_diffu_homN(v_prog,v_prog_old,h_prog,B_topography,delt,speye(nedges), speye(nedges),dissipation,tol,N);
    
    disp([i*delt,k])
    
    v_prog_old  = v_prog;
    h_prog      = 0*H +  Hs;
    v_prog      = Vs;
    
    %% update diagnostics
    h_cell = h_prog;
    
    mass(i+1)           =   sum(tri_A.*h_cell);
    K                     =  kinE_loc_irreg*(v_prog.*v_prog)  ;
    energy_kin_u(i+1)     =  sum(tri_A.*K.*h_cell);
    energy_pot(i+1)       =  sum(tri_A.*(0.5*g*(h_cell + B_topography).^2));
    energy_tot(i+1)      =    energy_pot(i+1) + energy_kin_u(i+1);
    
    
end

output = fullfile(myfolder,'final');
save(output);

end
