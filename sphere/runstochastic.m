function runstochastic(runNR)

global nu Ai2v Av2c Dnor Dtan Div_tri Curl_hex nor_x nor_y nor_z tang_x tang_y tang_z rec_u rec_v rec_w nedges edge_length_tri intersect_x intersect_y intersect_z Idneig Ac2i R
global tke

n=10242;

testcase='Galewsky';

% skale noise to same strength as POD
load('POD_10242.mat','S');
tke=sum(S.^2)/2;
clear S;


load([num2str(n),'.mat'])
load(['newops-',num2str(n),'.mat'])
load(['neighMatrix_',num2str(n),'.mat'])

disp(['Testcase: ',testcase]);
run([testcase,'_testcase.m'])

t = 20*(60*60*24);
dt=50;
steps = round(t/dt);

%% compute viscosity parameter

u = Ac2i*(rec_u*Vk);
v = Ac2i*(rec_v*Vk);

hex_u = Ai2v*u;
tri_u = Av2c*hex_u;
ux = (Dnor*tri_u(:)).*nor_x(:)+(Dtan*hex_u(:)).*tang_x(:);
uy = (Dnor*tri_u(:)).*nor_y(:)+(Dtan*hex_u(:)).*tang_y(:);
hex_u = Ai2v*v;
tri_u = Av2c*hex_u;
vx = (Dnor*tri_u(:)).*nor_x(:)+(Dtan*hex_u(:)).*tang_x(:);
vy = (Dnor*tri_u(:)).*nor_y(:)+(Dtan*hex_u(:)).*tang_y(:);

RMS = sqrt( mean(sqrt(((ux+vy)/2).^2+((uy+vx)/2).^2).^2) );
nu=40*RMS*mean(edge_length_hex).^4;
nu=nu/3;


%% initialize diagnostics
Ekin = zeros(steps,1);
Epot = zeros(steps,1);
E_tot= zeros(steps,1);
Mass = zeros(steps,1);



%% other parameters
plt_stp=5000;
sec2day = 1.15741e-5;
idxInM = 1;

hex2tri=0;
Vtk=0;

load(['newops-',num2str(n),'.mat'])
load(['neighMatrix_',num2str(n),'.mat'])

%% start time-stepping
disp('start timestepping:')

for tt=1:steps
    
    
    [Vs,Hs,k,dd,ee] = timestep_stochastic(Dij,Hinit,Vtk,Div_hex,hex2tri,tt,Vk,Vk_old,Hk,Hk_0,Bk,dt,dissipation,tol,Div_tri,Ac2i,Dnor,VAR_Curl_vxmis,VAR_Curl_vxpls,Djim,Dijm,Dijp,Djip,WeightIp,WeightIm,WeightJp,WeightJm,edge_length_hex,VAR_SWE_Ekin,g,f_pls,f_mis);
    

    Vk_old  = Vk;
    Hk_old  = Hk;
    Hk      = Hs;
    Vk      = Vs;
    
    
    %% update diagnostics
    Ekin(tt) = sum(tri_A.*((kinE_loc_irreg*(Vk.^2)).*Hk));
    Epot(tt) = sum(tri_A.*(0.5*g*(Hk + Bk).^2));
    E_tot(tt)    = Ekin(tt) + Epot(tt);
    Mass(tt) = sum(tri_A.*Hk);
    %%
    
    
    if(mod(tt*dt,60*60*6)==0) % save solution every 6 hours
        save(['SVD/run-',num2str(runNR),'/',num2str(round(dt*tt/(60*60*6))),'.mat'],'Hk','Vk')
        
    end
    
    
end

end
