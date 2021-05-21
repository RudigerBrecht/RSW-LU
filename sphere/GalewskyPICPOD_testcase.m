
CA_H=[10000 13000];

%% Define test case parameters
g = 9.80616;% m/s^2
Omega = 7.292*1e-5;% 1/s

h0 = 10000;% m

alpha = 0;

%% timestep
t  = 15*(24*60*60);%15 days in sec
dt = 50;% sec
steps = t/dt;

%%
u=uvelocity(nedges,e_lat,2*pi/14,5*pi/14);
v=0*u;
%% PIC
load(['POD_',num2str(n),'.mat']);
gamma = bsxfun(@times, randn(numel(S),1), S); 
sdBdt = reshape(U*gamma, [nedges 3]);
%%
f_pls=zeros(nhex,1);
f_mis=zeros(nhex,1);
ff_mis=zeros(nedges,1);
ff_pls=zeros(nedges,1);

for edge = 1:nedges

    x=intersect_x(edge);
    y=intersect_y(edge);
    z=intersect_z(edge);
    
    P=1/R^2*[R^2-x^2 , -y*x   , -x*z  ;...
               -x*y  ,R^2-y^2 ,-y*z   ;...
               -x*z  ,-z*y    , R^2-z^2 ];
    
    sdBdt(edge,:)=P*sdBdt(edge,:)'; % project noise field to the tangent plane of the sphere.
    
    
end
%%
u=u+sdBdt(:,1)/dt; % add noise to the initial velocity field. 
%%

for edge = 1:nedges    
     
    x=intersect_x(edge);
    y=intersect_y(edge);
    z=intersect_z(edge);


      A=[cos(e_lat(edge))*cos(e_lon(edge)),cos(e_lat(edge))*sin(e_lon(edge)),sin(e_lat(edge));...
        -sin(e_lat(edge))*cos(e_lon(edge)),-sin(e_lat(edge))*sin(e_lon(edge)),cos(e_lat(edge));...
        -sin(e_lon(edge)),cos(e_lon(edge)),0];

    
    P=1/R^2*[R^2-x^2 , -y*x , -x*z  ;...
        -x*y  ,R^2-y^2 ,-y*z   ;...
        -x*z  ,-z*y  , R^2-z^2 ];
    
    V=P*(inv(A)*[0;v(edge);u(edge)]);
    VV(edge,:)=V;
    Vk(edge) = dot(V,[nor_x(edge),nor_y(edge),nor_z(edge)]);
    
   %% find the two neighbouring dual cells for the edge. Then compute the Coriolis force value at the dual cells.    
    [~,~,~,~,~,~,~,~,~,~,~,~,~,~,Mis,Pls]=getiijj(edge,tri_for_e,edges_for_tri,orientation_nor_tri,orientation_tang_tri,hex_for_e,hex_A_section,edges_for_hex,tri_for_hex);
    
   x=hex_center_x(Mis);
   y=hex_center_y(Mis);
   z=hex_center_z(Mis);
   
   f_mis(edge)=2*Omega*(-cos(hex_lon(Mis)).*cos(hex_lat(Mis))*sin(alpha)+sin(hex_lat(Mis))*cos(alpha));

   x=hex_center_x(Pls);
   y=hex_center_y(Pls);
   z=hex_center_z(Pls);
   
   f_pls(edge)=2*Omega*(-cos(hex_lon(Pls)).*cos(hex_lat(Pls))*sin(alpha)+sin(hex_lat(Pls))*cos(alpha));

    
end

Vk=Vk';

%% compute the Coriolis force at the edges and on each dual cell.  
f_edge=2*Omega*(-cos(e_lon).*cos(e_lat)*sin(alpha)+sin(e_lat)*cos(alpha));
f_hex=2*Omega*(-cos(hex_lon).*cos(hex_lat)*sin(alpha)+sin(hex_lat)*cos(alpha));


%% initialize the water depth
Hk_0=h0;

h0=10157.9;
dtheta=abs(max(tri_lat)-min(tri_lat))/ntris;
hh=zeros(1,ntris);
utri=Av2c*(Ai2v*u);
for jj=1:ntris
    phi=(tri_lat(jj));
    if phi>0
        idp=tri_lat<=phi;% & tri_lat>0;
    else
        idp=phi>=tri_lat & tri_lat<0;
    end
    hh(jj)=h0-1/g*sum(dtheta*R*utri(idp).*(2*Omega*sin(tri_lat(idp))+(tan(tri_lat(idp))).*utri(idp)/R));
end
hh=hh';
Hk=hh+120*cos(tri_lat).*exp(-(tri_lon/(1/3)).^2-((pi/4-tri_lat)/(1/15)).^2);

%% bottom topography
Bk     = zeros(size(Hk));



%%%%%%%%%%%%%%%%%%%%
dissipation = 0;
tol         = 1E-10; %1E-12
%%%%%%%%%%%%%%%%%%%%

Hk_old    = Hk;
Vk_old    = Vk;
Hinit = Hk;
Vinit = Vk;


function u = uvelocity(nedges,e_lat,t0,t1)

u=zeros(nedges,1);




for ii=1:nedges
    
   t=e_lat(ii);
    
   if(t>=t0 && t<=t1)
      
       u(ii)=(80/exp(-4/(t1-t0)^2))*exp(1/((t-t0)*(t-t1)));
       
   end
    
end

end
