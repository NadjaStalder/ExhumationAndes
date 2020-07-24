MODULE definitions 
      
      type parm
         integer n,m_max,napt,nzirc,narar,dummy,contr,iterM
         double precision edot_mean,deltat,xL,sigma2,sigmaD,t_total,lat1,lat2,lon1
         double precision xmu
         double precision bot,lon2,dlon,kappa,hp,Ts,Tb,zl,angle,aspect
         character :: run*5,topofile*100
      end type

      type matr
         integer,dimension(:),allocatable ::B,nsystems
         double precision,dimension(:),allocatable :: isys,edot,edot_pr,edot_pr2,zc,zcp,eps,elev,eps_dat,eps_res,work,syn_age
         integer,dimension(:),allocatable :: iblock,idum_block
         double precision,dimension(:),allocatable :: x_true,y_true,x_dum_true,y_dum_true,elev_true 
         double precision,dimension(:),allocatable :: BB,edot_dat,edot_pr_dum,Cmm,eps_dum,x_dum,y_dum,misfits
         double precision,dimension(:),allocatable :: R,kernel,sf
         double precision,dimension(:),allocatable :: ta,a_error,x,y,cov_eps,depth,depth1,zz,tsteps,tsteps_sum
         double precision,dimension(:,:),allocatable :: A,G,cov,cpost,Y1,Y2,Y3,H,II,distm,ages,Y9
      end type

end MODULE definitions

