       program temp
       integer i,j,nsteps,niter,printstep
       real*8 l,u0,v0,unp1k,un,vn,udelta,tdelta,fdelta,tprint,g,
     & tend,uf(300000),vf(300000),energy(300000),time(300000)
       open(11,file='shape',status='unknown')
c      ************************************************
       l=2.0d0
       g=9.81d0
       u0=3.00d0
       v0=0.00d0
       tend=20.0d0
       niter=10
       tdelta=0.01d0
       tprint=0.1d0
c      prints after every tprint times
       udelta=0.0d0
       nsteps=int(tend/tdelta)
       printstep=nint(tprint/tdelta)
       if(printstep.lt.1)printstep=1
       if(nint(1.0d0*nsteps/printstep).gt.300000)then
        write(11,*) 'Increase size of arrays uf, vf etc.'
        stop
       endif
       print*,'nsteps,niter ',nsteps,niter
       time(1)=0.0d0
       un=u0
       vn=v0
       uf(1)=u0
       vf(1)=v0
       energy(1)=v0**2-3*g*dcos(u0)/l
       unp1k=u0
       do 1 i=2,nint(1.0d0*nsteps/printstep)+1
        time(i)=printstep*(i-1)*tdelta
1      continue
       do 100 i=1,nsteps
         if(unp1k.eq.un)unp1k=0.9d0*un
       do 10 j=1,niter
       fdelta=vn+(un-unp1k)/tdelta+3.0d0*g*tdelta*
     & (dcos(unp1k)-dcos(un))/(4.0d0*(unp1k-un)*l)
       udelta=fdelta/(1.0d0/tdelta+3.0d0*g*tdelta*(dcos(unp1k)-
     & dcos(un)+(unp1k-un)*dsin(unp1k))/(4.0d0*l*(unp1k-un)**2))
       unp1k=unp1k+udelta
       if(dabs(udelta).lt.1.0d-10)go to 20
10     continue
       print*,'Failed to converge'
       stop
20     continue
       vn=2.0d0*(unp1k-un)/tdelta-vn
       if(i.ge.printstep)then
       if(mod(i,printstep).eq.0)then
       uf(nint(1.0d0*i/printstep)+1)=unp1k
       vf(nint(1.0d0*i/printstep)+1)=vn
       energy(nint(1.0d0*i/printstep)+1)=vn**2-3.0d0*g*dcos(unp1k)/l
       endif
       endif
       un=unp1k
100    continue

       do 200 i=1,nint(1.0d0*nsteps/printstep)+1
       write(11,500) time(i),uf(i),vf(i),energy(i)
200    continue
500    format(1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7)
       end
