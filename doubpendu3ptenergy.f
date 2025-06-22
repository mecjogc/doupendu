       program temp
       integer i,j,nsteps,niter,printstep
       real*8 g,L1,L2,m1,m2,th10,th20,v10,v20,th1n,th2n,th1np1k,th2np1k,
     & v1n,v2n,v1np1k,v2np1k,th1delta,th2delta,v1delta,v2delta,
     & tdelta,fdelta1,fdelta2,fdelta3,fdelta4,pi,
     & xi1,tprint,m11,m12,m13,m14,m21,m22,m23,m24,m31,m32,m33,m34,
     & m41,m42,m43,m44,detm,tend,
     & uf1(300000),uf2(300000),vf1(300000),
     & vf2(300000),energy(300000),time(300000)
       open(11,file='shape',status='unknown')
c      ************************************************
       pi=4.0d0*datan(1.0d0)
       xi1=dsqrt(0.6d0)
       g=9.8d0
       m1=1.0d0
       m2=1.0d0
       L1=1.0d0
       L2=1.0d0
       th10=pi/2.0d0
       th20=0.00d0
       v10=0.00d0
       v20=0.00d0
       tend=20.0d0
       niter=10
       tdelta=0.0000001d0
       tprint=0.1d0
c      prints after every tprint times
       th1delta=0.0d0
       th2delta=0.0d0
       v1delta=0.0d0
       v2delta=0.0d0
       nsteps=nint(tend/tdelta)
       printstep=nint(tprint/tdelta)
       if(printstep.lt.1)printstep=1
       if(nint(1.0d0*nsteps/printstep).gt.300000)then
        write(11,*) 'Increase size of arrays uf, vf etc.'
        stop
       endif
       print*,'nsteps,niter ',nsteps,niter
       time(1)=0.0d0
       th1n=th10
       th2n=th20
       v1n=v10
       v2n=v20
       uf1(1)=th10
       uf2(1)=th20
       vf1(1)=v10
       vf2(1)=v20
       energy(1)=1.0d0/6.0d0*(m1*L1**2*v10**2+
     & m2*L2**2*v20**2+3.0d0*m2*L1*v10*(L1*v10+
     & L2*v20*dcos(th10-th20))
     & -3.0d0*(m1*g*L1*dcos(th10)
     & +m2*g*L2*dcos(th20)+2.0d0*m2*g*L1*dcos(th10)))
       th1np1k=th10
       th2np1k=th20
       v1np1k=v10
       v2np1k=v20
       do 1 i=2,nint(1.0d0*nsteps/printstep)+1
        time(i)=printstep*(i-1)*tdelta
1      continue
       do 100 i=1,nsteps
         if(th1np1k.eq.th1n)th1np1k=0.999d0*th1n-0.001d0
         if(th2np1k.eq.th2n)th2np1k=0.999d0*th2n-0.001d0
       do 10 j=1,niter
      m11= 1.0d0/tdelta + (L1*(-2*g*(m1 + 2*m2)*Cos(th1n) + 
     &  2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &  L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &   v1np1k*v2np1k*(Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k))))*
     & ((L1*(v1n + v1np1k)*(-2*g*(m1 + 2*m2)*Cos(th1n) + 
     &       2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &       L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &          v1np1k*v2np1k*(Cos(th1n - th2np1k) - 
     &             Cos(th1np1k - th2np1k)))))/(8.*(th1n - th1np1k)) + 
     &  (L2*m2*(v2n + v2np1k)*(L1*v1np1k*v2np1k*
     &        (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &       L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k)) + 
     &       2*g*(-Cos(th2n) + Cos(th2np1k))))/(8.*(th2n - th2np1k)) + 
     &  (L1*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &       3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &          v2n*Cos(th1np1k - th2np1k)))*
     &     ((-4*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k)/2.) + 
     &            m2*(3*g*Sin((th1n + th1np1k - 2*th2n - 2*th2np1k)/
     &                  2.) + 
     &               L2*(v2n + v2np1k)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &               (3*L1*(-v1n - v1np1k)**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k))/4.)))/
     &        (3.*L1*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((-3*(g*(4*m1 + 5*m2)*
     &                  Sin((th1n + th1np1k - th1n*xi1 + th1np1k*xi1)/
     &                   2.) + 
     &                 m2*(3*g*
     &                   Sin((th1n + th1np1k - 2*th2n - 2*th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + 2*th2n*xi1 - 
     &                  2*th2np1k*xi1)/2.) + 
     &                   L2*(v2n + v2np1k - v2n*xi1 + v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                  th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 + xi1) - v1np1k*(1 + xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4.)))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))) - (3*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k + th1n*xi1 - 
     &                   th1np1k*xi1)/2.) + 
     &                 m2*(L2*(v2n + v2np1k + v2n*xi1 - v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 - xi1) - v1np1k*(1 - xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4. + 
     &                   3*g*Sin((th1n + th1np1k - 2*th2n - 
     &                   2*th2np1k + th1n*xi1 - th1np1k*xi1 - 
     &                   2*th2n*xi1 + 2*th2np1k*xi1)/2.))))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                  th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & )))))/18.))/12. + (L2*m2*(2*L2*(v2n + v2np1k) + 
     &       3*L1*v1n*Cos(th1n - th2n) + 
     &       3*L1*v1np1k*Cos(th1np1k - th2np1k))*
     &     ((4*(3*L1*(m1 + 3*m2)*(v1n + v1np1k)**2*
     &             Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &            (9*L2*m2*(-v2n - v2np1k)**2*
     &               Sin(th1n + th1np1k - th2n - th2np1k))/4. + 
     &            9*g*(m1 + 2*m2)*
     &             Sin(th1n + th1np1k - th2n/2. - th2np1k/2.) - 
     &            3*g*(m1 + 6*m2)*Sin((th2n + th2np1k)/2.)))/
     &        (9.*L2*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((9*g*(m1 + 2*m2)*
     &                Sin(th1n + (th2n*(-1 - xi1))/2. + 
     &                  th1np1k*(1 - xi1) - (th2np1k*(1 - xi1))/2. + 
     &                  th1n*xi1) - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k + th2n*xi1 - th2np1k*xi1)/2.)
     & + 3*L1*(m1 + 3*m2)*(v1n + v1np1k + v1n*xi1 - v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 - xi1) - v2np1k*(1 - xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4.)/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                 th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))) + (3*L1*(m1 + 3*m2)*(v1n + v1np1k - v1n*xi1 + v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 + xi1) - v2np1k*(1 + xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4. - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k - th2n*xi1 + th2np1k*xi1)/2.)
     & + 9*g*(m1 + 2*m2)*Sin(th1n + (th2n*(-1 + xi1))/2. - th1n*xi1 + 
     &                  th1np1k*(1 + xi1) - (th2np1k*(1 + xi1))/2.))/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)
     & ))))/18.))/12.))/
     & (4.*(th1n - th1np1k)**2*((L1**2*
     &     (-2*g*(m1 + 2*m2)*Cos(th1n) + 2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &        L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &           v1np1k*v2np1k*
     &            (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k))))**2)/
     &   (16.*(th1n - th1np1k)**2) + 
     &  (L2**2*m2**2*(2*L2*(v2n + v2np1k) + 3*L1*v1n*Cos(th1n - th2n) + 
     &        3*L1*v1np1k*Cos(th1np1k - th2np1k))**2)/144. + 
     &  (L1**2*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &        3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &           v2n*Cos(th1np1k - th2np1k)))**2)/144. + 
     &  (L2**2*m2**2*(L1*v1np1k*v2np1k*
     &         (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &        L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k))+
     &    2*g*(-Cos(th2n)+Cos(th2np1k)))**2)/(16.*(th2n - th2np1k)**2))
     & ) + (L1*(-2*g*(m1 + 2*m2)*Sin(th1np1k) + 
     &  L2*m2*(v1n*v2n*Sin(th1np1k - th2n) + 
     &     v1np1k*v2np1k*Sin(th1np1k - th2np1k)))*
     & ((L1*(v1n + v1np1k)*(-2*g*(m1 + 2*m2)*Cos(th1n) + 
     &       2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &       L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &          v1np1k*v2np1k*(Cos(th1n - th2np1k) - 
     &             Cos(th1np1k - th2np1k)))))/(8.*(th1n - th1np1k)) + 
     &  (L2*m2*(v2n + v2np1k)*(L1*v1np1k*v2np1k*
     &        (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &       L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k)) +
     &       2*g*(-Cos(th2n) + Cos(th2np1k))))/(8.*(th2n - th2np1k)) + 
     &  (L1*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &       3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &          v2n*Cos(th1np1k - th2np1k)))*
     &     ((-4*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k)/2.) + 
     &            m2*(3*g*Sin((th1n + th1np1k - 2*th2n - 2*th2np1k)/
     &                  2.) + 
     &               L2*(v2n + v2np1k)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &               (3*L1*(-v1n - v1np1k)**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k))/4.)))/
     &        (3.*L1*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((-3*(g*(4*m1 + 5*m2)*
     &                  Sin((th1n + th1np1k - th1n*xi1 + th1np1k*xi1)/
     &                   2.) + 
     &                 m2*(3*g*
     &                   Sin((th1n + th1np1k - 2*th2n - 2*th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + 2*th2n*xi1 - 
     &                  2*th2np1k*xi1)/2.) + 
     &                   L2*(v2n + v2np1k - v2n*xi1 + v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                  th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 + xi1) - v1np1k*(1 + xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4.)))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))) - (3*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k + th1n*xi1 - 
     &                   th1np1k*xi1)/2.) + 
     &                 m2*(L2*(v2n + v2np1k + v2n*xi1 - v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 - xi1) - v1np1k*(1 - xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4. + 
     &                   3*g*Sin((th1n + th1np1k - 2*th2n - 
     &                   2*th2np1k + th1n*xi1 - th1np1k*xi1 - 
     &                   2*th2n*xi1 + 2*th2np1k*xi1)/2.))))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                 th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & )))))/18.))/12. + (L2*m2*(2*L2*(v2n + v2np1k) + 
     &       3*L1*v1n*Cos(th1n - th2n) + 
     &       3*L1*v1np1k*Cos(th1np1k - th2np1k))*
     &     ((4*(3*L1*(m1 + 3*m2)*(v1n + v1np1k)**2*
     &             Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &            (9*L2*m2*(-v2n - v2np1k)**2*
     &               Sin(th1n + th1np1k - th2n - th2np1k))/4. + 
     &            9*g*(m1 + 2*m2)*
     &             Sin(th1n + th1np1k - th2n/2. - th2np1k/2.) - 
     &            3*g*(m1 + 6*m2)*Sin((th2n + th2np1k)/2.)))/
     &        (9.*L2*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((9*g*(m1 + 2*m2)*
     &                Sin(th1n + (th2n*(-1 - xi1))/2. + 
     &                  th1np1k*(1 - xi1) - (th2np1k*(1 - xi1))/2. + 
     &                  th1n*xi1) - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k + th2n*xi1 - th2np1k*xi1)/2.)
     & + 3*L1*(m1 + 3*m2)*(v1n + v1np1k + v1n*xi1 - v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 - xi1) - v2np1k*(1 - xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4.)/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                  th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))) + (3*L1*(m1 + 3*m2)*(v1n + v1np1k - v1n*xi1 + v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 + xi1) - v2np1k*(1 + xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4. - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k - th2n*xi1 + th2np1k*xi1)/2.) 
     & + 9*g*(m1 + 2*m2)*Sin(th1n + (th2n*(-1 + xi1))/2. - th1n*xi1 + 
     &                  th1np1k*(1 + xi1) - (th2np1k*(1 + xi1))/2.))/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                 th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)
     & ))))/18.))/12.))/
     &    (4.*(th1n - th1np1k)*((L1**2*
     &     (-2*g*(m1 + 2*m2)*Cos(th1n) + 2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &        L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &           v1np1k*v2np1k*
     &            (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k))))**2)/
     &   (16.*(th1n - th1np1k)**2) + 
     &  (L2**2*m2**2*(2*L2*(v2n + v2np1k) + 3*L1*v1n*Cos(th1n - th2n) + 
     &        3*L1*v1np1k*Cos(th1np1k - th2np1k))**2)/144. + 
     &  (L1**2*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &        3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &           v2n*Cos(th1np1k - th2np1k)))**2)/144. + 
     &  (L2**2*m2**2*(L1*v1np1k*v2np1k*
     &         (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &    L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k)) + 
     &  2*g*(-Cos(th2n) + Cos(th2np1k)))**2)/(16.*(th2n - th2np1k)**2))
     & )-(L1*(-2*g*(m1+2*m2)*Cos(th1n) + 2*g*(m1 + 2*m2)*Cos(th1np1k)+ 
     &  L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &  v1np1k*v2np1k*(Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k))))*
     & ((L1**2*(-2*g*(m1 + 2*m2)*Cos(th1n) + 
     &        2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &        L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &           v1np1k*v2np1k*
     &            (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k))))**2)/
     &   (8.*(th1n - th1np1k)**3) - 
     &  (L1*L2**2*m2**2*v1np1k*
     &     (2*L2*(v2n + v2np1k) + 3*L1*v1n*Cos(th1n - th2n) + 
     &    3*L1*v1np1k*Cos(th1np1k - th2np1k))*Sin(th1np1k - th2np1k))/
     &   24. - (L1**2*L2*m2*v2n*
     &     (2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &       3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &       v2n*Cos(th1np1k - th2np1k)))*Sin(th1np1k - th2np1k))/24.+
     &  (L1*L2**2*m2**2*v1n*v2n*
     &     (L1*v1np1k*v2np1k*(Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &       L1*v1n*v2n*(Cos(th1np1k - th2n) - 
     &          Cos(th1np1k - th2np1k)) + 
     &       2*g*(-Cos(th2n) + Cos(th2np1k)))*
     &     (-Sin(th1np1k - th2n) + Sin(th1np1k - th2np1k)))/
     &   (8.*(th2n - th2np1k)**2) + 
     &  (L1**2*(-2*g*(m1 + 2*m2)*Cos(th1n) + 
     &       2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &       L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &          v1np1k*v2np1k*
     &           (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k))))*
     &     (-2*g*(m1 + 2*m2)*Sin(th1np1k) + 
     &       L2*m2*(v1n*v2n*Sin(th1np1k - th2n) + 
     &          v1np1k*v2np1k*Sin(th1np1k - th2np1k))))/
     &   (8.*(th1n - th1np1k)**2))*
     & ((L1*(v1n + v1np1k)*(-2*g*(m1 + 2*m2)*Cos(th1n) + 
     &       2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &       L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &          v1np1k*v2np1k*(Cos(th1n - th2np1k) - 
     &             Cos(th1np1k - th2np1k)))))/(8.*(th1n - th1np1k)) + 
     &  (L2*m2*(v2n + v2np1k)*(L1*v1np1k*v2np1k*
     &        (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &       L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k)) +
     &       2*g*(-Cos(th2n) + Cos(th2np1k))))/(8.*(th2n - th2np1k)) + 
     &  (L1*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &       3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &          v2n*Cos(th1np1k - th2np1k)))*
     &     ((-4*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k)/2.) + 
     &            m2*(3*g*Sin((th1n + th1np1k - 2*th2n - 2*th2np1k)/
     &                  2.) + 
     &               L2*(v2n + v2np1k)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &               (3*L1*(-v1n - v1np1k)**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k))/4.)))/
     &        (3.*L1*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((-3*(g*(4*m1 + 5*m2)*
     &                  Sin((th1n + th1np1k - th1n*xi1 + th1np1k*xi1)/
     &                   2.) + 
     &                 m2*(3*g*
     &                   Sin((th1n + th1np1k - 2*th2n - 2*th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + 2*th2n*xi1 - 
     &                  2*th2np1k*xi1)/2.) + 
     &                   L2*(v2n + v2np1k - v2n*xi1 + v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                  th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 + xi1) - v1np1k*(1 + xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4.)))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))) - (3*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k + th1n*xi1 - 
     &                   th1np1k*xi1)/2.) + 
     &                 m2*(L2*(v2n + v2np1k + v2n*xi1 - v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 - xi1) - v1np1k*(1 - xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4. + 
     &                   3*g*Sin((th1n + th1np1k - 2*th2n - 
     &                   2*th2np1k + th1n*xi1 - th1np1k*xi1 - 
     &                   2*th2n*xi1 + 2*th2np1k*xi1)/2.))))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                 th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & )))))/18.))/12. + (L2*m2*(2*L2*(v2n + v2np1k) + 
     &       3*L1*v1n*Cos(th1n - th2n) + 
     &       3*L1*v1np1k*Cos(th1np1k - th2np1k))*
     &     ((4*(3*L1*(m1 + 3*m2)*(v1n + v1np1k)**2*
     &             Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &            (9*L2*m2*(-v2n - v2np1k)**2*
     &               Sin(th1n + th1np1k - th2n - th2np1k))/4. + 
     &            9*g*(m1 + 2*m2)*
     &             Sin(th1n + th1np1k - th2n/2. - th2np1k/2.) - 
     &            3*g*(m1 + 6*m2)*Sin((th2n + th2np1k)/2.)))/
     &        (9.*L2*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((9*g*(m1 + 2*m2)*
     &                Sin(th1n + (th2n*(-1 - xi1))/2. + 
     &                  th1np1k*(1 - xi1) - (th2np1k*(1 - xi1))/2. + 
     &                  th1n*xi1) - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k + th2n*xi1 - th2np1k*xi1)/2.)
     & + 3*L1*(m1 + 3*m2)*(v1n + v1np1k + v1n*xi1 - v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 - xi1) - v2np1k*(1 - xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4.)/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                 th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))) + (3*L1*(m1 + 3*m2)*(v1n + v1np1k - v1n*xi1 + v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 + xi1) - v2np1k*(1 + xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4. - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k - th2n*xi1 + th2np1k*xi1)/2.) 
     & + 9*g*(m1 + 2*m2)*Sin(th1n + (th2n*(-1 + xi1))/2. - th1n*xi1 + 
     &                  th1np1k*(1 + xi1) - (th2np1k*(1 + xi1))/2.))/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)
     & ))))/18.))/12.))/
     &    (4.*(th1n - th1np1k)*((L1**2*
     &      (-2*g*(m1 + 2*m2)*Cos(th1n) + 2*g*(m1 + 2*m2)*Cos(th1np1k)+
     &         L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n))+
     &            v1np1k*v2np1k*
     &             (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k))))**2)/
     &    (16.*(th1n - th1np1k)**2) + 
     &   (L2**2*m2**2*(2*L2*(v2n + v2np1k) + 3*L1*v1n*Cos(th1n - th2n) +
     &         3*L1*v1np1k*Cos(th1np1k - th2np1k))**2)/144. + 
     &   (L1**2*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &         3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &            v2n*Cos(th1np1k - th2np1k)))**2)/144. + 
     &   (L2**2*m2**2*(L1*v1np1k*v2np1k*
     &          (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &         L1*v1n*v2n*(Cos(th1np1k - th2n) - 
     &            Cos(th1np1k - th2np1k)) + 
     &         2*g*(-Cos(th2n) + Cos(th2np1k)))**2)/
     &    (16.*(th2n - th2np1k)**2))**2) + 
     & (L1*(-2*g*(m1 + 2*m2)*Cos(th1n) + 2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &  L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &  v1np1k*v2np1k*(Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k))))*
     & ((L1*(v1n + v1np1k)*(-2*g*(m1 + 2*m2)*Cos(th1n) + 
     &       2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &       L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &          v1np1k*v2np1k*(Cos(th1n - th2np1k) - 
     &             Cos(th1np1k - th2np1k)))))/(8.*(th1n - th1np1k)**2) +
     &  (L1*L2*m2*v1n*v2n*(v2n + v2np1k)*
     &     (-Sin(th1np1k - th2n) + Sin(th1np1k - th2np1k)))/
     &   (8.*(th2n - th2np1k)) + 
     &  (L1*(v1n + v1np1k)*(-2*g*(m1 + 2*m2)*Sin(th1np1k) + 
     &       L2*m2*(v1n*v2n*Sin(th1np1k - th2n) + 
     &          v1np1k*v2np1k*Sin(th1np1k - th2np1k))))/
     &   (8.*(th1n - th1np1k)) - 
     &  (L1*L2*m2*v2n*Sin(th1np1k - th2np1k)*
     &     ((-4*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k)/2.) + 
     &            m2*(3*g*Sin((th1n + th1np1k - 2*th2n - 2*th2np1k)/
     &                  2.) + 
     &               L2*(v2n + v2np1k)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &               (3*L1*(-v1n - v1np1k)**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k))/4.)))/
     &        (3.*L1*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((-3*(g*(4*m1 + 5*m2)*
     &                  Sin((th1n + th1np1k - th1n*xi1 + th1np1k*xi1)/
     &                   2.) + 
     &                 m2*(3*g*
     &                   Sin((th1n + th1np1k - 2*th2n - 2*th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + 2*th2n*xi1 - 
     &                   2*th2np1k*xi1)/2.) + 
     &                   L2*(v2n + v2np1k - v2n*xi1 + v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 + xi1) - v1np1k*(1 + xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4.)))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                 th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))) - (3*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k + th1n*xi1 - 
     &                   th1np1k*xi1)/2.) + 
     &                 m2*(L2*(v2n + v2np1k + v2n*xi1 - v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &                    (3*L1*(v1n*(-1 - xi1) - v1np1k*(1 - xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4. + 
     &                    3*g*
     &                    Sin((th1n + th1np1k - 2*th2n - 2*th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - 2*th2n*xi1 + 
     &                   2*th2np1k*xi1)/2.))))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1)
     & ))))/18.))/4. + (L1*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &       3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &          v2n*Cos(th1np1k - th2np1k)))*
     &     ((-4*((g*(4*m1 + 5*m2)*Cos((th1n + th1np1k)/2.))/2. + 
     &            m2*((3*g*Cos((th1n + th1np1k - 2*th2n - 2*th2np1k)/
     &                   2.))/2. + 
     &               (L2*(v2n + v2np1k)**2*
     &                  Cos((th1n + th1np1k - th2n - th2np1k)/2.))/2. +
     &               (3*L1*(-v1n - v1np1k)**2*
     &                  Cos(th1n + th1np1k - th2n - th2np1k))/4.)))/
     &        (3.*L1*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (12*m2*Sin(th1n + th1np1k - th2n - th2np1k)*
     &          (g*(4*m1 + 5*m2)*Sin((th1n + th1np1k)/2.) + 
     &            m2*(3*g*Sin((th1n + th1np1k - 2*th2n - 2*th2np1k)/
     &                  2.) + 
     &               L2*(v2n + v2np1k)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &               (3*L1*(-v1n - v1np1k)**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k))/4.)))/
     &        (L1*(8*m1 + 15*m2 - 
     &             9*m2*Cos(th1n + th1np1k - th2n - th2np1k))**2) + 
     &       (5*((-3*((g*(4*m1 + 5*m2)*(1 + xi1)*
     &                   Cos((th1n + th1np1k - th1n*xi1 + 
     &                  th1np1k*xi1)/2.))/2. + 
     &                 m2*((3*g*(1 + xi1)*
     &                   Cos((th1n + th1np1k - 2*th2n - 2*th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + 2*th2n*xi1 - 
     &                  2*th2np1k*xi1)/2.))/2. + 
     &                   (L2*(1 + xi1)*
     &                   (v2n + v2np1k - v2n*xi1 + v2np1k*xi1)**2*
     &                   Cos((th1n + th1np1k - th2n - th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                  th2np1k*xi1)/2.))/2. + 
     &                   (3*L1*(1 + xi1)*
     &                   (v1n*(-1 + xi1) - v1np1k*(1 + xi1))**2*
     &                   Cos(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4.)))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                 th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))) - (3*((g*(4*m1 + 5*m2)*(1 - xi1)*
     &                   Cos((th1n + th1np1k + th1n*xi1 - 
     &                  th1np1k*xi1)/2.))/2. + 
     &                 m2*((L2*(1 - xi1)*
     &                   (v2n + v2np1k + v2n*xi1 - v2np1k*xi1)**2*
     &                   Cos((th1n + th1np1k - th2n - th2np1k + 
     &                  th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                  th2np1k*xi1)/2.))/2. + 
     &                   (3*L1*(v1n*(-1 - xi1) - v1np1k*(1 - xi1))**2*
     &                   (1 - xi1)*
     &                   Cos(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4. + 
     &                   (3*g*(1 - xi1)*
     &                   Cos((th1n + th1np1k - 2*th2n - 2*th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - 2*th2n*xi1 + 
     &                   2*th2np1k*xi1)/2.))/2.)))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                 th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))) + (27*m2*(1 + xi1)*Sin(th1n + th1np1k - th2n - th2np1k - 
     &                 th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)*
     &               (g*(4*m1 + 5*m2)*
     &                  Sin((th1n + th1np1k - th1n*xi1 + th1np1k*xi1)/
     &                   2.) + 
     &                 m2*(3*g*
     &                   Sin((th1n + th1np1k - 2*th2n - 2*th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + 2*th2n*xi1 - 
     &                   2*th2np1k*xi1)/2.) + 
     &                   L2*(v2n + v2np1k - v2n*xi1 + v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 + xi1) - v1np1k*(1 + xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4.)))/
     &             (L1*(8*m1 + 15*m2 - 
     &                  9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))**2) + (27*m2*(1 - xi1)*Sin(th1n + th1np1k - th2n - th2np1k + 
     &                 th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1)*
     &               (g*(4*m1 + 5*m2)*
     &                  Sin((th1n + th1np1k + th1n*xi1 - th1np1k*xi1)/
     &                   2.) + 
     &                 m2*(L2*(v2n + v2np1k + v2n*xi1 - v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &                    (3*L1*(v1n*(-1 - xi1) - v1np1k*(1 - xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4. + 
     &                    3*g*
     &                    Sin((th1n + th1np1k - 2*th2n - 2*th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - 2*th2n*xi1 + 
     &                   2*th2np1k*xi1)/2.))))/
     &             (L1*(8*m1 + 15*m2 - 
     &                  9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                 th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))**2)))/18.))/12. - (L1*L2*m2*v1np1k*Sin(th1np1k - th2np1k)*
     &     ((4*(3*L1*(m1 + 3*m2)*(v1n + v1np1k)**2*
     &             Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &            (9*L2*m2*(-v2n - v2np1k)**2*
     &               Sin(th1n + th1np1k - th2n - th2np1k))/4. + 
     &            9*g*(m1 + 2*m2)*
     &             Sin(th1n + th1np1k - th2n/2. - th2np1k/2.) - 
     &            3*g*(m1 + 6*m2)*Sin((th2n + th2np1k)/2.)))/
     &        (9.*L2*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((9*g*(m1 + 2*m2)*
     &                Sin(th1n + (th2n*(-1 - xi1))/2. + 
     &                  th1np1k*(1 - xi1) - (th2np1k*(1 - xi1))/2. + 
     &                  th1n*xi1) - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k + th2n*xi1 - th2np1k*xi1)/2.)
     & + 3*L1*(m1 + 3*m2)*(v1n + v1np1k + v1n*xi1 - v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 - xi1) - v2np1k*(1 - xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4.)/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                 th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))) + (3*L1*(m1 + 3*m2)*(v1n + v1np1k - v1n*xi1 + v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 + xi1) - v2np1k*(1 + xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4. - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k - th2n*xi1 + th2np1k*xi1)/2.)
     & + 9*g*(m1 + 2*m2)*Sin(th1n + (th2n*(-1 + xi1))/2. - th1n*xi1 + 
     &                  th1np1k*(1 + xi1) - (th2np1k*(1 + xi1))/2.))/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)
     & ))))/18.))/4.+(L2*m2*(2*L2*(v2n+v2np1k)+3*L1*v1n*Cos(th1n-th2n)+
     &       3*L1*v1np1k*Cos(th1np1k - th2np1k))*
     &     ((4*((3*L1*(m1 + 3*m2)*(v1n + v1np1k)**2*
     &               Cos((th1n + th1np1k - th2n - th2np1k)/2.))/2. + 
     &            (9*L2*m2*(-v2n - v2np1k)**2*
     &               Cos(th1n + th1np1k - th2n - th2np1k))/4. + 
     &            9*g*(m1 + 2*m2)*
     &             Cos(th1n + th1np1k - th2n/2. - th2np1k/2.)))/
     &        (9.*L2*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) - 
     &       (4*m2*Sin(th1n + th1np1k - th2n - th2np1k)*
     &          (3*L1*(m1 + 3*m2)*(v1n + v1np1k)**2*
     &             Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &            (9*L2*m2*(-v2n - v2np1k)**2*
     &               Sin(th1n + th1np1k - th2n - th2np1k))/4. + 
     &            9*g*(m1 + 2*m2)*
     &             Sin(th1n + th1np1k - th2n/2. - th2np1k/2.) - 
     &            3*g*(m1 + 6*m2)*Sin((th2n + th2np1k)/2.)))/
     &        (L2*(8*m1 + 15*m2 - 
     &             9*m2*Cos(th1n + th1np1k - th2n - th2np1k))**2) + 
     &       (5*((9*g*(m1 + 2*m2)*(1 - xi1)*
     &                Cos(th1n + (th2n*(-1 - xi1))/2. + 
     &                  th1np1k*(1 - xi1) - (th2np1k*(1 - xi1))/2. + 
     &                  th1n*xi1) + 
     &               (3*L1*(m1 + 3*m2)*(1 - xi1)*
     &                  (v1n + v1np1k + v1n*xi1 - v1np1k*xi1)**2*
     &                  Cos((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.))/2. + 
     &               (9*L2*m2*(v2n*(-1 - xi1) - v2np1k*(1 - xi1))**2*
     &                  (1 - xi1)*
     &                  Cos(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))/4.)/(L2*(8*m1 + 15*m2 - 9*m2*
     &                  Cos(th1n + th1np1k - th2n - th2np1k + th1n*xi1 -
     &                    th1np1k*xi1 - th2n*xi1 + th2np1k*xi1))) + 
     &            ((3*L1*(m1 + 3*m2)*(1 + xi1)*
     &                  (v1n + v1np1k - v1n*xi1 + v1np1k*xi1)**2*
     &                  Cos((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.))/2. + 
     &               (9*L2*m2*(1 + xi1)*
     &                  (v2n*(-1 + xi1) - v2np1k*(1 + xi1))**2*
     &                  Cos(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4. + 
     &               9*g*(m1 + 2*m2)*(1 + xi1)*
     &                Cos(th1n + (th2n*(-1 + xi1))/2. - th1n*xi1 + 
     &                  th1np1k*(1 + xi1) - (th2np1k*(1 + xi1))/2.))/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)
     & )) - (9*m2*(1-xi1)*Sin(th1n+th1np1k-th2n - th2np1k + th1n*xi1-
     &                 th1np1k*xi1 - th2n*xi1 + th2np1k*xi1)*
     &               (9*g*(m1 + 2*m2)*
     &                  Sin(th1n + (th2n*(-1 - xi1))/2. + 
     &                   th1np1k*(1 - xi1) - (th2np1k*(1 - xi1))/2. + 
     &                   th1n*xi1) - 
     &                 3*g*(m1 + 6*m2)*
     &                  Sin((th2n + th2np1k + th2n*xi1 - th2np1k*xi1)/
     &                   2.) + 
     &                 3*L1*(m1 + 3*m2)*
     &                  (v1n + v1np1k + v1n*xi1 - v1np1k*xi1)**2*
     &                  Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &                 (9*L2*m2*(v2n*(-1 - xi1) - v2np1k*(1 - xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))/4.))/(L2*(8*m1 + 15*m2 - 9*m2*
     &                   Cos(th1n + th1np1k - th2n - th2np1k + 
     &                  th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))**2) - (9*m2*(1 + xi1)*Sin(th1n + th1np1k - th2n - th2np1k - 
     &                 th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)*
     &               (3*L1*(m1 + 3*m2)*
     &                  (v1n + v1np1k - v1n*xi1 + v1np1k*xi1)**2*
     &                  Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & )/2.) + (9*L2*m2*(v2n*(-1 + xi1) - v2np1k*(1 + xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))/4. - 3*g*(m1 + 6*m2)*Sin((th2n+th2np1k-th2n*xi1+th2np1k*xi1)/
     &                   2.) + 
     &                 9*g*(m1 + 2*m2)*
     &                  Sin(th1n + (th2n*(-1 + xi1))/2. - th1n*xi1 + 
     &                    th1np1k*(1 + xi1) - (th2np1k*(1 + xi1))/2.)))/
     &             (L2*(8*m1 + 15*m2 - 
     &                  9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)
     & )**2)))/18.))/12.))/
     &    (4.*(th1n - th1np1k)*((L1**2*
     &     (-2*g*(m1 + 2*m2)*Cos(th1n) + 2*g*(m1 + 2*m2)*Cos(th1np1k) +
     &        L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) +
     &           v1np1k*v2np1k*
     &            (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k))))**2)/
     &   (16.*(th1n - th1np1k)**2) + 
     &  (L2**2*m2**2*(2*L2*(v2n + v2np1k) + 3*L1*v1n*Cos(th1n - th2n) + 
     &        3*L1*v1np1k*Cos(th1np1k - th2np1k))**2)/144. + 
     &  (L1**2*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &        3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &           v2n*Cos(th1np1k - th2np1k)))**2)/144. + 
     &  (L2**2*m2**2*(L1*v1np1k*v2np1k*
     &         (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &       L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k)) +
     &  2*g*(-Cos(th2n) + Cos(th2np1k)))**2)/(16.*(th2n - th2np1k)**2)))

      m12=(L1*L2*m2*v1np1k*v2np1k*(Sin(th1n-th2np1k)
     & -Sin(th1np1k-th2np1k))*
     & ((L1*(v1n + v1np1k)*(-2*g*(m1 + 2*m2)*Cos(th1n) + 
     &       2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &       L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &          v1np1k*v2np1k*(Cos(th1n - th2np1k) - 
     &             Cos(th1np1k - th2np1k)))))/(8.*(th1n - th1np1k)) + 
     &  (L2*m2*(v2n + v2np1k)*(L1*v1np1k*v2np1k*
     &        (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &       L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k)) +
     &       2*g*(-Cos(th2n) + Cos(th2np1k))))/(8.*(th2n - th2np1k)) + 
     &  (L1*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &       3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &          v2n*Cos(th1np1k - th2np1k)))*
     &     ((-4*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k)/2.) + 
     &            m2*(3*g*Sin((th1n + th1np1k - 2*th2n - 2*th2np1k)/
     &                  2.) + 
     &               L2*(v2n + v2np1k)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &               (3*L1*(-v1n - v1np1k)**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k))/4.)))/
     &        (3.*L1*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((-3*(g*(4*m1 + 5*m2)*
     &                  Sin((th1n + th1np1k - th1n*xi1 + th1np1k*xi1)/
     &                   2.) + 
     &                 m2*(3*g*
     &                   Sin((th1n + th1np1k - 2*th2n - 2*th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + 2*th2n*xi1 - 
     &                  2*th2np1k*xi1)/2.) + 
     &                   L2*(v2n + v2np1k - v2n*xi1 + v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                  th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 + xi1) - v1np1k*(1 + xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4.)))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))) - (3*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k + th1n*xi1 - 
     &                   th1np1k*xi1)/2.) + 
     &                 m2*(L2*(v2n + v2np1k + v2n*xi1 - v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 - xi1) - v1np1k*(1 - xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4. + 
     &                   3*g*Sin((th1n + th1np1k - 2*th2n - 
     &                   2*th2np1k + th1n*xi1 - th1np1k*xi1 - 
     &                   2*th2n*xi1 + 2*th2np1k*xi1)/2.))))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & )))))/18.))/12. + (L2*m2*(2*L2*(v2n + v2np1k) + 
     &       3*L1*v1n*Cos(th1n - th2n) + 
     &       3*L1*v1np1k*Cos(th1np1k - th2np1k))*
     &     ((4*(3*L1*(m1 + 3*m2)*(v1n + v1np1k)**2*
     &             Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &            (9*L2*m2*(-v2n - v2np1k)**2*
     &               Sin(th1n + th1np1k - th2n - th2np1k))/4. + 
     &            9*g*(m1 + 2*m2)*
     &             Sin(th1n + th1np1k - th2n/2. - th2np1k/2.) - 
     &            3*g*(m1 + 6*m2)*Sin((th2n + th2np1k)/2.)))/
     &        (9.*L2*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((9*g*(m1 + 2*m2)*
     &                Sin(th1n + (th2n*(-1 - xi1))/2. + 
     &                  th1np1k*(1 - xi1) - (th2np1k*(1 - xi1))/2. + 
     &                  th1n*xi1) - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k + th2n*xi1 - th2np1k*xi1)/2.)
     & + 3*L1*(m1 + 3*m2)*(v1n + v1np1k + v1n*xi1 - v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 - xi1) - v2np1k*(1 - xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4.)/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                  th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))) + (3*L1*(m1 + 3*m2)*(v1n + v1np1k - v1n*xi1 + v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 + xi1) - v2np1k*(1 + xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4. - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k - th2n*xi1 + th2np1k*xi1)/2.)
     & + 9*g*(m1 + 2*m2)*Sin(th1n + (th2n*(-1 + xi1))/2. - th1n*xi1 + 
     &                  th1np1k*(1 + xi1) - (th2np1k*(1 + xi1))/2.))/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)
     & ))))/18.))/12.))/
     &    (4.*(th1n - th1np1k)*((L1**2*
     &     (-2*g*(m1 + 2*m2)*Cos(th1n) + 2*g*(m1 + 2*m2)*Cos(th1np1k) +
     &        L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) +
     &           v1np1k*v2np1k*
     &            (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k))))**2)/
     &   (16.*(th1n - th1np1k)**2) + 
     &  (L2**2*m2**2*(2*L2*(v2n + v2np1k) + 3*L1*v1n*Cos(th1n - th2n) +
     &        3*L1*v1np1k*Cos(th1np1k - th2np1k))**2)/144. + 
     &  (L1**2*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &        3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &           v2n*Cos(th1np1k - th2np1k)))**2)/144. + 
     &  (L2**2*m2**2*(L1*v1np1k*v2np1k*
     &         (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &        L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k-th2np1k))+
     & 2*g*(-Cos(th2n) + Cos(th2np1k)))**2)/(16.*(th2n - th2np1k)**2))
     & ) - (L1*(-2*g*(m1 + 2*m2)*Cos(th1n) + 2*g*(m1+2*m2)*Cos(th1np1k)+
     &  L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &     v1np1k*v2np1k*(Cos(th1n - th2np1k) - Cos(th1np1k-th2np1k))))*
     & ((L2**2*m2**2*(L1*v1np1k*v2np1k*
     &         (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &        L1*v1n*v2n*(Cos(th1np1k - th2n) - 
     &           Cos(th1np1k - th2np1k)) + 
     &        2*g*(-Cos(th2n) + Cos(th2np1k)))**2)/
     &   (8.*(th2n - th2np1k)**3) + 
     &  (L1**2*L2*m2*v1np1k*v2np1k*
     &     (-2*g*(m1 + 2*m2)*Cos(th1n) + 2*g*(m1 + 2*m2)*Cos(th1np1k) +
     &       L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) +
     &          v1np1k*v2np1k*
     &           (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k))))*
     &     (Sin(th1n - th2np1k) - Sin(th1np1k - th2np1k)))/
     &   (8.*(th1n - th1np1k)**2) + 
     &  (L1*L2**2*m2**2*v1np1k*
     &     (2*L2*(v2n + v2np1k) + 3*L1*v1n*Cos(th1n - th2n) + 
     &    3*L1*v1np1k*Cos(th1np1k - th2np1k))*Sin(th1np1k - th2np1k))/
     &   24. + (L1**2*L2*m2*v2n*
     &     (2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &       3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &        v2n*Cos(th1np1k - th2np1k)))*Sin(th1np1k - th2np1k))/24.+
     &  (L2**2*m2**2*(L1*v1np1k*v2np1k*
     &        (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &       L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k)) +
     &       2*g*(-Cos(th2n) + Cos(th2np1k)))*
     &     (-(L1*v1np1k*v2np1k*Sin(th1n - th2np1k)) - 
     &       L1*v1n*v2n*Sin(th1np1k - th2np1k) - 2*g*Sin(th2np1k)))/
     &   (8.*(th2n - th2np1k)**2))*
     & ((L1*(v1n + v1np1k)*(-2*g*(m1 + 2*m2)*Cos(th1n) + 
     &       2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &       L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &          v1np1k*v2np1k*(Cos(th1n - th2np1k) - 
     &             Cos(th1np1k - th2np1k)))))/(8.*(th1n - th1np1k)) + 
     &  (L2*m2*(v2n + v2np1k)*(L1*v1np1k*v2np1k*
     &        (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &       L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k)) +
     &       2*g*(-Cos(th2n) + Cos(th2np1k))))/(8.*(th2n - th2np1k)) + 
     &  (L1*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &       3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &          v2n*Cos(th1np1k - th2np1k)))*
     &     ((-4*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k)/2.) + 
     &            m2*(3*g*Sin((th1n + th1np1k - 2*th2n - 2*th2np1k)/
     &                  2.) + 
     &               L2*(v2n + v2np1k)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &               (3*L1*(-v1n - v1np1k)**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k))/4.)))/
     &        (3.*L1*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((-3*(g*(4*m1 + 5*m2)*
     &                  Sin((th1n + th1np1k - th1n*xi1 + th1np1k*xi1)/
     &                   2.) + 
     &                 m2*(3*g*
     &                   Sin((th1n + th1np1k - 2*th2n - 2*th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + 2*th2n*xi1 - 
     &                  2*th2np1k*xi1)/2.) + 
     &                   L2*(v2n + v2np1k - v2n*xi1 + v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                  th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 + xi1) - v1np1k*(1 + xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4.)))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))) - (3*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k + th1n*xi1 - 
     &                   th1np1k*xi1)/2.) + 
     &                 m2*(L2*(v2n + v2np1k + v2n*xi1 - v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 - xi1) - v1np1k*(1 - xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4. + 
     &                   3*g*Sin((th1n + th1np1k - 2*th2n - 
     &                   2*th2np1k + th1n*xi1 - th1np1k*xi1 - 
     &                   2*th2n*xi1 + 2*th2np1k*xi1)/2.))))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &             th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & )))))/18.))/12. + (L2*m2*(2*L2*(v2n + v2np1k) + 
     &       3*L1*v1n*Cos(th1n - th2n) + 
     &       3*L1*v1np1k*Cos(th1np1k - th2np1k))*
     &     ((4*(3*L1*(m1 + 3*m2)*(v1n + v1np1k)**2*
     &             Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &            (9*L2*m2*(-v2n - v2np1k)**2*
     &               Sin(th1n + th1np1k - th2n - th2np1k))/4. + 
     &            9*g*(m1 + 2*m2)*
     &             Sin(th1n + th1np1k - th2n/2. - th2np1k/2.) - 
     &            3*g*(m1 + 6*m2)*Sin((th2n + th2np1k)/2.)))/
     &        (9.*L2*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((9*g*(m1 + 2*m2)*
     &                Sin(th1n + (th2n*(-1 - xi1))/2. + 
     &                  th1np1k*(1 - xi1) - (th2np1k*(1 - xi1))/2. + 
     &                  th1n*xi1) - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k + th2n*xi1 - th2np1k*xi1)/2.)
     & + 3*L1*(m1 + 3*m2)*(v1n + v1np1k + v1n*xi1 - v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 - xi1) - v2np1k*(1 - xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4.)/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &               th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))) + (3*L1*(m1 + 3*m2)*(v1n + v1np1k - v1n*xi1 + v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 + xi1) - v2np1k*(1 + xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4. - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k - th2n*xi1 + th2np1k*xi1)/2.)
     & + 9*g*(m1 + 2*m2)*Sin(th1n + (th2n*(-1 + xi1))/2. - th1n*xi1 + 
     &                  th1np1k*(1 + xi1) - (th2np1k*(1 + xi1))/2.))/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)
     & ))))/18.))/12.))/
     &    (4.*(th1n - th1np1k)*((L1**2*
     &      (-2*g*(m1 + 2*m2)*Cos(th1n) + 2*g*(m1 + 2*m2)*Cos(th1np1k) +
     &         L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) +
     &            v1np1k*v2np1k*
     &             (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k))))**2)/
     &    (16.*(th1n - th1np1k)**2) + 
     &   (L2**2*m2**2*(2*L2*(v2n + v2np1k) + 3*L1*v1n*Cos(th1n - th2n) +
     &         3*L1*v1np1k*Cos(th1np1k - th2np1k))**2)/144. + 
     &   (L1**2*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &         3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &            v2n*Cos(th1np1k - th2np1k)))**2)/144. + 
     &   (L2**2*m2**2*(L1*v1np1k*v2np1k*
     &          (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &         L1*v1n*v2n*(Cos(th1np1k - th2n) - 
     &            Cos(th1np1k - th2np1k)) + 
     &         2*g*(-Cos(th2n) + Cos(th2np1k)))**2)/
     &    (16.*(th2n - th2np1k)**2))**2) + 
     & (L1*(-2*g*(m1 + 2*m2)*Cos(th1n) + 2*g*(m1 + 2*m2)*Cos(th1np1k) +
     &  L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &  v1np1k*v2np1k*(Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k))))*
     & ((L2*m2*(v2n + v2np1k)*(L1*v1np1k*v2np1k*
     &        (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &    L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k)) + 
     &    2*g*(-Cos(th2n) + Cos(th2np1k))))/(8.*(th2n - th2np1k)**2) + 
     &  (L1*L2*m2*v1np1k*(v1n + v1np1k)*v2np1k*
     &     (Sin(th1n - th2np1k) - Sin(th1np1k - th2np1k)))/
     &   (8.*(th1n - th1np1k)) + 
     &  (L2*m2*(v2n + v2np1k)*(-(L1*v1np1k*v2np1k*Sin(th1n - th2np1k))-
     &       L1*v1n*v2n*Sin(th1np1k - th2np1k) - 2*g*Sin(th2np1k)))/
     &   (8.*(th2n - th2np1k)) + 
     &  (L1*L2*m2*v2n*Sin(th1np1k - th2np1k)*
     &     ((-4*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k)/2.) + 
     &            m2*(3*g*Sin((th1n + th1np1k - 2*th2n - 2*th2np1k)/
     &                  2.) + 
     &               L2*(v2n + v2np1k)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &               (3*L1*(-v1n - v1np1k)**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k))/4.)))/
     &        (3.*L1*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((-3*(g*(4*m1 + 5*m2)*
     &                  Sin((th1n + th1np1k - th1n*xi1 + th1np1k*xi1)/
     &                   2.) + 
     &                 m2*(3*g*
     &                   Sin((th1n + th1np1k - 2*th2n - 2*th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + 2*th2n*xi1 - 
     &                   2*th2np1k*xi1)/2.) + 
     &                   L2*(v2n + v2np1k - v2n*xi1 + v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 + xi1) - v1np1k*(1 + xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4.)))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))) - (3*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k + th1n*xi1 - 
     &                   th1np1k*xi1)/2.) + 
     &                 m2*(L2*(v2n + v2np1k + v2n*xi1 - v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &                    (3*L1*(v1n*(-1 - xi1) - v1np1k*(1 - xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4. + 
     &                    3*g*
     &                    Sin((th1n + th1np1k - 2*th2n - 2*th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - 2*th2n*xi1 + 
     &                   2*th2np1k*xi1)/2.))))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                  th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1)
     & ))))/18.))/4. + (L1*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &       3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &          v2n*Cos(th1np1k - th2np1k)))*
     &     ((-4*m2*(-3*g*Cos((th1n + th1np1k - 2*th2n - 2*th2np1k)/2.)-
     &            (L2*(v2n + v2np1k)**2*
     &               Cos((th1n + th1np1k - th2n - th2np1k)/2.))/2. - 
     &            (3*L1*(-v1n - v1np1k)**2*
     &               Cos(th1n + th1np1k - th2n - th2np1k))/4.))/
     &        (3.*L1*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) - 
     &       (12*m2*Sin(th1n + th1np1k - th2n - th2np1k)*
     &          (g*(4*m1 + 5*m2)*Sin((th1n + th1np1k)/2.) + 
     &            m2*(3*g*Sin((th1n + th1np1k - 2*th2n - 2*th2np1k)/
     &                  2.) + 
     &               L2*(v2n + v2np1k)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &               (3*L1*(-v1n - v1np1k)**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k))/4.)))/
     &        (L1*(8*m1 + 15*m2 - 
     &             9*m2*Cos(th1n + th1np1k - th2n - th2np1k))**2) + 
     &       (5*((-3*m2*((3*g*(-2 - 2*xi1)*
     &                   Cos((th1n + th1np1k - 2*th2n - 2*th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + 2*th2n*xi1 - 
     &                  2*th2np1k*xi1)/2.))/2. + 
     &                 (L2*(-1 - xi1)*
     &                   (v2n + v2np1k - v2n*xi1 + v2np1k*xi1)**2*
     &                   Cos((th1n + th1np1k - th2n - th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                  th2np1k*xi1)/2.))/2. + 
     &                 (3*L1*(-1 - xi1)*
     &                   (v1n*(-1 + xi1) - v1np1k*(1 + xi1))**2*
     &                   Cos(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4.))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & )))-(3*m2*((L2*(-1+xi1)*(v2n+v2np1k + v2n*xi1 - v2np1k*xi1)**2*
     &                   Cos((th1n + th1np1k - th2n - th2np1k + 
     &                  th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                  th2np1k*xi1)/2.))/2. + 
     &                 (3*L1*(v1n*(-1 - xi1) - v1np1k*(1 - xi1))**2*
     &                   (-1 + xi1)*
     &                   Cos(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4. + 
     &                 (3*g*(-2 + 2*xi1)*
     &                   Cos((th1n + th1np1k - 2*th2n - 2*th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - 2*th2n*xi1 + 
     &                   2*th2np1k*xi1)/2.))/2.))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                 th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))) + (27*m2*(-1 - xi1)*Sin(th1n + th1np1k - th2n - th2np1k - 
     &                 th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)*
     &               (g*(4*m1 + 5*m2)*
     &                  Sin((th1n + th1np1k - th1n*xi1 + th1np1k*xi1)/
     &                   2.) + 
     &                 m2*(3*g*
     &                   Sin((th1n + th1np1k - 2*th2n - 2*th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + 2*th2n*xi1 - 
     &                   2*th2np1k*xi1)/2.) + 
     &                   L2*(v2n + v2np1k - v2n*xi1 + v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 + xi1) - v1np1k*(1 + xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4.)))/
     &             (L1*(8*m1 + 15*m2 - 
     &                  9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))**2) + (27*m2*(-1 + xi1)*Sin(th1n + th1np1k - th2n - th2np1k + 
     &                 th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1)*
     &               (g*(4*m1 + 5*m2)*
     &                  Sin((th1n + th1np1k + th1n*xi1 - th1np1k*xi1)/
     &                   2.) + 
     &                 m2*(L2*(v2n + v2np1k + v2n*xi1 - v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &                    (3*L1*(v1n*(-1 - xi1) - v1np1k*(1 - xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4. + 
     &                    3*g*
     &                    Sin((th1n + th1np1k - 2*th2n - 2*th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - 2*th2n*xi1 + 
     &                   2*th2np1k*xi1)/2.))))/
     &             (L1*(8*m1 + 15*m2 - 
     &                  9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                  th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))**2)))/18.))/12. + (L1*L2*m2*v1np1k*Sin(th1np1k - th2np1k)*
     &     ((4*(3*L1*(m1 + 3*m2)*(v1n + v1np1k)**2*
     &             Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &            (9*L2*m2*(-v2n - v2np1k)**2*
     &               Sin(th1n + th1np1k - th2n - th2np1k))/4. + 
     &            9*g*(m1 + 2*m2)*
     &             Sin(th1n + th1np1k - th2n/2. - th2np1k/2.) - 
     &            3*g*(m1 + 6*m2)*Sin((th2n + th2np1k)/2.)))/
     &        (9.*L2*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((9*g*(m1 + 2*m2)*
     &                Sin(th1n + (th2n*(-1 - xi1))/2. + 
     &                  th1np1k*(1 - xi1) - (th2np1k*(1 - xi1))/2. + 
     &                  th1n*xi1) - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k + th2n*xi1 - th2np1k*xi1)/2.)
     & + 3*L1*(m1 + 3*m2)*(v1n + v1np1k + v1n*xi1 - v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 - xi1) - v2np1k*(1 - xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4.)/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                  th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))) + (3*L1*(m1 + 3*m2)*(v1n + v1np1k - v1n*xi1 + v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 + xi1) - v2np1k*(1 + xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4. - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k - th2n*xi1 + th2np1k*xi1)/2.)
     & + 9*g*(m1 + 2*m2)*Sin(th1n + (th2n*(-1 + xi1))/2. - th1n*xi1 + 
     &                  th1np1k*(1 + xi1) - (th2np1k*(1 + xi1))/2.))/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                 th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)
     & ))))/18.))/4.+(L2*m2*(2*L2*(v2n+v2np1k)+3*L1*v1n*Cos(th1n-th2n)+
     &       3*L1*v1np1k*Cos(th1np1k - th2np1k))*
     &     ((4*((-3*L1*(m1 + 3*m2)*(v1n + v1np1k)**2*
     &               Cos((th1n + th1np1k - th2n - th2np1k)/2.))/2. - 
     &            (9*L2*m2*(-v2n - v2np1k)**2*
     &               Cos(th1n + th1np1k - th2n - th2np1k))/4. - 
     &            (9*g*(m1 + 2*m2)*
     &               Cos(th1n + th1np1k - th2n/2. - th2np1k/2.))/2. - 
     &            (3*g*(m1 + 6*m2)*Cos((th2n + th2np1k)/2.))/2.))/
     &        (9.*L2*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (4*m2*Sin(th1n + th1np1k - th2n - th2np1k)*
     &          (3*L1*(m1 + 3*m2)*(v1n + v1np1k)**2*
     &             Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &            (9*L2*m2*(-v2n - v2np1k)**2*
     &               Sin(th1n + th1np1k - th2n - th2np1k))/4. + 
     &            9*g*(m1 + 2*m2)*
     &             Sin(th1n + th1np1k - th2n/2. - th2np1k/2.) - 
     &            3*g*(m1 + 6*m2)*Sin((th2n + th2np1k)/2.)))/
     &        (L2*(8*m1 + 15*m2 - 
     &             9*m2*Cos(th1n + th1np1k - th2n - th2np1k))**2) + 
     &       (5*(((9*g*(m1 + 2*m2)*(-1 + xi1)*
     &                  Cos(th1n + (th2n*(-1 - xi1))/2. + 
     &                   th1np1k*(1 - xi1) - (th2np1k*(1 - xi1))/2. + 
     &                   th1n*xi1))/2. - 
     &               (3*g*(m1 + 6*m2)*(1 - xi1)*
     &                  Cos((th2n + th2np1k + th2n*xi1 - th2np1k*xi1)/
     &                   2.))/2. + 
     &               (3*L1*(m1 + 3*m2)*(-1 + xi1)*
     &                  (v1n + v1np1k + v1n*xi1 - v1np1k*xi1)**2*
     &                  Cos((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.))/2. + 
     &               (9*L2*m2*(v2n*(-1 - xi1) - v2np1k*(1 - xi1))**2*
     &                  (-1 + xi1)*
     &                  Cos(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))/4.)/(L2*(8*m1 + 15*m2 - 9*m2*
     &                  Cos(th1n + th1np1k - th2n - th2np1k + th1n*xi1 - 
     &                    th1np1k*xi1 - th2n*xi1 + th2np1k*xi1))) + 
     &            ((3*L1*(m1 + 3*m2)*(-1 - xi1)*
     &                  (v1n + v1np1k - v1n*xi1 + v1np1k*xi1)**2*
     &                  Cos((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.))/2. + 
     &               (9*L2*m2*(-1 - xi1)*
     &                  (v2n*(-1 + xi1) - v2np1k*(1 + xi1))**2*
     &                  Cos(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4. - 
     &               (3*g*(m1 + 6*m2)*(1 + xi1)*
     &                  Cos((th2n + th2np1k - th2n*xi1 + th2np1k*xi1)/
     &                   2.))/2. + 
     &               (9*g*(m1 + 2*m2)*(-1 - xi1)*
     &                  Cos(th1n + (th2n*(-1 + xi1))/2. - th1n*xi1 + 
     &              th1np1k*(1 + xi1) - (th2np1k*(1 + xi1))/2.))/2.)/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)
     & ))-(9*m2*(-1+xi1)*Sin(th1n+th1np1k-th2n - th2np1k + th1n*xi1 -
     &                 th1np1k*xi1 - th2n*xi1 + th2np1k*xi1)*
     &               (9*g*(m1 + 2*m2)*
     &                  Sin(th1n + (th2n*(-1 - xi1))/2. + 
     &                   th1np1k*(1 - xi1) - (th2np1k*(1 - xi1))/2. + 
     &                   th1n*xi1) - 
     &                 3*g*(m1 + 6*m2)*
     &                  Sin((th2n + th2np1k + th2n*xi1 - th2np1k*xi1)/
     &                   2.) + 
     &                 3*L1*(m1 + 3*m2)*
     &                  (v1n + v1np1k + v1n*xi1 - v1np1k*xi1)**2*
     &                  Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &                 (9*L2*m2*(v2n*(-1 - xi1) - v2np1k*(1 - xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))/4.))/(L2*(8*m1 + 15*m2 - 9*m2*
     &                   Cos(th1n + th1np1k - th2n - th2np1k + 
     &                th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))**2) - (9*m2*(-1 - xi1)*Sin(th1n + th1np1k - th2n - th2np1k - 
     &                 th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)*
     &               (3*L1*(m1 + 3*m2)*
     &                  (v1n + v1np1k - v1n*xi1 + v1np1k*xi1)**2*
     &                  Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & )/2.) + (9*L2*m2*(v2n*(-1 + xi1) - v2np1k*(1 + xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))/4.-3*g*(m1+6*m2)*Sin((th2n+th2np1k - th2n*xi1 + th2np1k*xi1)/
     &                   2.) + 
     &                 9*g*(m1 + 2*m2)*
     &                  Sin(th1n + (th2n*(-1 + xi1))/2. - th1n*xi1 + 
     &                    th1np1k*(1 + xi1) - (th2np1k*(1 + xi1))/2.)))/
     &             (L2*(8*m1 + 15*m2 - 
     &                  9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)
     & )**2)))/18.))/12.))/
     &    (4.*(th1n - th1np1k)*((L1**2*
     &     (-2*g*(m1 + 2*m2)*Cos(th1n) + 2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &        L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &           v1np1k*v2np1k*
     &            (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k))))**2)/
     &   (16.*(th1n - th1np1k)**2) + 
     &  (L2**2*m2**2*(2*L2*(v2n + v2np1k) + 3*L1*v1n*Cos(th1n - th2n) + 
     &        3*L1*v1np1k*Cos(th1np1k - th2np1k))**2)/144. + 
     &  (L1**2*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &        3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &           v2n*Cos(th1np1k - th2np1k)))**2)/144. + 
     &  (L2**2*m2**2*(L1*v1np1k*v2np1k*
     &         (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &   L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k)) + 
     &  2*g*(-Cos(th2n) + Cos(th2np1k)))**2)/(16.*(th2n - th2np1k)**2)))

      m13=-0.5+(L1*(-2*g*(m1+2*m2)*Cos(th1n)+2*g*(m1+2*m2)*Cos(th1np1k)+
     &  L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &   v1np1k*v2np1k*(Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k))))*
     & ((L1*L2*m2*v2np1k*(v2n + v2np1k)*
     &  (Cos(th1n - th2n) - Cos(th1n - th2np1k)))/(8.*(th2n - th2np1k))
     & +(L1*(-2*g*(m1+2*m2)*Cos(th1n)+2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &       L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &          v1np1k*v2np1k*(Cos(th1n - th2np1k) - 
     &             Cos(th1np1k - th2np1k)))))/(8.*(th1n - th1np1k)) + 
     &  (L1*L2*m2*(v1n + v1np1k)*v2np1k*
     &     (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k)))/
     &   (8.*(th1n - th1np1k)) + 
     &  (L2*m2*(2*L2*(v2n + v2np1k) + 3*L1*v1n*Cos(th1n - th2n) + 
     &       3*L1*v1np1k*Cos(th1np1k - th2np1k))*
     &     ((8*L1*(m1 + 3*m2)*(v1n + v1np1k)*
     &          Sin((th1n + th1np1k - th2n - th2np1k)/2.))/
     &        (3.*L2*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((6*L1*(m1 + 3*m2)*(1 + xi1)*
     &               (v1n + v1np1k - v1n*xi1 + v1np1k*xi1)*
     &               Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.))/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))) + (6*L1*(m1 + 3*m2)*(1 - xi1)*
     &               (v1n + v1np1k + v1n*xi1 - v1np1k*xi1)*
     &               Sin((th1n + th1np1k - th2n - th2np1k + th1n*xi1 - 
     &                   th1np1k*xi1 - th2n*xi1 + th2np1k*xi1)/2.))/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & )))))/18.))/12. + (L1*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &       3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &          v2n*Cos(th1np1k - th2np1k)))*
     &     ((2*m2*(-v1n - v1np1k)*Sin(th1n + th1np1k - th2n - th2np1k))/
     &        (8*m1 + 15*m2 - 9*m2*Cos(th1n + th1np1k - th2n - th2np1k))
     & + (5*((-9*m2*(-1 - xi1)*(v1n*(-1 + xi1) - v1np1k*(1 + xi1))*
     &               Sin(th1n + th1np1k - th2n - th2np1k - th1n*xi1 + 
     &                 th1np1k*xi1 + th2n*xi1 - th2np1k*xi1))/
     &             (2.*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))) - (9*m2*(v1n*(-1 - xi1) - v1np1k*(1 - xi1))*(-1 + xi1)*
     &               Sin(th1n + th1np1k - th2n - th2np1k + th1n*xi1 - 
     &                 th1np1k*xi1 - th2n*xi1 + th2np1k*xi1))/
     &             (2.*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                 th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & )))))/18.))/12. + (L1**2*(m1 + 3*m2)*
     &     ((-4*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k)/2.) + 
     &            m2*(3*g*Sin((th1n + th1np1k - 2*th2n - 2*th2np1k)/
     &                  2.) + 
     &               L2*(v2n + v2np1k)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &               (3*L1*(-v1n - v1np1k)**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k))/4.)))/
     &        (3.*L1*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((-3*(g*(4*m1 + 5*m2)*
     &                  Sin((th1n + th1np1k - th1n*xi1 + th1np1k*xi1)/
     &                   2.) + 
     &                 m2*(3*g*
     &                   Sin((th1n + th1np1k - 2*th2n - 2*th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + 2*th2n*xi1 - 
     &                  2*th2np1k*xi1)/2.) + 
     &                   L2*(v2n + v2np1k - v2n*xi1 + v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                  th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 + xi1) - v1np1k*(1 + xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4.)))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))) - (3*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k + th1n*xi1 - 
     &                   th1np1k*xi1)/2.) + 
     &                 m2*(L2*(v2n + v2np1k + v2n*xi1 - v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 - xi1) - v1np1k*(1 - xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4. + 
     &                   3*g*Sin((th1n + th1np1k - 2*th2n - 
     &                   2*th2np1k + th1n*xi1 - th1np1k*xi1 - 
     &                   2*th2n*xi1 + 2*th2np1k*xi1)/2.))))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                 th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & )))))/18.))/6. + (L1*L2*m2*Cos(th1np1k - th2np1k)*
     &     ((4*(3*L1*(m1 + 3*m2)*(v1n + v1np1k)**2*
     &             Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &            (9*L2*m2*(-v2n - v2np1k)**2*
     &               Sin(th1n + th1np1k - th2n - th2np1k))/4. + 
     &            9*g*(m1 + 2*m2)*
     &             Sin(th1n + th1np1k - th2n/2. - th2np1k/2.) - 
     &            3*g*(m1 + 6*m2)*Sin((th2n + th2np1k)/2.)))/
     &        (9.*L2*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((9*g*(m1 + 2*m2)*
     &                Sin(th1n + (th2n*(-1 - xi1))/2. + 
     &                  th1np1k*(1 - xi1) - (th2np1k*(1 - xi1))/2. + 
     &                  th1n*xi1) - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k + th2n*xi1 - th2np1k*xi1)/2.)
     & + 3*L1*(m1 + 3*m2)*(v1n + v1np1k + v1n*xi1 - v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 - xi1) - v2np1k*(1 - xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4.)/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                  th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))) + (3*L1*(m1 + 3*m2)*(v1n + v1np1k - v1n*xi1 + v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 + xi1) - v2np1k*(1 + xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4. - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k - th2n*xi1 + th2np1k*xi1)/2.)
     & + 9*g*(m1 + 2*m2)*Sin(th1n + (th2n*(-1 + xi1))/2. - th1n*xi1 + 
     &                  th1np1k*(1 + xi1) - (th2np1k*(1 + xi1))/2.))/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)
     & ))))/18.))/4.))/
     & (4.*(th1n - th1np1k)*((L1**2*
     &     (-2*g*(m1 + 2*m2)*Cos(th1n) + 2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &        L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &           v1np1k*v2np1k*
     &            (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k))))**2)/
     &   (16.*(th1n - th1np1k)**2) + 
     &  (L2**2*m2**2*(2*L2*(v2n + v2np1k) + 3*L1*v1n*Cos(th1n - th2n) + 
     &        3*L1*v1np1k*Cos(th1np1k - th2np1k))**2)/144. + 
     &  (L1**2*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &        3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &           v2n*Cos(th1np1k - th2np1k)))**2)/144. + 
     &  (L2**2*m2**2*(L1*v1np1k*v2np1k*
     &         (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &        L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k))+
     &   2*g*(-Cos(th2n) + Cos(th2np1k)))**2)/(16.*(th2n - th2np1k)**2))
     & ) - (L1*(-2*g*(m1 + 2*m2)*Cos(th1n)+2*g*(m1+2*m2)*Cos(th1np1k)+ 
     &  L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &  v1np1k*v2np1k*(Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k))))*
     & ((L1**2*L2*m2*v2np1k*(-2*g*(m1 + 2*m2)*Cos(th1n) + 
     &       2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &       L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &          v1np1k*v2np1k*
     &           (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k))))*
     &     (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k)))/
     &   (8.*(th1n - th1np1k)**2) + 
     &  (L1*L2**2*m2**2*Cos(th1np1k - th2np1k)*
     &     (2*L2*(v2n + v2np1k) + 3*L1*v1n*Cos(th1n - th2n) + 
     &       3*L1*v1np1k*Cos(th1np1k - th2np1k)))/24. + 
     &  (L1**3*(m1 + 3*m2)*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &       3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &          v2n*Cos(th1np1k - th2np1k))))/36. + 
     &  (L1*L2**2*m2**2*v2np1k*(Cos(th1n - th2n) - Cos(th1n - th2np1k))*
     &     (L1*v1np1k*v2np1k*(Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &       L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k)) + 
     &     2*g*(-Cos(th2n) + Cos(th2np1k))))/(8.*(th2n - th2np1k)**2))*
     & ((L1*(v1n + v1np1k)*(-2*g*(m1 + 2*m2)*Cos(th1n) + 
     &       2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &       L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &          v1np1k*v2np1k*(Cos(th1n - th2np1k) - 
     &             Cos(th1np1k - th2np1k)))))/(8.*(th1n - th1np1k)) + 
     &  (L2*m2*(v2n + v2np1k)*(L1*v1np1k*v2np1k*
     &        (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &    L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k)) + 
     &    2*g*(-Cos(th2n) + Cos(th2np1k))))/(8.*(th2n - th2np1k)) + 
     &  (L1*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &       3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &          v2n*Cos(th1np1k - th2np1k)))*
     &     ((-4*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k)/2.) + 
     &            m2*(3*g*Sin((th1n + th1np1k - 2*th2n - 2*th2np1k)/
     &                  2.) + 
     &               L2*(v2n + v2np1k)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &               (3*L1*(-v1n - v1np1k)**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k))/4.)))/
     &        (3.*L1*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((-3*(g*(4*m1 + 5*m2)*
     &                  Sin((th1n + th1np1k - th1n*xi1 + th1np1k*xi1)/
     &                   2.) + 
     &                 m2*(3*g*
     &                   Sin((th1n + th1np1k - 2*th2n - 2*th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + 2*th2n*xi1 - 
     &                  2*th2np1k*xi1)/2.) + 
     &                   L2*(v2n + v2np1k - v2n*xi1 + v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                  th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 + xi1) - v1np1k*(1 + xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4.)))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))) - (3*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k + th1n*xi1 - 
     &                   th1np1k*xi1)/2.) + 
     &                 m2*(L2*(v2n + v2np1k + v2n*xi1 - v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 - xi1) - v1np1k*(1 - xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4. + 
     &                   3*g*Sin((th1n + th1np1k - 2*th2n - 
     &                   2*th2np1k + th1n*xi1 - th1np1k*xi1 - 
     &                   2*th2n*xi1 + 2*th2np1k*xi1)/2.))))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & )))))/18.))/12. + (L2*m2*(2*L2*(v2n + v2np1k) + 
     &       3*L1*v1n*Cos(th1n - th2n) + 
     &       3*L1*v1np1k*Cos(th1np1k - th2np1k))*
     &     ((4*(3*L1*(m1 + 3*m2)*(v1n + v1np1k)**2*
     &             Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &            (9*L2*m2*(-v2n - v2np1k)**2*
     &               Sin(th1n + th1np1k - th2n - th2np1k))/4. + 
     &            9*g*(m1 + 2*m2)*
     &             Sin(th1n + th1np1k - th2n/2. - th2np1k/2.) - 
     &            3*g*(m1 + 6*m2)*Sin((th2n + th2np1k)/2.)))/
     &        (9.*L2*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((9*g*(m1 + 2*m2)*
     &                Sin(th1n + (th2n*(-1 - xi1))/2. + 
     &                  th1np1k*(1 - xi1) - (th2np1k*(1 - xi1))/2. + 
     &                  th1n*xi1) - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k + th2n*xi1 - th2np1k*xi1)/2.)
     & + 3*L1*(m1 + 3*m2)*(v1n + v1np1k + v1n*xi1 - v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 - xi1) - v2np1k*(1 - xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4.)/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                 th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))) + (3*L1*(m1 + 3*m2)*(v1n + v1np1k - v1n*xi1 + v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 + xi1) - v2np1k*(1 + xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4. - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k - th2n*xi1 + th2np1k*xi1)/2.)
     & + 9*g*(m1 + 2*m2)*Sin(th1n + (th2n*(-1 + xi1))/2. - th1n*xi1 + 
     &                  th1np1k*(1 + xi1) - (th2np1k*(1 + xi1))/2.))/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)
     & ))))/18.))/12.))/
     &    (4.*(th1n - th1np1k)*((L1**2*
     &      (-2*g*(m1 + 2*m2)*Cos(th1n) + 2*g*(m1 + 2*m2)*Cos(th1np1k) +
     &         L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) +
     &            v1np1k*v2np1k*
     &             (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k))))**2)/
     &    (16.*(th1n - th1np1k)**2) + 
     &   (L2**2*m2**2*(2*L2*(v2n + v2np1k) + 3*L1*v1n*Cos(th1n - th2n) +
     &         3*L1*v1np1k*Cos(th1np1k - th2np1k))**2)/144. + 
     &   (L1**2*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &         3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &            v2n*Cos(th1np1k - th2np1k)))**2)/144. + 
     &   (L2**2*m2**2*(L1*v1np1k*v2np1k*
     &          (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &         L1*v1n*v2n*(Cos(th1np1k - th2n) - 
     &            Cos(th1np1k - th2np1k)) + 
     &         2*g*(-Cos(th2n) + Cos(th2np1k)))**2)/
     &    (16.*(th2n - th2np1k)**2))**2) + 
     &  (L1*L2*m2*v2np1k*(Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k))*
     & ((L1*(v1n + v1np1k)*(-2*g*(m1 + 2*m2)*Cos(th1n) + 
     &       2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &       L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &          v1np1k*v2np1k*(Cos(th1n - th2np1k) - 
     &             Cos(th1np1k - th2np1k)))))/(8.*(th1n - th1np1k)) + 
     &  (L2*m2*(v2n + v2np1k)*(L1*v1np1k*v2np1k*
     &        (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &       L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k)) +
     &       2*g*(-Cos(th2n) + Cos(th2np1k))))/(8.*(th2n - th2np1k)) + 
     &  (L1*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &       3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &          v2n*Cos(th1np1k - th2np1k)))*
     &     ((-4*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k)/2.) + 
     &            m2*(3*g*Sin((th1n + th1np1k - 2*th2n - 2*th2np1k)/
     &                  2.) + 
     &               L2*(v2n + v2np1k)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &               (3*L1*(-v1n - v1np1k)**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k))/4.)))/
     &        (3.*L1*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((-3*(g*(4*m1 + 5*m2)*
     &                  Sin((th1n + th1np1k - th1n*xi1 + th1np1k*xi1)/
     &                   2.) + 
     &                 m2*(3*g*
     &                   Sin((th1n + th1np1k - 2*th2n - 2*th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + 2*th2n*xi1 - 
     &                   2*th2np1k*xi1)/2.) + 
     &                   L2*(v2n + v2np1k - v2n*xi1 + v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 + xi1) - v1np1k*(1 + xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4.)))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                 th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))) - (3*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k + th1n*xi1 - 
     &                   th1np1k*xi1)/2.) + 
     &                 m2*(L2*(v2n + v2np1k + v2n*xi1 - v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &                    (3*L1*(v1n*(-1 - xi1) - v1np1k*(1 - xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4. + 
     &                    3*g*
     &                    Sin((th1n + th1np1k - 2*th2n - 2*th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - 2*th2n*xi1 + 
     &                   2*th2np1k*xi1)/2.))))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1)
     & ))))/18.))/12. + (L2*m2*(2*L2*(v2n + v2np1k) + 
     &   3*L1*v1n*Cos(th1n - th2n) + 3*L1*v1np1k*Cos(th1np1k - th2np1k)
     & )*((4*(3*L1*(m1 + 3*m2)*(v1n + v1np1k)**2*
     &             Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &            (9*L2*m2*(-v2n - v2np1k)**2*
     &               Sin(th1n + th1np1k - th2n - th2np1k))/4. + 
     &            9*g*(m1 + 2*m2)*
     &             Sin(th1n + th1np1k - th2n/2. - th2np1k/2.) - 
     &            3*g*(m1 + 6*m2)*Sin((th2n + th2np1k)/2.)))/
     &        (9.*L2*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((9*g*(m1 + 2*m2)*
     &                Sin(th1n + (th2n*(-1 - xi1))/2. + 
     &                  th1np1k*(1 - xi1) - (th2np1k*(1 - xi1))/2. + 
     &                  th1n*xi1) - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k + th2n*xi1 - th2np1k*xi1)/2.)
     & + 3*L1*(m1 + 3*m2)*(v1n + v1np1k + v1n*xi1 - v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 - xi1) - v2np1k*(1 - xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))/4.)/(L2*(8*m1 + 15*m2 - 9*m2*
     &                  Cos(th1n + th1np1k - th2n - th2np1k + th1n*xi1 - 
     &                    th1np1k*xi1 - th2n*xi1 + th2np1k*xi1))) + 
     &            (3*L1*(m1 + 3*m2)*
     &                (v1n + v1np1k - v1n*xi1 + v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k - th1n*xi1 +
     &                   th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 + xi1) - v2np1k*(1 + xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))/4.-3*g*(m1+6*m2)*Sin((th2n+th2np1k - th2n*xi1 + th2np1k*xi1)/
     &                  2.) + 9*g*(m1 + 2*m2)*
     &                Sin(th1n + (th2n*(-1 + xi1))/2. - th1n*xi1 + 
     &                  th1np1k*(1 + xi1) - (th2np1k*(1 + xi1))/2.))/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &               th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1))
     & )))/18.))/12.))/
     &    (4.*(th1n - th1np1k)*((L1**2*
     &     (-2*g*(m1 + 2*m2)*Cos(th1n) + 2*g*(m1 + 2*m2)*Cos(th1np1k) +
     &        L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) +
     &           v1np1k*v2np1k*
     &            (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k))))**2)/
     &   (16.*(th1n - th1np1k)**2) + 
     &  (L2**2*m2**2*(2*L2*(v2n + v2np1k) + 3*L1*v1n*Cos(th1n - th2n) +
     &        3*L1*v1np1k*Cos(th1np1k - th2np1k))**2)/144. + 
     &  (L1**2*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &        3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &           v2n*Cos(th1np1k - th2np1k)))**2)/144. + 
     &  (L2**2*m2**2*(L1*v1np1k*v2np1k*
     &         (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     & L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k)) + 
     & 2*g*(-Cos(th2n) + Cos(th2np1k)))**2)/(16.*(th2n - th2np1k)**2)))

      m14=(L1*(-2*g*(m1 + 2*m2)*Cos(th1n)+2*g*(m1+2*m2)*Cos(th1np1k)+
     &  L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &   v1np1k*v2np1k*(Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k))))*
     & ((L1*L2*m2*v1np1k*(v2n + v2np1k)*
     &   (Cos(th1n - th2n) - Cos(th1n - th2np1k)))/(8.*(th2n - th2np1k))
     & + (L1*L2*m2*v1np1k*(v1n + v1np1k)*
     &     (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k)))/
     &   (8.*(th1n - th1np1k)) + 
     &  (L2*m2*(L1*v1np1k*v2np1k*
     &        (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &       L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k)) +
     &       2*g*(-Cos(th2n) + Cos(th2np1k))))/(8.*(th2n - th2np1k)) + 
     &  (L1*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &       3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &          v2n*Cos(th1np1k - th2np1k)))*
     &     ((-8*L2*m2*(v2n + v2np1k)*
     &          Sin((th1n + th1np1k - th2n - th2np1k)/2.))/
     &        (3.*L1*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((-6*L2*m2*(1 + xi1)*
     &               (v2n + v2np1k - v2n*xi1 + v2np1k*xi1)*
     &               Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))) - (6*L2*m2*(1 - xi1)*(v2n + v2np1k + v2n*xi1 - v2np1k*xi1)*
     &               Sin((th1n + th1np1k - th2n - th2np1k + th1n*xi1 - 
     &                   th1np1k*xi1 - th2n*xi1 + th2np1k*xi1)/2.))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                  th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & )))))/18.))/12. + (L2*m2*(2*L2*(v2n + v2np1k) + 
     &       3*L1*v1n*Cos(th1n - th2n) + 
     &       3*L1*v1np1k*Cos(th1np1k - th2np1k))*
     &   ((-2*m2*(-v2n - v2np1k)*Sin(th1n + th1np1k - th2n - th2np1k))/
     &        (8*m1 + 15*m2 - 9*m2*Cos(th1n + th1np1k - th2n - th2np1k))
     & + (5*((9*m2*(-1 - xi1)*(v2n*(-1 + xi1) - v2np1k*(1 + xi1))*
     &               Sin(th1n + th1np1k - th2n - th2np1k - th1n*xi1 + 
     &                 th1np1k*xi1 + th2n*xi1 - th2np1k*xi1))/
     &             (2.*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))) + (9*m2*(v2n*(-1 - xi1) - v2np1k*(1 - xi1))*(-1 + xi1)*
     &               Sin(th1n + th1np1k - th2n - th2np1k + th1n*xi1 - 
     &                 th1np1k*xi1 - th2n*xi1 + th2np1k*xi1))/
     &             (2.*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                  th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & )))))/18.))/12. + (L1*L2*m2*Cos(th1n - th2n)*
     &     ((-4*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k)/2.) + 
     &            m2*(3*g*Sin((th1n + th1np1k - 2*th2n - 2*th2np1k)/
     &                  2.) + 
     &               L2*(v2n + v2np1k)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &               (3*L1*(-v1n - v1np1k)**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k))/4.)))/
     &        (3.*L1*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((-3*(g*(4*m1 + 5*m2)*
     &                  Sin((th1n + th1np1k - th1n*xi1 + th1np1k*xi1)/
     &                   2.) + 
     &                 m2*(3*g*
     &                   Sin((th1n + th1np1k - 2*th2n - 2*th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + 2*th2n*xi1 - 
     &                  2*th2np1k*xi1)/2.) + 
     &                   L2*(v2n + v2np1k - v2n*xi1 + v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                  th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 + xi1) - v1np1k*(1 + xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4.)))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                 th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))) - (3*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k + th1n*xi1 - 
     &                   th1np1k*xi1)/2.) + 
     &                 m2*(L2*(v2n + v2np1k + v2n*xi1 - v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 - xi1) - v1np1k*(1 - xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4. + 
     &                   3*g*Sin((th1n + th1np1k - 2*th2n - 
     &                   2*th2np1k + th1n*xi1 - th1np1k*xi1 - 
     &                   2*th2n*xi1 + 2*th2np1k*xi1)/2.))))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                  th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & )))))/18.))/4. + (L2**2*m2*((4*
     &          (3*L1*(m1 + 3*m2)*(v1n + v1np1k)**2*
     &             Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &            (9*L2*m2*(-v2n - v2np1k)**2*
     &               Sin(th1n + th1np1k - th2n - th2np1k))/4. + 
     &            9*g*(m1 + 2*m2)*
     &             Sin(th1n + th1np1k - th2n/2. - th2np1k/2.) - 
     &            3*g*(m1 + 6*m2)*Sin((th2n + th2np1k)/2.)))/
     &        (9.*L2*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((9*g*(m1 + 2*m2)*
     &                Sin(th1n + (th2n*(-1 - xi1))/2. + 
     &                  th1np1k*(1 - xi1) - (th2np1k*(1 - xi1))/2. + 
     &                  th1n*xi1) - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k + th2n*xi1 - th2np1k*xi1)/2.)
     & + 3*L1*(m1 + 3*m2)*(v1n + v1np1k + v1n*xi1 - v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 - xi1) - v2np1k*(1 - xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4.)/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                  th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))) + (3*L1*(m1 + 3*m2)*(v1n + v1np1k - v1n*xi1 + v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 + xi1) - v2np1k*(1 + xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4. - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k - th2n*xi1 + th2np1k*xi1)/2.)
     & + 9*g*(m1 + 2*m2)*Sin(th1n + (th2n*(-1 + xi1))/2. - th1n*xi1 + 
     &                  th1np1k*(1 + xi1) - (th2np1k*(1 + xi1))/2.))/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                 th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)
     & ))))/18.))/6.))/
     & (4.*(th1n - th1np1k)*((L1**2*
     &     (-2*g*(m1 + 2*m2)*Cos(th1n) + 2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &        L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &           v1np1k*v2np1k*
     &            (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k))))**2)/
     &   (16.*(th1n - th1np1k)**2) + 
     &  (L2**2*m2**2*(2*L2*(v2n + v2np1k) + 3*L1*v1n*Cos(th1n - th2n) + 
     &        3*L1*v1np1k*Cos(th1np1k - th2np1k))**2)/144. + 
     &  (L1**2*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &        3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &           v2n*Cos(th1np1k - th2np1k)))**2)/144. + 
     &  (L2**2*m2**2*(L1*v1np1k*v2np1k*
     &         (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &     L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k)) + 
     & 2*g*(-Cos(th2n) + Cos(th2np1k)))**2)/(16.*(th2n - th2np1k)**2))
     & ) - (L1*(-2*g*(m1 + 2*m2)*Cos(th1n)+2*g*(m1+2*m2)*Cos(th1np1k)+ 
     &  L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &   v1np1k*v2np1k*(Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k))))*
     & ((L1**2*L2*m2*v1np1k*(-2*g*(m1 + 2*m2)*Cos(th1n) + 
     &       2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &       L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &          v1np1k*v2np1k*
     &           (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k))))*
     &     (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k)))/
     &   (8.*(th1n - th1np1k)**2) + 
     &  (L2**3*m2**2*(2*L2*(v2n + v2np1k) + 3*L1*v1n*Cos(th1n - th2n) +
     &       3*L1*v1np1k*Cos(th1np1k - th2np1k)))/36. + 
     &  (L1**2*L2*m2*Cos(th1n - th2n)*
     &     (2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &       3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &          v2n*Cos(th1np1k - th2np1k))))/24. + 
     &  (L1*L2**2*m2**2*v1np1k*(Cos(th1n - th2n) - Cos(th1n - th2np1k))*
     &     (L1*v1np1k*v2np1k*(Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &       L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k)) + 
     &     2*g*(-Cos(th2n) + Cos(th2np1k))))/(8.*(th2n - th2np1k)**2))*
     & ((L1*(v1n + v1np1k)*(-2*g*(m1 + 2*m2)*Cos(th1n) + 
     &       2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &       L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &          v1np1k*v2np1k*(Cos(th1n - th2np1k) - 
     &             Cos(th1np1k - th2np1k)))))/(8.*(th1n - th1np1k)) + 
     &  (L2*m2*(v2n + v2np1k)*(L1*v1np1k*v2np1k*
     &        (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &       L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k)) +
     &       2*g*(-Cos(th2n) + Cos(th2np1k))))/(8.*(th2n - th2np1k)) + 
     &  (L1*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &       3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &          v2n*Cos(th1np1k - th2np1k)))*
     &     ((-4*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k)/2.) + 
     &            m2*(3*g*Sin((th1n + th1np1k - 2*th2n - 2*th2np1k)/
     &                  2.) + 
     &               L2*(v2n + v2np1k)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &               (3*L1*(-v1n - v1np1k)**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k))/4.)))/
     &        (3.*L1*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((-3*(g*(4*m1 + 5*m2)*
     &                  Sin((th1n + th1np1k - th1n*xi1 + th1np1k*xi1)/
     &                   2.) + 
     &                 m2*(3*g*
     &                   Sin((th1n + th1np1k - 2*th2n - 2*th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + 2*th2n*xi1 - 
     &                  2*th2np1k*xi1)/2.) + 
     &                   L2*(v2n + v2np1k - v2n*xi1 + v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                  th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 + xi1) - v1np1k*(1 + xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4.)))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))) - (3*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k + th1n*xi1 - 
     &                   th1np1k*xi1)/2.) + 
     &                 m2*(L2*(v2n + v2np1k + v2n*xi1 - v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 - xi1) - v1np1k*(1 - xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4. + 
     &                   3*g*Sin((th1n + th1np1k - 2*th2n - 
     &                   2*th2np1k + th1n*xi1 - th1np1k*xi1 - 
     &                   2*th2n*xi1 + 2*th2np1k*xi1)/2.))))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & )))))/18.))/12. + (L2*m2*(2*L2*(v2n + v2np1k) + 
     &       3*L1*v1n*Cos(th1n - th2n) + 
     &       3*L1*v1np1k*Cos(th1np1k - th2np1k))*
     &     ((4*(3*L1*(m1 + 3*m2)*(v1n + v1np1k)**2*
     &             Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &            (9*L2*m2*(-v2n - v2np1k)**2*
     &               Sin(th1n + th1np1k - th2n - th2np1k))/4. + 
     &            9*g*(m1 + 2*m2)*
     &             Sin(th1n + th1np1k - th2n/2. - th2np1k/2.) - 
     &            3*g*(m1 + 6*m2)*Sin((th2n + th2np1k)/2.)))/
     &        (9.*L2*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((9*g*(m1 + 2*m2)*
     &                Sin(th1n + (th2n*(-1 - xi1))/2. + 
     &                  th1np1k*(1 - xi1) - (th2np1k*(1 - xi1))/2. + 
     &                  th1n*xi1) - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k + th2n*xi1 - th2np1k*xi1)/2.)
     & + 3*L1*(m1 + 3*m2)*(v1n + v1np1k + v1n*xi1 - v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 - xi1) - v2np1k*(1 - xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4.)/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))) + (3*L1*(m1 + 3*m2)*(v1n + v1np1k - v1n*xi1 + v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 + xi1) - v2np1k*(1 + xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4. - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k - th2n*xi1 + th2np1k*xi1)/2.)
     & + 9*g*(m1 + 2*m2)*Sin(th1n + (th2n*(-1 + xi1))/2. - th1n*xi1 + 
     &                  th1np1k*(1 + xi1) - (th2np1k*(1 + xi1))/2.))/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                 th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)
     & ))))/18.))/12.))/
     & (4.*(th1n - th1np1k)*((L1**2*
     &      (-2*g*(m1 + 2*m2)*Cos(th1n) + 2*g*(m1 + 2*m2)*Cos(th1np1k) +
     &         L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) +
     &            v1np1k*v2np1k*
     &             (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k))))**2)/
     &    (16.*(th1n - th1np1k)**2) + 
     &   (L2**2*m2**2*(2*L2*(v2n + v2np1k) + 3*L1*v1n*Cos(th1n - th2n) +
     &         3*L1*v1np1k*Cos(th1np1k - th2np1k))**2)/144. + 
     &   (L1**2*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &         3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &            v2n*Cos(th1np1k - th2np1k)))**2)/144. + 
     &   (L2**2*m2**2*(L1*v1np1k*v2np1k*
     &          (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &         L1*v1n*v2n*(Cos(th1np1k - th2n) - 
     &            Cos(th1np1k - th2np1k)) + 
     &         2*g*(-Cos(th2n) + Cos(th2np1k)))**2)/
     &    (16.*(th2n - th2np1k)**2))**2) + 
     & (L1*L2*m2*v1np1k*(Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k))*
     & ((L1*(v1n + v1np1k)*(-2*g*(m1 + 2*m2)*Cos(th1n) + 
     &       2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &       L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &          v1np1k*v2np1k*(Cos(th1n - th2np1k) - 
     &             Cos(th1np1k - th2np1k)))))/(8.*(th1n - th1np1k)) + 
     &  (L2*m2*(v2n + v2np1k)*(L1*v1np1k*v2np1k*
     &        (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &       L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k)) +
     &       2*g*(-Cos(th2n) + Cos(th2np1k))))/(8.*(th2n - th2np1k)) + 
     &  (L1*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &       3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &          v2n*Cos(th1np1k - th2np1k)))*
     &     ((-4*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k)/2.) + 
     &            m2*(3*g*Sin((th1n + th1np1k - 2*th2n - 2*th2np1k)/
     &                  2.) + 
     &               L2*(v2n + v2np1k)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &               (3*L1*(-v1n - v1np1k)**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k))/4.)))/
     &        (3.*L1*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((-3*(g*(4*m1 + 5*m2)*
     &                  Sin((th1n + th1np1k - th1n*xi1 + th1np1k*xi1)/
     &                   2.) + 
     &                 m2*(3*g*
     &                   Sin((th1n + th1np1k - 2*th2n - 2*th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + 2*th2n*xi1 - 
     &                   2*th2np1k*xi1)/2.) + 
     &                   L2*(v2n + v2np1k - v2n*xi1 + v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 + xi1) - v1np1k*(1 + xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4.)))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                 th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))) - (3*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k + th1n*xi1 - 
     &                   th1np1k*xi1)/2.) + 
     &                 m2*(L2*(v2n + v2np1k + v2n*xi1 - v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &                    (3*L1*(v1n*(-1 - xi1) - v1np1k*(1 - xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4. + 
     &                    3*g*
     &                    Sin((th1n + th1np1k - 2*th2n - 2*th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - 2*th2n*xi1 + 
     &                   2*th2np1k*xi1)/2.))))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &         th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1)
     & ))))/18.))/12. + (L2*m2*(2*L2*(v2n + v2np1k) + 
     &  3*L1*v1n*Cos(th1n - th2n) + 3*L1*v1np1k*Cos(th1np1k - th2np1k)
     & )*((4*(3*L1*(m1 + 3*m2)*(v1n + v1np1k)**2*
     &             Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &            (9*L2*m2*(-v2n - v2np1k)**2*
     &               Sin(th1n + th1np1k - th2n - th2np1k))/4. + 
     &            9*g*(m1 + 2*m2)*
     &             Sin(th1n + th1np1k - th2n/2. - th2np1k/2.) - 
     &            3*g*(m1 + 6*m2)*Sin((th2n + th2np1k)/2.)))/
     &        (9.*L2*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((9*g*(m1 + 2*m2)*
     &                Sin(th1n + (th2n*(-1 - xi1))/2. + 
     &                  th1np1k*(1 - xi1) - (th2np1k*(1 - xi1))/2. + 
     &                  th1n*xi1) - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k + th2n*xi1 - th2np1k*xi1)/2.)
     & + 3*L1*(m1 + 3*m2)*(v1n + v1np1k + v1n*xi1 - v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 - xi1) - v2np1k*(1 - xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))/4.)/(L2*(8*m1 + 15*m2 - 9*m2*
     &                  Cos(th1n + th1np1k - th2n - th2np1k + th1n*xi1 -
     &                    th1np1k*xi1 - th2n*xi1 + th2np1k*xi1))) + 
     &            (3*L1*(m1 + 3*m2)*
     &                (v1n + v1np1k - v1n*xi1 + v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k - th1n*xi1 + 
     &                   th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 + xi1) - v2np1k*(1 + xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))/4. - 3*g*(m1 + 6*m2)*Sin((th2n+th2np1k-th2n*xi1+th2np1k*xi1)/
     &                  2.) + 9*g*(m1 + 2*m2)*
     &                Sin(th1n + (th2n*(-1 + xi1))/2. - th1n*xi1 + 
     &                  th1np1k*(1 + xi1) - (th2np1k*(1 + xi1))/2.))/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &             th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1))
     & )))/18.))/12.))/
     & (4.*(th1n - th1np1k)*((L1**2*
     &     (-2*g*(m1 + 2*m2)*Cos(th1n) + 2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &        L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &           v1np1k*v2np1k*
     &            (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k))))**2)/
     &   (16.*(th1n - th1np1k)**2) + 
     &  (L2**2*m2**2*(2*L2*(v2n + v2np1k) + 3*L1*v1n*Cos(th1n - th2n) + 
     &        3*L1*v1np1k*Cos(th1np1k - th2np1k))**2)/144. + 
     &  (L1**2*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &        3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &           v2n*Cos(th1np1k - th2np1k)))**2)/144. + 
     &  (L2**2*m2**2*(L1*v1np1k*v2np1k*
     &         (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &        L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k))+
     &  2*g*(-Cos(th2n) + Cos(th2np1k)))**2)/(16.*(th2n - th2np1k)**2)))
 
      m21=(L1*L2*m2*v1n*v2n*(-Sin(th1np1k-th2n)+Sin(th1np1k - th2np1k))*
     & ((L1*(v1n + v1np1k)*(-2*g*(m1 + 2*m2)*Cos(th1n) + 
     &       2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &       L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &          v1np1k*v2np1k*(Cos(th1n - th2np1k) - 
     &             Cos(th1np1k - th2np1k)))))/(8.*(th1n - th1np1k)) + 
     &  (L2*m2*(v2n + v2np1k)*(L1*v1np1k*v2np1k*
     &        (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &       L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k)) +
     &       2*g*(-Cos(th2n) + Cos(th2np1k))))/(8.*(th2n - th2np1k)) + 
     &  (L1*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &       3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &          v2n*Cos(th1np1k - th2np1k)))*
     &     ((-4*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k)/2.) + 
     &            m2*(3*g*Sin((th1n + th1np1k - 2*th2n - 2*th2np1k)/
     &                  2.) + 
     &               L2*(v2n + v2np1k)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &               (3*L1*(-v1n - v1np1k)**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k))/4.)))/
     &        (3.*L1*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((-3*(g*(4*m1 + 5*m2)*
     &                  Sin((th1n + th1np1k - th1n*xi1 + th1np1k*xi1)/
     &                   2.) + 
     &                 m2*(3*g*
     &                   Sin((th1n + th1np1k - 2*th2n - 2*th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + 2*th2n*xi1 - 
     &                  2*th2np1k*xi1)/2.) + 
     &                   L2*(v2n + v2np1k - v2n*xi1 + v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                  th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 + xi1) - v1np1k*(1 + xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4.)))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))) - (3*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k + th1n*xi1 - 
     &                   th1np1k*xi1)/2.) + 
     &                 m2*(L2*(v2n + v2np1k + v2n*xi1 - v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 - xi1) - v1np1k*(1 - xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4. + 
     &                   3*g*Sin((th1n + th1np1k - 2*th2n - 
     &                   2*th2np1k + th1n*xi1 - th1np1k*xi1 - 
     &                   2*th2n*xi1 + 2*th2np1k*xi1)/2.))))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                  th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & )))))/18.))/12. + (L2*m2*(2*L2*(v2n + v2np1k) + 
     &       3*L1*v1n*Cos(th1n - th2n) + 
     &       3*L1*v1np1k*Cos(th1np1k - th2np1k))*
     &     ((4*(3*L1*(m1 + 3*m2)*(v1n + v1np1k)**2*
     &             Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &            (9*L2*m2*(-v2n - v2np1k)**2*
     &               Sin(th1n + th1np1k - th2n - th2np1k))/4. + 
     &            9*g*(m1 + 2*m2)*
     &             Sin(th1n + th1np1k - th2n/2. - th2np1k/2.) - 
     &            3*g*(m1 + 6*m2)*Sin((th2n + th2np1k)/2.)))/
     &        (9.*L2*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((9*g*(m1 + 2*m2)*
     &                Sin(th1n + (th2n*(-1 - xi1))/2. + 
     &                  th1np1k*(1 - xi1) - (th2np1k*(1 - xi1))/2. + 
     &                  th1n*xi1) - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k + th2n*xi1 - th2np1k*xi1)/2.)
     & + 3*L1*(m1 + 3*m2)*(v1n + v1np1k + v1n*xi1 - v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 - xi1) - v2np1k*(1 - xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4.)/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                  th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))) + (3*L1*(m1 + 3*m2)*(v1n + v1np1k - v1n*xi1 + v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 + xi1) - v2np1k*(1 + xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4. - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k - th2n*xi1 + th2np1k*xi1)/2.)
     & + 9*g*(m1 + 2*m2)*Sin(th1n + (th2n*(-1 + xi1))/2. - th1n*xi1 + 
     &                  th1np1k*(1 + xi1) - (th2np1k*(1 + xi1))/2.))/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                 th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)
     & ))))/18.))/12.))/
     &    (4.*(th2n - th2np1k)*((L1**2*
     &     (-2*g*(m1 + 2*m2)*Cos(th1n) + 2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &        L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &           v1np1k*v2np1k*
     &            (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k))))**2)/
     &   (16.*(th1n - th1np1k)**2) + 
     &  (L2**2*m2**2*(2*L2*(v2n + v2np1k) + 3*L1*v1n*Cos(th1n - th2n) + 
     &        3*L1*v1np1k*Cos(th1np1k - th2np1k))**2)/144. + 
     &  (L1**2*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &        3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &           v2n*Cos(th1np1k - th2np1k)))**2)/144. + 
     &  (L2**2*m2**2*(L1*v1np1k*v2np1k*
     &         (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &      L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k))+
     &  2*g*(-Cos(th2n) + Cos(th2np1k)))**2)/(16.*(th2n - th2np1k)**2))
     & )-(L2*m2*(L1*v1np1k*v2np1k*(Cos(th1n-th2n)-Cos(th1n - th2np1k))+
     &  L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k)) + 
     &  2*g*(-Cos(th2n) + Cos(th2np1k)))*
     & ((L1**2*(-2*g*(m1 + 2*m2)*Cos(th1n) + 
     &        2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &        L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) +
     &           v1np1k*v2np1k*
     &            (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k))))**2)/
     &   (8.*(th1n - th1np1k)**3) - 
     &  (L1*L2**2*m2**2*v1np1k*
     &     (2*L2*(v2n + v2np1k) + 3*L1*v1n*Cos(th1n - th2n) + 
     &     3*L1*v1np1k*Cos(th1np1k - th2np1k))*Sin(th1np1k - th2np1k))/
     &   24. - (L1**2*L2*m2*v2n*
     &     (2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &       3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &      v2n*Cos(th1np1k - th2np1k)))*Sin(th1np1k - th2np1k))/24. + 
     &  (L1*L2**2*m2**2*v1n*v2n*
     &     (L1*v1np1k*v2np1k*(Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &       L1*v1n*v2n*(Cos(th1np1k - th2n) - 
     &          Cos(th1np1k - th2np1k)) + 
     &       2*g*(-Cos(th2n) + Cos(th2np1k)))*
     &     (-Sin(th1np1k - th2n) + Sin(th1np1k - th2np1k)))/
     &   (8.*(th2n - th2np1k)**2) + 
     &  (L1**2*(-2*g*(m1 + 2*m2)*Cos(th1n) + 
     &       2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &       L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &          v1np1k*v2np1k*
     &           (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k))))*
     &     (-2*g*(m1 + 2*m2)*Sin(th1np1k) + 
     &       L2*m2*(v1n*v2n*Sin(th1np1k - th2n) + 
     &          v1np1k*v2np1k*Sin(th1np1k - th2np1k))))/
     &   (8.*(th1n - th1np1k)**2))*
     & ((L1*(v1n + v1np1k)*(-2*g*(m1 + 2*m2)*Cos(th1n) + 
     &       2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &       L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &          v1np1k*v2np1k*(Cos(th1n - th2np1k) - 
     &             Cos(th1np1k - th2np1k)))))/(8.*(th1n - th1np1k)) + 
     &  (L2*m2*(v2n + v2np1k)*(L1*v1np1k*v2np1k*
     &        (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &       L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k)) +
     &       2*g*(-Cos(th2n) + Cos(th2np1k))))/(8.*(th2n - th2np1k)) + 
     &  (L1*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &       3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &          v2n*Cos(th1np1k - th2np1k)))*
     &     ((-4*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k)/2.) + 
     &            m2*(3*g*Sin((th1n + th1np1k - 2*th2n - 2*th2np1k)/
     &                  2.) + 
     &               L2*(v2n + v2np1k)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &               (3*L1*(-v1n - v1np1k)**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k))/4.)))/
     &        (3.*L1*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((-3*(g*(4*m1 + 5*m2)*
     &                  Sin((th1n + th1np1k - th1n*xi1 + th1np1k*xi1)/
     &                   2.) + 
     &                 m2*(3*g*
     &                   Sin((th1n + th1np1k - 2*th2n - 2*th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + 2*th2n*xi1 - 
     &                  2*th2np1k*xi1)/2.) + 
     &                   L2*(v2n + v2np1k - v2n*xi1 + v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                  th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 + xi1) - v1np1k*(1 + xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4.)))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))) - (3*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k + th1n*xi1 - 
     &                   th1np1k*xi1)/2.) + 
     &                 m2*(L2*(v2n + v2np1k + v2n*xi1 - v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 - xi1) - v1np1k*(1 - xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4. + 
     &                   3*g*Sin((th1n + th1np1k - 2*th2n - 
     &                   2*th2np1k + th1n*xi1 - th1np1k*xi1 - 
     &                   2*th2n*xi1 + 2*th2np1k*xi1)/2.))))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                 th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & )))))/18.))/12. + (L2*m2*(2*L2*(v2n + v2np1k) + 
     &       3*L1*v1n*Cos(th1n - th2n) + 
     &       3*L1*v1np1k*Cos(th1np1k - th2np1k))*
     &     ((4*(3*L1*(m1 + 3*m2)*(v1n + v1np1k)**2*
     &             Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &            (9*L2*m2*(-v2n - v2np1k)**2*
     &               Sin(th1n + th1np1k - th2n - th2np1k))/4. + 
     &            9*g*(m1 + 2*m2)*
     &             Sin(th1n + th1np1k - th2n/2. - th2np1k/2.) - 
     &            3*g*(m1 + 6*m2)*Sin((th2n + th2np1k)/2.)))/
     &        (9.*L2*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((9*g*(m1 + 2*m2)*
     &                Sin(th1n + (th2n*(-1 - xi1))/2. + 
     &                  th1np1k*(1 - xi1) - (th2np1k*(1 - xi1))/2. + 
     &                  th1n*xi1) - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k + th2n*xi1 - th2np1k*xi1)/2.)
     & + 3*L1*(m1 + 3*m2)*(v1n + v1np1k + v1n*xi1 - v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 - xi1) - v2np1k*(1 - xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4.)/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                  th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))) + (3*L1*(m1 + 3*m2)*(v1n + v1np1k - v1n*xi1 + v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 + xi1) - v2np1k*(1 + xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4. - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k - th2n*xi1 + th2np1k*xi1)/2.)
     & + 9*g*(m1 + 2*m2)*Sin(th1n + (th2n*(-1 + xi1))/2. - th1n*xi1 + 
     &                  th1np1k*(1 + xi1) - (th2np1k*(1 + xi1))/2.))/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                 th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)
     & ))))/18.))/12.))/
     &    (4.*(th2n - th2np1k)*((L1**2*
     &      (-2*g*(m1 + 2*m2)*Cos(th1n) + 2*g*(m1 + 2*m2)*Cos(th1np1k) +
     &         L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) +
     &            v1np1k*v2np1k*
     &             (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k))))**2)/
     &    (16.*(th1n - th1np1k)**2) + 
     &   (L2**2*m2**2*(2*L2*(v2n + v2np1k) + 3*L1*v1n*Cos(th1n - th2n) +
     &         3*L1*v1np1k*Cos(th1np1k - th2np1k))**2)/144. + 
     &   (L1**2*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &         3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &            v2n*Cos(th1np1k - th2np1k)))**2)/144. + 
     &   (L2**2*m2**2*(L1*v1np1k*v2np1k*
     &          (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &         L1*v1n*v2n*(Cos(th1np1k - th2n) - 
     &            Cos(th1np1k - th2np1k)) + 
     &         2*g*(-Cos(th2n) + Cos(th2np1k)))**2)/
     &    (16.*(th2n - th2np1k)**2))**2) + 
     & (L2*m2*(L1*v1np1k*v2np1k*(Cos(th1n-th2n)-Cos(th1n - th2np1k)) + 
     &  L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k)) + 
     &  2*g*(-Cos(th2n) + Cos(th2np1k)))*
     & ((L1*(v1n + v1np1k)*(-2*g*(m1 + 2*m2)*Cos(th1n) + 
     &       2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &       L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &          v1np1k*v2np1k*(Cos(th1n - th2np1k) - 
     &             Cos(th1np1k - th2np1k)))))/(8.*(th1n - th1np1k)**2) +
     &  (L1*L2*m2*v1n*v2n*(v2n + v2np1k)*
     &     (-Sin(th1np1k - th2n) + Sin(th1np1k - th2np1k)))/
     &   (8.*(th2n - th2np1k)) + 
     &  (L1*(v1n + v1np1k)*(-2*g*(m1 + 2*m2)*Sin(th1np1k) + 
     &       L2*m2*(v1n*v2n*Sin(th1np1k - th2n) + 
     &          v1np1k*v2np1k*Sin(th1np1k - th2np1k))))/
     &   (8.*(th1n - th1np1k)) - 
     &  (L1*L2*m2*v2n*Sin(th1np1k - th2np1k)*
     &     ((-4*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k)/2.) + 
     &            m2*(3*g*Sin((th1n + th1np1k - 2*th2n - 2*th2np1k)/
     &                  2.) + 
     &               L2*(v2n + v2np1k)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &               (3*L1*(-v1n - v1np1k)**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k))/4.)))/
     &        (3.*L1*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((-3*(g*(4*m1 + 5*m2)*
     &                  Sin((th1n + th1np1k - th1n*xi1 + th1np1k*xi1)/
     &                   2.) + 
     &                 m2*(3*g*
     &                   Sin((th1n + th1np1k - 2*th2n - 2*th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + 2*th2n*xi1 - 
     &                   2*th2np1k*xi1)/2.) + 
     &                   L2*(v2n + v2np1k - v2n*xi1 + v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 + xi1) - v1np1k*(1 + xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4.)))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))) - (3*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k + th1n*xi1 - 
     &                   th1np1k*xi1)/2.) + 
     &                 m2*(L2*(v2n + v2np1k + v2n*xi1 - v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &                    (3*L1*(v1n*(-1 - xi1) - v1np1k*(1 - xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4. + 
     &                    3*g*
     &                    Sin((th1n + th1np1k - 2*th2n - 2*th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - 2*th2n*xi1 + 
     &                   2*th2np1k*xi1)/2.))))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                 th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1)
     & ))))/18.))/4. + (L1*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &       3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &          v2n*Cos(th1np1k - th2np1k)))*
     &     ((-4*((g*(4*m1 + 5*m2)*Cos((th1n + th1np1k)/2.))/2. + 
     &            m2*((3*g*Cos((th1n + th1np1k - 2*th2n - 2*th2np1k)/
     &                   2.))/2. + 
     &               (L2*(v2n + v2np1k)**2*
     &                  Cos((th1n + th1np1k - th2n - th2np1k)/2.))/2. + 
     &               (3*L1*(-v1n - v1np1k)**2*
     &                  Cos(th1n + th1np1k - th2n - th2np1k))/4.)))/
     &        (3.*L1*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (12*m2*Sin(th1n + th1np1k - th2n - th2np1k)*
     &          (g*(4*m1 + 5*m2)*Sin((th1n + th1np1k)/2.) + 
     &            m2*(3*g*Sin((th1n + th1np1k - 2*th2n - 2*th2np1k)/
     &                  2.) + 
     &               L2*(v2n + v2np1k)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &               (3*L1*(-v1n - v1np1k)**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k))/4.)))/
     &        (L1*(8*m1 + 15*m2 - 
     &             9*m2*Cos(th1n + th1np1k - th2n - th2np1k))**2) + 
     &       (5*((-3*((g*(4*m1 + 5*m2)*(1 + xi1)*
     &                   Cos((th1n + th1np1k - th1n*xi1 + 
     &                  th1np1k*xi1)/2.))/2. + 
     &                 m2*((3*g*(1 + xi1)*
     &                   Cos((th1n + th1np1k - 2*th2n - 2*th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + 2*th2n*xi1 - 
     &                  2*th2np1k*xi1)/2.))/2. + 
     &                   (L2*(1 + xi1)*
     &                   (v2n + v2np1k - v2n*xi1 + v2np1k*xi1)**2*
     &                   Cos((th1n + th1np1k - th2n - th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                  th2np1k*xi1)/2.))/2. + 
     &                   (3*L1*(1 + xi1)*
     &                   (v1n*(-1 + xi1) - v1np1k*(1 + xi1))**2*
     &                   Cos(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4.)))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))) - (3*((g*(4*m1 + 5*m2)*(1 - xi1)*
     &                   Cos((th1n + th1np1k + th1n*xi1 - 
     &                  th1np1k*xi1)/2.))/2. + 
     &                 m2*((L2*(1 - xi1)*
     &                   (v2n + v2np1k + v2n*xi1 - v2np1k*xi1)**2*
     &                   Cos((th1n + th1np1k - th2n - th2np1k + 
     &                  th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                  th2np1k*xi1)/2.))/2. + 
     &                   (3*L1*(v1n*(-1 - xi1) - v1np1k*(1 - xi1))**2*
     &                   (1 - xi1)*
     &                   Cos(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4. + 
     &                   (3*g*(1 - xi1)*
     &                   Cos((th1n + th1np1k - 2*th2n - 2*th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - 2*th2n*xi1 + 
     &                   2*th2np1k*xi1)/2.))/2.)))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                 th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))) + (27*m2*(1 + xi1)*Sin(th1n + th1np1k - th2n - th2np1k - 
     &                 th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)*
     &               (g*(4*m1 + 5*m2)*
     &                  Sin((th1n + th1np1k - th1n*xi1 + th1np1k*xi1)/
     &                   2.) + 
     &                 m2*(3*g*
     &                   Sin((th1n + th1np1k - 2*th2n - 2*th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + 2*th2n*xi1 - 
     &                   2*th2np1k*xi1)/2.) + 
     &                   L2*(v2n + v2np1k - v2n*xi1 + v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 + xi1) - v1np1k*(1 + xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4.)))/
     &             (L1*(8*m1 + 15*m2 - 
     &                  9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))**2) + (27*m2*(1 - xi1)*Sin(th1n + th1np1k - th2n - th2np1k + 
     &                 th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1)*
     &               (g*(4*m1 + 5*m2)*
     &                  Sin((th1n + th1np1k + th1n*xi1 - th1np1k*xi1)/
     &                   2.) + 
     &                 m2*(L2*(v2n + v2np1k + v2n*xi1 - v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &                    (3*L1*(v1n*(-1 - xi1) - v1np1k*(1 - xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4. + 
     &                    3*g*
     &                    Sin((th1n + th1np1k - 2*th2n - 2*th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - 2*th2n*xi1 + 
     &                   2*th2np1k*xi1)/2.))))/
     &             (L1*(8*m1 + 15*m2 - 
     &                  9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                 th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))**2)))/18.))/12. - (L1*L2*m2*v1np1k*Sin(th1np1k - th2np1k)*
     &     ((4*(3*L1*(m1 + 3*m2)*(v1n + v1np1k)**2*
     &             Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &            (9*L2*m2*(-v2n - v2np1k)**2*
     &               Sin(th1n + th1np1k - th2n - th2np1k))/4. + 
     &            9*g*(m1 + 2*m2)*
     &             Sin(th1n + th1np1k - th2n/2. - th2np1k/2.) - 
     &            3*g*(m1 + 6*m2)*Sin((th2n + th2np1k)/2.)))/
     &        (9.*L2*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((9*g*(m1 + 2*m2)*
     &                Sin(th1n + (th2n*(-1 - xi1))/2. + 
     &                  th1np1k*(1 - xi1) - (th2np1k*(1 - xi1))/2. + 
     &                  th1n*xi1) - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k + th2n*xi1 - th2np1k*xi1)/2.)
     & + 3*L1*(m1 + 3*m2)*(v1n + v1np1k + v1n*xi1 - v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 - xi1) - v2np1k*(1 - xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4.)/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                  th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))) + (3*L1*(m1 + 3*m2)*(v1n + v1np1k - v1n*xi1 + v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 + xi1) - v2np1k*(1 + xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4. - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k - th2n*xi1 + th2np1k*xi1)/2.)
     & + 9*g*(m1 + 2*m2)*Sin(th1n + (th2n*(-1 + xi1))/2. - th1n*xi1 + 
     &                  th1np1k*(1 + xi1) - (th2np1k*(1 + xi1))/2.))/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                 th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)
     & ))))/18.))/4.+(L2*m2*(2*L2*(v2n+v2np1k)+3*L1*v1n*Cos(th1n-th2n)+
     &       3*L1*v1np1k*Cos(th1np1k - th2np1k))*
     &     ((4*((3*L1*(m1 + 3*m2)*(v1n + v1np1k)**2*
     &               Cos((th1n + th1np1k - th2n - th2np1k)/2.))/2. + 
     &            (9*L2*m2*(-v2n - v2np1k)**2*
     &               Cos(th1n + th1np1k - th2n - th2np1k))/4. + 
     &            9*g*(m1 + 2*m2)*
     &             Cos(th1n + th1np1k - th2n/2. - th2np1k/2.)))/
     &        (9.*L2*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) - 
     &       (4*m2*Sin(th1n + th1np1k - th2n - th2np1k)*
     &          (3*L1*(m1 + 3*m2)*(v1n + v1np1k)**2*
     &             Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &            (9*L2*m2*(-v2n - v2np1k)**2*
     &               Sin(th1n + th1np1k - th2n - th2np1k))/4. + 
     &            9*g*(m1 + 2*m2)*
     &             Sin(th1n + th1np1k - th2n/2. - th2np1k/2.) - 
     &            3*g*(m1 + 6*m2)*Sin((th2n + th2np1k)/2.)))/
     &        (L2*(8*m1 + 15*m2 - 
     &             9*m2*Cos(th1n + th1np1k - th2n - th2np1k))**2) + 
     &       (5*((9*g*(m1 + 2*m2)*(1 - xi1)*
     &                Cos(th1n + (th2n*(-1 - xi1))/2. + 
     &                  th1np1k*(1 - xi1) - (th2np1k*(1 - xi1))/2. + 
     &                  th1n*xi1) + 
     &               (3*L1*(m1 + 3*m2)*(1 - xi1)*
     &                  (v1n + v1np1k + v1n*xi1 - v1np1k*xi1)**2*
     &                  Cos((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.))/2. + 
     &               (9*L2*m2*(v2n*(-1 - xi1) - v2np1k*(1 - xi1))**2*
     &                  (1 - xi1)*
     &                  Cos(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))/4.)/(L2*(8*m1 + 15*m2 - 9*m2*
     &                  Cos(th1n + th1np1k - th2n - th2np1k + th1n*xi1 -
     &                    th1np1k*xi1 - th2n*xi1 + th2np1k*xi1))) + 
     &            ((3*L1*(m1 + 3*m2)*(1 + xi1)*
     &                  (v1n + v1np1k - v1n*xi1 + v1np1k*xi1)**2*
     &                  Cos((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.))/2. + 
     &               (9*L2*m2*(1 + xi1)*
     &                  (v2n*(-1 + xi1) - v2np1k*(1 + xi1))**2*
     &                  Cos(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4. + 
     &               9*g*(m1 + 2*m2)*(1 + xi1)*
     &                Cos(th1n + (th2n*(-1 + xi1))/2. - th1n*xi1 + 
     &                  th1np1k*(1 + xi1) - (th2np1k*(1 + xi1))/2.))/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)
     & ))-(9*m2*(1-xi1)*Sin(th1n + th1np1k - th2n - th2np1k + th1n*xi1-
     &                 th1np1k*xi1 - th2n*xi1 + th2np1k*xi1)*
     &               (9*g*(m1 + 2*m2)*
     &                  Sin(th1n + (th2n*(-1 - xi1))/2. + 
     &                   th1np1k*(1 - xi1) - (th2np1k*(1 - xi1))/2. + 
     &                   th1n*xi1) - 
     &                 3*g*(m1 + 6*m2)*
     &                  Sin((th2n + th2np1k + th2n*xi1 - th2np1k*xi1)/
     &                   2.) + 
     &                 3*L1*(m1 + 3*m2)*
     &                  (v1n + v1np1k + v1n*xi1 - v1np1k*xi1)**2*
     &                  Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &                 (9*L2*m2*(v2n*(-1 - xi1) - v2np1k*(1 - xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))/4.))/(L2*(8*m1 + 15*m2 - 9*m2*
     &                   Cos(th1n + th1np1k - th2n - th2np1k + 
     &                th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))**2) - (9*m2*(1 + xi1)*Sin(th1n + th1np1k - th2n - th2np1k - 
     &                 th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)*
     &               (3*L1*(m1 + 3*m2)*
     &                  (v1n + v1np1k - v1n*xi1 + v1np1k*xi1)**2*
     &                  Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & )/2.) + (9*L2*m2*(v2n*(-1 + xi1) - v2np1k*(1 + xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))/4.-3*g*(m1+6*m2)*Sin((th2n + th2np1k - th2n*xi1+th2np1k*xi1)/
     &                   2.) + 
     &                 9*g*(m1 + 2*m2)*
     &                  Sin(th1n + (th2n*(-1 + xi1))/2. - th1n*xi1 + 
     &                    th1np1k*(1 + xi1) - (th2np1k*(1 + xi1))/2.)))/
     &             (L2*(8*m1 + 15*m2 - 
     &                  9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)
     & )**2)))/18.))/12.))/
     & (4.*(th2n - th2np1k)*((L1**2*
     &     (-2*g*(m1 + 2*m2)*Cos(th1n) + 2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &        L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &           v1np1k*v2np1k*
     &            (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k))))**2)/
     &   (16.*(th1n - th1np1k)**2) + 
     &  (L2**2*m2**2*(2*L2*(v2n + v2np1k) + 3*L1*v1n*Cos(th1n - th2n) +
     &        3*L1*v1np1k*Cos(th1np1k - th2np1k))**2)/144. + 
     &  (L1**2*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &        3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &           v2n*Cos(th1np1k - th2np1k)))**2)/144. + 
     &  (L2**2*m2**2*(L1*v1np1k*v2np1k*
     &         (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &        L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k))+
     &  2*g*(-Cos(th2n) + Cos(th2np1k)))**2)/(16.*(th2n - th2np1k)**2)))

      m22= 1.0d0/tdelta + (L2*m2*(L1*v1np1k*v2np1k*
     &   (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &  L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k)) + 
     &  2*g*(-Cos(th2n) + Cos(th2np1k)))*
     & ((L1*(v1n + v1np1k)*(-2*g*(m1 + 2*m2)*Cos(th1n) + 
     &       2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &       L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &          v1np1k*v2np1k*(Cos(th1n - th2np1k) - 
     &             Cos(th1np1k - th2np1k)))))/(8.*(th1n - th1np1k)) + 
     &  (L2*m2*(v2n + v2np1k)*(L1*v1np1k*v2np1k*
     &        (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &       L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k)) +
     &       2*g*(-Cos(th2n) + Cos(th2np1k))))/(8.*(th2n - th2np1k)) + 
     &  (L1*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &       3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &          v2n*Cos(th1np1k - th2np1k)))*
     &     ((-4*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k)/2.) + 
     &            m2*(3*g*Sin((th1n + th1np1k - 2*th2n - 2*th2np1k)/
     &                  2.) + 
     &               L2*(v2n + v2np1k)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &               (3*L1*(-v1n - v1np1k)**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k))/4.)))/
     &        (3.*L1*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((-3*(g*(4*m1 + 5*m2)*
     &                  Sin((th1n + th1np1k - th1n*xi1 + th1np1k*xi1)/
     &                   2.) + 
     &                 m2*(3*g*
     &                   Sin((th1n + th1np1k - 2*th2n - 2*th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + 2*th2n*xi1 - 
     &                  2*th2np1k*xi1)/2.) + 
     &                   L2*(v2n + v2np1k - v2n*xi1 + v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                  th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 + xi1) - v1np1k*(1 + xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4.)))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))) - (3*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k + th1n*xi1 - 
     &                   th1np1k*xi1)/2.) + 
     &                 m2*(L2*(v2n + v2np1k + v2n*xi1 - v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 - xi1) - v1np1k*(1 - xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4. + 
     &                   3*g*Sin((th1n + th1np1k - 2*th2n - 
     &                   2*th2np1k + th1n*xi1 - th1np1k*xi1 - 
     &                   2*th2n*xi1 + 2*th2np1k*xi1)/2.))))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                  th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & )))))/18.))/12. + (L2*m2*(2*L2*(v2n + v2np1k) + 
     &       3*L1*v1n*Cos(th1n - th2n) + 
     &       3*L1*v1np1k*Cos(th1np1k - th2np1k))*
     &     ((4*(3*L1*(m1 + 3*m2)*(v1n + v1np1k)**2*
     &             Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &            (9*L2*m2*(-v2n - v2np1k)**2*
     &               Sin(th1n + th1np1k - th2n - th2np1k))/4. + 
     &            9*g*(m1 + 2*m2)*
     &             Sin(th1n + th1np1k - th2n/2. - th2np1k/2.) - 
     &            3*g*(m1 + 6*m2)*Sin((th2n + th2np1k)/2.)))/
     &        (9.*L2*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((9*g*(m1 + 2*m2)*
     &                Sin(th1n + (th2n*(-1 - xi1))/2. + 
     &                  th1np1k*(1 - xi1) - (th2np1k*(1 - xi1))/2. + 
     &                  th1n*xi1) - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k + th2n*xi1 - th2np1k*xi1)/2.)
     & + 3*L1*(m1 + 3*m2)*(v1n + v1np1k + v1n*xi1 - v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 - xi1) - v2np1k*(1 - xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4.)/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                  th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))) + (3*L1*(m1 + 3*m2)*(v1n + v1np1k - v1n*xi1 + v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 + xi1) - v2np1k*(1 + xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4. - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k - th2n*xi1 + th2np1k*xi1)/2.)
     & + 9*g*(m1 + 2*m2)*Sin(th1n + (th2n*(-1 + xi1))/2. - th1n*xi1 + 
     &                  th1np1k*(1 + xi1) - (th2np1k*(1 + xi1))/2.))/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)
     & ))))/18.))/12.))/
     & (4.*(th2n - th2np1k)**2*((L1**2*
     &     (-2*g*(m1 + 2*m2)*Cos(th1n) + 2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &        L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &           v1np1k*v2np1k*
     &            (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k))))**2)/
     &   (16.*(th1n - th1np1k)**2) + 
     &  (L2**2*m2**2*(2*L2*(v2n + v2np1k) + 3*L1*v1n*Cos(th1n - th2n) + 
     &        3*L1*v1np1k*Cos(th1np1k - th2np1k))**2)/144. + 
     &  (L1**2*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &        3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &           v2n*Cos(th1np1k - th2np1k)))**2)/144. + 
     &  (L2**2*m2**2*(L1*v1np1k*v2np1k*
     &         (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &        L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k))+
     & 2*g*(-Cos(th2n) + Cos(th2np1k)))**2)/(16.*(th2n - th2np1k)**2))
     & ) + (L2*m2*(-(L1*v1np1k*v2np1k*Sin(th1n - th2np1k)) - 
     &  L1*v1n*v2n*Sin(th1np1k - th2np1k) - 2*g*Sin(th2np1k))*
     & ((L1*(v1n + v1np1k)*(-2*g*(m1 + 2*m2)*Cos(th1n) + 
     &       2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &       L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &          v1np1k*v2np1k*(Cos(th1n - th2np1k) - 
     &             Cos(th1np1k - th2np1k)))))/(8.*(th1n - th1np1k)) + 
     &  (L2*m2*(v2n + v2np1k)*(L1*v1np1k*v2np1k*
     &        (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &       L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k)) +
     &       2*g*(-Cos(th2n) + Cos(th2np1k))))/(8.*(th2n - th2np1k)) + 
     &  (L1*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &       3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &          v2n*Cos(th1np1k - th2np1k)))*
     &     ((-4*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k)/2.) + 
     &            m2*(3*g*Sin((th1n + th1np1k - 2*th2n - 2*th2np1k)/
     &                  2.) + 
     &               L2*(v2n + v2np1k)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &               (3*L1*(-v1n - v1np1k)**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k))/4.)))/
     &        (3.*L1*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((-3*(g*(4*m1 + 5*m2)*
     &                  Sin((th1n + th1np1k - th1n*xi1 + th1np1k*xi1)/
     &                   2.) + 
     &                 m2*(3*g*
     &                   Sin((th1n + th1np1k - 2*th2n - 2*th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + 2*th2n*xi1 - 
     &                  2*th2np1k*xi1)/2.) + 
     &                   L2*(v2n + v2np1k - v2n*xi1 + v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                  th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 + xi1) - v1np1k*(1 + xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4.)))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))) - (3*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k + th1n*xi1 - 
     &                   th1np1k*xi1)/2.) + 
     &                 m2*(L2*(v2n + v2np1k + v2n*xi1 - v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 - xi1) - v1np1k*(1 - xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4. + 
     &                   3*g*Sin((th1n + th1np1k - 2*th2n - 
     &                   2*th2np1k + th1n*xi1 - th1np1k*xi1 - 
     &                   2*th2n*xi1 + 2*th2np1k*xi1)/2.))))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                 th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & )))))/18.))/12. + (L2*m2*(2*L2*(v2n + v2np1k) + 
     &       3*L1*v1n*Cos(th1n - th2n) + 
     &       3*L1*v1np1k*Cos(th1np1k - th2np1k))*
     &     ((4*(3*L1*(m1 + 3*m2)*(v1n + v1np1k)**2*
     &             Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &            (9*L2*m2*(-v2n - v2np1k)**2*
     &               Sin(th1n + th1np1k - th2n - th2np1k))/4. + 
     &            9*g*(m1 + 2*m2)*
     &             Sin(th1n + th1np1k - th2n/2. - th2np1k/2.) - 
     &            3*g*(m1 + 6*m2)*Sin((th2n + th2np1k)/2.)))/
     &        (9.*L2*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((9*g*(m1 + 2*m2)*
     &                Sin(th1n + (th2n*(-1 - xi1))/2. + 
     &                  th1np1k*(1 - xi1) - (th2np1k*(1 - xi1))/2. + 
     &                  th1n*xi1) - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k + th2n*xi1 - th2np1k*xi1)/2.)
     & + 3*L1*(m1 + 3*m2)*(v1n + v1np1k + v1n*xi1 - v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 - xi1) - v2np1k*(1 - xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4.)/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                  th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))) + (3*L1*(m1 + 3*m2)*(v1n + v1np1k - v1n*xi1 + v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 + xi1) - v2np1k*(1 + xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4. - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k - th2n*xi1 + th2np1k*xi1)/2.)
     & + 9*g*(m1 + 2*m2)*Sin(th1n + (th2n*(-1 + xi1))/2. - th1n*xi1 + 
     &                  th1np1k*(1 + xi1) - (th2np1k*(1 + xi1))/2.))/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                 th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)
     & ))))/18.))/12.))/
     &    (4.*(th2n - th2np1k)*((L1**2*
     &     (-2*g*(m1 + 2*m2)*Cos(th1n) + 2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &        L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &           v1np1k*v2np1k*
     &            (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k))))**2)/
     &   (16.*(th1n - th1np1k)**2) + 
     &  (L2**2*m2**2*(2*L2*(v2n + v2np1k) + 3*L1*v1n*Cos(th1n - th2n) +
     &        3*L1*v1np1k*Cos(th1np1k - th2np1k))**2)/144. + 
     &  (L1**2*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &        3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &           v2n*Cos(th1np1k - th2np1k)))**2)/144. + 
     &  (L2**2*m2**2*(L1*v1np1k*v2np1k*
     &         (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &        L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k))+
     &  2*g*(-Cos(th2n) + Cos(th2np1k)))**2)/(16.*(th2n - th2np1k)**2))
     & )-(L2*m2*(L1*v1np1k*v2np1k*(Cos(th1n-th2n)-Cos(th1n-th2np1k))+
     &  L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k)) + 
     &  2*g*(-Cos(th2n) + Cos(th2np1k)))*
     & ((L2**2*m2**2*(L1*v1np1k*v2np1k*
     &         (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &        L1*v1n*v2n*(Cos(th1np1k - th2n) - 
     &           Cos(th1np1k - th2np1k)) + 
     &        2*g*(-Cos(th2n) + Cos(th2np1k)))**2)/
     &   (8.*(th2n - th2np1k)**3) + 
     &  (L1**2*L2*m2*v1np1k*v2np1k*
     &     (-2*g*(m1 + 2*m2)*Cos(th1n) + 2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &       L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &          v1np1k*v2np1k*
     &           (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k))))*
     &     (Sin(th1n - th2np1k) - Sin(th1np1k - th2np1k)))/
     &   (8.*(th1n - th1np1k)**2) + 
     &  (L1*L2**2*m2**2*v1np1k*
     &     (2*L2*(v2n + v2np1k) + 3*L1*v1n*Cos(th1n - th2n) + 
     &     3*L1*v1np1k*Cos(th1np1k - th2np1k))*Sin(th1np1k - th2np1k))/
     &   24. + (L1**2*L2*m2*v2n*
     &     (2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &       3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &        v2n*Cos(th1np1k - th2np1k)))*Sin(th1np1k - th2np1k))/24.+
     &  (L2**2*m2**2*(L1*v1np1k*v2np1k*
     &        (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &     L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k)) + 
     &       2*g*(-Cos(th2n) + Cos(th2np1k)))*
     &     (-(L1*v1np1k*v2np1k*Sin(th1n - th2np1k)) - 
     &       L1*v1n*v2n*Sin(th1np1k - th2np1k) - 2*g*Sin(th2np1k)))/
     &   (8.*(th2n - th2np1k)**2))*
     & ((L1*(v1n + v1np1k)*(-2*g*(m1 + 2*m2)*Cos(th1n) + 
     &       2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &       L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) +
     &          v1np1k*v2np1k*(Cos(th1n - th2np1k) - 
     &             Cos(th1np1k - th2np1k)))))/(8.*(th1n - th1np1k)) + 
     &  (L2*m2*(v2n + v2np1k)*(L1*v1np1k*v2np1k*
     &        (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &       L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k))+
     &       2*g*(-Cos(th2n) + Cos(th2np1k))))/(8.*(th2n - th2np1k)) + 
     &  (L1*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &       3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &          v2n*Cos(th1np1k - th2np1k)))*
     &     ((-4*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k)/2.) + 
     &            m2*(3*g*Sin((th1n + th1np1k - 2*th2n - 2*th2np1k)/
     &                  2.) + 
     &               L2*(v2n + v2np1k)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &               (3*L1*(-v1n - v1np1k)**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k))/4.)))/
     &        (3.*L1*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((-3*(g*(4*m1 + 5*m2)*
     &                  Sin((th1n + th1np1k - th1n*xi1 + th1np1k*xi1)/
     &                   2.) + 
     &                 m2*(3*g*
     &                   Sin((th1n + th1np1k - 2*th2n - 2*th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + 2*th2n*xi1 - 
     &                  2*th2np1k*xi1)/2.) + 
     &                   L2*(v2n + v2np1k - v2n*xi1 + v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                  th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 + xi1) - v1np1k*(1 + xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4.)))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))) - (3*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k + th1n*xi1 - 
     &                   th1np1k*xi1)/2.) + 
     &                 m2*(L2*(v2n + v2np1k + v2n*xi1 - v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 - xi1) - v1np1k*(1 - xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4. + 
     &                   3*g*Sin((th1n + th1np1k - 2*th2n - 
     &                   2*th2np1k + th1n*xi1 - th1np1k*xi1 - 
     &                   2*th2n*xi1 + 2*th2np1k*xi1)/2.))))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                 th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & )))))/18.))/12. + (L2*m2*(2*L2*(v2n + v2np1k) + 
     &       3*L1*v1n*Cos(th1n - th2n) + 
     &       3*L1*v1np1k*Cos(th1np1k - th2np1k))*
     &     ((4*(3*L1*(m1 + 3*m2)*(v1n + v1np1k)**2*
     &             Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &            (9*L2*m2*(-v2n - v2np1k)**2*
     &               Sin(th1n + th1np1k - th2n - th2np1k))/4. + 
     &            9*g*(m1 + 2*m2)*
     &             Sin(th1n + th1np1k - th2n/2. - th2np1k/2.) - 
     &            3*g*(m1 + 6*m2)*Sin((th2n + th2np1k)/2.)))/
     &        (9.*L2*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((9*g*(m1 + 2*m2)*
     &                Sin(th1n + (th2n*(-1 - xi1))/2. + 
     &                  th1np1k*(1 - xi1) - (th2np1k*(1 - xi1))/2. + 
     &                  th1n*xi1) - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k + th2n*xi1 - th2np1k*xi1)/2.)
     & + 3*L1*(m1 + 3*m2)*(v1n + v1np1k + v1n*xi1 - v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 - xi1) - v2np1k*(1 - xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4.)/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                  th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))) + (3*L1*(m1 + 3*m2)*(v1n + v1np1k - v1n*xi1 + v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 + xi1) - v2np1k*(1 + xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4. - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k - th2n*xi1 + th2np1k*xi1)/2.)
     & + 9*g*(m1 + 2*m2)*Sin(th1n + (th2n*(-1 + xi1))/2. - th1n*xi1 + 
     &                  th1np1k*(1 + xi1) - (th2np1k*(1 + xi1))/2.))/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                 th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)
     & ))))/18.))/12.))/
     &    (4.*(th2n - th2np1k)*((L1**2*
     &      (-2*g*(m1 + 2*m2)*Cos(th1n) + 2*g*(m1 + 2*m2)*Cos(th1np1k) +
     &         L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) +
     &            v1np1k*v2np1k*
     &             (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k))))**2)/
     &    (16.*(th1n - th1np1k)**2) + 
     &   (L2**2*m2**2*(2*L2*(v2n + v2np1k) + 3*L1*v1n*Cos(th1n - th2n) +
     &         3*L1*v1np1k*Cos(th1np1k - th2np1k))**2)/144. + 
     &   (L1**2*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &         3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &            v2n*Cos(th1np1k - th2np1k)))**2)/144. + 
     &   (L2**2*m2**2*(L1*v1np1k*v2np1k*
     &          (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &         L1*v1n*v2n*(Cos(th1np1k - th2n) - 
     &            Cos(th1np1k - th2np1k)) + 
     &         2*g*(-Cos(th2n) + Cos(th2np1k)))**2)/
     &    (16.*(th2n - th2np1k)**2))**2) + 
     & (L2*m2*(L1*v1np1k*v2np1k*(Cos(th1n-th2n) - Cos(th1n - th2np1k))+
     &  L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k)) + 
     &  2*g*(-Cos(th2n) + Cos(th2np1k)))*
     & ((L2*m2*(v2n + v2np1k)*(L1*v1np1k*v2np1k*
     &        (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &       L1*v1n*v2n*(Cos(th1np1k - th2n)-Cos(th1np1k-th2np1k)) + 
     &    2*g*(-Cos(th2n) + Cos(th2np1k))))/(8.*(th2n-th2np1k)**2) + 
     &  (L1*L2*m2*v1np1k*(v1n + v1np1k)*v2np1k*
     &     (Sin(th1n - th2np1k) - Sin(th1np1k - th2np1k)))/
     &   (8.*(th1n - th1np1k)) + 
     &  (L2*m2*(v2n + v2np1k)*(-(L1*v1np1k*v2np1k*Sin(th1n - th2np1k)) -
     &       L1*v1n*v2n*Sin(th1np1k - th2np1k) - 2*g*Sin(th2np1k)))/
     &   (8.*(th2n - th2np1k)) + 
     &  (L1*L2*m2*v2n*Sin(th1np1k - th2np1k)*
     &     ((-4*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k)/2.) + 
     &            m2*(3*g*Sin((th1n + th1np1k - 2*th2n - 2*th2np1k)/
     &                  2.) + 
     &               L2*(v2n + v2np1k)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &               (3*L1*(-v1n - v1np1k)**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k))/4.)))/
     &        (3.*L1*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((-3*(g*(4*m1 + 5*m2)*
     &                  Sin((th1n + th1np1k - th1n*xi1 + th1np1k*xi1)/
     &                   2.) + 
     &                 m2*(3*g*
     &                   Sin((th1n + th1np1k - 2*th2n - 2*th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + 2*th2n*xi1 - 
     &                   2*th2np1k*xi1)/2.) + 
     &                   L2*(v2n + v2np1k - v2n*xi1 + v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 + xi1) - v1np1k*(1 + xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4.)))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                 th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))) - (3*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k + th1n*xi1 - 
     &                   th1np1k*xi1)/2.) + 
     &                 m2*(L2*(v2n + v2np1k + v2n*xi1 - v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &                    (3*L1*(v1n*(-1 - xi1) - v1np1k*(1 - xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4. + 
     &                    3*g*
     &                    Sin((th1n + th1np1k - 2*th2n - 2*th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - 2*th2n*xi1 + 
     &                   2*th2np1k*xi1)/2.))))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1)
     & ))))/18.))/4. + (L1*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &       3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &          v2n*Cos(th1np1k - th2np1k)))*
     &     ((-4*m2*(-3*g*Cos((th1n + th1np1k - 2*th2n - 2*th2np1k)/2.) -
     &            (L2*(v2n + v2np1k)**2*
     &               Cos((th1n + th1np1k - th2n - th2np1k)/2.))/2. - 
     &            (3*L1*(-v1n - v1np1k)**2*
     &               Cos(th1n + th1np1k - th2n - th2np1k))/4.))/
     &        (3.*L1*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) - 
     &       (12*m2*Sin(th1n + th1np1k - th2n - th2np1k)*
     &          (g*(4*m1 + 5*m2)*Sin((th1n + th1np1k)/2.) + 
     &            m2*(3*g*Sin((th1n + th1np1k - 2*th2n - 2*th2np1k)/
     &                  2.) + 
     &               L2*(v2n + v2np1k)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &               (3*L1*(-v1n - v1np1k)**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k))/4.)))/
     &        (L1*(8*m1 + 15*m2 - 
     &             9*m2*Cos(th1n + th1np1k - th2n - th2np1k))**2) + 
     &       (5*((-3*m2*((3*g*(-2 - 2*xi1)*
     &                   Cos((th1n + th1np1k - 2*th2n - 2*th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + 2*th2n*xi1 - 
     &                  2*th2np1k*xi1)/2.))/2. + 
     &                 (L2*(-1 - xi1)*
     &                   (v2n + v2np1k - v2n*xi1 + v2np1k*xi1)**2*
     &                   Cos((th1n + th1np1k - th2n - th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                  th2np1k*xi1)/2.))/2. + 
     &                 (3*L1*(-1 - xi1)*
     &                   (v1n*(-1 + xi1) - v1np1k*(1 + xi1))**2*
     &                   Cos(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4.))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                 th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & )))-(3*m2*((L2*(-1+xi1)*(v2n+v2np1k + v2n*xi1 - v2np1k*xi1)**2*
     &                   Cos((th1n + th1np1k - th2n - th2np1k + 
     &                  th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                  th2np1k*xi1)/2.))/2. + 
     &                 (3*L1*(v1n*(-1 - xi1) - v1np1k*(1 - xi1))**2*
     &                   (-1 + xi1)*
     &                   Cos(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4. + 
     &                 (3*g*(-2 + 2*xi1)*
     &                   Cos((th1n + th1np1k - 2*th2n - 2*th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - 2*th2n*xi1 + 
     &                   2*th2np1k*xi1)/2.))/2.))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                 th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))) + (27*m2*(-1 - xi1)*Sin(th1n + th1np1k - th2n - th2np1k - 
     &                 th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)*
     &               (g*(4*m1 + 5*m2)*
     &                  Sin((th1n + th1np1k - th1n*xi1 + th1np1k*xi1)/
     &                   2.) + 
     &                 m2*(3*g*
     &                   Sin((th1n + th1np1k - 2*th2n - 2*th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + 2*th2n*xi1 - 
     &                   2*th2np1k*xi1)/2.) + 
     &                   L2*(v2n + v2np1k - v2n*xi1 + v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 + xi1) - v1np1k*(1 + xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4.)))/
     &             (L1*(8*m1 + 15*m2 - 
     &                  9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))**2) + (27*m2*(-1 + xi1)*Sin(th1n + th1np1k - th2n - th2np1k + 
     &                 th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1)*
     &               (g*(4*m1 + 5*m2)*
     &                  Sin((th1n + th1np1k + th1n*xi1 - th1np1k*xi1)/
     &                   2.) + 
     &                 m2*(L2*(v2n + v2np1k + v2n*xi1 - v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &                    (3*L1*(v1n*(-1 - xi1) - v1np1k*(1 - xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4. + 
     &                    3*g*
     &                    Sin((th1n + th1np1k - 2*th2n - 2*th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - 2*th2n*xi1 + 
     &                   2*th2np1k*xi1)/2.))))/
     &             (L1*(8*m1 + 15*m2 - 
     &                  9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                 th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))**2)))/18.))/12. + (L1*L2*m2*v1np1k*Sin(th1np1k - th2np1k)*
     &     ((4*(3*L1*(m1 + 3*m2)*(v1n + v1np1k)**2*
     &             Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &            (9*L2*m2*(-v2n - v2np1k)**2*
     &               Sin(th1n + th1np1k - th2n - th2np1k))/4. + 
     &            9*g*(m1 + 2*m2)*
     &             Sin(th1n + th1np1k - th2n/2. - th2np1k/2.) - 
     &            3*g*(m1 + 6*m2)*Sin((th2n + th2np1k)/2.)))/
     &        (9.*L2*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((9*g*(m1 + 2*m2)*
     &                Sin(th1n + (th2n*(-1 - xi1))/2. + 
     &                  th1np1k*(1 - xi1) - (th2np1k*(1 - xi1))/2. + 
     &                  th1n*xi1) - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k + th2n*xi1 - th2np1k*xi1)/2.)
     & + 3*L1*(m1 + 3*m2)*(v1n + v1np1k + v1n*xi1 - v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 - xi1) - v2np1k*(1 - xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4.)/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                 th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))) + (3*L1*(m1 + 3*m2)*(v1n + v1np1k - v1n*xi1 + v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 + xi1) - v2np1k*(1 + xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4. - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k - th2n*xi1 + th2np1k*xi1)/2.)
     & + 9*g*(m1 + 2*m2)*Sin(th1n + (th2n*(-1 + xi1))/2. - th1n*xi1 + 
     &                  th1np1k*(1 + xi1) - (th2np1k*(1 + xi1))/2.))/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                 th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)
     & ))))/18.))/4.+(L2*m2*(2*L2*(v2n+v2np1k)+3*L1*v1n*Cos(th1n-th2n)+
     &       3*L1*v1np1k*Cos(th1np1k - th2np1k))*
     &     ((4*((-3*L1*(m1 + 3*m2)*(v1n + v1np1k)**2*
     &               Cos((th1n + th1np1k - th2n - th2np1k)/2.))/2. - 
     &            (9*L2*m2*(-v2n - v2np1k)**2*
     &               Cos(th1n + th1np1k - th2n - th2np1k))/4. - 
     &            (9*g*(m1 + 2*m2)*
     &               Cos(th1n + th1np1k - th2n/2. - th2np1k/2.))/2. - 
     &            (3*g*(m1 + 6*m2)*Cos((th2n + th2np1k)/2.))/2.))/
     &        (9.*L2*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (4*m2*Sin(th1n + th1np1k - th2n - th2np1k)*
     &          (3*L1*(m1 + 3*m2)*(v1n + v1np1k)**2*
     &             Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &            (9*L2*m2*(-v2n - v2np1k)**2*
     &               Sin(th1n + th1np1k - th2n - th2np1k))/4. + 
     &            9*g*(m1 + 2*m2)*
     &             Sin(th1n + th1np1k - th2n/2. - th2np1k/2.) - 
     &            3*g*(m1 + 6*m2)*Sin((th2n + th2np1k)/2.)))/
     &        (L2*(8*m1 + 15*m2 - 
     &             9*m2*Cos(th1n + th1np1k - th2n - th2np1k))**2) + 
     &       (5*(((9*g*(m1 + 2*m2)*(-1 + xi1)*
     &                  Cos(th1n + (th2n*(-1 - xi1))/2. + 
     &                   th1np1k*(1 - xi1) - (th2np1k*(1 - xi1))/2. + 
     &                   th1n*xi1))/2. - 
     &               (3*g*(m1 + 6*m2)*(1 - xi1)*
     &                  Cos((th2n + th2np1k + th2n*xi1 - th2np1k*xi1)/
     &                   2.))/2. + 
     &               (3*L1*(m1 + 3*m2)*(-1 + xi1)*
     &                  (v1n + v1np1k + v1n*xi1 - v1np1k*xi1)**2*
     &                  Cos((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.))/2. + 
     &               (9*L2*m2*(v2n*(-1 - xi1) - v2np1k*(1 - xi1))**2*
     &                  (-1 + xi1)*
     &                  Cos(th1n + th1np1k - th2n - th2np1k + 
     &                 th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))/4.)/(L2*(8*m1 + 15*m2 - 9*m2*
     &                  Cos(th1n + th1np1k - th2n - th2np1k + th1n*xi1 -
     &                    th1np1k*xi1 - th2n*xi1 + th2np1k*xi1))) + 
     &            ((3*L1*(m1 + 3*m2)*(-1 - xi1)*
     &                  (v1n + v1np1k - v1n*xi1 + v1np1k*xi1)**2*
     &                  Cos((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.))/2. + 
     &               (9*L2*m2*(-1 - xi1)*
     &                  (v2n*(-1 + xi1) - v2np1k*(1 + xi1))**2*
     &                  Cos(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4. - 
     &               (3*g*(m1 + 6*m2)*(1 + xi1)*
     &                  Cos((th2n + th2np1k - th2n*xi1 + th2np1k*xi1)/
     &                   2.))/2. + 
     &               (9*g*(m1 + 2*m2)*(-1 - xi1)*
     &                  Cos(th1n + (th2n*(-1 + xi1))/2. - th1n*xi1 + 
     &                th1np1k*(1 + xi1) - (th2np1k*(1 + xi1))/2.))/2.)/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                 th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)
     & ))-(9*m2*(-1+xi1)*Sin(th1n+th1np1k - th2n - th2np1k + th1n*xi1- 
     &                 th1np1k*xi1 - th2n*xi1 + th2np1k*xi1)*
     &               (9*g*(m1 + 2*m2)*
     &                  Sin(th1n + (th2n*(-1 - xi1))/2. + 
     &                   th1np1k*(1 - xi1) - (th2np1k*(1 - xi1))/2. + 
     &                   th1n*xi1) - 
     &                 3*g*(m1 + 6*m2)*
     &                  Sin((th2n + th2np1k + th2n*xi1 - th2np1k*xi1)/
     &                   2.) + 
     &                 3*L1*(m1 + 3*m2)*
     &                  (v1n + v1np1k + v1n*xi1 - v1np1k*xi1)**2*
     &                  Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &                 (9*L2*m2*(v2n*(-1 - xi1) - v2np1k*(1 - xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))/4.))/(L2*(8*m1 + 15*m2 - 9*m2*
     &                   Cos(th1n + th1np1k - th2n - th2np1k + 
     &                  th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))**2) - (9*m2*(-1 - xi1)*Sin(th1n + th1np1k - th2n - th2np1k - 
     &                 th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)*
     &               (3*L1*(m1 + 3*m2)*
     &                  (v1n + v1np1k - v1n*xi1 + v1np1k*xi1)**2*
     &                  Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & )/2.) + (9*L2*m2*(v2n*(-1 + xi1) - v2np1k*(1 + xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))/4.-3*g*(m1+6*m2)*Sin((th2n+th2np1k - th2n*xi1 + th2np1k*xi1)/
     &                   2.) + 
     &                 9*g*(m1 + 2*m2)*
     &                  Sin(th1n + (th2n*(-1 + xi1))/2. - th1n*xi1 + 
     &                    th1np1k*(1 + xi1) - (th2np1k*(1 + xi1))/2.)))/
     &             (L2*(8*m1 + 15*m2 - 
     &                  9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &               th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)
     & )**2)))/18.))/12.))/
     & (4.*(th2n - th2np1k)*((L1**2*
     &     (-2*g*(m1 + 2*m2)*Cos(th1n) + 2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &        L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &           v1np1k*v2np1k*
     &            (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k))))**2)/
     &   (16.*(th1n - th1np1k)**2) + 
     &  (L2**2*m2**2*(2*L2*(v2n + v2np1k) + 3*L1*v1n*Cos(th1n - th2n) + 
     &        3*L1*v1np1k*Cos(th1np1k - th2np1k))**2)/144. + 
     &  (L1**2*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &        3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &           v2n*Cos(th1np1k - th2np1k)))**2)/144. + 
     &  (L2**2*m2**2*(L1*v1np1k*v2np1k*
     &         (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &      L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k)) + 
     & 2*g*(-Cos(th2n) + Cos(th2np1k)))**2)/(16.*(th2n - th2np1k)**2)))

      m23=(L2*m2*(L1*v1np1k*v2np1k*(Cos(th1n-th2n)-Cos(th1n-th2np1k))+
     &  L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k)) + 
     &  2*g*(-Cos(th2n) + Cos(th2np1k)))*
     & ((L1*L2*m2*v2np1k*(v2n + v2np1k)*
     &   (Cos(th1n - th2n) - Cos(th1n - th2np1k)))/(8.*(th2n - th2np1k))
     &  + (L1*(-2*g*(m1 + 2*m2)*Cos(th1n)+2*g*(m1+2*m2)*Cos(th1np1k)+
     &       L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &          v1np1k*v2np1k*(Cos(th1n - th2np1k) - 
     &             Cos(th1np1k - th2np1k)))))/(8.*(th1n - th1np1k)) + 
     &  (L1*L2*m2*(v1n + v1np1k)*v2np1k*
     &     (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k)))/
     &   (8.*(th1n - th1np1k)) + 
     &  (L2*m2*(2*L2*(v2n + v2np1k) + 3*L1*v1n*Cos(th1n - th2n) + 
     &       3*L1*v1np1k*Cos(th1np1k - th2np1k))*
     &     ((8*L1*(m1 + 3*m2)*(v1n + v1np1k)*
     &          Sin((th1n + th1np1k - th2n - th2np1k)/2.))/
     &        (3.*L2*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((6*L1*(m1 + 3*m2)*(1 + xi1)*
     &               (v1n + v1np1k - v1n*xi1 + v1np1k*xi1)*
     &               Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.))/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))) + (6*L1*(m1 + 3*m2)*(1 - xi1)*
     &               (v1n + v1np1k + v1n*xi1 - v1np1k*xi1)*
     &               Sin((th1n + th1np1k - th2n - th2np1k + th1n*xi1 - 
     &                   th1np1k*xi1 - th2n*xi1 + th2np1k*xi1)/2.))/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                  th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & )))))/18.))/12. + (L1*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &       3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &          v2n*Cos(th1np1k - th2np1k)))*
     &     ((2*m2*(-v1n - v1np1k)*Sin(th1n + th1np1k - th2n - th2np1k))/
     &        (8*m1 + 15*m2 - 9*m2*Cos(th1n + th1np1k - th2n - th2np1k))
     & + (5*((-9*m2*(-1 - xi1)*(v1n*(-1 + xi1) - v1np1k*(1 + xi1))*
     &               Sin(th1n + th1np1k - th2n - th2np1k - th1n*xi1 + 
     &                 th1np1k*xi1 + th2n*xi1 - th2np1k*xi1))/
     &             (2.*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))) - (9*m2*(v1n*(-1 - xi1) - v1np1k*(1 - xi1))*(-1 + xi1)*
     &               Sin(th1n + th1np1k - th2n - th2np1k + th1n*xi1 - 
     &                 th1np1k*xi1 - th2n*xi1 + th2np1k*xi1))/
     &             (2.*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                  th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & )))))/18.))/12. + (L1**2*(m1 + 3*m2)*
     &     ((-4*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k)/2.) + 
     &            m2*(3*g*Sin((th1n + th1np1k - 2*th2n - 2*th2np1k)/
     &                  2.) + 
     &               L2*(v2n + v2np1k)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &               (3*L1*(-v1n - v1np1k)**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k))/4.)))/
     &        (3.*L1*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((-3*(g*(4*m1 + 5*m2)*
     &                  Sin((th1n + th1np1k - th1n*xi1 + th1np1k*xi1)/
     &                   2.) + 
     &                 m2*(3*g*
     &                   Sin((th1n + th1np1k - 2*th2n - 2*th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + 2*th2n*xi1 - 
     &                  2*th2np1k*xi1)/2.) + 
     &                   L2*(v2n + v2np1k - v2n*xi1 + v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                  th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 + xi1) - v1np1k*(1 + xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4.)))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))) - (3*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k + th1n*xi1 - 
     &                   th1np1k*xi1)/2.) + 
     &                 m2*(L2*(v2n + v2np1k + v2n*xi1 - v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 - xi1) - v1np1k*(1 - xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4. + 
     &                   3*g*Sin((th1n + th1np1k - 2*th2n - 
     &                   2*th2np1k + th1n*xi1 - th1np1k*xi1 - 
     &                   2*th2n*xi1 + 2*th2np1k*xi1)/2.))))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                  th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & )))))/18.))/6. + (L1*L2*m2*Cos(th1np1k - th2np1k)*
     &     ((4*(3*L1*(m1 + 3*m2)*(v1n + v1np1k)**2*
     &             Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &            (9*L2*m2*(-v2n - v2np1k)**2*
     &               Sin(th1n + th1np1k - th2n - th2np1k))/4. + 
     &            9*g*(m1 + 2*m2)*
     &             Sin(th1n + th1np1k - th2n/2. - th2np1k/2.) - 
     &            3*g*(m1 + 6*m2)*Sin((th2n + th2np1k)/2.)))/
     &        (9.*L2*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((9*g*(m1 + 2*m2)*
     &                Sin(th1n + (th2n*(-1 - xi1))/2. + 
     &                  th1np1k*(1 - xi1) - (th2np1k*(1 - xi1))/2. + 
     &                  th1n*xi1) - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k + th2n*xi1 - th2np1k*xi1)/2.)
     & + 3*L1*(m1 + 3*m2)*(v1n + v1np1k + v1n*xi1 - v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 - xi1) - v2np1k*(1 - xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4.)/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                  th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))) + (3*L1*(m1 + 3*m2)*(v1n + v1np1k - v1n*xi1 + v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 + xi1) - v2np1k*(1 + xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4. - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k - th2n*xi1 + th2np1k*xi1)/2.)
     & + 9*g*(m1 + 2*m2)*Sin(th1n + (th2n*(-1 + xi1))/2. - th1n*xi1 + 
     &                  th1np1k*(1 + xi1) - (th2np1k*(1 + xi1))/2.))/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                 th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)
     & ))))/18.))/4.))/
     & (4.*(th2n - th2np1k)*((L1**2*
     &     (-2*g*(m1 + 2*m2)*Cos(th1n) + 2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &        L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &           v1np1k*v2np1k*
     &            (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k))))**2)/
     &   (16.*(th1n - th1np1k)**2) + 
     &  (L2**2*m2**2*(2*L2*(v2n + v2np1k) + 3*L1*v1n*Cos(th1n - th2n) + 
     &        3*L1*v1np1k*Cos(th1np1k - th2np1k))**2)/144. + 
     &  (L1**2*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &        3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &           v2n*Cos(th1np1k - th2np1k)))**2)/144. + 
     &  (L2**2*m2**2*(L1*v1np1k*v2np1k*
     &         (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &        L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k-th2np1k))+ 
     &  2*g*(-Cos(th2n) + Cos(th2np1k)))**2)/(16.*(th2n - th2np1k)**2))
     & )-(L2*m2*(L1*v1np1k*v2np1k*(Cos(th1n-th2n)-Cos(th1n - th2np1k))+
     &  L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k)) + 
     &  2*g*(-Cos(th2n) + Cos(th2np1k)))*
     & ((L1**2*L2*m2*v2np1k*(-2*g*(m1 + 2*m2)*Cos(th1n) + 
     &       2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &       L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &          v1np1k*v2np1k*
     &           (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k))))*
     &     (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k)))/
     &   (8.*(th1n - th1np1k)**2) + 
     &  (L1*L2**2*m2**2*Cos(th1np1k - th2np1k)*
     &     (2*L2*(v2n + v2np1k) + 3*L1*v1n*Cos(th1n - th2n) + 
     &       3*L1*v1np1k*Cos(th1np1k - th2np1k)))/24. + 
     &  (L1**3*(m1 + 3*m2)*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &       3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &          v2n*Cos(th1np1k - th2np1k))))/36. + 
     &  (L1*L2**2*m2**2*v2np1k*(Cos(th1n - th2n) - Cos(th1n - th2np1k))*
     &     (L1*v1np1k*v2np1k*(Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &       L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k)) + 
     &      2*g*(-Cos(th2n) + Cos(th2np1k))))/(8.*(th2n - th2np1k)**2))*
     & ((L1*(v1n + v1np1k)*(-2*g*(m1 + 2*m2)*Cos(th1n) + 
     &       2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &       L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &          v1np1k*v2np1k*(Cos(th1n - th2np1k) - 
     &             Cos(th1np1k - th2np1k)))))/(8.*(th1n - th1np1k)) + 
     &  (L2*m2*(v2n + v2np1k)*(L1*v1np1k*v2np1k*
     &        (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &     L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k)) +
     &       2*g*(-Cos(th2n) + Cos(th2np1k))))/(8.*(th2n - th2np1k)) + 
     &  (L1*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &       3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &          v2n*Cos(th1np1k - th2np1k)))*
     &     ((-4*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k)/2.) + 
     &            m2*(3*g*Sin((th1n + th1np1k - 2*th2n - 2*th2np1k)/
     &                  2.) + 
     &               L2*(v2n + v2np1k)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &               (3*L1*(-v1n - v1np1k)**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k))/4.)))/
     &        (3.*L1*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((-3*(g*(4*m1 + 5*m2)*
     &                  Sin((th1n + th1np1k - th1n*xi1 + th1np1k*xi1)/
     &                   2.) + 
     &                 m2*(3*g*
     &                   Sin((th1n + th1np1k - 2*th2n - 2*th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + 2*th2n*xi1 - 
     &                  2*th2np1k*xi1)/2.) + 
     &                   L2*(v2n + v2np1k - v2n*xi1 + v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                  th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 + xi1) - v1np1k*(1 + xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4.)))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))) - (3*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k + th1n*xi1 - 
     &                   th1np1k*xi1)/2.) + 
     &                 m2*(L2*(v2n + v2np1k + v2n*xi1 - v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 - xi1) - v1np1k*(1 - xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4. + 
     &                   3*g*Sin((th1n + th1np1k - 2*th2n - 
     &                   2*th2np1k + th1n*xi1 - th1np1k*xi1 - 
     &                   2*th2n*xi1 + 2*th2np1k*xi1)/2.))))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                 th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & )))))/18.))/12. + (L2*m2*(2*L2*(v2n + v2np1k) + 
     &       3*L1*v1n*Cos(th1n - th2n) + 
     &       3*L1*v1np1k*Cos(th1np1k - th2np1k))*
     &     ((4*(3*L1*(m1 + 3*m2)*(v1n + v1np1k)**2*
     &             Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &            (9*L2*m2*(-v2n - v2np1k)**2*
     &               Sin(th1n + th1np1k - th2n - th2np1k))/4. + 
     &            9*g*(m1 + 2*m2)*
     &             Sin(th1n + th1np1k - th2n/2. - th2np1k/2.) - 
     &            3*g*(m1 + 6*m2)*Sin((th2n + th2np1k)/2.)))/
     &        (9.*L2*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((9*g*(m1 + 2*m2)*
     &                Sin(th1n + (th2n*(-1 - xi1))/2. + 
     &                  th1np1k*(1 - xi1) - (th2np1k*(1 - xi1))/2. + 
     &                  th1n*xi1) - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k + th2n*xi1 - th2np1k*xi1)/2.)
     & + 3*L1*(m1 + 3*m2)*(v1n + v1np1k + v1n*xi1 - v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 - xi1) - v2np1k*(1 - xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4.)/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                  th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))) + (3*L1*(m1 + 3*m2)*(v1n + v1np1k - v1n*xi1 + v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 + xi1) - v2np1k*(1 + xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4. - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k - th2n*xi1 + th2np1k*xi1)/2.)
     & + 9*g*(m1 + 2*m2)*Sin(th1n + (th2n*(-1 + xi1))/2. - th1n*xi1 + 
     &                  th1np1k*(1 + xi1) - (th2np1k*(1 + xi1))/2.))/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)
     & ))))/18.))/12.))/
     &    (4.*(th2n - th2np1k)*((L1**2*
     &      (-2*g*(m1 + 2*m2)*Cos(th1n) + 2*g*(m1 + 2*m2)*Cos(th1np1k) +
     &         L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) +
     &            v1np1k*v2np1k*
     &             (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k))))**2)/
     &    (16.*(th1n - th1np1k)**2) + 
     &   (L2**2*m2**2*(2*L2*(v2n + v2np1k) + 3*L1*v1n*Cos(th1n - th2n) +
     &         3*L1*v1np1k*Cos(th1np1k - th2np1k))**2)/144. + 
     &   (L1**2*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &         3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &            v2n*Cos(th1np1k - th2np1k)))**2)/144. + 
     &   (L2**2*m2**2*(L1*v1np1k*v2np1k*
     &          (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &         L1*v1n*v2n*(Cos(th1np1k - th2n) - 
     &            Cos(th1np1k - th2np1k)) + 
     &         2*g*(-Cos(th2n) + Cos(th2np1k)))**2)/
     &    (16.*(th2n - th2np1k)**2))**2) + 
     &   (L1*L2*m2*v2np1k*(Cos(th1n - th2n) - Cos(th1n - th2np1k))*
     & ((L1*(v1n + v1np1k)*(-2*g*(m1 + 2*m2)*Cos(th1n) + 
     &       2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &       L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &          v1np1k*v2np1k*(Cos(th1n - th2np1k) - 
     &             Cos(th1np1k - th2np1k)))))/(8.*(th1n - th1np1k)) + 
     &  (L2*m2*(v2n + v2np1k)*(L1*v1np1k*v2np1k*
     &        (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &       L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k)) +
     &       2*g*(-Cos(th2n) + Cos(th2np1k))))/(8.*(th2n - th2np1k)) + 
     &  (L1*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &       3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &          v2n*Cos(th1np1k - th2np1k)))*
     &     ((-4*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k)/2.) + 
     &            m2*(3*g*Sin((th1n + th1np1k - 2*th2n - 2*th2np1k)/
     &                  2.) + 
     &               L2*(v2n + v2np1k)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &               (3*L1*(-v1n - v1np1k)**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k))/4.)))/
     &        (3.*L1*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((-3*(g*(4*m1 + 5*m2)*
     &                  Sin((th1n + th1np1k - th1n*xi1 + th1np1k*xi1)/
     &                   2.) + 
     &                 m2*(3*g*
     &                   Sin((th1n + th1np1k - 2*th2n - 2*th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + 2*th2n*xi1 - 
     &                   2*th2np1k*xi1)/2.) + 
     &                   L2*(v2n + v2np1k - v2n*xi1 + v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 + xi1) - v1np1k*(1 + xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4.)))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                 th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))) - (3*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k + th1n*xi1 - 
     &                   th1np1k*xi1)/2.) + 
     &                 m2*(L2*(v2n + v2np1k + v2n*xi1 - v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &                    (3*L1*(v1n*(-1 - xi1) - v1np1k*(1 - xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4. + 
     &                    3*g*
     &                    Sin((th1n + th1np1k - 2*th2n - 2*th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - 2*th2n*xi1 + 
     &                   2*th2np1k*xi1)/2.))))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                 th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1)
     & ))))/18.))/12. + (L2*m2*(2*L2*(v2n + v2np1k) + 
     &  3*L1*v1n*Cos(th1n - th2n) + 3*L1*v1np1k*Cos(th1np1k - th2np1k)
     & )*((4*(3*L1*(m1 + 3*m2)*(v1n + v1np1k)**2*
     &             Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &            (9*L2*m2*(-v2n - v2np1k)**2*
     &               Sin(th1n + th1np1k - th2n - th2np1k))/4. + 
     &            9*g*(m1 + 2*m2)*
     &             Sin(th1n + th1np1k - th2n/2. - th2np1k/2.) - 
     &            3*g*(m1 + 6*m2)*Sin((th2n + th2np1k)/2.)))/
     &        (9.*L2*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((9*g*(m1 + 2*m2)*
     &                Sin(th1n + (th2n*(-1 - xi1))/2. + 
     &                  th1np1k*(1 - xi1) - (th2np1k*(1 - xi1))/2. + 
     &                  th1n*xi1) - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k + th2n*xi1 - th2np1k*xi1)/2.)
     & + 3*L1*(m1 + 3*m2)*(v1n + v1np1k + v1n*xi1 - v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 - xi1) - v2np1k*(1 - xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))/4.)/(L2*(8*m1 + 15*m2 - 9*m2*
     &                  Cos(th1n + th1np1k - th2n - th2np1k + th1n*xi1 - 
     &                    th1np1k*xi1 - th2n*xi1 + th2np1k*xi1))) + 
     &            (3*L1*(m1 + 3*m2)*
     &                (v1n + v1np1k - v1n*xi1 + v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k - th1n*xi1 + 
     &                   th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 + xi1) - v2np1k*(1 + xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))/4.-3*g*(m1+6*m2)*Sin((th2n + th2np1k - th2n*xi1+th2np1k*xi1)/
     &                  2.) + 9*g*(m1 + 2*m2)*
     &                Sin(th1n + (th2n*(-1 + xi1))/2. - th1n*xi1 + 
     &                  th1np1k*(1 + xi1) - (th2np1k*(1 + xi1))/2.))/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &              th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1))
     & )))/18.))/12.))/
     &    (4.*(th2n - th2np1k)*((L1**2*
     &     (-2*g*(m1 + 2*m2)*Cos(th1n) + 2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &        L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &           v1np1k*v2np1k*
     &            (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k))))**2)/
     &   (16.*(th1n - th1np1k)**2) + 
     &  (L2**2*m2**2*(2*L2*(v2n + v2np1k) + 3*L1*v1n*Cos(th1n - th2n) + 
     &        3*L1*v1np1k*Cos(th1np1k - th2np1k))**2)/144. + 
     &  (L1**2*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &        3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &           v2n*Cos(th1np1k - th2np1k)))**2)/144. + 
     &  (L2**2*m2**2*(L1*v1np1k*v2np1k*
     &         (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &   L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k)) + 
     &  2*g*(-Cos(th2n) + Cos(th2np1k)))**2)/(16.*(th2n - th2np1k)**2)))

      m24= -0.5 + (L2*m2*(L1*v1np1k*v2np1k*
     &   (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &  L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k)) + 
     &  2*g*(-Cos(th2n) + Cos(th2np1k)))*
     & ((L1*L2*m2*v1np1k*(v2n + v2np1k)*
     &   (Cos(th1n - th2n) - Cos(th1n - th2np1k)))/(8.*(th2n - th2np1k))
     & + (L1*L2*m2*v1np1k*(v1n + v1np1k)*
     &     (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k)))/
     &   (8.*(th1n - th1np1k)) + 
     &  (L2*m2*(L1*v1np1k*v2np1k*
     &        (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &       L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k)) +
     &       2*g*(-Cos(th2n) + Cos(th2np1k))))/(8.*(th2n - th2np1k)) + 
     &  (L1*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &       3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &          v2n*Cos(th1np1k - th2np1k)))*
     &     ((-8*L2*m2*(v2n + v2np1k)*
     &          Sin((th1n + th1np1k - th2n - th2np1k)/2.))/
     &        (3.*L1*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((-6*L2*m2*(1 + xi1)*
     &               (v2n + v2np1k - v2n*xi1 + v2np1k*xi1)*
     &               Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))) - (6*L2*m2*(1 - xi1)*(v2n + v2np1k + v2n*xi1 - v2np1k*xi1)*
     &               Sin((th1n + th1np1k - th2n - th2np1k + th1n*xi1 - 
     &                   th1np1k*xi1 - th2n*xi1 + th2np1k*xi1)/2.))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                  th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & )))))/18.))/12. + (L2*m2*(2*L2*(v2n + v2np1k) + 
     &       3*L1*v1n*Cos(th1n - th2n) + 
     &       3*L1*v1np1k*Cos(th1np1k - th2np1k))*
     &   ((-2*m2*(-v2n - v2np1k)*Sin(th1n + th1np1k - th2n - th2np1k))/
     &        (8*m1 + 15*m2 - 9*m2*Cos(th1n + th1np1k - th2n - th2np1k))
     & + (5*((9*m2*(-1 - xi1)*(v2n*(-1 + xi1) - v2np1k*(1 + xi1))*
     &               Sin(th1n + th1np1k - th2n - th2np1k - th1n*xi1 + 
     &                 th1np1k*xi1 + th2n*xi1 - th2np1k*xi1))/
     &             (2.*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))) + (9*m2*(v2n*(-1 - xi1) - v2np1k*(1 - xi1))*(-1 + xi1)*
     &               Sin(th1n + th1np1k - th2n - th2np1k + th1n*xi1 - 
     &                 th1np1k*xi1 - th2n*xi1 + th2np1k*xi1))/
     &             (2.*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                  th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & )))))/18.))/12. + (L1*L2*m2*Cos(th1n - th2n)*
     &     ((-4*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k)/2.) + 
     &            m2*(3*g*Sin((th1n + th1np1k - 2*th2n - 2*th2np1k)/
     &                  2.) + 
     &               L2*(v2n + v2np1k)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &               (3*L1*(-v1n - v1np1k)**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k))/4.)))/
     &        (3.*L1*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((-3*(g*(4*m1 + 5*m2)*
     &                  Sin((th1n + th1np1k - th1n*xi1 + th1np1k*xi1)/
     &                   2.) + 
     &                 m2*(3*g*
     &                   Sin((th1n + th1np1k - 2*th2n - 2*th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + 2*th2n*xi1 - 
     &                  2*th2np1k*xi1)/2.) + 
     &                   L2*(v2n + v2np1k - v2n*xi1 + v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                  th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 + xi1) - v1np1k*(1 + xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4.)))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))) - (3*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k + th1n*xi1 - 
     &                   th1np1k*xi1)/2.) + 
     &                 m2*(L2*(v2n + v2np1k + v2n*xi1 - v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 - xi1) - v1np1k*(1 - xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4. + 
     &                   3*g*Sin((th1n + th1np1k - 2*th2n - 
     &                   2*th2np1k + th1n*xi1 - th1np1k*xi1 - 
     &                   2*th2n*xi1 + 2*th2np1k*xi1)/2.))))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                 th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & )))))/18.))/4. + (L2**2*m2*((4*
     &          (3*L1*(m1 + 3*m2)*(v1n + v1np1k)**2*
     &             Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &            (9*L2*m2*(-v2n - v2np1k)**2*
     &               Sin(th1n + th1np1k - th2n - th2np1k))/4. + 
     &            9*g*(m1 + 2*m2)*
     &             Sin(th1n + th1np1k - th2n/2. - th2np1k/2.) - 
     &            3*g*(m1 + 6*m2)*Sin((th2n + th2np1k)/2.)))/
     &        (9.*L2*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((9*g*(m1 + 2*m2)*
     &                Sin(th1n + (th2n*(-1 - xi1))/2. + 
     &                  th1np1k*(1 - xi1) - (th2np1k*(1 - xi1))/2. + 
     &                  th1n*xi1) - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k + th2n*xi1 - th2np1k*xi1)/2.)
     & + 3*L1*(m1 + 3*m2)*(v1n + v1np1k + v1n*xi1 - v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 - xi1) - v2np1k*(1 - xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4.)/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                  th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))) + (3*L1*(m1 + 3*m2)*(v1n + v1np1k - v1n*xi1 + v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 + xi1) - v2np1k*(1 + xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4. - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k - th2n*xi1 + th2np1k*xi1)/2.)
     & + 9*g*(m1 + 2*m2)*Sin(th1n + (th2n*(-1 + xi1))/2. - th1n*xi1 + 
     &                  th1np1k*(1 + xi1) - (th2np1k*(1 + xi1))/2.))/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                 th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)
     & ))))/18.))/6.))/
     & (4.*(th2n - th2np1k)*((L1**2*
     &     (-2*g*(m1 + 2*m2)*Cos(th1n) + 2*g*(m1 + 2*m2)*Cos(th1np1k) +
     &        L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) +
     &           v1np1k*v2np1k*
     &            (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k))))**2)/
     &   (16.*(th1n - th1np1k)**2) + 
     &  (L2**2*m2**2*(2*L2*(v2n + v2np1k) + 3*L1*v1n*Cos(th1n - th2n) + 
     &        3*L1*v1np1k*Cos(th1np1k - th2np1k))**2)/144. + 
     &  (L1**2*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &        3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &           v2n*Cos(th1np1k - th2np1k)))**2)/144. + 
     &  (L2**2*m2**2*(L1*v1np1k*v2np1k*
     &         (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &        L1*v1n*v2n*(Cos(th1np1k - th2n)-Cos(th1np1k-th2np1k)) + 
     &  2*g*(-Cos(th2n) + Cos(th2np1k)))**2)/(16.*(th2n - th2np1k)**2))
     & )-(L2*m2*(L1*v1np1k*v2np1k*(Cos(th1n-th2n)-Cos(th1n - th2np1k)) +
     &  L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k)) + 
     &  2*g*(-Cos(th2n) + Cos(th2np1k)))*
     & ((L1**2*L2*m2*v1np1k*(-2*g*(m1 + 2*m2)*Cos(th1n) + 
     &       2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &       L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &          v1np1k*v2np1k*
     &           (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k))))*
     &     (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k)))/
     &   (8.*(th1n - th1np1k)**2) + 
     &  (L2**3*m2**2*(2*L2*(v2n + v2np1k) + 3*L1*v1n*Cos(th1n - th2n) + 
     &       3*L1*v1np1k*Cos(th1np1k - th2np1k)))/36. + 
     &  (L1**2*L2*m2*Cos(th1n - th2n)*
     &     (2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &       3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &          v2n*Cos(th1np1k - th2np1k))))/24. + 
     &  (L1*L2**2*m2**2*v1np1k*(Cos(th1n - th2n) - Cos(th1n - th2np1k))*
     &     (L1*v1np1k*v2np1k*(Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &       L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k)) + 
     &      2*g*(-Cos(th2n) + Cos(th2np1k))))/(8.*(th2n - th2np1k)**2))*
     & ((L1*(v1n + v1np1k)*(-2*g*(m1 + 2*m2)*Cos(th1n) + 
     &       2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &       L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &          v1np1k*v2np1k*(Cos(th1n - th2np1k) - 
     &             Cos(th1np1k - th2np1k)))))/(8.*(th1n - th1np1k)) + 
     &  (L2*m2*(v2n + v2np1k)*(L1*v1np1k*v2np1k*
     &        (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &       L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k)) +
     &       2*g*(-Cos(th2n) + Cos(th2np1k))))/(8.*(th2n - th2np1k)) + 
     &  (L1*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &       3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &          v2n*Cos(th1np1k - th2np1k)))*
     &     ((-4*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k)/2.) + 
     &            m2*(3*g*Sin((th1n + th1np1k - 2*th2n - 2*th2np1k)/
     &                  2.) + 
     &               L2*(v2n + v2np1k)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &               (3*L1*(-v1n - v1np1k)**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k))/4.)))/
     &        (3.*L1*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((-3*(g*(4*m1 + 5*m2)*
     &                  Sin((th1n + th1np1k - th1n*xi1 + th1np1k*xi1)/
     &                   2.) + 
     &                 m2*(3*g*
     &                   Sin((th1n + th1np1k - 2*th2n - 2*th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + 2*th2n*xi1 - 
     &                  2*th2np1k*xi1)/2.) + 
     &                   L2*(v2n + v2np1k - v2n*xi1 + v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                  th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 + xi1) - v1np1k*(1 + xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4.)))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))) - (3*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k + th1n*xi1 - 
     &                   th1np1k*xi1)/2.) + 
     &                 m2*(L2*(v2n + v2np1k + v2n*xi1 - v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 - xi1) - v1np1k*(1 - xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4. + 
     &                   3*g*Sin((th1n + th1np1k - 2*th2n - 
     &                   2*th2np1k + th1n*xi1 - th1np1k*xi1 - 
     &                   2*th2n*xi1 + 2*th2np1k*xi1)/2.))))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                  th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & )))))/18.))/12. + (L2*m2*(2*L2*(v2n + v2np1k) + 
     &       3*L1*v1n*Cos(th1n - th2n) + 
     &       3*L1*v1np1k*Cos(th1np1k - th2np1k))*
     &     ((4*(3*L1*(m1 + 3*m2)*(v1n + v1np1k)**2*
     &             Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &            (9*L2*m2*(-v2n - v2np1k)**2*
     &               Sin(th1n + th1np1k - th2n - th2np1k))/4. + 
     &            9*g*(m1 + 2*m2)*
     &             Sin(th1n + th1np1k - th2n/2. - th2np1k/2.) - 
     &            3*g*(m1 + 6*m2)*Sin((th2n + th2np1k)/2.)))/
     &        (9.*L2*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((9*g*(m1 + 2*m2)*
     &                Sin(th1n + (th2n*(-1 - xi1))/2. + 
     &                  th1np1k*(1 - xi1) - (th2np1k*(1 - xi1))/2. + 
     &                  th1n*xi1) - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k + th2n*xi1 - th2np1k*xi1)/2.)
     & + 3*L1*(m1 + 3*m2)*(v1n + v1np1k + v1n*xi1 - v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 - xi1) - v2np1k*(1 - xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4.)/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                  th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))) + (3*L1*(m1 + 3*m2)*(v1n + v1np1k - v1n*xi1 + v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 + xi1) - v2np1k*(1 + xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4. - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k - th2n*xi1 + th2np1k*xi1)/2.)
     & + 9*g*(m1 + 2*m2)*Sin(th1n + (th2n*(-1 + xi1))/2. - th1n*xi1 + 
     &                  th1np1k*(1 + xi1) - (th2np1k*(1 + xi1))/2.))/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                 th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)
     & ))))/18.))/12.))/
     & (4.*(th2n - th2np1k)*((L1**2*
     &      (-2*g*(m1 + 2*m2)*Cos(th1n) + 2*g*(m1 + 2*m2)*Cos(th1np1k) +
     &         L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) +
     &            v1np1k*v2np1k*
     &             (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k))))**2)/
     &    (16.*(th1n - th1np1k)**2) + 
     &   (L2**2*m2**2*(2*L2*(v2n + v2np1k) + 3*L1*v1n*Cos(th1n - th2n) +
     &         3*L1*v1np1k*Cos(th1np1k - th2np1k))**2)/144. + 
     &   (L1**2*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &         3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &            v2n*Cos(th1np1k - th2np1k)))**2)/144. + 
     &   (L2**2*m2**2*(L1*v1np1k*v2np1k*
     &          (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &         L1*v1n*v2n*(Cos(th1np1k - th2n) - 
     &            Cos(th1np1k - th2np1k)) + 
     &         2*g*(-Cos(th2n) + Cos(th2np1k)))**2)/
     &    (16.*(th2n - th2np1k)**2))**2) + 
     & (L1*L2*m2*v1np1k*(Cos(th1n - th2n) - Cos(th1n - th2np1k))*
     & ((L1*(v1n + v1np1k)*(-2*g*(m1 + 2*m2)*Cos(th1n) + 
     &       2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &       L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &          v1np1k*v2np1k*(Cos(th1n - th2np1k) - 
     &             Cos(th1np1k - th2np1k)))))/(8.*(th1n - th1np1k)) + 
     &  (L2*m2*(v2n + v2np1k)*(L1*v1np1k*v2np1k*
     &        (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &       L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k)) +
     &       2*g*(-Cos(th2n) + Cos(th2np1k))))/(8.*(th2n - th2np1k)) + 
     &  (L1*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &       3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &          v2n*Cos(th1np1k - th2np1k)))*
     &     ((-4*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k)/2.) + 
     &            m2*(3*g*Sin((th1n + th1np1k - 2*th2n - 2*th2np1k)/
     &                  2.) + 
     &               L2*(v2n + v2np1k)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &               (3*L1*(-v1n - v1np1k)**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k))/4.)))/
     &        (3.*L1*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((-3*(g*(4*m1 + 5*m2)*
     &                  Sin((th1n + th1np1k - th1n*xi1 + th1np1k*xi1)/
     &                   2.) + 
     &                 m2*(3*g*
     &                   Sin((th1n + th1np1k - 2*th2n - 2*th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + 2*th2n*xi1 - 
     &                   2*th2np1k*xi1)/2.) + 
     &                   L2*(v2n + v2np1k - v2n*xi1 + v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 + xi1) - v1np1k*(1 + xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4.)))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                 th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))) - (3*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k + th1n*xi1 - 
     &                   th1np1k*xi1)/2.) + 
     &                 m2*(L2*(v2n + v2np1k + v2n*xi1 - v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &                    (3*L1*(v1n*(-1 - xi1) - v1np1k*(1 - xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4. + 
     &                    3*g*
     &                    Sin((th1n + th1np1k - 2*th2n - 2*th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - 2*th2n*xi1 + 
     &                   2*th2np1k*xi1)/2.))))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                 th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1)
     & ))))/18.))/12. + (L2*m2*(2*L2*(v2n + v2np1k) + 
     &   3*L1*v1n*Cos(th1n - th2n) + 3*L1*v1np1k*Cos(th1np1k - th2np1k)
     & )*((4*(3*L1*(m1 + 3*m2)*(v1n + v1np1k)**2*
     &             Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &            (9*L2*m2*(-v2n - v2np1k)**2*
     &               Sin(th1n + th1np1k - th2n - th2np1k))/4. + 
     &            9*g*(m1 + 2*m2)*
     &             Sin(th1n + th1np1k - th2n/2. - th2np1k/2.) - 
     &            3*g*(m1 + 6*m2)*Sin((th2n + th2np1k)/2.)))/
     &        (9.*L2*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((9*g*(m1 + 2*m2)*
     &                Sin(th1n + (th2n*(-1 - xi1))/2. + 
     &                  th1np1k*(1 - xi1) - (th2np1k*(1 - xi1))/2. + 
     &                  th1n*xi1) - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k + th2n*xi1 - th2np1k*xi1)/2.)
     & + 3*L1*(m1 + 3*m2)*(v1n + v1np1k + v1n*xi1 - v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 - xi1) - v2np1k*(1 - xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))/4.)/(L2*(8*m1 + 15*m2 - 9*m2*
     &                  Cos(th1n + th1np1k - th2n - th2np1k + th1n*xi1 - 
     &                    th1np1k*xi1 - th2n*xi1 + th2np1k*xi1))) + 
     &            (3*L1*(m1 + 3*m2)*
     &                (v1n + v1np1k - v1n*xi1 + v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k - th1n*xi1 + 
     &                   th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 + xi1) - v2np1k*(1 + xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))/4.-3*g*(m1+6*m2)*Sin((th2n + th2np1k - th2n*xi1+th2np1k*xi1)/
     &                  2.) + 9*g*(m1 + 2*m2)*
     &                Sin(th1n + (th2n*(-1 + xi1))/2. - th1n*xi1 + 
     &                  th1np1k*(1 + xi1) - (th2np1k*(1 + xi1))/2.))/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &              th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1))
     & )))/18.))/12.))/
     & (4.*(th2n - th2np1k)*((L1**2*
     &     (-2*g*(m1 + 2*m2)*Cos(th1n) + 2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &        L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &           v1np1k*v2np1k*
     &            (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k))))**2)/
     &   (16.*(th1n - th1np1k)**2) + 
     &  (L2**2*m2**2*(2*L2*(v2n + v2np1k) + 3*L1*v1n*Cos(th1n - th2n) + 
     &        3*L1*v1np1k*Cos(th1np1k - th2np1k))**2)/144. + 
     &  (L1**2*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &        3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &           v2n*Cos(th1np1k - th2np1k)))**2)/144. + 
     &  (L2**2*m2**2*(L1*v1np1k*v2np1k*
     &         (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &  L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k)) + 
     &  2*g*(-Cos(th2n) + Cos(th2np1k)))**2)/(16.*(th2n - th2np1k)**2)))

      m31=(4*((g*(4*m1 + 5*m2)*Cos((th1n + th1np1k)/2.))/2. + 
     &  m2*((3*g*Cos((th1n + th1np1k - 2*th2n - 2*th2np1k)/2.))/2. + 
     &     (L2*(v2n + v2np1k)**2*
     &        Cos((th1n + th1np1k - th2n - th2np1k)/2.))/2. + 
     &     (3*L1*(-v1n - v1np1k)**2*Cos(th1n + th1np1k-th2n-th2np1k))/4.
     & )))/(3.*L1*(8*m1+15*m2-9*m2*Cos(th1n + th1np1k-th2n-th2np1k)))- 
     &    (12*m2*Sin(th1n + th1np1k - th2n - th2np1k)*
     & (g*(4*m1 + 5*m2)*Sin((th1n + th1np1k)/2.) + 
     &  m2*(3*g*Sin((th1n + th1np1k - 2*th2n - 2*th2np1k)/2.) + 
     &  L2*(v2n + v2np1k)**2*Sin((th1n + th1np1k - th2n - th2np1k)/2.) +
     &(3*L1*(-v1n - v1np1k)**2*Sin(th1n + th1np1k - th2n - th2np1k))/4.
     & )))/(L1*(8*m1+15*m2-9*m2*Cos(th1n+th1np1k-th2n-th2np1k))**2) -
     &   (5*((-3*((g*(4*m1 + 5*m2)*(1 + xi1)*
     &          Cos((th1n + th1np1k - th1n*xi1 + th1np1k*xi1)/2.))/2. + 
     &       m2*((3*g*(1 + xi1)*
     &             Cos((th1n + th1np1k - 2*th2n - 2*th2np1k - 
     &                 th1n*xi1 + th1np1k*xi1 + 2*th2n*xi1 - 
     &                 2*th2np1k*xi1)/2.))/2. + 
     &          (L2*(1 + xi1)*(v2n + v2np1k - v2n*xi1 + v2np1k*xi1)**2*
     &             Cos((th1n + th1np1k - th2n - th2np1k - th1n*xi1 + 
     &                 th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)/2.))/2. + 
     &          (3*L1*(1 + xi1)*(v1n*(-1 + xi1) - v1np1k*(1 + xi1))**2*
     &             Cos(th1n + th1np1k - th2n - th2np1k - th1n*xi1 + 
     &               th1np1k*xi1 + th2n*xi1 - th2np1k*xi1))/4.)))/
     &   (L1*(8*m1 + 15*m2 - 9*m2*
     &        Cos(th1n + th1np1k - th2n - th2np1k - th1n*xi1 + 
     &          th1np1k*xi1 + th2n*xi1 - th2np1k*xi1))) - 
     &  (3*((g*(4*m1 + 5*m2)*(1 - xi1)*
     &          Cos((th1n + th1np1k + th1n*xi1 - th1np1k*xi1)/2.))/2. + 
     &       m2*((L2*(1 - xi1)*(v2n + v2np1k + v2n*xi1 - v2np1k*xi1)**2*
     &             Cos((th1n + th1np1k - th2n - th2np1k + th1n*xi1 - 
     &                 th1np1k*xi1 - th2n*xi1 + th2np1k*xi1)/2.))/2. + 
     &          (3*L1*(v1n*(-1 - xi1) - v1np1k*(1 - xi1))**2*(1 - xi1)*
     &             Cos(th1n + th1np1k - th2n - th2np1k + th1n*xi1 - 
     &               th1np1k*xi1 - th2n*xi1 + th2np1k*xi1))/4. + 
     &          (3*g*(1 - xi1)*
     &             Cos((th1n + th1np1k - 2*th2n - 2*th2np1k + 
     &                 th1n*xi1 - th1np1k*xi1 - 2*th2n*xi1 + 
     &                 2*th2np1k*xi1)/2.))/2.)))/
     &   (L1*(8*m1 + 15*m2 - 9*m2*
     &        Cos(th1n + th1np1k - th2n - th2np1k + th1n*xi1 - 
     &          th1np1k*xi1 - th2n*xi1 + th2np1k*xi1))) + 
     &  (27*m2*(1 + xi1)*Sin(th1n + th1np1k - th2n - th2np1k - 
     &       th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)*
     &     (g*(4*m1 + 5*m2)*Sin((th1n + th1np1k - th1n*xi1 + 
     &            th1np1k*xi1)/2.) + 
     &       m2*(3*g*Sin((th1n + th1np1k - 2*th2n - 2*th2np1k - 
     &               th1n*xi1 + th1np1k*xi1 + 2*th2n*xi1 - 2*th2np1k*xi1
     & )/2.) + L2*(v2n + v2np1k - v2n*xi1 + v2np1k*xi1)**2*
     &           Sin((th1n + th1np1k - th2n - th2np1k - th1n*xi1 + 
     &               th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)/2.) + 
     &          (3*L1*(v1n*(-1 + xi1) - v1np1k*(1 + xi1))**2*
     &             Sin(th1n + th1np1k - th2n - th2np1k - th1n*xi1 + 
     &               th1np1k*xi1 + th2n*xi1 - th2np1k*xi1))/4.)))/
     &   (L1*(8*m1 + 15*m2 - 9*m2*
     &         Cos(th1n + th1np1k - th2n - th2np1k - th1n*xi1 + 
     &           th1np1k*xi1 + th2n*xi1 - th2np1k*xi1))**2) + 
     &  (27*m2*(1 - xi1)*Sin(th1n + th1np1k - th2n - th2np1k + th1n*xi1-
     &       th1np1k*xi1 - th2n*xi1 + th2np1k*xi1)*
     &  (g*(4*m1 + 5*m2)*Sin((th1n + th1np1k + th1n*xi1 - th1np1k*xi1)/
     &          2.) + m2*(L2*(v2n + v2np1k + v2n*xi1 - v2np1k*xi1)**2*
     &           Sin((th1n + th1np1k - th2n - th2np1k + th1n*xi1 - 
     &               th1np1k*xi1 - th2n*xi1 + th2np1k*xi1)/2.) + 
     &          (3*L1*(v1n*(-1 - xi1) - v1np1k*(1 - xi1))**2*
     &             Sin(th1n + th1np1k - th2n - th2np1k + th1n*xi1 - 
     &               th1np1k*xi1 - th2n*xi1 + th2np1k*xi1))/4. + 
     &          3*g*Sin((th1n + th1np1k - 2*th2n - 2*th2np1k + th1n*xi1-
     &               th1np1k*xi1 - 2*th2n*xi1 + 2*th2np1k*xi1)/2.))))/
     &   (L1*(8*m1 + 15*m2 - 9*m2*
     &         Cos(th1n + th1np1k - th2n - th2np1k + th1n*xi1 - 
     &           th1np1k*xi1 - th2n*xi1 + th2np1k*xi1))**2)))/18. - 
     &   (L1*L2*m2*v2n*Sin(th1np1k - th2np1k)*
     & ((L1*(v1n + v1np1k)*(-2*g*(m1 + 2*m2)*Cos(th1n) + 
     &       2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &       L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &          v1np1k*v2np1k*(Cos(th1n - th2np1k) - 
     &             Cos(th1np1k - th2np1k)))))/(8.*(th1n - th1np1k)) + 
     &  (L2*m2*(v2n + v2np1k)*(L1*v1np1k*v2np1k*
     &        (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &       L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k)) +
     &       2*g*(-Cos(th2n) + Cos(th2np1k))))/(8.*(th2n - th2np1k)) + 
     &  (L1*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &       3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &          v2n*Cos(th1np1k - th2np1k)))*
     &     ((-4*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k)/2.) + 
     &            m2*(3*g*Sin((th1n + th1np1k - 2*th2n - 2*th2np1k)/
     &                  2.) + 
     &               L2*(v2n + v2np1k)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &               (3*L1*(-v1n - v1np1k)**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k))/4.)))/
     &        (3.*L1*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((-3*(g*(4*m1 + 5*m2)*
     &                  Sin((th1n + th1np1k - th1n*xi1 + th1np1k*xi1)/
     &                   2.) + 
     &                 m2*(3*g*
     &                   Sin((th1n + th1np1k - 2*th2n - 2*th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + 2*th2n*xi1 - 
     &                  2*th2np1k*xi1)/2.) + 
     &                   L2*(v2n + v2np1k - v2n*xi1 + v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                  th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 + xi1) - v1np1k*(1 + xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4.)))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))) - (3*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k + th1n*xi1 - 
     &                   th1np1k*xi1)/2.) + 
     &                 m2*(L2*(v2n + v2np1k + v2n*xi1 - v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 - xi1) - v1np1k*(1 - xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4. + 
     &                   3*g*Sin((th1n + th1np1k - 2*th2n - 
     &                   2*th2np1k + th1n*xi1 - th1np1k*xi1 - 
     &                   2*th2n*xi1 + 2*th2np1k*xi1)/2.))))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                 th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & )))))/18.))/12. + (L2*m2*(2*L2*(v2n + v2np1k) + 
     &       3*L1*v1n*Cos(th1n - th2n) + 
     &       3*L1*v1np1k*Cos(th1np1k - th2np1k))*
     &     ((4*(3*L1*(m1 + 3*m2)*(v1n + v1np1k)**2*
     &             Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &            (9*L2*m2*(-v2n - v2np1k)**2*
     &               Sin(th1n + th1np1k - th2n - th2np1k))/4. + 
     &            9*g*(m1 + 2*m2)*
     &             Sin(th1n + th1np1k - th2n/2. - th2np1k/2.) - 
     &            3*g*(m1 + 6*m2)*Sin((th2n + th2np1k)/2.)))/
     &        (9.*L2*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((9*g*(m1 + 2*m2)*
     &                Sin(th1n + (th2n*(-1 - xi1))/2. + 
     &                  th1np1k*(1 - xi1) - (th2np1k*(1 - xi1))/2. + 
     &                  th1n*xi1) - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k + th2n*xi1 - th2np1k*xi1)/2.)
     & + 3*L1*(m1 + 3*m2)*(v1n + v1np1k + v1n*xi1 - v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 - xi1) - v2np1k*(1 - xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4.)/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                 th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))) + (3*L1*(m1 + 3*m2)*(v1n + v1np1k - v1n*xi1 + v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 + xi1) - v2np1k*(1 + xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4. - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k - th2n*xi1 + th2np1k*xi1)/2.)
     & + 9*g*(m1 + 2*m2)*Sin(th1n + (th2n*(-1 + xi1))/2. - th1n*xi1 + 
     &                  th1np1k*(1 + xi1) - (th2np1k*(1 + xi1))/2.))/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)
     & ))))/18.))/12.))/
     &    (4.*((L1**2*(-2*g*(m1 + 2*m2)*Cos(th1n) + 
     &        2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &        L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &           v1np1k*v2np1k*
     &            (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k))))**2)/
     &   (16.*(th1n - th1np1k)**2) + 
     &  (L2**2*m2**2*(2*L2*(v2n + v2np1k) + 3*L1*v1n*Cos(th1n - th2n) + 
     &        3*L1*v1np1k*Cos(th1np1k - th2np1k))**2)/144. + 
     &  (L1**2*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &        3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &           v2n*Cos(th1np1k - th2np1k)))**2)/144. + 
     &  (L2**2*m2**2*(L1*v1np1k*v2np1k*
     &         (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &        L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k))+
     &  2*g*(-Cos(th2n) + Cos(th2np1k)))**2)/(16.*(th2n - th2np1k)**2))
     & ) - (L1*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &  3*L2*m2*(v2np1k*Cos(th1n - th2n) + v2n*Cos(th1np1k - th2np1k)))*
     & ((L1**2*(-2*g*(m1 + 2*m2)*Cos(th1n) + 
     &        2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &        L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &           v1np1k*v2np1k*
     &            (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k))))**2)/
     &   (8.*(th1n - th1np1k)**3) - 
     &  (L1*L2**2*m2**2*v1np1k*
     &     (2*L2*(v2n + v2np1k) + 3*L1*v1n*Cos(th1n - th2n) + 
     &     3*L1*v1np1k*Cos(th1np1k - th2np1k))*Sin(th1np1k - th2np1k))/
     &   24. - (L1**2*L2*m2*v2n*
     &     (2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &       3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &      v2n*Cos(th1np1k - th2np1k)))*Sin(th1np1k - th2np1k))/24. + 
     &  (L1*L2**2*m2**2*v1n*v2n*
     &     (L1*v1np1k*v2np1k*(Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &       L1*v1n*v2n*(Cos(th1np1k - th2n) - 
     &          Cos(th1np1k - th2np1k)) + 
     &       2*g*(-Cos(th2n) + Cos(th2np1k)))*
     &     (-Sin(th1np1k - th2n) + Sin(th1np1k - th2np1k)))/
     &   (8.*(th2n - th2np1k)**2) + 
     &  (L1**2*(-2*g*(m1 + 2*m2)*Cos(th1n) + 
     &       2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &       L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &          v1np1k*v2np1k*
     &           (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k))))*
     &     (-2*g*(m1 + 2*m2)*Sin(th1np1k) + 
     &       L2*m2*(v1n*v2n*Sin(th1np1k - th2n) + 
     &          v1np1k*v2np1k*Sin(th1np1k - th2np1k))))/
     &   (8.*(th1n - th1np1k)**2))*
     & ((L1*(v1n + v1np1k)*(-2*g*(m1 + 2*m2)*Cos(th1n) + 
     &       2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &       L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &          v1np1k*v2np1k*(Cos(th1n - th2np1k) - 
     &             Cos(th1np1k - th2np1k)))))/(8.*(th1n - th1np1k)) + 
     &  (L2*m2*(v2n + v2np1k)*(L1*v1np1k*v2np1k*
     &        (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &       L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k)) +
     &       2*g*(-Cos(th2n) + Cos(th2np1k))))/(8.*(th2n - th2np1k)) + 
     &  (L1*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &       3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &          v2n*Cos(th1np1k - th2np1k)))*
     &     ((-4*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k)/2.) + 
     &            m2*(3*g*Sin((th1n + th1np1k - 2*th2n - 2*th2np1k)/
     &                  2.) + 
     &               L2*(v2n + v2np1k)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &               (3*L1*(-v1n - v1np1k)**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k))/4.)))/
     &        (3.*L1*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((-3*(g*(4*m1 + 5*m2)*
     &                  Sin((th1n + th1np1k - th1n*xi1 + th1np1k*xi1)/
     &                   2.) + 
     &                 m2*(3*g*
     &                   Sin((th1n + th1np1k - 2*th2n - 2*th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + 2*th2n*xi1 - 
     &                  2*th2np1k*xi1)/2.) + 
     &                   L2*(v2n + v2np1k - v2n*xi1 + v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                  th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 + xi1) - v1np1k*(1 + xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4.)))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))) - (3*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k + th1n*xi1 - 
     &                   th1np1k*xi1)/2.) + 
     &                 m2*(L2*(v2n + v2np1k + v2n*xi1 - v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 - xi1) - v1np1k*(1 - xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4. + 
     &                   3*g*Sin((th1n + th1np1k - 2*th2n - 
     &                   2*th2np1k + th1n*xi1 - th1np1k*xi1 - 
     &                   2*th2n*xi1 + 2*th2np1k*xi1)/2.))))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                 th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & )))))/18.))/12. + (L2*m2*(2*L2*(v2n + v2np1k) + 
     &       3*L1*v1n*Cos(th1n - th2n) + 
     &       3*L1*v1np1k*Cos(th1np1k - th2np1k))*
     &     ((4*(3*L1*(m1 + 3*m2)*(v1n + v1np1k)**2*
     &             Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &            (9*L2*m2*(-v2n - v2np1k)**2*
     &               Sin(th1n + th1np1k - th2n - th2np1k))/4. + 
     &            9*g*(m1 + 2*m2)*
     &             Sin(th1n + th1np1k - th2n/2. - th2np1k/2.) - 
     &            3*g*(m1 + 6*m2)*Sin((th2n + th2np1k)/2.)))/
     &        (9.*L2*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((9*g*(m1 + 2*m2)*
     &                Sin(th1n + (th2n*(-1 - xi1))/2. + 
     &                  th1np1k*(1 - xi1) - (th2np1k*(1 - xi1))/2. + 
     &                  th1n*xi1) - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k + th2n*xi1 - th2np1k*xi1)/2.)
     & + 3*L1*(m1 + 3*m2)*(v1n + v1np1k + v1n*xi1 - v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 - xi1) - v2np1k*(1 - xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4.)/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                 th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))) + (3*L1*(m1 + 3*m2)*(v1n + v1np1k - v1n*xi1 + v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 + xi1) - v2np1k*(1 + xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4. - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k - th2n*xi1 + th2np1k*xi1)/2.)
     & + 9*g*(m1 + 2*m2)*Sin(th1n + (th2n*(-1 + xi1))/2. - th1n*xi1 + 
     &                  th1np1k*(1 + xi1) - (th2np1k*(1 + xi1))/2.))/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                 th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)
     & ))))/18.))/12.))/
     &    (12.*((L1**2*(-2*g*(m1 + 2*m2)*Cos(th1n) + 
     &         2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &         L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) +
     &            v1np1k*v2np1k*
     &             (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k))))**2)/
     &    (16.*(th1n - th1np1k)**2) + 
     &   (L2**2*m2**2*(2*L2*(v2n + v2np1k) + 3*L1*v1n*Cos(th1n - th2n) +
     &         3*L1*v1np1k*Cos(th1np1k - th2np1k))**2)/144. + 
     &   (L1**2*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &         3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &            v2n*Cos(th1np1k - th2np1k)))**2)/144. + 
     &   (L2**2*m2**2*(L1*v1np1k*v2np1k*
     &          (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &         L1*v1n*v2n*(Cos(th1np1k - th2n) - 
     &            Cos(th1np1k - th2np1k)) + 
     &         2*g*(-Cos(th2n) + Cos(th2np1k)))**2)/
     &    (16.*(th2n - th2np1k)**2))**2) + 
     &   (L1*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &  3*L2*m2*(v2np1k*Cos(th1n - th2n) + v2n*Cos(th1np1k - th2np1k)))*
     & ((L1*(v1n + v1np1k)*(-2*g*(m1 + 2*m2)*Cos(th1n) + 
     &       2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &       L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &          v1np1k*v2np1k*(Cos(th1n - th2np1k) - 
     &             Cos(th1np1k - th2np1k)))))/(8.*(th1n - th1np1k)**2) + 
     &  (L1*L2*m2*v1n*v2n*(v2n + v2np1k)*
     &     (-Sin(th1np1k - th2n) + Sin(th1np1k - th2np1k)))/
     &   (8.*(th2n - th2np1k)) + 
     &  (L1*(v1n + v1np1k)*(-2*g*(m1 + 2*m2)*Sin(th1np1k) + 
     &       L2*m2*(v1n*v2n*Sin(th1np1k - th2n) + 
     &          v1np1k*v2np1k*Sin(th1np1k - th2np1k))))/
     &   (8.*(th1n - th1np1k)) - 
     &  (L1*L2*m2*v2n*Sin(th1np1k - th2np1k)*
     &     ((-4*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k)/2.) + 
     &            m2*(3*g*Sin((th1n + th1np1k - 2*th2n - 2*th2np1k)/
     &                  2.) + 
     &               L2*(v2n + v2np1k)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &               (3*L1*(-v1n - v1np1k)**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k))/4.)))/
     &        (3.*L1*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((-3*(g*(4*m1 + 5*m2)*
     &                  Sin((th1n + th1np1k - th1n*xi1 + th1np1k*xi1)/
     &                   2.) + 
     &                 m2*(3*g*
     &                   Sin((th1n + th1np1k - 2*th2n - 2*th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + 2*th2n*xi1 - 
     &                   2*th2np1k*xi1)/2.) + 
     &                   L2*(v2n + v2np1k - v2n*xi1 + v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 + xi1) - v1np1k*(1 + xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4.)))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                 th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))) - (3*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k + th1n*xi1 - 
     &                   th1np1k*xi1)/2.) + 
     &                 m2*(L2*(v2n + v2np1k + v2n*xi1 - v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &                    (3*L1*(v1n*(-1 - xi1) - v1np1k*(1 - xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4. + 
     &                    3*g*
     &                    Sin((th1n + th1np1k - 2*th2n - 2*th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - 2*th2n*xi1 + 
     &                   2*th2np1k*xi1)/2.))))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                 th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1)
     & ))))/18.))/4. + (L1*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &       3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &          v2n*Cos(th1np1k - th2np1k)))*
     &     ((-4*((g*(4*m1 + 5*m2)*Cos((th1n + th1np1k)/2.))/2. + 
     &            m2*((3*g*Cos((th1n + th1np1k - 2*th2n - 2*th2np1k)/
     &                   2.))/2. + 
     &               (L2*(v2n + v2np1k)**2*
     &                  Cos((th1n + th1np1k - th2n - th2np1k)/2.))/2. +
     &               (3*L1*(-v1n - v1np1k)**2*
     &                  Cos(th1n + th1np1k - th2n - th2np1k))/4.)))/
     &        (3.*L1*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (12*m2*Sin(th1n + th1np1k - th2n - th2np1k)*
     &          (g*(4*m1 + 5*m2)*Sin((th1n + th1np1k)/2.) + 
     &            m2*(3*g*Sin((th1n + th1np1k - 2*th2n - 2*th2np1k)/
     &                  2.) + 
     &               L2*(v2n + v2np1k)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &               (3*L1*(-v1n - v1np1k)**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k))/4.)))/
     &        (L1*(8*m1 + 15*m2 - 
     &             9*m2*Cos(th1n + th1np1k - th2n - th2np1k))**2) + 
     &       (5*((-3*((g*(4*m1 + 5*m2)*(1 + xi1)*
     &                   Cos((th1n + th1np1k - th1n*xi1 + 
     &                  th1np1k*xi1)/2.))/2. + 
     &                 m2*((3*g*(1 + xi1)*
     &                   Cos((th1n + th1np1k - 2*th2n - 2*th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + 2*th2n*xi1 - 
     &                  2*th2np1k*xi1)/2.))/2. + 
     &                   (L2*(1 + xi1)*
     &                   (v2n + v2np1k - v2n*xi1 + v2np1k*xi1)**2*
     &                   Cos((th1n + th1np1k - th2n - th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                  th2np1k*xi1)/2.))/2. + 
     &                   (3*L1*(1 + xi1)*
     &                   (v1n*(-1 + xi1) - v1np1k*(1 + xi1))**2*
     &                   Cos(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4.)))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                 th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))) - (3*((g*(4*m1 + 5*m2)*(1 - xi1)*
     &                   Cos((th1n + th1np1k + th1n*xi1 - 
     &                  th1np1k*xi1)/2.))/2. + 
     &                 m2*((L2*(1 - xi1)*
     &                   (v2n + v2np1k + v2n*xi1 - v2np1k*xi1)**2*
     &                   Cos((th1n + th1np1k - th2n - th2np1k + 
     &                  th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                  th2np1k*xi1)/2.))/2. + 
     &                   (3*L1*(v1n*(-1 - xi1) - v1np1k*(1 - xi1))**2*
     &                   (1 - xi1)*
     &                   Cos(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4. + 
     &                   (3*g*(1 - xi1)*
     &                   Cos((th1n + th1np1k - 2*th2n - 2*th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - 2*th2n*xi1 + 
     &                   2*th2np1k*xi1)/2.))/2.)))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                 th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))) + (27*m2*(1 + xi1)*Sin(th1n + th1np1k - th2n - th2np1k - 
     &                 th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)*
     &               (g*(4*m1 + 5*m2)*
     &                  Sin((th1n + th1np1k - th1n*xi1 + th1np1k*xi1)/
     &                   2.) + 
     &                 m2*(3*g*
     &                   Sin((th1n + th1np1k - 2*th2n - 2*th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + 2*th2n*xi1 - 
     &                   2*th2np1k*xi1)/2.) + 
     &                   L2*(v2n + v2np1k - v2n*xi1 + v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 + xi1) - v1np1k*(1 + xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4.)))/
     &             (L1*(8*m1 + 15*m2 - 
     &                  9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))**2) + (27*m2*(1 - xi1)*Sin(th1n + th1np1k - th2n - th2np1k + 
     &                 th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1)*
     &               (g*(4*m1 + 5*m2)*
     &                  Sin((th1n + th1np1k + th1n*xi1 - th1np1k*xi1)/
     &                   2.) + 
     &                 m2*(L2*(v2n + v2np1k + v2n*xi1 - v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &                    (3*L1*(v1n*(-1 - xi1) - v1np1k*(1 - xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4. + 
     &                    3*g*
     &                    Sin((th1n + th1np1k - 2*th2n - 2*th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - 2*th2n*xi1 + 
     &                   2*th2np1k*xi1)/2.))))/
     &             (L1*(8*m1 + 15*m2 - 
     &                  9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                 th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))**2)))/18.))/12. - (L1*L2*m2*v1np1k*Sin(th1np1k - th2np1k)*
     &     ((4*(3*L1*(m1 + 3*m2)*(v1n + v1np1k)**2*
     &             Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &            (9*L2*m2*(-v2n - v2np1k)**2*
     &               Sin(th1n + th1np1k - th2n - th2np1k))/4. + 
     &            9*g*(m1 + 2*m2)*
     &             Sin(th1n + th1np1k - th2n/2. - th2np1k/2.) - 
     &            3*g*(m1 + 6*m2)*Sin((th2n + th2np1k)/2.)))/
     &        (9.*L2*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((9*g*(m1 + 2*m2)*
     &                Sin(th1n + (th2n*(-1 - xi1))/2. + 
     &                  th1np1k*(1 - xi1) - (th2np1k*(1 - xi1))/2. + 
     &                  th1n*xi1) - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k + th2n*xi1 - th2np1k*xi1)/2.)
     & + 3*L1*(m1 + 3*m2)*(v1n + v1np1k + v1n*xi1 - v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 - xi1) - v2np1k*(1 - xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4.)/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                 th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))) + (3*L1*(m1 + 3*m2)*(v1n + v1np1k - v1n*xi1 + v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 + xi1) - v2np1k*(1 + xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4. - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k - th2n*xi1 + th2np1k*xi1)/2.)
     & + 9*g*(m1 + 2*m2)*Sin(th1n + (th2n*(-1 + xi1))/2. - th1n*xi1 + 
     &                  th1np1k*(1 + xi1) - (th2np1k*(1 + xi1))/2.))/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                 th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)
     & ))))/18.))/4.+(L2*m2*(2*L2*(v2n+v2np1k)+3*L1*v1n*Cos(th1n-th2n)+ 
     &       3*L1*v1np1k*Cos(th1np1k - th2np1k))*
     &     ((4*((3*L1*(m1 + 3*m2)*(v1n + v1np1k)**2*
     &               Cos((th1n + th1np1k - th2n - th2np1k)/2.))/2. + 
     &            (9*L2*m2*(-v2n - v2np1k)**2*
     &               Cos(th1n + th1np1k - th2n - th2np1k))/4. + 
     &            9*g*(m1 + 2*m2)*
     &             Cos(th1n + th1np1k - th2n/2. - th2np1k/2.)))/
     &        (9.*L2*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) - 
     &       (4*m2*Sin(th1n + th1np1k - th2n - th2np1k)*
     &          (3*L1*(m1 + 3*m2)*(v1n + v1np1k)**2*
     &             Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &            (9*L2*m2*(-v2n - v2np1k)**2*
     &               Sin(th1n + th1np1k - th2n - th2np1k))/4. + 
     &            9*g*(m1 + 2*m2)*
     &             Sin(th1n + th1np1k - th2n/2. - th2np1k/2.) - 
     &            3*g*(m1 + 6*m2)*Sin((th2n + th2np1k)/2.)))/
     &        (L2*(8*m1 + 15*m2 - 
     &             9*m2*Cos(th1n + th1np1k - th2n - th2np1k))**2) + 
     &       (5*((9*g*(m1 + 2*m2)*(1 - xi1)*
     &                Cos(th1n + (th2n*(-1 - xi1))/2. + 
     &                  th1np1k*(1 - xi1) - (th2np1k*(1 - xi1))/2. + 
     &                  th1n*xi1) + 
     &               (3*L1*(m1 + 3*m2)*(1 - xi1)*
     &                  (v1n + v1np1k + v1n*xi1 - v1np1k*xi1)**2*
     &                  Cos((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.))/2. + 
     &               (9*L2*m2*(v2n*(-1 - xi1) - v2np1k*(1 - xi1))**2*
     &                  (1 - xi1)*
     &                  Cos(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))/4.)/(L2*(8*m1 + 15*m2 - 9*m2*
     &                  Cos(th1n + th1np1k - th2n - th2np1k + th1n*xi1 - 
     &                    th1np1k*xi1 - th2n*xi1 + th2np1k*xi1))) + 
     &            ((3*L1*(m1 + 3*m2)*(1 + xi1)*
     &                  (v1n + v1np1k - v1n*xi1 + v1np1k*xi1)**2*
     &                  Cos((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.))/2. + 
     &               (9*L2*m2*(1 + xi1)*
     &                  (v2n*(-1 + xi1) - v2np1k*(1 + xi1))**2*
     &                  Cos(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4. + 
     &               9*g*(m1 + 2*m2)*(1 + xi1)*
     &                Cos(th1n + (th2n*(-1 + xi1))/2. - th1n*xi1 + 
     &                  th1np1k*(1 + xi1) - (th2np1k*(1 + xi1))/2.))/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)
     & ))-(9*m2*(1-xi1)*Sin(th1n + th1np1k - th2n - th2np1k + th1n*xi1-
     &                 th1np1k*xi1 - th2n*xi1 + th2np1k*xi1)*
     &               (9*g*(m1 + 2*m2)*
     &                  Sin(th1n + (th2n*(-1 - xi1))/2. + 
     &                   th1np1k*(1 - xi1) - (th2np1k*(1 - xi1))/2. + 
     &                   th1n*xi1) - 
     &                 3*g*(m1 + 6*m2)*
     &                  Sin((th2n + th2np1k + th2n*xi1 - th2np1k*xi1)/
     &                   2.) + 
     &                 3*L1*(m1 + 3*m2)*
     &                  (v1n + v1np1k + v1n*xi1 - v1np1k*xi1)**2*
     &                  Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &                 (9*L2*m2*(v2n*(-1 - xi1) - v2np1k*(1 - xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))/4.))/(L2*(8*m1 + 15*m2 - 9*m2*
     &                   Cos(th1n + th1np1k - th2n - th2np1k + 
     &                 th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))**2) - (9*m2*(1 + xi1)*Sin(th1n + th1np1k - th2n - th2np1k - 
     &                 th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)*
     &               (3*L1*(m1 + 3*m2)*
     &                  (v1n + v1np1k - v1n*xi1 + v1np1k*xi1)**2*
     &                  Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & )/2.) + (9*L2*m2*(v2n*(-1 + xi1) - v2np1k*(1 + xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))/4.-3*g*(m1+6*m2)*Sin((th2n+th2np1k - th2n*xi1 + th2np1k*xi1)/
     &                   2.) + 
     &                 9*g*(m1 + 2*m2)*
     &                  Sin(th1n + (th2n*(-1 + xi1))/2. - th1n*xi1 + 
     &                    th1np1k*(1 + xi1) - (th2np1k*(1 + xi1))/2.)))/
     &             (L2*(8*m1 + 15*m2 - 
     &                  9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)
     & )**2)))/18.))/12.))/
     &    (12.*((L1**2*(-2*g*(m1 + 2*m2)*Cos(th1n) + 
     &        2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &        L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &           v1np1k*v2np1k*
     &            (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k))))**2)/
     &   (16.*(th1n - th1np1k)**2) + 
     &  (L2**2*m2**2*(2*L2*(v2n + v2np1k) + 3*L1*v1n*Cos(th1n - th2n) + 
     &        3*L1*v1np1k*Cos(th1np1k - th2np1k))**2)/144. + 
     &  (L1**2*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &        3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &           v2n*Cos(th1np1k - th2np1k)))**2)/144. + 
     &  (L2**2*m2**2*(L1*v1np1k*v2np1k*
     &         (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &        L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k))+
     &  2*g*(-Cos(th2n) + Cos(th2np1k)))**2)/(16.*(th2n - th2np1k)**2)))

      m32=(4*m2*(-3*g*Cos((th1n + th1np1k - 2*th2n - 2*th2np1k)/2.) - 
     & (L2*(v2n + v2np1k)**2*Cos((th1n + th1np1k - th2n - th2np1k)/2.))/
     &   2. - (3*L1*(-v1n - v1np1k)**2*
     &     Cos(th1n + th1np1k - th2n - th2np1k))/4.))/
     & (3.*L1*(8*m1+15*m2-9*m2*Cos(th1n + th1np1k - th2n - th2np1k)))+
     & (12*m2*Sin(th1n + th1np1k - th2n - th2np1k)*
     &(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k)/2.) + 
     &  m2*(3*g*Sin((th1n + th1np1k - 2*th2n - 2*th2np1k)/2.) + 
     & L2*(v2n + v2np1k)**2*Sin((th1n + th1np1k - th2n - th2np1k)/2.)+
     &(3*L1*(-v1n - v1np1k)**2*Sin(th1n + th1np1k - th2n - th2np1k))/4.
     & )))/(L1*(8*m1+15*m2-9*m2*Cos(th1n+th1np1k - th2n - th2np1k))**2)-
     &   (5*((-3*m2*((3*g*(-2 - 2*xi1)*
     &          Cos((th1n + th1np1k - 2*th2n - 2*th2np1k - th1n*xi1 + 
     &              th1np1k*xi1 + 2*th2n*xi1 - 2*th2np1k*xi1)/2.))/2. + 
     &       (L2*(-1 - xi1)*(v2n + v2np1k - v2n*xi1 + v2np1k*xi1)**2*
     &          Cos((th1n + th1np1k - th2n - th2np1k - th1n*xi1 + 
     &              th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)/2.))/2. + 
     &       (3*L1*(-1 - xi1)*(v1n*(-1 + xi1) - v1np1k*(1 + xi1))**2*
     &          Cos(th1n + th1np1k - th2n - th2np1k - th1n*xi1 + 
     &            th1np1k*xi1 + th2n*xi1 - th2np1k*xi1))/4.))/
     &   (L1*(8*m1 + 15*m2 - 9*m2*
     &        Cos(th1n + th1np1k - th2n - th2np1k - th1n*xi1 + 
     &          th1np1k*xi1 + th2n*xi1 - th2np1k*xi1))) - 
     &  (3*m2*((L2*(-1 + xi1)*(v2n + v2np1k + v2n*xi1 - v2np1k*xi1)**2*
     &          Cos((th1n + th1np1k - th2n - th2np1k + th1n*xi1 - 
     &              th1np1k*xi1 - th2n*xi1 + th2np1k*xi1)/2.))/2. + 
     &       (3*L1*(v1n*(-1 - xi1) - v1np1k*(1 - xi1))**2*(-1 + xi1)*
     &          Cos(th1n + th1np1k - th2n - th2np1k + th1n*xi1 - 
     &            th1np1k*xi1 - th2n*xi1 + th2np1k*xi1))/4. + 
     &       (3*g*(-2 + 2*xi1)*
     &          Cos((th1n + th1np1k - 2*th2n - 2*th2np1k + th1n*xi1 - 
     &              th1np1k*xi1 - 2*th2n*xi1 + 2*th2np1k*xi1)/2.))/2.))/
     &   (L1*(8*m1 + 15*m2 - 9*m2*
     &        Cos(th1n + th1np1k - th2n - th2np1k + th1n*xi1 - 
     &          th1np1k*xi1 - th2n*xi1 + th2np1k*xi1))) + 
     &  (27*m2*(-1 - xi1)*Sin(th1n + th1np1k - th2n - th2np1k - 
     &       th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)*
     &     (g*(4*m1 + 5*m2)*Sin((th1n + th1np1k - th1n*xi1 + 
     &            th1np1k*xi1)/2.) + 
     &       m2*(3*g*Sin((th1n + th1np1k - 2*th2n - 2*th2np1k - 
     &               th1n*xi1 + th1np1k*xi1 + 2*th2n*xi1 - 2*th2np1k*xi1
     & )/2.) + L2*(v2n + v2np1k - v2n*xi1 + v2np1k*xi1)**2*
     &           Sin((th1n + th1np1k - th2n - th2np1k - th1n*xi1 + 
     &               th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)/2.) + 
     &          (3*L1*(v1n*(-1 + xi1) - v1np1k*(1 + xi1))**2*
     &             Sin(th1n + th1np1k - th2n - th2np1k - th1n*xi1 + 
     &               th1np1k*xi1 + th2n*xi1 - th2np1k*xi1))/4.)))/
     &   (L1*(8*m1 + 15*m2 - 9*m2*
     &         Cos(th1n + th1np1k - th2n - th2np1k - th1n*xi1 + 
     &           th1np1k*xi1 + th2n*xi1 - th2np1k*xi1))**2) + 
     &  (27*m2*(-1 + xi1)*Sin(th1n + th1np1k - th2n - th2np1k + 
     &       th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1)*
     &   (g*(4*m1 + 5*m2)*Sin((th1n + th1np1k + th1n*xi1 - th1np1k*xi1)/
     &          2.) + m2*(L2*(v2n + v2np1k + v2n*xi1 - v2np1k*xi1)**2*
     &           Sin((th1n + th1np1k - th2n - th2np1k + th1n*xi1 - 
     &               th1np1k*xi1 - th2n*xi1 + th2np1k*xi1)/2.) + 
     &          (3*L1*(v1n*(-1 - xi1) - v1np1k*(1 - xi1))**2*
     &             Sin(th1n + th1np1k - th2n - th2np1k + th1n*xi1 - 
     &               th1np1k*xi1 - th2n*xi1 + th2np1k*xi1))/4. + 
     &          3*g*Sin((th1n + th1np1k - 2*th2n - 2*th2np1k + th1n*xi1-
     &               th1np1k*xi1 - 2*th2n*xi1 + 2*th2np1k*xi1)/2.))))/
     &   (L1*(8*m1 + 15*m2 - 9*m2*
     &         Cos(th1n + th1np1k - th2n - th2np1k + th1n*xi1 - 
     &           th1np1k*xi1 - th2n*xi1 + th2np1k*xi1))**2)))/18. + 
     &   (L1*L2*m2*v2n*Sin(th1np1k - th2np1k)*
     & ((L1*(v1n + v1np1k)*(-2*g*(m1 + 2*m2)*Cos(th1n) + 
     &       2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &       L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &          v1np1k*v2np1k*(Cos(th1n - th2np1k) - 
     &             Cos(th1np1k - th2np1k)))))/(8.*(th1n - th1np1k)) + 
     &  (L2*m2*(v2n + v2np1k)*(L1*v1np1k*v2np1k*
     &        (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &       L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k)) +
     &       2*g*(-Cos(th2n) + Cos(th2np1k))))/(8.*(th2n - th2np1k)) + 
     &  (L1*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &       3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &          v2n*Cos(th1np1k - th2np1k)))*
     &     ((-4*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k)/2.) + 
     &            m2*(3*g*Sin((th1n + th1np1k - 2*th2n - 2*th2np1k)/
     &                  2.) + 
     &               L2*(v2n + v2np1k)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &               (3*L1*(-v1n - v1np1k)**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k))/4.)))/
     &        (3.*L1*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((-3*(g*(4*m1 + 5*m2)*
     &                  Sin((th1n + th1np1k - th1n*xi1 + th1np1k*xi1)/
     &                   2.) + 
     &                 m2*(3*g*
     &                   Sin((th1n + th1np1k - 2*th2n - 2*th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + 2*th2n*xi1 - 
     &                  2*th2np1k*xi1)/2.) + 
     &                   L2*(v2n + v2np1k - v2n*xi1 + v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                  th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 + xi1) - v1np1k*(1 + xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4.)))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))) - (3*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k + th1n*xi1 - 
     &                   th1np1k*xi1)/2.) + 
     &                 m2*(L2*(v2n + v2np1k + v2n*xi1 - v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 - xi1) - v1np1k*(1 - xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4. + 
     &                   3*g*Sin((th1n + th1np1k - 2*th2n - 
     &                   2*th2np1k + th1n*xi1 - th1np1k*xi1 - 
     &                   2*th2n*xi1 + 2*th2np1k*xi1)/2.))))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                  th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & )))))/18.))/12. + (L2*m2*(2*L2*(v2n + v2np1k) + 
     &       3*L1*v1n*Cos(th1n - th2n) + 
     &       3*L1*v1np1k*Cos(th1np1k - th2np1k))*
     &     ((4*(3*L1*(m1 + 3*m2)*(v1n + v1np1k)**2*
     &             Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &            (9*L2*m2*(-v2n - v2np1k)**2*
     &               Sin(th1n + th1np1k - th2n - th2np1k))/4. + 
     &            9*g*(m1 + 2*m2)*
     &             Sin(th1n + th1np1k - th2n/2. - th2np1k/2.) - 
     &            3*g*(m1 + 6*m2)*Sin((th2n + th2np1k)/2.)))/
     &        (9.*L2*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((9*g*(m1 + 2*m2)*
     &                Sin(th1n + (th2n*(-1 - xi1))/2. + 
     &                  th1np1k*(1 - xi1) - (th2np1k*(1 - xi1))/2. + 
     &                  th1n*xi1) - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k + th2n*xi1 - th2np1k*xi1)/2.)
     & + 3*L1*(m1 + 3*m2)*(v1n + v1np1k + v1n*xi1 - v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 - xi1) - v2np1k*(1 - xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4.)/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                 th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))) + (3*L1*(m1 + 3*m2)*(v1n + v1np1k - v1n*xi1 + v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 + xi1) - v2np1k*(1 + xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4. - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k - th2n*xi1 + th2np1k*xi1)/2.)
     & + 9*g*(m1 + 2*m2)*Sin(th1n + (th2n*(-1 + xi1))/2. - th1n*xi1 + 
     &                  th1np1k*(1 + xi1) - (th2np1k*(1 + xi1))/2.))/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)
     & ))))/18.))/12.))/
     &    (4.*((L1**2*(-2*g*(m1 + 2*m2)*Cos(th1n) + 
     &        2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &        L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &           v1np1k*v2np1k*
     &            (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k))))**2)/
     &   (16.*(th1n - th1np1k)**2) + 
     &  (L2**2*m2**2*(2*L2*(v2n + v2np1k) + 3*L1*v1n*Cos(th1n - th2n) +
     &        3*L1*v1np1k*Cos(th1np1k - th2np1k))**2)/144. + 
     &  (L1**2*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &        3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &           v2n*Cos(th1np1k - th2np1k)))**2)/144. + 
     &  (L2**2*m2**2*(L1*v1np1k*v2np1k*
     &         (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &        L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k))+
     &  2*g*(-Cos(th2n) + Cos(th2np1k)))**2)/(16.*(th2n - th2np1k)**2))
     & ) - (L1*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &  3*L2*m2*(v2np1k*Cos(th1n - th2n) + v2n*Cos(th1np1k - th2np1k)))*
     & ((L2**2*m2**2*(L1*v1np1k*v2np1k*
     &         (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &        L1*v1n*v2n*(Cos(th1np1k - th2n) - 
     &           Cos(th1np1k - th2np1k)) + 
     &        2*g*(-Cos(th2n) + Cos(th2np1k)))**2)/
     &   (8.*(th2n - th2np1k)**3) + 
     &  (L1**2*L2*m2*v1np1k*v2np1k*
     &     (-2*g*(m1 + 2*m2)*Cos(th1n) + 2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &       L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &          v1np1k*v2np1k*
     &           (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k))))*
     &     (Sin(th1n - th2np1k) - Sin(th1np1k - th2np1k)))/
     &   (8.*(th1n - th1np1k)**2) + 
     &  (L1*L2**2*m2**2*v1np1k*
     &     (2*L2*(v2n + v2np1k) + 3*L1*v1n*Cos(th1n - th2n) + 
     &     3*L1*v1np1k*Cos(th1np1k - th2np1k))*Sin(th1np1k - th2np1k))/
     &   24. + (L1**2*L2*m2*v2n*
     &     (2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &       3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &      v2n*Cos(th1np1k - th2np1k)))*Sin(th1np1k - th2np1k))/24.+
     &  (L2**2*m2**2*(L1*v1np1k*v2np1k*
     &        (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &       L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k)) +
     &       2*g*(-Cos(th2n) + Cos(th2np1k)))*
     &     (-(L1*v1np1k*v2np1k*Sin(th1n - th2np1k)) - 
     &       L1*v1n*v2n*Sin(th1np1k - th2np1k) - 2*g*Sin(th2np1k)))/
     &   (8.*(th2n - th2np1k)**2))*
     & ((L1*(v1n + v1np1k)*(-2*g*(m1 + 2*m2)*Cos(th1n) + 
     &       2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &       L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &          v1np1k*v2np1k*(Cos(th1n - th2np1k) - 
     &             Cos(th1np1k - th2np1k)))))/(8.*(th1n - th1np1k)) + 
     &  (L2*m2*(v2n + v2np1k)*(L1*v1np1k*v2np1k*
     &        (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &       L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k)) +
     &       2*g*(-Cos(th2n) + Cos(th2np1k))))/(8.*(th2n - th2np1k)) + 
     &  (L1*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &       3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &          v2n*Cos(th1np1k - th2np1k)))*
     &     ((-4*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k)/2.) + 
     &            m2*(3*g*Sin((th1n + th1np1k - 2*th2n - 2*th2np1k)/
     &                  2.) + 
     &               L2*(v2n + v2np1k)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &               (3*L1*(-v1n - v1np1k)**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k))/4.)))/
     &        (3.*L1*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((-3*(g*(4*m1 + 5*m2)*
     &                  Sin((th1n + th1np1k - th1n*xi1 + th1np1k*xi1)/
     &                   2.) + 
     &                 m2*(3*g*
     &                   Sin((th1n + th1np1k - 2*th2n - 2*th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + 2*th2n*xi1 - 
     &                  2*th2np1k*xi1)/2.) + 
     &                   L2*(v2n + v2np1k - v2n*xi1 + v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                  th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 + xi1) - v1np1k*(1 + xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4.)))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))) - (3*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k + th1n*xi1 - 
     &                   th1np1k*xi1)/2.) + 
     &                 m2*(L2*(v2n + v2np1k + v2n*xi1 - v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 - xi1) - v1np1k*(1 - xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4. + 
     &                   3*g*Sin((th1n + th1np1k - 2*th2n - 
     &                   2*th2np1k + th1n*xi1 - th1np1k*xi1 - 
     &                   2*th2n*xi1 + 2*th2np1k*xi1)/2.))))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                 th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & )))))/18.))/12. + (L2*m2*(2*L2*(v2n + v2np1k) + 
     &       3*L1*v1n*Cos(th1n - th2n) + 
     &       3*L1*v1np1k*Cos(th1np1k - th2np1k))*
     &     ((4*(3*L1*(m1 + 3*m2)*(v1n + v1np1k)**2*
     &             Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &            (9*L2*m2*(-v2n - v2np1k)**2*
     &               Sin(th1n + th1np1k - th2n - th2np1k))/4. + 
     &            9*g*(m1 + 2*m2)*
     &             Sin(th1n + th1np1k - th2n/2. - th2np1k/2.) - 
     &            3*g*(m1 + 6*m2)*Sin((th2n + th2np1k)/2.)))/
     &        (9.*L2*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((9*g*(m1 + 2*m2)*
     &                Sin(th1n + (th2n*(-1 - xi1))/2. + 
     &                  th1np1k*(1 - xi1) - (th2np1k*(1 - xi1))/2. + 
     &                  th1n*xi1) - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k + th2n*xi1 - th2np1k*xi1)/2.)
     & + 3*L1*(m1 + 3*m2)*(v1n + v1np1k + v1n*xi1 - v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 - xi1) - v2np1k*(1 - xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4.)/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))) + (3*L1*(m1 + 3*m2)*(v1n + v1np1k - v1n*xi1 + v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 + xi1) - v2np1k*(1 + xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4. - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k - th2n*xi1 + th2np1k*xi1)/2.)
     & + 9*g*(m1 + 2*m2)*Sin(th1n + (th2n*(-1 + xi1))/2. - th1n*xi1 + 
     &                  th1np1k*(1 + xi1) - (th2np1k*(1 + xi1))/2.))/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &              th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)
     & ))))/18.))/12.))/
     &    (12.*((L1**2*(-2*g*(m1 + 2*m2)*Cos(th1n) + 
     &         2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &         L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &            v1np1k*v2np1k*
     &             (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k))))**2)/
     &    (16.*(th1n - th1np1k)**2) + 
     &   (L2**2*m2**2*(2*L2*(v2n + v2np1k) + 3*L1*v1n*Cos(th1n - th2n) + 
     &         3*L1*v1np1k*Cos(th1np1k - th2np1k))**2)/144. + 
     &   (L1**2*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &         3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &            v2n*Cos(th1np1k - th2np1k)))**2)/144. + 
     &   (L2**2*m2**2*(L1*v1np1k*v2np1k*
     &          (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &         L1*v1n*v2n*(Cos(th1np1k - th2n) - 
     &            Cos(th1np1k - th2np1k)) + 
     &         2*g*(-Cos(th2n) + Cos(th2np1k)))**2)/
     &    (16.*(th2n - th2np1k)**2))**2) + 
     &   (L1*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &  3*L2*m2*(v2np1k*Cos(th1n - th2n) + v2n*Cos(th1np1k - th2np1k)))*
     & ((L2*m2*(v2n + v2np1k)*(L1*v1np1k*v2np1k*
     &        (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &       L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k)) + 
     &       2*g*(-Cos(th2n) + Cos(th2np1k))))/(8.*(th2n - th2np1k)**2)+
     &  (L1*L2*m2*v1np1k*(v1n + v1np1k)*v2np1k*
     &     (Sin(th1n - th2np1k) - Sin(th1np1k - th2np1k)))/
     &   (8.*(th1n - th1np1k)) + 
     &  (L2*m2*(v2n + v2np1k)*(-(L1*v1np1k*v2np1k*Sin(th1n - th2np1k)) - 
     &       L1*v1n*v2n*Sin(th1np1k - th2np1k) - 2*g*Sin(th2np1k)))/
     &   (8.*(th2n - th2np1k)) + 
     &  (L1*L2*m2*v2n*Sin(th1np1k - th2np1k)*
     &     ((-4*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k)/2.) + 
     &            m2*(3*g*Sin((th1n + th1np1k - 2*th2n - 2*th2np1k)/
     &                  2.) + 
     &               L2*(v2n + v2np1k)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &               (3*L1*(-v1n - v1np1k)**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k))/4.)))/
     &        (3.*L1*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((-3*(g*(4*m1 + 5*m2)*
     &                  Sin((th1n + th1np1k - th1n*xi1 + th1np1k*xi1)/
     &                   2.) + 
     &                 m2*(3*g*
     &                   Sin((th1n + th1np1k - 2*th2n - 2*th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + 2*th2n*xi1 - 
     &                   2*th2np1k*xi1)/2.) + 
     &                   L2*(v2n + v2np1k - v2n*xi1 + v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 + xi1) - v1np1k*(1 + xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4.)))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))) - (3*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k + th1n*xi1 - 
     &                   th1np1k*xi1)/2.) + 
     &                 m2*(L2*(v2n + v2np1k + v2n*xi1 - v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &                    (3*L1*(v1n*(-1 - xi1) - v1np1k*(1 - xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4. + 
     &                    3*g*
     &                    Sin((th1n + th1np1k - 2*th2n - 2*th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - 2*th2n*xi1 + 
     &                   2*th2np1k*xi1)/2.))))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1)
     & ))))/18.))/4. + (L1*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &       3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &          v2n*Cos(th1np1k - th2np1k)))*
     &     ((-4*m2*(-3*g*Cos((th1n + th1np1k - 2*th2n - 2*th2np1k)/2.) -
     &            (L2*(v2n + v2np1k)**2*
     &               Cos((th1n + th1np1k - th2n - th2np1k)/2.))/2. - 
     &            (3*L1*(-v1n - v1np1k)**2*
     &               Cos(th1n + th1np1k - th2n - th2np1k))/4.))/
     &        (3.*L1*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) - 
     &       (12*m2*Sin(th1n + th1np1k - th2n - th2np1k)*
     &          (g*(4*m1 + 5*m2)*Sin((th1n + th1np1k)/2.) + 
     &            m2*(3*g*Sin((th1n + th1np1k - 2*th2n - 2*th2np1k)/
     &                  2.) + 
     &               L2*(v2n + v2np1k)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &               (3*L1*(-v1n - v1np1k)**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k))/4.)))/
     &        (L1*(8*m1 + 15*m2 - 
     &             9*m2*Cos(th1n + th1np1k - th2n - th2np1k))**2) + 
     &       (5*((-3*m2*((3*g*(-2 - 2*xi1)*
     &                   Cos((th1n + th1np1k - 2*th2n - 2*th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + 2*th2n*xi1 - 
     &                  2*th2np1k*xi1)/2.))/2. + 
     &                 (L2*(-1 - xi1)*
     &                   (v2n + v2np1k - v2n*xi1 + v2np1k*xi1)**2*
     &                   Cos((th1n + th1np1k - th2n - th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                  th2np1k*xi1)/2.))/2. + 
     &                 (3*L1*(-1 - xi1)*
     &                   (v1n*(-1 + xi1) - v1np1k*(1 + xi1))**2*
     &                   Cos(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4.))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                 th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & )))-(3*m2*((L2*(-1+xi1)*(v2n+v2np1k + v2n*xi1 - v2np1k*xi1)**2*
     &                   Cos((th1n + th1np1k - th2n - th2np1k + 
     &                  th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                  th2np1k*xi1)/2.))/2. + 
     &                 (3*L1*(v1n*(-1 - xi1) - v1np1k*(1 - xi1))**2*
     &                   (-1 + xi1)*
     &                   Cos(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4. + 
     &                 (3*g*(-2 + 2*xi1)*
     &                   Cos((th1n + th1np1k - 2*th2n - 2*th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - 2*th2n*xi1 + 
     &                   2*th2np1k*xi1)/2.))/2.))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                 th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))) + (27*m2*(-1 - xi1)*Sin(th1n + th1np1k - th2n - th2np1k - 
     &                 th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)*
     &               (g*(4*m1 + 5*m2)*
     &                  Sin((th1n + th1np1k - th1n*xi1 + th1np1k*xi1)/
     &                   2.) + 
     &                 m2*(3*g*
     &                   Sin((th1n + th1np1k - 2*th2n - 2*th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + 2*th2n*xi1 - 
     &                   2*th2np1k*xi1)/2.) + 
     &                   L2*(v2n + v2np1k - v2n*xi1 + v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 + xi1) - v1np1k*(1 + xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4.)))/
     &             (L1*(8*m1 + 15*m2 - 
     &                  9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))**2) + (27*m2*(-1 + xi1)*Sin(th1n + th1np1k - th2n - th2np1k + 
     &                 th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1)*
     &               (g*(4*m1 + 5*m2)*
     &                  Sin((th1n + th1np1k + th1n*xi1 - th1np1k*xi1)/
     &                   2.) + 
     &                 m2*(L2*(v2n + v2np1k + v2n*xi1 - v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &                    (3*L1*(v1n*(-1 - xi1) - v1np1k*(1 - xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4. + 
     &                    3*g*
     &                    Sin((th1n + th1np1k - 2*th2n - 2*th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - 2*th2n*xi1 + 
     &                   2*th2np1k*xi1)/2.))))/
     &             (L1*(8*m1 + 15*m2 - 
     &                  9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                 th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))**2)))/18.))/12. + (L1*L2*m2*v1np1k*Sin(th1np1k - th2np1k)*
     &     ((4*(3*L1*(m1 + 3*m2)*(v1n + v1np1k)**2*
     &             Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &            (9*L2*m2*(-v2n - v2np1k)**2*
     &               Sin(th1n + th1np1k - th2n - th2np1k))/4. + 
     &            9*g*(m1 + 2*m2)*
     &             Sin(th1n + th1np1k - th2n/2. - th2np1k/2.) - 
     &            3*g*(m1 + 6*m2)*Sin((th2n + th2np1k)/2.)))/
     &        (9.*L2*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((9*g*(m1 + 2*m2)*
     &                Sin(th1n + (th2n*(-1 - xi1))/2. + 
     &                  th1np1k*(1 - xi1) - (th2np1k*(1 - xi1))/2. + 
     &                  th1n*xi1) - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k + th2n*xi1 - th2np1k*xi1)/2.)
     & + 3*L1*(m1 + 3*m2)*(v1n + v1np1k + v1n*xi1 - v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 - xi1) - v2np1k*(1 - xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4.)/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                 th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))) + (3*L1*(m1 + 3*m2)*(v1n + v1np1k - v1n*xi1 + v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 + xi1) - v2np1k*(1 + xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4. - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k - th2n*xi1 + th2np1k*xi1)/2.)
     & + 9*g*(m1 + 2*m2)*Sin(th1n + (th2n*(-1 + xi1))/2. - th1n*xi1 + 
     &                  th1np1k*(1 + xi1) - (th2np1k*(1 + xi1))/2.))/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                 th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)
     & ))))/18.))/4.+(L2*m2*(2*L2*(v2n+v2np1k)+3*L1*v1n*Cos(th1n-th2n)+
     &       3*L1*v1np1k*Cos(th1np1k - th2np1k))*
     &     ((4*((-3*L1*(m1 + 3*m2)*(v1n + v1np1k)**2*
     &               Cos((th1n + th1np1k - th2n - th2np1k)/2.))/2. - 
     &            (9*L2*m2*(-v2n - v2np1k)**2*
     &               Cos(th1n + th1np1k - th2n - th2np1k))/4. - 
     &            (9*g*(m1 + 2*m2)*
     &               Cos(th1n + th1np1k - th2n/2. - th2np1k/2.))/2. - 
     &            (3*g*(m1 + 6*m2)*Cos((th2n + th2np1k)/2.))/2.))/
     &        (9.*L2*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (4*m2*Sin(th1n + th1np1k - th2n - th2np1k)*
     &          (3*L1*(m1 + 3*m2)*(v1n + v1np1k)**2*
     &             Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &            (9*L2*m2*(-v2n - v2np1k)**2*
     &               Sin(th1n + th1np1k - th2n - th2np1k))/4. + 
     &            9*g*(m1 + 2*m2)*
     &             Sin(th1n + th1np1k - th2n/2. - th2np1k/2.) - 
     &            3*g*(m1 + 6*m2)*Sin((th2n + th2np1k)/2.)))/
     &        (L2*(8*m1 + 15*m2 - 
     &             9*m2*Cos(th1n + th1np1k - th2n - th2np1k))**2) + 
     &       (5*(((9*g*(m1 + 2*m2)*(-1 + xi1)*
     &                  Cos(th1n + (th2n*(-1 - xi1))/2. + 
     &                   th1np1k*(1 - xi1) - (th2np1k*(1 - xi1))/2. + 
     &                   th1n*xi1))/2. - 
     &               (3*g*(m1 + 6*m2)*(1 - xi1)*
     &                  Cos((th2n + th2np1k + th2n*xi1 - th2np1k*xi1)/
     &                   2.))/2. + 
     &               (3*L1*(m1 + 3*m2)*(-1 + xi1)*
     &                  (v1n + v1np1k + v1n*xi1 - v1np1k*xi1)**2*
     &                  Cos((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.))/2. + 
     &               (9*L2*m2*(v2n*(-1 - xi1) - v2np1k*(1 - xi1))**2*
     &                  (-1 + xi1)*
     &                  Cos(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))/4.)/(L2*(8*m1 + 15*m2 - 9*m2*
     &                  Cos(th1n + th1np1k - th2n - th2np1k + th1n*xi1 -
     &                    th1np1k*xi1 - th2n*xi1 + th2np1k*xi1))) + 
     &            ((3*L1*(m1 + 3*m2)*(-1 - xi1)*
     &                  (v1n + v1np1k - v1n*xi1 + v1np1k*xi1)**2*
     &                  Cos((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.))/2. + 
     &               (9*L2*m2*(-1 - xi1)*
     &                  (v2n*(-1 + xi1) - v2np1k*(1 + xi1))**2*
     &                  Cos(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4. - 
     &               (3*g*(m1 + 6*m2)*(1 + xi1)*
     &                  Cos((th2n + th2np1k - th2n*xi1 + th2np1k*xi1)/
     &                   2.))/2. + 
     &               (9*g*(m1 + 2*m2)*(-1 - xi1)*
     &                  Cos(th1n + (th2n*(-1 + xi1))/2. - th1n*xi1 + 
     &                th1np1k*(1 + xi1) - (th2np1k*(1 + xi1))/2.))/2.)/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &               th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)
     & ))-(9*m2*(-1+xi1)*Sin(th1n + th1np1k - th2n - th2np1k + th1n*xi1-
     &                 th1np1k*xi1 - th2n*xi1 + th2np1k*xi1)*
     &               (9*g*(m1 + 2*m2)*
     &                  Sin(th1n + (th2n*(-1 - xi1))/2. + 
     &                   th1np1k*(1 - xi1) - (th2np1k*(1 - xi1))/2. + 
     &                   th1n*xi1) - 
     &                 3*g*(m1 + 6*m2)*
     &                  Sin((th2n + th2np1k + th2n*xi1 - th2np1k*xi1)/
     &                   2.) + 
     &                 3*L1*(m1 + 3*m2)*
     &                  (v1n + v1np1k + v1n*xi1 - v1np1k*xi1)**2*
     &                  Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &                 (9*L2*m2*(v2n*(-1 - xi1) - v2np1k*(1 - xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))/4.))/(L2*(8*m1 + 15*m2 - 9*m2*
     &                   Cos(th1n + th1np1k - th2n - th2np1k + 
     &                th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))**2) - (9*m2*(-1 - xi1)*Sin(th1n + th1np1k - th2n - th2np1k - 
     &                 th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)*
     &               (3*L1*(m1 + 3*m2)*
     &                  (v1n + v1np1k - v1n*xi1 + v1np1k*xi1)**2*
     &                  Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & )/2.) + (9*L2*m2*(v2n*(-1 + xi1) - v2np1k*(1 + xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))/4.-3*g*(m1+6*m2)*Sin((th2n+th2np1k - th2n*xi1 + th2np1k*xi1)/
     &                   2.) + 
     &                 9*g*(m1 + 2*m2)*
     &                  Sin(th1n + (th2n*(-1 + xi1))/2. - th1n*xi1 + 
     &                    th1np1k*(1 + xi1) - (th2np1k*(1 + xi1))/2.)))/
     &             (L2*(8*m1 + 15*m2 - 
     &                  9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &               th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)
     & )**2)))/18.))/12.))/
     &    (12.*((L1**2*(-2*g*(m1 + 2*m2)*Cos(th1n) + 
     &        2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &        L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &           v1np1k*v2np1k*
     &            (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k))))**2)/
     &   (16.*(th1n - th1np1k)**2) + 
     &  (L2**2*m2**2*(2*L2*(v2n + v2np1k) + 3*L1*v1n*Cos(th1n - th2n) +
     &        3*L1*v1np1k*Cos(th1np1k - th2np1k))**2)/144. + 
     &  (L1**2*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &        3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &           v2n*Cos(th1np1k - th2np1k)))**2)/144. + 
     &  (L2**2*m2**2*(L1*v1np1k*v2np1k*
     &         (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &        L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k))+
     &  2*g*(-Cos(th2n) + Cos(th2np1k)))**2)/(16.*(th2n - th2np1k)**2)))

       m33=1.0d0/tdelta-(2*m2*(-v1n-v1np1k)*
     & Sin(th1n+th1np1k-th2n-th2np1k))/
     &    (8*m1 + 15*m2 - 9*m2*Cos(th1n + th1np1k - th2n - th2np1k)) - 
     &   (5*((-9*m2*(-1 - xi1)*(v1n*(-1 + xi1) - v1np1k*(1 + xi1))*
     &     Sin(th1n + th1np1k - th2n - th2np1k - th1n*xi1 + th1np1k*xi1+
     &       th2n*xi1 - th2np1k*xi1))/
     &   (2.*(8*m1 + 15*m2 - 9*m2*
     &        Cos(th1n + th1np1k - th2n - th2np1k - th1n*xi1 + 
     &          th1np1k*xi1 + th2n*xi1 - th2np1k*xi1))) - 
     &  (9*m2*(v1n*(-1 - xi1) - v1np1k*(1 - xi1))*(-1 + xi1)*
     &     Sin(th1n + th1np1k - th2n - th2np1k + th1n*xi1 - th1np1k*xi1-
     &       th2n*xi1 + th2np1k*xi1))/
     &   (2.*(8*m1 + 15*m2 - 9*m2*
     &        Cos(th1n + th1np1k - th2n - th2np1k + th1n*xi1 - 
     &          th1np1k*xi1 - th2n*xi1 + th2np1k*xi1)))))/18. + 
     &   (L1*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &  3*L2*m2*(v2np1k*Cos(th1n - th2n) + v2n*Cos(th1np1k - th2np1k)))*
     &((L1*L2*m2*v2np1k*(v2n + v2np1k)*
     & (Cos(th1n - th2n) - Cos(th1n - th2np1k)))/(8.*(th2n - th2np1k))
     & + (L1*(-2*g*(m1 + 2*m2)*Cos(th1n) + 2*g*(m1 + 2*m2)*Cos(th1np1k)+
     &       L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &          v1np1k*v2np1k*(Cos(th1n - th2np1k) - 
     &             Cos(th1np1k - th2np1k)))))/(8.*(th1n - th1np1k)) + 
     &  (L1*L2*m2*(v1n + v1np1k)*v2np1k*
     &     (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k)))/
     &   (8.*(th1n - th1np1k)) + 
     &  (L2*m2*(2*L2*(v2n + v2np1k) + 3*L1*v1n*Cos(th1n - th2n) + 
     &       3*L1*v1np1k*Cos(th1np1k - th2np1k))*
     &     ((8*L1*(m1 + 3*m2)*(v1n + v1np1k)*
     &          Sin((th1n + th1np1k - th2n - th2np1k)/2.))/
     &        (3.*L2*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((6*L1*(m1 + 3*m2)*(1 + xi1)*
     &               (v1n + v1np1k - v1n*xi1 + v1np1k*xi1)*
     &               Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.))/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))) + (6*L1*(m1 + 3*m2)*(1 - xi1)*
     &               (v1n + v1np1k + v1n*xi1 - v1np1k*xi1)*
     &               Sin((th1n + th1np1k - th2n - th2np1k + th1n*xi1 - 
     &                   th1np1k*xi1 - th2n*xi1 + th2np1k*xi1)/2.))/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & )))))/18.))/12. + (L1*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &       3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &          v2n*Cos(th1np1k - th2np1k)))*
     &     ((2*m2*(-v1n - v1np1k)*Sin(th1n + th1np1k - th2n - th2np1k))/
     &        (8*m1 + 15*m2 - 9*m2*Cos(th1n + th1np1k - th2n - th2np1k))
     & + (5*((-9*m2*(-1 - xi1)*(v1n*(-1 + xi1) - v1np1k*(1 + xi1))*
     &               Sin(th1n + th1np1k - th2n - th2np1k - th1n*xi1 + 
     &                 th1np1k*xi1 + th2n*xi1 - th2np1k*xi1))/
     &             (2.*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))) - (9*m2*(v1n*(-1 - xi1) - v1np1k*(1 - xi1))*(-1 + xi1)*
     &               Sin(th1n + th1np1k - th2n - th2np1k + th1n*xi1 - 
     &                 th1np1k*xi1 - th2n*xi1 + th2np1k*xi1))/
     &             (2.*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                  th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & )))))/18.))/12. + (L1**2*(m1 + 3*m2)*
     &     ((-4*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k)/2.) + 
     &            m2*(3*g*Sin((th1n + th1np1k - 2*th2n - 2*th2np1k)/
     &                  2.) + 
     &               L2*(v2n + v2np1k)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &               (3*L1*(-v1n - v1np1k)**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k))/4.)))/
     &        (3.*L1*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((-3*(g*(4*m1 + 5*m2)*
     &                  Sin((th1n + th1np1k - th1n*xi1 + th1np1k*xi1)/
     &                   2.) + 
     &                 m2*(3*g*
     &                   Sin((th1n + th1np1k - 2*th2n - 2*th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + 2*th2n*xi1 - 
     &                  2*th2np1k*xi1)/2.) + 
     &                   L2*(v2n + v2np1k - v2n*xi1 + v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                  th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 + xi1) - v1np1k*(1 + xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4.)))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))) - (3*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k + th1n*xi1 - 
     &                   th1np1k*xi1)/2.) + 
     &                 m2*(L2*(v2n + v2np1k + v2n*xi1 - v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 - xi1) - v1np1k*(1 - xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4. + 
     &                   3*g*Sin((th1n + th1np1k - 2*th2n - 
     &                   2*th2np1k + th1n*xi1 - th1np1k*xi1 - 
     &                   2*th2n*xi1 + 2*th2np1k*xi1)/2.))))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & )))))/18.))/6. + (L1*L2*m2*Cos(th1np1k - th2np1k)*
     &     ((4*(3*L1*(m1 + 3*m2)*(v1n + v1np1k)**2*
     &             Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &            (9*L2*m2*(-v2n - v2np1k)**2*
     &               Sin(th1n + th1np1k - th2n - th2np1k))/4. + 
     &            9*g*(m1 + 2*m2)*
     &             Sin(th1n + th1np1k - th2n/2. - th2np1k/2.) - 
     &            3*g*(m1 + 6*m2)*Sin((th2n + th2np1k)/2.)))/
     &        (9.*L2*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((9*g*(m1 + 2*m2)*
     &                Sin(th1n + (th2n*(-1 - xi1))/2. + 
     &                  th1np1k*(1 - xi1) - (th2np1k*(1 - xi1))/2. + 
     &                  th1n*xi1) - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k + th2n*xi1 - th2np1k*xi1)/2.)
     & + 3*L1*(m1 + 3*m2)*(v1n + v1np1k + v1n*xi1 - v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 - xi1) - v2np1k*(1 - xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4.)/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                  th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))) + (3*L1*(m1 + 3*m2)*(v1n + v1np1k - v1n*xi1 + v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 + xi1) - v2np1k*(1 + xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4. - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k - th2n*xi1 + th2np1k*xi1)/2.)
     & + 9*g*(m1 + 2*m2)*Sin(th1n + (th2n*(-1 + xi1))/2. - th1n*xi1 + 
     &                  th1np1k*(1 + xi1) - (th2np1k*(1 + xi1))/2.))/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &               th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)
     & ))))/18.))/4.))/
     &    (12.*((L1**2*(-2*g*(m1 + 2*m2)*Cos(th1n) + 
     &        2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &        L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &           v1np1k*v2np1k*
     &            (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k))))**2)/
     &   (16.*(th1n - th1np1k)**2) + 
     &  (L2**2*m2**2*(2*L2*(v2n + v2np1k) + 3*L1*v1n*Cos(th1n - th2n) + 
     &        3*L1*v1np1k*Cos(th1np1k - th2np1k))**2)/144. + 
     &  (L1**2*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &        3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &           v2n*Cos(th1np1k - th2np1k)))**2)/144. + 
     &  (L2**2*m2**2*(L1*v1np1k*v2np1k*
     &         (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &        L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k))+
     &  2*g*(-Cos(th2n) + Cos(th2np1k)))**2)/(16.*(th2n - th2np1k)**2))
     & ) - (L1*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &  3*L2*m2*(v2np1k*Cos(th1n - th2n) + v2n*Cos(th1np1k - th2np1k)))*
     & ((L1**2*L2*m2*v2np1k*(-2*g*(m1 + 2*m2)*Cos(th1n) + 
     &       2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &       L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &          v1np1k*v2np1k*
     &           (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k))))*
     &     (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k)))/
     &   (8.*(th1n - th1np1k)**2) + 
     &  (L1*L2**2*m2**2*Cos(th1np1k - th2np1k)*
     &     (2*L2*(v2n + v2np1k) + 3*L1*v1n*Cos(th1n - th2n) + 
     &       3*L1*v1np1k*Cos(th1np1k - th2np1k)))/24. + 
     &  (L1**3*(m1 + 3*m2)*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &       3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &          v2n*Cos(th1np1k - th2np1k))))/36. + 
     &  (L1*L2**2*m2**2*v2np1k*(Cos(th1n - th2n) - Cos(th1n - th2np1k))*
     &     (L1*v1np1k*v2np1k*(Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &       L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k)) +
     &     2*g*(-Cos(th2n) + Cos(th2np1k))))/(8.*(th2n - th2np1k)**2))*
     & ((L1*(v1n + v1np1k)*(-2*g*(m1 + 2*m2)*Cos(th1n) + 
     &       2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &       L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &          v1np1k*v2np1k*(Cos(th1n - th2np1k) - 
     &             Cos(th1np1k - th2np1k)))))/(8.*(th1n - th1np1k)) + 
     &  (L2*m2*(v2n + v2np1k)*(L1*v1np1k*v2np1k*
     &        (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &       L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k)) +
     &       2*g*(-Cos(th2n) + Cos(th2np1k))))/(8.*(th2n - th2np1k)) + 
     &  (L1*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &       3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &          v2n*Cos(th1np1k - th2np1k)))*
     &     ((-4*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k)/2.) + 
     &            m2*(3*g*Sin((th1n + th1np1k - 2*th2n - 2*th2np1k)/
     &                  2.) + 
     &               L2*(v2n + v2np1k)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &               (3*L1*(-v1n - v1np1k)**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k))/4.)))/
     &        (3.*L1*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((-3*(g*(4*m1 + 5*m2)*
     &                  Sin((th1n + th1np1k - th1n*xi1 + th1np1k*xi1)/
     &                   2.) + 
     &                 m2*(3*g*
     &                   Sin((th1n + th1np1k - 2*th2n - 2*th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + 2*th2n*xi1 - 
     &                  2*th2np1k*xi1)/2.) + 
     &                   L2*(v2n + v2np1k - v2n*xi1 + v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                  th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 + xi1) - v1np1k*(1 + xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4.)))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))) - (3*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k + th1n*xi1 - 
     &                   th1np1k*xi1)/2.) + 
     &                 m2*(L2*(v2n + v2np1k + v2n*xi1 - v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 - xi1) - v1np1k*(1 - xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4. + 
     &                   3*g*Sin((th1n + th1np1k - 2*th2n - 
     &                   2*th2np1k + th1n*xi1 - th1np1k*xi1 - 
     &                   2*th2n*xi1 + 2*th2np1k*xi1)/2.))))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                 th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & )))))/18.))/12. + (L2*m2*(2*L2*(v2n + v2np1k) + 
     &       3*L1*v1n*Cos(th1n - th2n) + 
     &       3*L1*v1np1k*Cos(th1np1k - th2np1k))*
     &     ((4*(3*L1*(m1 + 3*m2)*(v1n + v1np1k)**2*
     &             Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &            (9*L2*m2*(-v2n - v2np1k)**2*
     &               Sin(th1n + th1np1k - th2n - th2np1k))/4. + 
     &            9*g*(m1 + 2*m2)*
     &             Sin(th1n + th1np1k - th2n/2. - th2np1k/2.) - 
     &            3*g*(m1 + 6*m2)*Sin((th2n + th2np1k)/2.)))/
     &        (9.*L2*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((9*g*(m1 + 2*m2)*
     &                Sin(th1n + (th2n*(-1 - xi1))/2. + 
     &                  th1np1k*(1 - xi1) - (th2np1k*(1 - xi1))/2. + 
     &                  th1n*xi1) - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k + th2n*xi1 - th2np1k*xi1)/2.)
     & + 3*L1*(m1 + 3*m2)*(v1n + v1np1k + v1n*xi1 - v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 - xi1) - v2np1k*(1 - xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4.)/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                 th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))) + (3*L1*(m1 + 3*m2)*(v1n + v1np1k - v1n*xi1 + v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 + xi1) - v2np1k*(1 + xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4. - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k - th2n*xi1 + th2np1k*xi1)/2.)
     & + 9*g*(m1 + 2*m2)*Sin(th1n + (th2n*(-1 + xi1))/2. - th1n*xi1 + 
     &                  th1np1k*(1 + xi1) - (th2np1k*(1 + xi1))/2.))/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                 th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)
     & ))))/18.))/12.))/
     &    (12.*((L1**2*(-2*g*(m1 + 2*m2)*Cos(th1n) + 
     &         2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &         L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) +
     &            v1np1k*v2np1k*
     &             (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k))))**2)/
     &    (16.*(th1n - th1np1k)**2) + 
     &   (L2**2*m2**2*(2*L2*(v2n + v2np1k) + 3*L1*v1n*Cos(th1n - th2n) +
     &         3*L1*v1np1k*Cos(th1np1k - th2np1k))**2)/144. + 
     &   (L1**2*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &         3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &            v2n*Cos(th1np1k - th2np1k)))**2)/144. + 
     &   (L2**2*m2**2*(L1*v1np1k*v2np1k*
     &          (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &         L1*v1n*v2n*(Cos(th1np1k - th2n) - 
     &            Cos(th1np1k - th2np1k)) + 
     &         2*g*(-Cos(th2n) + Cos(th2np1k)))**2)/
     &    (16.*(th2n - th2np1k)**2))**2) + 
     &   (L1**2*(m1 + 3*m2)*((L1*(v1n + v1np1k)*
     &     (-2*g*(m1 + 2*m2)*Cos(th1n) + 2*g*(m1 + 2*m2)*Cos(th1np1k) +
     &       L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &          v1np1k*v2np1k*(Cos(th1n - th2np1k) - 
     &             Cos(th1np1k - th2np1k)))))/(8.*(th1n - th1np1k)) + 
     &  (L2*m2*(v2n + v2np1k)*(L1*v1np1k*v2np1k*
     &        (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &       L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k)) +
     &       2*g*(-Cos(th2n) + Cos(th2np1k))))/(8.*(th2n - th2np1k)) + 
     &  (L1*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &       3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &          v2n*Cos(th1np1k - th2np1k)))*
     &     ((-4*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k)/2.) + 
     &            m2*(3*g*Sin((th1n + th1np1k - 2*th2n - 2*th2np1k)/
     &                  2.) + 
     &               L2*(v2n + v2np1k)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &               (3*L1*(-v1n - v1np1k)**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k))/4.)))/
     &        (3.*L1*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((-3*(g*(4*m1 + 5*m2)*
     &                  Sin((th1n + th1np1k - th1n*xi1 + th1np1k*xi1)/
     &                   2.) + 
     &                 m2*(3*g*
     &                   Sin((th1n + th1np1k - 2*th2n - 2*th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + 2*th2n*xi1 - 
     &                   2*th2np1k*xi1)/2.) + 
     &                   L2*(v2n + v2np1k - v2n*xi1 + v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 + xi1) - v1np1k*(1 + xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4.)))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))) - (3*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k + th1n*xi1 - 
     &                   th1np1k*xi1)/2.) + 
     &                 m2*(L2*(v2n + v2np1k + v2n*xi1 - v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &                    (3*L1*(v1n*(-1 - xi1) - v1np1k*(1 - xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4. + 
     &                    3*g*
     &                    Sin((th1n + th1np1k - 2*th2n - 2*th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - 2*th2n*xi1 + 
     &                   2*th2np1k*xi1)/2.))))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                 th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1)
     & ))))/18.))/12. + (L2*m2*(2*L2*(v2n + v2np1k) + 
     & 3*L1*v1n*Cos(th1n - th2n) + 3*L1*v1np1k*Cos(th1np1k - th2np1k)
     & )*((4*(3*L1*(m1 + 3*m2)*(v1n + v1np1k)**2*
     &             Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &            (9*L2*m2*(-v2n - v2np1k)**2*
     &               Sin(th1n + th1np1k - th2n - th2np1k))/4. + 
     &            9*g*(m1 + 2*m2)*
     &             Sin(th1n + th1np1k - th2n/2. - th2np1k/2.) - 
     &            3*g*(m1 + 6*m2)*Sin((th2n + th2np1k)/2.)))/
     &        (9.*L2*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((9*g*(m1 + 2*m2)*
     &                Sin(th1n + (th2n*(-1 - xi1))/2. + 
     &                  th1np1k*(1 - xi1) - (th2np1k*(1 - xi1))/2. + 
     &                  th1n*xi1) - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k + th2n*xi1 - th2np1k*xi1)/2.)
     & + 3*L1*(m1 + 3*m2)*(v1n + v1np1k + v1n*xi1 - v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 - xi1) - v2np1k*(1 - xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))/4.)/(L2*(8*m1 + 15*m2 - 9*m2*
     &                  Cos(th1n + th1np1k - th2n - th2np1k + th1n*xi1 -
     &                    th1np1k*xi1 - th2n*xi1 + th2np1k*xi1))) + 
     &            (3*L1*(m1 + 3*m2)*
     &                (v1n + v1np1k - v1n*xi1 + v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k - th1n*xi1 + 
     &                   th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 + xi1) - v2np1k*(1 + xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))/4.-3*g*(m1+6*m2)*Sin((th2n + th2np1k - th2n*xi1+th2np1k*xi1)/
     &                  2.) + 9*g*(m1 + 2*m2)*
     &                Sin(th1n + (th2n*(-1 + xi1))/2. - th1n*xi1 + 
     &                  th1np1k*(1 + xi1) - (th2np1k*(1 + xi1))/2.))/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &           th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1))
     & )))/18.))/12.))/
     &    (6.*((L1**2*(-2*g*(m1 + 2*m2)*Cos(th1n) + 
     &        2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &        L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &           v1np1k*v2np1k*
     &            (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k))))**2)/
     &   (16.*(th1n - th1np1k)**2) + 
     &  (L2**2*m2**2*(2*L2*(v2n + v2np1k) + 3*L1*v1n*Cos(th1n - th2n) +
     &        3*L1*v1np1k*Cos(th1np1k - th2np1k))**2)/144. + 
     &  (L1**2*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &        3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &           v2n*Cos(th1np1k - th2np1k)))**2)/144. + 
     &  (L2**2*m2**2*(L1*v1np1k*v2np1k*
     &   (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &  L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k)) + 
     & 2*g*(-Cos(th2n) + Cos(th2np1k)))**2)/(16.*(th2n - th2np1k)**2)))

      m34=(8*L2*m2*(v2n + v2np1k)*Sin((th1n+th1np1k-th2n-th2np1k)/2.))/
     & (3.*L1*(8*m1+15*m2-9*m2*Cos(th1n + th1np1k - th2n - th2np1k)))-
     & (5*((-6*L2*m2*(1 + xi1)*(v2n + v2np1k - v2n*xi1 + v2np1k*xi1)*
     &     Sin((th1n + th1np1k - th2n - th2np1k - th1n*xi1 + 
     &         th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)/2.))/
     &   (L1*(8*m1 + 15*m2 - 9*m2*
     &        Cos(th1n + th1np1k - th2n - th2np1k - th1n*xi1 + 
     &          th1np1k*xi1 + th2n*xi1 - th2np1k*xi1))) - 
     &  (6*L2*m2*(1 - xi1)*(v2n + v2np1k + v2n*xi1 - v2np1k*xi1)*
     &     Sin((th1n + th1np1k - th2n - th2np1k + th1n*xi1 - 
     &         th1np1k*xi1 - th2n*xi1 + th2np1k*xi1)/2.))/
     &   (L1*(8*m1 + 15*m2 - 9*m2*
     &        Cos(th1n + th1np1k - th2n - th2np1k + th1n*xi1 - 
     &          th1np1k*xi1 - th2n*xi1 + th2np1k*xi1)))))/18. + 
     & (L1*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &  3*L2*m2*(v2np1k*Cos(th1n - th2n) + v2n*Cos(th1np1k - th2np1k)))*
     & ((L1*L2*m2*v1np1k*(v2n + v2np1k)*
     &   (Cos(th1n - th2n) - Cos(th1n - th2np1k)))/(8.*(th2n - th2np1k))
     & + (L1*L2*m2*v1np1k*(v1n + v1np1k)*
     &     (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k)))/
     &   (8.*(th1n - th1np1k)) + 
     &  (L2*m2*(L1*v1np1k*v2np1k*
     &        (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &       L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k))+
     &       2*g*(-Cos(th2n) + Cos(th2np1k))))/(8.*(th2n - th2np1k)) + 
     &  (L1*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &       3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &          v2n*Cos(th1np1k - th2np1k)))*
     &     ((-8*L2*m2*(v2n + v2np1k)*
     &          Sin((th1n + th1np1k - th2n - th2np1k)/2.))/
     &        (3.*L1*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((-6*L2*m2*(1 + xi1)*
     &               (v2n + v2np1k - v2n*xi1 + v2np1k*xi1)*
     &               Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))) - (6*L2*m2*(1 - xi1)*(v2n + v2np1k + v2n*xi1 - v2np1k*xi1)*
     &               Sin((th1n + th1np1k - th2n - th2np1k + th1n*xi1 - 
     &                   th1np1k*xi1 - th2n*xi1 + th2np1k*xi1)/2.))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                 th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & )))))/18.))/12. + (L2*m2*(2*L2*(v2n + v2np1k) + 
     &       3*L1*v1n*Cos(th1n - th2n) + 
     &       3*L1*v1np1k*Cos(th1np1k - th2np1k))*
     &   ((-2*m2*(-v2n - v2np1k)*Sin(th1n + th1np1k - th2n - th2np1k))/
     &        (8*m1 + 15*m2 - 9*m2*Cos(th1n + th1np1k - th2n - th2np1k))
     & + (5*((9*m2*(-1 - xi1)*(v2n*(-1 + xi1) - v2np1k*(1 + xi1))*
     &               Sin(th1n + th1np1k - th2n - th2np1k - th1n*xi1 + 
     &                 th1np1k*xi1 + th2n*xi1 - th2np1k*xi1))/
     &             (2.*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))) + (9*m2*(v2n*(-1 - xi1) - v2np1k*(1 - xi1))*(-1 + xi1)*
     &               Sin(th1n + th1np1k - th2n - th2np1k + th1n*xi1 - 
     &                 th1np1k*xi1 - th2n*xi1 + th2np1k*xi1))/
     &             (2.*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                 th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & )))))/18.))/12. + (L1*L2*m2*Cos(th1n - th2n)*
     &     ((-4*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k)/2.) + 
     &            m2*(3*g*Sin((th1n + th1np1k - 2*th2n - 2*th2np1k)/
     &                  2.) + 
     &               L2*(v2n + v2np1k)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &               (3*L1*(-v1n - v1np1k)**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k))/4.)))/
     &        (3.*L1*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((-3*(g*(4*m1 + 5*m2)*
     &                  Sin((th1n + th1np1k - th1n*xi1 + th1np1k*xi1)/
     &                   2.) + 
     &                 m2*(3*g*
     &                   Sin((th1n + th1np1k - 2*th2n - 2*th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + 2*th2n*xi1 - 
     &                  2*th2np1k*xi1)/2.) + 
     &                   L2*(v2n + v2np1k - v2n*xi1 + v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                  th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 + xi1) - v1np1k*(1 + xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4.)))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))) - (3*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k + th1n*xi1 - 
     &                   th1np1k*xi1)/2.) + 
     &                 m2*(L2*(v2n + v2np1k + v2n*xi1 - v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 - xi1) - v1np1k*(1 - xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4. + 
     &                   3*g*Sin((th1n + th1np1k - 2*th2n - 
     &                   2*th2np1k + th1n*xi1 - th1np1k*xi1 - 
     &                   2*th2n*xi1 + 2*th2np1k*xi1)/2.))))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                 th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & )))))/18.))/4. + (L2**2*m2*((4*
     &          (3*L1*(m1 + 3*m2)*(v1n + v1np1k)**2*
     &             Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &            (9*L2*m2*(-v2n - v2np1k)**2*
     &               Sin(th1n + th1np1k - th2n - th2np1k))/4. + 
     &            9*g*(m1 + 2*m2)*
     &             Sin(th1n + th1np1k - th2n/2. - th2np1k/2.) - 
     &            3*g*(m1 + 6*m2)*Sin((th2n + th2np1k)/2.)))/
     &        (9.*L2*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((9*g*(m1 + 2*m2)*
     &                Sin(th1n + (th2n*(-1 - xi1))/2. + 
     &                  th1np1k*(1 - xi1) - (th2np1k*(1 - xi1))/2. + 
     &                  th1n*xi1) - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k + th2n*xi1 - th2np1k*xi1)/2.)
     & + 3*L1*(m1 + 3*m2)*(v1n + v1np1k + v1n*xi1 - v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 - xi1) - v2np1k*(1 - xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4.)/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                  th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))) + (3*L1*(m1 + 3*m2)*(v1n + v1np1k - v1n*xi1 + v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 + xi1) - v2np1k*(1 + xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4. - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k - th2n*xi1 + th2np1k*xi1)/2.)
     & + 9*g*(m1 + 2*m2)*Sin(th1n + (th2n*(-1 + xi1))/2. - th1n*xi1 + 
     &                  th1np1k*(1 + xi1) - (th2np1k*(1 + xi1))/2.))/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &            th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)
     & ))))/18.))/6.))/
     &    (12.*((L1**2*(-2*g*(m1 + 2*m2)*Cos(th1n) + 
     &        2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &        L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) +
     &           v1np1k*v2np1k*
     &            (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k))))**2)/
     &   (16.*(th1n - th1np1k)**2) + 
     &  (L2**2*m2**2*(2*L2*(v2n + v2np1k) + 3*L1*v1n*Cos(th1n - th2n) +
     &        3*L1*v1np1k*Cos(th1np1k - th2np1k))**2)/144. + 
     &  (L1**2*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &        3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &           v2n*Cos(th1np1k - th2np1k)))**2)/144. + 
     &  (L2**2*m2**2*(L1*v1np1k*v2np1k*
     &         (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &        L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k))+
     &  2*g*(-Cos(th2n) + Cos(th2np1k)))**2)/(16.*(th2n - th2np1k)**2))
     & ) - (L1*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &  3*L2*m2*(v2np1k*Cos(th1n - th2n) + v2n*Cos(th1np1k - th2np1k)))*
     &((L1**2*L2*m2*v1np1k*(-2*g*(m1 + 2*m2)*Cos(th1n) + 
     &       2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &       L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &          v1np1k*v2np1k*
     &           (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k))))*
     &     (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k)))/
     &   (8.*(th1n - th1np1k)**2) + 
     &  (L2**3*m2**2*(2*L2*(v2n + v2np1k) + 3*L1*v1n*Cos(th1n - th2n) + 
     &       3*L1*v1np1k*Cos(th1np1k - th2np1k)))/36. + 
     &  (L1**2*L2*m2*Cos(th1n - th2n)*
     &     (2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &       3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &          v2n*Cos(th1np1k - th2np1k))))/24. + 
     &  (L1*L2**2*m2**2*v1np1k*(Cos(th1n - th2n) - Cos(th1n - th2np1k))*
     &     (L1*v1np1k*v2np1k*(Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &       L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k)) + 
     &     2*g*(-Cos(th2n) + Cos(th2np1k))))/(8.*(th2n - th2np1k)**2))*
     & ((L1*(v1n + v1np1k)*(-2*g*(m1 + 2*m2)*Cos(th1n) + 
     &       2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &       L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &          v1np1k*v2np1k*(Cos(th1n - th2np1k) - 
     &             Cos(th1np1k - th2np1k)))))/(8.*(th1n - th1np1k)) + 
     &  (L2*m2*(v2n + v2np1k)*(L1*v1np1k*v2np1k*
     &        (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &       L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k)) +
     &       2*g*(-Cos(th2n) + Cos(th2np1k))))/(8.*(th2n - th2np1k)) + 
     &  (L1*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &       3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &          v2n*Cos(th1np1k - th2np1k)))*
     &     ((-4*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k)/2.) + 
     &            m2*(3*g*Sin((th1n + th1np1k - 2*th2n - 2*th2np1k)/
     &                  2.) + 
     &               L2*(v2n + v2np1k)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &               (3*L1*(-v1n - v1np1k)**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k))/4.)))/
     &        (3.*L1*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((-3*(g*(4*m1 + 5*m2)*
     &                  Sin((th1n + th1np1k - th1n*xi1 + th1np1k*xi1)/
     &                   2.) + 
     &                 m2*(3*g*
     &                   Sin((th1n + th1np1k - 2*th2n - 2*th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + 2*th2n*xi1 - 
     &                  2*th2np1k*xi1)/2.) + 
     &                   L2*(v2n + v2np1k - v2n*xi1 + v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                  th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 + xi1) - v1np1k*(1 + xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4.)))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))) - (3*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k + th1n*xi1 - 
     &                   th1np1k*xi1)/2.) + 
     &                 m2*(L2*(v2n + v2np1k + v2n*xi1 - v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 - xi1) - v1np1k*(1 - xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4. + 
     &                   3*g*Sin((th1n + th1np1k - 2*th2n - 
     &                   2*th2np1k + th1n*xi1 - th1np1k*xi1 - 
     &                   2*th2n*xi1 + 2*th2np1k*xi1)/2.))))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                  th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & )))))/18.))/12. + (L2*m2*(2*L2*(v2n + v2np1k) + 
     &       3*L1*v1n*Cos(th1n - th2n) + 
     &       3*L1*v1np1k*Cos(th1np1k - th2np1k))*
     &     ((4*(3*L1*(m1 + 3*m2)*(v1n + v1np1k)**2*
     &             Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &            (9*L2*m2*(-v2n - v2np1k)**2*
     &               Sin(th1n + th1np1k - th2n - th2np1k))/4. + 
     &            9*g*(m1 + 2*m2)*
     &             Sin(th1n + th1np1k - th2n/2. - th2np1k/2.) - 
     &            3*g*(m1 + 6*m2)*Sin((th2n + th2np1k)/2.)))/
     &        (9.*L2*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((9*g*(m1 + 2*m2)*
     &                Sin(th1n + (th2n*(-1 - xi1))/2. + 
     &                  th1np1k*(1 - xi1) - (th2np1k*(1 - xi1))/2. + 
     &                  th1n*xi1) - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k + th2n*xi1 - th2np1k*xi1)/2.)
     & + 3*L1*(m1 + 3*m2)*(v1n + v1np1k + v1n*xi1 - v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 - xi1) - v2np1k*(1 - xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4.)/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                  th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))) + (3*L1*(m1 + 3*m2)*(v1n + v1np1k - v1n*xi1 + v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 + xi1) - v2np1k*(1 + xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4. - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k - th2n*xi1 + th2np1k*xi1)/2.)
     & + 9*g*(m1 + 2*m2)*Sin(th1n + (th2n*(-1 + xi1))/2. - th1n*xi1 + 
     &                  th1np1k*(1 + xi1) - (th2np1k*(1 + xi1))/2.))/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &           th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)
     & ))))/18.))/12.))/
     &    (12.*((L1**2*(-2*g*(m1 + 2*m2)*Cos(th1n) + 
     &         2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &         L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) +
     &            v1np1k*v2np1k*
     &             (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k))))**2)/
     &    (16.*(th1n - th1np1k)**2) + 
     &   (L2**2*m2**2*(2*L2*(v2n + v2np1k) + 3*L1*v1n*Cos(th1n - th2n) +
     &         3*L1*v1np1k*Cos(th1np1k - th2np1k))**2)/144. + 
     &   (L1**2*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &         3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &            v2n*Cos(th1np1k - th2np1k)))**2)/144. + 
     &   (L2**2*m2**2*(L1*v1np1k*v2np1k*
     &          (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &         L1*v1n*v2n*(Cos(th1np1k - th2n) - 
     &            Cos(th1np1k - th2np1k)) + 
     &         2*g*(-Cos(th2n) + Cos(th2np1k)))**2)/
     &    (16.*(th2n - th2np1k)**2))**2) + 
     &   (L1*L2*m2*Cos(th1n - th2n)*((L1*(v1n + v1np1k)*
     &     (-2*g*(m1 + 2*m2)*Cos(th1n) + 2*g*(m1 + 2*m2)*Cos(th1np1k) +
     &       L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &          v1np1k*v2np1k*(Cos(th1n - th2np1k) - 
     &             Cos(th1np1k - th2np1k)))))/(8.*(th1n - th1np1k)) + 
     &  (L2*m2*(v2n + v2np1k)*(L1*v1np1k*v2np1k*
     &        (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &       L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k)) +
     &       2*g*(-Cos(th2n) + Cos(th2np1k))))/(8.*(th2n - th2np1k)) + 
     &  (L1*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &       3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &          v2n*Cos(th1np1k - th2np1k)))*
     &     ((-4*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k)/2.) + 
     &            m2*(3*g*Sin((th1n + th1np1k - 2*th2n - 2*th2np1k)/
     &                  2.) + 
     &               L2*(v2n + v2np1k)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &               (3*L1*(-v1n - v1np1k)**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k))/4.)))/
     &        (3.*L1*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((-3*(g*(4*m1 + 5*m2)*
     &                  Sin((th1n + th1np1k - th1n*xi1 + th1np1k*xi1)/
     &                   2.) + 
     &                 m2*(3*g*
     &                   Sin((th1n + th1np1k - 2*th2n - 2*th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + 2*th2n*xi1 - 
     &                   2*th2np1k*xi1)/2.) + 
     &                   L2*(v2n + v2np1k - v2n*xi1 + v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 + xi1) - v1np1k*(1 + xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4.)))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))) - (3*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k + th1n*xi1 - 
     &                   th1np1k*xi1)/2.) + 
     &                 m2*(L2*(v2n + v2np1k + v2n*xi1 - v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &                    (3*L1*(v1n*(-1 - xi1) - v1np1k*(1 - xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4. + 
     &                    3*g*
     &                    Sin((th1n + th1np1k - 2*th2n - 2*th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - 2*th2n*xi1 + 
     &                   2*th2np1k*xi1)/2.))))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &               th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1)
     & ))))/18.))/12. + (L2*m2*(2*L2*(v2n + v2np1k) + 
     &  3*L1*v1n*Cos(th1n - th2n) + 3*L1*v1np1k*Cos(th1np1k - th2np1k)
     & )*((4*(3*L1*(m1 + 3*m2)*(v1n + v1np1k)**2*
     &             Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &            (9*L2*m2*(-v2n - v2np1k)**2*
     &               Sin(th1n + th1np1k - th2n - th2np1k))/4. + 
     &            9*g*(m1 + 2*m2)*
     &             Sin(th1n + th1np1k - th2n/2. - th2np1k/2.) - 
     &            3*g*(m1 + 6*m2)*Sin((th2n + th2np1k)/2.)))/
     &        (9.*L2*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((9*g*(m1 + 2*m2)*
     &                Sin(th1n + (th2n*(-1 - xi1))/2. + 
     &                  th1np1k*(1 - xi1) - (th2np1k*(1 - xi1))/2. + 
     &                  th1n*xi1) - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k + th2n*xi1 - th2np1k*xi1)/2.)
     & + 3*L1*(m1 + 3*m2)*(v1n + v1np1k + v1n*xi1 - v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 - xi1) - v2np1k*(1 - xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))/4.)/(L2*(8*m1 + 15*m2 - 9*m2*
     &                  Cos(th1n + th1np1k - th2n - th2np1k + th1n*xi1 -
     &                    th1np1k*xi1 - th2n*xi1 + th2np1k*xi1))) + 
     &            (3*L1*(m1 + 3*m2)*
     &                (v1n + v1np1k - v1n*xi1 + v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k - th1n*xi1 +
     &                   th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 + xi1) - v2np1k*(1 + xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))/4.-3*g*(m1+6*m2)*Sin((th2n+th2np1k - th2n*xi1 + th2np1k*xi1)/
     &                  2.) + 9*g*(m1 + 2*m2)*
     &                Sin(th1n + (th2n*(-1 + xi1))/2. - th1n*xi1 + 
     &                  th1np1k*(1 + xi1) - (th2np1k*(1 + xi1))/2.))/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &       th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1))
     & )))/18.))/12.))/
     &    (4.*((L1**2*(-2*g*(m1 + 2*m2)*Cos(th1n) + 
     &        2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &        L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) +
     &           v1np1k*v2np1k*
     &            (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k))))**2)/
     &   (16.*(th1n - th1np1k)**2) + 
     &  (L2**2*m2**2*(2*L2*(v2n + v2np1k) + 3*L1*v1n*Cos(th1n - th2n) +
     &        3*L1*v1np1k*Cos(th1np1k - th2np1k))**2)/144. + 
     &  (L1**2*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &        3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &           v2n*Cos(th1np1k - th2np1k)))**2)/144. + 
     &  (L2**2*m2**2*(L1*v1np1k*v2np1k*
     &         (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &        L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k))+
     & 2*g*(-Cos(th2n) + Cos(th2np1k)))**2)/(16.*(th2n - th2np1k)**2)))

      m41= (-4*((3*L1*(m1 + 3*m2)*(v1n + v1np1k)**2*
     &     Cos((th1n + th1np1k - th2n - th2np1k)/2.))/2. + 
     & (9*L2*m2*(-v2n - v2np1k)**2*Cos(th1n + th1np1k-th2n-th2np1k))/
     &  4. + 9*g*(m1 + 2*m2)*Cos(th1n + th1np1k - th2n/2.-th2np1k/2.)))/
     & (9.*L2*(8*m1 + 15*m2 - 9*m2*Cos(th1n + th1np1k-th2n-th2np1k)))+
     & (4*m2*Sin(th1n + th1np1k - th2n - th2np1k)*
     & (3*L1*(m1 + 3*m2)*(v1n + v1np1k)**2*
     &   Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     & (9*L2*m2*(-v2n-v2np1k)**2*Sin(th1n + th1np1k - th2n - th2np1k))/
     &   4. + 9*g*(m1 + 2*m2)*Sin(th1n + th1np1k-th2n/2.-th2np1k/2.)-
     &  3*g*(m1 + 6*m2)*Sin((th2n + th2np1k)/2.)))/
     & (L2*(8*m1 + 15*m2 - 9*m2*Cos(th1n + th1np1k-th2n-th2np1k))**2)-
     & (5*((9*g*(m1 + 2*m2)*(1 - xi1)*
     &      Cos(th1n + (th2n*(-1 - xi1))/2. + th1np1k*(1 - xi1) - 
     &        (th2np1k*(1 - xi1))/2. + th1n*xi1) + 
     &     (3*L1*(m1 + 3*m2)*(1 - xi1)*
     &        (v1n + v1np1k + v1n*xi1 - v1np1k*xi1)**2*
     &        Cos((th1n + th1np1k - th2n - th2np1k + th1n*xi1 - 
     &            th1np1k*xi1 - th2n*xi1 + th2np1k*xi1)/2.))/2. + 
     &     (9*L2*m2*(v2n*(-1 - xi1) - v2np1k*(1 - xi1))**2*(1 - xi1)*
     &        Cos(th1n + th1np1k - th2n - th2np1k + th1n*xi1 - 
     &          th1np1k*xi1 - th2n*xi1 + th2np1k*xi1))/4.)/
     &   (L2*(8*m1 + 15*m2 - 9*m2*
     &        Cos(th1n + th1np1k - th2n - th2np1k + th1n*xi1 - 
     &          th1np1k*xi1 - th2n*xi1 + th2np1k*xi1))) + 
     &  ((3*L1*(m1 + 3*m2)*(1 + xi1)*
     &        (v1n + v1np1k - v1n*xi1 + v1np1k*xi1)**2*
     &        Cos((th1n + th1np1k - th2n - th2np1k - th1n*xi1 + 
     &            th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)/2.))/2. + 
     &     (9*L2*m2*(1 + xi1)*(v2n*(-1 + xi1) - v2np1k*(1 + xi1))**2*
     &        Cos(th1n + th1np1k - th2n - th2np1k - th1n*xi1 + 
     &          th1np1k*xi1 + th2n*xi1 - th2np1k*xi1))/4. + 
     &     9*g*(m1 + 2*m2)*(1 + xi1)*
     &      Cos(th1n + (th2n*(-1 + xi1))/2. - th1n*xi1 + 
     &        th1np1k*(1 + xi1) - (th2np1k*(1 + xi1))/2.))/
     &   (L2*(8*m1 + 15*m2 - 9*m2*
     &        Cos(th1n + th1np1k - th2n - th2np1k - th1n*xi1 + 
     &          th1np1k*xi1 + th2n*xi1 - th2np1k*xi1))) - 
     &  (9*m2*(1 - xi1)*Sin(th1n + th1np1k - th2n - th2np1k + th1n*xi1 -
     &       th1np1k*xi1 - th2n*xi1 + th2np1k*xi1)*
     &     (9*g*(m1 + 2*m2)*Sin(th1n + (th2n*(-1 - xi1))/2. + 
     &          th1np1k*(1 - xi1) - (th2np1k*(1 - xi1))/2. + th1n*xi1) -
     &       3*g*(m1 + 6*m2)*Sin((th2n + th2np1k + th2n*xi1 - 
     &            th2np1k*xi1)/2.) + 
     &       3*L1*(m1 + 3*m2)*(v1n + v1np1k + v1n*xi1 - v1np1k*xi1)**2*
     &        Sin((th1n + th1np1k - th2n - th2np1k + th1n*xi1 - 
     &            th1np1k*xi1 - th2n*xi1 + th2np1k*xi1)/2.) + 
     &       (9*L2*m2*(v2n*(-1 - xi1) - v2np1k*(1 - xi1))**2*
     &          Sin(th1n + th1np1k - th2n - th2np1k + th1n*xi1 - 
     &            th1np1k*xi1 - th2n*xi1 + th2np1k*xi1))/4.))/
     &   (L2*(8*m1 + 15*m2 - 9*m2*
     &         Cos(th1n + th1np1k - th2n - th2np1k + th1n*xi1 - 
     &           th1np1k*xi1 - th2n*xi1 + th2np1k*xi1))**2) - 
     &  (9*m2*(1 + xi1)*Sin(th1n + th1np1k - th2n - th2np1k - th1n*xi1 +
     &       th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)*
     &     (3*L1*(m1 + 3*m2)*(v1n + v1np1k - v1n*xi1 + v1np1k*xi1)**2*
     &        Sin((th1n + th1np1k - th2n - th2np1k - th1n*xi1 + 
     &            th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)/2.) + 
     &       (9*L2*m2*(v2n*(-1 + xi1) - v2np1k*(1 + xi1))**2*
     &          Sin(th1n + th1np1k - th2n - th2np1k - th1n*xi1 + 
     &            th1np1k*xi1 + th2n*xi1 - th2np1k*xi1))/4. - 
     &       3*g*(m1 + 6*m2)*Sin((th2n + th2np1k - th2n*xi1 + 
     &            th2np1k*xi1)/2.) + 
     &       9*g*(m1 + 2*m2)*Sin(th1n + (th2n*(-1 + xi1))/2. - th1n*xi1+
     &          th1np1k*(1 + xi1) - (th2np1k*(1 + xi1))/2.)))/
     &   (L2*(8*m1 + 15*m2 - 9*m2*
     &         Cos(th1n + th1np1k - th2n - th2np1k - th1n*xi1 + 
     &           th1np1k*xi1 + th2n*xi1 - th2np1k*xi1))**2)))/18. - 
     & (L1*L2*m2*v1np1k*Sin(th1np1k - th2np1k)*
     & ((L1*(v1n + v1np1k)*(-2*g*(m1 + 2*m2)*Cos(th1n) + 
     &       2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &       L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &          v1np1k*v2np1k*(Cos(th1n - th2np1k) - 
     &             Cos(th1np1k - th2np1k)))))/(8.*(th1n - th1np1k)) + 
     &  (L2*m2*(v2n + v2np1k)*(L1*v1np1k*v2np1k*
     &        (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &       L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k))+
     &       2*g*(-Cos(th2n) + Cos(th2np1k))))/(8.*(th2n - th2np1k)) + 
     &  (L1*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &       3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &          v2n*Cos(th1np1k - th2np1k)))*
     &     ((-4*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k)/2.) + 
     &            m2*(3*g*Sin((th1n + th1np1k - 2*th2n - 2*th2np1k)/
     &                  2.) + 
     &               L2*(v2n + v2np1k)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &               (3*L1*(-v1n - v1np1k)**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k))/4.)))/
     &        (3.*L1*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((-3*(g*(4*m1 + 5*m2)*
     &                  Sin((th1n + th1np1k - th1n*xi1 + th1np1k*xi1)/
     &                   2.) + 
     &                 m2*(3*g*
     &                   Sin((th1n + th1np1k - 2*th2n - 2*th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + 2*th2n*xi1 - 
     &                  2*th2np1k*xi1)/2.) + 
     &                   L2*(v2n + v2np1k - v2n*xi1 + v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                  th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 + xi1) - v1np1k*(1 + xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4.)))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))) - (3*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k + th1n*xi1 - 
     &                   th1np1k*xi1)/2.) + 
     &                 m2*(L2*(v2n + v2np1k + v2n*xi1 - v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 - xi1) - v1np1k*(1 - xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4. + 
     &                   3*g*Sin((th1n + th1np1k - 2*th2n - 
     &                   2*th2np1k + th1n*xi1 - th1np1k*xi1 - 
     &                   2*th2n*xi1 + 2*th2np1k*xi1)/2.))))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                 th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & )))))/18.))/12. + (L2*m2*(2*L2*(v2n + v2np1k) + 
     &       3*L1*v1n*Cos(th1n - th2n) + 
     &       3*L1*v1np1k*Cos(th1np1k - th2np1k))*
     &     ((4*(3*L1*(m1 + 3*m2)*(v1n + v1np1k)**2*
     &             Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &            (9*L2*m2*(-v2n - v2np1k)**2*
     &               Sin(th1n + th1np1k - th2n - th2np1k))/4. + 
     &            9*g*(m1 + 2*m2)*
     &             Sin(th1n + th1np1k - th2n/2. - th2np1k/2.) - 
     &            3*g*(m1 + 6*m2)*Sin((th2n + th2np1k)/2.)))/
     &        (9.*L2*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((9*g*(m1 + 2*m2)*
     &                Sin(th1n + (th2n*(-1 - xi1))/2. + 
     &                  th1np1k*(1 - xi1) - (th2np1k*(1 - xi1))/2. + 
     &                  th1n*xi1) - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k + th2n*xi1 - th2np1k*xi1)/2.)
     & + 3*L1*(m1 + 3*m2)*(v1n + v1np1k + v1n*xi1 - v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 - xi1) - v2np1k*(1 - xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4.)/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                  th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))) + (3*L1*(m1 + 3*m2)*(v1n + v1np1k - v1n*xi1 + v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 + xi1) - v2np1k*(1 + xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4. - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k - th2n*xi1 + th2np1k*xi1)/2.)
     & + 9*g*(m1 + 2*m2)*Sin(th1n + (th2n*(-1 + xi1))/2. - th1n*xi1 + 
     &                  th1np1k*(1 + xi1) - (th2np1k*(1 + xi1))/2.))/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)
     & ))))/18.))/12.))/
     &    (4.*((L1**2*(-2*g*(m1 + 2*m2)*Cos(th1n) + 
     &        2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &        L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &           v1np1k*v2np1k*
     &            (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k))))**2)/
     &   (16.*(th1n - th1np1k)**2) + 
     &  (L2**2*m2**2*(2*L2*(v2n + v2np1k) + 3*L1*v1n*Cos(th1n - th2n) + 
     &        3*L1*v1np1k*Cos(th1np1k - th2np1k))**2)/144. + 
     &  (L1**2*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &        3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &           v2n*Cos(th1np1k - th2np1k)))**2)/144. + 
     &  (L2**2*m2**2*(L1*v1np1k*v2np1k*
     &         (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &        L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k))+
     &  2*g*(-Cos(th2n) + Cos(th2np1k)))**2)/(16.*(th2n - th2np1k)**2))
     & ) - (L2*m2*(2*L2*(v2n + v2np1k) + 3*L1*v1n*Cos(th1n - th2n) + 
     &  3*L1*v1np1k*Cos(th1np1k - th2np1k))*
     & ((L1**2*(-2*g*(m1 + 2*m2)*Cos(th1n) + 
     &        2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &        L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &           v1np1k*v2np1k*
     &            (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k))))**2)/
     &   (8.*(th1n - th1np1k)**3) - 
     &  (L1*L2**2*m2**2*v1np1k*
     &     (2*L2*(v2n + v2np1k) + 3*L1*v1n*Cos(th1n - th2n) + 
     &    3*L1*v1np1k*Cos(th1np1k - th2np1k))*Sin(th1np1k - th2np1k))/
     &   24. - (L1**2*L2*m2*v2n*
     &     (2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &       3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &          v2n*Cos(th1np1k - th2np1k)))*Sin(th1np1k-th2np1k))/24.+
     &  (L1*L2**2*m2**2*v1n*v2n*
     &     (L1*v1np1k*v2np1k*(Cos(th1n - th2n) - Cos(th1n - th2np1k)) +
     &       L1*v1n*v2n*(Cos(th1np1k - th2n) - 
     &          Cos(th1np1k - th2np1k)) + 
     &       2*g*(-Cos(th2n) + Cos(th2np1k)))*
     &     (-Sin(th1np1k - th2n) + Sin(th1np1k - th2np1k)))/
     &   (8.*(th2n - th2np1k)**2) + 
     &  (L1**2*(-2*g*(m1 + 2*m2)*Cos(th1n) + 
     &       2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &       L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &          v1np1k*v2np1k*
     &           (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k))))*
     &     (-2*g*(m1 + 2*m2)*Sin(th1np1k) + 
     &       L2*m2*(v1n*v2n*Sin(th1np1k - th2n) + 
     &          v1np1k*v2np1k*Sin(th1np1k - th2np1k))))/
     &   (8.*(th1n - th1np1k)**2))*
     & ((L1*(v1n + v1np1k)*(-2*g*(m1 + 2*m2)*Cos(th1n) + 
     &       2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &       L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &          v1np1k*v2np1k*(Cos(th1n - th2np1k) - 
     &             Cos(th1np1k - th2np1k)))))/(8.*(th1n - th1np1k)) + 
     &  (L2*m2*(v2n + v2np1k)*(L1*v1np1k*v2np1k*
     &        (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &       L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k))+
     &       2*g*(-Cos(th2n) + Cos(th2np1k))))/(8.*(th2n - th2np1k)) + 
     &  (L1*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &       3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &          v2n*Cos(th1np1k - th2np1k)))*
     &     ((-4*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k)/2.) + 
     &            m2*(3*g*Sin((th1n + th1np1k - 2*th2n - 2*th2np1k)/
     &                  2.) + 
     &               L2*(v2n + v2np1k)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &               (3*L1*(-v1n - v1np1k)**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k))/4.)))/
     &        (3.*L1*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((-3*(g*(4*m1 + 5*m2)*
     &                  Sin((th1n + th1np1k - th1n*xi1 + th1np1k*xi1)/
     &                   2.) + 
     &                 m2*(3*g*
     &                   Sin((th1n + th1np1k - 2*th2n - 2*th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + 2*th2n*xi1 - 
     &                  2*th2np1k*xi1)/2.) + 
     &                   L2*(v2n + v2np1k - v2n*xi1 + v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                  th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 + xi1) - v1np1k*(1 + xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4.)))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))) - (3*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k + th1n*xi1 - 
     &                   th1np1k*xi1)/2.) + 
     &                 m2*(L2*(v2n + v2np1k + v2n*xi1 - v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 - xi1) - v1np1k*(1 - xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4. + 
     &                   3*g*Sin((th1n + th1np1k - 2*th2n - 
     &                   2*th2np1k + th1n*xi1 - th1np1k*xi1 - 
     &                   2*th2n*xi1 + 2*th2np1k*xi1)/2.))))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & )))))/18.))/12. + (L2*m2*(2*L2*(v2n + v2np1k) + 
     &       3*L1*v1n*Cos(th1n - th2n) + 
     &       3*L1*v1np1k*Cos(th1np1k - th2np1k))*
     &     ((4*(3*L1*(m1 + 3*m2)*(v1n + v1np1k)**2*
     &             Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &            (9*L2*m2*(-v2n - v2np1k)**2*
     &               Sin(th1n + th1np1k - th2n - th2np1k))/4. + 
     &            9*g*(m1 + 2*m2)*
     &             Sin(th1n + th1np1k - th2n/2. - th2np1k/2.) - 
     &            3*g*(m1 + 6*m2)*Sin((th2n + th2np1k)/2.)))/
     &        (9.*L2*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((9*g*(m1 + 2*m2)*
     &                Sin(th1n + (th2n*(-1 - xi1))/2. + 
     &                  th1np1k*(1 - xi1) - (th2np1k*(1 - xi1))/2. + 
     &                  th1n*xi1) - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k + th2n*xi1 - th2np1k*xi1)/2.)
     & + 3*L1*(m1 + 3*m2)*(v1n + v1np1k + v1n*xi1 - v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 - xi1) - v2np1k*(1 - xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4.)/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                  th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))) + (3*L1*(m1 + 3*m2)*(v1n + v1np1k - v1n*xi1 + v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 + xi1) - v2np1k*(1 + xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4. - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k - th2n*xi1 + th2np1k*xi1)/2.)
     & + 9*g*(m1 + 2*m2)*Sin(th1n + (th2n*(-1 + xi1))/2. - th1n*xi1 + 
     &                  th1np1k*(1 + xi1) - (th2np1k*(1 + xi1))/2.))/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                 th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)
     & ))))/18.))/12.))/
     &    (12.*((L1**2*(-2*g*(m1 + 2*m2)*Cos(th1n) + 
     &         2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &         L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) +
     &            v1np1k*v2np1k*
     &             (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k))))**2)/
     &    (16.*(th1n - th1np1k)**2) + 
     &   (L2**2*m2**2*(2*L2*(v2n + v2np1k) + 3*L1*v1n*Cos(th1n - th2n) +
     &         3*L1*v1np1k*Cos(th1np1k - th2np1k))**2)/144. + 
     &   (L1**2*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &         3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &            v2n*Cos(th1np1k - th2np1k)))**2)/144. + 
     &   (L2**2*m2**2*(L1*v1np1k*v2np1k*
     &          (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &         L1*v1n*v2n*(Cos(th1np1k - th2n) - 
     &            Cos(th1np1k - th2np1k)) + 
     &         2*g*(-Cos(th2n) + Cos(th2np1k)))**2)/
     &    (16.*(th2n - th2np1k)**2))**2) + 
     &   (L2*m2*(2*L2*(v2n + v2np1k) + 3*L1*v1n*Cos(th1n - th2n) + 
     &  3*L1*v1np1k*Cos(th1np1k - th2np1k))*
     & ((L1*(v1n + v1np1k)*(-2*g*(m1 + 2*m2)*Cos(th1n) + 
     &       2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &       L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &          v1np1k*v2np1k*(Cos(th1n - th2np1k) - 
     &             Cos(th1np1k - th2np1k)))))/(8.*(th1n - th1np1k)**2) +
     &  (L1*L2*m2*v1n*v2n*(v2n + v2np1k)*
     &     (-Sin(th1np1k - th2n) + Sin(th1np1k - th2np1k)))/
     &   (8.*(th2n - th2np1k)) + 
     &  (L1*(v1n + v1np1k)*(-2*g*(m1 + 2*m2)*Sin(th1np1k) + 
     &       L2*m2*(v1n*v2n*Sin(th1np1k - th2n) + 
     &          v1np1k*v2np1k*Sin(th1np1k - th2np1k))))/
     &   (8.*(th1n - th1np1k)) - 
     &  (L1*L2*m2*v2n*Sin(th1np1k - th2np1k)*
     &     ((-4*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k)/2.) + 
     &            m2*(3*g*Sin((th1n + th1np1k - 2*th2n - 2*th2np1k)/
     &                  2.) + 
     &               L2*(v2n + v2np1k)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &               (3*L1*(-v1n - v1np1k)**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k))/4.)))/
     &        (3.*L1*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((-3*(g*(4*m1 + 5*m2)*
     &                  Sin((th1n + th1np1k - th1n*xi1 + th1np1k*xi1)/
     &                   2.) + 
     &                 m2*(3*g*
     &                   Sin((th1n + th1np1k - 2*th2n - 2*th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + 2*th2n*xi1 - 
     &                   2*th2np1k*xi1)/2.) + 
     &                   L2*(v2n + v2np1k - v2n*xi1 + v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 + xi1) - v1np1k*(1 + xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4.)))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))) - (3*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k + th1n*xi1 - 
     &                   th1np1k*xi1)/2.) + 
     &                 m2*(L2*(v2n + v2np1k + v2n*xi1 - v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &                    (3*L1*(v1n*(-1 - xi1) - v1np1k*(1 - xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4. + 
     &                    3*g*
     &                    Sin((th1n + th1np1k - 2*th2n - 2*th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - 2*th2n*xi1 + 
     &                   2*th2np1k*xi1)/2.))))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                  th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1)
     & ))))/18.))/4. + (L1*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &       3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &          v2n*Cos(th1np1k - th2np1k)))*
     &     ((-4*((g*(4*m1 + 5*m2)*Cos((th1n + th1np1k)/2.))/2. + 
     &            m2*((3*g*Cos((th1n + th1np1k - 2*th2n - 2*th2np1k)/
     &                   2.))/2. + 
     &               (L2*(v2n + v2np1k)**2*
     &                  Cos((th1n + th1np1k - th2n - th2np1k)/2.))/2. +
     &               (3*L1*(-v1n - v1np1k)**2*
     &                  Cos(th1n + th1np1k - th2n - th2np1k))/4.)))/
     &        (3.*L1*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (12*m2*Sin(th1n + th1np1k - th2n - th2np1k)*
     &          (g*(4*m1 + 5*m2)*Sin((th1n + th1np1k)/2.) + 
     &            m2*(3*g*Sin((th1n + th1np1k - 2*th2n - 2*th2np1k)/
     &                  2.) + 
     &               L2*(v2n + v2np1k)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &               (3*L1*(-v1n - v1np1k)**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k))/4.)))/
     &        (L1*(8*m1 + 15*m2 - 
     &             9*m2*Cos(th1n + th1np1k - th2n - th2np1k))**2) + 
     &       (5*((-3*((g*(4*m1 + 5*m2)*(1 + xi1)*
     &                   Cos((th1n + th1np1k - th1n*xi1 + 
     &                  th1np1k*xi1)/2.))/2. + 
     &                 m2*((3*g*(1 + xi1)*
     &                   Cos((th1n + th1np1k - 2*th2n - 2*th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + 2*th2n*xi1 - 
     &                  2*th2np1k*xi1)/2.))/2. + 
     &                   (L2*(1 + xi1)*
     &                   (v2n + v2np1k - v2n*xi1 + v2np1k*xi1)**2*
     &                   Cos((th1n + th1np1k - th2n - th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                  th2np1k*xi1)/2.))/2. + 
     &                   (3*L1*(1 + xi1)*
     &                   (v1n*(-1 + xi1) - v1np1k*(1 + xi1))**2*
     &                   Cos(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4.)))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))) - (3*((g*(4*m1 + 5*m2)*(1 - xi1)*
     &                   Cos((th1n + th1np1k + th1n*xi1 - 
     &                  th1np1k*xi1)/2.))/2. + 
     &                 m2*((L2*(1 - xi1)*
     &                   (v2n + v2np1k + v2n*xi1 - v2np1k*xi1)**2*
     &                   Cos((th1n + th1np1k - th2n - th2np1k + 
     &                  th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                  th2np1k*xi1)/2.))/2. + 
     &                   (3*L1*(v1n*(-1 - xi1) - v1np1k*(1 - xi1))**2*
     &                   (1 - xi1)*
     &                   Cos(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4. + 
     &                   (3*g*(1 - xi1)*
     &                   Cos((th1n + th1np1k - 2*th2n - 2*th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - 2*th2n*xi1 + 
     &                   2*th2np1k*xi1)/2.))/2.)))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                  th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))) + (27*m2*(1 + xi1)*Sin(th1n + th1np1k - th2n - th2np1k - 
     &                 th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)*
     &               (g*(4*m1 + 5*m2)*
     &                  Sin((th1n + th1np1k - th1n*xi1 + th1np1k*xi1)/
     &                   2.) + 
     &                 m2*(3*g*
     &                   Sin((th1n + th1np1k - 2*th2n - 2*th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + 2*th2n*xi1 - 
     &                   2*th2np1k*xi1)/2.) + 
     &                   L2*(v2n + v2np1k - v2n*xi1 + v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 + xi1) - v1np1k*(1 + xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4.)))/
     &             (L1*(8*m1 + 15*m2 - 
     &                  9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))**2) + (27*m2*(1 - xi1)*Sin(th1n + th1np1k - th2n - th2np1k + 
     &                 th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1)*
     &               (g*(4*m1 + 5*m2)*
     &                  Sin((th1n + th1np1k + th1n*xi1 - th1np1k*xi1)/
     &                   2.) + 
     &                 m2*(L2*(v2n + v2np1k + v2n*xi1 - v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &                    (3*L1*(v1n*(-1 - xi1) - v1np1k*(1 - xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4. + 
     &                    3*g*
     &                    Sin((th1n + th1np1k - 2*th2n - 2*th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - 2*th2n*xi1 + 
     &                   2*th2np1k*xi1)/2.))))/
     &             (L1*(8*m1 + 15*m2 - 
     &                  9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                  th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))**2)))/18.))/12. - (L1*L2*m2*v1np1k*Sin(th1np1k - th2np1k)*
     &     ((4*(3*L1*(m1 + 3*m2)*(v1n + v1np1k)**2*
     &             Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &            (9*L2*m2*(-v2n - v2np1k)**2*
     &               Sin(th1n + th1np1k - th2n - th2np1k))/4. + 
     &            9*g*(m1 + 2*m2)*
     &             Sin(th1n + th1np1k - th2n/2. - th2np1k/2.) - 
     &            3*g*(m1 + 6*m2)*Sin((th2n + th2np1k)/2.)))/
     &        (9.*L2*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((9*g*(m1 + 2*m2)*
     &                Sin(th1n + (th2n*(-1 - xi1))/2. + 
     &                  th1np1k*(1 - xi1) - (th2np1k*(1 - xi1))/2. +
     &                  th1n*xi1) - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k + th2n*xi1 - th2np1k*xi1)/2.)
     & + 3*L1*(m1 + 3*m2)*(v1n + v1np1k + v1n*xi1 - v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 - xi1) - v2np1k*(1 - xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4.)/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                  th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))) + (3*L1*(m1 + 3*m2)*(v1n + v1np1k - v1n*xi1 + v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 + xi1) - v2np1k*(1 + xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4. - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k - th2n*xi1 + th2np1k*xi1)/2.)
     & + 9*g*(m1 + 2*m2)*Sin(th1n + (th2n*(-1 + xi1))/2. - th1n*xi1 + 
     &                  th1np1k*(1 + xi1) - (th2np1k*(1 + xi1))/2.))/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                 th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)
     & ))))/18.))/4.+(L2*m2*(2*L2*(v2n+ v2np1k)+3*L1*v1n*Cos(th1n-th2n)+
     &       3*L1*v1np1k*Cos(th1np1k - th2np1k))*
     &     ((4*((3*L1*(m1 + 3*m2)*(v1n + v1np1k)**2*
     &               Cos((th1n + th1np1k - th2n - th2np1k)/2.))/2. + 
     &            (9*L2*m2*(-v2n - v2np1k)**2*
     &               Cos(th1n + th1np1k - th2n - th2np1k))/4. + 
     &            9*g*(m1 + 2*m2)*
     &             Cos(th1n + th1np1k - th2n/2. - th2np1k/2.)))/
     &        (9.*L2*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) - 
     &       (4*m2*Sin(th1n + th1np1k - th2n - th2np1k)*
     &          (3*L1*(m1 + 3*m2)*(v1n + v1np1k)**2*
     &             Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &            (9*L2*m2*(-v2n - v2np1k)**2*
     &               Sin(th1n + th1np1k - th2n - th2np1k))/4. + 
     &            9*g*(m1 + 2*m2)*
     &             Sin(th1n + th1np1k - th2n/2. - th2np1k/2.) - 
     &            3*g*(m1 + 6*m2)*Sin((th2n + th2np1k)/2.)))/
     &        (L2*(8*m1 + 15*m2 - 
     &             9*m2*Cos(th1n + th1np1k - th2n - th2np1k))**2) + 
     &       (5*((9*g*(m1 + 2*m2)*(1 - xi1)*
     &                Cos(th1n + (th2n*(-1 - xi1))/2. + 
     &                  th1np1k*(1 - xi1) - (th2np1k*(1 - xi1))/2. + 
     &                  th1n*xi1) + 
     &               (3*L1*(m1 + 3*m2)*(1 - xi1)*
     &                  (v1n + v1np1k + v1n*xi1 - v1np1k*xi1)**2*
     &                  Cos((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.))/2. + 
     &               (9*L2*m2*(v2n*(-1 - xi1) - v2np1k*(1 - xi1))**2*
     &                  (1 - xi1)*
     &                  Cos(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))/4.)/(L2*(8*m1 + 15*m2 - 9*m2*
     &                  Cos(th1n + th1np1k - th2n - th2np1k + th1n*xi1 - 
     &                    th1np1k*xi1 - th2n*xi1 + th2np1k*xi1))) + 
     &            ((3*L1*(m1 + 3*m2)*(1 + xi1)*
     &                  (v1n + v1np1k - v1n*xi1 + v1np1k*xi1)**2*
     &                  Cos((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.))/2. + 
     &               (9*L2*m2*(1 + xi1)*
     &                  (v2n*(-1 + xi1) - v2np1k*(1 + xi1))**2*
     &                  Cos(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4. + 
     &               9*g*(m1 + 2*m2)*(1 + xi1)*
     &                Cos(th1n + (th2n*(-1 + xi1))/2. - th1n*xi1 + 
     &                  th1np1k*(1 + xi1) - (th2np1k*(1 + xi1))/2.))/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                 th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)
     & ))-(9*m2*(1-xi1)*Sin(th1n + th1np1k - th2n - th2np1k + th1n*xi1-
     &                 th1np1k*xi1 - th2n*xi1 + th2np1k*xi1)*
     &               (9*g*(m1 + 2*m2)*
     &                  Sin(th1n + (th2n*(-1 - xi1))/2. + 
     &                   th1np1k*(1 - xi1) - (th2np1k*(1 - xi1))/2. + 
     &                   th1n*xi1) - 
     &                 3*g*(m1 + 6*m2)*
     &                  Sin((th2n + th2np1k + th2n*xi1 - th2np1k*xi1)/
     &                   2.) + 
     &                 3*L1*(m1 + 3*m2)*
     &                  (v1n + v1np1k + v1n*xi1 - v1np1k*xi1)**2*
     &                  Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &                 (9*L2*m2*(v2n*(-1 - xi1) - v2np1k*(1 - xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))/4.))/(L2*(8*m1 + 15*m2 - 9*m2*
     &                   Cos(th1n + th1np1k - th2n - th2np1k + 
     &                th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))**2) - (9*m2*(1 + xi1)*Sin(th1n + th1np1k - th2n - th2np1k - 
     &                 th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)*
     &               (3*L1*(m1 + 3*m2)*
     &                  (v1n + v1np1k - v1n*xi1 + v1np1k*xi1)**2*
     &                  Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & )/2.) + (9*L2*m2*(v2n*(-1 + xi1) - v2np1k*(1 + xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))/4.-3*g*(m1+6*m2)*Sin((th2n+th2np1k - th2n*xi1 + th2np1k*xi1)/
     &                   2.) + 
     &                 9*g*(m1 + 2*m2)*
     &                  Sin(th1n + (th2n*(-1 + xi1))/2. - th1n*xi1 + 
     &                    th1np1k*(1 + xi1) - (th2np1k*(1 + xi1))/2.)))/
     &             (L2*(8*m1 + 15*m2 - 
     &                  9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)
     & )**2)))/18.))/12.))/
     &    (12.*((L1**2*(-2*g*(m1 + 2*m2)*Cos(th1n) + 
     &        2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &        L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) +
     &           v1np1k*v2np1k*
     &            (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k))))**2)/
     &   (16.*(th1n - th1np1k)**2) + 
     &  (L2**2*m2**2*(2*L2*(v2n + v2np1k) + 3*L1*v1n*Cos(th1n - th2n) +
     &        3*L1*v1np1k*Cos(th1np1k - th2np1k))**2)/144. + 
     &  (L1**2*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &        3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &           v2n*Cos(th1np1k - th2np1k)))**2)/144. + 
     &  (L2**2*m2**2*(L1*v1np1k*v2np1k*
     &         (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &   L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k)) + 
     & 2*g*(-Cos(th2n) + Cos(th2np1k)))**2)/(16.*(th2n - th2np1k)**2)))

      m42=(-4*((-3*L1*(m1 + 3*m2)*(v1n + v1np1k)**2*
     &     Cos((th1n + th1np1k - th2n - th2np1k)/2.))/2. - 
     &(9*L2*m2*(-v2n - v2np1k)**2*Cos(th1n + th1np1k - th2n - th2np1k))/
     &4. - (9*g*(m1 + 2*m2)*Cos(th1n + th1np1k - th2n/2. - th2np1k/2.))/
     &   2. - (3*g*(m1 + 6*m2)*Cos((th2n + th2np1k)/2.))/2.))/
     & (9.*L2*(8*m1+15*m2 - 9*m2*Cos(th1n + th1np1k - th2n - th2np1k)))-
     & (4*m2*Sin(th1n + th1np1k - th2n - th2np1k)*
     & (3*L1*(m1 + 3*m2)*(v1n + v1np1k)**2*
     &   Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &(9*L2*m2*(-v2n - v2np1k)**2*Sin(th1n + th1np1k - th2n - th2np1k))/
     & 4. + 9*g*(m1 + 2*m2)*Sin(th1n + th1np1k - th2n/2. - th2np1k/2.)-
     &  3*g*(m1 + 6*m2)*Sin((th2n + th2np1k)/2.)))/
     &(L2*(8*m1 + 15*m2 - 9*m2*Cos(th1n + th1np1k-th2n - th2np1k))**2)-
     & (5*(((9*g*(m1 + 2*m2)*(-1 + xi1)*
     &        Cos(th1n + (th2n*(-1 - xi1))/2. + th1np1k*(1 - xi1) - 
     &          (th2np1k*(1 - xi1))/2. + th1n*xi1))/2. - 
     &     (3*g*(m1 + 6*m2)*(1 - xi1)*
     &        Cos((th2n + th2np1k + th2n*xi1 - th2np1k*xi1)/2.))/2. + 
     &     (3*L1*(m1 + 3*m2)*(-1 + xi1)*
     &        (v1n + v1np1k + v1n*xi1 - v1np1k*xi1)**2*
     &        Cos((th1n + th1np1k - th2n - th2np1k + th1n*xi1 - 
     &            th1np1k*xi1 - th2n*xi1 + th2np1k*xi1)/2.))/2. + 
     &     (9*L2*m2*(v2n*(-1 - xi1) - v2np1k*(1 - xi1))**2*(-1 + xi1)*
     &        Cos(th1n + th1np1k - th2n - th2np1k + th1n*xi1 - 
     &          th1np1k*xi1 - th2n*xi1 + th2np1k*xi1))/4.)/
     &   (L2*(8*m1 + 15*m2 - 9*m2*
     &        Cos(th1n + th1np1k - th2n - th2np1k + th1n*xi1 - 
     &          th1np1k*xi1 - th2n*xi1 + th2np1k*xi1))) + 
     &  ((3*L1*(m1 + 3*m2)*(-1 - xi1)*
     &        (v1n + v1np1k - v1n*xi1 + v1np1k*xi1)**2*
     &        Cos((th1n + th1np1k - th2n - th2np1k - th1n*xi1 + 
     &            th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)/2.))/2. + 
     &     (9*L2*m2*(-1 - xi1)*(v2n*(-1 + xi1) - v2np1k*(1 + xi1))**2*
     &        Cos(th1n + th1np1k - th2n - th2np1k - th1n*xi1 + 
     &          th1np1k*xi1 + th2n*xi1 - th2np1k*xi1))/4. - 
     &     (3*g*(m1 + 6*m2)*(1 + xi1)*
     &        Cos((th2n + th2np1k - th2n*xi1 + th2np1k*xi1)/2.))/2. + 
     &     (9*g*(m1 + 2*m2)*(-1 - xi1)*
     &        Cos(th1n + (th2n*(-1 + xi1))/2. - th1n*xi1 + 
     &          th1np1k*(1 + xi1) - (th2np1k*(1 + xi1))/2.))/2.)/
     &   (L2*(8*m1 + 15*m2 - 9*m2*
     &        Cos(th1n + th1np1k - th2n - th2np1k - th1n*xi1 + 
     &          th1np1k*xi1 + th2n*xi1 - th2np1k*xi1))) - 
     &  (9*m2*(-1 + xi1)*Sin(th1n + th1np1k - th2n - th2np1k + 
     &       th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1)*
     &     (9*g*(m1 + 2*m2)*Sin(th1n + (th2n*(-1 - xi1))/2. + 
     &          th1np1k*(1 - xi1) - (th2np1k*(1 - xi1))/2. + th1n*xi1) - 
     &       3*g*(m1 + 6*m2)*Sin((th2n + th2np1k + th2n*xi1 - 
     &            th2np1k*xi1)/2.) + 
     &       3*L1*(m1 + 3*m2)*(v1n + v1np1k + v1n*xi1 - v1np1k*xi1)**2*
     &        Sin((th1n + th1np1k - th2n - th2np1k + th1n*xi1 - 
     &            th1np1k*xi1 - th2n*xi1 + th2np1k*xi1)/2.) + 
     &       (9*L2*m2*(v2n*(-1 - xi1) - v2np1k*(1 - xi1))**2*
     &          Sin(th1n + th1np1k - th2n - th2np1k + th1n*xi1 - 
     &            th1np1k*xi1 - th2n*xi1 + th2np1k*xi1))/4.))/
     &   (L2*(8*m1 + 15*m2 - 9*m2*
     &         Cos(th1n + th1np1k - th2n - th2np1k + th1n*xi1 - 
     &           th1np1k*xi1 - th2n*xi1 + th2np1k*xi1))**2) - 
     &  (9*m2*(-1 - xi1)*Sin(th1n + th1np1k - th2n - th2np1k - th1n*xi1+
     &       th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)*
     &     (3*L1*(m1 + 3*m2)*(v1n + v1np1k - v1n*xi1 + v1np1k*xi1)**2*
     &        Sin((th1n + th1np1k - th2n - th2np1k - th1n*xi1 + 
     &            th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)/2.) + 
     &       (9*L2*m2*(v2n*(-1 + xi1) - v2np1k*(1 + xi1))**2*
     &          Sin(th1n + th1np1k - th2n - th2np1k - th1n*xi1 + 
     &            th1np1k*xi1 + th2n*xi1 - th2np1k*xi1))/4. - 
     &       3*g*(m1 + 6*m2)*Sin((th2n + th2np1k - th2n*xi1 + 
     &            th2np1k*xi1)/2.) + 
     &       9*g*(m1 + 2*m2)*Sin(th1n + (th2n*(-1 + xi1))/2. - th1n*xi1+
     &          th1np1k*(1 + xi1) - (th2np1k*(1 + xi1))/2.)))/
     &   (L2*(8*m1 + 15*m2 - 9*m2*
     &         Cos(th1n + th1np1k - th2n - th2np1k - th1n*xi1 + 
     &           th1np1k*xi1 + th2n*xi1 - th2np1k*xi1))**2)))/18. + 
     & (L1*L2*m2*v1np1k*Sin(th1np1k - th2np1k)*
     & ((L1*(v1n + v1np1k)*(-2*g*(m1 + 2*m2)*Cos(th1n) + 
     &       2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &       L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &          v1np1k*v2np1k*(Cos(th1n - th2np1k) - 
     &             Cos(th1np1k - th2np1k)))))/(8.*(th1n - th1np1k)) + 
     &  (L2*m2*(v2n + v2np1k)*(L1*v1np1k*v2np1k*
     &        (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &       L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k)) +
     &       2*g*(-Cos(th2n) + Cos(th2np1k))))/(8.*(th2n - th2np1k)) + 
     &  (L1*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &       3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &          v2n*Cos(th1np1k - th2np1k)))*
     &     ((-4*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k)/2.) + 
     &            m2*(3*g*Sin((th1n + th1np1k - 2*th2n - 2*th2np1k)/
     &                  2.) + 
     &               L2*(v2n + v2np1k)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &               (3*L1*(-v1n - v1np1k)**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k))/4.)))/
     &        (3.*L1*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((-3*(g*(4*m1 + 5*m2)*
     &                  Sin((th1n + th1np1k - th1n*xi1 + th1np1k*xi1)/
     &                   2.) + 
     &                 m2*(3*g*
     &                   Sin((th1n + th1np1k - 2*th2n - 2*th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + 2*th2n*xi1 - 
     &                  2*th2np1k*xi1)/2.) + 
     &                   L2*(v2n + v2np1k - v2n*xi1 + v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                  th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 + xi1) - v1np1k*(1 + xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4.)))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))) - (3*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k + th1n*xi1 - 
     &                   th1np1k*xi1)/2.) + 
     &                 m2*(L2*(v2n + v2np1k + v2n*xi1 - v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 - xi1) - v1np1k*(1 - xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4. + 
     &                   3*g*Sin((th1n + th1np1k - 2*th2n - 
     &                   2*th2np1k + th1n*xi1 - th1np1k*xi1 - 
     &                   2*th2n*xi1 + 2*th2np1k*xi1)/2.))))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                  th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & )))))/18.))/12. + (L2*m2*(2*L2*(v2n + v2np1k) + 
     &       3*L1*v1n*Cos(th1n - th2n) + 
     &       3*L1*v1np1k*Cos(th1np1k - th2np1k))*
     &     ((4*(3*L1*(m1 + 3*m2)*(v1n + v1np1k)**2*
     &             Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &            (9*L2*m2*(-v2n - v2np1k)**2*
     &               Sin(th1n + th1np1k - th2n - th2np1k))/4. + 
     &            9*g*(m1 + 2*m2)*
     &             Sin(th1n + th1np1k - th2n/2. - th2np1k/2.) - 
     &            3*g*(m1 + 6*m2)*Sin((th2n + th2np1k)/2.)))/
     &        (9.*L2*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((9*g*(m1 + 2*m2)*
     &                Sin(th1n + (th2n*(-1 - xi1))/2. + 
     &                  th1np1k*(1 - xi1) - (th2np1k*(1 - xi1))/2. + 
     &                  th1n*xi1) - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k + th2n*xi1 - th2np1k*xi1)/2.)
     & + 3*L1*(m1 + 3*m2)*(v1n + v1np1k + v1n*xi1 - v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 - xi1) - v2np1k*(1 - xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4.)/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                  th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))) + (3*L1*(m1 + 3*m2)*(v1n + v1np1k - v1n*xi1 + v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 + xi1) - v2np1k*(1 + xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4. - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k - th2n*xi1 + th2np1k*xi1)/2.)
     & + 9*g*(m1 + 2*m2)*Sin(th1n + (th2n*(-1 + xi1))/2. - th1n*xi1 + 
     &                  th1np1k*(1 + xi1) - (th2np1k*(1 + xi1))/2.))/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)
     & ))))/18.))/12.))/
     &    (4.*((L1**2*(-2*g*(m1 + 2*m2)*Cos(th1n) + 
     &        2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &        L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) +
     &           v1np1k*v2np1k*
     &            (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k))))**2)/
     &   (16.*(th1n - th1np1k)**2) + 
     &  (L2**2*m2**2*(2*L2*(v2n + v2np1k) + 3*L1*v1n*Cos(th1n - th2n) +
     &        3*L1*v1np1k*Cos(th1np1k - th2np1k))**2)/144. + 
     &  (L1**2*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &        3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &           v2n*Cos(th1np1k - th2np1k)))**2)/144. + 
     &  (L2**2*m2**2*(L1*v1np1k*v2np1k*
     &         (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &        L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k))+
     &  2*g*(-Cos(th2n) + Cos(th2np1k)))**2)/(16.*(th2n - th2np1k)**2))
     & ) - (L2*m2*(2*L2*(v2n + v2np1k) + 3*L1*v1n*Cos(th1n - th2n) + 
     &  3*L1*v1np1k*Cos(th1np1k - th2np1k))*
     & ((L2**2*m2**2*(L1*v1np1k*v2np1k*
     &         (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &        L1*v1n*v2n*(Cos(th1np1k - th2n) - 
     &           Cos(th1np1k - th2np1k)) + 
     &        2*g*(-Cos(th2n) + Cos(th2np1k)))**2)/
     &   (8.*(th2n - th2np1k)**3) + 
     &  (L1**2*L2*m2*v1np1k*v2np1k*
     &     (-2*g*(m1 + 2*m2)*Cos(th1n) + 2*g*(m1 + 2*m2)*Cos(th1np1k) +
     &       L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &          v1np1k*v2np1k*
     &           (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k))))*
     &     (Sin(th1n - th2np1k) - Sin(th1np1k - th2np1k)))/
     &   (8.*(th1n - th1np1k)**2) + 
     &  (L1*L2**2*m2**2*v1np1k*
     &     (2*L2*(v2n + v2np1k) + 3*L1*v1n*Cos(th1n - th2n) + 
     &     3*L1*v1np1k*Cos(th1np1k - th2np1k))*Sin(th1np1k - th2np1k))/
     &   24. + (L1**2*L2*m2*v2n*
     &     (2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &       3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &          v2n*Cos(th1np1k - th2np1k)))*Sin(th1np1k-th2np1k))/24.+
     &  (L2**2*m2**2*(L1*v1np1k*v2np1k*
     &        (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &       L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k))+
     &       2*g*(-Cos(th2n) + Cos(th2np1k)))*
     &     (-(L1*v1np1k*v2np1k*Sin(th1n - th2np1k)) - 
     &       L1*v1n*v2n*Sin(th1np1k - th2np1k) - 2*g*Sin(th2np1k)))/
     &   (8.*(th2n - th2np1k)**2))*
     & ((L1*(v1n + v1np1k)*(-2*g*(m1 + 2*m2)*Cos(th1n) + 
     &       2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &       L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &          v1np1k*v2np1k*(Cos(th1n - th2np1k) - 
     &             Cos(th1np1k - th2np1k)))))/(8.*(th1n - th1np1k)) + 
     &  (L2*m2*(v2n + v2np1k)*(L1*v1np1k*v2np1k*
     &        (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &       L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k)) +
     &       2*g*(-Cos(th2n) + Cos(th2np1k))))/(8.*(th2n - th2np1k)) + 
     &  (L1*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &       3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &          v2n*Cos(th1np1k - th2np1k)))*
     &     ((-4*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k)/2.) + 
     &            m2*(3*g*Sin((th1n + th1np1k - 2*th2n - 2*th2np1k)/
     &                  2.) + 
     &               L2*(v2n + v2np1k)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &               (3*L1*(-v1n - v1np1k)**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k))/4.)))/
     &        (3.*L1*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((-3*(g*(4*m1 + 5*m2)*
     &                  Sin((th1n + th1np1k - th1n*xi1 + th1np1k*xi1)/
     &                   2.) + 
     &                 m2*(3*g*
     &                   Sin((th1n + th1np1k - 2*th2n - 2*th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + 2*th2n*xi1 - 
     &                  2*th2np1k*xi1)/2.) + 
     &                   L2*(v2n + v2np1k - v2n*xi1 + v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                  th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 + xi1) - v1np1k*(1 + xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4.)))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))) - (3*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k + th1n*xi1 - 
     &                   th1np1k*xi1)/2.) + 
     &                 m2*(L2*(v2n + v2np1k + v2n*xi1 - v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 - xi1) - v1np1k*(1 - xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4. + 
     &                   3*g*Sin((th1n + th1np1k - 2*th2n - 
     &                   2*th2np1k + th1n*xi1 - th1np1k*xi1 - 
     &                   2*th2n*xi1 + 2*th2np1k*xi1)/2.))))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                 th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & )))))/18.))/12. + (L2*m2*(2*L2*(v2n + v2np1k) + 
     &       3*L1*v1n*Cos(th1n - th2n) + 
     &       3*L1*v1np1k*Cos(th1np1k - th2np1k))*
     &     ((4*(3*L1*(m1 + 3*m2)*(v1n + v1np1k)**2*
     &             Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &            (9*L2*m2*(-v2n - v2np1k)**2*
     &               Sin(th1n + th1np1k - th2n - th2np1k))/4. + 
     &            9*g*(m1 + 2*m2)*
     &             Sin(th1n + th1np1k - th2n/2. - th2np1k/2.) - 
     &            3*g*(m1 + 6*m2)*Sin((th2n + th2np1k)/2.)))/
     &        (9.*L2*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((9*g*(m1 + 2*m2)*
     &                Sin(th1n + (th2n*(-1 - xi1))/2. + 
     &                  th1np1k*(1 - xi1) - (th2np1k*(1 - xi1))/2. + 
     &                  th1n*xi1) - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k + th2n*xi1 - th2np1k*xi1)/2.)
     & + 3*L1*(m1 + 3*m2)*(v1n + v1np1k + v1n*xi1 - v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 - xi1) - v2np1k*(1 - xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4.)/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                  th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))) + (3*L1*(m1 + 3*m2)*(v1n + v1np1k - v1n*xi1 + v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 + xi1) - v2np1k*(1 + xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4. - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k - th2n*xi1 + th2np1k*xi1)/2.)
     & + 9*g*(m1 + 2*m2)*Sin(th1n + (th2n*(-1 + xi1))/2. - th1n*xi1 + 
     &                  th1np1k*(1 + xi1) - (th2np1k*(1 + xi1))/2.))/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                 th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)
     & ))))/18.))/12.))/
     &    (12.*((L1**2*(-2*g*(m1 + 2*m2)*Cos(th1n) + 
     &         2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &         L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) +
     &            v1np1k*v2np1k*
     &             (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k))))**2)/
     &    (16.*(th1n - th1np1k)**2) + 
     &   (L2**2*m2**2*(2*L2*(v2n + v2np1k) + 3*L1*v1n*Cos(th1n - th2n) + 
     &         3*L1*v1np1k*Cos(th1np1k - th2np1k))**2)/144. + 
     &   (L1**2*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &         3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &            v2n*Cos(th1np1k - th2np1k)))**2)/144. + 
     &   (L2**2*m2**2*(L1*v1np1k*v2np1k*
     &          (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &         L1*v1n*v2n*(Cos(th1np1k - th2n) - 
     &            Cos(th1np1k - th2np1k)) + 
     &         2*g*(-Cos(th2n) + Cos(th2np1k)))**2)/
     &    (16.*(th2n - th2np1k)**2))**2) + 
     &   (L2*m2*(2*L2*(v2n + v2np1k) + 3*L1*v1n*Cos(th1n - th2n) + 
     &  3*L1*v1np1k*Cos(th1np1k - th2np1k))*
     &      ((L2*m2*(v2n + v2np1k)*(L1*v1np1k*v2np1k*
     &        (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &       L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k)) + 
     &       2*g*(-Cos(th2n) + Cos(th2np1k))))/(8.*(th2n - th2np1k)**2)+
     &  (L1*L2*m2*v1np1k*(v1n + v1np1k)*v2np1k*
     &     (Sin(th1n - th2np1k) - Sin(th1np1k - th2np1k)))/
     &   (8.*(th1n - th1np1k)) + 
     &  (L2*m2*(v2n + v2np1k)*(-(L1*v1np1k*v2np1k*Sin(th1n - th2np1k)) - 
     &       L1*v1n*v2n*Sin(th1np1k - th2np1k) - 2*g*Sin(th2np1k)))/
     &   (8.*(th2n - th2np1k)) + 
     &  (L1*L2*m2*v2n*Sin(th1np1k - th2np1k)*
     &     ((-4*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k)/2.) + 
     &            m2*(3*g*Sin((th1n + th1np1k - 2*th2n - 2*th2np1k)/
     &                  2.) + 
     &               L2*(v2n + v2np1k)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &               (3*L1*(-v1n - v1np1k)**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k))/4.)))/
     &        (3.*L1*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((-3*(g*(4*m1 + 5*m2)*
     &                  Sin((th1n + th1np1k - th1n*xi1 + th1np1k*xi1)/
     &                   2.) + 
     &                 m2*(3*g*
     &                   Sin((th1n + th1np1k - 2*th2n - 2*th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + 2*th2n*xi1 - 
     &                   2*th2np1k*xi1)/2.) + 
     &                   L2*(v2n + v2np1k - v2n*xi1 + v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 + xi1) - v1np1k*(1 + xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4.)))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))) - (3*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k + th1n*xi1 - 
     &                   th1np1k*xi1)/2.) + 
     &                 m2*(L2*(v2n + v2np1k + v2n*xi1 - v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &                    (3*L1*(v1n*(-1 - xi1) - v1np1k*(1 - xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4. + 
     &                    3*g*
     &                    Sin((th1n + th1np1k - 2*th2n - 2*th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - 2*th2n*xi1 + 
     &                   2*th2np1k*xi1)/2.))))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                 th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1)
     & ))))/18.))/4. + (L1*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &       3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &          v2n*Cos(th1np1k - th2np1k)))*
     &     ((-4*m2*(-3*g*Cos((th1n + th1np1k - 2*th2n - 2*th2np1k)/2.) -
     &            (L2*(v2n + v2np1k)**2*
     &               Cos((th1n + th1np1k - th2n - th2np1k)/2.))/2. - 
     &            (3*L1*(-v1n - v1np1k)**2*
     &               Cos(th1n + th1np1k - th2n - th2np1k))/4.))/
     &        (3.*L1*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) - 
     &       (12*m2*Sin(th1n + th1np1k - th2n - th2np1k)*
     &          (g*(4*m1 + 5*m2)*Sin((th1n + th1np1k)/2.) + 
     &            m2*(3*g*Sin((th1n + th1np1k - 2*th2n - 2*th2np1k)/
     &                  2.) + 
     &               L2*(v2n + v2np1k)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &               (3*L1*(-v1n - v1np1k)**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k))/4.)))/
     &        (L1*(8*m1 + 15*m2 - 
     &             9*m2*Cos(th1n + th1np1k - th2n - th2np1k))**2) + 
     &       (5*((-3*m2*((3*g*(-2 - 2*xi1)*
     &                   Cos((th1n + th1np1k - 2*th2n - 2*th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + 2*th2n*xi1 - 
     &                  2*th2np1k*xi1)/2.))/2. + 
     &                 (L2*(-1 - xi1)*
     &                   (v2n + v2np1k - v2n*xi1 + v2np1k*xi1)**2*
     &                   Cos((th1n + th1np1k - th2n - th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                  th2np1k*xi1)/2.))/2. + 
     &                 (3*L1*(-1 - xi1)*
     &                   (v1n*(-1 + xi1) - v1np1k*(1 + xi1))**2*
     &                   Cos(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4.))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                 th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & )))-(3*m2*((L2*(-1+xi1)*(v2n+v2np1k + v2n*xi1 - v2np1k*xi1)**2*
     &                   Cos((th1n + th1np1k - th2n - th2np1k + 
     &                  th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                  th2np1k*xi1)/2.))/2. + 
     &                 (3*L1*(v1n*(-1 - xi1) - v1np1k*(1 - xi1))**2*
     &                   (-1 + xi1)*
     &                   Cos(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4. + 
     &                 (3*g*(-2 + 2*xi1)*
     &                   Cos((th1n + th1np1k - 2*th2n - 2*th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - 2*th2n*xi1 + 
     &                   2*th2np1k*xi1)/2.))/2.))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                 th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))) + (27*m2*(-1 - xi1)*Sin(th1n + th1np1k - th2n - th2np1k - 
     &                 th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)*
     &               (g*(4*m1 + 5*m2)*
     &                  Sin((th1n + th1np1k - th1n*xi1 + th1np1k*xi1)/
     &                   2.) + 
     &                 m2*(3*g*
     &                   Sin((th1n + th1np1k - 2*th2n - 2*th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + 2*th2n*xi1 - 
     &                   2*th2np1k*xi1)/2.) + 
     &                   L2*(v2n + v2np1k - v2n*xi1 + v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 + xi1) - v1np1k*(1 + xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4.)))/
     &             (L1*(8*m1 + 15*m2 - 
     &                  9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))**2) + (27*m2*(-1 + xi1)*Sin(th1n + th1np1k - th2n - th2np1k + 
     &                 th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1)*
     &               (g*(4*m1 + 5*m2)*
     &                  Sin((th1n + th1np1k + th1n*xi1 - th1np1k*xi1)/
     &                   2.) + 
     &                 m2*(L2*(v2n + v2np1k + v2n*xi1 - v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &                    (3*L1*(v1n*(-1 - xi1) - v1np1k*(1 - xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4. + 
     &                    3*g*
     &                    Sin((th1n + th1np1k - 2*th2n - 2*th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - 2*th2n*xi1 + 
     &                   2*th2np1k*xi1)/2.))))/
     &             (L1*(8*m1 + 15*m2 - 
     &                  9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                  th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))**2)))/18.))/12. + (L1*L2*m2*v1np1k*Sin(th1np1k - th2np1k)*
     &     ((4*(3*L1*(m1 + 3*m2)*(v1n + v1np1k)**2*
     &             Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &            (9*L2*m2*(-v2n - v2np1k)**2*
     &               Sin(th1n + th1np1k - th2n - th2np1k))/4. + 
     &            9*g*(m1 + 2*m2)*
     &             Sin(th1n + th1np1k - th2n/2. - th2np1k/2.) - 
     &            3*g*(m1 + 6*m2)*Sin((th2n + th2np1k)/2.)))/
     &        (9.*L2*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((9*g*(m1 + 2*m2)*
     &                Sin(th1n + (th2n*(-1 - xi1))/2. + 
     &                  th1np1k*(1 - xi1) - (th2np1k*(1 - xi1))/2. + 
     &                  th1n*xi1) - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k + th2n*xi1 - th2np1k*xi1)/2.)
     & + 3*L1*(m1 + 3*m2)*(v1n + v1np1k + v1n*xi1 - v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 - xi1) - v2np1k*(1 - xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4.)/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                 th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))) + (3*L1*(m1 + 3*m2)*(v1n + v1np1k - v1n*xi1 + v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 + xi1) - v2np1k*(1 + xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4. - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k - th2n*xi1 + th2np1k*xi1)/2.)
     & + 9*g*(m1 + 2*m2)*Sin(th1n + (th2n*(-1 + xi1))/2. - th1n*xi1 + 
     &                  th1np1k*(1 + xi1) - (th2np1k*(1 + xi1))/2.))/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)
     & ))))/18.))/4.+ (L2*m2*(2*L2*(v2n+v2np1k)+3*L1*v1n*Cos(th1n-th2n)+
     &       3*L1*v1np1k*Cos(th1np1k - th2np1k))*
     &     ((4*((-3*L1*(m1 + 3*m2)*(v1n + v1np1k)**2*
     &               Cos((th1n + th1np1k - th2n - th2np1k)/2.))/2. - 
     &            (9*L2*m2*(-v2n - v2np1k)**2*
     &               Cos(th1n + th1np1k - th2n - th2np1k))/4. - 
     &            (9*g*(m1 + 2*m2)*
     &               Cos(th1n + th1np1k - th2n/2. - th2np1k/2.))/2. - 
     &            (3*g*(m1 + 6*m2)*Cos((th2n + th2np1k)/2.))/2.))/
     &        (9.*L2*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (4*m2*Sin(th1n + th1np1k - th2n - th2np1k)*
     &          (3*L1*(m1 + 3*m2)*(v1n + v1np1k)**2*
     &             Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &            (9*L2*m2*(-v2n - v2np1k)**2*
     &               Sin(th1n + th1np1k - th2n - th2np1k))/4. + 
     &            9*g*(m1 + 2*m2)*
     &             Sin(th1n + th1np1k - th2n/2. - th2np1k/2.) - 
     &            3*g*(m1 + 6*m2)*Sin((th2n + th2np1k)/2.)))/
     &        (L2*(8*m1 + 15*m2 - 
     &             9*m2*Cos(th1n + th1np1k - th2n - th2np1k))**2) + 
     &       (5*(((9*g*(m1 + 2*m2)*(-1 + xi1)*
     &                  Cos(th1n + (th2n*(-1 - xi1))/2. + 
     &                   th1np1k*(1 - xi1) - (th2np1k*(1 - xi1))/2. +
     &                   th1n*xi1))/2. - 
     &               (3*g*(m1 + 6*m2)*(1 - xi1)*
     &                  Cos((th2n + th2np1k + th2n*xi1 - th2np1k*xi1)/
     &                   2.))/2. + 
     &               (3*L1*(m1 + 3*m2)*(-1 + xi1)*
     &                  (v1n + v1np1k + v1n*xi1 - v1np1k*xi1)**2*
     &                  Cos((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.))/2. + 
     &               (9*L2*m2*(v2n*(-1 - xi1) - v2np1k*(1 - xi1))**2*
     &                  (-1 + xi1)*
     &                  Cos(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))/4.)/(L2*(8*m1 + 15*m2 - 9*m2*
     &                  Cos(th1n + th1np1k - th2n - th2np1k + th1n*xi1 -
     &                    th1np1k*xi1 - th2n*xi1 + th2np1k*xi1))) + 
     &            ((3*L1*(m1 + 3*m2)*(-1 - xi1)*
     &                  (v1n + v1np1k - v1n*xi1 + v1np1k*xi1)**2*
     &                  Cos((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.))/2. + 
     &               (9*L2*m2*(-1 - xi1)*
     &                  (v2n*(-1 + xi1) - v2np1k*(1 + xi1))**2*
     &                  Cos(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4. - 
     &               (3*g*(m1 + 6*m2)*(1 + xi1)*
     &                  Cos((th2n + th2np1k - th2n*xi1 + th2np1k*xi1)/
     &                   2.))/2. + 
     &               (9*g*(m1 + 2*m2)*(-1 - xi1)*
     &                  Cos(th1n + (th2n*(-1 + xi1))/2. - th1n*xi1 + 
     &               th1np1k*(1 + xi1) - (th2np1k*(1 + xi1))/2.))/2.)/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                 th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)
     & ))-(9*m2*(-1+xi1)*Sin(th1n+th1np1k - th2n - th2np1k + th1n*xi1-
     &                 th1np1k*xi1 - th2n*xi1 + th2np1k*xi1)*
     &               (9*g*(m1 + 2*m2)*
     &                  Sin(th1n + (th2n*(-1 - xi1))/2. + 
     &                   th1np1k*(1 - xi1) - (th2np1k*(1 - xi1))/2. + 
     &                   th1n*xi1) - 
     &                 3*g*(m1 + 6*m2)*
     &                  Sin((th2n + th2np1k + th2n*xi1 - th2np1k*xi1)/
     &                   2.) + 
     &                 3*L1*(m1 + 3*m2)*
     &                  (v1n + v1np1k + v1n*xi1 - v1np1k*xi1)**2*
     &                  Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &                 (9*L2*m2*(v2n*(-1 - xi1) - v2np1k*(1 - xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))/4.))/(L2*(8*m1 + 15*m2 - 9*m2*
     &                   Cos(th1n + th1np1k - th2n - th2np1k + 
     &                  th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))**2) - (9*m2*(-1 - xi1)*Sin(th1n + th1np1k - th2n - th2np1k - 
     &                 th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)*
     &               (3*L1*(m1 + 3*m2)*
     &                  (v1n + v1np1k - v1n*xi1 + v1np1k*xi1)**2*
     &                  Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & )/2.) + (9*L2*m2*(v2n*(-1 + xi1) - v2np1k*(1 + xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))/4.-3*g*(m1+6*m2)*Sin((th2n+th2np1k - th2n*xi1 + th2np1k*xi1)/
     &                   2.) + 
     &                 9*g*(m1 + 2*m2)*
     &                  Sin(th1n + (th2n*(-1 + xi1))/2. - th1n*xi1 + 
     &                    th1np1k*(1 + xi1) - (th2np1k*(1 + xi1))/2.)))/
     &             (L2*(8*m1 + 15*m2 - 
     &                  9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)
     & )**2)))/18.))/12.))/
     &    (12.*((L1**2*(-2*g*(m1 + 2*m2)*Cos(th1n) + 
     &        2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &        L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &           v1np1k*v2np1k*
     &            (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k))))**2)/
     &   (16.*(th1n - th1np1k)**2) + 
     &  (L2**2*m2**2*(2*L2*(v2n + v2np1k) + 3*L1*v1n*Cos(th1n - th2n) +
     &        3*L1*v1np1k*Cos(th1np1k - th2np1k))**2)/144. + 
     &  (L1**2*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &        3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &           v2n*Cos(th1np1k - th2np1k)))**2)/144. + 
     &  (L2**2*m2**2*(L1*v1np1k*v2np1k*
     &         (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &        L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k))+
     & 2*g*(-Cos(th2n) + Cos(th2np1k)))**2)/(16.*(th2n - th2np1k)**2)))

      m43= (-8*L1*(m1 + 3*m2)*(v1n + v1np1k)*
     &       Sin((th1n + th1np1k - th2n - th2np1k)/2.))/
     & (3.*L2*(8*m1+15*m2-9*m2*Cos(th1n+th1np1k - th2n - th2np1k))) -
     & (5*((6*L1*(m1+3*m2)*(1+xi1)*(v1n+v1np1k - v1n*xi1 + v1np1k*xi1)*
     &     Sin((th1n + th1np1k - th2n - th2np1k - th1n*xi1 + 
     &         th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)/2.))/
     &   (L2*(8*m1 + 15*m2 - 9*m2*
     &        Cos(th1n + th1np1k - th2n - th2np1k - th1n*xi1 + 
     &          th1np1k*xi1 + th2n*xi1 - th2np1k*xi1))) + 
     &(6*L1*(m1 + 3*m2)*(1 - xi1)*(v1n + v1np1k + v1n*xi1 - v1np1k*xi1)*
     &     Sin((th1n + th1np1k - th2n - th2np1k + th1n*xi1 - 
     &         th1np1k*xi1 - th2n*xi1 + th2np1k*xi1)/2.))/
     &   (L2*(8*m1 + 15*m2 - 9*m2*
     &        Cos(th1n + th1np1k - th2n - th2np1k + th1n*xi1 - 
     &          th1np1k*xi1 - th2n*xi1 + th2np1k*xi1)))))/18. + 
     & (L2*m2*(2*L2*(v2n + v2np1k) + 3*L1*v1n*Cos(th1n - th2n) + 
     &  3*L1*v1np1k*Cos(th1np1k - th2np1k))*
     & ((L1*L2*m2*v2np1k*(v2n + v2np1k)*
     &   (Cos(th1n - th2n) - Cos(th1n - th2np1k)))/(8.*(th2n - th2np1k))
     & + (L1*(-2*g*(m1 + 2*m2)*Cos(th1n) + 2*g*(m1 + 2*m2)*Cos(th1np1k)+
     &       L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &          v1np1k*v2np1k*(Cos(th1n - th2np1k) - 
     &             Cos(th1np1k - th2np1k)))))/(8.*(th1n - th1np1k)) + 
     &  (L1*L2*m2*(v1n + v1np1k)*v2np1k*
     &     (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k)))/
     &   (8.*(th1n - th1np1k)) + 
     &  (L2*m2*(2*L2*(v2n + v2np1k) + 3*L1*v1n*Cos(th1n - th2n) + 
     &       3*L1*v1np1k*Cos(th1np1k - th2np1k))*
     &     ((8*L1*(m1 + 3*m2)*(v1n + v1np1k)*
     &          Sin((th1n + th1np1k - th2n - th2np1k)/2.))/
     &        (3.*L2*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((6*L1*(m1 + 3*m2)*(1 + xi1)*
     &               (v1n + v1np1k - v1n*xi1 + v1np1k*xi1)*
     &               Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.))/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))) + (6*L1*(m1 + 3*m2)*(1 - xi1)*
     &               (v1n + v1np1k + v1n*xi1 - v1np1k*xi1)*
     &               Sin((th1n + th1np1k - th2n - th2np1k + th1n*xi1 - 
     &                   th1np1k*xi1 - th2n*xi1 + th2np1k*xi1)/2.))/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                  th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & )))))/18.))/12. + (L1*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &       3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &          v2n*Cos(th1np1k - th2np1k)))*
     &     ((2*m2*(-v1n - v1np1k)*Sin(th1n + th1np1k - th2n - th2np1k))/
     &        (8*m1 + 15*m2 - 9*m2*Cos(th1n + th1np1k - th2n - th2np1k))
     & + (5*((-9*m2*(-1 - xi1)*(v1n*(-1 + xi1) - v1np1k*(1 + xi1))*
     &               Sin(th1n + th1np1k - th2n - th2np1k - th1n*xi1 + 
     &                 th1np1k*xi1 + th2n*xi1 - th2np1k*xi1))/
     &             (2.*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))) - (9*m2*(v1n*(-1 - xi1) - v1np1k*(1 - xi1))*(-1 + xi1)*
     &               Sin(th1n + th1np1k - th2n - th2np1k + th1n*xi1 - 
     &                 th1np1k*xi1 - th2n*xi1 + th2np1k*xi1))/
     &             (2.*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                  th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & )))))/18.))/12. + (L1**2*(m1 + 3*m2)*
     &     ((-4*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k)/2.) + 
     &            m2*(3*g*Sin((th1n + th1np1k - 2*th2n - 2*th2np1k)/
     &                  2.) + 
     &               L2*(v2n + v2np1k)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &               (3*L1*(-v1n - v1np1k)**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k))/4.)))/
     &        (3.*L1*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((-3*(g*(4*m1 + 5*m2)*
     &                  Sin((th1n + th1np1k - th1n*xi1 + th1np1k*xi1)/
     &                   2.) + 
     &                 m2*(3*g*
     &                   Sin((th1n + th1np1k - 2*th2n - 2*th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + 2*th2n*xi1 - 
     &                  2*th2np1k*xi1)/2.) + 
     &                   L2*(v2n + v2np1k - v2n*xi1 + v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                  th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 + xi1) - v1np1k*(1 + xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4.)))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))) - (3*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k + th1n*xi1 - 
     &                   th1np1k*xi1)/2.) + 
     &                 m2*(L2*(v2n + v2np1k + v2n*xi1 - v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 - xi1) - v1np1k*(1 - xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4. + 
     &                   3*g*Sin((th1n + th1np1k - 2*th2n - 
     &                   2*th2np1k + th1n*xi1 - th1np1k*xi1 - 
     &                   2*th2n*xi1 + 2*th2np1k*xi1)/2.))))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                  th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & )))))/18.))/6. + (L1*L2*m2*Cos(th1np1k - th2np1k)*
     &     ((4*(3*L1*(m1 + 3*m2)*(v1n + v1np1k)**2*
     &             Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &            (9*L2*m2*(-v2n - v2np1k)**2*
     &               Sin(th1n + th1np1k - th2n - th2np1k))/4. + 
     &            9*g*(m1 + 2*m2)*
     &             Sin(th1n + th1np1k - th2n/2. - th2np1k/2.) - 
     &            3*g*(m1 + 6*m2)*Sin((th2n + th2np1k)/2.)))/
     &        (9.*L2*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((9*g*(m1 + 2*m2)*
     &                Sin(th1n + (th2n*(-1 - xi1))/2. + 
     &                  th1np1k*(1 - xi1) - (th2np1k*(1 - xi1))/2. + 
     &                  th1n*xi1) - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k + th2n*xi1 - th2np1k*xi1)/2.)
     & + 3*L1*(m1 + 3*m2)*(v1n + v1np1k + v1n*xi1 - v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 - xi1) - v2np1k*(1 - xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4.)/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                  th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))) + (3*L1*(m1 + 3*m2)*(v1n + v1np1k - v1n*xi1 + v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 + xi1) - v2np1k*(1 + xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4. - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k - th2n*xi1 + th2np1k*xi1)/2.)
     & + 9*g*(m1 + 2*m2)*Sin(th1n + (th2n*(-1 + xi1))/2. - th1n*xi1 + 
     &                  th1np1k*(1 + xi1) - (th2np1k*(1 + xi1))/2.))/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)
     & ))))/18.))/4.))/
     & (12.*((L1**2*(-2*g*(m1 + 2*m2)*Cos(th1n) + 
     &        2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &        L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &           v1np1k*v2np1k*
     &            (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k))))**2)/
     &   (16.*(th1n - th1np1k)**2) + 
     &  (L2**2*m2**2*(2*L2*(v2n + v2np1k) + 3*L1*v1n*Cos(th1n - th2n) + 
     &        3*L1*v1np1k*Cos(th1np1k - th2np1k))**2)/144. + 
     &  (L1**2*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &        3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &           v2n*Cos(th1np1k - th2np1k)))**2)/144. + 
     &  (L2**2*m2**2*(L1*v1np1k*v2np1k*
     &         (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &        L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k))+
     &  2*g*(-Cos(th2n) + Cos(th2np1k)))**2)/(16.*(th2n - th2np1k)**2))
     & ) - (L2*m2*(2*L2*(v2n + v2np1k) + 3*L1*v1n*Cos(th1n - th2n) + 
     &  3*L1*v1np1k*Cos(th1np1k - th2np1k))*
     & ((L1**2*L2*m2*v2np1k*(-2*g*(m1 + 2*m2)*Cos(th1n) + 
     &       2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &       L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &          v1np1k*v2np1k*
     &           (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k))))*
     &     (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k)))/
     &   (8.*(th1n - th1np1k)**2) + 
     &  (L1*L2**2*m2**2*Cos(th1np1k - th2np1k)*
     &     (2*L2*(v2n + v2np1k) + 3*L1*v1n*Cos(th1n - th2n) + 
     &       3*L1*v1np1k*Cos(th1np1k - th2np1k)))/24. + 
     &  (L1**3*(m1 + 3*m2)*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &       3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &          v2n*Cos(th1np1k - th2np1k))))/36. + 
     &  (L1*L2**2*m2**2*v2np1k*(Cos(th1n - th2n) - Cos(th1n - th2np1k))*
     &     (L1*v1np1k*v2np1k*(Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &       L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k)) + 
     &      2*g*(-Cos(th2n) + Cos(th2np1k))))/(8.*(th2n - th2np1k)**2))*
     &  ((L1*(v1n + v1np1k)*(-2*g*(m1 + 2*m2)*Cos(th1n) + 
     &       2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &       L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &          v1np1k*v2np1k*(Cos(th1n - th2np1k) - 
     &             Cos(th1np1k - th2np1k)))))/(8.*(th1n - th1np1k)) + 
     &  (L2*m2*(v2n + v2np1k)*(L1*v1np1k*v2np1k*
     &        (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &       L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k)) +
     &       2*g*(-Cos(th2n) + Cos(th2np1k))))/(8.*(th2n - th2np1k)) + 
     &  (L1*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &       3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &          v2n*Cos(th1np1k - th2np1k)))*
     &     ((-4*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k)/2.) + 
     &            m2*(3*g*Sin((th1n + th1np1k - 2*th2n - 2*th2np1k)/
     &                  2.) + 
     &               L2*(v2n + v2np1k)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &               (3*L1*(-v1n - v1np1k)**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k))/4.)))/
     &        (3.*L1*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((-3*(g*(4*m1 + 5*m2)*
     &                  Sin((th1n + th1np1k - th1n*xi1 + th1np1k*xi1)/
     &                   2.) + 
     &                 m2*(3*g*
     &                   Sin((th1n + th1np1k - 2*th2n - 2*th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + 2*th2n*xi1 - 
     &                  2*th2np1k*xi1)/2.) + 
     &                   L2*(v2n + v2np1k - v2n*xi1 + v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                  th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 + xi1) - v1np1k*(1 + xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4.)))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))) - (3*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k + th1n*xi1 - 
     &                   th1np1k*xi1)/2.) + 
     &                 m2*(L2*(v2n + v2np1k + v2n*xi1 - v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 - xi1) - v1np1k*(1 - xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4. + 
     &                   3*g*Sin((th1n + th1np1k - 2*th2n - 
     &                   2*th2np1k + th1n*xi1 - th1np1k*xi1 - 
     &                   2*th2n*xi1 + 2*th2np1k*xi1)/2.))))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                  th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & )))))/18.))/12. + (L2*m2*(2*L2*(v2n + v2np1k) + 
     &       3*L1*v1n*Cos(th1n - th2n) + 
     &       3*L1*v1np1k*Cos(th1np1k - th2np1k))*
     &     ((4*(3*L1*(m1 + 3*m2)*(v1n + v1np1k)**2*
     &             Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &            (9*L2*m2*(-v2n - v2np1k)**2*
     &               Sin(th1n + th1np1k - th2n - th2np1k))/4. + 
     &            9*g*(m1 + 2*m2)*
     &             Sin(th1n + th1np1k - th2n/2. - th2np1k/2.) - 
     &            3*g*(m1 + 6*m2)*Sin((th2n + th2np1k)/2.)))/
     &        (9.*L2*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((9*g*(m1 + 2*m2)*
     &                Sin(th1n + (th2n*(-1 - xi1))/2. + 
     &                  th1np1k*(1 - xi1) - (th2np1k*(1 - xi1))/2. + 
     &                  th1n*xi1) - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k + th2n*xi1 - th2np1k*xi1)/2.)
     & + 3*L1*(m1 + 3*m2)*(v1n + v1np1k + v1n*xi1 - v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 - xi1) - v2np1k*(1 - xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4.)/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                  th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))) + (3*L1*(m1 + 3*m2)*(v1n + v1np1k - v1n*xi1 + v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 + xi1) - v2np1k*(1 + xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4. - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k - th2n*xi1 + th2np1k*xi1)/2.)
     & + 9*g*(m1 + 2*m2)*Sin(th1n + (th2n*(-1 + xi1))/2. - th1n*xi1 + 
     &                  th1np1k*(1 + xi1) - (th2np1k*(1 + xi1))/2.))/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)
     & ))))/18.))/12.))/
     &    (12.*((L1**2*(-2*g*(m1 + 2*m2)*Cos(th1n) + 
     &         2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &         L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) +
     &            v1np1k*v2np1k*
     &             (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k))))**2)/
     &    (16.*(th1n - th1np1k)**2) + 
     &   (L2**2*m2**2*(2*L2*(v2n + v2np1k) + 3*L1*v1n*Cos(th1n - th2n) +
     &         3*L1*v1np1k*Cos(th1np1k - th2np1k))**2)/144. + 
     &   (L1**2*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &         3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &            v2n*Cos(th1np1k - th2np1k)))**2)/144. + 
     &   (L2**2*m2**2*(L1*v1np1k*v2np1k*
     &          (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &         L1*v1n*v2n*(Cos(th1np1k - th2n) - 
     &            Cos(th1np1k - th2np1k)) + 
     &         2*g*(-Cos(th2n) + Cos(th2np1k)))**2)/
     &    (16.*(th2n - th2np1k)**2))**2) + 
     &   (L1*L2*m2*Cos(th1np1k - th2np1k)*
     & ((L1*(v1n + v1np1k)*(-2*g*(m1 + 2*m2)*Cos(th1n) + 
     &       2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &       L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &          v1np1k*v2np1k*(Cos(th1n - th2np1k) - 
     &             Cos(th1np1k - th2np1k)))))/(8.*(th1n - th1np1k)) + 
     &  (L2*m2*(v2n + v2np1k)*(L1*v1np1k*v2np1k*
     &        (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &       L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k)) +
     &       2*g*(-Cos(th2n) + Cos(th2np1k))))/(8.*(th2n - th2np1k)) + 
     &  (L1*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &       3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &          v2n*Cos(th1np1k - th2np1k)))*
     &     ((-4*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k)/2.) + 
     &            m2*(3*g*Sin((th1n + th1np1k - 2*th2n - 2*th2np1k)/
     &                  2.) + 
     &               L2*(v2n + v2np1k)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &               (3*L1*(-v1n - v1np1k)**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k))/4.)))/
     &        (3.*L1*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((-3*(g*(4*m1 + 5*m2)*
     &                  Sin((th1n + th1np1k - th1n*xi1 + th1np1k*xi1)/
     &                   2.) + 
     &                 m2*(3*g*
     &                   Sin((th1n + th1np1k - 2*th2n - 2*th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + 2*th2n*xi1 - 
     &                   2*th2np1k*xi1)/2.) + 
     &                   L2*(v2n + v2np1k - v2n*xi1 + v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 + xi1) - v1np1k*(1 + xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4.)))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))) - (3*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k + th1n*xi1 - 
     &                   th1np1k*xi1)/2.) + 
     &                 m2*(L2*(v2n + v2np1k + v2n*xi1 - v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &                    (3*L1*(v1n*(-1 - xi1) - v1np1k*(1 - xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4. + 
     &                    3*g*
     &                    Sin((th1n + th1np1k - 2*th2n - 2*th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - 2*th2n*xi1 + 
     &                   2*th2np1k*xi1)/2.))))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                 th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1)
     & ))))/18.))/12. + (L2*m2*(2*L2*(v2n + v2np1k) + 
     &  3*L1*v1n*Cos(th1n - th2n) + 3*L1*v1np1k*Cos(th1np1k - th2np1k)
     & )*((4*(3*L1*(m1 + 3*m2)*(v1n + v1np1k)**2*
     &             Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &            (9*L2*m2*(-v2n - v2np1k)**2*
     &               Sin(th1n + th1np1k - th2n - th2np1k))/4. + 
     &            9*g*(m1 + 2*m2)*
     &             Sin(th1n + th1np1k - th2n/2. - th2np1k/2.) - 
     &            3*g*(m1 + 6*m2)*Sin((th2n + th2np1k)/2.)))/
     &        (9.*L2*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((9*g*(m1 + 2*m2)*
     &                Sin(th1n + (th2n*(-1 - xi1))/2. + 
     &                  th1np1k*(1 - xi1) - (th2np1k*(1 - xi1))/2. + 
     &                  th1n*xi1) - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k + th2n*xi1 - th2np1k*xi1)/2.)
     & + 3*L1*(m1 + 3*m2)*(v1n + v1np1k + v1n*xi1 - v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 - xi1) - v2np1k*(1 - xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))/4.)/(L2*(8*m1 + 15*m2 - 9*m2*
     &                  Cos(th1n + th1np1k - th2n - th2np1k + th1n*xi1 - 
     &                    th1np1k*xi1 - th2n*xi1 + th2np1k*xi1))) + 
     &            (3*L1*(m1 + 3*m2)*
     &                (v1n + v1np1k - v1n*xi1 + v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k - th1n*xi1 + 
     &                   th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 + xi1) - v2np1k*(1 + xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))/4.-3*g*(m1+6*m2)*Sin((th2n+th2np1k - th2n*xi1 + th2np1k*xi1)/
     &                  2.) + 9*g*(m1 + 2*m2)*
     &                Sin(th1n + (th2n*(-1 + xi1))/2. - th1n*xi1 + 
     &                  th1np1k*(1 + xi1) - (th2np1k*(1 + xi1))/2.))/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1))
     & )))/18.))/12.))/
     &    (4.*((L1**2*(-2*g*(m1 + 2*m2)*Cos(th1n) + 
     &        2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &        L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &           v1np1k*v2np1k*
     &            (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k))))**2)/
     &   (16.*(th1n - th1np1k)**2) + 
     &  (L2**2*m2**2*(2*L2*(v2n + v2np1k) + 3*L1*v1n*Cos(th1n - th2n) + 
     &        3*L1*v1np1k*Cos(th1np1k - th2np1k))**2)/144. + 
     &  (L1**2*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &        3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &           v2n*Cos(th1np1k - th2np1k)))**2)/144. + 
     &  (L2**2*m2**2*(L1*v1np1k*v2np1k*
     &         (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &  L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k)) + 
     & 2*g*(-Cos(th2n) + Cos(th2np1k)))**2)/(16.*(th2n - th2np1k)**2)))

       m44=1.0d0/tdelta+(2*m2*(-v2n-v2np1k)*
     & Sin(th1n+th1np1k-th2n-th2np1k))/
     & (8*m1 + 15*m2 - 9*m2*Cos(th1n + th1np1k - th2n - th2np1k)) - 
     & (5*((9*m2*(-1 - xi1)*(v2n*(-1 + xi1) - v2np1k*(1 + xi1))*
     &     Sin(th1n + th1np1k - th2n - th2np1k - th1n*xi1 + th1np1k*xi1+
     &       th2n*xi1 - th2np1k*xi1))/
     &   (2.*(8*m1 + 15*m2 - 9*m2*
     &        Cos(th1n + th1np1k - th2n - th2np1k - th1n*xi1 + 
     &          th1np1k*xi1 + th2n*xi1 - th2np1k*xi1))) + 
     &  (9*m2*(v2n*(-1 - xi1) - v2np1k*(1 - xi1))*(-1 + xi1)*
     &     Sin(th1n + th1np1k - th2n - th2np1k + th1n*xi1 - th1np1k*xi1-
     &       th2n*xi1 + th2np1k*xi1))/
     &   (2.*(8*m1 + 15*m2 - 9*m2*
     &        Cos(th1n + th1np1k - th2n - th2np1k + th1n*xi1 - 
     &          th1np1k*xi1 - th2n*xi1 + th2np1k*xi1)))))/18. + 
     & (L2*m2*(2*L2*(v2n + v2np1k) + 3*L1*v1n*Cos(th1n - th2n) + 
     &  3*L1*v1np1k*Cos(th1np1k - th2np1k))*
     & ((L1*L2*m2*v1np1k*(v2n + v2np1k)*
     &   (Cos(th1n - th2n) - Cos(th1n - th2np1k)))/(8.*(th2n - th2np1k))
     & + (L1*L2*m2*v1np1k*(v1n + v1np1k)*
     &     (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k)))/
     &   (8.*(th1n - th1np1k)) + 
     &  (L2*m2*(L1*v1np1k*v2np1k*
     &        (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &       L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k)) +
     &       2*g*(-Cos(th2n) + Cos(th2np1k))))/(8.*(th2n - th2np1k)) + 
     &  (L1*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &       3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &          v2n*Cos(th1np1k - th2np1k)))*
     &     ((-8*L2*m2*(v2n + v2np1k)*
     &          Sin((th1n + th1np1k - th2n - th2np1k)/2.))/
     &        (3.*L1*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((-6*L2*m2*(1 + xi1)*
     &               (v2n + v2np1k - v2n*xi1 + v2np1k*xi1)*
     &               Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))) - (6*L2*m2*(1 - xi1)*(v2n + v2np1k + v2n*xi1 - v2np1k*xi1)*
     &               Sin((th1n + th1np1k - th2n - th2np1k + th1n*xi1 - 
     &                   th1np1k*xi1 - th2n*xi1 + th2np1k*xi1)/2.))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                 th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & )))))/18.))/12. + (L2*m2*(2*L2*(v2n + v2np1k) + 
     &       3*L1*v1n*Cos(th1n - th2n) + 
     &       3*L1*v1np1k*Cos(th1np1k - th2np1k))*
     &   ((-2*m2*(-v2n - v2np1k)*Sin(th1n + th1np1k - th2n - th2np1k))/
     &        (8*m1 + 15*m2 - 9*m2*Cos(th1n + th1np1k - th2n - th2np1k))
     & + (5*((9*m2*(-1 - xi1)*(v2n*(-1 + xi1) - v2np1k*(1 + xi1))*
     &               Sin(th1n + th1np1k - th2n - th2np1k - th1n*xi1 + 
     &                 th1np1k*xi1 + th2n*xi1 - th2np1k*xi1))/
     &             (2.*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))) + (9*m2*(v2n*(-1 - xi1) - v2np1k*(1 - xi1))*(-1 + xi1)*
     &               Sin(th1n + th1np1k - th2n - th2np1k + th1n*xi1 - 
     &                 th1np1k*xi1 - th2n*xi1 + th2np1k*xi1))/
     &             (2.*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                 th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & )))))/18.))/12. + (L1*L2*m2*Cos(th1n - th2n)*
     &     ((-4*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k)/2.) + 
     &            m2*(3*g*Sin((th1n + th1np1k - 2*th2n - 2*th2np1k)/
     &                  2.) + 
     &               L2*(v2n + v2np1k)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &               (3*L1*(-v1n - v1np1k)**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k))/4.)))/
     &        (3.*L1*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((-3*(g*(4*m1 + 5*m2)*
     &                  Sin((th1n + th1np1k - th1n*xi1 + th1np1k*xi1)/
     &                   2.) + 
     &                 m2*(3*g*
     &                   Sin((th1n + th1np1k - 2*th2n - 2*th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + 2*th2n*xi1 - 
     &                  2*th2np1k*xi1)/2.) + 
     &                   L2*(v2n + v2np1k - v2n*xi1 + v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                  th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 + xi1) - v1np1k*(1 + xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4.)))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))) - (3*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k + th1n*xi1 - 
     &                   th1np1k*xi1)/2.) + 
     &                 m2*(L2*(v2n + v2np1k + v2n*xi1 - v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 - xi1) - v1np1k*(1 - xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4. + 
     &                   3*g*Sin((th1n + th1np1k - 2*th2n - 
     &                   2*th2np1k + th1n*xi1 - th1np1k*xi1 - 
     &                   2*th2n*xi1 + 2*th2np1k*xi1)/2.))))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                  th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & )))))/18.))/4. + (L2**2*m2*((4*
     &          (3*L1*(m1 + 3*m2)*(v1n + v1np1k)**2*
     &             Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &            (9*L2*m2*(-v2n - v2np1k)**2*
     &               Sin(th1n + th1np1k - th2n - th2np1k))/4. + 
     &            9*g*(m1 + 2*m2)*
     &             Sin(th1n + th1np1k - th2n/2. - th2np1k/2.) - 
     &            3*g*(m1 + 6*m2)*Sin((th2n + th2np1k)/2.)))/
     &        (9.*L2*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((9*g*(m1 + 2*m2)*
     &                Sin(th1n + (th2n*(-1 - xi1))/2. + 
     &                  th1np1k*(1 - xi1) - (th2np1k*(1 - xi1))/2. + 
     &                  th1n*xi1) - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k + th2n*xi1 - th2np1k*xi1)/2.)
     & + 3*L1*(m1 + 3*m2)*(v1n + v1np1k + v1n*xi1 - v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 - xi1) - v2np1k*(1 - xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4.)/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                  th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))) + (3*L1*(m1 + 3*m2)*(v1n + v1np1k - v1n*xi1 + v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 + xi1) - v2np1k*(1 + xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4. - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k - th2n*xi1 + th2np1k*xi1)/2.)
     & + 9*g*(m1 + 2*m2)*Sin(th1n + (th2n*(-1 + xi1))/2. - th1n*xi1 + 
     &                  th1np1k*(1 + xi1) - (th2np1k*(1 + xi1))/2.))/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                 th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)
     & ))))/18.))/6.))/
     &    (12.*((L1**2*(-2*g*(m1 + 2*m2)*Cos(th1n) + 
     &        2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &        L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) +
     &           v1np1k*v2np1k*
     &            (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k))))**2)/
     &   (16.*(th1n - th1np1k)**2) + 
     &  (L2**2*m2**2*(2*L2*(v2n + v2np1k) + 3*L1*v1n*Cos(th1n - th2n) +
     &        3*L1*v1np1k*Cos(th1np1k - th2np1k))**2)/144. + 
     &  (L1**2*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &        3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &           v2n*Cos(th1np1k - th2np1k)))**2)/144. + 
     &  (L2**2*m2**2*(L1*v1np1k*v2np1k*
     &         (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &        L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k))+
     &  2*g*(-Cos(th2n) + Cos(th2np1k)))**2)/(16.*(th2n - th2np1k)**2))
     & ) - (L2*m2*(2*L2*(v2n + v2np1k) + 3*L1*v1n*Cos(th1n - th2n) + 
     &  3*L1*v1np1k*Cos(th1np1k - th2np1k))*
     & ((L1**2*L2*m2*v1np1k*(-2*g*(m1 + 2*m2)*Cos(th1n) + 
     &       2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &       L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &          v1np1k*v2np1k*
     &           (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k))))*
     &     (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k)))/
     &   (8.*(th1n - th1np1k)**2) + 
     &  (L2**3*m2**2*(2*L2*(v2n + v2np1k) + 3*L1*v1n*Cos(th1n - th2n) + 
     &       3*L1*v1np1k*Cos(th1np1k - th2np1k)))/36. + 
     &  (L1**2*L2*m2*Cos(th1n - th2n)*
     &     (2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &       3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &          v2n*Cos(th1np1k - th2np1k))))/24. + 
     &  (L1*L2**2*m2**2*v1np1k*(Cos(th1n - th2n) - Cos(th1n - th2np1k))*
     &     (L1*v1np1k*v2np1k*(Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &       L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k)) +
     &     2*g*(-Cos(th2n) + Cos(th2np1k))))/(8.*(th2n - th2np1k)**2))*
     & ((L1*(v1n + v1np1k)*(-2*g*(m1 + 2*m2)*Cos(th1n) + 
     &       2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &       L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &          v1np1k*v2np1k*(Cos(th1n - th2np1k) - 
     &             Cos(th1np1k - th2np1k)))))/(8.*(th1n - th1np1k)) + 
     &  (L2*m2*(v2n + v2np1k)*(L1*v1np1k*v2np1k*
     &        (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &       L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k)) +
     &       2*g*(-Cos(th2n) + Cos(th2np1k))))/(8.*(th2n - th2np1k)) + 
     &  (L1*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &       3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &          v2n*Cos(th1np1k - th2np1k)))*
     &     ((-4*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k)/2.) + 
     &            m2*(3*g*Sin((th1n + th1np1k - 2*th2n - 2*th2np1k)/
     &                  2.) + 
     &               L2*(v2n + v2np1k)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &               (3*L1*(-v1n - v1np1k)**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k))/4.)))/
     &        (3.*L1*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((-3*(g*(4*m1 + 5*m2)*
     &                  Sin((th1n + th1np1k - th1n*xi1 + th1np1k*xi1)/
     &                   2.) + 
     &                 m2*(3*g*
     &                   Sin((th1n + th1np1k - 2*th2n - 2*th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + 2*th2n*xi1 - 
     &                  2*th2np1k*xi1)/2.) + 
     &                   L2*(v2n + v2np1k - v2n*xi1 + v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                  th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 + xi1) - v1np1k*(1 + xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4.)))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))) - (3*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k + th1n*xi1 - 
     &                   th1np1k*xi1)/2.) + 
     &                 m2*(L2*(v2n + v2np1k + v2n*xi1 - v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 - xi1) - v1np1k*(1 - xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4. + 
     &                   3*g*Sin((th1n + th1np1k - 2*th2n - 
     &                   2*th2np1k + th1n*xi1 - th1np1k*xi1 - 
     &                   2*th2n*xi1 + 2*th2np1k*xi1)/2.))))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                 th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & )))))/18.))/12. + (L2*m2*(2*L2*(v2n + v2np1k) + 
     &       3*L1*v1n*Cos(th1n - th2n) + 
     &       3*L1*v1np1k*Cos(th1np1k - th2np1k))*
     &     ((4*(3*L1*(m1 + 3*m2)*(v1n + v1np1k)**2*
     &             Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &            (9*L2*m2*(-v2n - v2np1k)**2*
     &               Sin(th1n + th1np1k - th2n - th2np1k))/4. + 
     &            9*g*(m1 + 2*m2)*
     &             Sin(th1n + th1np1k - th2n/2. - th2np1k/2.) - 
     &            3*g*(m1 + 6*m2)*Sin((th2n + th2np1k)/2.)))/
     &        (9.*L2*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((9*g*(m1 + 2*m2)*
     &                Sin(th1n + (th2n*(-1 - xi1))/2. + 
     &                  th1np1k*(1 - xi1) - (th2np1k*(1 - xi1))/2. + 
     &                  th1n*xi1) - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k + th2n*xi1 - th2np1k*xi1)/2.)
     & + 3*L1*(m1 + 3*m2)*(v1n + v1np1k + v1n*xi1 - v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 - xi1) - v2np1k*(1 - xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4.)/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                 th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))) + (3*L1*(m1 + 3*m2)*(v1n + v1np1k - v1n*xi1 + v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 + xi1) - v2np1k*(1 + xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4. - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k - th2n*xi1 + th2np1k*xi1)/2.)
     & + 9*g*(m1 + 2*m2)*Sin(th1n + (th2n*(-1 + xi1))/2. - th1n*xi1 + 
     &                  th1np1k*(1 + xi1) - (th2np1k*(1 + xi1))/2.))/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)
     & ))))/18.))/12.))/
     &    (12.*((L1**2*(-2*g*(m1 + 2*m2)*Cos(th1n) + 
     &         2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &         L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) +
     &            v1np1k*v2np1k*
     &             (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k))))**2)/
     &    (16.*(th1n - th1np1k)**2) + 
     &   (L2**2*m2**2*(2*L2*(v2n + v2np1k) + 3*L1*v1n*Cos(th1n - th2n) +
     &         3*L1*v1np1k*Cos(th1np1k - th2np1k))**2)/144. + 
     &   (L1**2*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &         3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &            v2n*Cos(th1np1k - th2np1k)))**2)/144. + 
     &   (L2**2*m2**2*(L1*v1np1k*v2np1k*
     &          (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &         L1*v1n*v2n*(Cos(th1np1k - th2n) - 
     &            Cos(th1np1k - th2np1k)) + 
     &         2*g*(-Cos(th2n) + Cos(th2np1k)))**2)/
     &    (16.*(th2n - th2np1k)**2))**2) + 
     &   (L2**2*m2*((L1*(v1n + v1np1k)*
     &     (-2*g*(m1 + 2*m2)*Cos(th1n) + 2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &       L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &          v1np1k*v2np1k*(Cos(th1n - th2np1k) - 
     &             Cos(th1np1k - th2np1k)))))/(8.*(th1n - th1np1k)) + 
     &  (L2*m2*(v2n + v2np1k)*(L1*v1np1k*v2np1k*
     &        (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &       L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k)) +
     &       2*g*(-Cos(th2n) + Cos(th2np1k))))/(8.*(th2n - th2np1k)) + 
     &  (L1*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &       3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &          v2n*Cos(th1np1k - th2np1k)))*
     &     ((-4*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k)/2.) + 
     &            m2*(3*g*Sin((th1n + th1np1k - 2*th2n - 2*th2np1k)/
     &                  2.) + 
     &               L2*(v2n + v2np1k)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &               (3*L1*(-v1n - v1np1k)**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k))/4.)))/
     &        (3.*L1*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((-3*(g*(4*m1 + 5*m2)*
     &                  Sin((th1n + th1np1k - th1n*xi1 + th1np1k*xi1)/
     &                   2.) + 
     &                 m2*(3*g*
     &                   Sin((th1n + th1np1k - 2*th2n - 2*th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + 2*th2n*xi1 - 
     &                   2*th2np1k*xi1)/2.) + 
     &                   L2*(v2n + v2np1k - v2n*xi1 + v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 + xi1) - v1np1k*(1 + xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4.)))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                 th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))) - (3*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k + th1n*xi1 - 
     &                   th1np1k*xi1)/2.) + 
     &                 m2*(L2*(v2n + v2np1k + v2n*xi1 - v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &                    (3*L1*(v1n*(-1 - xi1) - v1np1k*(1 - xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4. + 
     &                    3*g*
     &                    Sin((th1n + th1np1k - 2*th2n - 2*th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - 2*th2n*xi1 + 
     &                   2*th2np1k*xi1)/2.))))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                 th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1)
     & ))))/18.))/12. + (L2*m2*(2*L2*(v2n + v2np1k) + 
     &    3*L1*v1n*Cos(th1n - th2n) + 3*L1*v1np1k*Cos(th1np1k - th2np1k)
     & )*((4*(3*L1*(m1 + 3*m2)*(v1n + v1np1k)**2*
     &             Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &            (9*L2*m2*(-v2n - v2np1k)**2*
     &               Sin(th1n + th1np1k - th2n - th2np1k))/4. + 
     &            9*g*(m1 + 2*m2)*
     &             Sin(th1n + th1np1k - th2n/2. - th2np1k/2.) - 
     &            3*g*(m1 + 6*m2)*Sin((th2n + th2np1k)/2.)))/
     &        (9.*L2*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((9*g*(m1 + 2*m2)*
     &                Sin(th1n + (th2n*(-1 - xi1))/2. + 
     &                  th1np1k*(1 - xi1) - (th2np1k*(1 - xi1))/2. + 
     &                  th1n*xi1) - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k + th2n*xi1 - th2np1k*xi1)/2.)
     & + 3*L1*(m1 + 3*m2)*(v1n + v1np1k + v1n*xi1 - v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 - xi1) - v2np1k*(1 - xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))/4.)/(L2*(8*m1 + 15*m2 - 9*m2*
     &                  Cos(th1n + th1np1k - th2n - th2np1k + th1n*xi1 - 
     &                    th1np1k*xi1 - th2n*xi1 + th2np1k*xi1))) + 
     &            (3*L1*(m1 + 3*m2)*
     &                (v1n + v1np1k - v1n*xi1 + v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k - th1n*xi1 + 
     &                   th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 + xi1) - v2np1k*(1 + xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))/4.-3*g*(m1+6*m2)*Sin((th2n + th2np1k -th2n*xi1+th2np1k*xi1)/
     &                  2.) + 9*g*(m1 + 2*m2)*
     &                Sin(th1n + (th2n*(-1 + xi1))/2. - th1n*xi1 + 
     &                  th1np1k*(1 + xi1) - (th2np1k*(1 + xi1))/2.))/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &               th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1))
     & )))/18.))/12.))/
     & (6.*((L1**2*(-2*g*(m1 + 2*m2)*Cos(th1n) + 
     &        2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &        L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) +
     &           v1np1k*v2np1k*
     &            (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k))))**2)/
     &   (16.*(th1n - th1np1k)**2) + 
     &  (L2**2*m2**2*(2*L2*(v2n + v2np1k) + 3*L1*v1n*Cos(th1n - th2n) +
     &        3*L1*v1np1k*Cos(th1np1k - th2np1k))**2)/144. + 
     &  (L1**2*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &        3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &           v2n*Cos(th1np1k - th2np1k)))**2)/144. + 
     &  (L2**2*m2**2*(L1*v1np1k*v2np1k*
     &         (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &  L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k)) + 
     &  2*g*(-Cos(th2n) + Cos(th2np1k)))**2)/(16.*(th2n - th2np1k)**2)))

        detm=m12*m24*m33*m41 - m12*m23*m34*m41 - m11*m24*m33*m42 
     & + m11*m23*m34*m42 - 
     &  m12*m24*m31*m43 + m11*m24*m32*m43 + m12*m21*m34*m43 
     & - m11*m22*m34*m43 + 
     &  m14*(m23*m32*m41 - m22*m33*m41 - m23*m31*m42 + m21*m33*m42 + 
     &     m22*m31*m43 - m21*m32*m43) + 
     &  (m12*m23*m31 - m11*m23*m32 - m12*m21*m33 + m11*m22*m33)*m44 + 
     &  m13*(-(m24*m32*m41) + m22*m34*m41 + m24*m31*m42 - m21*m34*m42 -
     &  m22*m31*m44 + m21*m32*m44)

      fdelta1= -((-th1n + th1np1k)/tdelta) + (v1n + v1np1k)/2. - 
     & (L1*(-2*g*(m1 + 2*m2)*Cos(th1n) + 2*g*(m1 + 2*m2)*Cos(th1np1k) +
     &  L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &   v1np1k*v2np1k*(Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k))))*
     & ((L1*(v1n + v1np1k)*(-2*g*(m1 + 2*m2)*Cos(th1n) + 
     &       2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &       L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &          v1np1k*v2np1k*(Cos(th1n - th2np1k) - 
     &             Cos(th1np1k - th2np1k)))))/(8.*(th1n - th1np1k)) + 
     &  (L2*m2*(v2n + v2np1k)*(L1*v1np1k*v2np1k*
     &        (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &       L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k)) +
     &       2*g*(-Cos(th2n) + Cos(th2np1k))))/(8.*(th2n - th2np1k)) + 
     &  (L1*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &       3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &          v2n*Cos(th1np1k - th2np1k)))*
     &     ((-4*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k)/2.) + 
     &            m2*(3*g*Sin((th1n + th1np1k - 2*th2n - 2*th2np1k)/
     &                  2.) + 
     &               L2*(v2n + v2np1k)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &               (3*L1*(-v1n - v1np1k)**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k))/4.)))/
     &        (3.*L1*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((-3*(g*(4*m1 + 5*m2)*
     &                  Sin((th1n + th1np1k - th1n*xi1 + th1np1k*xi1)/
     &                   2.) + 
     &                 m2*(3*g*
     &                   Sin((th1n + th1np1k - 2*th2n - 2*th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + 2*th2n*xi1 - 
     &                   2*th2np1k*xi1)/2.) + 
     &                   L2*(v2n + v2np1k - v2n*xi1 + v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 + xi1) - v1np1k*(1 + xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4.)))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                  th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))) - (3*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k + th1n*xi1 - 
     &                   th1np1k*xi1)/2.) + 
     &                 m2*(L2*(v2n + v2np1k + v2n*xi1 - v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &                    (3*L1*(v1n*(-1 - xi1) - v1np1k*(1 - xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4. + 
     &                    3*g*
     &                    Sin((th1n + th1np1k - 2*th2n - 2*th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - 2*th2n*xi1 + 
     &                   2*th2np1k*xi1)/2.))))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                 th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1)
     & ))))/18.))/12. + (L2*m2*(2*L2*(v2n + v2np1k) + 
     &   3*L1*v1n*Cos(th1n - th2n) + 3*L1*v1np1k*Cos(th1np1k - th2np1k)
     & )*((4*(3*L1*(m1 + 3*m2)*(v1n + v1np1k)**2*
     &             Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &            (9*L2*m2*(-v2n - v2np1k)**2*
     &               Sin(th1n + th1np1k - th2n - th2np1k))/4. + 
     &            9*g*(m1 + 2*m2)*
     &             Sin(th1n + th1np1k - th2n/2. - th2np1k/2.) - 
     &            3*g*(m1 + 6*m2)*Sin((th2n + th2np1k)/2.)))/
     &        (9.*L2*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((9*g*(m1 + 2*m2)*
     &                Sin(th1n + (th2n*(-1 - xi1))/2. + 
     &                  th1np1k*(1 - xi1) - (th2np1k*(1 - xi1))/2. + 
     &                  th1n*xi1) - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k + th2n*xi1 - th2np1k*xi1)/2.)
     & + 3*L1*(m1 + 3*m2)*(v1n + v1np1k + v1n*xi1 - v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 - xi1) - v2np1k*(1 - xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))/4.)/(L2*(8*m1 + 15*m2 - 9*m2*
     &                  Cos(th1n + th1np1k - th2n - th2np1k + th1n*xi1 - 
     &                    th1np1k*xi1 - th2n*xi1 + th2np1k*xi1))) + 
     &            (3*L1*(m1 + 3*m2)*
     &                (v1n + v1np1k - v1n*xi1 + v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k - th1n*xi1 + 
     &                   th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 + xi1) - v2np1k*(1 + xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))/4.-3*g*(m1+6*m2)*Sin((th2n+th2np1k - th2n*xi1 + th2np1k*xi1)/
     &                  2.) + 9*g*(m1 + 2*m2)*
     &                Sin(th1n + (th2n*(-1 + xi1))/2. - th1n*xi1 + 
     &                  th1np1k*(1 + xi1) - (th2np1k*(1 + xi1))/2.))/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1))
     & )))/18.))/12.))/
     &    (4.*(th1n - th1np1k)*((L1**2*
     &     (-2*g*(m1 + 2*m2)*Cos(th1n) + 2*g*(m1 + 2*m2)*Cos(th1np1k) +
     &        L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) +
     &           v1np1k*v2np1k*
     &            (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k))))**2)/
     &   (16.*(th1n - th1np1k)**2) + 
     &  (L2**2*m2**2*(2*L2*(v2n + v2np1k) + 3*L1*v1n*Cos(th1n - th2n) +
     &        3*L1*v1np1k*Cos(th1np1k - th2np1k))**2)/144. + 
     &  (L1**2*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &        3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &           v2n*Cos(th1np1k - th2np1k)))**2)/144. + 
     &  (L2**2*m2**2*(L1*v1np1k*v2np1k*
     &         (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &   L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k)) + 
     &  2*g*(-Cos(th2n) + Cos(th2np1k)))**2)/(16.*(th2n - th2np1k)**2)))

        fdelta2=-((-th2n + th2np1k)/tdelta) + (v2n + v2np1k)/2. - 
     & (L2*m2*(L1*v1np1k*v2np1k*(Cos(th1n - th2n) - Cos(th1n-th2np1k))+
     &  L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k)) + 
     &  2*g*(-Cos(th2n) + Cos(th2np1k)))*
     & ((L1*(v1n + v1np1k)*(-2*g*(m1 + 2*m2)*Cos(th1n) + 
     &       2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &       L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &          v1np1k*v2np1k*(Cos(th1n - th2np1k) - 
     &             Cos(th1np1k - th2np1k)))))/(8.*(th1n - th1np1k)) + 
     &  (L2*m2*(v2n + v2np1k)*(L1*v1np1k*v2np1k*
     &        (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &       L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k)) +
     &       2*g*(-Cos(th2n) + Cos(th2np1k))))/(8.*(th2n - th2np1k)) + 
     &  (L1*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &       3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &          v2n*Cos(th1np1k - th2np1k)))*
     &     ((-4*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k)/2.) + 
     &            m2*(3*g*Sin((th1n + th1np1k - 2*th2n - 2*th2np1k)/
     &                  2.) + 
     &               L2*(v2n + v2np1k)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &               (3*L1*(-v1n - v1np1k)**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k))/4.)))/
     &        (3.*L1*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((-3*(g*(4*m1 + 5*m2)*
     &                  Sin((th1n + th1np1k - th1n*xi1 + th1np1k*xi1)/
     &                   2.) + 
     &                 m2*(3*g*
     &                   Sin((th1n + th1np1k - 2*th2n - 2*th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + 2*th2n*xi1 - 
     &                   2*th2np1k*xi1)/2.) + 
     &                   L2*(v2n + v2np1k - v2n*xi1 + v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 + xi1) - v1np1k*(1 + xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4.)))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))) - (3*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k + th1n*xi1 - 
     &                   th1np1k*xi1)/2.) + 
     &                 m2*(L2*(v2n + v2np1k + v2n*xi1 - v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &                    (3*L1*(v1n*(-1 - xi1) - v1np1k*(1 - xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4. + 
     &                    3*g*
     &                    Sin((th1n + th1np1k - 2*th2n - 2*th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - 2*th2n*xi1 + 
     &                   2*th2np1k*xi1)/2.))))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1)
     & ))))/18.))/12. + (L2*m2*(2*L2*(v2n + v2np1k) + 
     &   3*L1*v1n*Cos(th1n - th2n) + 3*L1*v1np1k*Cos(th1np1k - th2np1k)
     & )*((4*(3*L1*(m1 + 3*m2)*(v1n + v1np1k)**2*
     &             Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &            (9*L2*m2*(-v2n - v2np1k)**2*
     &               Sin(th1n + th1np1k - th2n - th2np1k))/4. + 
     &            9*g*(m1 + 2*m2)*
     &             Sin(th1n + th1np1k - th2n/2. - th2np1k/2.) - 
     &            3*g*(m1 + 6*m2)*Sin((th2n + th2np1k)/2.)))/
     &        (9.*L2*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((9*g*(m1 + 2*m2)*
     &                Sin(th1n + (th2n*(-1 - xi1))/2. + 
     &                  th1np1k*(1 - xi1) - (th2np1k*(1 - xi1))/2. + 
     &                  th1n*xi1) - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k + th2n*xi1 - th2np1k*xi1)/2.) 
     & + 3*L1*(m1 + 3*m2)*(v1n + v1np1k + v1n*xi1 - v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 - xi1) - v2np1k*(1 - xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))/4.)/(L2*(8*m1 + 15*m2 - 9*m2*
     &                  Cos(th1n + th1np1k - th2n - th2np1k + th1n*xi1 - 
     &                    th1np1k*xi1 - th2n*xi1 + th2np1k*xi1))) + 
     &            (3*L1*(m1 + 3*m2)*
     &                (v1n + v1np1k - v1n*xi1 + v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k - th1n*xi1 + 
     &                   th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 + xi1) - v2np1k*(1 + xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))/4.-3*g*(m1+6*m2)*Sin((th2n+th2np1k - th2n*xi1 + th2np1k*xi1)/
     &                  2.) + 9*g*(m1 + 2*m2)*
     &                Sin(th1n + (th2n*(-1 + xi1))/2. - th1n*xi1 + 
     &                  th1np1k*(1 + xi1) - (th2np1k*(1 + xi1))/2.))/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &               th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1))
     & )))/18.))/12.))/
     & (4.*(th2n - th2np1k)*((L1**2*
     &     (-2*g*(m1 + 2*m2)*Cos(th1n) + 2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &        L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &           v1np1k*v2np1k*
     &            (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k))))**2)/
     &   (16.*(th1n - th1np1k)**2) + 
     &  (L2**2*m2**2*(2*L2*(v2n + v2np1k) + 3*L1*v1n*Cos(th1n - th2n) + 
     &        3*L1*v1np1k*Cos(th1np1k - th2np1k))**2)/144. + 
     &  (L1**2*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &        3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &           v2n*Cos(th1np1k - th2np1k)))**2)/144. + 
     &  (L2**2*m2**2*(L1*v1np1k*v2np1k*
     &         (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &   L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k)) + 
     &  2*g*(-Cos(th2n) + Cos(th2np1k)))**2)/(16.*(th2n - th2np1k)**2)))

       fdelta3=-((-v1n + v1np1k)/tdelta) - (4*
     & (g*(4*m1 + 5*m2)*Sin((th1n + th1np1k)/2.) + 
     &  m2*(3*g*Sin((th1n + th1np1k - 2*th2n - 2*th2np1k)/2.) + 
     & L2*(v2n + v2np1k)**2*Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     & (3*L1*(-v1n - v1np1k)**2*Sin(th1n + th1np1k - th2n - th2np1k))/4.
     & )))/(3.*L1*(8*m1+15*m2-9*m2*Cos(th1n+th1np1k-th2n- th2np1k)))+
     &   (5*((-3*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k - th1n*xi1 + 
     &            th1np1k*xi1)/2.) + 
     &       m2*(3*g*Sin((th1n + th1np1k - 2*th2n - 2*th2np1k - 
     &               th1n*xi1 + th1np1k*xi1 + 2*th2n*xi1 - 2*th2np1k*xi1
     & )/2.) + L2*(v2n + v2np1k - v2n*xi1 + v2np1k*xi1)**2*
     &           Sin((th1n + th1np1k - th2n - th2np1k - th1n*xi1 + 
     &               th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)/2.) + 
     &          (3*L1*(v1n*(-1 + xi1) - v1np1k*(1 + xi1))**2*
     &             Sin(th1n + th1np1k - th2n - th2np1k - th1n*xi1 + 
     &               th1np1k*xi1 + th2n*xi1 - th2np1k*xi1))/4.)))/
     &   (L1*(8*m1 + 15*m2 - 9*m2*
     &        Cos(th1n + th1np1k - th2n - th2np1k - th1n*xi1 + 
     &          th1np1k*xi1 + th2n*xi1 - th2np1k*xi1))) - 
     &(3*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k + th1n*xi1 - th1np1k*xi1)/
     &          2.) + m2*(L2*(v2n + v2np1k + v2n*xi1 - v2np1k*xi1)**2*
     &           Sin((th1n + th1np1k - th2n - th2np1k + th1n*xi1 - 
     &               th1np1k*xi1 - th2n*xi1 + th2np1k*xi1)/2.) + 
     &          (3*L1*(v1n*(-1 - xi1) - v1np1k*(1 - xi1))**2*
     &             Sin(th1n + th1np1k - th2n - th2np1k + th1n*xi1 - 
     &               th1np1k*xi1 - th2n*xi1 + th2np1k*xi1))/4. + 
     &          3*g*Sin((th1n + th1np1k - 2*th2n - 2*th2np1k + th1n*xi1-
     &               th1np1k*xi1 - 2*th2n*xi1 + 2*th2np1k*xi1)/2.))))/
     &   (L1*(8*m1 + 15*m2 - 9*m2*
     &        Cos(th1n + th1np1k - th2n - th2np1k + th1n*xi1 - 
     &          th1np1k*xi1 - th2n*xi1 + th2np1k*xi1)))))/18. - 
     &   (L1*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &  3*L2*m2*(v2np1k*Cos(th1n - th2n) + v2n*Cos(th1np1k - th2np1k)))*
     & ((L1*(v1n + v1np1k)*(-2*g*(m1 + 2*m2)*Cos(th1n) + 
     &       2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &       L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &          v1np1k*v2np1k*(Cos(th1n - th2np1k) - 
     &             Cos(th1np1k - th2np1k)))))/(8.*(th1n - th1np1k)) + 
     &  (L2*m2*(v2n + v2np1k)*(L1*v1np1k*v2np1k*
     &        (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &       L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k)) +
     &       2*g*(-Cos(th2n) + Cos(th2np1k))))/(8.*(th2n - th2np1k)) + 
     &  (L1*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &       3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &          v2n*Cos(th1np1k - th2np1k)))*
     &     ((-4*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k)/2.) + 
     &            m2*(3*g*Sin((th1n + th1np1k - 2*th2n - 2*th2np1k)/
     &                  2.) + 
     &               L2*(v2n + v2np1k)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &               (3*L1*(-v1n - v1np1k)**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k))/4.)))/
     &        (3.*L1*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((-3*(g*(4*m1 + 5*m2)*
     &                  Sin((th1n + th1np1k - th1n*xi1 + th1np1k*xi1)/
     &                   2.) + 
     &                 m2*(3*g*
     &                   Sin((th1n + th1np1k - 2*th2n - 2*th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + 2*th2n*xi1 - 
     &                   2*th2np1k*xi1)/2.) + 
     &                   L2*(v2n + v2np1k - v2n*xi1 + v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 + xi1) - v1np1k*(1 + xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4.)))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                 th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))) - (3*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k + th1n*xi1 - 
     &                   th1np1k*xi1)/2.) + 
     &                 m2*(L2*(v2n + v2np1k + v2n*xi1 - v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &                    (3*L1*(v1n*(-1 - xi1) - v1np1k*(1 - xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4. + 
     &                    3*g*
     &                    Sin((th1n + th1np1k - 2*th2n - 2*th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - 2*th2n*xi1 + 
     &                   2*th2np1k*xi1)/2.))))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &               th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1)
     & ))))/18.))/12. + (L2*m2*(2*L2*(v2n + v2np1k) + 
     &  3*L1*v1n*Cos(th1n - th2n) + 3*L1*v1np1k*Cos(th1np1k - th2np1k)
     & )*((4*(3*L1*(m1 + 3*m2)*(v1n + v1np1k)**2*
     &             Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &            (9*L2*m2*(-v2n - v2np1k)**2*
     &               Sin(th1n + th1np1k - th2n - th2np1k))/4. + 
     &            9*g*(m1 + 2*m2)*
     &             Sin(th1n + th1np1k - th2n/2. - th2np1k/2.) - 
     &            3*g*(m1 + 6*m2)*Sin((th2n + th2np1k)/2.)))/
     &        (9.*L2*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((9*g*(m1 + 2*m2)*
     &                Sin(th1n + (th2n*(-1 - xi1))/2. + 
     &                  th1np1k*(1 - xi1) - (th2np1k*(1 - xi1))/2. + 
     &                  th1n*xi1) - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k + th2n*xi1 - th2np1k*xi1)/2.)
     & + 3*L1*(m1 + 3*m2)*(v1n + v1np1k + v1n*xi1 - v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 - xi1) - v2np1k*(1 - xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))/4.)/(L2*(8*m1 + 15*m2 - 9*m2*
     &                  Cos(th1n + th1np1k - th2n - th2np1k + th1n*xi1 - 
     &                    th1np1k*xi1 - th2n*xi1 + th2np1k*xi1))) + 
     &            (3*L1*(m1 + 3*m2)*
     &                (v1n + v1np1k - v1n*xi1 + v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k - th1n*xi1 + 
     &                   th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 + xi1) - v2np1k*(1 + xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))/4.-3*g*(m1+6*m2)*Sin((th2n + th2np1k - th2n*xi1+th2np1k*xi1)/
     &                  2.) + 9*g*(m1 + 2*m2)*
     &                Sin(th1n + (th2n*(-1 + xi1))/2. - th1n*xi1 + 
     &                  th1np1k*(1 + xi1) - (th2np1k*(1 + xi1))/2.))/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &               th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1))
     & )))/18.))/12.))/
     &    (12.*((L1**2*(-2*g*(m1 + 2*m2)*Cos(th1n) + 
     &        2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &        L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) +
     &           v1np1k*v2np1k*
     &            (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k))))**2)/
     &   (16.*(th1n - th1np1k)**2) + 
     &  (L2**2*m2**2*(2*L2*(v2n + v2np1k) + 3*L1*v1n*Cos(th1n - th2n) +
     &        3*L1*v1np1k*Cos(th1np1k - th2np1k))**2)/144. + 
     &  (L1**2*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &        3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &           v2n*Cos(th1np1k - th2np1k)))**2)/144. + 
     &  (L2**2*m2**2*(L1*v1np1k*v2np1k*
     &         (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &     L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k)) + 
     & 2*g*(-Cos(th2n) + Cos(th2np1k)))**2)/(16.*(th2n - th2np1k)**2)))

       fdelta4=-((-v2n + v2np1k)/tdelta) + (4*
     & (3*L1*(m1 + 3*m2)*(v1n + v1np1k)**2*
     &   Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &(9*L2*m2*(-v2n - v2np1k)**2*Sin(th1n + th1np1k - th2n - th2np1k))/
     & 4. + 9*g*(m1 + 2*m2)*Sin(th1n + th1np1k - th2n/2. - th2np1k/2.)-
     &  3*g*(m1 + 6*m2)*Sin((th2n + th2np1k)/2.)))/
     & (9.*L2*(8*m1 + 15*m2 - 9*m2*Cos(th1n + th1np1k-th2n-th2np1k)))+
     & (5*((9*g*(m1 + 2*m2)*Sin(th1n + (th2n*(-1 - xi1))/2. + 
     &        th1np1k*(1 - xi1) - (th2np1k*(1 - xi1))/2. + th1n*xi1) - 
     &  3*g*(m1 + 6*m2)*Sin((th2n + th2np1k + th2n*xi1-th2np1k*xi1)/
     &        2.) + 3*L1*(m1 + 3*m2)*
     &      (v1n + v1np1k + v1n*xi1 - v1np1k*xi1)**2*
     &      Sin((th1n + th1np1k - th2n - th2np1k + th1n*xi1 - 
     &          th1np1k*xi1 - th2n*xi1 + th2np1k*xi1)/2.) + 
     &     (9*L2*m2*(v2n*(-1 - xi1) - v2np1k*(1 - xi1))**2*
     &        Sin(th1n + th1np1k - th2n - th2np1k + th1n*xi1 - 
     &          th1np1k*xi1 - th2n*xi1 + th2np1k*xi1))/4.)/
     &   (L2*(8*m1 + 15*m2 - 9*m2*
     &        Cos(th1n + th1np1k - th2n - th2np1k + th1n*xi1 - 
     &          th1np1k*xi1 - th2n*xi1 + th2np1k*xi1))) + 
     &  (3*L1*(m1 + 3*m2)*(v1n + v1np1k - v1n*xi1 + v1np1k*xi1)**2*
     &      Sin((th1n + th1np1k - th2n - th2np1k - th1n*xi1 + 
     &          th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)/2.) + 
     &     (9*L2*m2*(v2n*(-1 + xi1) - v2np1k*(1 + xi1))**2*
     &        Sin(th1n + th1np1k - th2n - th2np1k - th1n*xi1 + 
     &          th1np1k*xi1 + th2n*xi1 - th2np1k*xi1))/4. - 
     &  3*g*(m1 + 6*m2)*Sin((th2n + th2np1k - th2n*xi1 + th2np1k*xi1)/
     &        2.) + 9*g*(m1 + 2*m2)*
     &      Sin(th1n + (th2n*(-1 + xi1))/2. - th1n*xi1 + 
     &        th1np1k*(1 + xi1) - (th2np1k*(1 + xi1))/2.))/
     &   (L2*(8*m1 + 15*m2 - 9*m2*
     &        Cos(th1n + th1np1k - th2n - th2np1k - th1n*xi1 + 
     &          th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)))))/18. - 
     & (L2*m2*(2*L2*(v2n + v2np1k) + 3*L1*v1n*Cos(th1n - th2n) + 
     &  3*L1*v1np1k*Cos(th1np1k - th2np1k))*
     & ((L1*(v1n + v1np1k)*(-2*g*(m1 + 2*m2)*Cos(th1n) + 
     &       2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &       L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) +
     &          v1np1k*v2np1k*(Cos(th1n - th2np1k) - 
     &             Cos(th1np1k - th2np1k)))))/(8.*(th1n - th1np1k)) +
     &  (L2*m2*(v2n + v2np1k)*(L1*v1np1k*v2np1k*
     &        (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &       L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k)) +
     &       2*g*(-Cos(th2n) + Cos(th2np1k))))/(8.*(th2n - th2np1k)) + 
     &  (L1*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &       3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &          v2n*Cos(th1np1k - th2np1k)))*
     &     ((-4*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k)/2.) + 
     &            m2*(3*g*Sin((th1n + th1np1k - 2*th2n - 2*th2np1k)/
     &                  2.) + 
     &               L2*(v2n + v2np1k)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &               (3*L1*(-v1n - v1np1k)**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k))/4.)))/
     &        (3.*L1*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((-3*(g*(4*m1 + 5*m2)*
     &                  Sin((th1n + th1np1k - th1n*xi1 + th1np1k*xi1)/
     &                   2.) + 
     &                 m2*(3*g*
     &                   Sin((th1n + th1np1k - 2*th2n - 2*th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + 2*th2n*xi1 - 
     &                   2*th2np1k*xi1)/2.) + 
     &                   L2*(v2n + v2np1k - v2n*xi1 + v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1)/2.) + 
     &                   (3*L1*(v1n*(-1 + xi1) - v1np1k*(1 + xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - 
     &                   th2np1k*xi1))/4.)))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                 th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))) - (3*(g*(4*m1 + 5*m2)*Sin((th1n + th1np1k + th1n*xi1 - 
     &                   th1np1k*xi1)/2.) + 
     &                 m2*(L2*(v2n + v2np1k + v2n*xi1 - v2np1k*xi1)**2*
     &                   Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &                    (3*L1*(v1n*(-1 - xi1) - v1np1k*(1 - xi1))**2*
     &                   Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1))/4. + 
     &                    3*g*
     &                    Sin((th1n + th1np1k - 2*th2n - 2*th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - 2*th2n*xi1 + 
     &                   2*th2np1k*xi1)/2.))))/
     &             (L1*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k + 
     &                th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1)
     & ))))/18.))/12. + (L2*m2*(2*L2*(v2n + v2np1k) + 
     &    3*L1*v1n*Cos(th1n - th2n) + 3*L1*v1np1k*Cos(th1np1k - th2np1k)
     & )*((4*(3*L1*(m1 + 3*m2)*(v1n + v1np1k)**2*
     &             Sin((th1n + th1np1k - th2n - th2np1k)/2.) + 
     &            (9*L2*m2*(-v2n - v2np1k)**2*
     &               Sin(th1n + th1np1k - th2n - th2np1k))/4. + 
     &            9*g*(m1 + 2*m2)*
     &             Sin(th1n + th1np1k - th2n/2. - th2np1k/2.) - 
     &            3*g*(m1 + 6*m2)*Sin((th2n + th2np1k)/2.)))/
     &        (9.*L2*(8*m1 + 15*m2 - 
     &            9*m2*Cos(th1n + th1np1k - th2n - th2np1k))) + 
     &       (5*((9*g*(m1 + 2*m2)*
     &                Sin(th1n + (th2n*(-1 - xi1))/2. + 
     &                  th1np1k*(1 - xi1) - (th2np1k*(1 - xi1))/2. + 
     &                  th1n*xi1) - 
     &               3*g*(m1 + 6*m2)*
     &                Sin((th2n + th2np1k + th2n*xi1 - th2np1k*xi1)/2.)
     & + 3*L1*(m1 + 3*m2)*(v1n + v1np1k + v1n*xi1 - v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + 
     &                   th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 - xi1) - v2np1k*(1 - xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k + 
     &                   th1n*xi1 - th1np1k*xi1 - th2n*xi1 + th2np1k*xi1
     & ))/4.)/(L2*(8*m1 + 15*m2 - 9*m2*
     &                  Cos(th1n + th1np1k - th2n - th2np1k + th1n*xi1 - 
     &                    th1np1k*xi1 - th2n*xi1 + th2np1k*xi1))) + 
     &            (3*L1*(m1 + 3*m2)*
     &                (v1n + v1np1k - v1n*xi1 + v1np1k*xi1)**2*
     &                Sin((th1n + th1np1k - th2n - th2np1k - th1n*xi1 + 
     &                   th1np1k*xi1 + th2n*xi1 - th2np1k*xi1)/2.) + 
     &               (9*L2*m2*(v2n*(-1 + xi1) - v2np1k*(1 + xi1))**2*
     &                  Sin(th1n + th1np1k - th2n - th2np1k - 
     &                   th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1
     & ))/4.-3*g*(m1+6*m2)*Sin((th2n+th2np1k - th2n*xi1 + th2np1k*xi1)/
     &                  2.) + 9*g*(m1 + 2*m2)*
     &                Sin(th1n + (th2n*(-1 + xi1))/2. - th1n*xi1 + 
     &                  th1np1k*(1 + xi1) - (th2np1k*(1 + xi1))/2.))/
     &             (L2*(8*m1 + 15*m2 - 
     &                 9*m2*Cos(th1n + th1np1k - th2n - th2np1k - 
     &                th1n*xi1 + th1np1k*xi1 + th2n*xi1 - th2np1k*xi1))
     & )))/18.))/12.))/
     &    (12.*((L1**2*(-2*g*(m1 + 2*m2)*Cos(th1n) + 
     &        2*g*(m1 + 2*m2)*Cos(th1np1k) + 
     &        L2*m2*(v1n*v2n*(Cos(th1n - th2n) - Cos(th1np1k - th2n)) + 
     &           v1np1k*v2np1k*
     &            (Cos(th1n - th2np1k) - Cos(th1np1k - th2np1k))))**2)/
     &   (16.*(th1n - th1np1k)**2) + 
     &  (L2**2*m2**2*(2*L2*(v2n + v2np1k) + 3*L1*v1n*Cos(th1n - th2n) +
     &        3*L1*v1np1k*Cos(th1np1k - th2np1k))**2)/144. + 
     &  (L1**2*(2*L1*(m1 + 3*m2)*(v1n + v1np1k) + 
     &        3*L2*m2*(v2np1k*Cos(th1n - th2n) + 
     &           v2n*Cos(th1np1k - th2np1k)))**2)/144. + 
     &  (L2**2*m2**2*(L1*v1np1k*v2np1k*
     &         (Cos(th1n - th2n) - Cos(th1n - th2np1k)) + 
     &    L1*v1n*v2n*(Cos(th1np1k - th2n) - Cos(th1np1k - th2np1k)) + 
     &  2*g*(-Cos(th2n) + Cos(th2np1k)))**2)/(16.*(th2n - th2np1k)**2)))

       th1delta=(fdelta4*(m14*m23*m32 - m13*m24*m32 - m14*m22*m33 
     & + m12*m24*m33 + 
     & m13*m22*m34 - m12*m23*m34) + fdelta2*m14*m33*m42 - 
     &fdelta1*m24*m33*m42 - fdelta2*m13*m34*m42 + fdelta1*m23*m34*m42 -
     &fdelta2*m14*m32*m43 + fdelta1*m24*m32*m43 + fdelta2*m12*m34*m43 - 
     &fdelta1*m22*m34*m43 + (fdelta2*m13*m32 - fdelta1*m23*m32 - 
     & fdelta2*m12*m33 + fdelta1*m22*m33)*m44 + 
     &fdelta3*(-(m14*m23*m42) + m13*m24*m42 + m14*m22*m43 - m12*m24*m43-
     & m13*m22*m44 + m12*m23*m44))/detm
       th2delta=(fdelta4*(-(m14*m23*m31) + m13*m24*m31 + 
     & m14*m21*m33 - m11*m24*m33 - 
     &m13*m21*m34 + m11*m23*m34) - fdelta2*m14*m33*m41 + 
     &fdelta1*m24*m33*m41 + fdelta2*m13*m34*m41 - fdelta1*m23*m34*m41 +
     &fdelta2*m14*m31*m43 - fdelta1*m24*m31*m43 - fdelta2*m11*m34*m43 +
     &fdelta1*m21*m34*m43 - fdelta2*m13*m31*m44 + fdelta1*m23*m31*m44 +
     &fdelta2*m11*m33*m44 - fdelta1*m21*m33*m44 + 
     &fdelta3*(m14*m23*m41 - m13*m24*m41 - m14*m21*m43 + m11*m24*m43 +
     &m13*m21*m44 - m11*m23*m44))/detm
        v1delta=(fdelta4*(m14*m22*m31 - m12*m24*m31 - 
     &m14*m21*m32 + m11*m24*m32 + 
     &m12*m21*m34 - m11*m22*m34) + fdelta2*m14*m32*m41 - 
     &fdelta1*m24*m32*m41 - fdelta2*m12*m34*m41 + fdelta1*m22*m34*m41 - 
     &fdelta2*m14*m31*m42 + fdelta1*m24*m31*m42 + fdelta2*m11*m34*m42 - 
     &fdelta1*m21*m34*m42 + (fdelta2*m12*m31 - fdelta1*m22*m31 - 
     &fdelta2*m11*m32 + fdelta1*m21*m32)*m44 + 
     &fdelta3*(-(m14*m22*m41) + m12*m24*m41 + m14*m21*m42 - m11*m24*m42-
     &m12*m21*m44 + m11*m22*m44))/detm
        v2delta=(fdelta4*(-(m13*m22*m31) + m12*m23*m31 + 
     &m13*m21*m32 - m11*m23*m32 - 
     &m12*m21*m33 + m11*m22*m33) - fdelta2*m13*m32*m41 + 
     &fdelta1*m23*m32*m41 + fdelta2*m12*m33*m41 - fdelta1*m22*m33*m41 + 
     &fdelta2*m13*m31*m42 - fdelta1*m23*m31*m42 - fdelta2*m11*m33*m42 + 
     &fdelta1*m21*m33*m42 - fdelta2*m12*m31*m43 + fdelta1*m22*m31*m43 + 
     &fdelta2*m11*m32*m43 - fdelta1*m21*m32*m43 + 
     &fdelta3*(m13*m22*m41 - m12*m23*m41 - m13*m21*m42 + m11*m23*m42 + 
     &m12*m21*m43 - m11*m22*m43))/detm
       th1np1k=th1np1k+th1delta
       th2np1k=th2np1k+th2delta
       v1np1k=v1np1k+v1delta
       v2np1k=v2np1k+v2delta
       if(dabs(th1delta)+dabs(th2delta)+dabs(v1delta)+dabs(v2delta).lt.
     & 1.0d-10)go to 20
10     continue
       print*,'Failed to converge'
       stop
20     continue
       if(i.ge.printstep)then
       if(mod(i,printstep).eq.0)then
       uf1(nint(1.0d0*i/printstep)+1)=th1np1k
       uf2(nint(1.0d0*i/printstep)+1)=th2np1k
       vf1(nint(1.0d0*i/printstep)+1)=v1np1k
       vf2(nint(1.0d0*i/printstep)+1)=v2np1k
       energy(nint(1.0d0*i/printstep)+1)=1.0d0/6.0d0*(m1*L1**2*v1np1k**2
     & +m2*L2**2*v2np1k**2+3.0d0*m2*L1*v1np1k*(L1*v1np1k+
     & L2*v2np1k*dcos(th1np1k-th2np1k))
     & -3.0d0*(m1*g*L1*dcos(th1np1k)
     & +m2*g*L2*dcos(th2np1k)+2.0d0*m2*g*L1*dcos(th1np1k)))
       endif
       endif
       th1n=th1np1k
       th2n=th2np1k
       v1n=v1np1k
       v2n=v2np1k
100    continue

       do 200 i=1,nint(1.0d0*nsteps/printstep)+1
       write(11,500) time(i),uf1(i),uf2(i),vf1(i),vf2(i),energy(i)
200    continue
500    format(1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7)
       end
