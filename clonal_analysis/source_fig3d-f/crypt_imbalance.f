c---- non-neutral drift crypt dynamics
      implicit double precision (a-h,o-z)
c--------------------------------------------------------------------
c
c     program to determine least-square fit to the intestinal clonal data
c     based on biased drift dynamics using a 2d grid sampling of parameter
c     space. Programme samples though parameter space of the ISC 
c     loss/replacement rate (dltmin to dltmax) and degree of bias 
c     (delmin to delmax) for a given time-offset dt0 in units of weeks
c     
c     ntime - total number of time points
c     ngrid - grid dimension
c
c--------------------------------------------------------------------
      parameter (nmax=1000,nsize=100,ntime=4,ngrid=25,
     +   dltmin=2.0d0,
     1   dltmax=3.5d0,
     2   delmin=0.4d0,
     3   delmax=0.8d0,
     4   dt0=0.29d0)
      double precision draw(nmax,ntime),dpnexp(nsize,ntime),
     +   dpnop(nsize,ntime),dsq(ngrid,ngrid),dpnbar(nsize),
     1   dpnexo(nsize,ntime),dpnbro(nsize,ntime),
     2   dcnexp(0:nsize,ntime),dcnbar(0:nsize),
     3   dsq2(ngrid,ngrid),
     4   dcnexo(0:nsize,ntime),dcnbro(0:nsize,ntime),dtime(ntime)
      integer npoints(ntime)
c---- effective stem cell number
      ns=5
c---- raw data files - adjust for given data set
      open(unit=18,file='kras.dat')
c      open(unit=18,file='pi3k.dat')
c      open(unit=18,file='n1.dat')
c
c---- time points for clonal data in weeks without time offset
      dtime(1)=4.0d0/7.0d0
      dtime(2)=1.0d0
      dtime(3)=2.0d0
      dtime(4)=3.0d0
c      
c---- enter raw clonal data - draw - raw clonal data 
c                             npoints - number of clones
c
      i1=1
 20   read(18,*,end=40) (draw(i1,it),it=1,ntime)
         do 30 it=1,ntime
         if (draw(i1,it).ne.-1.0d0) npoints(it)=i1
 30      continue
      i1=i1+1
      goto 20
 40   close(unit=18)
c---- reset optimal value record for fitting parameters from clonal distribution
      dsqmin=1.0d9
      dltop=0.0d0
      delop=0.0d0
c---- and cumulative distribution
      dsqmin2=1.0d9
      dltop2=0.0d0
      delop2=0.0d0
      call expdat(ns,draw,dpnexp,dcnexp,npoints)
c---- loop through all trial values of loss/replacement rate and bias
         do 600 l1=1,ngrid
         dlt=dltmin+(dltmax-dltmin)*dble(l1-1)/dble(ngrid-1)
            do 550 l2=1,ngrid
            del=delmin+(delmax-delmin)*dble(l2-1)/dble(ngrid-1)
c---- reset the least-square parameters
            dsq(l1,l2)=0.0d0
            dsq2(l1,l2)=0.0d0
               do 300 it=1,ntime
               dlt1=dlt*(dtime(it)-dt0)
c---- get theoretical distribution for given parameter set
               call cdistd(ns,dlt1,del,dpnbar,dcnbar)
c---- determine least-square statistic for direct and cumulative distribution
                  do 200 l3=1,ns
                  if (dpnbar(l3).gt.1d-9) 
     +               dsq(l1,l2)=dsq(l1,l2)+
     1                  (dpnbar(l3)-dpnexp(l3,it))**2.0d0
                  if (dcnbar(l3).gt.1d-9) 
     +               dsq2(l1,l2)=dsq2(l1,l2)+
     1                  (dcnbar(l3)-dcnexp(l3,it))**2.0d0
 200              continue
 300           continue
c---- update least-square statistic for optimal distribution
            if (dsq(l1,l2).lt.dsqmin) then
               dsqmin=dsq(l1,l2)
               dltop=dlt
               delop=del
c---- store optimal experiment and model fit
                  do 320 it=1,ntime
                  dlt1=dlt*(dtime(it)-dt0)
                  call cdistd(ns,dlt1,del,dpnbar,dcnbar)
                     do 310 l3=1,ns
                     dpnexo(l3,it)=dpnexp(l3,it)
                     dpnbro(l3,it)=dpnbar(l3)
 310                 continue
 320              continue
            endif
c---- update least-square statistics for cumulative distribution
            if (dsq2(l1,l2).lt.dsqmin2) then
               dsqmin2=dsq2(l1,l2)
               dltop2=dlt
               delop2=del
c---- store optimal experiment and model
                  do 340 it=1,ntime
                  dlt1=dlt*(dtime(it)-dt0)
                  call cdistd(ns,dlt1,del,dpnbar,dcnbar)
                     do 330 l3=1,nsize
                     dcnexo(l3,it)=dcnexp(l3,it)
                     dcnbro(l3,it)=dcnbar(l3)
 330                 continue
 340             continue
           endif
 550       continue
 600    continue
c---- output optimal parameters
      write(*,'(2a5,f6.4,a5,f6.4,a5,f6.4,a5,f6.4)') 'dis ','lsq=',
     +   dsqmin,' lam=',dltop,' del=',delop,' lef=',dltop*(1.0d0-delop)
      write(*,'(2a5,f6.4,a5,f6.4,a5,f6.4,a5,f6.4)') 'cum ','lsq=',
     +     dsqmin2,'lam=',dltop2,'del=',delop2,' lef=',
     1     dltop*(1.0d0-delop2)
c---- output optimal distribution
      open(unit=19,file='dist_optimal.dat',status='replace')
         do 750 l1=1,ns
         write(19,'(i3,12f6.3)') l1,(dpnexo(l1,it),
     +      dsqrt(dpnbro(l1,it)*(1.0d0-dpnbro(l1,it))/dble(npoints(it)))
     1      ,dpnbro(l1,it),it=1,ntime)
 750     continue
c---- output optimal cumulative distribution
      open(unit=18,file='cdist_optimal.dat',status='replace')
      write(18,*) 0,1.0d0,1.0d0,1.0d0,1.0d0,1.0d0,1.0d0
         do 760 l1=1,ns
         write(18,*) l1,(dcnexo(l1,it),dcnbro(l1,it),it=1,ntime)
 760     continue
c---- output least-square measures across parameter space for mathematica
      open(unit=20,file='lsq.dat',status='replace')
      write(20,'(a19)') 'a=ListContourPlot[{'
         do 800 l1=1,ngrid-1
         write(20,'(a3,29(f10.5,a3),f10.5,a3)') '{',
     +      (dsq(l1,l2),',',l2=1,ngrid-1),
     1      dsq(l1,ngrid),'},'
 800     continue
      write(20,'(a3,29(f10.5,a3),f10.5,a3)') '{',
     +   (dsq(ngrid,l2),',',l2=1,ngrid-1),
     1   dsq(ngrid,ngrid),'}},'
      write(20,'(a120)') 'LabelStyle -> {FontFamily -> "Helvetica", 
     +   FontSize -> 14}, PlotLabel -> "least-square statistic",' 
      write(20,'(a120)') 'FrameLabel -> {"Drift bias, \[Delta]", 
     +   "loss/replacement rate, \[Lambda]"},'
      write(20,'(a120)') 'PlotLegends -> BarLegend[All, 20],'
      write(20,'(a31,f6.3,a3,f6.4,a5,f6.3,a3,f6.4,a3)') 
     +   'Contours -> 20, DataRange -> {{',
     1   delmin,',',delmax,'},{',dltmin,',',dltmax,'}}]'
      close(unit=18)
      close(unit=19)
      close(unit=20)
c---- output fit to averages clone sizes
      open(unit=19,file='average_optimal.dat',status='replace')
         do 820 it=1,ntime
         davexp=0.0d0
         davbar=0.0d0   
            do 810 l1=1,ns
            davexp=davexp+dble(l1)*dpnexp(l1,it)
            davbar=davbar+dble(l1)*dpnbro(l1,it)
 810        continue
            write(19,*) dtime(it),davexp,davexp/dsqrt(dble(npoints(it)))
     +         ,davbar
 820     continue
      close(unit=19)            
      stop
      end
c
c
c
c
c
      subroutine cdistd(Nstem,dlt,del,dpnbar,dcnbar)
      implicit double precision (a-h,o-z)
      parameter (nsize=100,dpi=dacos(0.0d0)*2.0d0)
      double precision dpn(0:nsize),dpnbar(nsize),dcnbar(0:nsize)
      if ((Nstem.lt.2).or.(Nstem.gt.nsize)) then
         write(*,*) 'Nstem error ',Nstem
         stop
      endif
c---- zero distribution
         do 100 l1=0,Nstem
         dpn(l1)=0.0d0
 100     continue
      dmu=dsqrt(1.0d0-del**2.0d0)
      dmul=dsqrt((1.0d0-del)/(1.0d0+del))
c---- loop over k sum values
         do 200 lk=1,Nstem-1
         dfk=2.0d0*(1.0d0/dmu-1.0d0)+
     +      (2.0d0*dsin(dpi*dble(lk)/2.0d0/dble(Nstem)))**2.0d0
         d1=(1.0d0-dexp(-dlt*dmu*dfk))
     +      *(dsin(dpi*dble(lk)/dble(Nstem)))**2.0d0/dfk
         dpn(0)=dpn(0)+2.0d0*d1/dble(Nstem)*dmul
         dpn(Nstem)=dpn(Nstem)+dble((-1)**(lk+1))*2.0d0*d1/dble(Nstem)*
     +      dmul**(dble(1-Nstem))
            do 150 n1=1,Nstem-1
            dpn(n1)=dpn(n1)+2.0d0*dsin(dpi*dble(lk)/dble(Nstem))
     +         *dsin(dpi*dble(lk*n1)/dble(Nstem))*dexp(-dlt*dmu*dfk)
     1         /dble(Nstem)*dmul**(dble(1.0d0-n1))
 150        continue
 200     continue
c---- obtain surviving distribution
         do 300 l1=1,Nstem
         dpnbar(l1)=dpn(l1)/(1.0d0-dpn(0))
 300     continue
c---- get cumulative distribution
      dcnbar(0)=1.0d0
         do 400 l1=1,Nstem
         dcnbar(l1)=dcnbar(l1-1)-dpnbar(l1)
 400     continue
c---- checksum
      if (dabs(dcnbar(Nstem)).gt.1.0d-7) then
         write(*,*) 'distribution err2',dcnbar(Nstem)
         stop
      endif
      return
      end
c
c
c
c     
c
      subroutine expdat(ns,draw,dpnexp,dcnexp,npoints)
c--------------------------------------------------------------------
c     determine experimental distributions
c
c     input
c     ns - effective stem cell number
c     draw - raw data
c     npoints - number of clones
c
c     output
c     dpnexp - clone size distribution
c     dcnexp - cumulative clone size distribution
c--------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      parameter (nmax=1000,nsize=100,ntime=4)
      double precision draw(nmax,ntime),dpnexp(nsize,ntime),
     +   dcnexp(0:nsize,ntime)
      integer npoints(ntime)
c
c---- zero out distribution
         do 200 it=1,ntime
            do 100 l1=1,ns
            dpnexp(l1,it)=0.0d0
 100        continue
c---- set clone counter
         ncl=0
c---- loop over clones and set size by effective cell number
            do 120 l1=1,npoints(it)
            nval=int(draw(l1,it)/360.0d0*dble(ns))+1
c            nval=int((draw(l1,it)+360.0d0/(2.0d0*dble(ns)))/360.0d0
c     +         *dble(ns))
c---- checksum
            if (nval.gt.ns) then
               write(*,*) 'Nstem error'
               stop
            endif
c---- input clones and update clone counter
            if (nval.gt.0) then
               dpnexp(nval,it)=dpnexp(nval,it)+1.0d0
               ncl=ncl+1
            endif
 120        continue
         npoints(it)=ncl
c---- define distribution and cumulative distribution
         dcnexp(0,it)=1.0d0
            do 130 l1=1,ns
c---- normalize distribution
            dpnexp(l1,it)=dpnexp(l1,it)/dble(npoints(it))
            dcnexp(l1,it)=dcnexp(l1-1,it)-dpnexp(l1,it)
 130        continue
c---- checksum
         if (dabs(dcnexp(ns,it)).gt.1.0d-7) then
            write(*,*) 'distribution err',dcnexp(ns,it),it
            stop
         endif
 200     continue
      return
      end
