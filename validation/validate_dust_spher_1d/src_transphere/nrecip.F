c     -------------------------------------------------------------
c               HUNT FOR LOCATION OF A NUMBER IN A TABLE
c     -------------------------------------------------------------
      subroutine hunt(xx,n,x,jlo)
      integer jlo,n
      doubleprecision x,xx(n)
      integer inc,jhi,jm
      logical ascnd
      ascnd=xx(n).gt.xx(1)
      if(jlo.le.0.or.jlo.gt.n)then
        jlo=0
        jhi=n+1
        goto 3
      endif
      inc=1
      if(x.ge.xx(jlo).eqv.ascnd)then
1       jhi=jlo+inc
        if(jhi.gt.n)then
          jhi=n+1
        else if(x.ge.xx(jhi).eqv.ascnd)then
          jlo=jhi
          inc=inc+inc
          goto 1
        endif
      else
        jhi=jlo
2       jlo=jhi-inc
        if(jlo.lt.1)then
          jlo=0
        else if(x.lt.xx(jlo).eqv.ascnd)then
          jhi=jlo
          inc=inc+inc
          goto 2
        endif
      endif
3     if(jhi-jlo.eq.1)return
      jm=(jhi+jlo)/2
      if(x.gt.xx(jm).eqv.ascnd)then
        jlo=jm
      else
        jhi=jm
      endif
      goto 3
      end
C  (C) Copr. 1986-92 Numerical Recipes Software =v1.9"217..


c     --------------------------------------------------------------
c                Brent's algorithm for root finding
c     --------------------------------------------------------------
      FUNCTION zbrent(func,x1,x2,tol)
      INTEGER ITMAX
      DOUBLEPRECISION zbrent,tol,x1,x2,func,EPSS
      EXTERNAL func
      PARAMETER (ITMAX=100,EPSS=3.d-8)
      INTEGER iter
      DOUBLEPRECISION a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
      a=x1
      b=x2
      fa=func(a)
      fb=func(b)
      if((fa.gt.0..and.fb.gt.0.).or.(fa.lt.0..and.fb.lt.0.))stop
     *'root must be bracketed for zbrent'
      c=b
      fc=fb
      do 11 iter=1,ITMAX
        if((fb.gt.0..and.fc.gt.0.).or.(fb.lt.0..and.fc.lt.0.))then
          c=a
          fc=fa
          d=b-a
          e=d
        endif
        if(abs(fc).lt.abs(fb)) then
          a=b
          b=c
          c=a
          fa=fb
          fb=fc
          fc=fa
        endif
        tol1=2.*EPSS*abs(b)+0.5*tol
        xm=.5*(c-b)
        if(abs(xm).le.tol1 .or. fb.eq.0.)then
          zbrent=b
          return
        endif
        if(abs(e).ge.tol1 .and. abs(fa).gt.abs(fb)) then
          s=fb/fa
          if(a.eq.c) then
            p=2.*xm*s
            q=1.-s
          else
            q=fa/fc
            r=fb/fc
            p=s*(2.*xm*q*(q-r)-(b-a)*(r-1.))
            q=(q-1.)*(r-1.)*(s-1.)
          endif
          if(p.gt.0.) q=-q
          p=abs(p)
          if(2.*p .lt. min(3.*xm*q-abs(tol1*q),abs(e*q))) then
            e=d
            d=p/q
          else
            d=xm
            e=d
          endif
        else
          d=xm
          e=d
        endif
        a=b
        fa=fb
        if(abs(d) .gt. tol1) then
          b=b+d
        else
          b=b+sign(tol1,xm)
        endif
        fb=func(b)
11    continue
      stop 'zbrent exceeding maximum iterations'
      zbrent=b
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software =v1.9"217..


c     --------------------------------------------------------------
c                      MAKE INDEXED FILENAME
c
c     This routine makes filenames of the type dusttemp_4.dat etc.
c     --------------------------------------------------------------
      subroutine make_indexed_filename(base,index,ext,filename)
      implicit none
      character*80 base
      integer index
      character*80 ext
      character*80 filename
      character*12 ch
c     
      if((index.lt.0).or.(index.ge.1000)) then
          write(*,*) 'ERROR in make_indexed_filename()'
          stop 729
      endif
      if(index.lt.10) then
          write(ch,11) index
 11       format(I1)
      elseif(index.lt.100) then
          write(ch,12) index
 12       format(I2)
      elseif(index.lt.1000) then
          write(ch,13) index
 13       format(I3)
      endif
      filename = base(1:lnblnk(base))//
     %           ch(1:lnblnk(ch))//
     %           ext(1:lnblnk(ext))
c
      end


c     --------------------------------------------------------------
c                   NUMERICAL RECIPES ROUTINE: TRIDAG
c     --------------------------------------------------------------
      SUBROUTINE tridag(a,b,c,r,u,n)
      INTEGER n,NMAX
      DOUBLEPRECISION a(n),b(n),c(n),r(n),u(n)
      PARAMETER (NMAX=500)
      INTEGER j
      DOUBLEPRECISION bet,gam(NMAX)
      if(b(1).eq.0.)stop 'tridag: rewrite equations'
      bet=b(1)
      u(1)=r(1)/bet
      do 11 j=2,n
        gam(j)=c(j-1)/bet
        bet=b(j)-a(j)*gam(j)
        if(bet.eq.0.)stop 'tridag failed'
        u(j)=(r(j)-a(j)*u(j-1))/bet
11    continue
      do 12 j=n-1,1,-1
        u(j)=u(j)-gam(j+1)*u(j+1)
12    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software =v1.9"217..
