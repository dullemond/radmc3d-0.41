c     ==============================================================
c         MODULE FOR INTEGRATION OF THE FORMAL TRANSFER EQUATION
c                       IN SPHERICAL COORDINATES
c
c     This module integrates the transfer equation for given opacity
c     and source function in spherical coordinates. One first needs
c     to call the subroutine create_rays() in order to set up the
c     set of rays and compute the crossings with the radial shells.
c     Then one can call do_transfer() to compute the transfer 
c     integral for all frequencies. This module uses internal
c     common blocks in common_transfer.h. 
c
c     MODIFICATION (14-03-02): From now on this module is completely
c         frequency-ignorant. So you call these routines for each
c         frequency bin separately. This is to save on memory.
c
c     DANGER: Note that the srcarr and alparr in do_transfer() are
c         arrays from 0:FRSIZE_R !! Not from 1:FRSIZE_R.
c     ==============================================================

#include "main.h"


c     --------------------------------------------------------------
c                       CREATE THE TANGENT RAYS
c
c     This subroutine sets up the commons used by the tangent ray
c     transfer integrator. It computes the crossings of the rays
c     with the radial shells, and the corresponding mu points and
c     integration weights. It must be called before the do_transfer
c     is called.
c     --------------------------------------------------------------
      subroutine create_rays(ncst,ncex,ncnr,nr,rr,rstar,itypemuwgt)
      implicit none
      integer ncst,ncex,ncnr,nr,itypemuwgt
      doubleprecision rstar,rr(FRSIZE_R)
c
#include "common_transfer.h"
c
      doubleprecision r,b,dr,dummy,s
      doubleprecision mui(0:FRSIZE_RAYS)
      integer is,iss,ir,iray,iraydum,ifound,imu,ihit,irlow
      integer irefine
c
c     Store numbers in commons
c
      rays_ncst = ncst
      rays_ncex = ncex
      rays_ncnr = ncnr
c
c     Check array sizes
c
      if(nr.gt.FRSIZE_R) then
          write(*,*) 'Weird error: R-grid too big'
          stop 13
      endif
      if(ncst+ncex+nr*ncnr.gt.FRSIZE_RAYS) then
          write(*,*)'ERROR: Too many rays...'
          stop 13
      endif
c
c     Check if star not too big
c
      if(rstar.ge.rr(1)) then
          write(*,*) 'Star .ge. R_in!!'
          stop 13
      endif
c
c     Reset stuff
c
      iray = 1
      if(rays_ncst.lt.1) stop 45645
c
c     First the rays hitting the star
c
      do iraydum=1,rays_ncst
          rays_b(iray) = (iraydum-0.5d0) * (rstar/(1.d0*rays_ncst))
          iray = iray + 1
      enddo
c
c     Then some extra rays
c
      do iraydum=1,rays_ncex
          rays_b(iray) = rstar + (iraydum-0.5d0) * 
     %             ((rr(1)-rstar)/(rays_ncex*1.d0))
          iray = iray + 1
      enddo
c
c     Then the tangent rays 
c
      do ir=1,nr-1 
          dr = rr(ir+1)-rr(ir)
          do iraydum=1,rays_ncnr
              rays_b(iray) = rr(ir)+dr*(iraydum-1)/(rays_ncnr)
              if(iraydum.eq.1) then
                  rays_itang(iray) = 1
              else
                  rays_itang(iray) = 0
              endif
              iray = iray + 1
          enddo
      enddo
c
c     Total nr of rays
c
      rays_nr = iray-1
c
c     Now find the corresponding R-grid crossings
c
      do iray=1,rays_nr
          is = 1
c
c         First from outwards inwards
c
          irlow = nr+1
          do ir=nr,1,-1
              if(rr(ir).gt.rays_b(iray)) then
                  oneray_ir(is) = ir
                  oneray_dr(is) = 0.d0
                  oneray_s(is)  = 
     %                 -sqrt(rr(ir)**2-rays_b(iray)**2)
                  oneray_node(is) = 1
                  irlow = ir
                  is = is + 1
              endif
          enddo
c
c         If star?
c
          if(rays_b(iray).le.rstar) then
              rays_istar(iray)  = 1
              oneray_ir(is)  = 0
              oneray_dr(is)  = 0.d0
              oneray_s(is)   = sqrt(rstar**2-rays_b(iray)**2)
              oneray_node(is) = 0
              is = is + 1
          else
              rays_istar(iray) = 0
          endif
c
c         Add point at s=0
c
          if((rays_itang(iray).eq.0).and.
     %       (rays_b(iray).gt.rr(1))) then
              oneray_ir(is)  = irlow-1
              oneray_dr(is)  = (rays_b(iray)-rr(irlow-1))/
     %                            (rr(irlow)-rr(irlow-1))
              oneray_s(is)   = 0.d0
              oneray_node(is) = 0
              is = is + 1
          endif
c
c         Now from inwards outwards
c
          do ir=1,nr
              if(rr(ir).ge.rays_b(iray)) then
                  oneray_ir(is) = ir
                  oneray_dr(is) = 0.d0
                  oneray_s(is)  = 
     %                 sqrt(rr(ir)**2-rays_b(iray)**2)
                  oneray_node(is) = 1
                  is = is + 1
              endif
          enddo
          rays_ns(iray) = is - 1
c
c         Now copy into real arrays
c
          if(rays_nrefs.gt.SREF) then
              write(*,*) 'Sorry, cannot refine more than ',SREF
              stop 13
          endif
          iss=1
          do is=1,rays_ns(iray)
              rays_ir(iss,iray)   = oneray_ir(is)
              rays_dr(iss,iray)   = oneray_dr(is)
              rays_s(iss,iray)    = oneray_s(is)
              rays_node(iss,iray) = oneray_node(is)
              iss = iss + 1
              if((is.lt.rays_ns(iray)).and.(rays_nrefs.gt.0)) then
                  if((oneray_ir(is).ge.1).and.
     %               (oneray_ir(is+1).ge.1).and.
     %               (oneray_s(is)*oneray_s(is+1).gt.0.d0)) then
                      do irefine=1,rays_nrefs
                          s = oneray_s(is)+irefine*
     %                         (oneray_s(is+1)-oneray_s(is))/
     %                         (rays_nrefs+1)
                          r = sqrt(rays_b(iray)**2+s**2)
                          rays_ir(iss,iray) = min(oneray_ir(is),
     %                                            oneray_ir(is+1))
                          ir = rays_ir(iss,iray)
                          rays_dr(iss,iray) = (r-rr(ir))/
     %                                     (rr(ir+1)-rr(ir))
                          rays_s(iss,iray)  = s
                          rays_node(iss,iray) = 0
                          if(abs(rays_dr(iss,iray)).gt.1.d0) then
                              write(*,*) oneray_ir(is),
     %                                            oneray_ir(is+1)
                              write(*,*) r,rr(ir),rr(ir+1)
                              write(*,*) s,oneray_s(is),oneray_s(is+1)
                              stop 23721
                          endif
                          iss = iss + 1
                      enddo
                  endif
              endif
          enddo
          rays_ns(iray) = iss - 1
c############################
          if(rays_ir(iss-1,iray).eq.0) stop 1453
c############################
c
c         Check
c
          do is=2,rays_ns(iray) 
              if(rays_s(is,iray).le.rays_s(is-1,iray)) stop 88421
          enddo
c
c         Just for safety reset rays_imu ---> -100000
c
          do is=1,rays_ns(iray)
              rays_imu(is,iray)=-10000
          enddo
c
      enddo
c
c     Now create the mu at every point in R
c
      do ir=1,nr
          ifound=0
          r = rr(ir)
          node_nmu(ir) = 1
          node_iray(0,ir)=-1
          node_mu(0,ir)=0.d0
          do iray=rays_nr,1,-1
              b = rays_b(iray)
              if(b.eq.r) then
                  node_mu(0,ir)=0.d0
                  node_iray(0,ir)=iray
                  ifound=ifound+1
                  ihit=0
                  do is=1,rays_ns(iray)
                      if((rays_ir(is,iray).eq.ir).and.
     %                   (rays_dr(is,iray).eq.0.d0).and.
     %                   (rays_node(is,iray).eq.1).and.
     %                   (ihit.eq.0)) then
                          rays_imu(is,iray)=0
                          ihit=1
                      endif
                  enddo
              elseif(b.lt.r) then
                  node_mu(node_nmu(ir),ir)=sqrt(1.d0-(b/r)**2)
                  node_mu(-node_nmu(ir),ir)=-sqrt(1.d0-(b/r)**2)
                  node_iray(node_nmu(ir),ir)=iray
                  node_iray(-node_nmu(ir),ir)=iray
                  ihit=0
                  do is=1,rays_ns(iray)
                      if((rays_ir(is,iray).eq.ir).and.
     %                   (rays_dr(is,iray).eq.0.d0).and.
     %                   (rays_node(is,iray).eq.1).and.
     %                     (ihit.eq.0)) then
                          rays_imu(is,iray)=-node_nmu(ir)
                          ihit=1
                      elseif((rays_ir(is,iray).eq.ir).and.
     %                   (rays_dr(is,iray).eq.0.d0).and.
     %                   (rays_node(is,iray).eq.1).and.
     %                     (ihit.eq.1)) then
                          rays_imu(is,iray)=node_nmu(ir)
                          ihit=2
                      endif
                  enddo
                  if(ihit.ne.2) stop 23521
                  node_nmu(ir)=node_nmu(ir)+1
              endif
          enddo
          node_nmu(ir) = node_nmu(ir)-1
          if((ifound.ne.1).and.(ir.lt.nr)) then
              write(*,*) 'Error: did not find unique mu=0 ray'
              write(*,*) '  (in principle not wrong, but check)'
              write(*,*) '  ir,iray=',ir,iray
              stop 7675
          endif
      enddo
c
c     Find mu integration weights
c
      call compute_mu_weights(itypemuwgt,nr,rr,rstar)
c     
c     Do a safety check
c
      do iray=1,rays_nr
          do is=1,rays_ns(iray)
              if((rays_imu(is,iray).eq.-10000).and.
     %        (rays_dr(is,iray).eq.0.d0).and.
     %        (rays_ir(is,iray).ne.0)) then
                  write(*,*) 'Internal error while creating rays'
                  write(*,*) 'iray,is,ns,ir=',iray,is,rays_ns(iray),
     %              rays_ir(is,iray)
                  stop 3245
              endif
          enddo
      enddo
c
      end



c     --------------------------------------------------------------
c                        COMPUTE THE MU-WEIGHTS
c     --------------------------------------------------------------
      subroutine compute_mu_weights(itype,nr,rr,rstar)
      implicit none
c
      integer itype,nr
      doubleprecision rr(FRSIZE_R),rstar
c
#include "common_transfer.h"
c
      doubleprecision r,b,dr,dummy,rm1,rm2,eps,rm
      doubleprecision mui(0:FRSIZE_RAYS)
      integer is,ir,iray,iraydum,ifound,imu,ihit
c
      do ir=1,nr
c
          r = rr(ir)
c
c         Find the mu-interfaces. First the regular ones
c
          if(itype.eq.1) then
              do imu=0,node_nmu(ir)-rays_ncst-1
                  mui(imu) = sqrt(1.d0-sqrt(
     %                 (1.d0-node_mu(imu,ir)**2)*
     %                 (1.d0-node_mu(imu+1,ir)**2)))
              enddo
          elseif(itype.eq.2) then
              do imu=0,node_nmu(ir)-rays_ncst-1
                  mui(imu) = (node_mu(imu,ir)+node_mu(imu+1,ir))/2.
              enddo
          elseif(itype.eq.3) then
              do imu=0,node_nmu(ir)-rays_ncst-1
                  mui(imu) = ((node_mu(imu,ir)**0.5+
     %                   node_mu(imu+1,ir)**0.5)/2.)**2
              enddo
          elseif(itype.eq.4) then
              do imu=4,node_nmu(ir)-rays_ncst-1
                  mui(imu) = sqrt(1.d0-sqrt(
     %                 (1.d0-node_mu(imu,ir)**2)*
     %                 (1.d0-node_mu(imu+1,ir)**2)))
              enddo
              do imu=0,3
                  mui(imu) = ((node_mu(imu,ir)**0.5+
     %                   node_mu(imu+1,ir)**0.5)/2.)**2
              enddo
          elseif(itype.eq.5) then
              do imu=0,node_nmu(ir)-rays_ncst-1
                  rm1      = sqrt(1.d0-node_mu(imu,ir)**2)
                  rm2      = sqrt(1.d0-node_mu(imu+1,ir)**2)
                  eps      = 0.95
                  rm       = exp(eps*log(rm1)+(1.d0-eps)*log(rm2))
                  mui(imu) = sqrt(1.d0-rm**2)
              enddo
          else
              write(*,*) 'ERROR: Do not know mu_wgt_itype=',
     %             itype
              stop 55
          endif
c
c         Then the star ones. Make sure that the interface of the 
c         outermost ray is exactly tangent to the stellar surface.
c
          if(rays_ncst.ge.1) then
              imu = node_nmu(ir)-rays_ncst
              mui(imu) = sqrt(1.d0-(rstar/r)**2)
              do imu=node_nmu(ir)-rays_ncst+1,node_nmu(ir)-1
                  mui(imu) = 0.5 * ( node_mu(imu,ir) + 
     %                               node_mu(imu+1,ir) )
              enddo
          endif
c
c         Now the first and last ones
c
          if(node_iray(0,ir).lt.0) then
              mui(0)  = 0.d0
          endif
          mui(node_nmu(ir))  = 1.d0
c
c         From the mu interfaces, find the mu weights
c
          do imu=1,node_nmu(ir)
              node_dmu(imu,ir) = mui(imu)-mui(imu-1)
              node_dmu(-imu,ir) = node_dmu(imu,ir)
          enddo
          if(node_iray(0,ir).gt.0) then
              node_dmu(0,ir)=2*mui(0)
              if(mui(0).eq.0.d0) then 
                  write(*,*) 'Weird error while making dmu'
                  stop 5676
              endif
          else
              node_dmu(0,ir)=0.d0
          endif
c
c         Test the mu-weights
c
          dummy=0.d0
          do imu=-node_nmu(ir),node_nmu(ir)
              if(node_dmu(imu,ir).lt.0.d0) then
                  write(*,*) 'Internal error while making dmu:'
                  write(*,*) '  Negative dmu at imu=',imu,' ir=',ir
                  stop 9321
              endif
              dummy = dummy + node_dmu(imu,ir)
          enddo
          if(abs(dummy-2.d0).gt.1.d-4) then
              write(*,*) 'Internal error while making dmu:'
              write(*,*) '  int dmu neq 2'
              write(*,*) '  Answer=',dummy
              write(*,*) '  at ir=',ir,' with nmu=',node_nmu(ir)
              write(*,*) '  mu= ',(node_mu(imu,ir),imu=-node_nmu(ir),
     %                            node_nmu(ir))
              write(*,*) '  mui=',(mui(imu),imu=0,node_nmu(ir))
              write(*,*) '  dmu=',(node_dmu(imu,ir),imu=-node_nmu(ir),
     %                            node_nmu(ir))
              stop 4564
          endif
      enddo
c
      end



c     --------------------------------------------------------------
c                  DO THE TRANSFER FOR ONE FREQUENCY
c
c     This is the main routine for this module: it integrates the
c     formal radiative transfer equation along all rays through the
c     spherically symmetric cloud. It does so for a given frequency.
c
c     DANGER: Note that the srcarr and alparr in do_transfer() are
c         arrays from 0:FRSIZE_R !! Not from 1:FRSIZE_R.
c     --------------------------------------------------------------
      subroutine do_transfer(nr,rr,starint,bcout,srcarr,alparr)
      implicit none
c
      integer ir,imu,iray,is,nr
      doubleprecision intens,rr(FRSIZE_R)
      doubleprecision starint,bcout
      doubleprecision srcarr(0:FRSIZE_R)
      doubleprecision alparr(0:FRSIZE_R)
c
#include "common_transfer.h"
c
c     Check if the tangent rays have been set up
c
      if(rays_nr.lt.2) then
          write(*,*) 'ERROR: Called do_transfer() before having'
          write(*,*) '    called setup_rays().'
          stop 13
      endif
c
c     Safety measure:
c
      do ir=1,nr
          do imu=-node_nmu(ir),node_nmu(ir)
              node_intensity(imu,ir)=-1d99
          enddo
          do imu=node_nmu(ir)+1,FRSIZE_RAYS
              node_intensity(imu,ir)=1d99
              node_intensity(-imu,ir)=1d99
          enddo
      enddo
c
c     Do transfer ray by ray
c
      do iray=1,rays_nr
c
c         Set outer BC
c
          intens = bcout
c
c         Now do the transfer for this ray
c     
          call rayint(nr,rr,srcarr,alparr,
     %                starint,iray,intens)
c
c         Put intensities into appropriate node arrays
c
          do is=1,rays_ns(iray)
              if(rays_dr(is,iray).eq.0.d0) then
                  ir  = rays_ir(is,iray)
                  imu = rays_imu(is,iray)
                  node_intensity(imu,ir) = 
     %                 rays_intensity(is,iray)
              endif
          enddo
      enddo
c
c     Check if some things have been skipped
c
      do ir=1,nr
          do imu=1,node_nmu(ir)
              if((node_intensity(imu,ir).lt.0.d0).or.
     %             (node_intensity(-imu,ir).lt.0.d0)) then
                  stop 45611
              endif
          enddo
          if((node_iray(0,ir).gt.0).and.
     %         (node_intensity(0,ir).lt.0.d0)) then
              stop 45612
          endif
          do imu=node_nmu(ir)+1,FRSIZE_RAYS
              if((node_intensity(imu,ir).lt.1.d98).or.
     %             (node_intensity(-imu,ir).lt.1.d98)) then
                  stop 45613
              endif
          enddo              
      enddo
c
      end


c     --------------------------------------------------------------
c                  PERFORM THE INTEGRAL ALONG THE RAY
c     --------------------------------------------------------------
      subroutine rayint(nr,r,srcarr,alparr,starint,iray,intens)
      implicit none
c
      integer iray,nr
      doubleprecision intens
      doubleprecision srcarr(0:FRSIZE_R)
      doubleprecision alparr(0:FRSIZE_R)
      doubleprecision r(FRSIZE_R),starint
c
#include "common_transfer.h"
c
      integer is,iorder
      doubleprecision ds,ds2,alpprv,srcprv,alpcur,srccur
      doubleprecision alpfut,srcfut,qdr_src_2,qdr_src_3
      doubleprecision dr,dr1
c
      alpprv = alparr(rays_ir(1,iray))
      srcprv = srcarr(rays_ir(1,iray))
      rays_intensity(1,iray)=intens
      if(iter_qdr.eq.1) then
c
c         First order quadrature
c
          do is=2,rays_ns(iray)
              ds = rays_s(is,iray)-rays_s(is-1,iray)
              if(rays_ir(is,iray).gt.0) then
                  call find_srcalp(nr,r,srcarr,alparr,
     %                 rays_ir(is,iray),
     %                 rays_dr(is,iray),srccur,alpcur)
                  if(((rays_ir(is-1,iray).gt.1).and.
     %                 (rays_ir(is,iray).ge.1)).or.
     %                 ((rays_ir(is-1,iray).ge.1).and.
     %                 (rays_ir(is,iray).gt.1)).or.
     %                 ((rays_ir(is-1,iray).eq.1).and.
     %                 (rays_ir(is,iray).eq.1).and.
     %                 ((rays_dr(is,iray).gt.0.d0).or.
     %                 (rays_dr(is-1,iray).gt.0.d0)))) then
                      intens = qdr_src_2(intens,srcprv,alpprv,
     %                     srccur,alpcur,ds)
                  endif
              else
                  intens = starint
              endif
              rays_intensity(is,iray) = intens
              alpprv=alpcur
              srcprv=srccur
          enddo
      elseif(iter_qdr.eq.2) then
c
c         Second order quadrature
c
          if(rays_ir(2,iray).gt.0) then
              call find_srcalp(nr,r,srcarr,alparr,
     %             rays_ir(2,iray),
     %             rays_dr(2,iray),srccur,alpcur)
          else
              stop 77721
          endif
          do is=2,rays_ns(iray)-1
              ds  = rays_s(is,iray)-rays_s(is-1,iray)
              ds2 = rays_s(is+1,iray)-rays_s(is,iray)
              if(rays_ir(is,iray).gt.0) then
                  if(rays_ir(is+1,iray).gt.0) then
                      call find_srcalp(nr,r,
     %                     srcarr,alparr,
     %                     rays_ir(is+1,iray),
     %                     rays_dr(is+1,iray),srcfut,alpfut)
                      iorder=2
                  else
                      iorder=1
                  endif
                  if(((rays_ir(is-1,iray).gt.1).and.
     %                 (rays_ir(is,iray).ge.1)).or.
     %                 ((rays_ir(is-1,iray).ge.1).and.
     %                 (rays_ir(is,iray).gt.1)).or.
     %                 ((rays_ir(is-1,iray).eq.1).and.
     %                 (rays_ir(is,iray).eq.1).and.
     %                 ((rays_dr(is,iray).gt.0.d0).or.
     %                 (rays_dr(is-1,iray).gt.0.d0)))) then
                      if(iorder.eq.2) then
                          intens = qdr_src_3(intens,srcprv,alpprv,
     %                         srccur,alpcur,srcfut,alpfut,
     %                         ds,ds2)
                      else
                          intens = qdr_src_2(intens,srcprv,alpprv,
     %                         srccur,alpcur,ds)
                      endif
                  endif
              else
                  intens = starint
              endif
              rays_intensity(is,iray) = intens
              alpprv=alpcur
              srcprv=srccur
              alpcur=alpfut
              srccur=srcfut
          enddo
          is = rays_ns(iray)
          ds = rays_s(is,iray)-rays_s(is-1,iray)
          intens = qdr_src_2(intens,srcprv,alpprv,
     %         srccur,alpcur,ds)
          rays_intensity(is,iray) = intens
      else
          stop 435
      endif
c
      end



c     --------------------------------------------------------------
c                     FIND SOURCE AND OPACITY
c
c     This routine returns the src and alp at a point on the grid, 
c     given by the index ir and the relative offset dr. It will 
c     automatically choose the right interpolation style.
c     --------------------------------------------------------------
      subroutine find_srcalp(nr,rr,srcarr,alparr,ir,dr,src,alp)
      implicit none
c
      integer ir,nr
      doubleprecision dr,src,alp,rr(FRSIZE_R)
      doubleprecision srcarr(0:FRSIZE_R)
      doubleprecision alparr(0:FRSIZE_R)
c
#include "common_transfer.h"
c
      doubleprecision dr1,xx(4),yy(4),zz(4),r
      integer ix
c
      if((ir.lt.1).or.(ir.gt.nr).or.
     %     ((ir.eq.1).and.(dr.lt.0.d0)).or.
     %     ((ir.eq.nr).and.(dr.gt.0.d0))) then
          write(*,*) ir,dr
          stop 82955
      endif
      if(dr.eq.0.d0) then
          alp = alparr(ir)
          src = srcarr(ir)
      else
          if((iter_qdr.eq.1).or.(ir.eq.1).or.(ir.eq.nr-1)) then
c
c             Linear interpolation
c
              dr1 = 1.d0-dr
              alp = dr1*alparr(ir)+dr*alparr(ir+1)
              src = dr1*srcarr(ir)+dr*srcarr(ir+1)
          else
c
c             Cubic interpolation
c
              r = rr(ir) + dr*(rr(ir+1)-rr(ir))
              do ix=1,4
                  xx(ix) = rr(ir+ix-2)
                  yy(ix) = srcarr(ir+ix-2)
                  zz(ix) = alparr(ir+ix-2)
              enddo
              call polint4(xx,yy,r,src)
              call polint4(xx,zz,r,alp)
          endif
      endif
c
      end


c     --------------------------------------------------------------
c      INTEGRATE THE FORMAL SOLUTION ALONG ONE ELEMENT OF SHORT CHAR
c
c     This function computes the formal integral over a segment of
c     a local characteristic. This will be done first order here.
c
c     This linear integration is not suitable for optically thick
c     systems. In that case the 3-point quadratic integrator,
c     qdr_src_3() should be used. 
c     --------------------------------------------------------------
      function qdr_src_2(int,js1,alp1,js2,alp2,ds)
      implicit none
      doubleprecision qdr_src_2
      doubleprecision int,js1,alp1,js2,alp2,ds
c
      doubleprecision e0,e1,a,b,dtau1,xp,src1,src2,theomax
c
      dtau1     = 0.5d0 * ( alp1 + alp2 ) * ds
      theomax   = 0.5d0 * ( js1 + js2 ) * ds
      if(dtau1.gt.1.d-6) then
          xp        = exp(-dtau1)
          e0        = 1.d0 - xp
          e1        = dtau1 - e0
          b         = e1 / dtau1
          a         = e0 - b
      else
          a         = 0.5d0 * dtau1
          b         = 0.5d0 * dtau1
          xp        = 1.d0 - dtau1 
      endif
      if(alp1.gt.0.d0) then
          src1 = js1 / alp1
      elseif(alp2.gt.0.d0) then
          src1 = js2 / alp2
      else
          src1 = 0.d0
      endif
      if(alp2.gt.0.d0) then
          src2 = js2 / alp2
      elseif(alp1.gt.0.d0) then
          src2 = js1 / alp1
      else
          src2 = 0.d0
      endif
      if(dtau1.gt.1e-9) then
         qdr_src_2 = a * src1 + b * src2
      else
         qdr_src_2 = 0.5d0 * ( js1 + js2 ) * ds
      endif
c
c     Verify if serious overshoot occur (say: factor 2 or more)
c
c#     ifndef SRCQDR_NO_ERRCHECK
c      if(qdr_src_2.gt.2.d0*theomax) then
c          transerr_overshoot =  transerr_overshoot + 1
c      endif
c#     endif
c
c     NEW: Bracket the total integrated emission between the
c          theoretical maximum and 0.
c
#     ifndef SRCQDR_NO_BRACKET
      qdr_src_2 = min(qdr_src_2,theomax)
#     endif
c
      qdr_src_2 = int * xp + qdr_src_2
      return
c
      write(*,*) 'Ehh, something wrong in dqr_src_2'
      stop 13
      return
      end




c     --------------------------------------------------------------
c      INTEGRATE THE FORMAL SOLUTION ALONG ONE ELEMENT OF SHORT CHAR
c
c     This function computes the formal integral over a segment of
c     a local characteristic. This will be done second order here.
c     This is necessary for optically thick systems, but hardy so
c     for optically thin systems.
c
c     WARNING: Be careful for optical depths close to 1e-6, because
c              this regime might be uncertain.....
c
c     --------------------------------------------------------------
      function qdr_src_3(int,js1,alp1,js2,alp2,js3,alp3,ds1,ds2)
      implicit none
      doubleprecision qdr_src_3
      doubleprecision int,js1,alp1,js2,alp2,js3,alp3,ds1,ds2
c
      doubleprecision e0,e1,e2,a,b,c,dtau1,dtau2,dtsum,xp
      doubleprecision src1,src2,src3,theomax,theomin
c
      dtau1     = 0.5d0 * ( alp1 + alp2 ) * ds1
      dtau2     = 0.5d0 * ( alp2 + alp3 ) * ds2
      dtsum     = dtau1 + dtau2
c
c     If the safety checks are on, then check something out
c
#     ifdef SAFETY_CHECKS_ACTIVE
      if(dtau1.gt.selfch_taumax) then
          selfch_taumax = dtau1
      endif
#     endif
c
c     Now compute the theoretical maximum and minimum of the
c     quadrature. These are used to bracket the Olson-Kunasz formula
c     which may sometimes go wrong.
c     I did this previously in an erronous way.  But since I've used
c     that (wrong) formula for a long time (until 6 July 1999) I rather
c     not remove it completely.  I still keep it an option, just in case
c     old runs may be incompatible with this new thing.
c
#     ifndef SRCQDR_OLDTHEOMAX
c
c     The new one
c
      theomax   = max(js1,js2) * ds1
      theomin   = 0.d0        ! min(js1,js2) * ds1 * exp(-dtau1)
#     else
c
c     The old one ( < 6 july 1999 )
c
      theomax   = 0.5d0 * ( js1 + js2 ) * ds1
      theomin   = 0.d0
#     endif
c
      if(dtau1.gt.1.d-6) then
          xp        = exp(-dtau1)
          e0        = 1.d0 - xp
          e1        = dtau1 - e0
          e2        = dtau1 * dtau1 - 2.d0 * e1
          if(dtau2.gt.1.d-6) then
              b  = ( dtsum * e1 - e2 ) / ( dtau1 * dtau2 )
              a  = e0 + ( e2 - ( dtsum + dtau1 ) * e1 ) 
     %                    / ( dtau1 * dtsum )
              c  = ( e2 - dtau1 * e1 ) / ( dtau2 * dtsum )
          else  
              b  = e1 / dtau1
              a  = e0 - b
              c  = 0.d0
          endif
          if(alp1.gt.0.d0) then
             src1 = js1 / alp1
          elseif(alp2.gt.0.d0) then
             src1 = js2 / alp2
          else
             src1 = 0.d0
          endif
          if(alp2.gt.0.d0) then
             src2 = js2 / alp2
          elseif(alp1.gt.0.d0) then
             src2 = js1 / alp1
          elseif(alp3.gt.0.d0) then
             src2 = js3 / alp3
          else
             src2 = 0.d0
          endif
          if(alp3.gt.0.d0) then
             src3 = js3 / alp3
          elseif(alp2.gt.0.d0) then
             src3 = js2 / alp2
          else
             src3 = 0.d0
          endif
          qdr_src_3 = a * src1 + b * src2 + c * src3
c
c         Verify if serious overshoot occur (say: factor 2 or more)
c
c#         ifndef SRCQDR_NO_ERRCHECK
c          if(qdr_src_3.gt.2.d0*theomax) then
c              transerr_overshoot =  transerr_overshoot + 1
c          endif
c#         endif
c
c         NEW: Bracket the total integrated emission between the
c              theoretical maximum and minimum.
c
#         ifndef SRCQDR_NO_BRACKET
          qdr_src_3 = min(qdr_src_3,theomax)
          qdr_src_3 = max(qdr_src_3,theomin)
#         endif
c
          qdr_src_3 = int * xp + qdr_src_3
c
      else
          qdr_src_3 = int + 0.5d0 * ( js1 + js2 ) * ds1
      endif
      return
c
      write(*,*) 'Ehh, something wrong in dqr_src_3'
      stop 13
      return
      end


c     --------------------------------------------------------------
c             POLINOMIAL INTERPOLATION FOR FOUR POINTS
c     --------------------------------------------------------------
      subroutine polint4(xa,ya,x,y)
      doubleprecision x,y,xa(4),ya(4)
      
      y = (((x-xa(2))*(x-xa(3))*(x-xa(4)))
     % /((xa(1)-xa(2))*(xa(1)-xa(3))*(xa(1)-xa(4))))*ya(1)  
 
     % + (((x-xa(1))*(x-xa(3))*(x-xa(4)))
     % /((xa(2)-xa(1))*(xa(2)-xa(3))*(xa(2)-xa(4))))*ya(2) 

     % + (((x-xa(1))*(x-xa(2))*(x-xa(4)))
     % /((xa(3)-xa(1))*(xa(3)-xa(2))*(xa(3)-xa(4))))*ya(3) 

     % + (((x-xa(1))*(x-xa(2))*(x-xa(3)))
     % /((xa(4)-xa(1))*(xa(4)-xa(2))*(xa(4)-xa(3))))*ya(4)

      return
      end
