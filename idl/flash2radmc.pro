;======================================================================
;                  FLASH-CODE TO RADMC-3D INTERFACE
;
; The FLASH code is an astrophysical (magneto-)hydrodynamics code:
;
;    http://flash.uchicago.edu/
;
; The subroutine flash2radmc takes FLASH output and converts it into
; input for RADMC-3D. But please keep in mind that this only produces
; the dust density distribution, assuming "gastodust" as the gas-to-
; dust ratio. The positions of stars and their properties, as well as
; opacities and such have to be separately made by the user. 
;
; NOTE: This routine requires the read_amr_hdf5.pro routine from the
;       flash code distribution.
;======================================================================
@read_amr_hdf5.pro
;----------------------------------------------------------------------
;    PRODUCE RADMC-3D INPUT DATA FILES FROM FLASH OUTPUT DATA FILE
;                 Bernd Voelkl, February 4, 2009
;             Cornelis Dullemond, February 9, 2009
;                         
;----------------------------------------------------------------------
pro flash2radmc,filename,gastodust=gastodust,writetemp=writetemp,form=form,$
                center=center
  ;;
  ;; Some constants
  ;;
  GG  = 6.672d-8                ; Gravitational constant
  mp  = 1.6726d-24              ; Mass of proton          [g]
  kk  = 1.3807d-16              ; Bolzmann's constant     [erg/K]
  cc  = 2.9979245800000d10      ; Light speed             [cm/s]
  ss  = 5.6703d-5               ; Stefan-Boltzmann const  [erg/cm^2/K^4/s]
  LS  = 3.8525d33               ; Solar luminosity        [erg/s]
  RS  = 6.96d10                 ; Solar radius            [cm]
  MS  = 1.99d33                 ; Solar mass              [g]
  TS  = 5.78d3                  ; Solar temperature       [K]
  AU  = 1.496d13                ; Astronomical Unit       [cm]
  pc  = 3.08572d18              ; Parsec                  [cm]
  ;;
  ;; Some defaults
  ;;
  if not keyword_set(gastodust) then gastodust=100.
  ;;
  ;; Read the HDF5 data file from FLASH
  ;;
  read_amr_hdf5, filename, TREE=tree, DATA=data, $ 
        PARAMETERS=params, STORED_VARS=varnames
  ;;
  ;; Find total number of blocks (branches and leafs)
  ;;
  nblocks    = n_elements(tree)
  ;;
  ;; Find the leaf blocks
  ;;
  leaf_nodes = where(tree[*].nodetype EQ 1)
  nleafs     = n_elements(leaf_nodes)
  ;;
  ;; Keys for data request
  ;; Standard FLASH keys, there are more
  ;;
  idens = (where( varnames eq 'dens' ))[0]
  itemp = (where( varnames eq 'pres' ))[0]
  ipres = (where( varnames eq 'pres' ))[0]
  ivelx = (where( varnames eq 'velx' ))[0]
  ively = (where( varnames eq 'vely' ))[0]
  ivelz = (where( varnames eq 'velz' ))[0]
  ;;
  ;; Get min, max refinement levels
  ;;
  lrefine_max = max(tree[leaf_nodes].lrefine)
  lrefine_min = min(tree[leaf_nodes].lrefine)
  ;;
  ;; Find the overall bounding box
  ;;
  xmin  = min(tree[*].bndbox[*,0])
  xmax  = max(tree[*].bndbox[*,0])
  ymin  = min(tree[*].bndbox[*,1])
  ymax  = max(tree[*].bndbox[*,1])
  zmin  = min(tree[*].bndbox[*,2])
  zmax  = max(tree[*].bndbox[*,2])
  ;;
  ;; How much do we have to shift to center
  ;;
  if keyword_set(center) then begin
     centerpos = 0.5 * [xmax+xmin,ymax+ymin,zmax+zmin] 
  endif else begin
     centerpos = [0.,0.,0.]
  endelse
  ;;
  ;; Make a memory of parent blocks
  ;;
  ichildblock  = intarr(3,lrefine_max+1)
  ilevel       = 0
  itheleafs    = intarr(nleafs)
  ;;
  ;; Open grid file for RADMC-3D
  ;;
  print,'Writing amr_grid.inp...'
  openw,1,'amr_grid.inp'
  printf,1,1               ; Format file
  printf,1,1               ; AMR Style
  printf,1,3               ; Coord system
  printf,1,0               ; No grid info
  printf,1,1,1,1           ; Include all dimensions
  printf,1,1,1,1           ; Base grid has 1x1x1 dims = 1 base cell
  printf,1,lrefine_max+3L,$                         ; Max refinement level, counted from 0
           nleafs*8L*8L*8L,$                        ; Nr of leafs for RADMC-3D
           nblocks*(1L+2L*2L*2L+4L*4L*4L+8L*8L*8L)  ; Nr of branches for RADMC-3D
  printf,1,xmin-centerpos[0],xmax-centerpos[0]       ; RADMC-3D Base grid is 1 cell
  printf,1,ymin-centerpos[1],ymax-centerpos[1]
  printf,1,zmin-centerpos[2],zmax-centerpos[2]
  ;;
  ;; First level is the base block
  ;;
  iiblock = intarr(lrefine_max+1)
  iiblock[0] = (where(tree[*].lrefine eq 1))[0]+1
  ;;
  ;; While in RADMC-3D we have in fact a base grid, we will follow
  ;; the Flash convention and just have one base cell and simply
  ;; refine it recursively
  ;;
  ileaf = 0
  while ilevel ge 0 do begin
     ;;
     ;; Find the current block from the parent
     ;;
     if ilevel gt 0 then begin 
        ich = ichildblock[2,ilevel-1]*4+ichildblock[1,ilevel-1]*2+ichildblock[0,ilevel-1]
        iiblock[ilevel] = tree[iiblock[ilevel-1]-1].gid[7+ich]
     endif
     ;;
     ;; Find out if the current block is refined or not
     ;;
     ;;;print,iiblock[ilevel],tree[iiblock[ilevel]-1].lrefine-1
     if tree[iiblock[ilevel]-1].nodetype eq 1 then begin
        ;;
        ;; Leaf block of 8x8x8 cells. 
        ;;
        itheleafs[ileaf] = iiblock[ilevel]
        ileaf = ileaf+1
        ;;
        ;; For RADMC-3D this is a further subtree of 3 levels.
        ;;
        printf,1,'1'
        for iz1=0,1 do begin
        for iy1=0,1 do begin
        for ix1=0,1 do begin
           printf,1,'1'
           for iz2=0,1 do begin
           for iy2=0,1 do begin
           for ix2=0,1 do begin
              printf,1,'1'
              for iz3=0,1 do begin
              for iy3=0,1 do begin
              for ix3=0,1 do begin
                 ix = ix1*4+ix2*2+ix3
                 iy = iy1*4+iy2*2+iy3
                 iz = iz1*4+iz2*2+iz3
                 printf,1,'0'
              endfor
              endfor
              endfor
           endfor
           endfor
           endfor
        endfor
        endfor
        endfor
        ;;
        ;; Update parent's child list
        ;;
        ret = 1
        while ilevel gt 0 and ret eq 1 do begin
           ret = 0
           ichildblock[0,ilevel-1] = ichildblock[0,ilevel-1] + 1
           if ichildblock[0,ilevel-1] eq 2 then begin
              ichildblock[0,ilevel-1] = 0
              ichildblock[1,ilevel-1] = ichildblock[1,ilevel-1] + 1
              if ichildblock[1,ilevel-1] eq 2 then begin
                 ichildblock[1,ilevel-1] = 0
                 ichildblock[2,ilevel-1] = ichildblock[2,ilevel-1] + 1
                 if ichildblock[2,ilevel-1] eq 2 then begin
                    ichildblock[2,ilevel-1] = 0
                    ilevel=ilevel-1
                    ret = 1
                 endif
              endif
           endif
        endwhile
     endif else begin
        ;;
        ;; Branch, so refine
        ;;
        printf,1,'1'
        ichildblock[*,ilevel] = 0
        ilevel = ilevel + 1
     endelse
     if ilevel eq 0 then ilevel = -1
  endwhile
  ;;
  ;; Now close this file amr_grid.inp
  ;;
  close,1
   ;;
  ;; Check
  ;;
  if ileaf ne nleafs then begin
     print,'ERROR: Nr of leafs wrong'
     stop
  endif
  ;;
  ;; Now we simply dump the density, divided by the gas-to-dust ratio
  ;;
  if keyword_set(form) then begin
     print,'Writing dust_density.inp...'
     openw,1,'dust_density.inp'
     printf,1,1                      ; Format number
     printf,1,nleafs*8L*8L*8L        ; Nr of leafs
     printf,1,1                      ; Nr of dust species
     for ileaf=0,nleafs-1 do begin
        for iz1=0,1 do begin
        for iy1=0,1 do begin
        for ix1=0,1 do begin
           for iz2=0,1 do begin
           for iy2=0,1 do begin
           for ix2=0,1 do begin
              for iz3=0,1 do begin
              for iy3=0,1 do begin
              for ix3=0,1 do begin
                 ix = ix1*4+ix2*2+ix3
                 iy = iy1*4+iy2*2+iy3
                 iz = iz1*4+iz2*2+iz3
                 printf,1,data[idens,itheleafs[ileaf]-1,ix,iy,iz]/gastodust
              endfor
              endfor
              endfor
           endfor
           endfor
           endfor
        endfor
        endfor
        endfor
     endfor
     close,1
  endif else begin
     dumdat=dblarr(8*8*8)
     print,'Writing dust_density.uinp...'
     openw,1,'dust_density.uinp',/f77_unformatted
     writeu,1,1LL,8LL*8LL*8LL*8LL       ; Format number , record length
     writeu,1,nleafs*8LL*8LL*8LL,1LL    ; Nr of leafs , Nr of dust species
     for ileaf=0,nleafs-1 do begin
        ii=0
        for iz1=0,1 do begin
        for iy1=0,1 do begin
        for ix1=0,1 do begin
           for iz2=0,1 do begin
           for iy2=0,1 do begin
           for ix2=0,1 do begin
              for iz3=0,1 do begin
              for iy3=0,1 do begin
              for ix3=0,1 do begin
                 ix = ix1*4+ix2*2+ix3
                 iy = iy1*4+iy2*2+iy3
                 iz = iz1*4+iz2*2+iz3
                 dumdat[ii] = data[idens,itheleafs[ileaf]-1,ix,iy,iz]/gastodust
                 ii = ii + 1
              endfor
              endfor
              endfor
           endfor
           endfor
           endfor
        endfor
        endfor
        endfor
        writeu,1,dumdat
     endfor
     close,1
  endelse
  ;;
  ;; Now if the user wants this, we can also immmediately dump the 
  ;; temperature as well, so that we do not have to calculate the
  ;; dust temperature with the Monte Carlo method.
  ;;
  if keyword_set(writetemp) then begin
     if keyword_set(form) then begin
        print,'Writing dust_temperature.dat...'
        openw,1,'dust_temperature.dat'
        printf,1,1                 ; Format number
        printf,1,nleafs*8L*8L*8L   ; Nr of leafs
        printf,1,1                 ; Nr of dust species
        for ileaf=0,nleafs-1 do begin
           for iz1=0,1 do begin
           for iy1=0,1 do begin
           for ix1=0,1 do begin
              for iz2=0,1 do begin
              for iy2=0,1 do begin
              for ix2=0,1 do begin
                 for iz3=0,1 do begin
                 for iy3=0,1 do begin
                 for ix3=0,1 do begin
                    ix = ix1*4+ix2*2+ix3
                    iy = iy1*4+iy2*2+iy3
                    iz = iz1*4+iz2*2+iz3
                    temp = 2.3*mp*(data[ipres,itheleafs[ileaf]-1,ix,iy,iz]/$
                           data[idens,itheleafs[ileaf]-1,ix,iy,iz])/kk
                    printf,1,temp
                 endfor
                 endfor
                 endfor
              endfor
              endfor
              endfor
           endfor
           endfor
           endfor
        endfor
        close,1
     endif else begin
        dumdat=dblarr(8*8*8)
        print,'Writing dust_temperature.udat...'
        openw,1,'dust_temperature.udat',/f77_unformatted
        writeu,1,1LL,8LL*8LL*8LL*8LL      ; Format number , record length
        writeu,1,nleafs*8LL*8LL*8LL,1LL   ; Nr of leafs , Nr of dust species
        for ileaf=0,nleafs-1 do begin
           ii = 0
           for iz1=0,1 do begin
           for iy1=0,1 do begin
           for ix1=0,1 do begin
              for iz2=0,1 do begin
              for iy2=0,1 do begin
              for ix2=0,1 do begin
                 for iz3=0,1 do begin
                 for iy3=0,1 do begin
                 for ix3=0,1 do begin
                    ix = ix1*4+ix2*2+ix3
                    iy = iy1*4+iy2*2+iy3
                    iz = iz1*4+iz2*2+iz3
                    dumdat[ii] = 2.3*mp*(data[ipres,itheleafs[ileaf]-1,ix,iy,iz]/$
                           data[idens,itheleafs[ileaf]-1,ix,iy,iz])/kk
                    ii = ii + 1
                 endfor
                 endfor
                 endfor
              endfor
              endfor
              endfor
           endfor
           endfor
           endfor
           writeu,1,dumdat
        endfor
        close,1
     endelse
  endif
end
