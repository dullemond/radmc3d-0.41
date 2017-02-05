pro amr2radmc,nout=nout,gastodust=gastodust
;===================================================================
; This routine convert RAMSES AMR structure to the RADMC-3D format
; Inoutput, this returns the files amr_grid.inp, dust_density.inp
; and dust_temperature.dat 
;
; Written by:	Benoit Commercon, 07/2011.
;===================================================================

strout='_00'+string(nout,format='(i3)')
if(nout lt 100) then strout='_000' +string(nout,format='(i2)')
if(nout lt 10 ) then strout='_0000'+string(nout,format='(i1)')

  
;; Some constants
;;
GG  = 6.672d-8                  ; Gravitational constant
mp  = 1.6726d-24                ; Mass of proton          [g]
kk  = 1.3807d-16                ; Bolzmann's constant     [erg/K]
cc  = 2.9979245800000d10        ; Light speed             [cm/s]
ss  = 5.6703d-5                 ; Stefan-Boltzmann const  [erg/cm^2/K^4/s]
LS  = 3.8525d33                 ; Solar luminosity        [erg/s]
RS  = 6.96d10                   ; Solar radius            [cm]
MS  = 1.99d33                   ; Solar mass              [g]
TS  = 5.78d3                    ; Solar temperature       [K]
AU  = 1.496d13                  ; Astronomical Unit       [cm]
pc  = 3.08572d18                ; Parsec                  [cm]
mu = 2.375                      ; Mean molecular weigth
;;
;; Some defaults
;;
if not keyword_set(gastodust) then gastodust=100.

;Read grid and hydro
rd_all,grid,hydro,nout=nout
rd_info,info,nout=nout
amr2cell,grid,hydro,ccell
help,ccell,/str

unit_tp = info.unit_l*info.unit_l*mu*mp/(kk*info.unit_t*info.unit_t)

nrcells = ccell.n
ngridtot = ((total(grid.ngrid)*8L+1L))
lmax = grid.nlevelmax

lun_d   = 2
lun_tp  = 3
lun_amr = 4

openw,lun_amr,/get_lun,'amr_grid'+strout+'.inp'
printf,lun_amr,'1'         ; Format
printf,lun_amr,'1'         ; AMR style!
printf,lun_amr,'1'         ; Coordinate system: cartesioan (<100)
printf,lun_amr,'0'        ; Grid info
printf,lun_amr,'1',' 1',' 1' ; Three dimensions are active
printf,lun_amr,'1',' 1',' 1' ; Grids at level 0 = 1 base cell
printf,lun_amr,lmax, $     ; Max refinement level, counted form 0
       nrcells,$           ; Number of leaf cells for RADMC
       ngridtot,FORMAT = '(I3,1X,I9,1X,I9)' ; Number of branches for RADMC
               
xc=grid.boxlen*info.unit_l/2.0d0

zero=0.0
printf,lun_amr,zero-xc,grid.boxlen*info.unit_l-xc ; Grid dimension x
printf,lun_amr,zero-xc,grid.boxlen*info.unit_l-xc ; Grid dimension y
printf,lun_amr,zero-xc,grid.boxlen*info.unit_l-xc ; Grid dimension z



openw,lun_d,/get_lun,'dust_density'+strout+'.inp'
printf,lun_d,'1'           ; Format
printf,lun_d,nrcells       ; Number of leaf cells
printf,lun_d,'1'           ; Number of dust specy

openw,lun_tp,/get_lun,'dust_temperature'+strout+'.dat'
printf,lun_tp,'1'          ; Format
printf,lun_tp,nrcells      ; Number of leaf cells
printf,lun_tp,'1'          ; Number of dust specy


ison=0 ; son index in the oct
ilevel=0
cpu_1st_arr=where((grid.ngrid(0,*)) gt 0)


cpu_1st = cpu_1st_arr(0)

; Look for the first cell at level 0

cpu_father=cpu_1st              ; 
cell=*(grid.level(0,cpu_1st))
son_cpu = cell.cpumap(0,ison)-1                      ; Get ison cpu
ind_son = (*grid.level(0,cpu_1st)).son(0,ison)       ; Get ison index
father_cell = *(grid.level(0,cpu_1st))
father_cell.indgrid = 1

son_coarse=0

cell=*(grid.level(ilevel+1,son_cpu))

cell.cpu_father(*)=cpu_father      
father_cell.cpu(*)  = cpu_father

(*grid.level(0,cpu_1st)).ison(*)=ison+1
(*grid.level(0,cpu_1st)).cpu(*)=cpu_1st

coarse_cell=father_cell
printf,lun_amr,'1'
printf,lun_amr,'1'

gson = where(cell.indgrid eq ind_son(0)) ; get son index in the arrays of 'cell' (0=1st element)
cell.cpu(gson(0)) = father_cell.cpumap(0,0)-1
(*grid.level(ilevel+1,son_cpu)).cpu(gson(0)) = father_cell.cpumap(0,0)
(*grid.level(ilevel+1,son_cpu)).cpu_father(gson(0))=cpu_father
(*grid.level(ilevel+1,son_cpu)).ison_father(gson(0))=ison
(*grid.level(ilevel+1,son_cpu)).indgrid_father(gson(0))=father_cell.indgrid 

cell=*(grid.level(ilevel+1,son_cpu))
gfather=0
father_cell = *(grid.level(0,cpu_1st))

;=================================================
; Store data for the first oct
;=================================================
for ison1=0,7 do begin         
   ind_son2 = father_cell.son(0,ison1) 
   if(ind_son2 gt 0) then begin

      cell_cpu2=father_cell.cpumap(0,ison1)-1
      cell2=*(grid.level(ilevel+1,cell_cpu2))

      gson2 = where(cell2.indgrid eq ind_son2(0))
      (*grid.level(ilevel+1,cell_cpu2(0))).cpu(gson2(0)) = father_cell.cpumap(0,ison1)
      (*grid.level(ilevel+1,cell_cpu2(0))).ison(gson2(0))=0
      if(ison1 lt 7)then (*grid.level(ilevel+1,cell_cpu2(0))).next(gson2(0))=father_cell.son(0,ison1+1)
      if(ison1 lt 7)then (*grid.level(ilevel+1,cell_cpu2(0))).cpu_next(gson2(0))=father_cell.cpumap(0,ison1+1)
      (*grid.level(ilevel+1,cell_cpu2(0))).cpu_father(gson2(0))=cpu_father
      (*grid.level(ilevel+1,cell_cpu2(0))).ison_father(gson2(0))=ison1
      (*grid.level(ilevel+1,cell_cpu2(0))).indgrid_father(gson2(0))=father_cell.indgrid(gfather(0))
   endif
endfor

;=================================================
; Ecrire ici au format radmc dans amr_grid.inp
;=================================================

ilevel = ilevel + 1
nbcell=0

ison=0

while(ilevel ge 0) do begin
  BIG_LOOP: nb=0
   
   if(cell.son(gson(0),ison) eq 0) then begin ; This is a leaf cell
      
      if(cell.ison(gson(0)) gt 8)then begin
         print,'PROBLEM ison > 8!!!, line 496'
         read,ooo
         stop
      endif
      
;=================================================
; Ecrire ici au format radmc dans amr_grid.inp
;=================================================
      
      nbcell=nbcell+1L
      
      icpu = cell.cpumap(gson(0),ison)-1
      
      dens  = (*hydro.levelh(ilevel,son_cpu)).u(gson(0),ison, 0)
      press = (*hydro.levelh(ilevel,son_cpu)).u(gson(0),ison,10)
      (*hydro.levelh(ilevel,son_cpu)).u(gson(0),ison, 0) = -1.0
      (*hydro.levelh(ilevel,son_cpu)).u(gson(0),ison, 10) = -1.0

      if(dens lt 0 or press lt 0)then stop ; Cell already checked... There is a bug...

      tp = press/dens
      
      dens = dens * info.unit_d
      tp   = tp   * unit_tp

      printf,lun_amr,'0';,nbcell
      printf,lun_d,dens/gastodust,FORMAT = '(e12.5)'
      printf,lun_tp,tp,FORMAT = '(e12.5)'

      LOOP1:  loop1=1
      
      if(ison lt 7) then begin ; Go to the neighboring cell, correponding to the next cell in the father oct         
         (*grid.level(ilevel,son_cpu)).ison(gson(0))=ison+1
         ison=ison+1
      endif else begin ; Go to the neighboring cell of father, correponding to the next cell in the opa oct ;else l. 536
         ind_next = cell.next(gson(0))
         cpu_next = cell.cpu_next(gson(0))-1
         if(ind_next le 0)then begin ; No next cell in the oct. Go to opa oct.
            ilevel=ilevel-1
            if(ilevel le 1) then begin ; if l693
               
               ;; LOOP_COARSE:            
               father_Cell=coarse_Cell
               cell=father_Cell
               son_coarse=son_coarse +1 
               ison =         son_coarse                  
               gfather = 0 
               gson=0
               if(ison gt 7) then GOTO,FIN
               
               ind_next = cell.son(gson(0),ison)
               cpu_next = cell.cpumap(gson(0),ison)-1
               cell = *(grid.level(ilevel,cpu_next))
               gson=where(cell.indgrid eq ind_next)
               cpu_father = cpu_next
               ilevel=1
               
               printf,lun_amr,'1' ; Format

               GOTO, LOOP2
            endif
            ind_father = cell.indgrid_father(gson(0))
            cpu_father = cell.cpu_father(gson(0))
            cell = *(grid.level(ilevel,cpu_father))
            gson = where(father_cell.indgrid eq ind_father)
            
            son_cpu=cpu_father
            ind_father = cell.indgrid_father(gson(0))
            cpu_father = cell.cpu_father(gson(0))
 
            father_cell = *(grid.level(ilevel-1,cpu_father))
            gfather = where(father_cell.indgrid eq ind_father)
            ison_father = cell.ison_father(gson(0))
            (*grid.level(ilevel-1,cpu_father)).ison(gfather(0))=ison_father
            father_cell = *(grid.level(ilevel-1,cpu_father))    
            
            ison = cell.ison(gson(0))

            GOTO, LOOP1

         endif else begin
            cell = *(grid.level(ilevel,cpu_next))
            gson = where(cell.indgrid eq ind_next)
            ison = (*grid.level(ilevel,cpu_next)).ison(gson(0))

            ind_father = cell.indgrid_father(gson(0))
            cpu_father = cell.cpu_father(gson(0))
            father_cell = *(grid.level(ilevel-1,cpu_father))
            gfather = where(father_cell.indgrid eq ind_father)

            ison_father = cell.ison_father(gson(0))
            (*grid.level(ilevel-1,cpu_father)).ison(gfather(0))=ison_father
            father_cell = *(grid.level(ilevel-1,cpu_father))     
            son_cpu=cpu_next
            
            printf,lun_amr,'1'  ; Format
 
        endelse


      endelse

    endif else begin             ; The cell is refined
       cpu_father=cell.cpu(gson(0))-1 ; 
       
       LOOP2:      printf,lun_amr,'1'
       
       father_cell=cell
       gfather=gson
       
       ilevel=ilevel+1
       
       ison=(*grid.level(ilevel-1,cpu_father)).ison(gfather(0))
       
       son_cpu = cell.cpumap(gson(0),ison)-1 ; Get ison cpu
       cell=*(grid.level(ilevel,son_cpu))
       
       gson = where(((cell.father mod grid.ngridmax)-grid.ncoarse)eq father_cell.indgrid(gfather(0)) and fix(cell.father/grid.ngridmax) eq ison,na) 
       
       if(son_cpu ne cpu_father or na eq 0)then begin      
          dx = 1./2.^(ilevel+1)
          ix=-1.0d0 
          iy=-1.0d0 
          iz=-1.0d0 
          
          if((ison mod 2) eq 1) then ix=1.0d0   
          if((ison mod 4) ge 2) then iy=1.0d0   
          if((ison ge 4)      ) then iz=1.0d0   
          
          xson = father_cell.xg(gfather(0),0)+ix*dx
          yson = father_cell.xg(gfather(0),1)+iy*dx
          zson = father_cell.xg(gfather(0),2)+iz*dx
          
          aa=where(cell.xg(*,0) eq xson and cell.xg(*,1) eq yson and cell.xg(*,2) eq zson,na)
          
          gson = aa
       endif 
       
       (*grid.level(ilevel,son_cpu)).cpu(gson(0))  = father_cell.cpumap(gfather(0),ison)
       (*grid.level(ilevel,son_cpu)).ison(gson(0)) = ison 
       if(ison lt 7)then (*grid.level(ilevel,son_cpu)).next(gson(0)) = father_cell.son(gfather(0),ison+1)
       if(ison lt 7)then (*grid.level(ilevel,son_cpu)).cpu_next(gson(0)) = father_cell.cpumap(gfather(0),ison+1)
       (*grid.level(ilevel,son_cpu)).cpu_father(gson(0))     = cpu_father
       (*grid.level(ilevel,son_cpu)).ison_father(gson(0))    = fix(cell.father(gson(0))/grid.ngridmax)
       (*grid.level(ilevel,son_cpu)).indgrid_father(gson(0)) = father_cell.indgrid(gfather(0))
       cell=*(grid.level(ilevel,son_cpu))

       if(son_cpu ne father_cell.cpumap(gfather(0),ison)-1)then stop
       
       
       ind_father = cell.indgrid_father(gson(0))
       cpu_father = cell.cpu_father(gson(0))
       father_cell=*(grid.level(ilevel-1,cpu_father))
       
       GFATHER = WHERE(FATHER_CELL.INDGRID EQ IND_FATHER,NA)
       if(na le 0)then begin    ;cpu_father ne son_cpu)then begin
          xx = cell.xg(*,0)
          yy = cell.xg(*,1)
          zz = cell.xg(*,2)
          
          dx = 1./2.^(ilevel+1)
          ix=1.0d0 
          iy=1.0d0 
          iz=1.0d0 
          
          if((ison mod 2) eq 1) then ix=-1.0d0   
          if((ison mod 4) ge 2) then iy=-1.0d0   
          if((ison ge 4)      ) then iz=-1.0d0   
          
          xfather = cell.xg(gson(0),0)+ix*dx
          yfather = cell.xg(gson(0),1)+iy*dx
          zfather = cell.xg(gson(0),2)+iz*dx
          
          aa=where(father_cell.xg(*,0) eq xfather and father_cell.xg(*,1) eq yfather and father_cell.xg(*,2) eq zfather,na)
          
          gfather = aa 
       endif              
       
       (*grid.level(ilevel-1,cpu_father)).ison(gfather(0)) = fix(cell.father(gson(0))/grid.ngridmax)
       
       ison_opa = (*grid.level(ilevel-1,cpu_father)).ison_father(gfather(0))
       (*grid.level(ilevel-1,cpu_father)).ison_father(gfather(0))    = ison_opa  + 1
       
       father_cell=*(grid.level(ilevel-1,cpu_father))
              
       ind_son    = cell.indgrid(gson(0)) 
       ind_father = cell.indgrid_father(gson(0))


;=================================================
; Store data for the rest of the oct
;=================================================
       ind_son_prev=-1
       for ison1=0,7 do begin         
          ind_son2 = father_cell.son(gfather(0),ison1) 
          if(ind_son2 gt 0) then begin
             
             cell_cpu2=father_cell.cpumap(gfather(0),ison1)-1
             cell2=*(grid.level(ilevel,cell_cpu2))
             
             gson2 = where(cell2.indgrid eq ind_son2(0))
             if(cpu_father ne cell_cpu2)then begin
                xx = cell2.xg(*,0)
                yy = cell2.xg(*,1)
                zz = cell2.xg(*,2)
                
                dx = 1./2.^(ilevel+1)
                ix=-1.0d0 
                iy=-1.0d0 
                iz=-1.0d0 
                
                if((ison1 mod 2) eq 1) then ix=1.0d0   
                if((ison1 mod 4) ge 2) then iy=1.0d0   
                if((ison1 ge 4)      ) then iz=1.0d0   
                
                xfather = father_cell.xg(gfather(0),0)+ix*dx
                yfather = father_cell.xg(gfather(0),1)+iy*dx
                zfather = father_cell.xg(gfather(0),2)+iz*dx
                
                aa=where(cell2.xg(*,0) eq xfather and cell2.xg(*,1) eq yfather and cell2.xg(*,2) eq zfather,na)
                
                gson2 = aa 
                
; Change next index in previous cell
                if(ison1 gt 0 and ind_son_prev gt 0)then (*grid.level(ilevel,cpu_prev(0))).next(gson_prev(0))=cell2.indgrid(gson2(0))                     
             endif              
 
             (*grid.level(ilevel,cell_cpu2(0))).cpu(gson2(0)) = father_cell.cpumap(gfather(0),ison1)
             (*grid.level(ilevel,cell_cpu2(0))).ison(gson2(0))=0
             if(ison1 lt 7)then (*grid.level(ilevel,cell_cpu2(0))).next(gson2(0))=father_cell.son(gfather(0),ison1+1)
             if(ison1 lt 7)then (*grid.level(ilevel,cell_cpu2(0))).cpu_next(gson2(0))=father_cell.cpumap(gfather(0),ison1+1)
             
             
             (*grid.level(ilevel,cell_cpu2(0))).cpu_father(gson2(0))=cpu_father
             (*grid.level(ilevel,cell_cpu2(0))).ison_father(gson2(0))=ison1
             (*grid.level(ilevel,cell_cpu2(0))).indgrid_father(gson2(0))=father_cell.indgrid(gfather(0))
             ind_son2    = cell2.indgrid(gson2(0)) 
             
             cpu_prev = cell_cpu2
             gson_prev = gson2 
          endif
          ind_son_prev = ind_son2
       endfor
       
       ison=0
       cell=*(grid.level(ilevel,son_cpu))

    endelse
        
 endwhile

FIN: close,lun_amr
close,lun_d
close,lun_tp
print,nbcell            
end
;+
; NAME:
;	RD_AMR
;
; PURPOSE:
;	This procedure reads the mesh structure from a RAMSES AMR file for 
;       interfacing with RADMC-3D.
;
; CATEGORY:
;	Input/Output.
;
; CALLING SEQUENCE:
;	RD_AMR, Grid, FILE=file, SWAP=swap, NCPU=ncpu, ICPU=icpu,
;	VERBOSE=verbose 
;
; OPTIONAL INPUTS:
;	FILE:   if set, input the scalar string containing the name of
;	        the file to be read. Otherwise, a PICKFILE widget is
;	        launched.  
;
;       SWAP:   if set, reverse the bit ordering (Little Endian versus
;               Big Endian)
;
;       ICPU:   first cpu file to be read. Default: 1.
;
;       NCPU:   number of cpu files to read, starting from
;               icpu. Default: all files from icpu to ncpu_max.  
;
; OUTPUTS:
;	Grid:   store the AMR tree in structure Grid.
;
; COMMON BLOCKS:
;       None.
;
; EXAMPLE:
;       To read on a SGI architecture a RAMSES AMR file created on a
;       COMPAQ Workstation, type:
;
;	        RD_AMR, Grid, file='amr_00001.out',/swap
;
;       If the file was generated on the same IEEE system, just type:
;
;               RD_AMR, Grid, file='amr_00001.out'
;
; MODIFICATION HISTORY:
; 	Written by:	Romain Teyssier, 01/01/2000.
;                       e-mail: Romain.Teyssier@cea.fr
;	Fevrier, 2001:	Comments and header added by Romain Teyssier.
;-
;###################################################
;###################################################
;###################################################
pro rd_amr, grid, file=file, swap=swap, icpu=icpu, verbose=verbose, nout=nout

IF N_PARAMS() NE 1 THEN BEGIN
    PRINT, 'Wrong number of arguments'
    DOC_LIBRARY,'rd_amr'
    RETURN
ENDIF

if not keyword_set(icpu) then icpu=0
if icpu eq 0 then jcpu=1 else jcpu=icpu

suffix=getcarnum(jcpu)
if not keyword_set(file) and not keyword_set(nout) then begin
    key='*amr*.out'+suffix(jcpu-1)
    file=DIALOG_PICKFILE(/READ,filter=key)
endif
if keyword_set(nout) then begin 
    suffnout=getcarnum(nout)
    file='output_'+suffnout(nout-1)+'/amr_'+suffnout(nout-1)+'.out'
endif
if not keyword_set(file) then return
base_offset=strpos(file,'.out')+4
file_base=strmid(file,0,base_offset)

; Free memory associated to grid
del_amr,grid

; Initialize header variables
ncpu_run=0L & ndim=0L & nx=0L & ny=0L & nz=0L
nlevelmax=0L & ngridmax=0L & nboundary=0L & ngridactual=0L & nstep=0L 
noutput=0L & boxlen=0.0d0 & t=0.0d0
iout=0L & ifout=0L
aexp=0.0d0 & hexp=0.0d0 & aexp_old=0.0d0 & epot_tot_int=0.0d0
epot_tot_old=0.0d0
const=0.0d0 & mass_tot_0=0.0d0 & rho_tot=0.0d0
omega_m=0.0d0 & omega_l=0.0d0 & omega_k=0.0d0 & omega_b=0.0d0 & h0=0.0d0
aexp_ini=0.0d0 & mass_sph=0.0d0
headf=0L & tailf=0L & numbf=0L & used_mem=0L & used_mem_tot=0L

; Read first file to get header
print,'Reading file ',trim(file_base)
file=trim(file_base+suffix(jcpu-1))
openr,1,file,/f77_unformatted,swap_endian=swap
readu,1,ncpu_run
readu,1,ndim
readu,1,nx,ny,nz
ncoarse=nx*ny*nz
readu,1,nlevelmax
readu,1,ngridmax
readu,1,nboundary
readu,1,ngridactual
readu,1,boxlen
readu,1,noutput,iout,ifout
tout=dblarr(noutput)
aout=dblarr(noutput)
readu,1,tout
readu,1,aout
readu,1,t
dtold=dblarr(nlevelmax)
dtnew=dblarr(nlevelmax)
readu,1,dtold
readu,1,dtnew
readu,1,nstep,nstep_coarse
readu,1,const,mass_tot_0,rho_tot
readu,1,omega_m,omega_l,omega_k,omega_b,h0,aexp_ini
readu,1,aexp,hexp,aexp_old,epot_tot_int,epot_tot_old
readu,1,mass_sph
close,1

; Write header to screen
print,'ncpu      =',ncpu_run
print,'ndim      =',ndim
print,'nlevelmax =',nlevelmax
print,'nstep     =',nstep
print,'boxlen    =',boxlen
print,'time      =',t
if(hexp gt 0.0)then begin ; detect cosmo run
print,'aexp      =',aexp
print,'omega_m   =',omega_m
print,'omega_l   =',omega_l
print,'omega_k   =',omega_k
print,'omega_b   =',omega_b
endif
if nboundary eq 0 then begin
    print,"Periodic boundary conditions"
endif

; Allocate arrays
ncpu=ncpu_run
icpumin=1L & icpumax=ncpu & listmax=ncpu
if icpu gt 0 then begin
    icpumin=icpu & icpumax=icpu & listmax=ncpu+nboundary
endif
suffix=getcarnum(ncpu_run)
ncell=2L^ndim
nnbor=2L*ndim
ngrid=LONARR(nlevelmax,listmax)
level=PTRARR(nlevelmax,listmax)
headl=LONARR(ncpu,nlevelmax)
taill=LONARR(ncpu,nlevelmax)
numbl=LONARR(ncpu,nlevelmax)
if(nboundary gt 0)then begin
    headb=LONARR(nboundary,nlevelmax)
    tailb=LONARR(nboundary,nlevelmax)
    numbb=LONARR(nboundary,nlevelmax)
endif
numbtot=LONARR(10,nlevelmax)
bound_key=DBLARR(ncpu+1)
son=LONARR(ncoarse)
nbor=LONARR(ncoarse)
next=LONARR(ncoarse)
prev=LONARR(ncoarse)
flag1=LONARR(ncoarse)
cpu_map=LONARR(ncoarse)
xbound=[0d0,0d0,0d0]
ordering='                       '
; Read AMR grids

; Loop over cpu files
ngridtot=0L
list=0L
for jcpu=icpumin,icpumax do begin

    file=trim(file_base+suffix(jcpu-1))
    if keyword_set(verbose) then print,'Reading file ',trim(file)
    openr,1,file,/f77_unformatted,swap_endian=swap
    
; Read header
    readu,1,ncpu_run
    readu,1,ndim
    readu,1,nx,ny,nz
    readu,1,nlevelmax
    readu,1,ngridmax
    readu,1,nboundary
    readu,1,ngridactual
    readu,1,boxlen
    readu,1,noutput,iout,ifout
    readu,1,tout
    readu,1,aout
    readu,1,t
    readu,1,dtold
    readu,1,dtnew
    readu,1,nstep,nstep_coarse
    readu,1,const,mass_tot_0,rho_tot
    readu,1,omega_m,omega_l,omega_k,omega_b,h0,aexp_ini
    readu,1,aexp,hexp,aexp_old,epot_tot_int,epot_tot_old
    readu,1,mass_sph
    readu,1,headl
    readu,1,taill
    readu,1,numbl
    readu,1,numbot
    if nboundary gt 0 then begin
        readu,1,headb
        readu,1,tailb
        readu,1,numbb
        xbound=[double(nx/2),double(ny/2),double(nz/2)]
    endif
    readu,1,headf,tailf,numbf,used_mem,used_mem_tot
    readu,1,ordering
    if(strcompress(ordering) eq 'bisection ')then begin
        readu,1
        readu,1
        readu,1
        readu,1
        readu,1
    endif else begin
        readu,1,bound_key
    endelse
    readu,1,son
    readu,1,flag1
    readu,1,cpu_map
    
; Read fine levels
    nlevel=0L & ilevel=0L & ng=0L
    kcpumin=1L & kcpumax=nboundary+ncpu
    for ilevel=0L,nlevelmax-1L do begin
        for kcpu=kcpumin,kcpumax do begin
            if(kcpu le ncpu)then begin
                ng=numbl(kcpu-1,ilevel)
            endif else begin
                ng=numbb(kcpu-ncpu-1,ilevel)
            endelse
            if icpu eq 0 then begin
                if (kcpu eq jcpu) then ngrid(ilevel,kcpu-1)=ng
            endif else begin
                ngrid(ilevel,kcpu-1)=ng
            endelse
            if(ng gt 0)then begin
                ngridtot=ngridtot+ng
                if keyword_set(verbose) then begin
                    print,ilevel+1,ng,kcpu $
                      ,format='("Level ",i2," has ",i6," grids in proc",i6)'
                endif
                nlevel=nlevel+1L
                mesh2={ilevel:ilevel,nc:ng,xg:dblarr(ng,ndim),son:lonarr(ng,ncell),cpumap:lonarr(ng,ncell),nbor:lonarr(ng,nnbor),indgrid:lonarr(ng),next:lonarr(ng),prev:lonarr(ng),father:lonarr(ng),ison:lonarr(ng),cpu_father:lonarr(ng),cpu:lonarr(ng),ison_father:lonarr(ng),indgrid_father:lonarr(ng),cpu_next:lonarr(ng)}
                ii=lonarr(ng)   ; Define level structure
                xx=dblarr(ng)
                readu,1,ii      ; Read grid index
                mesh2.indgrid=ii

                readu,1,ii      ; Read next index
                mesh2.next=ii

                readu,1,ii      ; Read prev index
                mesh2.prev=ii
                for idim=0,ndim-1 do begin
                    readu,1,xx  ; Read grid center
                    mesh2.xg(*,idim)=xx-xbound(idim)
                endfor
                readu,1,ii      ; Read father index
                mesh2.father=ii
                for idim=0,2*ndim-1 do begin
                    readu,1,ii  ; Read nbor index
                    mesh2.nbor(*,idim)=ii
                 endfor
                for idim=0,2^ndim-1 do begin
                    readu,1,ii  ; Read son index
                    mesh2.son(*,idim)=ii
                endfor
                for idim=0,2^ndim-1 do begin
                    readu,1,ii  ; Read cpu map
                    mesh2.cpumap(*,idim)=ii
                endfor
                for idim=0,2^ndim-1 do begin
                    readu,1,ii  ; Read refinement map
                endfor
                
                ; Add a vector for the son index in the father oct
                mesh2.ison=mesh2.indgrid
                mesh2.ison(*)=0
                mesh2.ison_father=mesh2.indgrid
                mesh2.ison_father(*)=0
                mesh2.indgrid_father=mesh2.indgrid
                mesh2.indgrid_father(*)=0
                mesh2.cpu_father=mesh2.indgrid
                mesh2.cpu_father(*)=0
                mesh2.cpu=mesh2.indgrid
                mesh2.cpu(*)=0
                mesh2.cpu_next=mesh2.indgrid
                mesh2.cpu_next(*)=-1
                mesh2.next=mesh2.indgrid
                mesh2.next(*)=-1
                xx=0d0
                if icpu eq 0 then begin
                    if(kcpu eq jcpu)then begin
                        pc=ptr_new(mesh2)
                        level(ilevel,jcpu-1)=pc
                    endif
                endif else begin
                    pc=ptr_new(mesh2)
                    level(ilevel,kcpu-1)=pc
                endelse
                list=list+1L
            endif
        endfor
    endfor    
    close,1
endfor 
ngrid=ngrid[0:nlevelmax-1,0L:listmax-1L]
level=level[0:nlevelmax-1,0L:listmax-1L]
ncoarse=nx*ny*nz
grid={ncpu:listmax,ndim:ndim,time:t,aexp:aexp,nlevelmax:nlevelmax,boxlen:boxlen $
      ,ngridtot:ngridtot,ngrid:ngrid,level:level,ngridmax:ngridmax,ncoarse:ncoarse}

mesh2=0

return

bad_luck:  print,'I/O Error, exiting...'
close,1
mesh2=0

return

end
;###################################################
;###################################################
;###################################################




