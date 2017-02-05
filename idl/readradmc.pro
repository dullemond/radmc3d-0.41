;==========================================================================
;==========================================================================
;==========================================================================
;           IDL PACKAGE TO READ AND PLOT RESULTS FROM RADMC-3D
;==========================================================================
;==========================================================================
;==========================================================================





;==========================================================================
;                        ROUTINES FOR IMAGES
;==========================================================================


;-----------------------------------------------------------------
;              READ THE RECTANGULAR TELESCOPE IMAGE
;-----------------------------------------------------------------
function readimage,ext=ext,file=file,iounit=iounit
nx=0
ny=0
nf=0
sizepix_x=0.d0
sizepix_y=0.d0
if not keyword_set(iounit) then funit=1 else funit=iounit
if not keyword_set(iounit) then begin
   ;;
   ;; Read from normal file, so make filename
   ;;
   if not keyword_set(file) then begin
      if n_elements(ext) eq 0 then begin
         file='image.out'
      endif else begin
         file='image_'+strcompress(string(ext),/remove_all)+'.out'
      endelse
   endif
   str=findfile(file)
   if(str(0) ne file) then begin
      print,'Sorry, cannot find ',file
      print,'Presumably radmc3d exited without succes.'
      print,'See above for possible error messages of radmc3d!'
      stop
   endif
   openr,funit,file
endif else begin
   ;;
   ;; Read from radmc3d directly via standard output
   ;; 
   ;; First ask the image from radmc3d
   ;;
   printf,funit,'writeimage'
   flush,funit
   ;;
   ;; Check if RADMC-3D has already responded
   ;;
   if eof(funit) then begin
      print,'Waiting for RADMC-3D...'
      wait,1
      if eof(funit) then begin
         print,'Aborting: RADMC-3D does not respond...'
         stop
      endif
   endif
   ;;
   ;; Now we can read...
   ;;
endelse
;
; Read the image
;
iformat=0
readf,funit,iformat
if iformat lt 1 or iformat gt 4 then begin
   print,'ERROR: File format of ',file,' not recognized.'
   stop
endif
if iformat eq 1 or iformat eq 3 then begin
   radian = (1 eq 0)
endif else begin
   radian = (1 eq 1)
endelse
if iformat eq 1 or iformat eq 2 then begin
   stokes = (1 eq 0)
endif else begin
   stokes = (1 eq 1)
endelse
readf,funit,nx,ny
readf,funit,nf
readf,funit,sizepix_x,sizepix_y
lambda=dblarr(nf)
readf,funit,lambda
if stokes then begin
   image=dblarr(4,nx,ny,nf)
endif else begin
   image=dblarr(nx,ny,nf)
endelse
readf,funit,image
if funit eq 1 then close,1
;
; If the image contains all four Stokes vector components,
; then it is useful to transpose the image array such that
; the Stokes index is the third index, so that the first
; two indices remain x and y
;
if stokes then begin
   if nf gt 1 then begin
      image = transpose(image,[1,2,0,3])
   endif else begin
      image = transpose(image,[1,2,0])
   endelse
endif
;
; Compute the flux in this image as seen at 1 pc
;
flux=dblarr(nf)
if stokes then begin
   for ix=0,nx-1 do for iy=0,ny-1 do flux=flux+image[ix,iy,0,*]
endif else begin
   for ix=0,nx-1 do for iy=0,ny-1 do flux=flux+image[ix,iy,*]
endelse
flux=flux*sizepix_x*sizepix_y
if not radian then begin
   pc=3.0857200d+18
   flux=flux/pc^2
endif
;
; ADDED 13.12.06:
; Compute the x- and y- coordinates
;
x=((dindgen(nx)+0.5)/(nx*1.d0)-0.5d0)*sizepix_x*nx
y=((dindgen(ny)+0.5)/(ny*1.d0)-0.5d0)*sizepix_y*ny
;
; Return all
;
return,{nx:nx,ny:ny,nrfr:nf,sizepix_x:sizepix_x,sizepix_y:sizepix_y,$
        image:image,flux:flux,x:x,y:y,lambda:lambda,radian:radian,stokes:stokes}
end



;--------------------------------------------------------------------------
;                 CALL RADMC-3D FOR IMAGE, READ AND PLOT
;--------------------------------------------------------------------------
pro doimage,incl=incl,phi=phi,npix=npix,sizecm=size,sizeau=sizeau,$
            sizepc=sizepc,zoomau=zoomau,zoompc=zoompc,$
            posang=posang,lines=lines,secondorder=secondorder,$
            ifreq=ifreq,lambda=lambda,nostar=nostar,$
            log=log,pc=pc,au=au,contour=contour,nlevels=nlevels,$
            noimage=noimage,pos=pos,thrcol=thrcol,thres=thres,maxlog=maxlog,$
            saturate=saturate,xticklen=xticklen,yticklen=yticklen,$
            xminor=xminor,yminor=yminor,plottau=plottau,lgrange=lgrange,$
            image=image,dcatch=dcatch,options=options,fixseed=fixseed,$
            nphot=nphot,_EXTRA=_extra
makeimage,incl=incl,phi=phi,npix=npix,sizecm=size,sizeau=sizeau,$
              sizepc=sizepc,posang=posang,lines=lines,$
              ifreq=ifreq,lambda=lambda,nostar=nostar,$
              zoomau=zoomau,zoompc=zoompc,options=options,$
              secondorder=secondorder,dcatch=dcatch,fixseed=fixseed,$
              nphot=nphot
image=readimage()
plotimage,image,log=log,pc=pc,au=au,contour=contour,nlevels=nlevels,$
        noimage=noimage,pos=pos,thrcol=thrcol,thres=thres,maxlog=maxlog,$
        saturate=saturate,xticklen=xticklen,yticklen=yticklen,$
        xminor=xminor,yminor=yminor,plottau=plottau,lgrange=lgrange,$
        _EXTRA=_extra
end


;--------------------------------------------------------------------------
;                   CALL RADMC-3D TO MAKE THE IMAGE
;--------------------------------------------------------------------------
pro makeimage,incl=incl,phi=phi,npix=npix,sizecm=sizecm,sizeau=sizeau,$
              sizepc=sizepc,posang=posang,nofluxcons=nofluxcons,dcatch=dcatch,$
              pointcm=pointcm,pointau=pointau,pointpc=pointpc,lines=lines,$
              ifreq=ifreq,lambda=lambda,nostar=nostar,nrimages=nrimages,$
              zoomau=zoomau,zoompc=zoompc,plottau=plottau,iounit=iounit,$
              linelist=linelist,iline=iline,ispec=ispec,vkms=vkms,$
              secondorder=secondorder,radmc3d=radmc3d,options=options,$
              fixseed=fixseed,nphot=nphot
AU  = 1.496d13       ; Astronomical Unit       [cm]
pc  = 3.0857200d18   ; Parsec                  [cm]
cc  = 2.9979245800000d10      ; Light speed             [cm/s]
;cc  = 2.9979d10      ; Light speed             [cm/s]
if not keyword_set(radmc3d) then begin
   tmp=file_search('./radmc3d',count=count1)
   if count1 eq 0 then begin
      radmc3d='radmc3d'
   endif else begin
      radmc3d='./radmc3d'
   endelse
endif
if n_elements(iline) gt 0 or n_elements(ispec) gt 0 or n_elements(vkms) gt 0  then begin
   linelamspec=1
   if n_elements(iline) eq 0 then iline=1
   if n_elements(ispec) eq 0 then ispec=1
   if n_elements(vkms) eq 0 then vkms=0.d0
endif else begin
   linelamspec=0
endelse
if n_elements(incl) eq 0 then incl=0.d0
if n_elements(phi) eq 0 then phi=0.d0
if n_elements(npix) eq 0 then npix=100
if n_elements(sizecm) ne 0 then sizeau=sizecm/AU
if n_elements(sizepc) ne 0 then sizeau=sizepc*pc/AU
if n_elements(zoompc) ne 0 then zoomau=zoompc*pc/AU
if n_elements(pointcm) ne 0 then pointau=pointcm/AU
if n_elements(pointpc) ne 0 then pointau=pointpc*pc/AU
if n_elements(pointau) ne 0 then begin
   if n_elements(pointau[*,0]) ne 3 then begin
      print,'ERROR in makeimage: point** keywords must be array of 3 numbers (x,y,z)'
      return
   endif
endif
if n_elements(pointau) eq 0 then pointau=[0.,0.,0.]
if n_elements(posang) eq 0 then posang=0.d0
if n_elements(plottau) eq 0 then plottau=0
;;
;; If we do not have iline, ispec, vkms (i.e. lines), then we must specify
;; the wavelength in the following manner: 
;;
if linelamspec eq 0 then begin
   if n_elements(ifreq) gt 0 and n_elements(lambda) gt 0 then begin
      print,'ERROR: Cannot specify both ifreq and lambda'
      stop
   endif
   if n_elements(ifreq) eq 0 and n_elements(lambda) eq 0 then begin
      print,'ERROR: Must specify either ifreq or lambda'
      stop
   endif
   if n_elements(ifreq) then begin
      nlam=n_elements(ifreq)
   endif else begin
      nlam=n_elements(lambda)
   endelse
   if nlam eq 0 then begin
      print,'ERROR: You must specify the wavelength or frequency-bin of the image'
      return
   endif
   ;;
   ;; Now check if we have multiple colors or a single color
   ;;
   if nlam eq 1 then begin
      spawn,'\rm color_inus.inp >& /dev/null'
      spawn,'\rm camera_wavelength_micron.inp >& /dev/null'
      if keyword_set(ifreq) then begin
         lamstr0 = 'ilambda '
         lamstr1 = strcompress(string(ifreq[0]),/remove_all)
      endif else begin
         lamstr0 = 'lambda '
         lamstr1 = strcompress(string(lambda[0],format='(E21.14)'),/remove_all)
      endelse
   endif else begin
      if n_elements(ifreq) ne 0 then begin
         lamstr0 = 'color'
         lamstr1 = ''
         openw,1,'color_inus.inp'
         printf,1,1
         printf,1,nlam
         printf,1,' '
         for ilam=0,nlam-1 do printf,1,ifreq[ilam]
         close,1
      endif else begin
         lamstr0 = 'loadlambda'
         lamstr1 = ''
         openw,1,'camera_wavelength_micron.inp'
         printf,1,nlam
         for ilam=0,nlam-1 do printf,1,lambda[ilam]
         close,1
      endelse
   endelse
endif else begin
   spawn,'\rm -f camera_wavelength_micron.inp >& /dev/null'
   spawn,'\rm -f color_inus.inp >& /dev/null'
   lamstr0 = ''
   lamstr1 = ''
endelse
;;
;; Now check if we should make multiple images
;;
if n_elements(pointau[0,*]) gt 2 and n_elements(nrimages) eq 0 then nrimages=n_elements(pointau[0,*])
if n_elements(sizeau) gt 2 and n_elements(nrimages) eq 0 then nrimages=n_elements(sizeau)
if n_elements(incl) gt 2 and n_elements(nrimages) eq 0 then nrimages=n_elements(incl)
if n_elements(phi) gt 2 and n_elements(nrimages) eq 0 then nrimages=n_elements(phi)
if n_elements(posang) gt 2 and n_elements(nrimages) eq 0 then nrimages=n_elements(posang)
if n_elements(nrimages) eq 0 then nrimages=1
if nrimages gt 1 then begin
   ;;
   ;; Temprary
   ;;
   if not keyword_set(sizeau) then begin
      if keyword_set(zoomau) then begin
         print,'ERROR: For now zoomau is not enabled for movies'
         stop
      endif else begin
         print,'ERROR: Need to specify sizeau or sizepc for movie'
         stop
      endelse
   endif
   ;;
   ;; Yes. Now make arrays for the observer's orientation
   ;;
   pt = dblarr(3,nrimages)
   hs = dblarr(nrimages)
   pa = dblarr(nrimages)
   th = dblarr(nrimages)
   ph = dblarr(nrimages)
   xx = dindgen(nrimages)/(nrimages-1.d0)
   if n_elements(pointau[0,*]) eq 1 then begin
      pt[0,*] = pointau[0,0]
      pt[1,*] = pointau[1,0]
      pt[2,*] = pointau[2,0]
   endif else begin
      if n_elements(pointau[0,*]) eq 2 then begin
         pt[0,*] = pointau[0,0] + xx*(pointau[0,1]-pointau[0,0])
         pt[1,*] = pointau[1,0] + xx*(pointau[1,1]-pointau[1,0])
         pt[2,*] = pointau[2,0] + xx*(pointau[2,1]-pointau[2,0])
      endif else begin
         if n_elements(pointau[0,*]) eq nrimages then begin
            pt[0,*] = pointau[0,*]
            pt[1,*] = pointau[1,*]
            pt[2,*] = pointau[2,*]
         endif else begin
            print,'ERROR: pointau/pc/cm has wrog number of elements'
            return
         endelse
      endelse
   endelse
   if n_elements(sizeau) ne 0 and n_elements(zoomau) ne 0 then begin
      print,'ERROR: Cannot specify size and zoom simultaneously'
      return
   endif
   if n_elements(sizeau) eq 1 then begin
      hs[*] = au * sizeau[0] / 2.
   endif else begin
      if n_elements(sizeau) eq 2 then begin
         hs[*] = au * ( sizeau[0] + xx*(sizeau[1]-sizeau[0]) ) / 2.
      endif else begin
         if n_elements(sizeau) eq nrimages then begin
            hs[*] = au * sizeau[*] / 2.
         endif else begin
            hs[*] = 0.d0    ; Telling radmc3d to use its own estimate
         endelse
      endelse
   endelse
   if n_elements(posang) eq 1 then begin
      pa[*] = posang[0]
   endif else begin
      if n_elements(posang) eq 2 then begin
         pa[*] = ( posang[0] + xx*(posang[1]-posang[0]) ) 
      endif else begin
         if n_elements(posang) eq nrimages then begin
            pa[*] = posang[*]
         endif else begin
            print,'ERROR: posang has wrong number of elements'
            return
         endelse
      endelse
   endelse
   if n_elements(incl) eq 1 then begin
      th[*] = incl[0]
   endif else begin
      if n_elements(incl) eq 2 then begin
         th[*] = ( incl[0] + xx*(incl[1]-incl[0]) )
      endif else begin
         if n_elements(incl) eq nrimages then begin
            th[*] = incl[*]
         endif else begin
            print,'ERROR: incl has wrong number of elements'
            return
         endelse
      endelse
   endelse
   if n_elements(phi) eq 1 then begin
      ph[*] = phi[0]
   endif else begin
      if n_elements(phi) eq 2 then begin
         ph[*] = ( phi[0] + xx*(phi[1]-phi[0]) )
      endif else begin
         if n_elements(phi) eq nrimages then begin
            ph[*] = phi[*]
         endif else begin
            print,'ERROR: phi has wrong number of elements'
            return
         endelse
      endelse
   endelse
   ;;
   ;; Now write the movie.inp file
   ;;
   openw,1,'movie.inp'
   printf,1,1
   printf,1,nrimages
   for i=0,nrimages-1 do printf,1,pt[0:2,i],hs[i],hs[i],pa[i],th[i],ph[i],1d99
   close,1
endif else begin
   spawn,'\rm movie.inp >& /dev/null'
   ;;
   ;; Do some checks 
   ;;
   if n_elements(sizeau) ne 0 and n_elements(zoomau) ne 0 then begin
      print,'ERROR: Cannot specify size and zoom simultaneously'
      return
   endif
   if n_elements(zoomau) ne 0 and n_elements(zoomau) ne 4 then begin
      print,'ERROR: Zoomau must have 4 elements.'
      return
   endif
endelse
;;
;; Now make the radmc3d command
;; NOTE: This command is always built, but only really used if the
;;       iounit is not set.
;;
if nrimages le 1 then begin
   command=radmc3d+" image npix "+strcompress(string(npix),/remove_all)+$
           " "+lamstr0+lamstr1+$
           " posang "+strcompress(string(posang),/remove_all)+$
           " incl "+strcompress(string(incl),/remove_all)+$
           " phi "+strcompress(string(phi),/remove_all)+$
           " pointau "+strcompress(string(pointau[0]),/remove_all)+" "+$
           strcompress(string(pointau[1]),/remove_all)+" "+$
           strcompress(string(pointau[2]),/remove_all)
   if keyword_set(sizeau) then command = command+" sizeau "+strcompress(string(sizeau),/remove_all)
   if keyword_set(zoomau) then command = command+" zoomau "+   $
              strcompress(string(zoomau[0]),/remove_all)+" "+  $
              strcompress(string(zoomau[1]),/remove_all)+" "+  $
              strcompress(string(zoomau[2]),/remove_all)+" "+  $
              strcompress(string(zoomau[3]),/remove_all)
   if keyword_set(nofluxcons) then command = command+" nofluxcons" else $
      command = command+" fluxcons"
   if keyword_set(nostar) then command = command+" nostar" else $
      command = command+" inclstar"
   if keyword_set(secondorder) then command = command+" secondorder"
   if keyword_set(fixseed) then command = command+" resetseed"
   if keyword_set(plottau) then command = command+" tracetau"
   if n_elements(ifreq) gt 1 then command = command+" loadcolor"
   if n_elements(lambda) gt 1 then command = command+" loadlambda"
   if keyword_set(nphot) then command = command + " nphot_scat "+  $
           strcompress(string(nphot),/remove_all)
   if keyword_set(lines) then command = command + " inclline"
   if keyword_set(linelist) then command = command + " linelist"
   if keyword_set(dcatch) then command = command + " doppcatch"
   if keyword_set(iline) then command = command + " iline "+  $
           strcompress(string(iline),/remove_all)
   if keyword_set(ispec) then command = command + " imolspec "+  $
           strcompress(string(ispec),/remove_all)
   if keyword_set(vkms) then command = command + " vkms "+  $
           strcompress(string(vkms),/remove_all)
   if (n_elements(options) gt 0) then begin
      no = n_elements(options)
      for io=0,no-1 do begin
         command = command + options[io] + " "
      endfor
   endif
endif else begin
   command=radmc3d+" movie npix "+strcompress(string(npix),/remove_all)+" "+lamstr0+lamstr1
endelse
print,command
;;
;; Are we calling radmc3d or are we communicating with radmc3d as a child process?
;;
if keyword_set(iounit) then begin
   ;;
   ;; Do a check
   ;;
   if nrimages gt 1 then begin
      print,'SORRY: Movies not supported yet in child mode.'
      stop
   endif
   ;;
   ;; Radmc3D is a child process, so communicate with it via the iounit 
   ;; file unit port
   ;;
   printf,iounit,'image'
   printf,iounit,'npix'
   printf,iounit,strcompress(string(npix),/remove_all)
   if lamstr0 ne '' then printf,iounit,lamstr0
   if lamstr1 ne '' then printf,iounit,lamstr1
   printf,iounit,'posang'
   printf,iounit,strcompress(string(posang),/remove_all)
   printf,iounit,'incl'
   printf,iounit,strcompress(string(incl),/remove_all)
   printf,iounit,'phi'
   printf,iounit,strcompress(string(phi),/remove_all)
   printf,iounit,'pointau'
   printf,iounit,strcompress(string(pointau[0]),/remove_all)
   printf,iounit,strcompress(string(pointau[1]),/remove_all)
   printf,iounit,strcompress(string(pointau[2]),/remove_all)
   if keyword_set(sizeau) then begin
      printf,iounit,'sizeau'
      printf,iounit,strcompress(string(sizeau),/remove_all)
   endif
   if keyword_set(zoomau) then begin
      printf,iounit,'zoomau'
      printf,iounit,strcompress(string(zoomau[0]),/remove_all)
      printf,iounit,strcompress(string(zoomau[1]),/remove_all)
      printf,iounit,strcompress(string(zoomau[2]),/remove_all)
      printf,iounit,strcompress(string(zoomau[3]),/remove_all)
   endif
   if keyword_set(nofluxcons) then begin
      printf,iounit,'nofluxcons'
   endif else begin
      printf,iounit,'fluxcons'
   endelse
   if keyword_set(nostar) then begin
      printf,iounit,'nostar'
   endif else begin
      printf,iounit,'inclstar'
   endelse
   if keyword_set(secondorder) then begin
      printf,iounit,'secondorder'
   endif
   if keyword_set(fixseed) then begin
      printf,iounit,'resetseed'
   endif
   if keyword_set(plottau) then begin
      printf,iounit,'tracetau'
   endif else begin
      printf,iounit,'tracenormal'
   endelse
   if keyword_set(nphot) then begin
      printf,iounit,'nphot_scat'
      printf,iounit,strcompress(string(nphot),/remove_all)
   endif
   ;;if n_elements(ifreq) gt 1 then printf,iounit,'loadcolor'
   ;;if n_elements(lambda) gt 1 then printf,iounit,'loadlambda'
   if keyword_set(lines) then printf,iounit,'inclline'
   if keyword_set(linelist) then printf,iounit,'linelist'
   if keyword_set(dcatch) then printf,iounit,'doppcatch'
   if keyword_set(iline) then begin
      printf,iounit,'iline'
      printf,iounit,strcompress(string(iline),/remove_all)
   endif 
   if keyword_set(ispec) then begin
      printf,iounit,'imolspec'
      printf,iounit,strcompress(string(ispec),/remove_all)
   endif
   if keyword_set(vkms) then begin
      printf,iounit,'vkms'
      printf,iounit,strcompress(string(vkms),/remove_all)
   endif
   if (n_elements(options) gt 0) then begin
      no = n_elements(options)
      for io=0,no-1 do begin
         printf,iounit,options[io]
      endfor
   endif
   printf,iounit,'enter'
   flush,iounit
endif else begin
   ;;
   ;; Radmc3D must be called from a shell
   ;;
   ;;
   ;; Call radmc3d
   ;;
   ;spawn,'../src/'+command
   spawn,command
endelse
;;
end

;--------------------------------------------------------------------------
;        CALL RADMC-3D TO MAKE THE IMAGE FOR A LOCAL OBSERVER
;--------------------------------------------------------------------------
pro makeimagelocalobs,npix=npix,posang=posang,nofluxcons=nofluxcons,$
              pointcm=pointcm,pointau=pointau,pointpc=pointpc,lines=lines,$
              locobscm=locobscm,locobsau=locobsau,locobspc=locobspc,sizeradian=sizeradian,$
              ifreq=ifreq,lambda=lambda,nostar=nostar,$
              plottau=plottau,iounit=iounit,secondorder=secondorder,$
              linelist=linelist,iline=iline,ispec=ispec,vkms=vkms,$
              radmc3d=radmc3d,options=options
AU  = 1.496d13       ; Astronomical Unit       [cm]
pc  = 3.0857200d18   ; Parsec                  [cm]
cc  = 2.9979245800000d10      ; Light speed             [cm/s]
;cc  = 2.9979d10      ; Light speed             [cm/s]
if not keyword_set(radmc3d) then begin
   tmp=file_search('./radmc3d',count=count1)
   if count1 eq 0 then begin
      radmc3d='radmc3d'
   endif else begin
      radmc3d='./radmc3d'
   endelse
endif
if n_elements(iline) gt 0 or n_elements(ispec) gt 0 or n_elements(vkms) gt 0  then begin
   linelamspec=1
   if n_elements(iline) eq 0 then iline=1
   if n_elements(ispec) eq 0 then ispec=1
   if n_elements(vkms) eq 0 then vkms=0.d0
endif else begin
   linelamspec=0
endelse
if n_elements(npix) eq 0 then npix=100
if n_elements(pointcm) ne 0 then pointau=pointcm/AU
if n_elements(pointpc) ne 0 then pointau=pointpc*pc/AU
if n_elements(pointau) ne 0 then begin
   if n_elements(pointau[*,0]) ne 3 then begin
      print,'ERROR in makeimage: point** keywords must be array of 3 numbers (x,y,z)'
      return
   endif
endif
if n_elements(pointau) eq 0 then pointau=[0.,0.,0.]
if n_elements(locobscm) ne 0 then locobsau=locobscm/AU
if n_elements(locobspc) ne 0 then locobsau=locobspc*pc/AU
if n_elements(locobsau) ne 0 then begin
   if n_elements(locobsau[*,0]) ne 3 then begin
      print,'ERROR in makeimage: locobs** keywords must be array of 3 numbers (x,y,z)'
      return
   endif
endif
if n_elements(locobsau) eq 0 then locobsau=[0.,-1000.,0.]
if n_elements(posang) eq 0 then posang=0.d0
if n_elements(sizeradian) eq 0 then sizeradian=0.8d0
if n_elements(plottau) eq 0 then plottau=0
;;
;; If we do not have iline, ispec, vkms (i.e. lines), then we must specify
;; the wavelength in the following manner: 
;;
if linelamspec eq 0 then begin
   if n_elements(ifreq) gt 0 and n_elements(lambda) gt 0 then begin
      print,'ERROR: Cannot specify both ifreq and lambda'
      stop
   endif
   if n_elements(ifreq) eq 0 and n_elements(lambda) eq 0 then begin
      print,'ERROR: Must specify either ifreq or lambda'
      stop
   endif
   if n_elements(ifreq) then begin
      nlam=n_elements(ifreq)
   endif else begin
      nlam=n_elements(lambda)
   endelse
   if nlam eq 0 then begin
      print,'ERROR: You must specify the wavelength or frequency-bin of the image'
      return
   endif
   ;;
   ;; Now check if we have multiple colors or a single color
   ;;
   if nlam eq 1 then begin
      spawn,'\rm color_inus.inp >& /dev/null'
      spawn,'\rm camera_wavelength_micron.inp >& /dev/null'
      if keyword_set(ifreq) then begin
         lamstr0 = 'ilambda '
         lamstr1 = strcompress(string(ifreq[0]),/remove_all)
      endif else begin
         lamstr0 = 'lambda '
         lamstr1 = strcompress(string(lambda[0],format='(E21.14)'),/remove_all)
      endelse
   endif else begin
      if n_elements(ifreq) ne 0 then begin
         lamstr0 = 'color'
         lamstr1 = ''
         openw,1,'color_inus.inp'
         printf,1,1
         printf,1,nlam
         printf,1,' '
         for ilam=0,nlam-1 do printf,1,ifreq[ilam]
         close,1
      endif else begin
         lamstr0 = 'loadlambda'
         lamstr1 = ''
         openw,1,'camera_wavelength_micron.inp'
         printf,1,nlam
         for ilam=0,nlam-1 do printf,1,lambda[ilam]
         close,1
      endelse
   endelse
endif else begin
   spawn,'\rm -f camera_wavelength_micron.inp >& /dev/null'
   spawn,'\rm -f color_inus.inp >& /dev/null'
   lamstr0 = ''
   lamstr1 = ''
endelse
;;
;; Now make the radmc3d command
;; NOTE: This command is always built, but only really used if the
;;       iounit is not set.
;;
command=radmc3d+" image npix "+strcompress(string(npix),/remove_all)+$
           " "+lamstr0+lamstr1+$
           " posang "+strcompress(string(posang),/remove_all)+$
           " pointau "+strcompress(string(pointau[0]),/remove_all)+" "+$
           strcompress(string(pointau[1]),/remove_all)+" "+$
           strcompress(string(pointau[2]),/remove_all)+$
           " locobsau "+strcompress(string(locobsau[0]),/remove_all)+" "+$
           strcompress(string(locobsau[1]),/remove_all)+" "+$
           strcompress(string(locobsau[2]),/remove_all)+$
           " sizeradian "+strcompress(string(sizeradian),/remove_all)
if keyword_set(nofluxcons) then command = command+" nofluxcons" else $
      command = command+" fluxcons"
if keyword_set(nostar) then command = command+" nostar" else $
      command = command+" inclstar"
if keyword_set(secondorder) then command = command+" secondorder"
if keyword_set(plottau) then command = command+" tracetau"
if n_elements(ifreq) gt 1 then command = command+" loadcolor"
if n_elements(lambda) gt 1 then command = command+" loadlambda"
if keyword_set(lines) then command = command + " inclline"
if keyword_set(linelist) then command = command + " linelist"
if keyword_set(iline) then command = command + " iline "+  $
           strcompress(string(iline),/remove_all)
if keyword_set(ispec) then command = command + " imolspec "+  $
           strcompress(string(ispec),/remove_all)
if keyword_set(vkms) then command = command + " vkms "+  $
           strcompress(string(vkms),/remove_all)
if (n_elements(options) gt 0) then begin
   no = n_elements(options)
   for io=0,no-1 do begin
      command = command + options[io] + " "
   endfor
endif
print,command
;;
;; Are we calling radmc3d or are we communicating with radmc3d as a child process?
;;
if keyword_set(iounit) then begin
   ;;
   ;; Radmc3D is a child process, so communicate with it via the iounit 
   ;; file unit port
   ;;
   printf,iounit,'image'
   printf,iounit,'npix'
   printf,iounit,strcompress(string(npix),/remove_all)
   if lamstr0 ne '' then printf,iounit,lamstr0
   if lamstr1 ne '' then printf,iounit,lamstr1
   printf,iounit,'posang'
   printf,iounit,strcompress(string(posang),/remove_all)
   printf,iounit,'pointau'
   printf,iounit,strcompress(string(pointau[0]),/remove_all)
   printf,iounit,strcompress(string(pointau[1]),/remove_all)
   printf,iounit,strcompress(string(pointau[2]),/remove_all)
   printf,iounit,'locobsau'
   printf,iounit,strcompress(string(locobsau[0]),/remove_all)
   printf,iounit,strcompress(string(locobsau[1]),/remove_all)
   printf,iounit,strcompress(string(locobsau[2]),/remove_all)
   printf,iounit,'sizeradian'
   printf,iounit,strcompress(string(sizeradian),/remove_all)
   if keyword_set(nofluxcons) then begin
      printf,iounit,'nofluxcons'
   endif else begin
      printf,iounit,'fluxcons'
   endelse
   if keyword_set(nostar) then begin
      printf,iounit,'nostar'
   endif else begin
      printf,iounit,'inclstar'
   endelse
   if keyword_set(secondorder) then begin
      printf,iounit,'secondorder'
   endif
   if keyword_set(plottau) then begin
      printf,iounit,'tracetau'
   endif else begin
      printf,iounit,'tracenormal'
   endelse
   ;;if n_elements(ifreq) gt 1 then printf,iounit,'loadcolor'
   ;;if n_elements(lambda) gt 1 then printf,iounit,'loadlambda'
   if keyword_set(lines) then printf,iounit,'inclline'
   if keyword_set(linelist) then printf,iounit,'linelist'
   if keyword_set(iline) then begin
      printf,iounit,'iline'
      printf,iounit,strcompress(string(iline),/remove_all)
   endif 
   if keyword_set(ispec) then begin
      printf,iounit,'imolspec'
      printf,iounit,strcompress(string(ispec),/remove_all)
   endif
   if keyword_set(vkms) then begin
      printf,iounit,'vkms'
      printf,iounit,strcompress(string(vkms),/remove_all)
   endif
   if (n_elements(options) gt 0) then begin
      no = n_elements(options)
      for io=0,no-1 do begin
         printf,iounit,options[io]
      endfor
   endif
   printf,iounit,'enter'
   flush,iounit
endif else begin
   ;;
   ;; Radmc3D must be called from a shell
   ;;
   ;;
   ;; Call radmc3d
   ;;
   ;spawn,'../src/'+command
   spawn,command
endelse
;;
end

;--------------------------------------------------------------------------
;                 ASK RADMC3D TO WRITE A GRID FILE
;--------------------------------------------------------------------------
pro radmc3d_write_gridfile,iounit=iounit,radmc3d=radmc3d
if not keyword_set(radmc3d) then begin
   tmp=file_search('./radmc3d',count=count1)
   if count1 eq 0 then begin
      radmc3d='radmc3d'
   endif else begin
      radmc3d='./radmc3d'
   endelse
endif
if keyword_set(iounit) then begin
   printf,iounit,'writegridfile'
   printf,iounit,'enter'
   printf,iounit,'respondwhenready'
   printf,iounit,'enter'
   flush,iounit
   idum=0
   readf,iounit,idum
endif else begin
   command = radmc3d+' writegridfile'
   spawn,command
endelse
end

;-----------------------------------------------------------------
;                   PLOT THE RECTANGULAR IMAGE
;
; This routine plots the image read in with readimage(). It can 
; plot grey-scale images or color-table monochromatic images. It
; can also plot color images if you provide an image with 3 
; channels (RGB). 
;
; ARGUMENTS:
;  a                  The image to be plotted
;  log                If set, then plot the log of the image
;  pc                 If set, then use parsec as unit of the axes
;  au                 If set, then use AU as unit of the axes
;  contour            If set, make contour plot
;  nlevels            Nlevels of the contour plot
;  noimage            If set, only plot axes
;  pos                If set, determines the position on the canvas
;  maxlog             If /log, then maxlog sets max range (e.g. maxlog=6)
;  jpg                If set, then write a JPG file: image.jpg.
;                     The value of jpg is the size factor. If jpg==3, then
;                     for each pixel in the image, a 3x3 set of pixels on
;                     the JPG image are reserved. If jpg=-3, then the same
;                     scaling is applied, but instead of a grey-scale
;                     image, a red-temperature color table is used.
;                     Looks nicer.
;  mincol, maxcol     Limit the use of the color range [0,255]. 
;                     So if maxcol=100 you get darker images.
;  lgrange            Set specifically the range of the image brightness
;                     in terms of the log of the image values.
;  ilam               If a is a multi-color image, ilam specfies
;                     which wavelength you wish to use.
;  coltune            If set to 1, then rescale the brightness of 
;                     all channels the same value, to get the best
;                     color depth. If set to a 3-element array, you
;                     can directly specify the weight of each color.
;                     In this way you can really fine-tune the colors.
;  dev                If set, then automatically set some basic 
;                     device settings such that the image looks as
;                     it should look on 'x' or 'ps' devices. For
;                     X it in fact uses dev as the wset number.
;  zoom               If set, use this 'zoom' array (see zoomau and
;                     zoompc in makeimage) for the coordinates.
;
;-----------------------------------------------------------------
pro plotimage,a,log=log,pc=pc,au=au,contour=contour,nlevels=nlevels,$
        noimage=noimage,position=position,maxlog=maxlog,jpg=jpg,saturate=saturate,$
        xticklen=xticklen,yticklen=yticklen,mincol=mincol,maxcol=maxcol,$
        xminor=xminor,yminor=yminor,lgrange=lgrange,filenr=filenr,$
        ilam=ilam,coltune=coltune,dev=dev,jim=jim,zoom=zoom,$
        xsize=xsize,ysize=ysize,cmx=cmx,cmn=cmn,_EXTRA=_extra
common colcodes,green,yellow,red,blue,black,white
common colors,r_orig,g_orig,b_orig,r_curr,g_curr,b_curr
;
; Check if Stokes is switched on
;
stokes = a.stokes
;
; Image characteristics
;
if stokes then begin
   nx=n_elements(a.image[*,0,0,0])
   ny=n_elements(a.image[0,*,0,0])
endif else begin
   nx=n_elements(a.image[*,0,0])
   ny=n_elements(a.image[0,*,0])
endelse
;
; Display settings
;
if !d.name eq 'X' then begin
   if n_elements(position) ne 4 then position=[70,70,400,400]
   position=position*1.d0
   if nx lt ny then begin
      position[2]=position[3]*(1.d0*nx)/(1.d0*ny)
   endif
   if ny lt nx then begin
      position[3]=position[2]*(1.d0*ny)/(1.d0*nx)
   endif
   if not keyword_set(xsize) then xsize=position[0]+position[2]+30
   if not keyword_set(ysize) then ysize=position[1]+position[3]+30
endif else begin
   if n_elements(position) ne 4 then position=[0.1,0.1,0.8,0.8]
endelse
if n_elements(mincol) eq 0 then mincol=0
if n_elements(maxcol) eq 0 then maxcol=255
color=0
;
; Color settings
;
if stokes then begin
   nlam=n_elements(a.image[0,0,0,*])
endif else begin
   nlam=n_elements(a.image[0,0,*])
endelse
if nlam eq 1 then begin
   ilam=0
   color=0
endif else begin
   if nlam eq 3 and n_elements(ilam) eq 0 then begin
      ilam=[0,1,2]
      color=1
   endif
   if n_elements(ilam) eq 0 then begin
      print,'ERROR: Image has multiple colors. Must specify ilam.'
      return
   endif
   if min(ilam) lt 0 or max(ilam) ge nlam then begin
      print,'ERROR: ilam value(s) out of range.'
      return
   endif
endelse
nilam = n_elements(ilam)
if nilam ne 1 and nilam ne 3 then begin
   print,'ERROR: Either B/W or color images.'
   return
endif
if nilam eq 1 then color=0
if nilam eq 3 then color=1
;;
;; Get the right frequency bins of the the multi-color image
;; For grey-scale images, just pick one of the colors
;;
if stokes then begin
   imagedata = a.image[*,*,0,ilam]
endif else begin
   imagedata = a.image[*,*,ilam]
endelse
;;
;; Saturate the image by some factor
;;
if keyword_set(saturate) then begin
   mxx       = max(imagedata)
   imagedata = imagedata < saturate*mxx
endif
;;
;; Log range
;;
if n_elements(lgrange) gt 0 then begin
   if n_elements(log) eq 0 then log=1
   logrange=lgrange
   sz=size(logrange)
   if sz[0] eq 1 then logrange = rebin(logrange,2,nilam)
endif
;;
;; Lin or log
;;
if keyword_set(log) then begin
   ;;
   ;; Log this image
   ;;
   imagedata = alog10(imagedata+1.d-200)
   ;;
   ;; Make sure to have a limited range
   ;;
   ;; For color images, note that the rescaling is not done for each
   ;; color separately, but for all at once.
   ;;
   mn = dblarr(nilam)
   mx = dblarr(nilam)
   if not keyword_set(lgrange) then begin
      ;;
      ;; Log range will be automatically determined
      ;;
      ;; NOTE: Here we will determine the min and max for each
      ;;       color separately. This is different from the old
      ;;       version.
      ;;
      for icol=0,nilam-1 do begin
         mx[icol] = max(imagedata[*,*,icol])
         mn[icol] = min(imagedata[*,*,icol])
      endfor
      if keyword_set(maxlog) gt 0.d0 then begin
         if max(abs(mx-mn)) gt maxlog then begin
            mn = mx - maxlog
         endif
      endif
      ;;
      ;; If the coltune is set 
      ;;
      if n_elements(coltune) ne 0 then begin
         mn = mn - alog10(coltune)
         mx = mx - alog10(coltune)
      endif
   endif else begin
      ;;
      ;; Log range is fixed
      ;;
      mn = transpose(logrange[0,*])
      mx = transpose(logrange[1,*])
   endelse
endif else begin
   ;;
   ;; Linear scaling. Be sure to stick to range.
   ;;
   mn = dblarr(nilam)
   mx = dblarr(nilam)
   for icol=0,nilam-1 do begin
      if not keyword_set(lgrange) then begin
         mx[icol] = max(imagedata[*,*,icol])
         mn[icol] = 0.d0
      endif else begin
         mn[icol] = 10^logrange[0,icol]
         mx[icol] = 10^logrange[1,icol]
      endelse
   endfor
   ;;
   ;; If the coltune is set 
   ;;
   if n_elements(coltune) ne 0 then begin
      mn = mn / coltune
      mx = mx / coltune
   endif
endelse
;;
;; Loop over colors
;;
for icol=0,nilam-1 do begin
   ;;
   ;; Limit to range
   ;;
   imagedata[*,*,icol] = imagedata[*,*,icol] > mn[icol]
   imagedata[*,*,icol] = imagedata[*,*,icol] < mx[icol]
   ;;
   ;; Convert to [0,1]
   ;;
   imagedata[*,*,icol] = ( imagedata[*,*,icol] - mn[icol] ) / ( mx[icol] - mn[icol] )
   ;;
   ;; Convert to [0,255], keeping into account the mincol and maxcol
   ;;
   imagedata[*,*,icol] = mincol + imagedata[*,*,icol] * ( maxcol - mincol )
endfor
;;
;; Interpreting the input
;;
if a.radian then begin
   distscale=1.d0
   distunit='radian'
endif else begin
   distscale=1.d0
   distunit='cm'
   if keyword_set(pc) then begin
      distscale = 3.24073473938d-19
      distunit  = 'pc'
   endif
   if keyword_set(au) then begin
      distscale = 6.68449197861d-14
      distunit  = 'AU'
   endif
endelse
if not keyword_set(nlevels) then nlevels=10
if not keyword_set(xticklen) then xticklen=-0.02
if not keyword_set(yticklen) then yticklen=-0.02
if not keyword_set(xminor) then xminor=-1
if not keyword_set(yminor) then yminor=-1
;
; Cast the image to byte
;
imagedata = byte(imagedata)
;
; Make the x and y coordinates
;
if n_elements(zoom) eq 4 then begin
   x=zoom[0]+(zoom[1]-zoom[0])*(findgen(a.nx)+0.5)/(a.nx*1.d0)
   y=zoom[2]+(zoom[3]-zoom[2])*(findgen(a.ny)+0.5)/(a.ny*1.d0)
endif else begin
   x=((findgen(a.nx)+0.5)/(a.nx*1.d0)-0.5d0)*a.sizepix_x*a.nx*distscale
   y=((findgen(a.ny)+0.5)/(a.ny*1.d0)-0.5d0)*a.sizepix_y*a.ny*distscale
endelse
;
; Now choose between the various output modes
;
if keyword_set(jpg) then begin
   ;;
   ;; Make image for jpeg. The abs(jpg) sets the zoom-factor
   ;; Negative jpg means use red temperature color table
   ;;
   if color eq 1 then begin
      zf      = byte(abs(jpg))
      qq      = bytarr(3,zf*nx,zf*ny)
      qq      = rebin(imagedata,3,zf*nx,zf*ny)
   endif else begin
      if jpg lt 0 then begin
         ;;
         ;; Use red temperature color table
         ;;
         i=indgen(256)
         r=floor(i*256./176.)<255
         g=floor((i-120.)*(255./(255.-120.)))>0
         b=floor((i-190.)*(255./(255.-190.)))>0
         ;;
         ;; cast the image to this color table
         ;;
         zf      = byte(abs(jpg))
         qq      = bytarr(3,zf*nx,zf*ny)
         q       = rebin(imagedata,zf*nx,zf*ny)
         qq[0,*,*] = r[q]
         qq[1,*,*] = g[q]
         qq[2,*,*] = b[q]
         ;;
      endif else begin
         zf      = byte(abs(jpg))
         qq      = bytarr(3,zf*nx,zf*ny)
         qq      = transpose(rebin(imagedata,zf*nx,zf*ny,3),[2,0,1])
      endelse
   endelse
   ;;
   ;; Write the JPG file or return image to keyword jim
   ;;
   if not keyword_set(jim) then begin
      if not keyword_set(filenr) then begin
         filename='image.jpg'
      endif else begin
         filename='image_'+filenr+'.jpg'
      endelse
      write_jpeg,filename,qq,true=1
   endif else begin
      jim=qq
      return
   endelse
   ;;
endif else begin
   ;;
   ;; Plot to the current plotting device
   ;;
   if !d.name eq 'X' then begin
      ;;
      ;; Current plotting device is X, so make image for the X-display
      ;;
      if n_elements(dev) eq 1 then begin
         window,dev,xsize=xsize,ysize=ysize
         if nilam eq 1 then device,decomposed=0
         if nilam eq 3 then device,decomposed=1
      endif
      erase
      if not keyword_set(noimage) then begin
         bimg = byte(congrid(imagedata,position(2),position(3),nilam))
         if nilam eq 1 then tv,bimg,position(0),position(1)
         if nilam eq 3 then tv,bimg,position(0),position(1),true=3
      endif
      plot,x,y,/noerase,position=[position(0),position(1),$
           position(0)+position(2),position(1)+position(3)],/device,/nodata,$
           xminor=xminor,yminor=yminor,xstyle=1,ystyle=1,$
           xticklen=xticklen,yticklen=yticklen,_EXTRA=_extra
      if keyword_set(contour) then begin
         contour,imagedata[*,*,0],x,y,nlevels=nlevels,position=[position(0),position(1),$
                 position(0)+position(2),position(1)+position(3)],/device,$
                 xminor=xminor,yminor=yminor,xstyle=1,ystyle=1,$
                 xticklen=xticklen,yticklen=yticklen,/noerase,_EXTRA=_extra
      endif
   endif else begin
      ;;
      ;; Make image for postscript
      ;;
      if keyword_set(dev) then begin
         device,/color,bits_per_pixel=24,xsize=xsize,ysize=ysize
      endif
      if not keyword_set(noimage) then begin
         bimg=byte(imagedata)
         tv,bimg,position(0),position(1),$
            xsize=position(2),ysize=position(3),/norm
      endif
      plot,x,y,/noerase,position=[position(0),position(1),$
           position(0)+position(2),position(1)+position(3)],/norm,/nodata,$
           xminor=xminor,yminor=yminor,xstyle=1,ystyle=1,$
           xticklen=xticklen,yticklen=yticklen,_EXTRA=_extra
      if keyword_set(contour) then begin
         contour,imagedata[*,*,0],x,y,nlevels=nlevels,position=[position(0),position(1),$
                 position(0)+position(2),position(1)+position(3)],/norm,$
                 xminor=xminor,yminor=yminor,xstyle=1,ystyle=1,$
                 xticklen=xticklen,yticklen=yticklen,/noerase,_EXTRA=_extra
      endif
   endelse
endelse
;;
;; Return the mx and mn if the user wants to have these
;;
if keyword_set(log) then begin
   cmn = 10^mn
   cmx = 10^mx
endif else begin
   cmn = mn
   cmx = mx   
endelse
end

;==========================================================================
;                          MOVIE SUBROUTINES
;==========================================================================

;--------------------------------------------------------------------------
;                             MAKE MOVIE
;--------------------------------------------------------------------------
pro makemovie,_EXTRA=_extra
  makeimage,_EXTRA=_extra
end


;--------------------------------------------------------------------------
;                         MAKE MOVIE FROM FRAMES
;--------------------------------------------------------------------------
pro linkmovie,tf=tf,sf=sf,jpg=jpg,img=img,clean=clean,log=log,$
              maxlog=maxlog,lgrange=lgrange,_EXTRA=_extra
;;
;; Default
;;
if not keyword_set(jpg) and keyword_set(sf) then jpg=sf
if not keyword_set(jpg) then jpg=1
if jpg eq 0 then begin
   print,'ERROR: jpg cannot be 0'
   return
endif
zf=abs(jpg)
if not keyword_set(tf) then tf=1
;;
;; Find the number of movie frames
;;
openr,1,'movie.info'
nim=0
readf,1,nim
close,1
;;
;; Read the first to get the basic data
;;
image = readimage(file='image_0001.out')
nx=image.nx
ny=image.ny
nrfr=image.nrfr
;;
;; Create or clear the movie frames directory, or initialize MPG movie
;;
if not keyword_set(img) then begin
   mpegID = MPEG_OPEN([zf*nx,zf*ny],bitrate=10485720.0)
endif else begin
   spawn,'mkdir MOVIEFRAMES >& /dev/null'
   spawn,'\rm -f MOVIEFRAMES/movieframe_*.jpg >& /dev/null'
   spawn,'\rm -f MOVIEFRAMES/movie.mpg >& /dev/null'
endelse
;;
;; Now make a loop over all images
;;
imarr=dblarr(nx,ny,nim*tf)
for im=1,nim do begin
   ;;
   ;; Read the image
   ;;
   print,'Reading image ',im
   stri = strcompress(string(im),/remove_all)
   if im lt 1000 then stri='0'+stri
   if im lt 100 then stri='0'+stri
   if im lt 10 then stri='0'+stri
   filename='image_'+stri+'.out'
   image = readimage(file=filename)
   if image.stokes then begin
      imarr[*,*,im-1]=image.image[*,*,0]
   endif else begin
      imarr[*,*,im-1]=image.image[*,*]
   endelse
endfor
;;
;; Determine the limits
;;
immax = max(imarr)
immin = min(imarr)
;;
;; If the user does not specify the log of the brightness scale, then
;; do some computations here.
;;
if not keyword_set(lgrange) then begin
   lgrange=alog10([immin,immax])
   if keyword_set(log) and keyword_set(maxlog) then begin
      lgrange[0]=lgrange[1]-maxlog
   endif
   if not keyword_set(log) then lgrange[0]=-99.
endif
;;
;; Now make the movie frames
;;
iframe=1
for im=1,nim do begin
   ;;
   ;; Plot the image to a bitmap array
   ;;
   jim=1
   image.image[*,*]=imarr[*,*,iframe-1]
   plotimage,image,jpg=jpg,jim=jim,log=log,maxlog=maxlog,lgrange=lgrange,$
             _EXTRA=_extra
   ;;
   ;; Check if it is color or not
   ;;
   sz=size(jim)
   if sz(0) eq 3 then begin
      true=1
   endif else begin
      true=0
   endelse
   ;;
   ;; Now write to movie
   ;;
   for ii=1,tf do begin
      strf = strcompress(string(iframe),/remove_all)
      if iframe lt 1000 then strf='0'+strf
      if iframe lt 100 then strf='0'+strf
      if iframe lt 10 then strf='0'+strf
      print,'Image ',stri,', movie frame ',strf
      if not keyword_set(img) then begin
         MPEG_PUT, mpegID, FRAME=iframe, IMAGE=jim
      endif else begin
         write_jpeg,'MOVIEFRAMES/movieframe_'+strf+'.jpg',jim,true=true
      endelse
      iframe=iframe+1
   endfor
endfor
nframes=iframe-1
;;
;; Now finalize the movie
;;
if not keyword_set(img) then begin
   MPEG_SAVE, mpegID, FILENAME='movie.mpg'
   MPEG_CLOSE, mpegID
endif else begin
   cd,current=current
   cd,'MOVIEFRAMES'
   openw,1,'script'
   printf,1,'SIZE movieframe_0001.jpg'
   ;;printf,1,'TEXT 5 -color yellow -font 18 <<EOT'
   ;;printf,1,'   Simple model of '
   ;;printf,1,'turbulent star formation'
   ;;printf,1,'EOT'
   ;;printf,1,'TEXT 2 -font 14 -color white <<EOT'
   ;;printf,1,'          (c) 2008'
   ;;printf,1,'       C.P. Dullemond'
   ;;printf,1,'EOT'
   ;;printf,1,'30x red'
   ;;printf,1,'CLEAR 5 2'
   printf,1,'movieframe_%04d.jpg 1 '+strcompress(string(nframes),/remove_all)
   printf,1,'MPEG PATTERN I'
   close,1
   spawn,'mpp -mpeg script'
   spawn,'\mv movie.mpg ../'
   if keyword_set(clean) then begin
      spawn,'\rm -f movieframe_*.jpg'
      spawn,'\rm -f MPP_frames/*.ppm'
   endif
   cd,current
endelse

spawn,'open movie.mpg >& /dev/null &'

end


;==========================================================================
;                        ROUTINES FOR SPECTRA
;==========================================================================

;-----------------------------------------------------------------
;                       READ A SPECTRUM
;-----------------------------------------------------------------
function readspectrum,ext=ext,file=file,iounit=iounit
nf=0
if not keyword_set(iounit) then funit=1 else funit=iounit
if not keyword_set(iounit) then begin
   ;;
   ;; Read from normal file, so make filename
   ;;
   if not keyword_set(file) then begin
      if n_elements(ext) eq 0 then begin
         file='spectrum.out'
      endif else begin
         file='spectrum_'+strcompress(string(ext),/remove_all)+'.out'
      endelse
   endif
   str=findfile(file)
   if(str(0) ne file) then begin
      print,'Sorry, cannot find ',file
      print,'Presumably radmc3d exited without succes.'
      print,'See above for possible error messages of radmc3d!'
      stop
   endif
   openr,funit,file
endif else begin
   ;;
   ;; Read from radmc3d directly via standard output
   ;; 
   ;; First ask the spectrum from radmc3d
   ;;
   printf,funit,'writespec'
   flush,funit
   ;;
   ;; Now we can read...
   ;;
endelse
;
; Read the image
;
iformat=0
readf,funit,iformat
if iformat ne 1 then begin
   print,'ERROR: File format of ',file,' not recognized.'
   stop
endif
readf,funit,nf
data=dblarr(2,nf)
readf,funit,data
lambda = transpose(data[0,*])
spectrum = transpose(data[1,*])
cc  = 2.9979245800000d10      ; Light speed             [cm/s]
freq = 1d4*cc/lambda
if funit eq 1 then close,1
;
; Return all
; [Use format of structure that originates from the older RADICAL code
; and the RADMC/RAYTRACE codes, but with the lambda entry added]
;
return,{spectrum:spectrum,nfr:nf,freq:freq,lambda:lambda,fr_units:-1,obs:0}
end


;--------------------------------------------------------------------------
;                   CALL RADMC-3D TO MAKE THE SPECTRUM
;--------------------------------------------------------------------------
pro makespectrum,incl=incl,phi=phi,sizecm=sizecm,sizeau=sizeau,$
              sizepc=sizepc,posang=posang,$
              pointcm=pointcm,pointau=pointau,pointpc=pointpc,$
              nostar=nostar,zoomau=zoomau,zoompc=zoompc,iounit=iounit
AU  = 1.496d13       ; Astronomical Unit       [cm]
pc  = 3.0857200d18   ; Parsec                  [cm]
cc  = 2.9979245800000d10      ; Light speed             [cm/s]
;cc  = 2.9979d10      ; Light speed             [cm/s]
if n_elements(incl) eq 0 then incl=0.d0
if n_elements(phi) eq 0 then phi=0.d0
if n_elements(sizecm) ne 0 then sizeau=sizecm/AU
if n_elements(sizepc) ne 0 then sizeau=sizepc*pc/AU
if n_elements(zoompc) ne 0 then zoomau=zoompc*pc/AU
if n_elements(pointcm) ne 0 then pointau=pointcm/AU
if n_elements(pointpc) ne 0 then pointau=pointpc*pc/AU
if n_elements(pointau) ne 0 then begin
   if n_elements(pointau[*,0]) ne 3 then begin
      print,'ERROR in makeimage: point** keywords must be array of 3 numbers (x,y,z)'
      return
   endif
endif
if n_elements(pointau) eq 0 then pointau=[0.,0.,0.]
;;
;; Do some checks 
;;
if n_elements(sizeau) ne 0 and n_elements(zoomau) ne 0 then begin
   print,'ERROR: Cannot specify size and zoom simultaneously'
   return
endif
if n_elements(zoomau) ne 0 and n_elements(zoomau) ne 4 then begin
   print,'ERROR: Zoomau must have 4 elements.'
   return
endif
;;
;; Now make the radmc3d command
;;
;;   command="radmc3d spectrum "+$
command="radmc3d sed"+$
        " incl "+strcompress(string(incl),/remove_all)+$
        " phi "+strcompress(string(phi),/remove_all)+$
        " pointau "+strcompress(string(pointau[0]),/remove_all)+" "+$
        strcompress(string(pointau[1]),/remove_all)+" "+$
        strcompress(string(pointau[2]),/remove_all)
if keyword_set(sizeau) then command = command+" sizeau "+strcompress(string(sizeau),/remove_all)
if keyword_set(zoomau) then command = command+" zoomau "+   $
        strcompress(string(zoomau[0]),/remove_all)+" "+  $
        strcompress(string(zoomau[1]),/remove_all)+" "+  $
        strcompress(string(zoomau[2]),/remove_all)+" "+  $
        strcompress(string(zoomau[3]),/remove_all)
if keyword_set(nostar) then command = command+" nostar" else $
                            command = command+" inclstar"
print,command
;;
;; Are we calling radmc3d or are we communicating with radmc3d as a child process?
;;
if keyword_set(iounit) then begin
   ;;
   ;; Radmc3D is a child process, so communicate with it via the iounit 
   ;; file unit port
   ;;
;   printf,iounit,'spectrum'
   printf,iounit,'sed'
   printf,iounit,'incl'
   printf,iounit,strcompress(string(incl),/remove_all)
   printf,iounit,'phi'
   printf,iounit,strcompress(string(phi),/remove_all)
   printf,iounit,'pointau'
   printf,iounit,strcompress(string(pointau[0]),/remove_all)
   printf,iounit,strcompress(string(pointau[1]),/remove_all)
   printf,iounit,strcompress(string(pointau[2]),/remove_all)
   if keyword_set(sizeau) then begin
      printf,iounit,'sizeau'
      printf,iounit,strcompress(string(sizeau),/remove_all)
   endif
   if keyword_set(zoomau) then begin
      printf,iounit,'zoomau'
      printf,iounit,strcompress(string(zoomau[0]),/remove_all)
      printf,iounit,strcompress(string(zoomau[1]),/remove_all)
      printf,iounit,strcompress(string(zoomau[2]),/remove_all)
      printf,iounit,strcompress(string(zoomau[3]),/remove_all)
   endif
   if keyword_set(nostar) then begin
      printf,iounit,'nostar'
   endif else begin
      printf,iounit,'inclstar'
   endelse
   printf,iounit,'enter'
   flush,iounit
endif else begin
   ;;
   ;; Radmc3D must be called from a shell
   ;;
   ;; Call radmc3d
   ;;
   ;spawn,'../src/'+command
   spawn,command
endelse
;;
end

;-----------------------------------------------------------------
;                       PLOT SPECTRUM 
; NOTE: This routine is nearly identical to the one of RADMC/RAYTRACE
;       and RADICAL.
;
; ARGUMENTS:
;   ev              = 1 --> frequency in electronvolt (default=Hz)
;   kev             = 1 --> frequency in kiloelectronvolt (default=Hz)
;   lnu             = 1 --> L_nu (default nu*L_nu)
;   rpc             = Distance of observer in units of parsec
;                   = 0 (or undef) --> Lum = Total Lum in erg/s
;                     instead of erg / s cm^2
;   lumtot          = 1 --> total luminosity (default=flux at 1 parsec)
;
;-----------------------------------------------------------------
pro plotspectrum,a,ev=ev,kev=kev,hz=hz,micron=micron,$
        lnu=lnu,fnu=fnu,nulnu=nulnu,nufnu=nufnu,rpc=rpc,$
        dpc=dpc,xlg=xlg,oplot=oplot,jy=jy,lsun=lsun,$
        itheta=itheta,ldat=ldat,ylin=ylin,lum=lum,$
        thick=thick,xthick=xthick,ythick=ythick,xlin=xlin,$
        charthick=charthick,charsize=charsize,_extra=_extra
if not keyword_set(itheta) then itheta=0
fr_units = a.fr_units
if keyword_set(micron) then begin
    fr_units = -1
endif
if keyword_set(hz) then begin
    fr_units = 0
endif
if keyword_set(ev) then begin
    fr_units = 1
endif
if keyword_set(kev) then begin
    fr_units = 2
endif
case fr_units of
    -1: begin
        cc  = 2.9979245800000d10 ; Light speed             [cm/s]
        xcoord = 1d4*cc / a.freq(*)
        xtitle = '!4k!X [!4l!Xm]'
    end
    0: begin
        xcoord = a.freq(*)
        xtitle = '!4m!X [Hz]'
    end
    1: begin
        xcoord = 4.13568842841d-15 * a.freq(*)
        xtitle = '!4m!X [eV]'
    end
    2: begin
        xcoord = 4.13568842841d-18 * a.freq(*)
        xtitle = '!4m!X [KeV]'
    end
endcase

if keyword_set(dpc) then rpc=dpc
;
; Plot nuFnu or Fnu (same with Lnu)? And what about Fnu vs Lnu?
;
sed=1
ylum=0
if keyword_set(jy) then begin
   sed=0
endif
if keyword_set(fnu) then begin
   ylum=0
   sed=0
endif
if keyword_set(lnu) then begin
   ylum=1
   sed=0
endif
if keyword_set(nufnu) then begin
   ylum=0
   sed=1
endif
if keyword_set(nulnu) then begin
   ylum=1
   sed=1
endif
if keyword_set(lum) then ylum=1 
if keyword_set(jy) then begin
   ylum=0
   ;sed=0
endif
if keyword_set(lsun) then begin
   ylum=1
   sed=1
endif
;
; If the distance is not given, then for an observed spectrum
; we can only plot the flux, while for a computed spectrum we
; can only plot the luminosity
;
if a.obs eq 0 then begin
   if not keyword_set(rpc) and keyword_set(jy) then begin
      print,"Cannot use Jy units of theoretical spectrum if dpc not known."
      return
   endif
   if not keyword_set(rpc) and ylum eq 0 then begin
      ylum=1
      print,"Do not know distance to source. Plotting Luminosity instead."
   endif
endif else begin
   if not keyword_set(rpc) and ylum eq 1 then begin
      ylum=0
      print,"Do not know distance to source. Plotting Flux instead."
   endif
endelse
;
; Which plot to make? Lum or flux?
;
if ylum eq 0 then begin
   ;
   ; Plot spectrum as flux at a certain distance
   ;
   if a.obs eq 0 then begin
      rpc=1.d0*rpc
      distfact = 1.d0/ (rpc^2)
   endif else begin
      distfact = 1.d0
   endelse
   if not keyword_set(jy) then begin
      if sed eq 0 then begin
         lumfact=1.d0
         ytitle='F!I!4m!X!N  [erg cm!E-2!N Hz!E-1!N s!E-1!N]'
      endif else begin
         lumfact=a.freq(*)
         ytitle='!4m!XF!I!4m!X!N [erg cm!E-2!N s!E-1!N]'
      endelse
   endif else begin
      if sed eq 0 then begin
         lumfact=1d+23
         ytitle='F!I!4m!X!N [Jy]'
      endif else begin
         lumfact=1d+23*a.freq(*)
         ytitle='!4m!XF!I!4m!X!N [JyHz]'
      endelse
   endelse
endif else begin
   ;
   ; Plot spectrum as luminosity
   ;
   if a.obs eq 0 then begin
      distfact = 1.1965280793d38 ; = 4*pi*(1 parsec)^2 = 1.19d38 cm^2
   endif else begin
      rpc=1.d0*rpc
      distfact = rpc^2 * 1.1965280793d38
   endelse
   if sed eq 0 then begin
      lumfact=1.d0
      ytitle='L!I!4m!X!N  [erg Hz!E-1!N s!E-1!N]'
   endif else begin
      if not keyword_set(lsun) then begin
         lumfact=a.freq(*)
         ytitle='!4m!XL!I!4m!X!N [erg s!E-1!N]'
      endif else begin
         lumfact=a.freq(*)*2.5956986d-34
;         ytitle='!4m!XL!I!4m!X!N [Lsun]'
         ytitle='!4m!XL!I!4m!X!N [L!D!9n!X!N]'
      endelse
   endelse
endelse
if not keyword_set(xlg) then begin
    if not keyword_set(xlin) then begin
        xlog=1
    endif else begin
        xlog=0
    endelse
endif else begin
    xcoord=alog10(xcoord)
    xtitle = 'log '+xtitle
    xlog=0
endelse
if not keyword_set(ylin) then begin
    ylog=1
endif else begin
    ylog=0
endelse
if not keyword_set(oplot) then begin
    plot,xcoord,distfact*lumfact*a.spectrum(*,itheta),$
         xtitle=xtitle,ytitle=ytitle,$
         xlog=xlog,ylog=ylog,thick=thick,xthick=xthick,ythick=ythick,$
        charthick=charthick,charsize=charsize,_extra=_extra
endif else begin
    oplot,xcoord,distfact*lumfact*a.spectrum(*,itheta),$
         thick=thick,_extra=_extra
endelse
end


;==========================================================================
;               ROUTINES FOR READING DENSITY AND TEMPERATURE
;==========================================================================

;--------------------------------------------------------------------------
;                    READ THE AMR GRID INFORMATION
;--------------------------------------------------------------------------
function read_amr_grid,basic=basic,mirror=mirror
  close,1
  if n_elements(mirror) eq 0 then mirror=0
  tmp=file_search('amr_grid.inp',count=count1)
  tmp=file_search('amr_grid.uinp',count=count2)
  tmp=file_search('amr_grid.binp',count=count3)
  count = count1+count2+count3
  if count ne 1 then begin
     print,'ERROR: Need 1 and only 1 file amr_grid.*inp...'
     stop
  endif
  nx=0LL
  ny=0LL
  nz=0LL
  nxmax=0LL
  nymax=0LL
  nzmax=0LL
  levelmax=0LL
  nrleafsmax=0LL
  nrbranchesmax=0LL
  nlayers=0LL
  octtree=0LL
  iparent=0LL
  ixyz=0LL
  nxyz=0LL
  nnxyz=0LL
  amrstyle=0LL
  coordsys=0LL
  gridinfo=0LL
  iformat=0LL
  iprecis=0LL
  incd=[0LL,0LL,0LL]
  layer_xi=0.d0
  layer_yi=0.d0
  layer_zi=0.d0
  layer_x=0.d0
  layer_y=0.d0
  layer_z=0.d0
  if count1 gt 0 then begin
     openr,1,'amr_grid.inp'
     ;; For now do things simple
     readf,1,iformat
     readf,1,amrstyle,coordsys,gridinfo
     readf,1,incd
     readf,1,nx,ny,nz
     case amrstyle of
        0: begin
           ;;
           ;; Regular grid
           ;;
        end
        1: begin
           ;;
           ;; Oct-tree style AMR
           ;;
           readf,1,levelmax,nrleafsmax,nrbranchesmax
        end
        10: begin
           ;;
           ;; Layer style AMR
           ;;
           readf,1,levelmax,nlayers
        end
     endcase
     nxmax=nx
     nymax=ny
     nzmax=nz
     xi=dblarr(nx+1)
     yi=dblarr(ny+1)
     zi=dblarr(nz+1)
     readf,1,xi
     readf,1,yi
     readf,1,zi
     if keyword_set(basic) then begin
        ;;
        ;; Only return the basic information of the grid
        ;;
        close,1
        if coordsys lt 100 then begin
           ;;
           ;; Cartesian coordinates
           ;;
           x = 0.5d0 * ( xi[0:nx-1] + xi[1:nx] )
           y = 0.5d0 * ( yi[0:ny-1] + yi[1:ny] )
           z = 0.5d0 * ( zi[0:nz-1] + zi[1:nz] )
           return,{x:x,y:y,z:z,xi:xi,yi:yi,zi:zi,nx:nx,ny:ny,nz:nz,mirror:mirror,$
             coordsys:coordsys,incx:incd[0],incy:incd[1],incz:incd[2],nlayers:nlayers,$
             gridstyle:amrstyle}
        endif else begin
           ;;
           ;; Spherical coordinates
           ;;
           ri = xi
           thetai = yi
           phii = zi
           nr = nx
           ntheta = ny
           nphi = nz
           if keyword_set(mirror) then begin
              thetai = [thetai,!pi-rotate(thetai[0:ny-1],2)]
              ntheta = ntheta * 2 ; I deliberately do not increase ny; only ntheta
           endif
           r = sqrt( ri[0:nr-1] * ri[1:nr] )
           theta = 0.5d0 * ( thetai[0:ntheta-1] + thetai[1:ntheta] )
           phi = 0.5d0 * ( phii[0:nphi-1] + phii[1:nphi] )
           return,{r:r,theta:theta,phi:phi,ri:ri,thetai:thetai,phii:phii,$
             nr:nr,ntheta:ntheta,nphi:nphi,coordsys:coordsys,$
             nx:nx,ny:ny,nz:nz,incx:incd[0],incy:incd[1],incz:incd[2],$
             mirror:mirror,nlayers:nlayers,gridstyle:amrstyle}
        endelse
     endif
     case amrstyle of
        0: begin
           ;;
           ;; Regular grid
           ;;
           ncell    = nx*ny*nz
           ncellinp = ncell
        end
        1: begin
           ;;
           ;; Oct-tree style AMR
           ;;
           octtree  = bytarr(nrbranchesmax)
           readf,1,octtree
           ncell    = n_elements(where(octtree eq 0))
           ncellinp = ncell
        end
        10: begin
           ;;
           ;; Layer style AMR
           ;;
           iparent = intarr(nlayers+1)*0L
           ixyz    = intarr(3,nlayers+1)*0L
           nxyz    = intarr(3,nlayers+1)*0L
           nnxyz   = intarr(3,nlayers+1)*0L
           nnxyz[*,0] = [nx,ny,nz]
           idat    = intarr(7)*0L
           ncell   = nx*ny*nz
           ncellinp= nx*ny*nz
           for i=1,nlayers do begin
              readf,1,idat
              iparent[i] = idat[0]
              ixyz[*,i]  = idat[1:3]
              nxyz[*,i]  = idat[4:6]
              nnxyz[*,i] = idat[4:6]
              if incd[0] eq 1 then nnxyz[0,i]=nnxyz[0,i]*2
              if incd[1] eq 1 then nnxyz[1,i]=nnxyz[1,i]*2
              if incd[2] eq 1 then nnxyz[2,i]=nnxyz[2,i]*2
              ncell      = ncell + nnxyz[0,i]*nnxyz[1,i]*nnxyz[2,i] - $
                                   nxyz[0,i]*nxyz[1,i]*nxyz[2,i]
              ncellinp   = ncellinp + nnxyz[0,i]*nnxyz[1,i]*nnxyz[2,i]
              nxmax = max([nxmax,nnxyz[0,i]])
              nymax = max([nymax,nnxyz[1,i]])
              nzmax = max([nzmax,nnxyz[2,i]])
           endfor
        end
     endcase
     close,1
  endif else begin
     if count2 gt 0 then begin
        openr,1,'amr_grid.uinp',/f77_unformatted
        ;; For now do things simple
        readu,1,iformat
        readu,1,amrstyle
        readu,1,coordsys
        readu,1,gridinfo
        readu,1,incd
        readu,1,nx,ny,nz
        case amrstyle of
           0: begin
              ;;
              ;; Regular grid
              ;;
           end
           1: begin
              ;;
              ;; Oct-tree style AMR
              ;;
              readu,1,levelmax,nrleafsmax,nrbranchesmax
           end
           10: begin
              ;;
              ;; Layer style AMR
              ;;
              readu,1,levelmax,nlayers
           end
        endcase
        nxmax=nx
        nymax=ny
        nzmax=nz
        xi=dblarr(nx+1)
        yi=dblarr(ny+1)
        zi=dblarr(nz+1)
        readu,1,xi
        readu,1,yi
        readu,1,zi
        if keyword_set(basic) then begin
           ;;
           ;; Only return the basic information of the grid
           ;;
           close,1
           if coordsys lt 100 then begin
              ;;
              ;; Cartesian coordinates
              ;;
              x = 0.5d0 * ( xi[0:nx-1] + xi[1:nx] )
              y = 0.5d0 * ( yi[0:ny-1] + yi[1:ny] )
              z = 0.5d0 * ( zi[0:nz-1] + zi[1:nz] )
              return,{x:x,y:y,z:z,xi:xi,yi:yi,zi:zi,nx:nx,ny:ny,nz:nz,mirror:mirror,$
                coordsys:coordsys,incx:incd[0],incy:incd[1],incz:incd[2],nlayers:nlayers,$
                gridstyle:amrstyle}
           endif else begin
              ;;
              ;; Spherical coordinates
              ;;
              ri = xi
              thetai = yi
              phii = zi
              nr = nx
              ntheta = ny
              nphi = nz
              if keyword_set(mirror) then begin
                 thetai = [thetai,!pi-rotate(thetai[0:ny-1],2)]
                 ntheta = ntheta * 2 ; I deliberately do not increase ny; only ntheta
              endif
              r = sqrt( ri[0:nr-1] * ri[1:nr] )
              theta = 0.5d0 * ( thetai[0:ntheta-1] + thetai[1:ntheta] )
              phi = 0.5d0 * ( phii[0:nphi-1] + phii[1:nphi] )
              return,{r:r,theta:theta,phi:phi,ri:ri,thetai:thetai,phii:phii,$
                nr:nr,ntheta:ntheta,nphi:nphi,coordsys:coordsys,$
                nx:nx,ny:ny,nz:nz,incx:incd[0],incy:incd[1],incz:incd[2],$
                mirror:mirror,nlayers:nlayers,gridstyle:amrstyle}
           endelse
        endif
        case amrstyle of
           0: begin
              ;;
              ;; Regular grid
              ;;
              ncell    = nx*ny*nz
              ncellinp = ncell
           end
           1: begin
              ;;
              ;; Oct-tree style AMR
              ;;
              print,'STOPPING: The reading of the oct-tree AMR in unformatted style is not yet ready...'
              stop
              octtree  = bytarr(nrbranchesmax)
              readu,1,octtree
              ncell    = n_elements(where(octtree eq 0))
              ncellinp = ncell
           end
           10: begin
              ;;
              ;; Layer style AMR
              ;;
              iparent = intarr(nlayers+1)*0L
              ixyz    = intarr(3,nlayers+1)*0L
              nxyz    = intarr(3,nlayers+1)*0L
              nnxyz   = intarr(3,nlayers+1)*0L
              nnxyz[*,0] = [nx,ny,nz]
              idat    = intarr(7)*0L
              ncell   = nx*ny*nz
              ncellinp= nx*ny*nz
              for i=1,nlayers do begin
                 readu,1,idat
                 iparent[i] = idat[0]
                 ixyz[*,i]  = idat[1:3]
                 nxyz[*,i]  = idat[4:6]
                 nnxyz[*,i] = idat[4:6]
                 if incd[0] eq 1 then nnxyz[0,i]=nnxyz[0,i]*2
                 if incd[1] eq 1 then nnxyz[1,i]=nnxyz[1,i]*2
                 if incd[2] eq 1 then nnxyz[2,i]=nnxyz[2,i]*2
                 ncell      = ncell + nnxyz[0,i]*nnxyz[1,i]*nnxyz[2,i] - $
                                      nxyz[0,i]*nxyz[1,i]*nxyz[2,i]
                 ncellinp   = ncellinp + nnxyz[0,i]*nnxyz[1,i]*nnxyz[2,i]
                 nxmax = max([nxmax,nnxyz[0,i]])
                 nymax = max([nymax,nnxyz[1,i]])
                 nzmax = max([nzmax,nnxyz[2,i]])
              endfor
           end
        endcase
        close,1
     endif else begin
        openr,1,'amr_grid.binp'
        ;; For now do things simple
        readu,1,iformat
        readu,1,amrstyle
        readu,1,coordsys
        readu,1,gridinfo
        readu,1,incd
        readu,1,nx,ny,nz
        case amrstyle of
           0: begin
              ;;
              ;; Regular grid
              ;;
           end
           1: begin
              ;;
              ;; Oct-tree style AMR
              ;;
              readu,1,levelmax,nrleafsmax,nrbranchesmax
           end
           10: begin
              ;;
              ;; Layer style AMR
              ;;
              readu,1,levelmax,nlayers
           end
        endcase
        nxmax=nx
        nymax=ny
        nzmax=nz
        xi=dblarr(nx+1)
        yi=dblarr(ny+1)
        zi=dblarr(nz+1)
        readu,1,xi
        readu,1,yi
        readu,1,zi
        if keyword_set(basic) then begin
           ;;
           ;; Only return the basic information of the grid
           ;;
           close,1
           if coordsys lt 100 then begin
              ;;
              ;; Cartesian coordinates
              ;;
              x = 0.5d0 * ( xi[0:nx-1] + xi[1:nx] )
              y = 0.5d0 * ( yi[0:ny-1] + yi[1:ny] )
              z = 0.5d0 * ( zi[0:nz-1] + zi[1:nz] )
              return,{x:x,y:y,z:z,xi:xi,yi:yi,zi:zi,nx:nx,ny:ny,nz:nz,mirror:mirror,$
                coordsys:coordsys,incx:incd[0],incy:incd[1],incz:incd[2],nlayers:nlayers,$
                gridstyle:amrstyle}
           endif else begin
              ;;
              ;; Spherical coordinates
              ;;
              ri = xi
              thetai = yi
              phii = zi
              nr = nx
              ntheta = ny
              nphi = nz
              if keyword_set(mirror) then begin
                 thetai = [thetai,!pi-rotate(thetai[0:ny-1],2)]
                 ntheta = ntheta * 2 ; I deliberately do not increase ny; only ntheta
              endif
              r = sqrt( ri[0:nr-1] * ri[1:nr] )
              theta = 0.5d0 * ( thetai[0:ntheta-1] + thetai[1:ntheta] )
              phi = 0.5d0 * ( phii[0:nphi-1] + phii[1:nphi] )
              return,{r:r,theta:theta,phi:phi,ri:ri,thetai:thetai,phii:phii,$
                nr:nr,ntheta:ntheta,nphi:nphi,coordsys:coordsys,$
                nx:nx,ny:ny,nz:nz,incx:incd[0],incy:incd[1],incz:incd[2],$
                mirror:mirror,nlayers:nlayers,gridstyle:amrstyle}
           endelse
        endif
        case amrstyle of
           0: begin
              ;;
              ;; Regular grid
              ;;
              ncell    = nx*ny*nz
              ncellinp = ncell
           end
           1: begin
              ;;
              ;; Oct-tree style AMR
              ;;
              print,'STOPPING: The reading of the oct-tree AMR in unformatted style is not yet ready...'
              stop
              octtree  = bytarr(nrbranchesmax)
              readu,1,octtree
              ncell    = n_elements(where(octtree eq 0))
              ncellinp = ncell
           end
           10: begin
              ;;
              ;; Layer style AMR
              ;;
              iparent = intarr(nlayers+1)*0L
              ixyz    = intarr(3,nlayers+1)*0L
              nxyz    = intarr(3,nlayers+1)*0L
              nnxyz   = intarr(3,nlayers+1)*0L
              nnxyz[*,0] = [nx,ny,nz]
              idat    = intarr(7)*0L
              ncell   = nx*ny*nz
              ncellinp= nx*ny*nz
              for i=1,nlayers do begin
                 readu,1,idat
                 iparent[i] = idat[0]
                 ixyz[*,i]  = idat[1:3]
                 nxyz[*,i]  = idat[4:6]
                 nnxyz[*,i] = idat[4:6]
                 if incd[0] eq 1 then nnxyz[0,i]=nnxyz[0,i]*2
                 if incd[1] eq 1 then nnxyz[1,i]=nnxyz[1,i]*2
                 if incd[2] eq 1 then nnxyz[2,i]=nnxyz[2,i]*2
                 ncell      = ncell + nnxyz[0,i]*nnxyz[1,i]*nnxyz[2,i] - $
                                      nxyz[0,i]*nxyz[1,i]*nxyz[2,i]
                 ncellinp   = ncellinp + nnxyz[0,i]*nnxyz[1,i]*nnxyz[2,i]
                 nxmax = max([nxmax,nnxyz[0,i]])
                 nymax = max([nymax,nnxyz[1,i]])
                 nzmax = max([nzmax,nnxyz[2,i]])
              endfor
           end
        endcase
        close,1
     endelse
  endelse
  ;;
  ;; Some post-processing
  ;;
  if amrstyle eq 10 then begin
     layer_xi = dblarr(nxmax+1,nlayers+1)
     layer_yi = dblarr(nymax+1,nlayers+1)
     layer_zi = dblarr(nzmax+1,nlayers+1)
     layer_x  = dblarr(nxmax,nlayers+1)
     layer_y  = dblarr(nymax,nlayers+1)
     layer_z  = dblarr(nzmax,nlayers+1)
     layer_xi[0:nx,0] = xi
     layer_yi[0:ny,0] = yi
     layer_zi[0:nz,0] = zi
     layer_x[0:nx-1,0]  = 0.5d0 * ( xi[0:nx-1] + xi[1:nx] )
     layer_y[0:ny-1,0]  = 0.5d0 * ( yi[0:ny-1] + yi[1:ny] )
     layer_z[0:nz-1,0]  = 0.5d0 * ( zi[0:nz-1] + zi[1:nz] )
     for i=1,nlayers do begin
        for k=0,nnxyz[0,i],2 do begin
           layer_xi[k,i] = layer_xi[ixyz[0,i]-1+k/2,iparent[i]] 
        endfor
        for k=0,nnxyz[1,i],2 do begin
           layer_yi[k,i] = layer_yi[ixyz[1,i]-1+k/2,iparent[i]] 
        endfor
        for k=0,nnxyz[2,i],2 do begin
           layer_zi[k,i] = layer_zi[ixyz[2,i]-1+k/2,iparent[i]] 
        endfor
        if coordsys lt 100 then begin
           for k=1,nnxyz[0,i]-1,2 do begin
              layer_xi[k,i] = 0.5d0 * ( layer_xi[ixyz[0,i]-1+(k-1)/2,iparent[i]] + $
                                        layer_xi[ixyz[0,i]-1+(k+1)/2,iparent[i]] )
           endfor
        endif else begin
           for k=1,nnxyz[0,i]-1,2 do begin
              layer_xi[k,i] =     sqrt( layer_xi[ixyz[0,i]-1+(k-1)/2,iparent[i]] * $
                                        layer_xi[ixyz[0,i]-1+(k+1)/2,iparent[i]] )
           endfor
        endelse
        for k=1,nnxyz[1,i]-1,2 do begin
           layer_yi[k,i] = 0.5d0 * ( layer_yi[ixyz[1,i]-1+(k-1)/2,iparent[i]] + $
                                     layer_yi[ixyz[1,i]-1+(k+1)/2,iparent[i]] )
        endfor
        for k=1,nnxyz[2,i]-1,2 do begin
           layer_zi[k,i] = 0.5d0 * ( layer_zi[ixyz[2,i]-1+(k-1)/2,iparent[i]] + $
                                     layer_zi[ixyz[2,i]-1+(k+1)/2,iparent[i]] )
        endfor
        if coordsys lt 100 then begin
           for k=0,nnxyz[0,i]-1 do begin
              layer_x[k,i] = 0.5d0 * ( layer_xi[k,i] + layer_xi[k+1,i] )
           endfor
        endif else begin
           for k=0,nnxyz[0,i]-1 do begin
              layer_x[k,i] = sqrt( layer_xi[k,i] * layer_xi[k+1,i] )
           endfor
        endelse
        for k=0,nnxyz[1,i]-1 do begin
           layer_y[k,i] = 0.5d0 * ( layer_yi[k,i] + layer_yi[k+1,i] )
        endfor
        for k=0,nnxyz[2,i]-1 do begin
           layer_z[k,i] = 0.5d0 * ( layer_zi[k,i] + layer_zi[k+1,i] )
        endfor
     endfor
  endif
  if coordsys lt 100 then begin
     ;;
     ;; Cartesian coordinates
     ;;
     x = 0.5d0 * ( xi[0:nx-1] + xi[1:nx] )
     y = 0.5d0 * ( yi[0:ny-1] + yi[1:ny] )
     z = 0.5d0 * ( zi[0:nz-1] + zi[1:nz] )
     return,{x:x,y:y,z:z,xi:xi,yi:yi,zi:zi,nx:nx,ny:ny,nz:nz,$
             nxmax:nxmax,nymax:nymax,nzmax:nzmax,mirror:mirror,$
             ncell:ncell,ncellinp:ncellinp,coordsys:coordsys,octtree:octtree,$
             incx:incd[0],incy:incd[1],incz:incd[2],nlayers:nlayers,$
             ixyz:ixyz,nxyz:nxyz,nnxyz:nnxyz,iparent:iparent,gridstyle:amrstyle,$
             layer_xi:layer_xi,layer_yi:layer_yi,layer_zi:layer_zi,$
             layer_x:layer_x,layer_y:layer_y,layer_z:layer_z}
  endif else begin
     ;;
     ;; Spherical coordinates
     ;;
     ri = xi
     thetai = yi
     phii = zi
     nr = nx
     ntheta = ny
     nphi = nz
     if keyword_set(mirror) then begin
        thetai = [thetai,!pi-rotate(thetai[0:ny-1],2)]
        ntheta = ntheta * 2 ; I deliberately do not increase ny; only ntheta
     endif
     r = sqrt( ri[0:nr-1] * ri[1:nr] )
     theta = 0.5d0 * ( thetai[0:ntheta-1] + thetai[1:ntheta] )
     phi = 0.5d0 * ( phii[0:nphi-1] + phii[1:nphi] )
     return,{r:r,theta:theta,phi:phi,ri:ri,thetai:thetai,phii:phii,$
             nr:nr,ntheta:ntheta,nphi:nphi,ncell:ncell,coordsys:coordsys,$
             nx:nx,ny:ny,nz:nz,incx:incd[0],incy:incd[1],incz:incd[2],$
             octtree:octtree,ncellinp:ncellinp,$
             nxmax:nxmax,nymax:nymax,nzmax:nzmax,mirror:mirror,nlayers:nlayers,$
             ixyz:ixyz,nxyz:nxyz,nnxyz:nnxyz,iparent:iparent,gridstyle:amrstyle,$
             layer_xi:layer_xi,layer_yi:layer_yi,layer_zi:layer_zi,$
             layer_x:layer_x,layer_y:layer_y,layer_z:layer_z}
  endelse
end

;--------------------------------------------------------------------------
;                      LAYER FILL ROUTINE HELPER
;--------------------------------------------------------------------------
pro layers_fill_helper,grid,data,nspec,specinner=specinner
for ilayer=grid.nlayers,1,-1 do begin
   iparent=grid.iparent[ilayer]
   if keyword_set(specinner) then begin
      data[*,grid.ixyz[0,ilayer]-1:grid.ixyz[0,ilayer]+grid.nxyz[0,ilayer]-2,$
        grid.ixyz[1,ilayer]-1:grid.ixyz[1,ilayer]+grid.nxyz[1,ilayer]-2,$
        grid.ixyz[2,ilayer]-1:grid.ixyz[2,ilayer]+grid.nxyz[2,ilayer]-2,$
        iparent] = rebin(data[*,0:grid.nnxyz[0,ilayer]-1,$
                                0:grid.nnxyz[1,ilayer]-1,$
                                0:grid.nnxyz[2,ilayer]-1,$
                                ilayer],nspec,grid.nxyz[0,ilayer],$
                                grid.nxyz[1,ilayer],grid.nxyz[2,ilayer])
   endif else begin
      data[grid.ixyz[0,ilayer]-1:grid.ixyz[0,ilayer]+grid.nxyz[0,ilayer]-2,$
        grid.ixyz[1,ilayer]-1:grid.ixyz[1,ilayer]+grid.nxyz[1,ilayer]-2,$
        grid.ixyz[2,ilayer]-1:grid.ixyz[2,ilayer]+grid.nxyz[2,ilayer]-2,$
        *,iparent] = rebin(data[0:grid.nnxyz[0,ilayer]-1,$
                                0:grid.nnxyz[1,ilayer]-1,$
                                0:grid.nnxyz[2,ilayer]-1,$
                                *,ilayer],grid.nxyz[0,ilayer],$
                                grid.nxyz[1,ilayer],grid.nxyz[2,ilayer],$
                                nspec)
   endelse
endfor
end

;--------------------------------------------------------------------------
;                LAYER FILL ROUTINE (TO FILL THE "HOLES")
;--------------------------------------------------------------------------
pro layers_fill,a,specinner=specinner
if n_elements(a.rho) gt 1 and a.ilayer eq -1 then begin
   data=a.rho
   layers_fill_helper,a.grid,data,a.nspec,specinner=specinner
   a.rho=data
endif
if n_elements(a.temp) gt 1 and a.ilayer eq -1 then begin
   data=a.temp
   layers_fill_helper,a.grid,data,a.nspec,specinner=specinner
   a.temp=data
endif
if n_elements(a.pop) gt 1 and a.ilayer eq -1 then begin
   data=a.pop
   layers_fill_helper,a.grid,data,a.nspec,specinner=specinner
   a.pop=data
endif
end


;--------------------------------------------------------------------------
;                    HELPER ROUTINE FOR READ DATA
;--------------------------------------------------------------------------
pro read_data_helper,grid,data,nspec,ilayer=ilayer,specinner=specinner
case grid.gridstyle of
   0: begin
      ;;
      ;; Regular grid
      ;;
      if keyword_set(specinner) then begin
         data=dblarr(nspec,grid.nx,grid.ny,grid.nz)
         readf,1,data
         if grid.mirror ne 0 and grid.coordsys eq 100 then begin
            data0 = data
            data  = dblarr(nspec,grid.nx,2*grid.ny,grid.nz)
            for iy=0,grid.ny-1 do data[*,*,iy,*] = data0[*,*,iy,*]
            for iy=grid.ny,2*grid.ny-1 do data[*,*,iy,*] = data0[*,*,2*grid.ny-iy-1,*]
         endif
      endif else begin
         data=dblarr(grid.nx,grid.ny,grid.nz,nspec)
         readf,1,data
         if grid.mirror ne 0 and grid.coordsys eq 100 then begin
            data0 = data
            data  = dblarr(grid.nx,2*grid.ny,grid.nz,nspec)
            for iy=0,grid.ny-1 do data[*,iy,*,*] = data0[*,iy,*,*]
            for iy=grid.ny,2*grid.ny-1 do data[*,iy,*,*] = data0[*,2*grid.ny-iy-1,*,*]
         endif
      endelse
   end
   1: begin
      ;;
      ;; Oct-tree AMR grid
      ;;
      ;; Mirror style not explicitly treated here (i.e. if mirror: then
      ;; still only the cells above the equatorial plane are included)
      ;;
      if keyword_set(specinner) then begin
         data=dblarr(nspec,grid.ncellinp)
         readf,1,data
      endif else begin
         data=dblarr(grid.ncellinp,nspec)
         readf,1,data
      endelse
   end
   10: begin
      ;;
      ;; Layer-style AMR grid
      ;;
      ;; Mirror style not explicitly treated here (i.e. if mirror: then
      ;; still only the cells above the equatorial plane are included)
      ;;
      if keyword_set(specinner) then begin
         if n_elements(ilayer) eq 1 then begin
            for ilr=0,ilayer do begin
               data = dblarr(nspec,grid.nnxyz[0,ilr],grid.nnxyz[1,ilr],grid.nnxyz[2,ilr])
               readf,1,data
            endfor
         endif else begin
            data = dblarr(nspec,grid.nxmax,grid.nymax,grid.nzmax,grid.nlayers+1)
            for ilr=0,grid.nlayers do begin
               data0 = dblarr(nspec,grid.nnxyz[0,ilr],grid.nnxyz[1,ilr],grid.nnxyz[2,ilr])
               readf,1,data0
               data[*,0:grid.nnxyz[0,ilr]-1,0:grid.nnxyz[1,ilr]-1,$
                    0:grid.nnxyz[2,ilr]-1,ilr] = data0
            endfor
         endelse
      endif else begin      
         if n_elements(ilayer) eq 1 then begin
            for ilr=0,ilayer do begin
               data = dblarr(grid.nnxyz[0,ilr],grid.nnxyz[1,ilr],grid.nnxyz[2,ilr],nspec)
               readf,1,data
            endfor
         endif else begin
            data = dblarr(grid.nxmax,grid.nymax,grid.nzmax,nspec,grid.nlayers+1)
            for ilr=0,grid.nlayers do begin
               data0 = dblarr(grid.nnxyz[0,ilr],grid.nnxyz[1,ilr],grid.nnxyz[2,ilr],nspec)
               readf,1,data0
               data[0:grid.nnxyz[0,ilr]-1,0:grid.nnxyz[1,ilr]-1,$
                    0:grid.nnxyz[2,ilr]-1,*,ilr] = data0
            endfor
         endelse
      endelse
   end
endcase
end

;--------------------------------------------------------------------------
;                       READ THE GRID-BASED DATA
;--------------------------------------------------------------------------
function read_data,ddens=ddens,dtemp=dtemp,meanint=meanint,mirror=mirror,ilayer=ilayer,fill=fill
  rho = 0.d0
  temp = 0.d0
  jnu = 0.d0
  nspec=0
  nlam=0
  close,1
  grid=read_amr_grid(mirror=mirror)
  if keyword_set(ddens) then begin
     openr,1,'dust_density.inp'
     iformat=0
     readf,1,iformat
     ncell=0L
     readf,1,ncell
     nspec=0L
     readf,1,nspec
     if ncell ne grid.ncellinp then stop
     read_data_helper,grid,rho,nspec,ilayer=ilayer
     close,1
  endif
  if keyword_set(dtemp) then begin
     openr,1,'dust_temperature.dat'
     iformat=0
     readf,1,iformat
     ncell=0L
     readf,1,ncell
     nspec=0L
     readf,1,nspec
     if ncell ne grid.ncellinp then stop
     read_data_helper,grid,temp,nspec,ilayer=ilayer
     close,1
  endif 
  if keyword_set(meanint) then begin
     openr,1,'mean_intensity.out'
     iformat=0
     readf,1,iformat
     if iformat ne 2 then stop
     ncell=0L
     readf,1,ncell
     nlam=0L
     readf,1,nlam
     freq=dblarr(nlam)
     readf,1,freq
     if ncell ne grid.ncellinp then stop
     read_data_helper,grid,jnu,nlam,ilayer=ilayer
     close,1
  endif 
  if n_elements(ilayer) eq 0 then ilayer=-1
  a = {grid:grid,rho:rho,temp:temp,jnu:jnu,ilayer:ilayer,nspec:nspec,nlam:nlam}
  if keyword_set(fill) and grid.gridstyle eq 10 then begin
     layers_fill,a
  endif
  return,a
end

;--------------------------------------------------------------------------
;                       READ LEVEL POPULATIONS
;--------------------------------------------------------------------------
function read_levelpop,molecule,mirror=mirror,ilayer=ilayer,fill=fill,ext=ext,alpha=alpha
  ;;
  ;; Read the populations
  ;; 
  close,1
  grid=read_amr_grid(mirror=mirror)
  filename='levelpop_'+molecule+'.dat'
  openr,1,filename
  iformat=0
  readf,1,iformat
  ncell=0L
  readf,1,ncell
  nlev=0L
  readf,1,nlev
  if ncell ne grid.ncellinp then stop
  ilev=intarr(nlev)
  readf,1,ilev
  pop=dblarr(nlev,ncell)
  read_data_helper,grid,pop,nlev,ilayer=ilayer,/specinner
  close,1
  if n_elements(ilayer) eq 0 then ilayer=-1
  ;;
  ;; If requested, compute also the excitation temperature
  ;;
  if keyword_set(ext) then begin
     filemol='molecule_'+molecule+'.inp'
     openr,1,filemol
     str=''
     readf,1,str
     readf,1,str
     readf,1,str
     mumol=0.d0
     readf,1,mumol
     readf,1,str
     nlv=0
     readf,1,nlv
     if nlv ne nlev then begin
        print,'ERROR in computing excitation temperatures: the file ',filename
        print,'    does not contain all the levels given in ',filemol
        stop
     endif
     readf,1,str
     levdata=dblarr(4,nlev)
     readf,1,levdata
     readf,1,str
     nlin=0
     readf,1,nlin
     readf,1,str
     lindata=dblarr(6,nlin)
     readf,1,lindata
     close,1
     ext=dblarr(nlin,ncell)
     for iline=1,nlin do begin
        iup=lindata[1,iline-1]
        ilo=lindata[2,iline-1] ;; Bugfix 22.02.16
        gup=levdata[2,iup-1]
        glo=levdata[2,ilo-1]
        de=6.6262d-27*lindata[4,iline-1]*1d9
        ext[iline-1,*]=-de/(1.3807d-16*alog(glo*pop[iup-1,*]/(gup*pop[ilo-1,*])))
     endfor
  endif else begin  
     ext=0
  endelse
  ;;
  ;; If requested, compute also the line center opacity alpha
  ;;
  if keyword_set(alpha) then begin
     filename='gas_temperature.inp'
     openr,1,filename
     iformat=0
     readf,1,iformat
     ncell=0L
     readf,1,ncell
     if ncell ne grid.ncellinp then stop
     tgas=dblarr(ncell,1)
     read_data_helper,grid,tgas,1,ilayer=ilayer
     close,1
     filename='microturbulence.inp'
     str=findfile(filename)
     if(str(0) eq filename) then begin
        openr,1,filename
        iformat=0
        readf,1,iformat
        ncell=0L
        readf,1,ncell
        if ncell ne grid.ncellinp then stop
        aturb=dblarr(ncell,1)
        read_data_helper,grid,aturb,1,ilayer=ilayer
        close,1
     endif else begin
        aturb = 0.d0
     endelse
     filemol='molecule_'+molecule+'.inp'
     openr,1,filemol
     str=''
     readf,1,str
     readf,1,str
     readf,1,str
     mumol=0.d0
     readf,1,mumol
     readf,1,str
     nlv=0
     readf,1,nlv
     if nlv ne nlev then begin
        print,'ERROR in computing line center opacities: the file ',filename
        print,'    does not contain all the levels given in ',filemol
        stop
     endif
     readf,1,str
     levdata=dblarr(4,nlev)
     readf,1,levdata
     readf,1,str
     nlin=0
     readf,1,nlin
     readf,1,str
     lindata=dblarr(6,nlin)
     readf,1,lindata
     close,1
     alpha=dblarr(nlin,ncell)
     for iline=1,nlin do begin
        iup=lindata[1,iline-1]
        ilo=lindata[2,iline-1] ;; Bugfix 22.02.16
        gup=levdata[2,iup-1]
        glo=levdata[2,ilo-1]
        nu0=lindata[4,iline-1]*1d9
        aud=lindata[3,iline-1]
        bud=aud*6.7818295d+46/nu0^3
        bdu=bud*(gup/glo)
        mmol=mumol*1.6726000d-24
        atot=sqrt(aturb^2+2*1.3807000d-16*tgas[*]/mmol)
        phi0=1.6913978d+10/(atot*nu0)
        de=6.6262d-27*nu0
        alpha[iline-1,*]=(de/12.566371d0)*(pop[ilo-1,*]*bdu-pop[iup-1,*]*bud)*phi0
     endfor
  endif else begin  
     alpha=0.d0
     tgas=0.d0
     aturb=0.d0
  endelse
  a = {grid:grid,pop:pop,ilayer:ilayer,nlev:nlev,ilev:ilev,molecule:molecule,ext:ext,alpha:alpha}
  if keyword_set(fill) and grid.gridstyle eq 10 then begin
     layers_fill,a,/specinner
  endif
  return,a
end

;--------------------------------------------------------------------------
;                       READ THE GRID-BASED DATA
;--------------------------------------------------------------------------
function readdata,_extra=_extra
return,read_data(_extra=_extra)
end

;==========================================================================
;               ROUTINES FOR READING SUBBOX DATA SETS
;==========================================================================


;-----------------------------------------------------------------
;                          MAKE SUBBOX 
;-----------------------------------------------------------------
pro makesubbox,variable,nxyz=nxyz,box=box,size=size,pos=pos,  $
               phi1=phi1,theta=theta,phi2=phi2,               $
               radmc3d=radmc3d,au=au,pc=pc,iounit=iounit
  AAU  = 1.496d13                ; Astronomical Unit       [cm]
  ppc  = 3.0857200d18            ; Parsec                  [cm]
  ccc  = 2.9979245800000d10      ; Light speed             [cm/s]
;  ccc  = 2.9979d10               ; Light speed             [cm/s]
  if not keyword_set(radmc3d) then begin
     tmp=file_search('./radmc3d',count=count1)
     if count1 eq 0 then begin
        radmc3d='radmc3d'
     endif else begin
        radmc3d='./radmc3d'
     endelse
  endif
  nx = 64
  ny = 64
  nz = 64
  if keyword_set(nxyz) then begin
     if n_elements(nxyz) eq 1 then begin
        nx=nxyz[0]
        ny=nxyz[0]
        nz=nxyz[0]
     endif else begin
        nx=nxyz[0]
        ny=nxyz[1]
        nz=nxyz[2]
     endelse
  endif
  if not keyword_set(phi1) then phi1=0.d0
  if not keyword_set(theta) then theta=0.d0
  if not keyword_set(phi2) then phi2=0.d0
  if keyword_set(box) then begin
     x0=box[0]*1.d0
     x1=box[1]*1.d0
     y0=box[2]*1.d0
     y1=box[3]*1.d0
     z0=box[4]*1.d0
     z1=box[5]*1.d0
  endif
  if keyword_set(size) then begin
     if n_elements(size) eq 1 then begin
        x0=-size[0]/2.d0
        x1=size[0]/2.d0
        y0=-size[0]/2.d0
        y1=size[0]/2.d0
        z0=-size[0]/2.d0
        z1=size[0]/2.d0
     endif else begin
        x0=-size[0]/2.d0
        x1=size[0]/2.d0
        y0=-size[1]/2.d0
        y1=size[1]/2.d0
        z0=-size[2]/2.d0
        z1=size[2]/2.d0
     endelse
     if keyword_set(pos) then begin
        x0=x0+pos[0]
        x1=x1+pos[0]
        y0=y0+pos[1]
        y1=y1+pos[1]
        z0=z0+pos[2]
        z1=z1+pos[2]
     endif
  endif
  if keyword_set(au) and n_elements(x0) ge 1 then begin
     x0=x0*aau
     x1=x1*aau
     y0=y0*aau
     y1=y1*aau
     z0=z0*aau
     z1=z1*aau
  endif
  if keyword_set(pc) and n_elements(x0) ge 1 then begin
     x0=x0*ppc
     x1=x1*ppc
     y0=y0*ppc
     y1=y1*ppc
     z0=z0*ppc
     z1=z1*ppc
  endif
  command=radmc3d+" subbox_"+variable+" subbox_nxyz "+  $
          strcompress(string(nx),/remove_all)+" "+      $
          strcompress(string(ny),/remove_all)+" "+      $
          strcompress(string(nz),/remove_all)+" "
  if n_elements(x0) ge 1 then begin
     command=command+" subbox_xyz01 "+                  $
          strcompress(string(x0),/remove_all)+" "+      $
          strcompress(string(x1),/remove_all)+" "+      $
          strcompress(string(y0),/remove_all)+" "+      $
          strcompress(string(y1),/remove_all)+" "+      $
          strcompress(string(z0),/remove_all)+" "+      $
          strcompress(string(z1),/remove_all)+" "
  endif
  command=command+" subbox_phi1 "+                      $
          strcompress(string(phi1),/remove_all)+" "+    $
          " subbox_theta "+                             $
          strcompress(string(theta),/remove_all)+" "+   $
          " subbox_phi2 "+                              $
          strcompress(string(phi2),/remove_all)
  print,command
  ;;
  ;; Are we calling radmc3d or are we communicating with radmc3d as a child process?
  ;;
  if keyword_set(iounit) then begin
     printf,iounit,"subbox_"+variable
     printf,iounit,"subbox_nxyz"
     printf,iounit,strcompress(string(nx),/remove_all)
     printf,iounit,strcompress(string(ny),/remove_all)
     printf,iounit,strcompress(string(nz),/remove_all)
     if n_elements(x0) ge 1 then begin
        printf,iounit,"subbox_xyz01"
        printf,iounit,strcompress(string(x0),/remove_all)
        printf,iounit,strcompress(string(x1),/remove_all)
        printf,iounit,strcompress(string(y0),/remove_all)
        printf,iounit,strcompress(string(y1),/remove_all)
        printf,iounit,strcompress(string(z0),/remove_all)
        printf,iounit,strcompress(string(z1),/remove_all)
     endif
     printf,iounit,"subbox_phi1"
     printf,iounit,strcompress(string(phi1),/remove_all)
     printf,iounit,"subbox_theta"
     printf,iounit,strcompress(string(theta),/remove_all)
     printf,iounit,"subbox_phi2"
     printf,iounit,strcompress(string(phi2),/remove_all)
     printf,iounit,'enter'
     flush,iounit
  endif else begin
     ;;
     ;; Radmc3D must be called from a shell
     ;;
     spawn,command
  endelse
end


;-----------------------------------------------------------------
;                          READ SUBBOX 
;-----------------------------------------------------------------
function readsubbox,filename
openr,1,filename
iformat=0
readf,1,iformat
if iformat gt 2 then begin
   print,'ERROR reading subbox. Format number not known.'
   stop
endif
nx=0
ny=0
nz=0
nv=0
readf,1,nx,ny,nz,nv
x0=0.d0
x1=0.d0
y0=0.d0
y1=0.d0
z0=0.d0
z1=0.d0
readf,1,x0,x1,y0,y1,z0,z1
phi1=0.d0
theta=0.d0
phi2=0.d0
readf,1,phi1,theta,phi2
id=intarr(nv)
if iformat gt 1 then begin
   readf,1,id
endif
data=dblarr(nx,ny,nz,nv)
readf,1,data
close,1
return,{nx:nx,ny:ny,nz:nz,x0:x0,x1:x1,y0:y0,y1:y1,z0:z0,z1:z1, $
        phi1:phi1,theta:theta,phi2:phi2,id:id,data:data}
end


;-----------------------------------------------------------------
;                     MAKE SUBBOX AND READ
;-----------------------------------------------------------------
function subbox,variable,_extra=_extra
  makesubbox,variable,_extra=_extra
  filename=variable+"_subbox.out"
  q=readsubbox(filename)
  return,q
end


;==========================================================================
;            ROUTINES FOR READING AND PLOTTING THE OPACITIES
;
; Note: For "historic reasons" the kappa_abs and the kappa_scat were written
; as cabs and csca (which are actually officially the notations of the
; cross sec for 1 particle). From now on they will be called kappa_abs
; and kappa_scat, but the cabs and csca will still be there for backward
; compatibility reasons.
;==========================================================================


;------------------------------------------------------;
; READ THE DUST OPACITY FILES                          ;
;------------------------------------------------------;
function readopac,spec=spec,used=used
close,1
if not keyword_set(spec) then spec=1
filename_kap = 'dustkappa_'+strcompress(string(spec),/remove_all)+'.inp'
kapfile = file_test(filename_kap)
filename_scatmat = 'dustkapscatmat_'+strcompress(string(spec),/remove_all)+'.inp'
scatmatfile = file_test(filename_scatmat)
if not kapfile and not scatmatfile then begin
   print,'ERROR: Neither ',filename_kap,' nor ',filename_scatmat,' present.'
   stop
endif
if kapfile and scatmatfile then begin
   print,'ERROR: Both ',filename_kap,' nor ',filename_scatmat,' present.'
   stop
endif
if kapfile then begin
   filename = filename_kap
   if keyword_set(used) then filename=filename+'.used'
   print,"Reading ",filename
   openr,1,filename
   icnt=0
   str=' '
   flag=0
   while flag eq 0 do begin
      readf,1,str
      i=0
      while strmid(str,i,1) eq ' ' do i=i+1
      if((strmid(str,i,1) eq '#') or (strmid(str,i,1) eq ';') or (strmid(str,i,1) eq '!')) then begin
         icnt=icnt+1
      endif else begin
         flag = 1
      endelse
   endwhile
   close,1
   openr,1,filename
   for i=1,icnt do readf,1,str
   iformat=0
   readf,1,iformat
   nf=0
   readf,1,nf
   case iformat of
      1: ncol=2
      2: ncol=3
      3: ncol=4
   endcase
   data=dblarr(ncol,nf)
   readf,1,data
   close,1
   lambda    = dblarr(nf)
   kappa_abs = dblarr(nf)
   kappa_scat= dblarr(nf)
   g         = dblarr(nf)
   lambda[*] = data[0,*]
   kappa_abs[*] = data[1,*]
   if iformat ge 2 then kappa_scat[*] = data[2,*]
   if iformat ge 3 then g[*] = data[3,*]
   angle=0
   zmat=0
   pmat=0
   mu=0
   na=0
endif else begin
   filename = filename_scatmat
   print,"Reading ",filename
   openr,1,filename
   icnt=0
   str=' '
   flag=0
   while flag eq 0 do begin
      readf,1,str
      i=0
      while strmid(str,i,1) eq ' ' do i=i+1
      if((strmid(str,i,1) eq '#') or (strmid(str,i,1) eq ';') or (strmid(str,i,1) eq '!')) then begin
         icnt=icnt+1
      endif else begin
         flag = 1
      endelse
   endwhile
   close,1
   openr,1,filename
   for i=1,icnt do readf,1,str
   iformat=0
   readf,1,iformat
   nf=0
   readf,1,nf
   na=0
   readf,1,na
   data=dblarr(4,nf)
   readf,1,data
   lambda     = dblarr(nf)
   kappa_abs  = dblarr(nf)
   kappa_scat = dblarr(nf)
   g          = dblarr(nf)
   lambda[*]     = data[0,*]
   kappa_abs[*]  = data[1,*]
   kappa_scat[*] = data[2,*]
   g[*]          = data[3,*]
   angle = dblarr(na)
   readf,1,angle
   mu=cos(angle*!dpi/180.d0)
   mu[0]=1.d0
   mu[na-1]=-1.d0
   ;;
   ;; Read the scattering matrix
   ;;
   zmat = dblarr(6,na,nf)
   readf,1,zmat
   close,1
   ;;
   ;; Now normalize to obtain the phase matrix
   ;;
   pmat = dblarr(6,na,nf)
   for inu=0,nf-1 do begin
      pmat[*,*,inu] = 4*!dpi * zmat[*,*,inu] / kappa_scat[inu]
   endfor
endelse
;
; Some postprocessing
;
cc  = 2.9979245800000d10        ; Light speed             [cm/s]
freq = 1d4*cc/lambda
;
return,{nf:nf,ns:1,na:na,freq:freq,wave:lambda,lambda:lambda,cabs:kappa_abs,csca:kappa_scat,$
        kappa_abs:kappa_abs,kappa_scat:kappa_scat,g:g,nrt:1,trange:0.d0,$
        angle:angle,mu:mu,zmat:zmat,pmat:pmat}
end

;------------------------------------------------------;
; PLOT THE DUST OPACITIES                              ;
;------------------------------------------------------;
pro plotopac,a,abs=abs,scat=scat,is=is,oplot=oplot,$
             ylin=ylin,irt=irt,mult=mult,xlin=xlin,$
             _extra=_extra
ipl=0
ylog=1
xlog=1
if not keyword_set(irt) then irt=0
if irt ge a.nrt then begin
   print,'irt too large'
   return
endif
xtitle='!4k!X [!4l!Xm]'
ytitle='!4j!X!D!4k!X!N [cm!U2!Ng!U-1!N]'
if keyword_set(xlin) then xlog=0
if keyword_set(ylin) then ylog=0
if not keyword_set(oplot) then oplot=0
if not keyword_set(mult) then mult=1.d0
if keyword_set(abs) then begin
    if n_elements(is) gt 0 then begin
        if is ge 0 then begin
	    if oplot eq 0 then begin
                plot,a.wave,mult*a.kappa_abs(*,is,irt),ylog=ylog,xlog=xlog,xtitle=xtitle,ytitle=ytitle,_extra=_extra
            endif else begin
                oplot,a.wave,mult*a.kappa_abs(*,is,irt),_extra=_extra
            endelse
        endif else begin
            dum = mult*a.kappa_abs(*,0,irt)
            for is=1,a.ns-1 do begin
                dum = dum + mult*a.kappa_abs(*,is,irt)
            endfor
	    if oplot eq 0 then begin
               plot,a.wave,dum,ylog=ylog,xlog=xlog,_extra=_extra
	    endif else begin
               oplot,a.wave,dum,_extra=_extra
            endelse
        endelse
    endif else begin
        if oplot eq 0 then begin
            plot,a.wave,mult*a.kappa_abs(*,0,irt),ylog=ylog,xlog=xlog,xtitle=xtitle,ytitle=ytitle,_extra=_extra
        endif else begin
            oplot,a.wave,mult*a.kappa_abs(*,0,irt),_extra=_extra
	endelse
        for is=1,a.ns-1 do begin
            oplot,a.wave,mult*a.kappa_abs(*,is,irt)
        endfor
    endelse
    ipl = 1
endif
if keyword_set(scat) then begin
    if n_elements(is) gt 0 then begin
        if is ge 0 then begin
            if oplot eq 0 then begin
                plot,a.wave,mult*a.kappa_scat(*,is,irt),ylog=ylog,xlog=xlog,xtitle=xtitle,ytitle=ytitle,_extra=_extra
            endif else begin
                oplot,a.wave,mult*a.kappa_scat(*,is,irt),_extra=_extra
	    endelse
        endif else begin
            dum = mult*a.kappa_scat(*,0,irt)
            for is=1,a.ns-1 do begin
                dum = dum + mult*a.kappa_scat(*,is,irt)
            endfor
            if oplot eq 0 then begin
                plot,a.wave,dum,ylog=ylog,xlog=xlog,xtitle=xtitle,ytitle=ytitle,_extra=_extra
            endif else begin
                oplot,a.wave,dum,_extra=_extra
	    endelse
        endelse
    endif else begin
        if oplot eq 0 then begin
            plot,a.wave,mult*a.kappa_scat(*,0,irt),ylog=ylog,xlog=xlog,xtitle=xtitle,ytitle=ytitle,_extra=_extra
        endif else begin
            oplot,a.wave,mult*a.kappa_scat(*,0,irt),_extra=_extra
        endelse
        for is=1,a.ns-1 do begin
            oplot,a.wave,mult*a.kappa_scat(*,is,irt)
        endfor
    endelse
    ipl = 1
endif
if ipl eq 0 then begin
    dum = mult*a.kappa_scat(*,0,irt) + mult*a.kappa_abs(*,0,irt)
    for is=1,a.ns-1 do begin
        dum = dum + mult*a.kappa_scat(*,is,irt) + mult*a.kappa_abs(*,is,irt)
    endfor
    if oplot eq 0 then begin
       plot,a.wave,dum,ylog=ylog,xlog=xlog,xtitle=xtitle,ytitle=ytitle,_extra=_extra
    endif else begin
       oplot,a.wave,dum,_extra=_extra
    endelse
endif
end

;-----------------------------------------------------------------
;                 FIND A FREQUENCY IN A TABLE
;-----------------------------------------------------------------
function find_freq,freq,nu,eps=eps
nf=n_elements(freq)
idx=-100
for inu=0,nf-2 do begin
   if (nu ge freq(inu) and nu le freq(inu+1)) or $
     (nu ge freq(inu+1) and nu le freq(inu)) then begin
      idx = inu
      eps = (nu-freq[inu]) / (freq[inu+1]-freq[inu])
   endif
endfor
if idx lt 0 then begin
   print,"ERROR: nu out of range"
   stop
endif
return,idx
end

;-----------------------------------------------------------------
;                 FIND OPACITY AT CERTAIN WAVELENGTH
;-----------------------------------------------------------------
function findopac,o,lambda=lambda,freq=freq,abs=abs,sca=sca
cc  = 2.9979245800000d10      ; Light speed             [cm/s]
if n_elements(abs) eq 0 and n_elements(sca) eq 0 then begin
   abs=1
   sca=1
endif else begin
   if n_elements(abs) eq 0 then abs=0
   if n_elements(sca) eq 0 then sca=0
endelse
if keyword_set(lambda) then freq=1d4*cc/lambda
if not keyword_set(freq) then stop
if freq ge max(o.freq) then return,0.d0
if freq le min(o.freq) then return,0.d0
idx=find_freq(o.freq,freq)
eps=(freq-o.freq[idx])/(o.freq[idx+1]-o.freq[idx])
kappa=0.d0
if abs gt 0 then kappa = kappa + $
   (1.d0-eps)*o.kappa_abs(idx)+eps*o.kappa_abs(idx+1)
if sca gt 0 then kappa = kappa + $
   (1.d0-eps)*o.kappa_scat(idx)+eps*o.kappa_scat(idx+1)
return,kappa
end

;==========================================================================
;                 ROUTINE FOR READING THE TAU SURFACE 3-D DATA
;==========================================================================

;--------------------------------------------------------------------------
;                        READ THE 3-D TAU SURFACE DATA
;--------------------------------------------------------------------------
function read_tausurf_3d
  openr,1,'tausurface_3d.out'
  iformat=0
  readf,1,iformat
  if iformat ne 1 then stop
  nx=0
  ny=0
  readf,1,nx,ny
  nf=0
  readf,1,nf
  lambda_micron=dblarr(nf)
  readf,1,lambda_micron
  data=dblarr(3,nx,ny)
  readf,1,data
  x=dblarr(nx,ny)
  y=dblarr(nx,ny)
  z=dblarr(nx,ny)
  x[*,*] = data[0,*,*]
  y[*,*] = data[1,*,*]
  z[*,*] = data[2,*,*]
  close,1
  return,{nx:nx,ny:ny,nf:nf,lambda_micron:lambda_micron,x:x,y:y,z:z}
end

;==========================================================================
;          ROUTINE FOR READING THE SUBPIXELING DIAGNOSTICS DATA
;==========================================================================

;-----------------------------------------------------------------
;           READING THE SUBPIXELING DIAGNOSTICS DATA
;-----------------------------------------------------------------
function readsubpixels
  readcol,'subpixeling_diagnostics.out',x,y,dx,dy
  return,{x:x,y:y,dx:dx,dy:dy}
end

;-----------------------------------------------------------------
;                       PLOT SUBPIXELS
;-----------------------------------------------------------------
pro plotsubpixels,a,level=level,minlev=minlev,maxlev=maxlev,_extra=_extra
  if n_elements(level) eq 1 then begin
     minlev=level
     maxlev=level
  endif else begin
     if not keyword_set(minlev) then minlev=0
     if not keyword_set(maxlev) then maxlev=10000
  endelse
  np=n_elements(a.x)
  dxmax=max(a.dx)
  thelevels=intarr(np)
  i=0LL
  for i=0LL,np-1LL do begin
     qx=floor(1.1*dxmax/a.dx[i])
     if qx lt 1 then stop
     l = 0
     while qx gt 1 do begin
        qx=qx/2
        l=l+1
     endwhile
     thelevels[i]=l
  endfor
  print,'Maximum level = ',max(thelevels)
  ii=where(thelevels ge minlev and thelevels le maxlev)
  if ii[0] ge 0 then begin
     x=a.x[ii]
     y=a.y[ii]
     plot,x,y,_extra=_extra
  endif else begin
     print,'Nothing to plot. Check minlev and maxlev keywords'
  endelse
end

;==========================================================================
;                      SOME MISCELLANEOUS ROUTINES
;==========================================================================

;------------------------------------------------------------------------
;                       FUNCTION: PLANCK
;------------------------------------------------------------------------
function bplanck,nnu,TT
nu=nnu*1.d0
T=TT*1.d0
;cc=2.9979d10
cc  = 2.9979245800000d10      ; Light speed             [cm/s]
hh=6.6262d-27
kk=1.3807e-16
n=n_elements(nu)
bpl=dblarr(n)
if TT eq 0.d0 then return,dblarr(n)
for i=0,n-1 do begin
  x=hh*nu[i]/(kk*T)
  ;print,x
  if x gt 100.0 then b=0.d0
  if x gt 1.d-3 then b=(2.d0*hh*nu[i]^3/cc^2)/(exp(x)-1.d0) $
    else b=2.0*nu[i]^2*kk*T/cc^2
  bpl[i]=b
endfor
return,bpl
end


;-----------------------------------------------------------------
;               CALCULATE THE PLANCK MEAN OPACITY
;-----------------------------------------------------------------
function planckmean_abs,opac,temp
if n_elements(temp) ne 1 then begin
   print,'ERROR: Can only compute planck mean for a single temperature'
   stop
endif
nu   = opac.freq
kabs = opac.kappa_abs
bpl  = bplanck(nu,temp)
dum1 = 0.d0
dum2 = 0.d0
for inu=0,n_elements(nu)-2 do begin
   dum1 = dum1 + 0.5d0*(kabs[inu]*bpl[inu]+kabs[inu+1]*bpl[inu+1]) * $
          (nu[inu+1]-nu[inu])
   dum2 = dum2 + 0.5d0*(bpl[inu]+bpl[inu+1]) * $
          (nu[inu+1]-nu[inu])
endfor
return,dum1/dum2
end


;-----------------------------------------------------------------
;            CALCULATE OPTICALLY THIN DUST TEMPERATURE
;-----------------------------------------------------------------
function optthindusttemp,opac,dist,tstar,rstar,tol=tol
if n_elements(tol) ne 1 then tol = 1d-6
nu     = opac.freq
tdust0 = tstar * sqrt(rstar/(2*dist))
kps    = planckmean_abs(opac,tstar)
kpd    = planckmean_abs(opac,tdust0)
iter   = 0
err    = 1d99
while err gt tol do begin
   kpdold = kpd
   tdust  = tdust0 * (kps/kpd)^0.25
   kpd    = planckmean_abs(opac,tdust)
   err    = abs(kpd/(kpdold+1d-99)-1.d0)
   iter   = iter+1
   if iter gt 100 then err=-1
endwhile
if err lt 0.d0 then print,'Temperature not converged'
return,tdust
end
