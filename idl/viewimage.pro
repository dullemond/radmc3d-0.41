@readradmc.pro
;========================================================================
;      WIDGET ROUTINE FOR VIEWING IMAGES FROM RADMC-3D CALCULATION
;========================================================================

;----------------------------------------------------------
; IMAGE VIEWER
;----------------------------------------------------------
pro viewimage,color=color,nochild=nochild,fluxunits=fluxunits,$
              pc=pc,au=au,small=small,verti=verti,lines=lines,$
              radmc3d=radmc3d,local=local,usertrans=usertrans
;
; Check if a local radmc3d executable is present
;
if not keyword_set(radmc3d) then begin
   z=findfile('radmc3d',count=radmchere)
   if radmchere ge 1 then begin
      radmc3d = './radmc3d'
   endif else begin
      radmc3d = 'radmc3d'
   endelse
endif
;
; First start radmc-3d as a child process
; NOTE: All radmc3d output is redirected to radmc3d.out, so we don't see this 
;       on the screen. Also it may take a while before all the data have
;       been read into radmc3d. 
;
if not keyword_set(nochild) then begin
   print,'Starting RADMC-3D... This may take some time, because of data reading...'
   print,'In the mean time the widget will already pop up, but will only activate'
   print,'when RADMC-3D is ready...'
   spawn,'nice '+radmc3d+' child',unit=iounit
   print,'RADMC-3D is now running as a child of IDL, communicating via unit ',$
          strcompress(string(iounit),/remove_all)
endif else begin
   iounit = 0
endelse
print,'NOTE: viewimage has the options:'
print,'  /pc      ---> Use parsec as scale'
print,'  /au      ---> Use AU as scale'
print,'  /color   ---> Make RGB color images'
print,'  /verti   ---> Put the controllers below the image'
print,'  /small   ---> Use a smaller drawing window for small displays'
print,'  /nochild ---> Do not spawn radmc3d as a child, but call each time'
print,'  /lines   ---> Include line options'
print,'  /local   ---> Allow (in addition to normal view) local observer view'
;
; Max values
;
im_nx_max        = 400
im_ny_max        = 400
npmax            = im_nx_max*2
;
; Default values
;
if not keyword_set(small) then begin
    draw_xsize       = 500
    draw_ysize       = 500
    margins          = [70,70,30,30]
endif else begin
    draw_xsize       = 400
    draw_ysize       = 400
    margins          = [50,50,20,20]
endelse
incldcatchbutton = 0
ctinit           = 3
xstxt            = 50
im_nx_default    = 100
im_incl_default  = 60.e0
im_phi_default   = 30.e0
im_orient_default= 0e0
nphot_default    = 100000LL
im_maxlog_default= 6.d0
lambda_default_r = 1000.d0
lambda_default_g = 10.d0
lambda_default_b = 1.d0
if n_elements(fluxunits) eq 0 then fluxunits = 1
if not keyword_set(color) then color=0
if keyword_set(background) then begin
    bg_temp_default  = 3.d0
    bg_dil_default   = 1.d0
endif else begin
    bg_temp_default  = 0.d0
    bg_dil_default   = 0.d0
    background       = 0
    wtx_bgtmp        = 0
    wtx_bgdil        = 0
endelse
if not keyword_set(lines) then lines=0
plot_log_default = 0
rotpts  = [0,0,0,0]*1d0
zoombox = [-1,-1,-1,-1]*1d99
zoomcur = [-1,-1,-1,-1]*1d99
zoomcen = [-1,-1]*1d99
if not keyword_set(au) then au=0
if not keyword_set(pc) then pc=0
if not keyword_set(bg) then bg=0
if not keyword_set(local) then local=0
;
; Find wavelength file and read it
;
nlam = 0
z=findfile('wavelength_micron.inp',count=count)
if(count gt 0) then begin
   openr,1,'wavelength_micron.inp'
   readf,1,nlam
   lambda = dblarr(nlam)
   readf,1,lambda
   close,1
   if nlam le 1 then stop
endif else begin
   z=findfile('frequency.inp',count=count)
   if(count gt 0) then begin
      openr,1,'frequency.inp'
      readf,1,nlam
      lambda = dblarr(nlam)
      readf,1,lambda
      close,1
      cc  = 2.9979245800000d10  ; Light speed             [cm/s]
      lambda = 1d4*cc/lambda
      if nlam le 1 then stop
   endif else begin
      stop
   endelse
endelse
;
; Find the index for the default wavelength
;
if lambda_default_r le max(lambda) and lambda_default_r ge min(lambda) then begin
   ilam_r = 1
   for inu=1,nlam-1 do begin
      if (lambda_default_r-lambda[inu-1])*(lambda_default_r-lambda[inu]) le 0 then begin
         ilam_r=inu
      endif
   endfor
endif else begin
   if(lambda[1] gt lambda[0]) then begin
      if lambda_default_r gt max(lambda) then ilam_r = nlam else ilam_r = 1
   endif else begin
      if lambda_default_r gt max(lambda) then ilam_r = 1 else ilam_r = nlam
   endelse
endelse
if lambda_default_g le max(lambda) and lambda_default_g ge min(lambda) then begin
   ilam_g = 1
   for inu=1,nlam-1 do begin
      if (lambda_default_g-lambda[inu-1])*(lambda_default_g-lambda[inu]) le 0 then begin
         ilam_g=inu
      endif
   endfor
endif else begin
   if(lambda[1] gt lambda[0]) then begin
      if lambda_default_g gt max(lambda) then ilam_g = nlam else ilam_g = 1
   endif else begin
      if lambda_default_g gt max(lambda) then ilam_g = 1 else ilam_g = nlam
   endelse
endelse
if lambda_default_b le max(lambda) and lambda_default_b ge min(lambda) then begin
   ilam_b = 1
   for inu=1,nlam-1 do begin
      if (lambda_default_b-lambda[inu-1])*(lambda_default_b-lambda[inu]) le 0 then begin
         ilam_b=inu
      endif
   endfor
endif else begin
   if(lambda[1] gt lambda[0]) then begin
      if lambda_default_b gt max(lambda) then ilam_b = nlam else ilam_b = 1
   endif else begin
      if lambda_default_b gt max(lambda) then ilam_b = 1 else ilam_b = nlam
   endelse
endelse
;
; Find radius file and read it
;
if radmchere ne 0 then begin
   ;;
   ;; If using a local radmc3d executable, then probably the user has
   ;; made a model internally. So ask radmc3d to write the grid file.
   ;;
   radmc3d_write_gridfile,iounit=iounit
endif
grid = read_amr_grid(/basic)
if grid.coordsys lt 100 then begin
   im_xsize_default = 3*max([abs(grid.xi),abs(grid.yi),abs(grid.zi)])
endif else begin
   im_xsize_default = 2*1.1*max(grid.ri)
endelse
if keyword_set(pc) then begin
    im_xsize_default = im_xsize_default / 3.08572d18
    if keyword_set(au) then stop ; Not both pc and au allowed
endif
if keyword_set(au) then begin
    im_xsize_default = im_xsize_default / 1.496d13
endif
;
; Create the backdrop pixmap
;
window,/free,/pixmap,xsize=draw_xsize,ysize=draw_ysize
bitmap=!d.window
;
; Creating the widgets
;
if not keyword_set(verti) then begin
    wbase      = WIDGET_BASE(column=2)
endif else begin
    wbase      = WIDGET_BASE(column=1)
endelse
wbase1     = WIDGET_BASE(wbase,column=1)
wbase2     = WIDGET_BASE(wbase,column=1)
bbase1     = WIDGET_BASE(wbase1,column=2)
but_quit   = WIDGET_BUTTON(bbase1,value='Quit Viewer',$
             xsize=255,uvalue='but_quit')            
but_wrtim  = WIDGET_BUTTON(bbase1,value='Write Image',$ 
             xsize=255,uvalue='but_wrtim') 
wdr_main   = WIDGET_DRAW(wbase1,xsize=draw_xsize,$
             ysize=draw_ysize,frame=10,uvalue='wdraw',/button_events)
;
; Create widgets for choice between lin and log
;
imbase2    = WIDGET_BASE(wbase2,column=2)
iibase2    = WIDGET_BASE(imbase2,column=1,frame=1)
ibase2     = WIDGET_BASE(iibase2,column=3)
wbg_mouse  = WIDGET_DROPLIST(ibase2,value=['rotate','zoom a','zoom b','zoom c'],$
             uvalue='wbg_mouse')
if not keyword_set(color) then begin
   wbg_con = CW_BGROUP(ibase2,['contour'],$
             column=1,uvalue='wbg_con',/nonexclusive)
endif else begin
   wbg_con = 0
endelse
wbg_lilg   = CW_BGROUP(ibase2,['lin'],$
             column=1,uvalue='wbg_lilg',/nonexclusive)
wbg_star   = CW_BGROUP(ibase2,['star'],$
             column=1,uvalue='wbg_star',set_value=1,/nonexclusive)
if not keyword_set(local) then begin
   wbg_prev = CW_BGROUP(ibase2,['preview'],$
              column=1,uvalue='wbg_prev',/nonexclusive,set_value=1)
   wbg_local = 0
endif else begin
   wbg_local = CW_BGROUP(ibase2,['loc obs'],$
               column=1,uvalue='wbg_local',/nonexclusive,set_value=0)
   wbg_prev  = 0
endelse
wbg_fixseed = CW_BGROUP(ibase2,['fxseed'],$
             column=1,uvalue='wbg_fixseed',/nonexclusive,set_value=0)
wbg_tau    = CW_BGROUP(ibase2,['tau'],$
             column=1,uvalue='wbg_tau',/nonexclusive,set_value=0)
wbg_second = CW_BGROUP(ibase2,['2nd-O'],$
               column=1,uvalue='wbg_second',/nonexclusive,set_value=0)
if not keyword_set(color) then begin
   wsl_ct  = WIDGET_SLIDER(ibase2,minimum=0,maximum=15,$  
             xsize=50,uvalue='wsl_ct',value=0,/drag)
endif else wsl_ct=0
z=findfile('lines.inp',count=linesinp)
if linesinp gt 0 then begin
   wbg_dcatch  = CW_BGROUP(ibase2,['dcatch'],$
               column=1,uvalue='wbg_dcatch',/nonexclusive,set_value=0)
   incldcatchbutton = 1
endif else begin
   wbg_dcatch=0
   incldcatchbutton = 0
endelse
imbase2b   = WIDGET_BASE(imbase2,column=1)
ibase3     = WIDGET_BASE(imbase2b,column=4,frame=1)
wlb_maxlog = WIDGET_LABEL(ibase3,xsize=xstxt,$
             ysize=30,uvalue='wlb_maxlog',$
             value='MaxLog')
wlb_nrcont = WIDGET_LABEL(ibase3,xsize=xstxt,$
             ysize=30,uvalue='wlb_nrcont',$
             value='Nr Cont')
wtx_maxlog = WIDGET_TEXT(ibase3,xsize=3,$    
             ysize=1,frame=0,uvalue='wtx_maxlog',$
             /editable)                       
wtx_nrcont = WIDGET_TEXT(ibase3,xsize=3,$    
             ysize=1,frame=0,uvalue='wtx_nrcont',$
             /editable)                       
wlb_satur  = WIDGET_LABEL(ibase3,xsize=xstxt,$
             ysize=30,uvalue='wlb_satur',$
             value='Saturate')
;wlb_posang = WIDGET_LABEL(ibase3,xsize=xstxt,$
;             ysize=30,uvalue='wlb_posang',$
;             value='PosAng')
wlb_nphot  = WIDGET_LABEL(ibase3,xsize=xstxt,$
             ysize=30,uvalue='wlb_nphot',$
             value='Nphot')
wtx_satur  = WIDGET_TEXT(ibase3,xsize=6,$    
             ysize=1,frame=0,uvalue='wtx_satur',$
             /editable)                       
;wtx_posang = WIDGET_TEXT(ibase3,xsize=6,$    
;             ysize=1,frame=0,uvalue='wtx_posang',$
;             /editable)                       
wtx_nphot  = WIDGET_TEXT(ibase3,xsize=8,$    
             ysize=1,frame=0,uvalue='wtx_nphot',$
             /editable)                       
but_spec   = WIDGET_BUTTON(imbase2b,value='Render Spectrum',$
             xsize=200,ysize=35,uvalue='but_spec')            
;
; Create positioning control widgets and render + unzoom button
;
tmbase1    = WIDGET_BASE(wbase2,column=3)
tbase0     = WIDGET_BASE(tmbase1,column=1,frame=1)
but_makeim = WIDGET_BUTTON(tbase0,value="Render Image",$
             xsize=150,ysize=40,uvalue='but_makeim')
but_unzoom = WIDGET_BUTTON(tbase0,value="Unzoom",$
             xsize=150,ysize=22,uvalue='but_unzoom')
tbase1     = WIDGET_BASE(tmbase1,column=2,frame=1)
wlb_nx     = WIDGET_LABEL(tbase1,xsize=xstxt,$
             ysize=30,uvalue='wlb_nx',$
             value='Npix')
wlb_sx     = WIDGET_LABEL(tbase1,xsize=xstxt,$
             ysize=30,uvalue='wlb_sx',$
             value='Size')
wtx_nx     = WIDGET_TEXT(tbase1,xsize=10,$    
             ysize=1,frame=0,uvalue='wtx_nx',$
             /editable)                       
wtx_sx     = WIDGET_TEXT(tbase1,xsize=10,$    
             ysize=1,frame=0,uvalue='wtx_sx',$
             /editable)                       
tbase4     = WIDGET_BASE(tmbase1,column=2,frame=1)
wlb_incl   = WIDGET_LABEL(tbase4,xsize=90,$
             ysize=30,uvalue='wlb_incl',$
             value='Inclination')
wlb_phi    = WIDGET_LABEL(tbase4,xsize=xstxt,$
             ysize=30,uvalue='wlb_phi',$
             value='Phi')
wtx_incl   = WIDGET_TEXT(tbase4,xsize=10,$    
             ysize=1,frame=0,uvalue='wtx_incl',$
             /editable)
wtx_phi    = WIDGET_TEXT(tbase4,xsize=10,$    
             ysize=1,frame=0,uvalue='wtx_phi',$
             /editable)
if keyword_set(local) then begin
   ptng_base = WIDGET_BASE(wbase2,column=2)
   ptng_b1 = WIDGET_BASE(ptng_base,row=2,frame=1)
   ptng_b2 = WIDGET_BASE(ptng_base,row=2,frame=1)
   wlb_pnt   = WIDGET_LABEL(ptng_b1,xsize=90,$
               ysize=30,uvalue='wlb_pnt',$
               value='Pointing:')
   wtx_pntx  = WIDGET_TEXT(ptng_b1,xsize=10,$    
               ysize=1,frame=0,uvalue='wtx_pntx',$
               /editable)
   wtx_pnty  = WIDGET_TEXT(ptng_b1,xsize=10,$    
               ysize=1,frame=0,uvalue='wtx_pnty',$
               /editable)
   wtx_pntz  = WIDGET_TEXT(ptng_b1,xsize=10,$    
               ysize=1,frame=0,uvalue='wtx_pntz',$
               /editable)
   wlb_obs   = WIDGET_LABEL(ptng_b1,xsize=90,$
               ysize=30,uvalue='wlb_obs',$
               value='Observer Pos:')
   wtx_obsx  = WIDGET_TEXT(ptng_b1,xsize=10,$    
               ysize=1,frame=0,uvalue='wtx_obsx',$
               /editable)
   wtx_obsy  = WIDGET_TEXT(ptng_b1,xsize=10,$    
               ysize=1,frame=0,uvalue='wtx_obsy',$
               /editable)
   wtx_obsz  = WIDGET_TEXT(ptng_b1,xsize=10,$    
               ysize=1,frame=0,uvalue='wtx_obsz',$
               /editable)
   wlb_view  = WIDGET_LABEL(ptng_b2,xsize=xstxt,$
               ysize=30,uvalue='wlb_view',$
               value='ViewAng:')
   wtx_view  = WIDGET_TEXT(ptng_b2,xsize=10,$    
               ysize=1,frame=0,uvalue='wtx_view',$
               /editable)
   wbg_relobs= CW_BGROUP(ptng_b2,['relobs'],$
             column=1,uvalue='wbg_relobs',set_value=0,/nonexclusive)
   wtx_scale = WIDGET_TEXT(ptng_b2,xsize=4,$    
               ysize=1,frame=0,uvalue='wtx_scale',$
               /editable)
endif else begin
   ptng_base = 0
   wtx_pntx  = 0
   wtx_pnty  = 0
   wtx_pntz  = 0
   wtx_obsx  = 0
   wtx_obsy  = 0
   wtx_obsz  = 0
   wtx_view  = 0
   wbg_relobs = 0
   wtx_scale = 0
endelse
;
; Slider for choosing wavelength
;
if keyword_set(color) then begin
   tmbase2   = WIDGET_BASE(wbase2,row=3)
   wslsx     = 300
   wbg_scl   = CW_BGROUP(ibase2,['absolute scale'],$
               column=1,uvalue='wbg_scl',/nonexclusive)
endif else begin
   tmbase2   = WIDGET_BASE(wbase2,row=1)
   wslsx     = 400
   wbg_scl   = 0
endelse
wsl_lam_r = WIDGET_SLIDER(tmbase2,minimum=1,maximum=nlam,$  
             xsize=wslsx,uvalue='wsl_lam_r',value=1,/drag)
wtx_lam_r = WIDGET_TEXT(tmbase2,xsize=11,/editable,$
             ysize=1,uvalue='wtx_lam_r')
if keyword_set(color) then begin
   wtx_coltune_r = WIDGET_TEXT(tmbase2,xsize=10,$
               /editable,ysize=1,uvalue='wtx_coltune_r')
   wsl_lam_g = WIDGET_SLIDER(tmbase2,minimum=1,maximum=nlam,$  
               xsize=wslsx,uvalue='wsl_lam_g',value=1,/drag)
   wtx_lam_g = WIDGET_TEXT(tmbase2,xsize=11,/editable,$
               ysize=1,uvalue='wtx_lam_g')
   wtx_coltune_g = WIDGET_TEXT(tmbase2,xsize=10,$
               /editable,ysize=1,uvalue='wtx_coltune_g')
   wsl_lam_b = WIDGET_SLIDER(tmbase2,minimum=1,maximum=nlam,$  
               xsize=wslsx,uvalue='wsl_lam_b',value=1,/drag)
   wtx_lam_b = WIDGET_TEXT(tmbase2,xsize=11,/editable,$
               ysize=1,uvalue='wtx_lam_b')
   wtx_coltune_b = WIDGET_TEXT(tmbase2,xsize=10,$
               /editable,ysize=1,uvalue='wtx_coltune_b')
endif else begin
   wsl_lam_g = 0
   wsl_lam_b = 0
   wtx_lam_g = 0
   wtx_lam_b = 0
   wtx_coltune_r = 0
   wtx_coltune_g = 0
   wtx_coltune_b = 0
endelse
if keyword_set(lines) then begin
   linebase   = WIDGET_BASE(wbase2,row=1)
   wlb_imol   = WIDGET_LABEL(linebase,xsize=40,$
                ysize=30,uvalue='wlb_imol',$
                value='iMolec')
   wtx_imol   = WIDGET_TEXT(linebase,xsize=4,/editable,$
                ysize=1,uvalue='wtx_imol')
   wlb_iline  = WIDGET_LABEL(linebase,xsize=40,$
                ysize=30,uvalue='wlb_iline',$
                value='iLine')
   wtx_iline  = WIDGET_TEXT(linebase,xsize=4,/editable,$
                ysize=1,uvalue='wtx_iline')
   wlb_vkms   = WIDGET_LABEL(linebase,xsize=70,$
                ysize=30,uvalue='wlb_vkms',$
                value='v [km/s]')
   wtx_vkms   = WIDGET_TEXT(linebase,xsize=11,/editable,$
                ysize=1,uvalue='wtx_vkms')
endif else begin
   wlb_imol   = 0
   wtx_imol   = 0 
   wlb_iline  = 0 
   wtx_iline  = 0 
   wlb_vkms   = 0 
   wtx_vkms   = 0 
endelse
vars=0
vals=0
nvars=0
if keyword_set(usertrans) then begin
   ;;
   ;; Find the transfer.inp.default file, which is essential for this
   ;; routine.
   ;;
   z=findfile('transfer.inp.default',count=cnt)
   if cnt eq 0 then begin
      vars=0
      vals=0
      nvars=0
   endif else begin
      openr,1,'transfer.inp.default'
      vars = [' ']
      vals = [' ']
      while not eof(1) do begin
         str=''
         readf,1,str
         ps=strpos(str,'=')
         if ps gt 0 then begin
            name = strcompress(strmid(str,0,ps-1),/remove_all)
            value = strcompress(strmid(str,ps+1,100),/remove_all)
            vars=[vars,name]
            vals=[vals,value]
         endif
      endwhile
      close,1
      nvars=n_elements(vars)
      vars=vars[1:nvars-1]
      vals=vals[1:nvars-1]
      nvars=nvars-1
   endelse
endif 
if nvars gt 0 then begin
   utbase   = WIDGET_BASE(wbase2,column=1)
   for ivar=0,nvars-1 do begin
      name       = strcompress(string(ivar),/remove_all)
      utvbase    = WIDGET_BASE(utbase,column=2)
      wlb_utr    = WIDGET_LABEL(utvbase,xsize=250,$
             ysize=30,uvalue='wlb_uvar_'+name,$
             value=vars[ivar])
      wtx_utr    = WIDGET_TEXT(utvbase,xsize=6,$    
             ysize=1,frame=0,uvalue='wtx_uvar_'+name,$
             /editable)
      if ivar eq 0 then begin
         wtx_utrans = [wtx_utr]
      endif else begin
         wtx_utrans = [wtx_utrans,wtx_utr]
      endelse
   endfor
endif else begin
   wtx_utrans = 0
endelse
;
; Put things on the screen
;
WIDGET_CONTROL,wbase,/realize
WIDGET_CONTROL,wdr_main,get_value=window
;
; Handle the coloring
;
if not keyword_set(color) then begin
   set_plot,'x'
   device,decomposed=0
   loadct,ctinit
endif else begin
   set_plot,'x'
   device,decomposed=1
endelse
;
; Set up a local image buffer
; 
if keyword_set(color) then nrcol=3 else nrcol=1
cmn   = dblarr(nrcol)
cmx   = dblarr(nrcol)
img   = dblarr(npmax,npmax,nrcol)
xx    = dblarr(npmax)
flux  = dblarr(nrcol)
lamb  = dblarr(nrcol)
radian= (1 eq 1)
image = {nx:0,ny:0,nrfr:nrcol,sizepix_x:0.d0,sizepix_y:0.d0,$
         image:img,flux:flux,x:xx,y:xx,lambda:lamb,radian:radian,stokes:(1 eq 0)}
;
; Now make a big structure with all info
;
state={iounit:iounit,$                  ; The unit number for data transfer to radmc3d
       window:window,$                  ; Main window of the widget
       wdr_main:wdr_main,$              ; Draw window
       wtx_sx:wtx_sx,$                  ; Entry: Size of image in model
       wtx_nx:wtx_nx,$                  ; Entry: Nr of pixels in a dimension
       wtx_maxlog:wtx_maxlog,$          ; Entry: Maximum range in 10log(brightness)
       wtx_incl:wtx_incl,$              ; Entry: Inclination (from z-axis)
       wtx_phi:wtx_phi,$                ; Entry: Rotation around z-axis
       ;wtx_posang:wtx_posang,$          ; Entry: Rotation in image plane
       wtx_nphot:wtx_nphot,$            ; Entry: Nr of photons in MC
       wsl_lam_r:wsl_lam_r,$            ; Slider: Wavelength
       wsl_lam_g:wsl_lam_g,$            ; Slider: Wavelength
       wsl_lam_b:wsl_lam_b,$            ; Slider: Wavelength
       wtx_lam_r:wtx_lam_r,$            ; Output: Wavelength in microns
       wtx_lam_g:wtx_lam_g,$            ; Output: Wavelength in microns
       wtx_lam_b:wtx_lam_b,$            ; Output: Wavelength in microns
       wtx_coltune_r:wtx_coltune_r,$    ; Input: Color tuning
       wtx_coltune_g:wtx_coltune_g,$    ; Input: Color tuning
       wtx_coltune_b:wtx_coltune_b,$    ; Input: Color tuning
       wtx_satur:wtx_satur,$            ; Input: Saturation of image
       wbg_lilg:wbg_lilg,$              ; Button: Lin or Log image
       wbg_second:wbg_second,$          ; Button: Second order integration
       wbg_dcatch:wbg_dcatch,$          ; Button: Doppler catching algor active
       wbg_prev:wbg_prev,$              ; Button: Preview image or flux-consistent image
       wbg_fixseed:wbg_fixseed,$        ; Button: Always reset the initial random seed
       wbg_tau:wbg_tau,$                ; Button: Image = optical depth 
       wbg_star:wbg_star,$              ; Button: Include/exclude star(s)
       wbg_con:wbg_con,$                ; Button: For contours
       wbg_mouse:wbg_mouse,$            ; Droplist: The function of the mouse in the draw panel
       wbg_scl:wbg_scl,$                ; Button: For auto or manual color scale
       mousefunction:0,$                ; For use with the mouse action droplist
       nlam:nlam,$                      ; Number of wavelengths
       zoom:0,$                         ; Zoom mode active?
       zoomcur:zoomcur,$                ; Current geometry of image
       zoombox:zoombox,$                ; Coordinates of zoom box
       zoomcen:zoomcen,$                ; Coordinates of zoom box
       lambda:lambda,$                  ; Wavelength array in micron
       color:color,$                    ; Switch for RGB color mode
       npmax:npmax,$                    ; Max size of image buffer
       image:image,$                    ; Local image buffer
       bitmap:bitmap,$                  ; The backdrop bitmap
       draw_xsize:draw_xsize,$          ; Draw window geometry x
       draw_ysize:draw_ysize,$          ; Draw window geometry y
       pc:pc,$                          ; Use parsec as units
       au:au,$                          ; Use AU as units
       margins:margins,$                ; Margins in pixels for display
       rotpts:rotpts,$                  ; Memory of the rotation points
       fluxunits:fluxunits,$            ; Unit of flux used 
       wtx_nrcont:wtx_nrcont,$          ; Widget for number of contours
       wsl_ct:wsl_ct,$                  ; Color table slider
       local:local,$                    ; Local observer feature available?
       wbg_local:wbg_local,$            ; Button for local observer
       wtx_pntx:wtx_pntx,$              ; Text entries for local observer mode
       wtx_pnty:wtx_pnty,$              ; Text entries for local observer mode
       wtx_pntz:wtx_pntz,$              ; Text entries for local observer mode
       wtx_obsx:wtx_obsx,$              ; Text entries for local observer mode
       wtx_obsy:wtx_obsy,$              ; Text entries for local observer mode
       wtx_obsz:wtx_obsz,$              ; Text entries for local observer mode
       wtx_view:wtx_view,$              ; Viewing width
       wbg_relobs:wbg_relobs,$          ; Do observer position relative or absolute wrt pointing
       wtx_scale:wtx_scale,$            ; Optional scaling factor for obs pos
       radmc3d:radmc3d,$                ; The command used to execute radmc3d
       vars:vars,$                      ; For user-defined transfer function
       vals:vals,$                      ; For user-defined transfer function
       nvars:nvars,$                    ; For user-defined transfer function
       wtx_utrans:wtx_utrans,$          ; For user-defined transfer function
       cmn:cmn,cmx:cmx,$                ; Memory of color range of last plot
       incldcatchbutton:incldcatchbutton,$; Include the doppler catching button?
       lines:lines,$                    ; Include lines or not
       doline:0,$                       ; Internal flag
       wtx_imol:wtx_imol,$              ; Widget selecting molecule
       wtx_iline:wtx_iline,$            ; Widget selecting line
       wtx_vkms:wtx_vkms}               ; Widget selecting velocity

;
; Give this structure to the main widget 
;
WIDGET_CONTROL,wbase,set_uvalue=state
;
; Set default values in widgets
;
WIDGET_CONTROL,wtx_nx,set_value=string(format='(I4)',im_nx_default)
WIDGET_CONTROL,wtx_sx,set_value=string(format='(E9.2)',im_xsize_default)
WIDGET_CONTROL,wtx_incl,set_value=string(format='(F9.2)',im_incl_default)
WIDGET_CONTROL,wtx_phi,set_value=string(format='(F9.2)',im_phi_default)
;WIDGET_CONTROL,wtx_posang,set_value=strcompress(string(format='(F9.2)',im_orient_default),/remove_all)
WIDGET_CONTROL,wtx_nphot,set_value=strcompress(string(format='(I9)',nphot_default),/remove_all)
WIDGET_CONTROL,wtx_maxlog,set_value=string(format='(F3.1)',im_maxlog_default)
WIDGET_CONTROL,wtx_satur,set_value='1.00'
WIDGET_CONTROL,wtx_nrcont,set_value='20'
if not keyword_set(color) then WIDGET_CONTROL,wsl_ct,set_value=ctinit
WIDGET_CONTROL,wsl_lam_r,set_value=ilam_r
WIDGET_CONTROL,wtx_lam_r,set_value=strcompress(string(format='(F9.2)',state.lambda[ilam_r-1]),/remove_all)+' mu'
if keyword_set(color) then begin
   WIDGET_CONTROL,wsl_lam_g,set_value=ilam_g
   WIDGET_CONTROL,wsl_lam_b,set_value=ilam_b
   WIDGET_CONTROL,wtx_lam_g,set_value=strcompress(string(format='(F9.2)',state.lambda[ilam_g-1]),/remove_all)+' mu'
   WIDGET_CONTROL,wtx_lam_b,set_value=strcompress(string(format='(F9.2)',state.lambda[ilam_b-1]),/remove_all)+' mu'
   WIDGET_CONTROL,wtx_coltune_r,set_value='1.00'
   WIDGET_CONTROL,wtx_coltune_g,set_value='1.00'
   WIDGET_CONTROL,wtx_coltune_b,set_value='1.00'
endif
if keyword_set(lines) then begin
   WIDGET_CONTROL,wtx_imol,set_value='1'
   WIDGET_CONTROL,wtx_iline,set_value='1'
   WIDGET_CONTROL,wtx_vkms,set_value='0.0'
endif
if keyword_set(local) then begin
   WIDGET_CONTROL,wtx_pntx,set_value='0.E+00'
   WIDGET_CONTROL,wtx_pnty,set_value='0.E+00'
   WIDGET_CONTROL,wtx_pntz,set_value='0.E+00'
   WIDGET_CONTROL,wtx_obsx,set_value='0.E+00'
   WIDGET_CONTROL,wtx_obsy,set_value= $
      strcompress(string(-2*im_xsize_default,format='(E10.3)'),/remove_all)
   WIDGET_CONTROL,wtx_obsz,set_value='0.E+00'
   WIDGET_CONTROL,wtx_view,set_value='0.8'
endif
if nvars gt 0 then begin
   for ivar=0,nvars-1 do begin
      WIDGET_CONTROL,wtx_utrans[ivar],set_value=strcompress(string(vals[ivar]),/remove_all)
   endfor
endif
;
; If lines included, make linelist
;
if state.lines eq 0 then begin
   linelist = 0
endif else begin
   linelist = 1
endelse
;
; Make a first plot
;
view_the_image,state,/makeimage,linelist=linelist
;
; Give this structure to the main widget 
;
WIDGET_CONTROL,wbase,set_uvalue=state
;
; Now activate the widget
;
xmanager,'viewimage',wbase,/NO_BLOCK
;
end


;----------------------------------------------------------
; WIDGET EVENT HANDLER
;----------------------------------------------------------
PRO viewimage_event,sevent
WIDGET_CONTROL,sevent.id,get_uvalue=uvalue
WIDGET_CONTROL,sevent.top,get_uvalue=state,/no_copy
notrecog=0
case uvalue of
   'but_quit': begin
      WIDGET_CONTROL,sevent.top,/destroy
      !except=1
      if state.iounit ne 0 then begin
         printf,state.iounit,'quit'
         close,state.iounit
         free_lun,state.iounit
      endif
      return
   end
   'but_makeim': begin
      view_the_image,state,/makeimage
   end
   'but_spec': begin
      if state.local eq 1 then begin
         WIDGET_CONTROL,state.wbg_local,get_value=txt
         txt=txt(0)
         localobs = txt + 0
      endif else begin
         localobs = 0
      endelse
      if localobs eq 0 then begin
         view_the_spectrum,state,/makespec
      endif
   end
   'wtx_incl': begin
      view_the_image,state,/makeimage
   end
   'wtx_phi': begin
      view_the_image,state,/makeimage
   end
;   'wtx_posang': begin
;      view_the_image,state,/makeimage
;   end
   'wtx_nphot': begin
      view_the_image,state,/makeimage
   end
   'wtx_nx': begin
      view_the_image,state,/makeimage
   end
   'wtx_sx': begin
      view_the_image,state,/makeimage
   end
   'but_unzoom': begin
      view_the_image,state,/makeimage,/unzoom
   end
   'but_wrtim': begin
      view_the_image,state,/ps
   end
   'wsl_lam_r': begin
      WIDGET_CONTROL,state.wsl_lam_r,get_value=txt
      txt=strcompress(string(txt(0)),/remove_all)
      ilam_r = txt + 0
      WIDGET_CONTROL,state.wtx_lam_r,set_value=strcompress(string(format='(F9.2)',state.lambda[ilam_r-1]),/remove_all)+' mu'
   end
   'wsl_lam_g': begin
      WIDGET_CONTROL,state.wsl_lam_g,get_value=txt
      txt=strcompress(string(txt(0)),/remove_all)
      ilam_g = txt + 0
      WIDGET_CONTROL,state.wtx_lam_g,set_value=strcompress(string(format='(F9.2)',state.lambda[ilam_g-1]),/remove_all)+' mu'
   end
   'wsl_lam_b': begin
      WIDGET_CONTROL,state.wsl_lam_b,get_value=txt
      txt=strcompress(string(txt(0)),/remove_all)
      ilam_b = txt + 0
      WIDGET_CONTROL,state.wtx_lam_b,set_value=strcompress(string(format='(F9.2)',state.lambda[ilam_b-1]),/remove_all)+' mu'
   end
   'wtx_lam_r': begin
      view_the_image,state,/makeimage
   end
   'wtx_lam_g': begin
      view_the_image,state,/makeimage
   end
   'wtx_lam_b': begin
      view_the_image,state,/makeimage
   end
   'wtx_imol': begin
      state.doline = 1
      view_the_image,state,/makeimage
   end
   'wtx_iline': begin
      state.doline = 1
      view_the_image,state,/makeimage
   end
   'wtx_vkms': begin
      state.doline = 1
      view_the_image,state,/makeimage
   end
   'wtx_maxlog': begin
      view_the_image,state
   end
   'wbg_lilg': begin
      view_the_image,state
   end
   'wbg_mouse': begin
      state.zoombox[*] = -1d99
      state.zoomcen[*] = -1d99
      state.mousefunction = sevent.index
      view_the_image,state
   end
   'wtx_angp': begin
      view_the_image,state
   end
   'wbg_star': begin
      view_the_image,state,/makeimage
   end
   'wbg_second': begin
      view_the_image,state,/makeimage
   end
   'wbg_dcatch': begin
      view_the_image,state,/makeimage
   end
   'wbg_prev': begin
      view_the_image,state,/makeimage
   end
   'wbg_fixseed': begin
      view_the_image,state,/makeimage
   end
   'wbg_tau': begin
      view_the_image,state,/makeimage
   end
   'wtx_satur': begin
      view_the_image,state
   end
   'wtx_coltune_r': begin
      view_the_image,state
   end
   'wtx_coltune_g': begin
      view_the_image,state
   end
   'wtx_coltune_b': begin
      view_the_image,state
   end
   'wbg_scl': begin
      WIDGET_CONTROL,state.wbg_scl,get_value=txt
      colorscale = txt + 0
      if colorscale eq 1 then begin
         txt = strcompress(string(state.cmx[0],format='(E10.3)'),/remove_all)
         WIDGET_CONTROL,state.wtx_coltune_r,set_value=txt
         txt = strcompress(string(state.cmx[1],format='(E10.3)'),/remove_all)
         WIDGET_CONTROL,state.wtx_coltune_g,set_value=txt
         txt = strcompress(string(state.cmx[2],format='(E10.3)'),/remove_all)
         WIDGET_CONTROL,state.wtx_coltune_b,set_value=txt
      endif else begin
         txt = '1.0'
         WIDGET_CONTROL,state.wtx_coltune_r,set_value=txt
         WIDGET_CONTROL,state.wtx_coltune_g,set_value=txt
         WIDGET_CONTROL,state.wtx_coltune_b,set_value=txt
      endelse
      view_the_image,state
   end
   'wbg_con': begin
      view_the_image,state
   end
   'wtx_nrcont': begin
      view_the_image,state
   end
   'wsl_ct': begin
      WIDGET_CONTROL,state.wsl_ct,get_value=txt
      txt=strcompress(string(txt(0)),/remove_all)
      ct = txt + 0
      loadct,ct,/silent
      view_the_image,state
   end
   'wtx_obsx': begin
      view_the_image,state,/makeimage
   end
   'wtx_obsy': begin
      view_the_image,state,/makeimage
   end
   'wtx_obsz': begin
      view_the_image,state,/makeimage
   end
   'wbg_local': begin
      view_the_image,state,/makeimage
   end
   'wbg_relobs': begin
      view_the_image,state,/makeimage
   end
   'wtx_pntx': begin
      view_the_image,state,/makeimage
   end
   'wtx_pnty': begin
      view_the_image,state,/makeimage
   end
   'wtx_pntz': begin
      view_the_image,state,/makeimage
   end
   'wtx_scale': begin
      WIDGET_CONTROL,state.wtx_scale,get_value=txt
      scale=txt[0]+0.d0
      if scale gt 0.d0 then begin         
         WIDGET_CONTROL,state.wtx_obsx,get_value=txt
         txt=txt(0)
         obsx = txt + 0.d0
         WIDGET_CONTROL,state.wtx_obsy,get_value=txt
         txt=txt(0)
         obsy = txt + 0.d0
         WIDGET_CONTROL,state.wtx_obsz,get_value=txt
         txt=txt(0)
         obsz = txt + 0.d0
         obsx = scale*obsx
         obsy = scale*obsy
         obsz = scale*obsz
         WIDGET_CONTROL,state.wtx_obsx,set_value=  $
                strcompress(string(obsx,format='(E10.3)'),/remove_all)
         WIDGET_CONTROL,state.wtx_obsy,set_value=  $
                strcompress(string(obsy,format='(E10.3)'),/remove_all)
         WIDGET_CONTROL,state.wtx_obsz,set_value=  $
                strcompress(string(obsz,format='(E10.3)'),/remove_all)
         view_the_image,state,/makeimage
      endif
   end
   'wtx_view': begin
      view_the_image,state,/makeimage
   end
   'wdraw': begin
      case state.mousefunction of
         0: begin
            ;;
            ;; Buttom "mouse rotate" is not pressed, so use the mouse motions
            ;; for rotation
            ;;
            if sevent.type eq 0 then begin
               state.rotpts[0] = (1.d0*sevent.x-state.margins[0])/(1.d0*state.draw_xsize-state.margins[0]-state.margins[2])
               state.rotpts[1] = (1.d0*sevent.y-state.margins[1])/(1.d0*state.draw_ysize-state.margins[1]-state.margins[3])
            endif
            if sevent.type eq 1 then begin
               state.rotpts[2] = (1.d0*sevent.x-state.margins[0])/(1.d0*state.draw_xsize-state.margins[0]-state.margins[2])
               state.rotpts[3] = (1.d0*sevent.y-state.margins[1])/(1.d0*state.draw_ysize-state.margins[1]-state.margins[3])
               WIDGET_CONTROL,state.wtx_incl,get_value=txt
               txt=strcompress(string(txt(0)),/remove_all)
               incl = ( txt + 0.e0 )
               WIDGET_CONTROL,state.wtx_phi,get_value=txt
               txt=strcompress(string(txt(0)),/remove_all)
               phi = ( txt + 0.e0 )
               incl = incl + ( state.rotpts[3] - state.rotpts[1] ) * 120.d0
               ;if incl lt 0.d0 then incl = 0.d0
               ;if incl gt 180.d0 then incl = 180.d0
               if sin(!dpi*incl/180.) ge 0 then begin
                  phi  = phi  + ( state.rotpts[2] - state.rotpts[0] ) * 120.d0
               endif else begin
                  phi  = phi  - ( state.rotpts[2] - state.rotpts[0] ) * 120.d0
               endelse
               WIDGET_CONTROL,state.wtx_incl,set_value=string(format='(F9.2)',incl)
               WIDGET_CONTROL,state.wtx_phi,set_value=string(format='(F9.2)',phi)
            endif
            ;;
            ;; Replot
            ;; 
            view_the_image,state,/makeimage
         end
         1: begin
            ;;
            ;; Use the mouse motions as zoom-in tool, center dragged
            ;;
            makeimage = 0
            if sevent.type eq 0 then begin
               ;;
               ;; Start of the box
               ;;
               widget_control, state.wdr_main, draw_motion_events = 1
               dum = (1.d0*sevent.x-state.margins[0])/(1.d0*state.draw_xsize-state.margins[0]-state.margins[2])
               state.zoombox[0] = state.zoomcur[0] + (state.zoomcur[1]-state.zoomcur[0])*dum
               dum = (1.d0*sevent.y-state.margins[1])/(1.d0*state.draw_ysize-state.margins[1]-state.margins[3])
               state.zoombox[2] = state.zoomcur[2] + (state.zoomcur[3]-state.zoomcur[2])*dum
            endif
            if sevent.type eq 1 or sevent.type eq 2 then begin
               if sevent.type eq 1 then begin
                  ;;
                  ;; End of the box
                  ;;
                  widget_control, state.wdr_main, draw_motion_events = 0
               endif            
               ;;
               ;; Plot the box
               ;;
               dumx = (1.d0*sevent.x-state.margins[0])/(1.d0*state.draw_xsize-state.margins[0]-state.margins[2])
               state.zoombox[1] = state.zoomcur[0] + (state.zoomcur[1]-state.zoomcur[0])*dumx
               dumy = (1.d0*sevent.y-state.margins[1])/(1.d0*state.draw_ysize-state.margins[1]-state.margins[3])
               state.zoombox[3] = state.zoomcur[2] + (state.zoomcur[3]-state.zoomcur[2])*dumy
               if sevent.type eq 1 and ( state.zoombox[1] eq state.zoombox[0] or state.zoombox[3] eq state.zoombox[2] ) then begin
                  ;;
                  ;; If endpoint but a box of size zero, then no zoom box
                  ;;
                  state.zoombox[*] = -1d99
                  ;;
                  ;; If this zero-size box is outside of the frame/image, then
                  ;; the user wants to move the panel to left/right/up/down.
                  ;;
                  ddd = 0.1
                  if dumx lt 0 then begin
                     dx               = state.zoomcur[1]-state.zoomcur[0]
                     state.zoomcur[0] = state.zoomcur[0] - ddd*dx
                     state.zoomcur[1] = state.zoomcur[1] - ddd*dx
                     makeimage = 1
                  endif
                  if dumx gt 1 then begin
                     dx               = state.zoomcur[1]-state.zoomcur[0]
                     state.zoomcur[0] = state.zoomcur[0] + ddd*dx
                     state.zoomcur[1] = state.zoomcur[1] + ddd*dx
                     makeimage = 1
                  endif
                  if dumy lt 0 then begin
                     dy               = state.zoomcur[3]-state.zoomcur[2]
                     state.zoomcur[2] = state.zoomcur[2] - ddd*dy
                     state.zoomcur[3] = state.zoomcur[3] - ddd*dy
                     makeimage = 1
                  endif
                  if dumy gt 1 then begin
                     dy               = state.zoomcur[3]-state.zoomcur[2]
                     state.zoomcur[2] = state.zoomcur[2] + ddd*dy
                     state.zoomcur[3] = state.zoomcur[3] + ddd*dy
                     makeimage = 1
                  endif
               endif else begin
                  ;;
                  ;; Make sure the box is square
                  ;;
                  dx = abs( state.zoombox[1]-state.zoombox[0] )
                  dy = abs( state.zoombox[3]-state.zoombox[2] )
                  if dx lt dy then begin
                     if dx eq 0.d0 then begin
                        state.zoombox[1] = state.zoombox[0] + dy
                     endif else begin
                        dum = (state.zoombox[1]-state.zoombox[0])*dy/dx
                        state.zoombox[1] = state.zoombox[0] + dum
                     endelse
                  endif else begin
                     if dy eq 0.d0 then begin
                        state.zoombox[3] = state.zoombox[2] + dx
                     endif else begin
                        dum = (state.zoombox[3]-state.zoombox[2])*dx/dy
                        state.zoombox[3] = state.zoombox[2] + dum
                     endelse
                  endelse
                  ;;
                  ;; Make sure the box is such that x1>x0 and y1>y0
                  ;; But only when the box is finished.
                  ;;
                  if sevent.type eq 1 then begin
                     if state.zoombox[1] lt state.zoombox[0] then begin
                        dum = state.zoombox[1] 
                        state.zoombox[1] = state.zoombox[0]
                        state.zoombox[0] = dum
                     endif
                     if state.zoombox[3] lt state.zoombox[2] then begin
                        dum = state.zoombox[3] 
                        state.zoombox[3] = state.zoombox[2]
                        state.zoombox[2] = dum
                     endif
                  endif
               endelse
            endif 
            ;;
            ;; Replot
            ;; 
            view_the_image,state,makeimage=makeimage
         end
         2: begin
            ;;
            ;; Use the mouse motions as zoom-in tool, center dragged
            ;;
            makeimage = 0
            if sevent.type eq 0 then begin
               ;;
               ;; Start of the box
               ;;
               widget_control, state.wdr_main, draw_motion_events = 1
               dum = (1.d0*sevent.x-state.margins[0])/(1.d0*state.draw_xsize-state.margins[0]-state.margins[2])
               state.zoomcen[0] = state.zoomcur[0] + (state.zoomcur[1]-state.zoomcur[0])*dum
               dum = (1.d0*sevent.y-state.margins[1])/(1.d0*state.draw_ysize-state.margins[1]-state.margins[3])
               state.zoomcen[1] = state.zoomcur[2] + (state.zoomcur[3]-state.zoomcur[2])*dum
            endif
            if sevent.type eq 1 or sevent.type eq 2 then begin
               if sevent.type eq 1 then begin
                  ;;
                  ;; End of the box
                  ;;
                  widget_control, state.wdr_main, draw_motion_events = 0
               endif            
               ;;
               ;; Plot the box
               ;;
               dumx = (1.d0*sevent.x-state.margins[0])/(1.d0*state.draw_xsize-state.margins[0]-state.margins[2])
               state.zoombox[1] = state.zoomcur[0] + (state.zoomcur[1]-state.zoomcur[0])*dumx
               dumy = (1.d0*sevent.y-state.margins[1])/(1.d0*state.draw_ysize-state.margins[1]-state.margins[3])
               state.zoombox[3] = state.zoomcur[2] + (state.zoomcur[3]-state.zoomcur[2])*dumy
               state.zoombox[0] = state.zoomcen[0] - (state.zoombox[1]-state.zoomcen[0])
               state.zoombox[2] = state.zoomcen[1] - (state.zoombox[3]-state.zoomcen[1])
               if sevent.type eq 1 and ( state.zoombox[1] eq state.zoombox[0] or state.zoombox[3] eq state.zoombox[2] ) then begin
                  ;;
                  ;; If endpoint but a box of size zero, then no zoom box
                  ;;
                  state.zoombox[*] = -1d99
                  ;;
                  ;; If this zero-size box is outside of the frame/image, then
                  ;; the user wants to move the panel to left/right/up/down.
                  ;;
                  ddd = 0.1
                  if dumx lt 0 then begin
                     dx               = state.zoomcur[1]-state.zoomcur[0]
                     state.zoomcur[0] = state.zoomcur[0] - ddd*dx
                     state.zoomcur[1] = state.zoomcur[1] - ddd*dx
                     makeimage = 1
                  endif
                  if dumx gt 1 then begin
                     dx               = state.zoomcur[1]-state.zoomcur[0]
                     state.zoomcur[0] = state.zoomcur[0] + ddd*dx
                     state.zoomcur[1] = state.zoomcur[1] + ddd*dx
                     makeimage = 1
                  endif
                  if dumy lt 0 then begin
                     dy               = state.zoomcur[3]-state.zoomcur[2]
                     state.zoomcur[2] = state.zoomcur[2] - ddd*dy
                     state.zoomcur[3] = state.zoomcur[3] - ddd*dy
                     makeimage = 1
                  endif
                  if dumy gt 1 then begin
                     dy               = state.zoomcur[3]-state.zoomcur[2]
                     state.zoomcur[2] = state.zoomcur[2] + ddd*dy
                     state.zoomcur[3] = state.zoomcur[3] + ddd*dy
                     makeimage = 1
                  endif
               endif else begin
                  ;;
                  ;; Make sure the box is square
                  ;;
                  dx = abs( state.zoombox[1]-state.zoombox[0] )
                  dy = abs( state.zoombox[3]-state.zoombox[2] )
                  dd = max([dx,dy])
                  state.zoombox[0] = state.zoomcen[0]-dd*0.5d0
                  state.zoombox[1] = state.zoomcen[0]+dd*0.5d0
                  state.zoombox[2] = state.zoomcen[1]-dd*0.5d0
                  state.zoombox[3] = state.zoomcen[1]+dd*0.5d0
               endelse
            endif 
            ;;
            ;; Replot
            ;; 
            view_the_image,state,makeimage=makeimage
         end
         3: begin
            ;;
            ;; Use the mouse motions as zoom-in tool, funny dragged
            ;; (this strange way was found by accident, but could be
            ;; useful!)
            ;;
            makeimage = 0
            if sevent.type eq 0 then begin
               ;;
               ;; Start of the box
               ;;
               widget_control, state.wdr_main, draw_motion_events = 1
               dum = (1.d0*sevent.x-state.margins[0])/(1.d0*state.draw_xsize-state.margins[0]-state.margins[2])
               state.zoomcen[0] = state.zoomcur[0] + (state.zoomcur[1]-state.zoomcur[0])*dum
               dum = (1.d0*sevent.y-state.margins[1])/(1.d0*state.draw_ysize-state.margins[1]-state.margins[3])
               state.zoomcen[1] = state.zoomcur[2] + (state.zoomcur[3]-state.zoomcur[2])*dum
            endif
            if sevent.type eq 1 or sevent.type eq 2 then begin
               if sevent.type eq 1 then begin
                  ;;
                  ;; End of the box
                  ;;
                  widget_control, state.wdr_main, draw_motion_events = 0
               endif            
               ;;
               ;; Plot the box
               ;;
               dumx = (1.d0*sevent.x-state.margins[0])/(1.d0*state.draw_xsize-state.margins[0]-state.margins[2])
               state.zoombox[1] = state.zoomcur[0] + (state.zoomcur[1]-state.zoomcur[0])*dumx
               dumy = (1.d0*sevent.y-state.margins[1])/(1.d0*state.draw_ysize-state.margins[1]-state.margins[3])
               state.zoombox[3] = state.zoomcur[2] + (state.zoomcur[3]-state.zoomcur[2])*dumy
               state.zoombox[0] = state.zoomcen[0] - (state.zoombox[1]-state.zoomcen[0])
               state.zoombox[2] = state.zoomcen[1] - (state.zoombox[3]-state.zoomcen[1])
               if sevent.type eq 1 and ( state.zoombox[1] eq state.zoombox[0] or state.zoombox[3] eq state.zoombox[2] ) then begin
                  ;;
                  ;; If endpoint but a box of size zero, then no zoom box
                  ;;
                  state.zoombox[*] = -1d99
                  ;;
                  ;; If this zero-size box is outside of the frame/image, then
                  ;; the user wants to move the panel to left/right/up/down.
                  ;;
                  ddd = 0.1
                  if dumx lt 0 then begin
                     dx               = state.zoomcur[1]-state.zoomcur[0]
                     state.zoomcur[0] = state.zoomcur[0] - ddd*dx
                     state.zoomcur[1] = state.zoomcur[1] - ddd*dx
                     makeimage = 1
                  endif
                  if dumx gt 1 then begin
                     dx               = state.zoomcur[1]-state.zoomcur[0]
                     state.zoomcur[0] = state.zoomcur[0] + ddd*dx
                     state.zoomcur[1] = state.zoomcur[1] + ddd*dx
                     makeimage = 1
                  endif
                  if dumy lt 0 then begin
                     dy               = state.zoomcur[3]-state.zoomcur[2]
                     state.zoomcur[2] = state.zoomcur[2] - ddd*dy
                     state.zoomcur[3] = state.zoomcur[3] - ddd*dy
                     makeimage = 1
                  endif
                  if dumy gt 1 then begin
                     dy               = state.zoomcur[3]-state.zoomcur[2]
                     state.zoomcur[2] = state.zoomcur[2] + ddd*dy
                     state.zoomcur[3] = state.zoomcur[3] + ddd*dy
                     makeimage = 1
                  endif
               endif else begin
                  ;;
                  ;; Make sure the box is square
                  ;;
                  dx = abs( state.zoombox[1]-state.zoombox[0] )
                  dy = abs( state.zoombox[3]-state.zoombox[2] )
                  if dx lt dy then begin
                     if dx eq 0.d0 then begin
                        state.zoombox[1] = state.zoombox[0] + dy
                     endif else begin
                        dum = (state.zoombox[1]-state.zoombox[0])*dy/dx
                        state.zoombox[1] = state.zoombox[0] + dum
                     endelse
                  endif else begin
                     if dy eq 0.d0 then begin
                        state.zoombox[3] = state.zoombox[2] + dx
                     endif else begin
                        dum = (state.zoombox[3]-state.zoombox[2])*dx/dy
                        state.zoombox[3] = state.zoombox[2] + dum
                     endelse
                  endelse
                  ;;
                  ;; Make sure the box is such that x1>x0 and y1>y0
                  ;; But only when the box is finished.
                  ;;
                  if sevent.type eq 1 then begin
                     if state.zoombox[1] lt state.zoombox[0] then begin
                        dum = state.zoombox[1] 
                        state.zoombox[1] = state.zoombox[0]
                        state.zoombox[0] = dum
                     endif
                     if state.zoombox[3] lt state.zoombox[2] then begin
                        dum = state.zoombox[3] 
                        state.zoombox[3] = state.zoombox[2]
                        state.zoombox[2] = dum
                     endif
                  endif
               endelse
            endif 
            ;;
            ;; Replot
            ;; 
            view_the_image,state,makeimage=makeimage
         end
      endcase
   end
   else: notrecog=1
endcase
;
; Check if it could be one of the usertrans variables
;
if notrecog eq 1 then begin
   if strmid(uvalue,0,8) eq 'wtx_uvar' and state.nvars gt 0 then begin
      ivar = strcompress(strmid(uvalue,9,8),/remove_all)+0
      if ivar lt 0 or ivar gt state.nvars-1 then begin
         print,'ERROR: ivar out of range'
         stop
      endif
      WIDGET_CONTROL,state.wtx_utrans[ivar],get_value=txt
      value=txt[0]+0.d0
      state.vals[ivar]=value
      view_the_image,state,/makeimage
   endif
endif
WIDGET_CONTROL,sevent.top,set_uvalue=state,/no_copy
end


;----------------------------------------------------------
; RENDER THE IMAGE AND VIEW IT ON SCREEN OR IN POSTSCRIPT
;----------------------------------------------------------
pro view_the_image,state,makeimage=makeimage,ps=ps,unzoom=unzoom,linelist=linelist
;
; Get the current values of the widgets
;
WIDGET_CONTROL,state.wtx_maxlog,get_value=txt
txt=txt(0)
maxlog = txt + 0.d0
WIDGET_CONTROL,state.wtx_nx,get_value=txt
txt=strcompress(string(txt(0)),/remove_all)
npix = txt + 0
if npix gt state.npmax then npix=state.npmax
WIDGET_CONTROL,state.wtx_sx,get_value=txt
txt=strcompress(string(txt(0)),/remove_all)
size = txt + 0.e0
WIDGET_CONTROL,state.wtx_incl,get_value=txt
txt=strcompress(string(txt(0)),/remove_all)
incl = ( txt + 0.e0 )
WIDGET_CONTROL,state.wtx_phi,get_value=txt
txt=strcompress(string(txt(0)),/remove_all)
phi = ( txt + 0.e0 )
;WIDGET_CONTROL,state.wtx_posang,get_value=txt
;txt=strcompress(string(txt(0)),/remove_all)
;posang = ( txt + 0.e0 ) 
posang = 0.e0
WIDGET_CONTROL,state.wtx_nphot,get_value=txt
txt=strcompress(string(txt(0)),/remove_all)
nphot = ( txt + 0LL ) 
;;WIDGET_CONTROL,state.wsl_lam_r,get_value=txt
;;txt=strcompress(string(txt(0)),/remove_all)
;;ilam_r = txt + 0
WIDGET_CONTROL,state.wtx_lam_r,get_value=txt
txt=strcompress(string(txt(0)),/remove_all)
lam_r = txt + 0.d0
if keyword_set(state.color) then begin
   ;;WIDGET_CONTROL,state.wsl_lam_g,get_value=txt
   ;;txt=strcompress(string(txt(0)),/remove_all)
   ;;ilam_g = txt + 0
   ;;WIDGET_CONTROL,state.wsl_lam_b,get_value=txt
   ;;txt=strcompress(string(txt(0)),/remove_all)
   ;;ilam_b = txt + 0
   WIDGET_CONTROL,state.wtx_lam_g,get_value=txt
   txt=strcompress(string(txt(0)),/remove_all)
   lam_g = txt + 0
   WIDGET_CONTROL,state.wtx_lam_b,get_value=txt
   txt=strcompress(string(txt(0)),/remove_all)
   lam_b = txt + 0
   ;;ilambda = [ilam_r,ilam_g,ilam_b]
   lambda = [lam_r,lam_g,lam_b]
   WIDGET_CONTROL,state.wbg_scl,get_value=txt
   colorscale = txt + 0
   if colorscale eq 0 then begin
      coltune = dblarr(3)
      WIDGET_CONTROL,state.wtx_coltune_r,get_value=txt
      txt=strcompress(string(txt(0)),/remove_all)
      coltune[0] = txt + 0.d0
      WIDGET_CONTROL,state.wtx_coltune_g,get_value=txt
      txt=strcompress(string(txt(0)),/remove_all)
      coltune[1] = txt + 0.d0
      WIDGET_CONTROL,state.wtx_coltune_b,get_value=txt
      txt=strcompress(string(txt(0)),/remove_all)
      coltune[2] = txt + 0.d0
   endif else begin
      lgrange = dblarr(2,3)
      WIDGET_CONTROL,state.wtx_coltune_r,get_value=txt
      txt=strcompress(string(txt(0)),/remove_all)
      lgrange[1,0] = alog10(txt+0.d0)
      WIDGET_CONTROL,state.wtx_coltune_g,get_value=txt
      txt=strcompress(string(txt(0)),/remove_all)
      lgrange[1,1] = alog10(txt+0.d0)
      WIDGET_CONTROL,state.wtx_coltune_b,get_value=txt
      txt=strcompress(string(txt(0)),/remove_all)
      lgrange[1,2] = alog10(txt+0.d0)
      lgrange[0,*] = lgrange[1,*] - maxlog
   endelse
endif else begin
   if state.lines ne 0 and state.doline eq 1 then begin
      WIDGET_CONTROL,state.wtx_imol,get_value=txt
      txt=strcompress(string(txt(0)),/remove_all)
      imol = txt + 0
      WIDGET_CONTROL,state.wtx_iline,get_value=txt
      txt=strcompress(string(txt(0)),/remove_all)
      iline = txt + 0
      WIDGET_CONTROL,state.wtx_vkms,get_value=txt
      txt=strcompress(string(txt(0)),/remove_all)
      vkms = txt + 0.d0
   endif else begin
      lambda = lam_r
   endelse
endelse
WIDGET_CONTROL,state.wbg_lilg,get_value=txt
txt=txt(0)
plot_log = 1 - txt 
WIDGET_CONTROL,state.wbg_second,get_value=txt
txt=txt(0)
secondorder = txt + 0
if state.incldcatchbutton eq 1 then begin
   WIDGET_CONTROL,state.wbg_dcatch,get_value=txt
   txt=txt(0)
   dcatch = txt + 0
endif else dcatch=0
WIDGET_CONTROL,state.wtx_satur,get_value=txt
txt=txt(0)
satur = txt + 0.d0
WIDGET_CONTROL,state.wbg_star,get_value=txt
txt=strcompress(string(txt(0)),/remove_all)
star = txt + 0
nostar = 1-star
if state.local eq 0 then begin
   WIDGET_CONTROL,state.wbg_prev,get_value=txt
   txt=strcompress(string(txt(0)),/remove_all)
   nofluxcons = txt + 0
endif else begin
   nofluxcons = 1
endelse
WIDGET_CONTROL,state.wbg_fixseed,get_value=txt
txt=strcompress(string(txt(0)),/remove_all)
fixseed = txt + 0
WIDGET_CONTROL,state.wbg_tau,get_value=txt
txt=strcompress(string(txt(0)),/remove_all)
plottau = txt + 0
if state.color eq 0 then begin
   WIDGET_CONTROL,state.wbg_con,get_value=txt
   txt=txt(0)
   plotcont = txt 
   WIDGET_CONTROL,state.wtx_nrcont,get_value=txt
   txt=txt(0)
   nlevels = txt
endif else begin
   plotcont = 0
   nlevels  = 0
endelse 
if state.local eq 1 then begin
   WIDGET_CONTROL,state.wbg_local,get_value=txt
   txt=txt(0)
   localobs = txt + 0
   if localobs eq 1 then begin
      state.zoom = 0
      state.zoombox = [-1d99,-1d99,-1d99,-1d99]
   endif
   WIDGET_CONTROL,state.wbg_relobs,get_value=txt
   txt=txt(0)
   relobs = txt + 0
   WIDGET_CONTROL,state.wtx_pntx,get_value=txt
   txt=txt(0)
   pntx = txt + 0.d0
   WIDGET_CONTROL,state.wtx_pnty,get_value=txt
   txt=txt(0)
   pnty = txt + 0.d0
   WIDGET_CONTROL,state.wtx_pntz,get_value=txt
   txt=txt(0)
   pntz = txt + 0.d0
   WIDGET_CONTROL,state.wtx_obsx,get_value=txt
   txt=txt(0)
   obsx = txt + 0.d0
   WIDGET_CONTROL,state.wtx_obsy,get_value=txt
   txt=txt(0)
   obsy = txt + 0.d0
   WIDGET_CONTROL,state.wtx_obsz,get_value=txt
   txt=txt(0)
   obsz = txt + 0.d0
   WIDGET_CONTROL,state.wtx_view,get_value=txt
   txt=txt(0)
   viewang = txt + 0.d0
   ;;WIDGET_CONTROL,state.wtx_posang,get_value=txt
   ;;txt=txt(0)
   ;;posang = txt + 0.d0
   posang = 0.d0
   WIDGET_CONTROL,state.wtx_nphot,get_value=txt
   txt=txt(0)
   nphot = txt + 0LL
   point  = [pntx,pnty,pntz]
   if relobs eq 0 then begin
      locobs = [obsx,obsy,obsz]
   endif else begin
      locobs = [pntx+obsx,pnty+obsy,pntz+obsz]
   endelse
endif else begin
   localobs = 0
   relobs = 0
   point  = [0.,0.,0.]
   locobs = [0.,0.,0.]
endelse
if state.nvars gt 0 then begin
   ;;
   ;; Read all the usertrans variable values from the widget
   ;;
   for ivar=0,state.nvars-1 do begin
      WIDGET_CONTROL,state.wtx_utrans[ivar],get_value=txt
      value=txt[0]+0.d0
      state.vals[ivar]=value
   endfor
   ;;
   ;; Write them to transfer.inp
   ;;
   openw,1,'transfer.inp'
   for ivar=0,state.nvars-1 do begin
      printf,1,state.vars[ivar],'=',state.vals[ivar]
   endfor
   close,1
   ;;
   ;; Ask RADMC-3D to read this file
   ;; 
   if state.iounit ne 0 then begin
      printf,state.iounit,'myaction'
      printf,state.iounit,'enter'
      flush,state.iounit
   endif
endif
;
; Convert size
;
unit = 1.d0
if state.pc ne 0 then begin
   unit = 3.08572d18
   if state.au ne 0 then stop    ; Not both au and pc possible
endif
if state.au ne 0 then begin
   unit = 1.496d13
endif
;
; If not in zoom-mode, then use size to make the box, else interpret
; the zoom coordinates
;
if state.zoom eq 0 then begin
   ;;
   ;; Normal unzoomed mode, so create the zoomcur from size
   ;;
   state.zoomcur=[-size,size,-size,size]/2.d0
   ;;
endif
;
; Line list?
;
if not keyword_set(linelist) then linelist=0
;
; Make sure that the wavelength stays within the bounds of
; the wavelength_micron.inp
;
if state.lines eq 0 or state.doline ne 1 then begin
   if lambda gt max(state.lambda) then begin
      if abs(lambda-max(state.lambda)) lt 1d-3*lambda then begin
         lambda = max(state.lambda)
      endif else begin
         print,'ERROR: Wavelength out of range: lambda too large'
         stop
      endelse
   endif
   if lambda lt min(state.lambda) then begin
      if abs(lambda-min(state.lambda)) lt 1d-3*lambda then begin
         lambda = min(state.lambda)
      endif else begin
         print,'ERROR: Wavelength out of range: lambda too small'
         stop
      endelse
   endif
endif
;
; Now render the image
;
if keyword_set(makeimage) then begin
   ;;
   ;; If there is a zoombox, then use this and put zoom to 1
   ;; Else unzoom
   ;;
   if min(state.zoombox) gt -1d90 then begin
      ;;
      ;; Yes, a box is defined, so zoom!
      ;;
      state.zoomcur = state.zoombox
      state.zoom = 1
      state.zoombox[*] = -1d99
   endif else begin
      ;;
      ;; No, a box is not defined
      ;;
      if keyword_set(unzoom) then begin
         ;;
         ;; Unzoom, i.e. go back to original geometry
         ;;
         state.zoom = 0
         state.zoomcur=[-size,size,-size,size]/2.d0
      endif
   endelse
   ;;
   ;; If plotting to the widget, then print a statement on the screen
   ;;
   if not keyword_set(ps) then begin
      oldwin=1
      wincur=!d.window 
      wset,state.bitmap
      xyouts,10,10,/device,'Please wait... RADMC-3D is working...'
      wset,state.window
      device,copy=[0,0,state.draw_xsize,state.draw_ysize,0,0,state.bitmap]
      if oldwin ne 0 then wset,wincur
   endif
   ;;
   ;; Call Radmc3d to compute the image
   ;;
   if localobs eq 0 then begin
      makeimage,iounit=state.iounit,incl=incl,phi=phi,npix=npix,posang=posang,$
             iline=iline,ispec=imol,vkms=vkms,secondorder=secondorder,$
             lambda=lambda,lines=state.lines,linelist=linelist,dcatch=dcatch,$
             nostar=nostar,zoomau=state.zoomcur*unit/1.496d13,nphot=nphot,$
             nofluxcons=nofluxcons,plottau=plottau,fixseed=fixseed,radmc3d=state.radmc3d
   endif else begin
      makeimagelocalobs,iounit=state.iounit,npix=npix,posang=posang,nofluxcons=nofluxcons,$
              pointau=point*unit/1.496d13,lines=state.lines,dcatch=dcatch,$
              locobsau=locobs*unit/1.496d13,sizeradian=viewang,nphot=nphot,$
              ifreq=ifreq,lambda=lambda,nostar=nostar,secondorder=secondorder,$
              plottau=plottau,linelist=linelist,iline=iline,ispec=ispec,vkms=vkms,$
              fixseed=fixseed,radmc3d=state.radmc3d
   endelse
   ;;
   ;; Read the image
   ;;
   image=readimage(iounit=state.iounit)
   if image.nx gt state.npmax or image.ny gt state.npmax then begin
      print,'IMAGE TOO LARGE...'
      if state.iounit gt 1 then printf,iounit,'quit'
      stop
   endif
   state.image.nx                         = image.nx
   state.image.ny                         = image.ny
   state.image.nrfr                       = image.nrfr
   state.image.sizepix_x                  = image.sizepix_x
   state.image.sizepix_y                  = image.sizepix_y
   state.image.flux[*]                    = image.flux[*]
   state.image.x[0:image.nx-1]            = image.x[*]
   state.image.y[0:image.ny-1]            = image.y[*]
   state.image.lambda[*]                  = image.lambda[*]
   state.image.image[0:image.nx-1,0:image.ny-1,*] = image.image[0:image.nx-1,0:image.ny-1,*]
   state.image.radian                     = image.radian
   state.image.stokes                     = image.stokes
   ;;
   ;; Refresh the wavelength indicator
   ;;
   if state.doline eq 1 then begin
      WIDGET_CONTROL,state.wtx_lam_r,set_value=string(format='(F10.4)',image.lambda[0])
   endif
   ;;
   ;; If RADMC-3D is linked via a tunnel, then write out the image to file,
   ;; so that the user can analyze it separately.
   ;;
   if keyword_set(state.iounit) then begin
      openw,1,'image.out'
      printf,1,1
      printf,1,image.nx,image.ny
      printf,1,image.nrfr
      printf,1,image.sizepix_x,image.sizepix_y
      printf,1,image.lambda,format='(E21.14)'
      printf,1,image.image,format='(E21.14,1X)'
      close,1
   endif
endif
;
; Make local version of image
;
image = {nx:state.image.nx,ny:state.image.ny,nrfr:state.image.nrfr,$
         sizepix_x:state.image.sizepix_x,sizepix_y:state.image.sizepix_y,$
         image:state.image.image[0:state.image.nx-1,0:state.image.ny-1,*],$
         flux:state.image.flux,x:state.image.x[0:state.image.nx-1],$
         y:state.image.y[0:state.image.ny-1],lambda:state.image.lambda,$
         radian:state.image.radian,stokes:state.image.stokes}
;
; Now check where to plot to
;
if localobs eq 0 then begin
   zoomcur = state.zoomcur
endif
if keyword_set(ps) then begin
   ;;
   ;; Plot to postscript
   ;;
   print,'Writing idl.ps image...'
   set_plot,'ps',/interpolate
   device,file='idl.ps',/portrait,/color,bits_per_pixel=24,font_size=fontsize,$
          xsize=16.0,ysize=16.0,yoffset=1
   plotimage,image,pc=state.pc,au=state.au,log=plot_log,maxlog=maxlog,$
           zoom=zoomcur,coltune=coltune,saturate=satur,contour=plotcont,$
           nlevels=nlevels,cmn=cmn,cmx=cmx,lgrange=lgrange
   device,/close
   set_plot,'x'
   print,'..done ploting to postscript'
   ;;
endif else begin
   ;;
   ;; Plot to window 
   ;;
   oldwin=1
   wincur=!d.window 
   wset,state.bitmap
   plotimage,image,pc=state.pc,au=state.au,log=plot_log,maxlog=maxlog,$
         position=[state.margins[0],state.margins[1],$
                   state.draw_xsize-state.margins[0]-state.margins[2],$
                   state.draw_ysize-state.margins[1]-state.margins[3]],$
         zoom=zoomcur,coltune=coltune,saturate=satur,contour=plotcont,$
           nlevels=nlevels,cmn=cmn,cmx=cmx,lgrange=lgrange
   ;;
   ;; Write a little idl script that the user may want to use to reproduce
   ;; the current image
   ;;
   if localobs eq 0 then begin
      openw,20,'image_script.pro'
      printf,20,'@readradmc.pro'
      printf,20,'makeimage',+$
             ',incl='+strcompress(string(incl),/remove_all)+$
             ',phi='+strcompress(string(phi),/remove_all)+$
             ',npix='+strcompress(string(npix),/remove_all)+$
             ',fixseed='+strcompress(string(fixseed),/remove_all)+$
             ',plottau='+strcompress(string(plottau),/remove_all)+',$'
      if keyword_set(lambda) then begin
         printf,20,'      posang='+strcompress(string(posang),/remove_all)+$
             ',lambda='+strcompress(string(lambda),/remove_all)+$
             ',nostar='+strcompress(string(nostar),/remove_all)+$
             ',nofluxcons='+strcompress(string(nofluxcons),/remove_all)+',$'
      endif else begin
         printf,20,'      posang='+strcompress(string(posang),/remove_all)+$
             ',ispec='+strcompress(string(imol),/remove_all)+$
             ',iline='+strcompress(string(iline),/remove_all)+$
             ',vkms='+strcompress(string(vkms),/remove_all)+$
             ',nostar='+strcompress(string(nostar),/remove_all)+$
             ',nofluxcons='+strcompress(string(nofluxcons),/remove_all)+',$'
      endelse
      if secondorder eq 1 then begin
         printf,20,'      secondorder=1,$'
      endif
      if dcatch eq 1 then begin
         printf,20,'      dcatch=1,$'
      endif
      printf,20,'      zoomau=['+strcompress(string(state.zoomcur[0]*unit/1.496d13),/remove_all)+$
             ','+strcompress(string(state.zoomcur[1]*unit/1.496d13),/remove_all)+$
             ','+strcompress(string(state.zoomcur[2]*unit/1.496d13),/remove_all)+$
             ','+strcompress(string(state.zoomcur[3]*unit/1.496d13),/remove_all)+$
             ',"'+state.radmc3d+'"]'
      printf,20,'a=readimage()'
      printf,20,'plotimage,a,pc='+strcompress(string(state.pc),/remove_all)+$
             ',au='+strcompress(string(state.au),/remove_all)+$
             ',log='+strcompress(string(plot_log),/remove_all)+$
             ',maxlog='+strcompress(string(maxlog),/remove_all)+$
             ',saturate='+strcompress(string(satur),/remove_all)+',$'
      if keyword_set(state.color) then begin
         if n_elements(coltune) gt 0 then begin
            printf,20,'      coltune='+strcompress(string(coltune),/remove_all)+',$'
         endif
         if n_elements(lgrange) gt 0 then begin
            printf,20,'      lgrange='+strcompress(string(lgrange),/remove_all)+',$'
         endif
      endif
      printf,20,'      zoom=['+strcompress(string(state.zoomcur[0]),/remove_all)+$
          ','+strcompress(string(state.zoomcur[1]),/remove_all)+$
          ','+strcompress(string(state.zoomcur[2]),/remove_all)+$
          ','+strcompress(string(state.zoomcur[3]),/remove_all)+']'
      printf,20,'end'
      close,20
   endif
   ;;
   ;; Plot a zoom box if necessary
   ;;
   if min(state.zoombox) gt -1d90 then begin
      x = [state.zoombox[0],state.zoombox[1],state.zoombox[1],state.zoombox[0],state.zoombox[0]]
      y = [state.zoombox[2],state.zoombox[2],state.zoombox[3],state.zoombox[3],state.zoombox[2]]
      oplot,x,y
   endif
   ;;
   ;; Write the overal image-integrated flux or the max optical depth
   ;;
   if not keyword_set(plottau) and localobs eq 0 then begin
      case state.fluxunits of
         0: begin   
            xyouts,250,10,/device,'Flux = '+string(image.flux,format='(E10.4)')+$
                   ' erg/cm^2/s/Hz @ 1 pc'
         end
         1: begin
            xyouts,280,10,/device,'Flux = '+string(1d23*image.flux,format='(E10.4)')+$
                   ' Jy @ 1 pc'
         end
      endcase
   endif else begin
      xyouts,250,10,/device,'Max Tau = '+string(max(image.image),format='(E10.4)')
   endelse
   ;;
   ;; Set everything on the screen
   ;;
   wset,state.window
   device,copy=[0,0,state.draw_xsize,state.draw_ysize,0,0,state.bitmap]
   if oldwin ne 0 then wset,wincur
   ;;
endelse
;;
state.doline = 0
;;
state.cmn = cmn
state.cmx = cmx
;;
end


;----------------------------------------------------------
; RENDER THE SPECTRUM AND VIEW IT ON SCREEN OR IN POSTSCRIPT
;----------------------------------------------------------
pro view_the_spectrum,state,makespec=makespec
;
; Get the current values of the widgets
;
WIDGET_CONTROL,state.wtx_incl,get_value=txt
txt=strcompress(string(txt(0)),/remove_all)
incl = ( txt + 0.e0 )
WIDGET_CONTROL,state.wtx_phi,get_value=txt
txt=strcompress(string(txt(0)),/remove_all)
phi = ( txt + 0.e0 )
WIDGET_CONTROL,state.wbg_star,get_value=txt
txt=strcompress(string(txt(0)),/remove_all)
star = txt + 0
nostar = 1-star
;
; Convert size
;
unit = 1.d0
if state.pc ne 0 then begin
   unit = 3.08572d18
   if state.au ne 0 then stop    ; Not both au and pc possible
endif
if state.au ne 0 then begin
   unit = 1.496d13
endif
;
; Now render the spectrum
;
if keyword_set(makespec) then begin
   ;;
   ;; If plotting to the widget, then print a statement on the screen
   ;;
   if not keyword_set(ps) then begin
      oldwin=1
      wincur=!d.window 
      wset,state.bitmap
      xyouts,10,25,/device,'Please be pretty patient... RADMC-3D is working on the spectrum...'
      xyouts,10,10,/device,'This takes generally a lot of time.'
      wset,state.window
      device,copy=[0,0,state.draw_xsize,state.draw_ysize,0,0,state.bitmap]
      if oldwin ne 0 then wset,wincur
   endif
   ;;
   ;; Call Radmc3d to compute the spectrum
   ;;
   makespectrum,iounit=state.iounit,incl=incl,phi=phi,$
             zoomau=state.zoomcur*unit/1.496d13,nostar=nostar
   ;;
   ;; Message
   ;;
   if nostar ne 0 then begin
      print,'WARNING: Spectrum made without the stars!! Click the star button'
      print,'         and redo spectrum to get the stars included.'
   endif
   ;;
   ;; Read the spectrum
   ;;
   spectrum=readspectrum(iounit=state.iounit)
   ;;
   ;; If RADMC-3D is linked via a tunnel, then write out the spectrum to file,
   ;; so that the user can analyze it separately.
   ;;
   if keyword_set(state.iounit) then begin
      openw,1,'spectrum.out'
      printf,1,1
      printf,1,spectrum.nfr
      printf,1,' '
      for inu=0,spectrum.nfr-1 do begin
         printf,1,spectrum.lambda[inu],spectrum.spectrum[inu]
      endfor
      close,1
   endif
endif
;;
;; Plot to window 
;;
oldwin=1
wincur=!d.window 
wset,state.bitmap
plotspectrum,spectrum,dpc=1
;;
;; Write a little idl script that the user may want to use to reproduce
;; the current image
;;
openw,20,'spectrum_script.pro'
printf,20,'@readradmc.pro'
printf,20,'makespectrum',+$
       ',incl='+strcompress(string(incl),/remove_all)+$
       ',phi='+strcompress(string(phi),/remove_all)+',$'
printf,20,'      zoomau=['+strcompress(string(state.zoomcur[0]*unit/1.496d13),/remove_all)+$
       ','+strcompress(string(state.zoomcur[1]*unit/1.496d13),/remove_all)+$
       ','+strcompress(string(state.zoomcur[2]*unit/1.496d13),/remove_all)+$
       ','+strcompress(string(state.zoomcur[3]*unit/1.496d13),/remove_all)+']'
printf,20,'a=readspectrum()'
printf,20,'plotspectrum,a,dpc=1'
printf,20,'end'
close,20
;;
;; Set plot onto the screen
;;
wset,state.window
device,copy=[0,0,state.draw_xsize,state.draw_ysize,0,0,state.bitmap]
if oldwin ne 0 then wset,wincur
;;
end



