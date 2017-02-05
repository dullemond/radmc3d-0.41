@../../idl/readradmc.pro
@read_params_radmc3d.pro
@create_ext_tab.pro
;------------------------------------------------------------------------
;                     MAKE OPTICAL DEPTHS
; Copyright C.P. Dullemond
;------------------------------------------------------------------------
function maketau,a,kappa
sz=size(a.rho)
nr=sz(1)
nt=sz(2)
nspec=n_elements(kappa)
if sz[0] le 2 then nspecc=1 else nspecc=sz[3]
if nspec ne nspecc then begin
    print,'ERROR: Nr of opacities must equal nr of dust components'
    stop
endif
taur=fltarr(nr,nt)
taut=fltarr(nr,nt)
for it=0,nt-1 do begin
    taur(0,it) = 0.d0
    for ir=1,nr-1 do begin
        dummy = 0.5*(a.rho(ir,it,*)+a.rho(ir-1,it,*))*kappa[*] * $
                (a.r(ir)-a.r(ir-1))
        taur(ir,it) = taur(ir-1,it) + total(dummy)
    endfor
endfor
for ir=0,nr-1 do begin
    taut(ir,0) = 0.d0
    for it=1,nt-1 do begin
        dummy = 0.5*(a.rho(ir,it,*)+a.rho(ir,it-1,*))*kappa[*] * $
                (a.theta(it)-a.theta(it-1))*a.r(ir)
        taut(ir,it) = taut(ir,it-1) + total(dummy)
    endfor
endfor
return,{taur:taur,taut:taut}
end

;**********************************************************************************************
; Function to create the optical depth structure using the 'maketau' function            
;     Attila Juhasz -  Heidelberg - 04.01.2010
;**********************************************************************************************

function get_tau_1dust, dens=dens, wav=wav, nr=nr, $
                        dustopac_fname=dustopac_fname, imirt=imirt

  ;; Number of the dust species the density of which should be used
  ;; to calculate the optical depth
    if not keyword_set(nr) then nr = 1

  ;; Check if we have theta-mirroring and convert the structure
  ;;  into the RADMC2D convention
    ii = where(!pi/2.-dens.grid.theta lt 0.)
    if (ii(0) ge 0) then begin
       nt = n_elements(dens.grid.theta)/2
       dens2 = {r      : dens.grid.r,$
                theta  : dens.grid.theta(0:nt-1),$
                nr     : dens.grid.nr,$
                nt     : nt,$
                nspec  : n_elements(dens.rho(0,0,*)),$
                rho    : dens.rho(*,0:nt-1,*)}
    endif else begin
       nt = n_elements(dens.grid.theta)
       dens2 = {r      : dens.grid.r,$
                theta  : dens.grid.theta(0:nt-1),$
                nr     : dens.grid.nr,$
                nt     : dens.grid.ntheta,$
                nspec  : n_elements(dens.rho(0,0,*)),$
                rho    : dens.rho}
    endelse
            
    r = dens2.r
    theta = dens2.theta
    
    ndust = n_elements(dens2.rho(0,0,*))
    rho = dens2.rho(*,*,nr-1)
;
; Reading the dust opacities
;
    res = findfile(dustopac_fname(0), count=count)
    
    openr, 1, res
    readf, 1, iformat
    readf, 1, nlam
    if iformat eq 1 then begin
       dum = [0d0, 0d0]
       opac = {wave:dblarr(nlam), cabs:dblarr(nlam), csca:dblarr(nlam)}
       for ilam=0, nlam-1 do begin
          readf, 1, dum
          opac.wave(ilam) = dum(0)
          opac.cabs(ilam) = dum(1)
       endfor
    endif
    if iformat ge 2 then begin
       dum = [0d0, 0d0, 0d0]
       opac = {wave:dblarr(nlam), cabs:dblarr(nlam), csca:dblarr(nlam)}
       for ilam=0, nlam-1 do begin
          readf, 1, dum
          opac.wave(ilam) = dum(0)
          opac.cabs(ilam) = dum(1)
          opac.csca(ilam) = dum(2)
       endfor
    endif
    close, 1
    
    if opac.cabs(0) ne 0 then begin
       kappa = 10.^(interpol(alog10(opac.cabs), alog10(opac.wave), alog10(wav))) +$
               10.^(interpol(alog10(opac.csca), alog10(opac.wave), alog10(wav)))
    endif else begin
       kappa = 10.^(interpol(alog10(opac.cabs), alog10(opac.wave), alog10(wav))) 
    endelse

    kappa = 10.^(interpol(alog10(opac.cabs), alog10(opac.wave), alog10(wav))) 
    
    aa = {r:r, theta:theta, rho:rho}
    
    tau = maketau(aa, kappa)

    if not keyword_set(imirt) then begin
       nnr = n_elements(dens.grid.r)
       nnt = n_elements(dens.grid.theta)
       taur = dblarr(nnr, nnt)
       taut = dblarr(nnr, nnt)
       
       for ir=0, nnr-1 do begin
          taur(ir,0:nt-1) = tau.taur(ir,*)
          taut(ir,0:nt-1) = tau.taut(ir,*)
          taur(ir,nt:2*nt-1) = reverse(reform(tau.taur(ir,*)))
          taut(ir,nt:2*nt-1) = reverse(reform(tau.taut(ir,*)))
       endfor

       theta = dens.grid.theta

       t = {taur:taur, taut:taut}
    endif else begin
       t = tau
    endelse


    return, {r:r, theta:theta, rho:rho, kappa:kappa, wav:wav, tau:t}
 end
;**********************************************************************************************
; Procedure to read out the values of the input widgets 
;     Attila Juhasz -  Heidelberg - 04.01.2010
;**********************************************************************************************
pro read_all_widget_values, state, params

;
; Thermal Monte-Carlo
;

  ii = where(strcompress(params.par.name, /remove_all) eq 'nphot')
   if ii(0) ge 0 then WIDGET_CONTROL, state.run_nr_phot_inp, get_value=dum
   params.par(ii(0)).value = dum
  
;
; Central star
;
  ii = where(strcompress(params.par.name,/remove_all) eq 'star_enable')
  WIDGET_CONTROL, state.star_source_select, get_value=dum
  if dum eq 0 then params.par(ii(0)).value = '1'
  if dum eq 1 then params.par(ii(0)).value = '0'
  

  ii = where(strcompress(params.par.name, /remove_all) eq 'tstar')
   if ii(0) ge 0 then WIDGET_CONTROL, state.star_source_temp_inp, get_value=dum
   params.par(ii(0)).value = dum
  ii = where(strcompress(params.par.name, /remove_all) eq 'rstar')
   if ii(0) ge 0 then WIDGET_CONTROL, state.star_source_rad_inp, get_value=dum
   params.par(ii(0)).value = dum
  ii = where(strcompress(params.par.name, /remove_all) eq 'mstar')
   if ii(0) ge 0 then WIDGET_CONTROL, state.star_source_mass_inp, get_value=dum
   params.par(ii(0)).value = dum
  
;
; Continuous starlike radiation source
;
  ii = where(strcompress(params.par.name,/remove_all) eq 'cstar_enable')
  WIDGET_CONTROL, state.cstar_source_select, get_value=dum
  if dum eq 0 then params.par(ii(0)).value = '1'
  if dum eq 1 then params.par(ii(0)).value = '0'

  
  ii = where(strcompress(params.par.name, /remove_all) eq 'cstar_rin')
   if ii(0) ge 0 then WIDGET_CONTROL, state.cstar_source_rin_inp, get_value=dum
   params.par(ii(0)).value = dum
  ii = where(strcompress(params.par.name, /remove_all) eq 'cstar_rout')
   if ii(0) ge 0 then WIDGET_CONTROL, state.cstar_source_rout_inp, get_value=dum
   params.par(ii(0)).value = dum
  ii = where(strcompress(params.par.name, /remove_all) eq 'cstar_trin')
   if ii(0) ge 0 then WIDGET_CONTROL, state.cstar_source_tin_inp, get_value=dum
   params.par(ii(0)).value = dum
  ii = where(strcompress(params.par.name, /remove_all) eq 'cstar_powex')
   if ii(0) ge 0 then WIDGET_CONTROL, state.cstar_source_powex_inp, get_value=dum
   params.par(ii(0)).value = dum

;
; External radiation source
;

  ii = where(strcompress(params.par.name,/remove_all) eq 'ext_enable')
  WIDGET_CONTROL, state.ext_source_select, get_value=dum
  if dum eq 0 then params.par(ii(0)).value = '1'
  if dum eq 1 then params.par(ii(0)).value = '0'

 
  ii = where(strcompress(params.par.name, /remove_all) eq 'ext_temp')
   if ii(0) ge 0 then WIDGET_CONTROL, state.ext_source_temp_inp, get_value=dum
   params.par(ii(0)).value = dum

  
;
; Circumstellar disk
;

  ii = where(strcompress(params.par.name,/remove_all) eq 'ppdisk_enable')
  WIDGET_CONTROL, state.ppdisk_select, get_value=dum
  if dum eq 0 then params.par(ii(0)).value = '1'
  if dum eq 1 then params.par(ii(0)).value = '0'

 
  ii = where(strcompress(params.par.name, /remove_all) eq 'ppdisk_rin')
   if ii(0) ge 0 then WIDGET_CONTROL, state.ppdisk_rin_inp, get_value=dum
   params.par(ii(0)).value = dum
  ii = where(strcompress(params.par.name, /remove_all) eq 'ppdisk_rout')
   if ii(0) ge 0 then WIDGET_CONTROL, state.ppdisk_rout_inp, get_value=dum
   params.par(ii(0)).value = dum
  ii = where(strcompress(params.par.name, /remove_all) eq 'ppdisk_sig0')
   if ii(0) ge 0 then WIDGET_CONTROL, state.ppdisk_sig0_inp, get_value=dum
   params.par(ii(0)).value = dum
  ii = where(strcompress(params.par.name, /remove_all) eq 'ppdisk_plsig1')
   if ii(0) ge 0 then WIDGET_CONTROL, state.ppdisk_plsig1_inp, get_value=dum
   params.par(ii(0)).value = dum
  ii = where(strcompress(params.par.name, /remove_all) eq 'ppdisk_hrdisk')
   if ii(0) ge 0 then WIDGET_CONTROL, state.ppdisk_hrdisk_inp, get_value=dum
   params.par(ii(0)).value = dum
  ii = where(strcompress(params.par.name, /remove_all) eq 'ppdisk_plh')
   if ii(0) ge 0 then WIDGET_CONTROL, state.ppdisk_plh_inp, get_value=dum
   params.par(ii(0)).value = dum
  ii = where(strcompress(params.par.name, /remove_all) eq 'ppdisk_dustopac_file')
   if ii(0) ge 0 then WIDGET_CONTROL, state.ppdisk_dustopac_inp, get_value=dum
   params.par(ii(0)).value = dum

;
; Circumstellar envelope
;
  ii = where(strcompress(params.par.name,/remove_all) eq 'envelope_enable')
  WIDGET_CONTROL, state.envelope_select, get_value=dum
  if dum eq 0 then params.par(ii(0)).value = '1'
  if dum eq 1 then params.par(ii(0)).value = '0'

  
  ii = where(strcompress(params.par.name, /remove_all) eq 'envelope_rin')
   if ii(0) ge 0 then WIDGET_CONTROL, state.envelope_rin_inp, get_value=dum
   params.par(ii(0)).value = dum
  ii = where(strcompress(params.par.name, /remove_all) eq 'envelope_rout')
   if ii(0) ge 0 then WIDGET_CONTROL, state.envelope_rout_inp, get_value=dum
   params.par(ii(0)).value = dum
  ii = where(strcompress(params.par.name, /remove_all) eq 'envelope_rho0')
   if ii(0) ge 0 then WIDGET_CONTROL, state.envelope_rho0_inp, get_value=dum
   params.par(ii(0)).value = dum
  ii = where(strcompress(params.par.name, /remove_all) eq 'envelope_plrho')
   if ii(0) ge 0 then WIDGET_CONTROL, state.envelope_plrho_inp, get_value=dum
   params.par(ii(0)).value = dum
  ii = where(strcompress(params.par.name, /remove_all) eq 'envelope_cavrad')
   if ii(0) ge 0 then WIDGET_CONTROL, state.envelope_cavrad_inp, get_value=dum
   params.par(ii(0)).value = dum
  ii = where(strcompress(params.par.name, /remove_all) eq 'envelope_cavrfact')
   if ii(0) ge 0 then WIDGET_CONTROL, state.envelope_cavrfact_inp, get_value=dum
   params.par(ii(0)).value = dum
  ii = where(strcompress(params.par.name, /remove_all) eq 'envelope_dustopac_file')
   if ii(0) ge 0 then WIDGET_CONTROL, state.envelope_dustopac_inp, get_value=dum
   params.par(ii(0)).value = dum
;
; Background density
;

  ii = where(strcompress(params.par.name,/remove_all) eq 'bg_enable')
  WIDGET_CONTROL, state.bg_select, get_value=dum
  if dum eq 0 then params.par(ii(0)).value = '1'
  if dum eq 1 then params.par(ii(0)).value = '0'


  ii = where(strcompress(params.par.name, /remove_all) eq 'bg_rho')
   if ii(0) ge 0 then WIDGET_CONTROL, state.bg_rho_inp, get_value=dum
   params.par(ii(0)).value = dum
  ii = where(strcompress(params.par.name, /remove_all) eq 'bg_dustopac_file')
   if ii(0) ge 0 then WIDGET_CONTROL, state.bg_dustopac_inp, get_value=dum
   params.par(ii(0)).value = dum
;
; Spatial grid
;
  
  ii = where(strcompress(params.par.name, /remove_all) eq 'grid_rin')
   if ii(0) ge 0 then WIDGET_CONTROL, state.spatial_grid_rin_inp, get_value=dum
   params.par(ii(0)).value = dum
  ii = where(strcompress(params.par.name, /remove_all) eq 'grid_rout')
   if ii(0) ge 0 then WIDGET_CONTROL, state.spatial_grid_rout_inp, get_value=dum
   params.par(ii(0)).value = dum
  ii = where(strcompress(params.par.name, /remove_all) eq 'grid_nr')
   if ii(0) ge 0 then WIDGET_CONTROL, state.spatial_grid_nr_inp, get_value=dum
   params.par(ii(0)).value = dum
  ii = where(strcompress(params.par.name, /remove_all) eq 'grid_tmax')
   if ii(0) ge 0 then WIDGET_CONTROL, state.spatial_grid_theta_max_inp, get_value=dum
   params.par(ii(0)).value = dum
  ii = where(strcompress(params.par.name, /remove_all) eq 'grid_nt')
   if ii(0) ge 0 then WIDGET_CONTROL, state.spatial_grid_nr_theta_inp, get_value=dum
   params.par(ii(0)).value = dum

;
; Wavelength grid
;
  
  WIDGET_CONTROL, state.wav_grid_lam_bound_inp, get_value=dum
  txt = strsplit(dum, ',', /extract)
  ii = where(strcompress(params.par.name, /remove_all) eq 'lambda1')
  params.par(ii(0)).value = txt(0)
  ii = where(strcompress(params.par.name, /remove_all) eq 'lambda2')
  params.par(ii(0)).value = txt(1)
  ii = where(strcompress(params.par.name, /remove_all) eq 'lambda3')
  params.par(ii(0)).value = txt(2)
  ii = where(strcompress(params.par.name, /remove_all) eq 'lambda4')
  params.par(ii(0)).value = txt(3)
  
  WIDGET_CONTROL, state.wav_grid_nr_inp, get_value=dum
  txt = strsplit(dum, ',', /extract)
  ii = where(strcompress(params.par.name, /remove_all) eq 'n12')
  params.par(ii(0)).value = txt(0)
  ii = where(strcompress(params.par.name, /remove_all) eq 'n23')
  params.par(ii(0)).value = txt(1)
  ii = where(strcompress(params.par.name, /remove_all) eq 'n34')
  params.par(ii(0)).value = txt(2)
  
end

;**********************************************************************************************
; Procedure to fill up the widgets with the values read from problem_params.pro 
;     Attila Juhasz -  Heidelberg - 04.01.2010
;**********************************************************************************************
pro fillup_all_widget_values, state, params

;
; Thermal Monte-Carlo
;

  ii = where(strcompress(params.par.name, /remove_all) eq 'nphot')
   if ii(0) ge 0 then WIDGET_CONTROL, state.run_nr_phot_inp, set_value=params.par(ii(0)).value
  
;
; Central star
;
  ii = where(strcompress(params.par.name,/remove_all) eq 'star_enable')
  if strcompress(params.par(ii(0)).value, /remove_all) eq '1' then begin
     WIDGET_CONTROL, state.star_source_select, set_value=0
     WIDGET_CONTROL, state.star_source_inp_base, sensitive=1
  endif
  if strcompress(params.par(ii(0)).value, /remove_all) eq '0' then begin
     WIDGET_CONTROL, state.star_source_select, set_value=1
     WIDGET_CONTROL, state.star_source_inp_base, sensitive=0
  endif

  ii = where(strcompress(params.par.name, /remove_all) eq 'tstar')
   if ii(0) ge 0 then WIDGET_CONTROL, state.star_source_temp_inp, set_value=params.par(ii(0)).value
  ii = where(strcompress(params.par.name, /remove_all) eq 'rstar')
   if ii(0) ge 0 then WIDGET_CONTROL, state.star_source_rad_inp, set_value=params.par(ii(0)).value
  ii = where(strcompress(params.par.name, /remove_all) eq 'mstar')
   if ii(0) ge 0 then WIDGET_CONTROL, state.star_source_mass_inp, set_value=params.par(ii(0)).value
  
;
; Continuous starlike radiation source
;
  ii = where(strcompress(params.par.name,/remove_all) eq 'cstar_enable')
  if strcompress(params.par(ii(0)).value, /remove_all) eq '1' then begin
     WIDGET_CONTROL, state.cstar_source_select, set_value=0
     WIDGET_CONTROL, state.cstar_source_inp_base, sensitive=1
  endif
  if strcompress(params.par(ii(0)).value, /remove_all) eq '0' then begin
     WIDGET_CONTROL, state.cstar_source_select, set_value=1
     WIDGET_CONTROL, state.cstar_source_inp_base, sensitive=0
  endif

  ii = where(strcompress(params.par.name, /remove_all) eq 'cstar_rin')
   if ii(0) ge 0 then WIDGET_CONTROL, state.cstar_source_rin_inp, set_value=params.par(ii(0)).value
  ii = where(strcompress(params.par.name, /remove_all) eq 'cstar_rout')
   if ii(0) ge 0 then WIDGET_CONTROL, state.cstar_source_rout_inp, set_value=params.par(ii(0)).value
  ii = where(strcompress(params.par.name, /remove_all) eq 'cstar_trin')
   if ii(0) ge 0 then WIDGET_CONTROL, state.cstar_source_tin_inp, set_value=params.par(ii(0)).value
  ii = where(strcompress(params.par.name, /remove_all) eq 'cstar_powex')
   if ii(0) ge 0 then WIDGET_CONTROL, state.cstar_source_powex_inp, set_value=params.par(ii(0)).value

;
; External radiation source
;
  ii = where(strcompress(params.par.name,/remove_all) eq 'ext_enable')
  if strcompress(params.par(ii(0)).value, /remove_all) eq '1' then begin
     WIDGET_CONTROL, state.ext_source_select, set_value=0
     WIDGET_CONTROL, state.ext_source_inp_base, sensitive=1
  endif
  if strcompress(params.par(ii(0)).value, /remove_all) eq '0' then begin
     WIDGET_CONTROL, state.ext_source_select, set_value=1
     WIDGET_CONTROL, state.ext_source_inp_base, sensitive=0
  endif

  ii = where(strcompress(params.par.name, /remove_all) eq 'ext_temp')
   if ii(0) ge 0 then WIDGET_CONTROL, state.ext_source_temp_inp, set_value=params.par(ii(0)).value
  
;
; Circumstellar disk
;
  ii = where(strcompress(params.par.name,/remove_all) eq 'ppdisk_enable')
  if strcompress(params.par(ii(0)).value, /remove_all) eq '1' then begin
     WIDGET_CONTROL, state.ppdisk_select, set_value=0
     WIDGET_CONTROL, state.ppdisk_inp_base, sensitive=1
  endif
  if strcompress(params.par(ii(0)).value, /remove_all) eq '0' then begin
     WIDGET_CONTROL, state.ppdisk_select, set_value=1
     WIDGET_CONTROL, state.ppdisk_inp_base, sensitive=0
  endif

  ii = where(strcompress(params.par.name, /remove_all) eq 'ppdisk_rin')
   if ii(0) ge 0 then WIDGET_CONTROL, state.ppdisk_rin_inp, set_value=params.par(ii(0)).value
  ii = where(strcompress(params.par.name, /remove_all) eq 'ppdisk_rout')
   if ii(0) ge 0 then WIDGET_CONTROL, state.ppdisk_rout_inp, set_value=params.par(ii(0)).value
  ii = where(strcompress(params.par.name, /remove_all) eq 'ppdisk_sig0')
   if ii(0) ge 0 then WIDGET_CONTROL, state.ppdisk_sig0_inp, set_value=params.par(ii(0)).value
  ii = where(strcompress(params.par.name, /remove_all) eq 'ppdisk_plsig1')
   if ii(0) ge 0 then WIDGET_CONTROL, state.ppdisk_plsig1_inp, set_value=params.par(ii(0)).value
  ii = where(strcompress(params.par.name, /remove_all) eq 'ppdisk_hrdisk')
   if ii(0) ge 0 then WIDGET_CONTROL, state.ppdisk_hrdisk_inp, set_value=params.par(ii(0)).value
  ii = where(strcompress(params.par.name, /remove_all) eq 'ppdisk_plh')
   if ii(0) ge 0 then WIDGET_CONTROL, state.ppdisk_plh_inp, set_value=params.par(ii(0)).value
  ii = where(strcompress(params.par.name, /remove_all) eq 'ppdisk_dustopac_file')
   if ii(0) ge 0 then WIDGET_CONTROL, state.ppdisk_dustopac_inp, set_value=params.par(ii(0)).value

;
; Circumstellar envelope
;
  ii = where(strcompress(params.par.name,/remove_all) eq 'envelope_enable')
  if strcompress(params.par(ii(0)).value, /remove_all) eq '1' then begin
     WIDGET_CONTROL, state.envelope_select, set_value=0
     WIDGET_CONTROL, state.envelope_inp_base, sensitive=1
  endif
  if strcompress(params.par(ii(0)).value, /remove_all) eq '0' then begin
     WIDGET_CONTROL, state.envelope_select, set_value=1
     WIDGET_CONTROL, state.envelope_inp_base, sensitive=0
  endif

  ii = where(strcompress(params.par.name, /remove_all) eq 'envelope_rin')
   if ii(0) ge 0 then WIDGET_CONTROL, state.envelope_rin_inp, set_value=params.par(ii(0)).value
  ii = where(strcompress(params.par.name, /remove_all) eq 'envelope_rout')
   if ii(0) ge 0 then WIDGET_CONTROL, state.envelope_rout_inp, set_value=params.par(ii(0)).value
  ii = where(strcompress(params.par.name, /remove_all) eq 'envelope_rho0')
   if ii(0) ge 0 then WIDGET_CONTROL, state.envelope_rho0_inp, set_value=params.par(ii(0)).value
  ii = where(strcompress(params.par.name, /remove_all) eq 'envelope_plrho')
   if ii(0) ge 0 then WIDGET_CONTROL, state.envelope_plrho_inp, set_value=params.par(ii(0)).value
  ii = where(strcompress(params.par.name, /remove_all) eq 'envelope_cavrad')
   if ii(0) ge 0 then WIDGET_CONTROL, state.envelope_cavrad_inp, set_value=params.par(ii(0)).value
  ii = where(strcompress(params.par.name, /remove_all) eq 'envelope_cavrfact')
   if ii(0) ge 0 then WIDGET_CONTROL, state.envelope_cavrfact_inp, set_value=params.par(ii(0)).value
  ii = where(strcompress(params.par.name, /remove_all) eq 'envelope_dustopac_file')
   if ii(0) ge 0 then WIDGET_CONTROL, state.envelope_dustopac_inp, set_value=params.par(ii(0)).value
;
; Background density
;
  ii = where(strcompress(params.par.name,/remove_all) eq 'bg_enable')
  if strcompress(params.par(ii(0)).value, /remove_all) eq '1' then begin
     WIDGET_CONTROL, state.bg_select, set_value=0
     WIDGET_CONTROL, state.bg_inp_base, sensitive=1
  endif
  if strcompress(params.par(ii(0)).value, /remove_all) eq '0' then begin
     WIDGET_CONTROL, state.bg_select, set_value=1
     WIDGET_CONTROL, state.bg_inp_base, sensitive=0
  endif

  ii = where(strcompress(params.par.name, /remove_all) eq 'bg_rho')
   if ii(0) ge 0 then WIDGET_CONTROL, state.bg_rho_inp, set_value=params.par(ii(0)).value
  ii = where(strcompress(params.par.name, /remove_all) eq 'bg_dustopac_file')
   if ii(0) ge 0 then WIDGET_CONTROL, state.bg_dustopac_inp, set_value=params.par(ii(0)).value
;
; Spatial grid
;
  
  ii = where(strcompress(params.par.name, /remove_all) eq 'grid_rin')
   if ii(0) ge 0 then WIDGET_CONTROL, state.spatial_grid_rin_inp, set_value=params.par(ii(0)).value
  ii = where(strcompress(params.par.name, /remove_all) eq 'grid_rout')
   if ii(0) ge 0 then WIDGET_CONTROL, state.spatial_grid_rout_inp, set_value=params.par(ii(0)).value
  ii = where(strcompress(params.par.name, /remove_all) eq 'grid_nr')
   if ii(0) ge 0 then WIDGET_CONTROL, state.spatial_grid_nr_inp, set_value=params.par(ii(0)).value
  ii = where(strcompress(params.par.name, /remove_all) eq 'grid_tmax')
   if ii(0) ge 0 then WIDGET_CONTROL, state.spatial_grid_theta_max_inp, set_value=params.par(ii(0)).value
  ii = where(strcompress(params.par.name, /remove_all) eq 'grid_nt')
   if ii(0) ge 0 then WIDGET_CONTROL, state.spatial_grid_nr_theta_inp, set_value=params.par(ii(0)).value
  
;
; Wavelength grid
;
  
  ii = where(strcompress(params.par.name, /remove_all) eq 'lambda1')
  txt = strcompress(params.par(ii(0)).value, /remove_all)
  ii = where(strcompress(params.par.name, /remove_all) eq 'lambda2')
  txt = txt +', '+ strcompress(params.par(ii(0)).value, /remove_all)
  ii = where(strcompress(params.par.name, /remove_all) eq 'lambda3')
  txt = txt +', '+ strcompress(params.par(ii(0)).value, /remove_all)
  ii = where(strcompress(params.par.name, /remove_all) eq 'lambda4')
  txt = txt +', '+ strcompress(params.par(ii(0)).value, /remove_all)
  WIDGET_CONTROL, state.wav_grid_lam_bound_inp, set_value=txt
 
  ii = where(strcompress(params.par.name, /remove_all) eq 'n12')
  txt = strcompress(params.par(ii(0)).value, /remove_all)
  ii = where(strcompress(params.par.name, /remove_all) eq 'n23')
  txt = txt +', '+ strcompress(params.par(ii(0)).value, /remove_all)
  ii = where(strcompress(params.par.name, /remove_all) eq 'n34')
  txt = txt +', '+ strcompress(params.par(ii(0)).value, /remove_all)
  WIDGET_CONTROL, state.wav_grid_nr_inp, set_value=txt

end

;**********************************************************************************************
; Procedure to re-plot the whole disk structure window (1D plot)
;     Attila Juhasz -  Heidelberg - 05.01.2010
;**********************************************************************************************
pro replot_1d_disk_struct, _EXTRA=_EXTRA

COMMON DISK_STRUCTURE, struct, tau, cont_line_var, cont_fill_var, xrange_2d, yrange_2d, xlog, ylog,$
                       plot_var, plot_crd, plot_crd_id, xrange_1d, yrange_1d, plot_dim, plot_crd_sys


au = 1.496d13

;
; Plot physical quantities as a function of radius
;
if plot_crd eq 0 then begin
;
; Make the titles of the axes
;
   xtitle = 'R [AU]'
   title  = '!4h!X = '+strcompress((!pi/2.-struct.grid.theta(plot_crd_id)), /remove_all)+' rad'
   
   case plot_var of 
;
; Plot the temperature
;
      0 : begin
         ytitle = 'Temperature [K]'
         plot, struct.grid.r/au, struct.temp(*, plot_crd_id), xrange=xrange_1d, /xs, xlog=xlog, ylog=ylog, $
               xtitle=xtitle, ytitle=ytitle, title=title, psym=-4, _EXTRA=_EXTRA
      end
;
; Plot the density
;
      1 : begin
         ytitle = 'Density [g/cm^3]'
         plot, struct.grid.r/au, struct.rho(*, plot_crd_id), xrange=xrange_1d, /xs, xlog=xlog, ylog=ylog, $
               xtitle=xtitle, ytitle=ytitle, title=title, psym=-4, _EXTRA=_EXTRA
      end
;
; Plot the radial optical depth
;
      2 : begin
         ytitle = 'Radial optical depth'
         plot, tau.r/au, tau.tau.taur(*, plot_crd_id), xrange=xrange_1d, /xs, xlog=xlog, ylog=ylog, $
               xtitle=xtitle, ytitle=ytitle, title=title, psym=-4, _EXTRA=_EXTRA
      end
;
; Plot the vertical optical depth
;
      3 : begin
         ytitle = 'Vertical optical depth'
         plot, tau.r/au, tau.tau.taut(*, plot_crd_id), xrange=xrange_1d, /xs, xlog=xlog, ylog=ylog, $
               xtitle=xtitle, ytitle=ytitle, title=title, psym=-4, _EXTRA=_EXTRA
      end
   endcase
endif else begin
;
; Plot physical quantities as a function of the meridional angle (theta)
;
;
; Make the titles of the axes
;
   xtitle = '!4p!X/2-!4h!X [radian] !9A!X z/R'
   title  = 'R = '+strcompress(struct.grid.r(plot_crd_id)/au, /remove_all)+' AU'
   case plot_var of 
;
; Plot the temperature
;
      0 : begin
         ytitle = 'Temperature [K]'
         plot, !pi/2.-struct.grid.theta, struct.temp(plot_crd_id, *), xrange=xrange_1d, /xs, xlog=xlog, ylog=ylog, $
               xtitle=xtitle, ytitle=ytitle, title=title, psym=-4, _EXTRA=_EXTRA
      end
;
; Plot the density
;
      1 : begin
         ytitle = 'Density [g/cm^3]'
         plot, !pi/2.-struct.grid.theta, struct.rho(plot_crd_id, *), xrange=xrange_1d, /xs, xlog=xlog, ylog=ylog, $
               xtitle=xtitle, ytitle=ytitle, title=title, psym=-4, _EXTRA=_EXTRA
      end
;
; Plot the radial optical depth
;
      2 : begin
         ytitle = 'Radial optical depth'
         plot, !pi/2.-tau.theta, tau.tau.taur(plot_crd_id, *), xrange=xrange_1d, /xs, xlog=xlog, ylog=ylog, $
               xtitle=xtitle, ytitle=ytitle, title=title, psym=-4, _EXTRA=_EXTRA
      end
;
; Plot the vertical optical depth
;
      3 : begin
         ytitle = 'Vertical optical depth'
         plot, !pi/2.-tau.theta, tau.tau.taut(plot_crd_id, *), xrange=xrange_1d, /xs, xlog=xlog, ylog=ylog, $
               xtitle=xtitle, ytitle=ytitle, title=title, psym=-4, _EXTRA=_EXTRA
      end
   endcase

endelse

end
;**********************************************************************************************
; Procedure to re-plot the whole disk structure window (2D contours)
;     Attila Juhasz -  Heidelberg - 05.01.2010
;**********************************************************************************************
pro replot_2d_disk_struct, _EXTRA=_EXTRA
  
COMMON DISK_STRUCTURE, struct, tau, cont_line_var, cont_fill_var, xrange_2d, yrange_2d, xlog, ylog,$
                       plot_var, plot_crd, plot_crd_id, xrange_1d, yrange_1d, plot_dim, plot_crd_sys

au = 1.496d13

;
; Check if the temperature is known (i.e. radmc3d has already run)
;

if n_elements(struct.temp) eq 1 then temp_def = 0
if n_elements(struct.temp) gt 1 then temp_def = 1

;
; Normal descartes contours (z/r vs r)
;
if plot_crd_sys eq 0 then begin
;
; Make the titles of the axes
;
   xtitle = 'R [AU]'
   ytitle = '!4p!X/2 - !4h!X !9A!X z/R'

;
; Plot the filled contour background
;
   case cont_fill_var of
      -1 : begin
         plot, [-10000], [-10000], xlog=xlog, ylog=ylog, xrange=xrange_2d, yrange=yrange_2d, /xs, /ys, $
               xtitle=xtitle, ytitle=ytitle
      end
      
      0 : begin
         loadct, 3
         col = findgen(50) * 255/49
         if temp_def eq 1 then begin
            contour, struct.temp, struct.grid.r/au, !pi/2.-struct.grid.theta, nlev=50, c_col=col, /fill, $
                     xlog=xlog, ylog=ylog, xrange=xrange_2d, yrange=yrange_2d, /xs, /ys, xtitle=xtitle, ytitle=ytitle, $
                     charsize=1.0, _EXTRA=_EXTRA
         endif else begin
            plot, [-10000], [-10000], xlog=xlog, ylog=ylog, xrange=xrange_2d, yrange=yrange_2d, /xs, /ys, $
                  xtitle=xtitle, ytitle=ytitle, charsize=1.0, _EXTRA=_EXTRA
         endelse
      end
      
      1 : begin
         loadct, 1
         col = findgen(50) * 255/49
         contour, alog10(struct.rho), struct.grid.r/au, !pi/2.-struct.grid.theta, nlev=50, c_col=col, /fill, $
                  xlog=xlog, ylog=ylog, xrange=xrange_2d, yrange=yrange_2d, /xs, /ys, xtitle=xtitle, ytitle=ytitle, $
                  charsize=1.0, _EXTRA=_EXTRA
      end
      
      2 : begin
         loadct, 8
         col = findgen(50) * 255/49
         contour, tau.tau.taur, tau.r/au, !pi/2.-tau.theta, nlev=50, /fill, $
                  xlog=xlog, ylog=ylog, xrange=xrange_2d, yrange=yrange_2d, /xs, /ys, xtitle=xtitle, ytitle=ytitle, $
                  charsize=1.0, _EXTRA=_EXTRA
      end
      
      3 : begin
         loadct, 7
         col = findgen(50) * 255/49
         contour, tau.tau.taut, tau.r/au, !pi/2.-tau.theta, nlev=50, /fill, $
                  xlog=xlog, ylog=ylog, xrange=xrange_2d, yrange=yrange_2d, /xs, /ys, xtitle=xtitle, ytitle=ytitle, $
                  charsize=1.0, _EXTRA=_EXTRA
      end
      
   endcase
   loadct, 39
   
;
; Plot the variables with contour lines
;

   if cont_line_var(0) eq 1 then begin
      if cont_fill_var eq 0 then begin
         col = fltarr(12) + 255
      endif else begin
         col = fltarr(12) + 254
      endelse
      
      if temp_def eq 1 then begin
         lev = dindgen(12) * (max(struct.temp) - min(struct.temp)) / 12d0  + min(struct.temp)
         contour, struct.temp, struct.grid.r/au, !pi/2.-struct.grid.theta, lev=lev, c_label=intarr(12)+1, $
                  c_col=col, xlog=xlog, ylog=ylog, xrange=xrange_2d, yrange=yrange_2d, /overplot, /xs, /ys
      endif
   endif
   
   if cont_line_var(1) eq 1 then begin
      if cont_fill_var eq 1 then begin
         col = fltarr(12) + 255
      endif else begin
         col = fltarr(12) + 70
      endelse
      
      lev = dindgen(12) * (max(alog10(struct.rho)) - min(alog10(struct.rho))) / 12d0  + min(alog10(struct.rho))
      contour, alog10(struct.rho), struct.grid.r/au, !pi/2.-struct.grid.theta, lev=lev, c_label=intarr(12)+1, $
               c_col=col, xlog=xlog, ylog=ylog, xrange=xrange_2d, yrange=yrange_2d, /overplot, /xs, /ys
   endif
   
   if cont_line_var(2) eq 1 then begin
      if cont_fill_var eq 2 then begin
         col = fltarr(12) + 255
      endif else begin
         col = fltarr(12) + 130
      endelse
      
      lev = [dindgen(12) * (max(tau.tau.taur) - min(tau.tau.taur)) / 12d0  + min(tau.tau.taur), 1.0]
      id = sort(lev)
      lev = lev(id)
      contour, tau.tau.taur, tau.r/au, !pi/2.-tau.theta, lev=lev, c_label=intarr(12)+1, $
               c_col=col, xlog=xlog, ylog=ylog, xrange=xrange_2d, yrange=yrange_2d, /overplot, /xs, /ys
   endif

   if cont_line_var(3) eq 1 then begin
      if cont_fill_var eq 3 then begin
         col = fltarr(12) + 255
      endif else begin
      col = fltarr(12) + 160
   endelse
      
        
      loadct, 2
      lev = [dindgen(12) * (max(tau.tau.taut) - min(tau.tau.taut)) / 12d0  + min(tau.tau.taut), 1.0]
      id = sort(lev)
      lev = lev(id)
      contour, tau.tau.taut, tau.r/au, !pi/2.-tau.theta, lev=lev, c_label=intarr(12)+1, $
               c_col=col, xlog=xlog, ylog=ylog, xrange=xrange_2d, yrange=yrange_2d, /overplot, /xs, /ys
      loadct, 39
   endif

endif

;
; Polar contours
;

if plot_crd_sys eq 1 then begin
;
; Make the titles of the axes
;
   xtitle = 'R [AU]'
   ytitle = 'Z [AU]'

;
; Plot the filled contour background
;
   case cont_fill_var of
      0 : begin
         loadct, 3
         col = findgen(50) * 255/49

         if temp_def eq 1 then begin
            polar_contour, transpose(struct.temp), !pi/2.-struct.grid.theta, struct.grid.r/au, nlev=50, c_col=col, /fill, $
                           xlog=xlog, ylog=ylog, xrange=xrange_2d, yrange=yrange_2d, /xs, /ys, xtitle=xtitle, ytitle=ytitle, $
                           charsize=1.0, /dither, _EXTRA=_EXTRA
         endif else begin
            plot, [-10000], [-10000], xlog=xlog, ylog=ylog, xrange=xrange_2d, yrange=yrange_2d, /xs, /ys, $
                  xtitle=xtitle, ytitle=ytitle, charsize=1.0, _EXTRA=_EXTRA
         endelse
      end
      
      1 : begin
         loadct, 1
         col = findgen(50) * 255/49
         polar_contour, alog10(transpose(struct.rho)), !pi/2.-struct.grid.theta, struct.grid.r/au, nlev=50, c_col=col, /fill, $
                  xlog=xlog, ylog=ylog, xrange=xrange_2d, yrange=yrange_2d, /xs, /ys, xtitle=xtitle, ytitle=ytitle, $
                  charsize=1.0, /dither, _EXTRA=_EXTRA
      end
      
      2 : begin
         loadct, 8
         col = findgen(50) * 255/49
         polar_contour, transpose(tau.tau.taur), !pi/2.-tau.theta, tau.r/au, nlev=50, /fill, $
                  xlog=xlog, ylog=ylog, xrange=xrange_2d, yrange=yrange_2d, /xs, /ys, xtitle=xtitle, ytitle=ytitle, $
                  charsize=1.0, /dither, _EXTRA=_EXTRA
      end
      
      3 : begin
         loadct, 7
         col = findgen(50) * 255/49
         polar_contour, transpose(tau.tau.taut), !pi/2.-tau.theta, tau.r/au, nlev=50, /fill, $
                  xlog=xlog, ylog=ylog, xrange=xrange_2d, yrange=yrange_2d, /xs, /ys, xtitle=xtitle, ytitle=ytitle, $
                  charsize=1.0, /dither, _EXTRA=_EXTRA
      end
      
   endcase
   loadct, 39
   
;
; Plot the variables with contour lines
;

   if cont_line_var(0) eq 1 then begin
      if cont_fill_var eq 0 then begin
         col = fltarr(12) + 255
      endif else begin
         col = fltarr(12) + 254
      endelse
      
      if temp_def eq 1 then begin
         lev = dindgen(12) * (max(struct.temp) - min(struct.temp)) / 12d0  + min(struct.temp)
         polar_contour, transpose(struct.temp), !pi/2.-struct.grid.theta, struct.grid.r/au, lev=lev, c_label=intarr(12)+1, $
                        c_col=col, xlog=xlog, ylog=ylog, xrange=xrange_2d, yrange=yrange_2d, /overplot, /xs, /ys, /dither
      endif
   endif
   
   if cont_line_var(1) eq 1 then begin
      if cont_fill_var eq 1 then begin
         col = fltarr(12) + 255
      endif else begin
         col = fltarr(12) + 70
      endelse
      
      lev = dindgen(12) * (max(alog10(struct.rho)) - min(alog10(struct.rho))) / 12d0  + min(alog10(struct.rho))
      polar_contour, transpose(alog10(struct.rho)), !pi/2.-struct.grid.theta, struct.grid.r/au, lev=lev, c_label=intarr(12)+1, $
               c_col=col, xlog=xlog, ylog=ylog, xrange=xrange_2d, yrange=yrange_2d, /overplot, /xs, /ys, /dither
   endif
   
   if cont_line_var(2) eq 1 then begin
      if cont_fill_var eq 2 then begin
         col = fltarr(12) + 255
      endif else begin
         col = fltarr(12) + 130
      endelse
      
      lev = [dindgen(12) * (max(tau.tau.taur) - min(tau.tau.taur)) / 12d0  + min(tau.tau.taur), 1.0]
      id = sort(lev)
      lev = lev(id)
      polar_contour, transpose(tau.tau.taur), !pi/2.-tau.theta, tau.r/au, lev=lev, c_label=intarr(12)+1, $
               c_col=col, xlog=xlog, ylog=ylog, xrange=xrange_2d, yrange=yrange_2d, /overplot, /xs, /ys, /dither
   endif

   if cont_line_var(3) eq 1 then begin
      if cont_fill_var eq 3 then begin
         col = fltarr(12) + 255
      endif else begin
      col = fltarr(12) + 160
   endelse
      
        
      loadct, 2
      lev = [dindgen(12) * (max(tau.tau.taut) - min(tau.tau.taut)) / 12d0  + min(tau.tau.taut), 1.0]
      id = sort(lev)
      lev = lev(id)
      polar_contour, transpose(tau.tau.taut), !pi/2.-tau.theta, tau.r/au, lev=lev, c_label=intarr(12)+1, $
               c_col=col, xlog=xlog, ylog=ylog, xrange=xrange_2d, yrange=yrange_2d, /overplot, /xs, /ys, /dither
      loadct, 39
   endif

endif
end

;**********************************************************************************************
; Procedure to plot mass absorption coefficients
;     Attila Juhasz -  Heidelberg - 07.01.2010
;**********************************************************************************************
pro plot_kappa, _EXTRA=_EXTRA
  
  COMMON SPECTRUM, spec, xunit, yunit, xrange_sed, yrange_sed, dist, plot_sed_kappa, opac
  
  plot, opac.wave, opac.cabs, /xl, /yl, xrange=xrange_sed, yrange=yrange, /xs, /ys, $
        psym=-4, symsize=0.3, xtitle='!4k!X [!4l!Xm]', ytitle='!4j!X [cm^2/g]', _EXTRA=_EXTRA
 
  oplot, opac.wave, opac.cabs, psym=-4, symsize=0.3, col=70

  ii = where(opac.csca gt 0.)
  
  if ii(0) ge 0 then begin
     oplot, opac.wave, opac.csca, psym=-4, symsize=0.3, col=254
     ;; Make a legend
     plots, [0.8, 0.9], [0.86, 0.86], psym=-4, symsize=0.3, col=70, /normal
     plots, [0.8, 0.9], [0.8, 0.8], psym=-4, symsize=0.3, col=254, /normal
     xyouts, 0.76, 0.85, '!4j!X!DABS!N', align=0.5, /normal, charsize=1.6
     xyouts, 0.76, 0.79, '!4j!X!DSCA!N', align=0.5, /normal, charsize=1.6
  endif else begin
     ;; Make a legend
     plots, [0.8, 0.9], [0.86, 0.86], psym=-4, symsize=0.3, col=70, /normal
     xyouts, 0.76, 0.85, '!4j!X!DABS!N', align=0.5, /normal, charsize=1.6
  endelse
end

;**********************************************************************************************
; Procedure to re-plot the SED window
;     Attila Juhasz -  Heidelberg - 07.01.2010
;**********************************************************************************************
pro replot_sed, _EXTRA=_EXTRA

    COMMON SPECTRUM, spec, xunit, yunit, xrange_sed, yrange_sed, dist, plot_sed_kappa, opac
   
;;
;; Choose the unit of the x-axis
;;
    case xunit of 
       ;; Micrometer
       0 : begin
          x = spec.lambda
          xtitle = '!4k!X [!4l!Xm]'
       end
       ;; Hertz
       1 : begin
          x = spec.freq
          xtitle = '!4m!X [Hz]'
       end
    endcase

;;
;; Choose the unit of the y-axis
;;
    case yunit of 
       ;; Fnu [erg/s/cm/cm/Hz]
       0 : begin
          y = spec.spectrum
          ytitle = 'F!D!4m!X!N [erg/s/cm/cm/Hz]'
       end
       ;; nuFnu [erg/s/cm/cm/Hz]
       1 : begin
          y = spec.spectrum * spec.freq
          ytitle = '!4m!XF!D!4m!X!N [erg/s/cm/cm]'
       end
       ;; Jy [1d-23*erg/s/cm/cm/Hz]
       2 : begin
          y = spec.spectrum * 1d23
          ytitle = 'F!D!4m!X!N [Jy]'
       end

    endcase
;;
;; Scale the SED according to the distance
;;
    dfact = 1./dist^2.
    
;;
;; Make the plot
;;

    plot, x, y*dfact, /xl, /yl, xtitle=xtitle, ytitle=ytitle, xrange=xrange_sed, yrange=yrange_sed, $
          /xs, /ys, _EXTRA=_EXTRA

end

;**********************************************************************************************
; Display help messages
;     Attila Juhasz -  Heidelberg - 05.01.2010
;**********************************************************************************************
pro display_help, help_id=help_id

  COMMON PARAMETERS, params
  
  help_main_base    = WIDGET_BASE(/column, title='HELP RADMC3D GUI V0.01', uname='main_window')
     
  help_text         = WIDGET_TEXT(help_main_base, uvalue='help_text', scr_xsize=480, scr_ysize=600, $
                                  /scroll, value=['EMPTY'], /align_center)
  help_close_base   = WIDGET_BASE(help_main_base, /row, /align_center)
  help_close_button = WIDGET_BUTTON(help_close_base, uvalue='help_close_button', $
                                    xsize=100, value='Close')

bspace = " "
pname_len0 = 30

;;
;; Help text for general help
;;
  if help_id eq 0 then begin
     help = [' RADMC3D GUI V0.01 - TEST VERSION ', $
             '',$
             'This Graphical User Interface provides an easy way to run and',$
             'analyze the results of RADMC3D. Some more description to come...']
             
     widget_control, help_text, set_value = get_started_help_text
  endif
;;
;; Help text for help on standard built-in parameters
;; 
  if help_id eq 1 then begin
     
     help = strarr(n_elements(params.par.name))
     for ipar=0, n_elements(params.par.name)-1 do begin
        diff = pname_len0 - strlen(params.par(ipar).name) 
        if diff gt 0 then begin
           separator = strarr(diff) + bspace
           separator = strjoin(separator, /single)
        endif else begin
           separator = ' '
        endelse
        help(ipar) = params.par(ipar).name + separator + $
                     params.par(ipar).comment
     endfor
     
;;
;; Help text on extra parameters
;;     
     if n_elements(params.ext) gt 1 or params.ext(0).group_id gt 0 then begin
        help_extra = [' ','EXTRA PARAMETERS: ',' ']
        
        for ipar=0, n_elements(params.ext)-1 do begin
           diff = pname_len0 - strlen(params.ext(ipar).name) 
           if diff gt 0 then begin
              separator = strarr(diff) + bspace
              separator = strjoin(separator, /single)
           endif else begin
              separator = ' '
           endelse
           help_extra = [help_extra, params.ext(ipar).name + separator + $
                         params.ext(ipar).comment]
        endfor
        
        help = [help, help_extra]
     endif

     
     help = ['DEFAULT PARAMETERS : ', help]
     
  endif
;;
;; Help text for natural constants
;; 
  if help_id eq 2 then begin
     help      = ['NATURAL CONSTANTS:', $
                  ' ',$
                  ' All input parameters should be in cgs units (see RADMC3D Manual)  ',$
                  ' Other units are, however, usually more convenient (e.g. AU for    ',$
                  ' linear scales of protoplanetary disks). Therefore, there are      ',$
                  ' natural constants defined to convert the input parameters in these',$
                  ' units to cgs units',$
                  ' ',$
                  ' The GUI input text widgets handle the input as strings, so one can ',$
                  ' define any arbitrary/convenient units/expressions', $
                  ' For example : The radius of the star can be given as 6.96d10, 1d0*RS,',$
                  ' or as 0.0046524064*AU.',$
                  ' ',$
                  ' DEFINED NATURAL CONSTANTS : ',$
                  ' ',$
                  'AU  = 1.49598d13     ; Astronomical Unit       [cm]',$
                  'pc  = 3.08572d18     ; Parsec                  [cm]',$
                  'MS  = 1.98892d33     ; Solar mass              [g]',$
                  'TS  = 5.78d3         ; Solar temperature       [K]',$
                  'LS  = 3.8525d33      ; Solar luminosity        [erg/s]',$
                  'RS  = 6.96d10        ; Solar radius            [cm]']

  endif
 
  widget_control, help_text, set_value = help
  help_state = {help_text:help_text, help_close_button:help_close_button}

  widget_control, help_main_base, /realize
  widget_control, help_main_base, set_uvalue=help_state
  xmanager, 'display_help', help_main_base, /no_block
end


pro display_help_event, help_event
  
  WIDGET_CONTROL, help_event.top, get_uvalue=help_state
  WIDGET_CONTROL, help_event.id, get_uvalue=uvalue

  case uvalue of
     'help_text' : begin
     end
     'help_close_button' : begin
        widget_control, help_event.top, /destroy
     end
  endcase

end

;**********************************************************************************************
; Graphical User Interface for RADMC3D V0.01
;
;     Attila Juhasz -  Heidelberg - 04.01.2010
;**********************************************************************************************

pro radmc3d_gui 
  
COMMON PARAMETERS, params
COMMON SPECTRUM, spec, xunit, yunit, xrange_sed, yrange_sed, dist, plot_sed_kappa, opac
COMMON DISK_STRUCTURE, struct, tau, cont_line_var, cont_fill_var, xrange_2d, yrange_2d, xlog, ylog,$
                       plot_var, plot_crd, plot_crd_id, xrange_1d, yrange_1d, plot_dim, plot_crd_sys

;;
;; Read the parameters from problem_params.pro
;;
params = read_params_radmc3d(read_file_name='problem_params.pro')

;----------------------------------------------------------------------------------------------
; Create the main window
;----------------------------------------------------------------------------------------------
main_window                  = WIDGET_BASE(/row, MBAR=menubar, title=' RADMC3D GUI V0.01 - TEST VERSION', $
                                           uname='main_window', ysize=725, xsize=1273, ypad=0)

main_left_window             = WIDGET_BASE(main_window, /column, xsize=1150)
main_right_window            = WIDGET_BASE(main_window, /column, xsize=120)
 
; Base for the drawing windows and their control buttons
drawer_base                  = WIDGET_BASE(main_left_window, xsize=1000, /row, ysize=430, frame=0)
; Base for the input (below the drawing windows)
struct_base                  = WIDGET_BASE(main_left_window, xsize=1000, /row, frame=0)

;----------------------------------------------------------------------------------------------
; Create the menu
;----------------------------------------------------------------------------------------------
file_menu                    = WIDGET_BUTTON(menubar, value='File', /Menu)
file_menu_new                = WIDGET_BUTTON(file_menu, value='  New model  ', uvalue='file_menu_new')
file_menu_open               = WIDGET_BUTTON(file_menu, value='  Open model ', uvalue='file_menu_open')
file_menu_save               = WIDGET_BUTTON(file_menu, value='  Save model ', uvalue='file_menu_save')
file_menu_quit               = WIDGET_BUTTON(file_menu, value='  Quit ', uvalue='file_menu_quit')

pref_menu                    = WIDGET_BUTTON(menubar, value='Preferences', /Menu)
pref_menu_paths              = WIDGET_BUTTON(pref_menu, value='   Set paths   ', uvalue='pref_menu_paths')

help_menu                    = WIDGET_BUTTON(menubar, value='Help', /Menu)
help_menu_general_help       = WIDGET_BUTTON(help_menu, value='  Get start           ', $
                                             uvalue='help_menu_general_help') 
help_menu_parameter_help     = WIDGET_BUTTON(help_menu, value='  Help on parameters  ', $
                                             uvalue='help_menu_parameter_help') 
help_menu_natconst_help      = WIDGET_BUTTON(help_menu, value='  Help on natural constants  ', $
                                             uvalue='help_menu_natconst_help') 
;----------------------------------------------------------------------------------------------
; Create drawing widgets for spectrum/images
;----------------------------------------------------------------------------------------------
drawer_left_base             = WIDGET_BASE(drawer_base, xsize=550, /column, frame=0, ypad=6)
spec_imag_drawer             = WIDGET_DRAW(drawer_left_base, xsize=550, ysize=350, $
                                           retain=2, uvalue='spec_imag_drawer', /drop_event)
;----------------------------------------------------------------------------------------------
; Buttons for the spectrum drawing widget
;----------------------------------------------------------------------------------------------
spec_imag_par_upper_row_base = widget_base(drawer_left_base, ysize=30, /row)
spec_imag_par_lower_row_base = widget_base(drawer_left_base, ysize=30, /row)



spec_imag_xunit_label        = WIDGET_LABEL(spec_imag_par_upper_row_base, value='Xunit:', $
                                            /align_left)
spec_imag_xunit_select       = WIDGET_COMBOBOX(spec_imag_par_upper_row_base, xsize=120, $
                                               uvalue='spec_imag_xunit_select',$
                                               value=['Micron', 'Hz'], frame=1)

spec_imag_xrange_inp         = CW_FIELD(spec_imag_par_upper_row_base, title='Xrange:', $
                                        uvalue='spec_imag_xrange_inp', value=' ', /row, $
                                        xsize=12, /return_events)

spec_imag_dist_inp           = CW_FIELD(spec_imag_par_upper_row_base, title='Dist [pc]:', $
                                        uvalue='spec_imag_dist_inp', value='140.', /row, $
                                        xsize=6, /return_events)
spec_imag_plot_spec_button   = WIDGET_BUTTON(spec_imag_par_upper_row_base, value='Plot SED',$
                                             uvalue='spec_imag_plot_spec_button', xsize=90, $
                                            tooltip='Plot the SED (and reset the x-/y-ranges)')

spec_imag_yunit_label        = WIDGET_LABEL(spec_imag_par_lower_row_base, value='Yunit:', $
                                            /align_left)
spec_imag_yunit_select       = WIDGET_COMBOBOX(spec_imag_par_lower_row_base, xsize=120, $
                                               uvalue='spec_imag_yunit_select',$
                                               value=['erg/s/cm/cm/Hz', 'erg/cm/cm/s', 'Jy'], frame=1)
spec_imag_yrange_inp         = CW_FIELD(spec_imag_par_lower_row_base, title='Yrange:', $
                                        uvalue='spec_imag_yrange_inp', value=' ', /row, $
                                        xsize=12, /return_events)
spec_imag_oplot_data_button  = WIDGET_BUTTON(spec_imag_par_lower_row_base, value='Oplot data',$
                                             uvalue='spec_imag_oplot_data_button', xsize=110, $
                                             tooltip='Overplot observational data')
spec_imag_make_ps_button     = WIDGET_BUTTON(spec_imag_par_lower_row_base, value='Postscript',$
                                             uvalue='spec_imag_make_ps_button', xsize=110, $
                                             tooltip='Plot the SED/Opacity window into a postscript file')
;----------------------------------------------------------------------------------------------
; Create drawing widgets to plot/contour the disk/envelope structure
;----------------------------------------------------------------------------------------------
drawer_right_base            = WIDGET_BASE(drawer_base, xsize=590, /column, frame=0)
; Base for the y-slider and the drawing widget
drawer_right_upper_row_base  = WIDGET_BASE(drawer_right_base, xsize=580, /row, frame=0)
; Base for the x-slider
drawer_right_lower_row_base  = WIDGET_BASE(drawer_right_base, xsize=580, /row, frame=0, xpad=37)

disk_struct_yslide           = WIDGET_SLIDER(drawer_right_upper_row_base, maximum=90.0, $
                                             minimum=1.0, uvalue='disk_struct_yslide', value=90.,$
                                             /vertical, /drag)
disk_struct_drawer           = WIDGET_DRAW(drawer_right_upper_row_base, xsize=550, ysize=350, $
                                           retain=2, uvalue='disk_struct_drawer')
disk_struct_xslide           = WIDGET_SLIDER(drawer_right_lower_row_base, maximum=100.0, $
                                             minimum=1.0, uvalue='disk_struct_xslide', $
                                             value=100., scr_xsize=555, /drag)

;----------------------------------------------------------------------------------------------
; Create buttons what to plot for the disk/envelope structure
;   (Density, temperature, optical depth, make ps file)
;----------------------------------------------------------------------------------------------
plot2d_label                 = WIDGET_LABEL(main_right_window, value='2D Plot')
plot2d_column_labels         = WIDGET_LABEL(main_right_window, value='Var. Lines Fill')
plot2d_base                  = WIDGET_BASE(main_right_window, xsize=110,  /row,$
                                           /align_left, frame=1)

plot2d_label_base            = WIDGET_BASE(plot2d_base, /column, ypad=2)
plot2d_temp_label            = WIDGET_LABEL(plot2d_label_base, value='Temp.', ysize=23)
plot2d_dens_label            = WIDGET_LABEL(plot2d_label_base, value='Dens.', ysize=23)
plot2d_taur_label            = WIDGET_LABEL(plot2d_label_base, value='Taur ', ysize=23)
plot2d_taut_label            = WIDGET_LABEL(plot2d_label_base, value='Taut ', ysize=23)
plot2d_clines_button         = CW_BGROUP(plot2d_base, [' ', ' ', ' ', ' '], /nonexclusive, /column, space=0, $
                                         uvalue='plot2d_clines_button', xsize=33, set_value=[0,1,0,0])
plot2d_cfill_button          = CW_BGROUP(plot2d_base, [' ', ' ', ' ', ' '], /exclusive, /column, space=0, $
                                         uvalue='plot2d_cfill_button', set_value=1)
plot_nd_coord_sys_button     = CW_BGROUP(main_right_window, ['Theta', 'Z'], /exclusive, /row, space=0, $
                                         uvalue='plot_nd_coord_sys_button', set_value=0)


plot1d_label                 = WIDGET_LABEL(main_right_window, value='1D Plot')
plot1d_base                  = WIDGET_BASE(main_right_window, xsize=110, ysize=110, /column, $
                                           /align_left, frame=1)
plot1d_coord_button          = CW_BGROUP(plot1d_base, ['R', 'Theta'], /exclusive, /row, space=0, $
                                         label_left='', uvalue='plot1d_coord_button', set_value=0)
plot1d_var_button            = CW_BGROUP(plot1d_base, ['Temp', 'Dens', 'Taur', 'Taut'], /exclusive,  $
                                         row=2, space=0, uvalue='plot1d_var_button')

plot_nd_tauwav_base          = WIDGET_BASE(main_right_window, /row)
plot_nd_tauwav_inp           = CW_FIELD(plot_nd_tauwav_base, title='Wav of tau:', uvalue='plot1d_tauwav_inp', $
                                        xsize=10, /column, value='0.55')
plot_nd_axis_scale_button    = CW_BGROUP(main_right_window, ['xlog', 'ylog'], /nonexclusive, /row, space=0, $
                                         uvalue='plot_nd_axis_scale_button', set_value=[1,0])

plot_nd_update_base          = WIDGET_BASE(main_right_window, /row)
plot_struct_update_button    = WIDGET_BUTTON(plot_nd_update_base, value='Update',$
                                             uvalue='plot_struct_update_button', xsize=110)
plot_nd_ps_base              = WIDGET_BASE(main_right_window, /row)
plot_struct_ps_button        = WIDGET_BUTTON(plot_nd_ps_base, value='Postscript',$
                                             uvalue='plot_struct_ps_button', xsize=110)

;----------------------------------------------------------------------------------------------
; Radiation sources
;----------------------------------------------------------------------------------------------
inp_par_base                 = WIDGET_BASE(main_left_window,/row, frame=0)
rad_source_base              = WIDGET_BASE(inp_par_base,/row, frame=1, title='Radiation Sources')
rad_source_tab_base          = WIDGET_TAB(rad_source_base, /align_top, xsize=300, ysize=245, $
                                         uvalue='rad_source_tab_base')

;----------------------------------------------------------------------------------------------
; Input parameters for the stellar radiation source
;----------------------------------------------------------------------------------------------
star_source_base             = WIDGET_BASE(rad_source_tab_base,/column, frame=1, sensitive=1, title='Central star')

star_source_select_base      = WIDGET_BASE(star_source_base, /row)
star_source_title            = WIDGET_LABEL(star_source_select_base, value='Central stellar source : ')
star_source_select           = CW_BGROUP(star_source_select_base, ['Enable', 'Disable'],$
                                         /row, /return_name, uvalue='star_source_select', /exclusive)
        
star_source_inp_base         = WIDGET_BASE(star_source_base, /column, frame=0, sensitive=0)

star_source_temp_inp         = CW_FIELD(star_source_inp_base, title='Temperature [K] : ', $
                                        uvalue='star_source_temp_inp', xsize=10)
star_source_rad_inp          = CW_FIELD(star_source_inp_base, title='Radius :          ', $
                                        uvalue='star_source_rad_inp', xsize=10)
star_source_mass_inp         = CW_FIELD(star_source_inp_base, title='Mass :            ', $
                                        uvalue='star_source_mass_inp', xsize=10)

star_source_spec_base        = WIDGET_BASE(star_source_inp_base, /row)
star_source_spec_label       = WIDGET_LABEL(star_source_spec_base, value='Stellar spectrum : ', $
                                           xsize=110, /align_left)

star_source_spec_select      = CW_BGROUP(star_source_spec_base, ['Blackbody', 'Kurucz'], /return_name,$
                                         uvalue='star_source_spec_select', /no_release, /exclusive, /row)
;----------------------------------------------------------------------------------------------
; Input parameters for the continuous starlike radiation source
;----------------------------------------------------------------------------------------------
cstar_source_base            = WIDGET_BASE(rad_source_tab_base,/column, frame=1, sensitive=1, $
                                           title='Cont. starlike')

cstar_source_select_base     = WIDGET_BASE(cstar_source_base, /row)
cstar_source_title           = WIDGET_LABEL(cstar_source_select_base, value='Cont. starlike source : ')
cstar_source_select          = CW_BGROUP(cstar_source_select_base, ['Enable', 'Disable'],$
                                         /row, /return_name, uvalue='cstar_source_select', /exclusive)
  
cstar_source_inp_base        = WIDGET_BASE(cstar_source_base, /column, frame=0, sensitive=0)

cstar_source_rin_inp         = CW_FIELD(cstar_source_inp_base, title='Inner Radius :           ', $
                                        uvalue='cstar_source_rin_inp', xsize=10)
cstar_source_rout_inp        = CW_FIELD(cstar_source_inp_base, title='Outer Radius :           ', $
                                        uvalue='cstar_source_rin_inp', xsize=10)
cstar_source_tin_inp         = CW_FIELD(cstar_source_inp_base, title='Temperature at Rin :     ', $
                                        uvalue='cstar_source_rin_inp', xsize=10)
cstar_source_powex_inp       = CW_FIELD(cstar_source_inp_base, title='Temp. dist. pow. exp. :  ', $
                                        uvalue='cstar_source_rin_inp', xsize=10)
;----------------------------------------------------------------------------------------------
; Input parameters for the external radiation source
;----------------------------------------------------------------------------------------------
ext_source_base              = WIDGET_BASE(rad_source_tab_base,/column, frame=1, sensitive=1, $
                                           title='External source')

ext_source_select_base       = WIDGET_BASE(ext_source_base, /row)
ext_source_title             = WIDGET_LABEL(ext_source_select_base, value='External source : ')
ext_source_select            = CW_BGROUP(ext_source_select_base, ['Enable', 'Disable'],$
                                         /row, /return_name, uvalue='ext_source_select', /exclusive)

ext_source_inp_base          = WIDGET_BASE(ext_source_base, /column, frame=0, sensitive=0)
ext_source_temp_inp          = CW_FIELD(ext_source_inp_base, title='Temperature [K] :      ', $
                                        uvalue='ext_source_temp_inp', xsize=10)
;----------------------------------------------------------------------------------------------
; Density structures
;----------------------------------------------------------------------------------------------
dens_struct_base             = WIDGET_BASE(inp_par_base,/row, frame=1, title='Density structures')
dens_struct_tab_base         = WIDGET_TAB(dens_struct_base, /align_top, xsize=400, ysize=245, $
                                         uvalue='dens_struct_tab_base')
;----------------------------------------------------------------------------------------------
; Input parameters for a circumstellar disk
;----------------------------------------------------------------------------------------------
ppdisk_base                  = WIDGET_BASE(dens_struct_tab_base, /column, frame=1, $
                                           title='Disk')

ppdisk_select_base           = WIDGET_BASE(ppdisk_base, /row)
ppdisk_title                 = WIDGET_LABEL(ppdisk_select_base, value='Circumstellar Disk : ')
ppdisk_select                = CW_BGROUP(ppdisk_select_base, ['Enable', 'Disable'],$
                                         /row, /return_name, uvalue='ppdisk_select', /exclusive)

ppdisk_inp_base              = WIDGET_BASE(ppdisk_base, /row, frame=0, sensitive=1)
ppdisk_inp_left_base         = WIDGET_BASE(ppdisk_inp_base, /column, frame=0, sensitive=1)
ppdisk_inp_right_base        = WIDGET_BASE(ppdisk_inp_base, /column, frame=0, sensitive=1)

ppdisk_rin_inp               = CW_FIELD(ppdisk_inp_left_base, title='Inner Radius :         ', $
                                        uvalue='ppdisk_rin_inp', xsize=8)
ppdisk_rout_inp              = CW_FIELD(ppdisk_inp_left_base, title='Outer Radius :         ', $
                                        uvalue='ppdisk_rout_inp', xsize=8)
ppdisk_sig0_inp              = CW_FIELD(ppdisk_inp_left_base, title= 'Surf. dens. at Rout :  ', $
                                        uvalue='ppdisk_sig0_inp', xsize=8)
ppdisk_plsig1_inp            = CW_FIELD(ppdisk_inp_left_base, title='Surf. dens. pow. exp. :', $
                                        uvalue='ppdisk_plsig1_inp', xsize=8)
ppdisk_hrdisk_inp            = CW_FIELD(ppdisk_inp_left_base, title='Hp/R at Rout :         ', $
                                        uvalue='ppdisk_hrdisk_inp', xsize=8)

ppdisk_plh_inp               = CW_FIELD(ppdisk_inp_right_base, title='Flare index : ', $
                                        uvalue='ppdisk_plh_inp', xsize=8)
ppdisk_dustopac_inp          = CW_FIELD(ppdisk_inp_right_base, title='Dust opacity file :       ', $
                                        uvalue='ppdisk_dustopac_inp', xsize=8, /column)
ppdisk_plotopac_button       = WIDGET_BUTTON(ppdisk_inp_right_base, value='Plot opacity', $
                                             uvalue='ppdisk_plotopac_button')
;----------------------------------------------------------------------------------------------
; Input parameters for a circumstellar envelope
;----------------------------------------------------------------------------------------------
envelope_base                = WIDGET_BASE(dens_struct_tab_base, /column, frame=1, $
                                           title='Envelope')

envelope_select_base         = WIDGET_BASE(envelope_base, /row)
envelope_title               = WIDGET_LABEL(envelope_select_base, value='Circumstellar Envelope : ')
envelope_select              = CW_BGROUP(envelope_select_base, ['Enable', 'Disable'],$
                                         /row, /return_name, uvalue='envelope_select', /exclusive)

envelope_inp_base            = WIDGET_BASE(envelope_base, /row, frame=0, sensitive=1)
envelope_inp_left_base       = WIDGET_BASE(envelope_inp_base, /column, frame=0, sensitive=1)
envelope_inp_right_base      = WIDGET_BASE(envelope_inp_base, /column, frame=0, sensitive=1)

envelope_rin_inp             = CW_FIELD(envelope_inp_left_base, title='Inner Radius :       ', $
                                        uvalue='envelope_rin_inp', xsize=6)
envelope_rout_inp            = CW_FIELD(envelope_inp_left_base, title='Outer Radius :       ', $
                                        uvalue='envelope_rout_inp', xsize=6)
envelope_rho0_inp            = CW_FIELD(envelope_inp_left_base, title='Vol. density at Rout:', $
                                        uvalue='envelope_rho0_inp', xsize=6)
envelope_plrho_inp           = CW_FIELD(envelope_inp_left_base, title='Density pow. exp. :  ', $
                                        uvalue='envelope_plrho_inp', xsize=6)
envelope_cavrad_inp          = CW_FIELD(envelope_inp_left_base, title='Cavity radius [deg] :', $
                                        uvalue='envelope_cavrad_inp', xsize=6)


envelope_cavrfact_inp        = CW_FIELD(envelope_inp_right_base, title='Cavity dens. reduc. :', $
                                        uvalue='envelope_cavrfact_inp', xsize=5)
envelope_dustopac_inp        = CW_FIELD(envelope_inp_right_base, title='Dust opacity file :       ', $
                                        uvalue='envelope_dustopac_inp', xsize=8, /column)
envelope_plotopac_button     = WIDGET_BUTTON(envelope_inp_right_base, value='Plot opacity', $
                                             uvalue='envelope_plotopac_button')

;----------------------------------------------------------------------------------------------
; Input parameters for the background
;----------------------------------------------------------------------------------------------
bg_base                      = WIDGET_BASE(dens_struct_tab_base, /column, frame=1, $
                                           title='Background')

bg_select_base               = WIDGET_BASE(bg_base, /row)
bg_title                     = WIDGET_LABEL(bg_select_base, value='Background : ')
bg_select                    = CW_BGROUP(bg_select_base, ['Enable', 'Disable'],$
                                         /row, /return_name, uvalue='bg_select', /exclusive)

bg_inp_base                  = WIDGET_BASE(bg_base, /row, frame=0, sensitive=1)
bg_inp_left_base             = WIDGET_BASE(bg_inp_base, /column, frame=0, sensitive=1)
bg_inp_right_base            = WIDGET_BASE(bg_inp_base, /column, frame=0, sensitive=1)

bg_rho_inp                   = CW_FIELD(bg_inp_left_base, title='Backg. dens. [g/cm^3] :  ', $
                                        uvalue='bg_rin_inp', xsize=6)

bg_dustopac_inp              = CW_FIELD(bg_inp_left_base, title='Dust opacity file :       ', $
                                        uvalue='bg_dustopac_inp', xsize=8, /column)
bg_plotopac_button           = WIDGET_BUTTON(bg_inp_left_base, value='Plot opacity', $
                                             uvalue='bg_plotopac_button')



;----------------------------------------------------------------------------------------------
; Build the extra tabs
;----------------------------------------------------------------------------------------------
dum = uniq(params.ext.group)
next = n_elements(dum)
if params.ext(0).group_id eq -1 then next = 0

;;---------------------------
if next ge 1 then begin
;;
;; Look for the variable 'extra_'+group_name+'_tab_name' 
;;
   gid = dum(0)
   ii = where(params.ext.group eq params.ext(gid).group)
   jj = where(strtrim(params.ext(ii).name,2) eq 'extra_'+params.ext(gid).group+'_tab_name')
   if jj(0) ge 0 then begin
      tab_title = strsplit(params.ext(ii(jj)).value, "'", /extract)
      tab_title = strjoin(tab_title, /single)
   endif else begin
;;
;; If no variable 'extra_'+group_name+'_tab_name' is present use the
;; group name for the tab name
;;
      tab_title = params.ext(gid).group
   endelse
   
   ext1_base                    = WIDGET_BASE(dens_struct_tab_base, /column, frame=1, $
                                              title=tab_title)
   
   ext1_select_base             = WIDGET_BASE(ext1_base, /row)
   ext1_title                   = WIDGET_LABEL(ext1_select_base, value=tab_title+' : ')
   ext1_select                  = CW_BGROUP(ext1_select_base, ['Enable', 'Disable'],$
                                            /row, /return_name, uvalue='ext1_select', $
                                            uname='extra_'+params.ext(gid).group+'_enable',$
                                            /exclusive, set_value=0)
   ext1_inp = create_ext_tab2(parent=ext1_base, params=params, group=params.ext(gid).group)

   ii = where(strcompress(params.ext.name, /remove_all) eq 'extra_'+params.ext(gid).group+'_enable')

   if ii(0) ge 0 then begin
      if fix(params.ext(ii(0)).value) eq 0 then begin
         widget_control, ext1_inp.par_row_base, sensitive=0
         widget_control, ext1_select, set_value=1
      endif
   endif else begin
      widget_control, ext1_inp.par_row_base, sensitive=0
      widget_control, ext1_select, set_value=1
   endelse
   
endif else begin
   ext1_select = -1L
   ext1_inp = -1L
endelse
;;---------------------------
if next ge 2 then begin
;;
;; Look for the variable 'extra_'+group_name+'_tab_name' 
;;
   gid = dum(1)
   ii = where(params.ext.group eq params.ext(gid).group)
   jj = where(strtrim(params.ext(ii).name,2) eq 'extra_'+params.ext(gid).group+'_tab_name')
   if jj(0) ge 0 then begin
      tab_title = strsplit(params.ext(ii(jj)).value, "'", /extract)
      tab_title = strjoin(tab_title, /single)
   endif else begin
;;
;; If no variable 'extra_'+group_name+'_tab_name' is present use the
;; group name for the tab name
;;
      tab_title = params.ext(gid).group
   endelse
   
   ext2_base                    = WIDGET_BASE(dens_struct_tab_base, /column, frame=1, $
                                              title=tab_title)
   
   ext2_select_base             = WIDGET_BASE(ext2_base, /row)
   ext2_title                   = WIDGET_LABEL(ext2_select_base, value=tab_title+' : ')
   ext2_select                  = CW_BGROUP(ext2_select_base, ['Enable', 'Disable'],$
                                            /row, /return_name, uvalue='ext2_select', $
                                            uname='extra_'+params.ext(gid).group+'_enable',$
                                            /exclusive, set_value=0)
   

   ext2_inp = create_ext_tab2(parent=ext2_base, params=params, group=params.ext(gid).group)

   ii = where(strcompress(params.ext.name, /remove_all) eq 'extra_'+params.ext(gid).group+'_enable')
   if ii(0) ge 0 then begin
      if fix(params.ext(ii(0)).value) eq 0 then begin
         widget_control, ext2_inp.par_row_base, sensitive=0
         widget_control, ext2_select, set_value=1
      endif
   endif else begin
      widget_control, ext2_inp.par_row_base, sensitive=0
      widget_control, ext2_select, set_value=1
   endelse
endif else begin
   ext2_select = -1L
   ext2_inp = -1L
endelse

;;---------------------------
if next ge 3 then begin
;;
;; Look for the variable 'extra_'+group_name+'_tab_name' 
;;
   gid = dum(2)
   ii = where(params.ext.group eq params.ext(gid).group)
   jj = where(strtrim(params.ext(ii).name,2) eq 'extra_'+params.ext(gid).group+'_tab_name')
   if jj(0) ge 0 then begin
      tab_title = strsplit(params.ext(ii(jj)).value, "'", /extract)
      tab_title = strjoin(tab_title, /single)
   endif else begin
;;
;; If no variable 'extra_'+group_name+'_tab_name' is present use the
;; group name for the tab name
;;
      tab_title = params.ext(gid).group
   endelse
   
   ext3_base                    = WIDGET_BASE(dens_struct_tab_base, /column, frame=1, $
                                              title=tab_title)
   
   ext3_select_base             = WIDGET_BASE(ext3_base, /row)
   ext3_title                   = WIDGET_LABEL(ext3_select_base, value=tab_title+' : ')
   ext3_select                  = CW_BGROUP(ext3_select_base, ['Enable', 'Disable'],$
                                            /row, /return_name, uvalue='ext3_select', $
                                            uname='extra_'+params.ext(gid).group+'_enable',$
                                            /exclusive, set_value=0)
   
   ext3_inp = create_ext_tab2(parent=ext3_base, params=params, group=params.ext(gid).group)

   ii = where(strcompress(params.ext.name, /remove_all) eq 'extra_'+params.ext(gid).group+'_enable')
   if ii(0) ge 0 then begin
      if fix(params.ext(ii(0)).value) eq 0 then begin
         widget_control, ext3_inp.par_row_base, sensitive=0
         widget_control, ext3_select, set_value=1
      endif
   endif else begin
      widget_control, ext3_inp.par_row_base, sensitive=0
      widget_control, ext3_select, set_value=1
   endelse
endif else begin
   ext3_select = -1L
   ext3_inp = -1L
endelse

;;---------------------------
if next ge 4 then begin
;;
;; Look for the variable 'extra_'+group_name+'_tab_name' 
;;
   gid = dum(3)
   ii = where(params.ext.group eq params.ext(gid).group)
   jj = where(strtrim(params.ext(ii).name,2) eq 'extra_'+params.ext(gid).group+'_tab_name')
   if jj(0) ge 0 then begin
      tab_title = strsplit(params.ext(ii(jj)).value, "'", /extract)
      tab_title = strjoin(tab_title, /single)
   endif else begin
;;
;; If no variable 'extra_'+group_name+'_tab_name' is present use the
;; group name for the tab name
;;
      tab_title = params.ext(gid).group
   endelse
   
   ext4_base                    = WIDGET_BASE(dens_struct_tab_base, /column, frame=1, $
                                              title=tab_title)
   
   ext4_select_base             = WIDGET_BASE(ext4_base, /row)
   ext4_title                   = WIDGET_LABEL(ext4_select_base, value=tab_title+' : ')
   ext4_select                  = CW_BGROUP(ext4_select_base, ['Enable', 'Disable'],$
                                            /row, /return_name, uvalue='ext4_select',  $
                                            uname='extra_'+params.ext(gid).group+'_enable',$
                                            /exclusive, set_value=0)
   
   ext4_inp = create_ext_tab2(parent=ext4_base, params=params, group=params.ext(gid).group)

   ii = where(strcompress(params.ext.name, /remove_all) eq 'extra_'+params.ext(gid).group+'_enable')
   if ii(0) ge 0 then begin
      if fix(params.ext(ii(0)).value) eq 0 then begin
         widget_control, ext4_inp.par_row_base, sensitive=0
         widget_control, ext4_select, set_value=1
      endif
   endif else begin
      widget_control, ext4_inp.par_row_base, sensitive=0
      widget_control, ext4_select, set_value=1
   endelse
endif else begin
   ext4_select = -1L
   ext4_inp = -1L
endelse

;----------------------------------------------------------------------------------------------
; Grid
;----------------------------------------------------------------------------------------------
grid_base                    = WIDGET_BASE(inp_par_base,/row, frame=1, title='Grid')
grid_tab_base                = WIDGET_TAB(grid_base, /align_top, ysize=245, $
                                         uvalue='grid_tab_base')

;----------------------------------------------------------------------------------------------
; Spatial grid parameters
;----------------------------------------------------------------------------------------------
inp_spatial_grid_base        = WIDGET_BASE(grid_tab_base, /column, frame=1, sensitive=1, title='Spatial grid')

spatial_grid_rin_inp         = CW_FIELD(inp_spatial_grid_base, title='R min :                ', $
                                        uvalue='spatial_grid_rin_inp', xsize=7)
spatial_grid_rout_inp        = CW_FIELD(inp_spatial_grid_base, title='R max :                ', $
                                        uvalue='spatial_grid_rout_inp', xsize=7)
spatial_grid_nr_inp          = CW_FIELD(inp_spatial_grid_base, title='Nr of R grid points :  ', $
                                        uvalue='spatial_grid_nr_inp', xsize=7)

spatial_grid_theta_max_inp   = CW_FIELD(inp_spatial_grid_base, title='Theta max [deg] :      ', $
                                        uvalue='spatial_grid_theta_max_inp', xsize=7)
spatial_grid_nr_theta_inp    = CW_FIELD(inp_spatial_grid_base,  title='Nr of theta gridp. :   ', $
                                        uvalue='spatial_grid_theta_max_inp', xsize=7)

;----------------------------------------------------------------------------------------------
; Wavelength grid parameters
;----------------------------------------------------------------------------------------------
inp_wav_grid_base            = WIDGET_BASE(grid_tab_base, /column, frame=1, sensitive=1, $
                                           title='Wavel. grid', uvalue='inp_wav_grid_base')

wav_grid_lam_bound_inp       = CW_FIELD(inp_wav_grid_base, title='Lambda bounds :              ', $
                                        uvalue='wav_grid_lam_bound_inp', xsize=13, /column)
wav_grid_nr_inp              = CW_FIELD(inp_wav_grid_base, title='Nr of wavelengts :           ', $
                                        uvalue='wav_grid_nr_inp', xsize=7, /column)
;----------------------------------------------------------------------------------------------
; Big Red Buttons
;----------------------------------------------------------------------------------------------
inp_run_base                 = WIDGET_BASE(inp_par_base, /column, frame=1, sensitive=1)
inp_run_title                = WIDGET_LABEL(inp_run_base, value='Run RADMC3D')

run_nr_phot_inp              = CW_FIELD(inp_run_base,  title='Nr of photons:', $
                                        uvalue='run_nr_phot_inp', xsize=8)
run_therm_mc_button          = WIDGET_BUTTON(inp_run_base,  value='Thermal MC', uvalue='run_therm_mc',$
                                             xsize=150, ysize=35, sensitive=1, $
                                             tooltip=' Run the thermal Monte-Carlo simulation'+$
                                             'to calculate the temperature structure')
run_scat_mc_button           = WIDGET_BUTTON(inp_run_base,  value='Scattering MC', uvalue='run_scat_mc',$
                                             xsize=150, ysize=40, sensitive=0, $
                                             tooltip=' Run a monochromatic Monte-Carlo simulation '+$
                                             'to calculate anizotropic scattering (currently disabled)')

run_incl_inp                 = CW_FIELD(inp_run_base, title='Incl. [deg] : ', uvalue='run_incl_inp', $
                                        xsize=8, value='45')

run_make_spectrum_button     = WIDGET_BUTTON(inp_run_base,  value='Make spectrum', uvalue='run_make_spectrum',$
                                             xsize=150, ysize=40, sensitive=1,$
                                             tooltip='Run RADMC3D to calculate a spectrum at '+$
                                             'the given inclination angle')
run_make_image_button        = WIDGET_BUTTON(inp_run_base,  value='Make image', uvalue='run_make_image',$
                                             xsize=150, ysize=40, sensitive=1,$
                                             tooltip='Run RADMC3D to calculate an image ... GUI '+$
                                             'is under developement')
;----------------------------------------------------------------------------------------------
; Realize the widget
;----------------------------------------------------------------------------------------------
WIDGET_CONTROL, main_window, /realize
WIDGET_CONTROL, spec_imag_drawer, get_value = spec_imag_drawer
WIDGET_CONTROL, disk_struct_drawer, get_value = disk_struct_drawer
device, decomposed=0

wset, spec_imag_drawer
xyouts, 0.5, 0.5, 'SED', /normal, align=0.5, charsize=3

wset, disk_struct_drawer
xyouts, 0.5, 0.55, 'Density / Temperature /', /normal, align=0.5, charsize=3
xyouts, 0.5, 0.45, 'Optical depth', /normal, align=0.5, charsize=3
;
; Make a big fat structure for all widgets I'll need to operate with
;
state = {main_window:main_window,$ ; Main window
         file_menu_new:file_menu_new, $ ; New model 
         file_menu_open:file_menu_open, $ ; Open an existing model
         file_menu_save:file_menu_save, $ ; Save the current model - Does it have a meaning??
         file_menu_quit:file_menu_quit, $ ; Quit the GUI
         pref_menu_paths:pref_menu_paths, $ ; Set some paths - Does it have a meaning?
         help_menu_general_help:help_menu_general_help, $ ; General help to get started
         help_menu_parameter_help:help_menu_parameter_help,$ ; Help on the parameters
         help_menu_natconst_help:help_menu_natconst_help,$ ; Help on natural constants 
         spec_imag_drawer:spec_imag_drawer,$ ; Drawing window for the SED
         spec_imag_xrange_inp:spec_imag_xrange_inp, $ ; X-range for the SED drawing window
         spec_imag_yrange_inp:spec_imag_yrange_inp,$ ; Y-range for the SED drawing window
         spec_imag_dist_inp:spec_imag_dist_inp,$ ; Distance of the source in pc [to scale the SED] 
         spec_imag_oplot_data_button:spec_imag_oplot_data_button, $   ; Overplot data to the SED
         spec_imag_plot_spec_button:spec_imag_plot_spec_button, $ ; Overplot the stellar spectrum
         spec_imag_make_ps_button:spec_imag_make_ps_button, $ ; Plot the SED window into a postscript file
         spec_imag_xunit_select:spec_imag_xunit_select, $ ; Unit of the X-axis in the SED drawing window
         spec_imag_yunit_select:spec_imag_yunit_select,$ ; Unit of the Y-axis in the SED drawing window
         disk_struct_drawer:disk_struct_drawer, $ ; Drawing window for the disk structure
         disk_struct_yslide:disk_struct_yslide, $ ; Slider for maximum of the Y-axis (for 'zoom'-ing purpose)
         disk_struct_xslide:disk_struct_xslide, $ ; Slider for maximum of the X-axis (for 'zoom'-ing purpose)
         plot2d_clines_button:plot2d_clines_button, $ ; Plot contour lines of physical quantities in 2D
         plot2d_cfill_button:plot2d_cfill_button,$ ; Plot filled contours of physical quantities in 2D
         plot1d_coord_button:plot1d_coord_button, $ ; Choose the type of the X-scale in 1D plot (R or theta)
         plot1d_var_button:plot1d_var_button, $ ; Plot physical quantities as a function of either R or theta
         plot_nd_tauwav_inp:plot_nd_tauwav_inp, $ ; Wavelength at which the optical depth is calculated
         plot_nd_axis_scale_button:plot_nd_axis_scale_button, $ ; Choose the axis scale type of the disk structure window (lin/log)
         plot_nd_coord_sys_button:plot_nd_coord_sys_button, $ ; Choose the coordinate system type (Descartes vs. polar)
         plot_struct_update_button:plot_struct_update_button, $ ; Update the disk structure window according to the actual parameters
         plot_struct_ps_button:plot_struct_ps_button,$ ; Plot the disk structure window into a postscript file
         rad_source_tab_base:rad_source_tab_base, $ ; Tab widget containing radiation sources
         star_source_select:star_source_select, $ ; Enable/Disable a central star as a radiation source
         star_source_temp_inp:star_source_temp_inp, $ ; Temperature of the central star
         star_source_rad_inp:star_source_rad_inp, $ ; Radius of the central star
         star_source_mass_inp:star_source_mass_inp, $ ; Mass of the central star
         star_source_spec_select:star_source_spec_select, $ ; Type of the stellar radiation field (blackbody / Kurucz-model)
         star_source_inp_base:star_source_inp_base, $ ; Base widget containing the central star parameters in rad_source_tab_base
         cstar_source_select:cstar_source_select, $ ; Enable/Disable continuous starlike radiation source
         cstar_source_rin_inp:cstar_source_rin_inp,$ ; Inner radius of the radiation source
         cstar_source_rout_inp:cstar_source_rout_inp, $ ; Outer radius of the radiation source
         cstar_source_tin_inp:cstar_source_tin_inp,$ ; Temperature at the inner radius
         cstar_source_powex_inp:cstar_source_powex_inp, $ ; Power exponent of the radial temperature distribution
         cstar_source_inp_base:cstar_source_inp_base, $ ; Base widget containing the continous starlike source parameters in rad_source_tab_base
         ext_source_select:ext_source_select, $ ; Enable/Disable external radiation field
         ext_source_inp_base:ext_source_inp_base, $ ; Base widget containing the external source parameters in rad_source_tab_base
         ext_source_temp_inp:ext_source_temp_inp, $ ; Temperature of the external radiation field (assumed to be a blackbody)
         dens_struct_tab_base:dens_struct_tab_base, $ ; Tab widget containing the density structure parameters
         ppdisk_select:ppdisk_select, $ ; Enable/Disable a circumstellar disk
         ppdisk_rin_inp:ppdisk_rin_inp, $ ; Inner radius of the disk
         ppdisk_rout_inp:ppdisk_rout_inp,$ ; Outer radius of the disk
         ppdisk_sig0_inp:ppdisk_sig0_inp, $ ; Surface density at the outer radius
         ppdisk_plsig1_inp:ppdisk_plsig1_inp, $ ; Power exponent of the radial surface density distribution
         ppdisk_hrdisk_inp:ppdisk_hrdisk_inp, $ ; Ratio of the pressure scale height to the radius at the outer radius
         ppdisk_plh_inp:ppdisk_plh_inp, $ ; Flare index
         ppdisk_dustopac_inp:ppdisk_dustopac_inp, $ ; Name of the dust opacity file used for the circumstellar disk
         ppdisk_plotopac_button:ppdisk_plotopac_button, $ ; Plot the dust opacity of the disk
         ppdisk_inp_base:ppdisk_inp_base, $ ; Base widget containing the disk parameters in dens_struct_tab_base
         envelope_select:envelope_select, $ ; Enable/Disable a spherical circumstellar envelope
         envelope_rin_inp:envelope_rin_inp, $ ; Inner radius of the envelope
         envelope_rout_inp:envelope_rout_inp,$ ; Outer radius of the envelope
         envelope_rho0_inp:envelope_rho0_inp, $ ; Mass of the envelope 
         envelope_plrho_inp:envelope_plrho_inp, $ ; Power exponent of the radial volume density distribution
         envelope_cavrad_inp:envelope_cavrad_inp, $ ; Radius of the cavity around the pole
         envelope_cavrfact_inp:envelope_cavrfact_inp, $ ; Reduction of the density in the polar cavity
         envelope_dustopac_inp:envelope_dustopac_inp,  $ ; Name of the dust opacity file used for the envelope
         envelope_plotopac_button:envelope_plotopac_button, $ ; Plot the opacity of the envelope
         envelope_inp_base:envelope_inp_base, $ ; Base widget containing the envelope parameters in dens_struct_tab_base
         bg_select:bg_select, $ ; Enable/Disable a constant background density into which the envelope/disk is enbedded
         bg_inp_base:bg_inp_base, $                               ; Base widget containing the background density parameters in dens_struct_tab_base
         bg_rho_inp:bg_rho_inp, $ ; Background density into which the envelope/disk is enbedded
         bg_dustopac_inp:bg_dustopac_inp, $ ; Name of the dust opacity file used for the background density
         bg_plotopac_button:bg_plotopac_button, $ ; Plot the opacity of the background
         grid_tab_base:grid_tab_base, $ ; Tab widget containing the spatial and the wavelength grid parameters
         spatial_grid_rin_inp:spatial_grid_rin_inp, $ ; Inner radius of the spatial grid
         spatial_grid_rout_inp:spatial_grid_rout_inp, $ ; Outer radius of the spatial grid
         spatial_grid_nr_inp:spatial_grid_nr_inp, $ ; Number of radial grid points
         spatial_grid_theta_max_inp:spatial_grid_theta_max_inp, $ ; Upper bound of the theta grid 
         spatial_grid_nr_theta_inp:spatial_grid_nr_theta_inp,$ ; Number of theta grid points
         wav_grid_lam_bound_inp:wav_grid_lam_bound_inp, $ ; Boundary values for the wavelength grid (the same style as RADMC2D)
         wav_grid_nr_inp:wav_grid_nr_inp, $ ; Number of grid points in each segments of the wavelength grid (the same style as RADMC2D)
         run_nr_phot_inp:run_nr_phot_inp, $ ; Number of photons used for the thermal Monte-Carlo simulation
         run_therm_mc_button:run_therm_mc_button, $ ; Run the thermal Monte-Carlo
         run_scat_mc_button:run_scat_mc_button, $ ; Run the single-wavelength Monte-Carlo for scattering (currently disabled)
         run_incl_inp:run_incl_inp, $ ; Inclination angel at which the spectrum/image should be calculated
         run_make_spectrum_button:run_make_spectrum_button, $ ; Run RADMC3D to make a spectrum
         run_make_image_button:run_make_image_button,$ ; Run RADMC3D to make an image
         ext1_select:ext1_select,$ ; 
         ext1_inp:ext1_inp, $      ; Extra tab #1 in the density setup field
         ext2_select:ext2_select,$ ; 
         ext2_inp:ext2_inp, $ ; Extra tab #2 in the density setup field
         ext3_select:ext3_select,$ ; 
         ext3_inp:ext3_inp, $ ; Extra tab #3 in the density setup field
         ext4_select:ext4_select,$ ; 
         ext4_inp:ext4_inp} ; Extra tab #4 in the density setup field
         
;----------------------------------------------------------------------------------------------
; Set some defaults in the GUI
;----------------------------------------------------------------------------------------------

;widget_control, state.ext1_inp, set_value={EXTRA_GAP_FUNC_INP:'ize'}

;; Fill up the widgets with the paremeter values read from the file
fillup_all_widget_values, state, params

;; Fill up the extra tabs with the paremeter values read from the file
dum = uniq(params.ext.group)
next = n_elements(dum)
if params.ext(0).group_id eq -1 then next = 0

if next gt 0 then fillup_ext_tabs, state.ext1_inp, params, group=params.ext(dum(0)).group
if next gt 1 then fillup_ext_tabs, state.ext2_inp, params, group=params.ext(dum(1)).group
if next gt 2 then fillup_ext_tabs, state.ext3_inp, params, group=params.ext(dum(2)).group
if next gt 3 then fillup_ext_tabs, state.ext4_inp, params, group=params.ext(dum(3)).group

;; Run the problem_setup.pro to set up the model
spawn, 'idl do_setup'
;; Check if RADMC3D has already be run and we have the temperature and
;; density distributions
dum = findfile('dust_temperature.dat', count=tcount)
dum = findfile('dust_density.inp', count=dcount)
dum = findfile('spectrum.out', count=scount)

;; Read the density and the temperature if these files are present
struct = read_data(/ddens)

;; Calculate the optical depth (TODO; currently I am only using
;; the disk's opacity file - later on each structure needs to
;;                           have its own opacity!)
widget_control, state.ppdisk_dustopac_inp, get_value=dustopac_fname
tau    = get_tau_1dust(dens=struct, wav=0.55, dustopac_fname=strcompress(dustopac_fname, /remove_all),$
                      /imirt)
;
; Set the disk structure x-, y-range sliders maximum value to the 
;   number of r, theta grid points
;
widget_control, state.disk_struct_xslide, set_slider_max=struct.grid.nr-1
widget_control, state.disk_struct_yslide, set_slider_max=struct.grid.ntheta-1
au = 1.496d13
;; Set the density contour lines to be plotted in the disk structure window
cont_line_var = [0, 1, 0, 0]
;; Set the density contours (filled) to be plotted in the disk
;; structure window
cont_fill_var = 1
;; Set the x-axis to be logarithmic and the y-axis to be linear
xlog=1
ylog=0
;; Set the temperature to be selected in the 1D plot (not plotted!)
plot_var = 0
;; Set the radius to be on the X-axis in 1D plots (not plotted)
plot_crd = 0
;; Set the x-slider to its maximum value
plot_crd_id = struct.grid.ntheta-1
;; Choose 2D contours to be plotted instead of 1D curves
plot_dim = 2
;; Set coordinates the 2D contour plots to be R vs theta (~z/R) 
plot_crd_sys = 0
;; Set the x- and y-ranges in 1D plots
xrange_1d = [min(struct.grid.r)/au, max(struct.grid.r)/au]
yrange_1d = [0., 1.]
;; Set the x- and y-ranges in 2D plots
xrange_2d = [min(struct.grid.r)/au, max(struct.grid.r)/au]
yrange_2d = [0., max(!pi/2.-struct.grid.theta)]
   
;; Now update the screen (if we have something to plot, of course...)
wset, state.disk_struct_drawer
replot_2d_disk_struct

;; Set the x-axis to be wavelength (micron) in the SED plot
xunit = 0
;; Set the y-axis to be nu*Fnu (erg/s/cm/cm) in the SED plot
yunit = 1
;; Set the distance by default to 140pc
dist = 140.


;; Switch to SED plot instead of Kappa plot in the left screen
plot_sed_kappa = 0 
if scount gt 0 then begin
;; Read the spectrum
   spec = readspectrum()

   xrange_sed = [min(spec.lambda), max(spec.lambda)]
   dist = 140.
   yrange_sed = [min(spec.spectrum/dist^2.*spec.freq)*2., max(spec.spectrum/dist^2.*spec.freq)]
;; Update the x- and y-range input widgets
   widget_control, state.spec_imag_xrange_inp, $
                   set_value = [strcompress(xrange_sed(0), /remove_all)+', '+strcompress(xrange_sed(1), /remove_all)]
   widget_control, state.spec_imag_yrange_inp, $
                   set_value = [strcompress(yrange_sed(0), /remove_all)+', '+strcompress(yrange_sed(1), /remove_all)]

;; Now update the screen (if we have something to plot, of course...)
   wset, state.spec_imag_drawer
   replot_sed
endif else begin
   xrange_sed = [0.1, 1000.]
   yrange_sed = [1d-18, 1d-8]
endelse

;----------------------------------------------------------------------------------------------
; Release the control
;----------------------------------------------------------------------------------------------
WIDGET_CONTROL, main_window, set_uvalue=state
xmanager, 'radmc3d_gui', main_window, /no_block
end

;**********************************************************************************************
; Event handler for RADMC3D GUI
;     Attila Juhasz -  Heidelberg - 05.01.2010
;**********************************************************************************************
pro radmc3d_gui_event, xrevent

  COMMON PARAMETERS, params
  COMMON SPECTRUM, spec, xunit, yunit, xrange_sed, yrange_sed, dist, plot_sed_kappa, opac
  COMMON DISK_STRUCTURE, struct, tau, cont_line_var, cont_fill_var, xrange_2d, yrange_2d, xlog, ylog,$
     plot_var, plot_crd, plot_crd_id, xrange_1d, yrange_1d, plot_dim, plot_crd_sys

  WIDGET_CONTROL, xrevent.top, get_uvalue=state     
  WIDGET_CONTROL, xrevent.id, get_uvalue=uvalue   

  AU  = 1.49598d13              ; Astronomical Unit       [cm]


 
  case uvalue of
     'file_menu_new' : begin
        dum = dialog_message('This feature is under developement (i.e. temporarily disabled)!', /error, $
                             title='ERROR')
     end

     'file_menu_open' : begin
        dum = dialog_message('This feature is under developement (i.e. temporarily disabled)!', /error, $
                             title='ERROR')
     end
     
     'file_menu_save' : begin
        dum = dialog_message('This feature is under developement (i.e. temporarily disabled)!', /error, $
                             title='ERROR')
     end

     'file_menu_quit' : begin
        widget_control, xrevent.top, /destroy
     end

     'pref_menu_paths' : begin
        dum = dialog_message('This feature is under developement (i.e. temporarily disabled)!', /error, $
                             title='ERROR')
     end

     'help_menu_general_help' : begin
        display_help, help_id=0        
     end

     'help_menu_parameter_help' : begin
        display_help, help_id=1
     end

     'help_menu_natconst_help' : begin
        display_help, help_id=2
     end
;; ------------------------------------
;; Set the x-range of the SED plot
;; ------------------------------------
     'spec_imag_xrange_inp' : begin
        ;; Update the xrange
        widget_control, state.spec_imag_xrange_inp, get_value=txt
        txt = strsplit(txt, ',', /extract)
        xrange_sed = double(txt)
        ;; Update the yrange
        widget_control, state.spec_imag_yrange_inp, get_value=txt
        txt = strsplit(txt, ',', /extract)
        yrange_sed = double(txt)
        ;; Now update the screen
        wset, state.spec_imag_drawer
        replot_sed
     end
;; ------------------------------------
;; Set the y-range of the SED plot
;; ------------------------------------
     'spec_imag_yrange_inp' : begin
        ;; Update the xrange
        widget_control, state.spec_imag_xrange_inp, get_value=txt
        txt = strsplit(txt, ',', /extract)
        xrange_sed = double(txt)
        ;; Update the xrange
        widget_control, state.spec_imag_yrange_inp, get_value=txt
        txt = strsplit(txt, ',', /extract)
        yrange_sed = double(txt)
        ;; Now update the screen
        wset, state.spec_imag_drawer
        replot_sed
     end
;; ------------------------------------
;; Set the distance of the source 
;; ------------------------------------
     'spec_imag_dist_inp' : begin
        widget_control, state.spec_imag_dist_inp, get_value=txt
        dist = double(txt(0))
        
        ;; Update the yrange
        print, yrange_sed
        case yunit of
           0 : yrange_sed = [min(spec.spectrum), max(spec.spectrum)] / dist^2.
           1 : yrange_sed = [min(spec.spectrum*spec.freq), max(spec.spectrum*spec.freq)] / dist^2. 
           2 : yrange_sed = [min(spec.spectrum), max(spec.spectrum)] / dist^2. * 1d23
        endcase
        
        ;; Now update the screen
        wset, state.spec_imag_drawer
        replot_sed
     end
;; ------------------------------------
;; Overplot observational  data (currently disabled)
;; ------------------------------------
     'spec_imag_oplot_data_button' : begin
        dum = dialog_message('This feature is under developement (i.e. temporarily disabled)!', /error, $
                             title='ERROR')
     end
;; ------------------------------------
;; Plot the SED
;; ------------------------------------
     'spec_imag_plot_spec_button' : begin
        ;; Switch to 'plot_SED' mode for the postscript plotting
        plot_sed_kappa = 0
        ;; Check if we have a spectrum.out file at all
        dum = findfile('spectrum.out', count=scount)	
        if scount eq 0 then begin
           dum = dialog_message('No spectrum.out file has been found! Either RADMC3D has not yet been run or it did not finish properly.', /error, $
                             title='ERROR')
           goto, spec_imag_plot_skip
        endif else begin
           ;; Read the spectrum
           spec = readspectrum()
        endelse
        ;; Re-set the xrange
        case xunit of
           0 : xrange_sed = [min(spec.lambda), max(spec.lambda)]
           1 : xrange_sed = [min(spec.freq), max(spec.freq)]
        endcase
        ;; Re-set the yrange
        case yunit of
           0 : yrange_sed = [min(spec.spectrum), max(spec.spectrum)] / dist^2.
           1 : yrange_sed = [min(spec.spectrum*spec.freq), max(spec.spectrum*spec.freq)] / dist^2. 
           2 : yrange_sed = [min(spec.spectrum), max(spec.spectrum)] / dist^2. * 1d23
        endcase
        ;; Update the input text widgets
        txt = strcompress(xrange_sed(0), /remove_all) + ', ' + strcompress(xrange_sed(1))
        widget_control, state.spec_imag_xrange_inp, set_value=txt
        txt = strcompress(yrange_sed(0), /remove_all) + ', ' + strcompress(yrange_sed(1))
        widget_control, state.spec_imag_yrange_inp, set_value=txt
        

        ;; Now update the screen        
        wset, state.spec_imag_drawer
        replot_sed

spec_imag_plot_skip:
     end
;; ------------------------------------
;; Plot the SED into a postscript file
;; ------------------------------------
     'spec_imag_make_ps_button' : begin
        
        if plot_sed_kappa eq 0 then begin
           ;; Check if we have a spectrum.out file at all
           dum = findfile('spectrum.out', count=scount)	
           if scount eq 0 then begin
              dum = dialog_message('No spectrum.out file has been found! Either RADMC3D has not yet been run or it did not finish properly.', /error, $
                                   title='ERROR')
              goto, spec_imag_ps_skip
           endif else begin
              ;; Read the spectrum
              spec = readspectrum()
           endelse
           ;; Pick the file name in which we want to plot the dust opacity
           ps_filename = DIALOG_PICKFILE(file='sed.eps', /write,$
                                      /overwrite_prompt)
           if strlen(ps_filename) ne 0 then begin
              if strmid(ps_filename, 0, 1, /reverse_offset) ne '/' then begin
                 

                 set_plot, 'ps'
                 device, file=ps_filename, xsize=28, ysize=20, /landscape, /color, bits=8, xoffset=0, yoffset=28, $
                         /encapsulated
                 
                 replot_sed, charsize=1.6, xthick=2, ythick=2, charthick=2
                 
                 device, /close
                 set_plot, 'x'
                 spawn, 'gv '+ps_filename+' &'
              endif
           endif
             
        endif else begin
           ;; Pick the file name in which we want to plot the dust opacity
           ps_filename = DIALOG_PICKFILE(file='dustkappa.eps', /write,$
                                      /overwrite_prompt)
           if strlen(ps_filename) ne 0 then begin
              if strmid(ps_filename, 0, 1, /reverse_offset) ne '/' then begin
                 

                 set_plot, 'ps'
                 device, file=ps_filename, xsize=28, ysize=20, /landscape, /color, bits=8, xoffset=0, yoffset=28, $
                         /encapsulated
                 
                 plot_kappa, charsize=1.6, xthick=2, ythick=2, charthick=2
                 
                 device, /close
                 set_plot, 'x'
                 spawn, 'gv '+ps_filename+' &'
              endif
           endif
        endelse

spec_imag_ps_skip:
     end
;; ------------------------------------
;; Select the unit of the X-axis
;; ------------------------------------
     'spec_imag_xunit_select' : begin
        plot_sed_kappa = 0
        xunit = xrevent.index
        case xunit of
           0 : xrange_sed = [min(spec.lambda), max(spec.lambda)]
           1 : xrange_sed = [min(spec.freq), max(spec.freq)]
        endcase
        ;; Update the x-range input widgets
        widget_control, state.spec_imag_xrange_inp, $
                        set_value = [strcompress(xrange_sed(0), /remove_all)+$
                                     ', '+strcompress(xrange_sed(1), /remove_all)]
        ;; Now update the screen
        wset, state.spec_imag_drawer
        replot_sed
     end
;; ------------------------------------
;; Select the unit of the Y-axis
;; ------------------------------------
     'spec_imag_yunit_select' : begin
        plot_sed_kappa = 0
        yunit = xrevent.index
        case yunit of
           0 : yrange_sed = [min(spec.spectrum), max(spec.spectrum)] / dist^2.
           1 : yrange_sed = [min(spec.spectrum*spec.freq), max(spec.spectrum*spec.freq)] / dist^2. 
           2 : yrange_sed = [min(spec.spectrum), max(spec.spectrum)] / dist^2. * 1d23
        endcase
        ;; Update the x-range input widgets
        widget_control, state.spec_imag_yrange_inp, $
                        set_value = [strcompress(yrange_sed(0), /remove_all)+$
                                     ', '+strcompress(yrange_sed(1), /remove_all)]
        ;; Now update the screen
        wset, state.spec_imag_drawer
        replot_sed
     end

;; ------------------------------------
;; Set the x-range of the plot using the slider
;; ------------------------------------
     'disk_struct_xslide' : begin
        grid_id = xrevent.value
        if plot_dim eq 1 then begin
           if plot_crd eq 0 then begin
              if grid_id gt struct.grid.nr-1 then grid_id = struct.grid.nr-1
              xrange_1d(1) = struct.grid.r(grid_id) / au
           endif
           if plot_crd eq 1 then begin
              if grid_id gt struct.grid.ntheta-1 then grid_id = struct.grid.ntheta-1
              xrange_1d(1) = !pi/2. - struct.grid.theta(struct.grid.ntheta-grid_id-1)
           endif
           wset, state.disk_struct_drawer
           replot_1d_disk_struct
        endif
        if plot_dim eq 2 then begin
           if grid_id gt struct.grid.nr-1 then grid_id = struct.grid.nr-1
           xrange_2d(1) = struct.grid.r(grid_id) / au
           wset, state.disk_struct_drawer
           replot_2d_disk_struct
        endif 
     end

;; ------------------------------------
;; Set the y-range of the plot using the slider
;; ------------------------------------
     'disk_struct_yslide' : begin
        grid_id = xrevent.value
        ;; If we have a 1D plot
        if plot_dim eq 1 then begin
           ;; Get the number of the grid point along which we plot the variables
           plot_crd_id = grid_id
           ;; Now update the screen
           wset, state.disk_struct_drawer
           replot_1d_disk_struct
        endif
        if plot_dim eq 2 then begin
           if plot_crd_sys eq 0 then yrange_2d(1) = !pi/2.-struct.grid.theta(struct.grid.ntheta-grid_id-1)
           if plot_crd_sys eq 1 then yrange_2d(1) = struct.grid.r(grid_id) / au
           wset, state.disk_struct_drawer
           replot_2d_disk_struct
        endif
     end
;; ------------------------------------
;; Plot contour lines in 2D 
;; ------------------------------------
     'plot2d_clines_button' : begin
        ;; If there is a 1D plot in the drawing widget and we switch
        ;; to a 2D plot we need to re-set the maximum value of the
        ;; sliders and the x- and y-range
        if plot_dim ne 2 then begin
           plot_dim = 2 
           widget_control, state.disk_struct_xslide, set_slider_max=struct.grid.nr-1
           widget_control, state.disk_struct_yslide, set_slider_max=struct.grid.ntheta-1
           xrange_2d = [min(struct.grid.r), max(struct.grid.r)]/au
           yrange_2d = [min(!pi/2-struct.grid.theta), max(!pi/2.-struct.grid.theta)]
        endif
        ;; Temperature contours
        if xrevent.value eq 0 then begin
           if xrevent.select eq 0 then cont_line_var[0] = 0 ; Turn off the contour lines
           if xrevent.select eq 1 then cont_line_var[0] = 1 ; Turn on the contour lines
        end
        ;; Density contours
        if xrevent.value eq 1 then begin
           if xrevent.select eq 0 then cont_line_var[1] = 0 ; Turn off the contour lines
           if xrevent.select eq 1 then cont_line_var[1] = 1 ; Turn on the contour lines
        end
        ;; Radial optical depth contours
        if xrevent.value eq 2 then begin
           if xrevent.select eq 0 then cont_line_var[2] = 0 ; Turn off the contour lines
           if xrevent.select eq 1 then cont_line_var[2] = 1 ; Turn on the contour lines
        end
        ;; Vertical optical depth contours
        if xrevent.value eq 3 then begin
           if xrevent.select eq 0 then cont_line_var[3] = 0 ; Turn off the contour lines
           if xrevent.select eq 1 then cont_line_var[3] = 1 ; Turn on the contour lines
        end
        
        ;; Now update the screen
        wset, state.disk_struct_drawer
        replot_2d_disk_struct  
        
     end
;; ------------------------------------
;; Plot filled contours in 2D
;; ------------------------------------
     'plot2d_cfill_button' : begin
        ;; If there is a 1D plot in the drawing widget and we switch
        ;; to a 2D plot we need to re-set the maximum value of the
        ;; sliders and the x- and y-range
        if plot_dim ne 2 then begin
           plot_dim = 2 
           widget_control, state.disk_struct_xslide, set_slider_max=struct.grid.nr-1
           widget_control, state.disk_struct_yslide, set_slider_max=struct.grid.ntheta-1
           xrange_2d = [min(struct.grid.r), max(struct.grid.r)]/au
           yrange_2d = [min(!pi/2-struct.grid.theta), max(!pi/2.-struct.grid.theta)]
        endif
        plot_dim = 2
        cont_fill_var = xrevent.value
        ;; Now update the screen
        wset, state.disk_struct_drawer
        replot_2d_disk_struct
     end
;; ------------------------------------
;; Choose what to plot on the x-axis in 1D plots
;; ------------------------------------
     'plot1d_coord_button' : begin
        plot_dim = 1
        plot_crd = xrevent.value
        ;; If there is a 2D plot in the drawing widget and we switch
        ;; to a 1D plot or we change the coordinate type on the
        ;; x-axis we need to re-set the maximum value of the
        ;; sliders and the x- and y-range
        if plot_crd eq 0 then begin
           widget_control, state.disk_struct_xslide, set_slider_max=struct.grid.nr-1
           widget_control, state.disk_struct_yslide, set_slider_max=struct.grid.ntheta-1
           xrange_1d = [min(struct.grid.r), max(struct.grid.r)]/au
        endif
        if plot_crd eq 1 then begin
           widget_control, state.disk_struct_xslide, set_slider_max=struct.grid.ntheta-1
           widget_control, state.disk_struct_yslide, set_slider_max=struct.grid.nr-1
           xrange_1d = [min(!pi/2-struct.grid.theta), max(!pi/2.-struct.grid.theta)]
        endif

        ;; Now update the screen
        wset, state.disk_struct_drawer
        replot_1d_disk_struct
     end
;; ------------------------------------
;; Choose what to plot on the y-axis in 1D plots
;; ------------------------------------
     'plot1d_var_button' : begin
        plot_dim = 1
        ;; If there is a 2D plot in the drawing widget and we switch
        ;; to a 1D plot or we change the coordinate type on the
        ;; x-axis we need to re-set the maximum value of the
        ;; sliders and the x- and y-range
        if plot_crd eq 0 then begin
           widget_control, state.disk_struct_xslide, set_slider_max=struct.grid.nr-1
           widget_control, state.disk_struct_yslide, set_slider_max=struct.grid.ntheta-1
           ;xrange_1d = [min(struct.grid.r), max(struct.grid.r)]/au
        endif
        if plot_crd eq 1 then begin
           widget_control, state.disk_struct_xslide, set_slider_max=struct.grid.ntheta-1
           widget_control, state.disk_struct_yslide, set_slider_max=struct.grid.nr-1
           ;xrange_1d = [min(!pi/2-struct.grid.theta), max(!pi/2.-struct.grid.theta)]
        endif

        plot_var = xrevent.value
        widget_control, state.disk_struct_yslide, get_value=plot_crd_id
        
        ;; Now update the screen
        wset, state.disk_struct_drawer
        replot_1d_disk_struct
     end

     'plot_nd_tauwav_inp' : begin
     end
;; ------------------------------------
;; Choose the scale of the axes (lin/log)
;; ------------------------------------
     'plot_nd_axis_scale_button' : begin
        if xrevent.value eq 0 then begin
           xlog = xrevent.select
        endif
        if xrevent.value eq 1 then begin
           ylog = xrevent.select
        endif
        
        ;; Now update the screen
        wset, state.disk_struct_drawer
        if plot_dim eq 1 then replot_1d_disk_struct
        if plot_dim eq 2 then replot_2d_disk_struct
     end
;; ------------------------------------
;; Choose the coordinate system (Descartes vs. polar)
;; ------------------------------------ 
     'plot_nd_coord_sys_button' : begin
        plot_crd_sys = xrevent.value
        
        if plot_crd_sys eq 0 then begin
           ;; Update the y-slider
           yrange_2d = [0., max(!pi/2.-struct.grid.theta)]
           widget_control, state.disk_struct_yslide, set_slider_max=struct.grid.ntheta-1
        endif 

        if plot_crd_sys eq 1 then begin
           ;; Use the same slider behaviour in both axes
           yrange_2d = [0, max(struct.grid.r/au)]
           widget_control, state.disk_struct_yslide, set_slider_max=struct.grid.nr-1
        endif
        ;; Now update the screen
        if plot_dim eq 2 then begin
           wset, state.disk_struct_drawer
           replot_2d_disk_struct
        endif

     end
;; ------------------------------------
;; Update the disk structure screen with 
;;  the actual values read from the widgets
;; ------------------------------------
     'plot_struct_update_button' : begin
        ;; Read out all input widgets
        read_all_widget_values, state, params
        ;; Read out all input widgets in the extra tabs

        dum = uniq(params.ext.group)
        next = n_elements(dum)
        if params.ext(0).group_id eq -1 then next = 0
        if next gt 0 then read_ext_tabs, state.ext1_inp, params, group=params.ext(dum(0)).group
        if next gt 1 then read_ext_tabs, state.ext2_inp, params, group=params.ext(dum(1)).group
        if next gt 2 then read_ext_tabs, state.ext3_inp, params, group=params.ext(dum(2)).group
        if next gt 3 then read_ext_tabs, state.ext4_inp, params, group=params.ext(dum(3)).group
        
        ;; Write the problem_params.pro file
        write_params_radmc3d, params, write_file_name='problem_params.pro'
        ;; Make the model density structure
        spawn, 'idl do_setup'
        
        ;; Run the extra function #1 if it is enabled 
        if state.ext1_select ne -1L then begin
           widget_control, state.ext1_select, get_value=dum
           ;; Check which group/function is associated with this tab
           if dum eq 0 then begin
              dum2 = widget_info(state.ext1_select, /uname)
              group = strsplit(dum2, '_', /extract)
              group = group(1)
              
              ;; Now get the function_name
              ii = where(strcompress(params.ext.name, /remove_all) eq 'extra_'+group+'_func')
              
              if ii(0) ge 0 then begin
                 if strmid(strcompress(params.ext(ii(0)).value, /remove_all), 0, 1) eq "'" then begin
                    pro_name = strsplit(params.ext(ii).value, "'", /extract)
                 endif else pro_name = strcompress(params.ext(ii(0)).value)
                 resolve_routine, pro_name, /compile_full_file
                 call_procedure, pro_name(0)
              endif
           endif
        endif
        
        if state.ext2_select ne -1L then begin
           ;; Run the extra function #2 if it is enabled 
           widget_control, state.ext2_select, get_value=dum
           ;; Check which group/function is associated with this tab
           if dum eq 0 then begin
              dum2 = widget_info(state.ext2_select, /uname)
              group = strsplit(dum2, '_', /extract)
              group = group(1)
              
              ;; Now get the function_name
              ii = where(strcompress(params.ext.name, /remove_all) eq 'extra_'+group+'_func')
              
              if ii(0) ge 0 then begin
                 if strmid(strcompress(params.ext(ii(0)).value, /remove_all), 0, 1) eq "'" then begin
                    pro_name = strsplit(params.ext(ii).value, "'", /extract)
                 endif else pro_name = strcompress(params.ext(ii(0)).value)
                 resolve_routine, pro_name, /compile_full_file
                 call_procedure, pro_name(0)
              endif
           endif
        endif
        
        if state.ext3_select ne -1L then begin
           ;; Run the extra function #3 if it is enabled 
           widget_control, state.ext3_select, get_value=dum
           ;; Check which group/function is associated with this tab
           if dum eq 0 then begin
              dum2 = widget_info(state.ext3_select, /uname)
              group = strsplit(dum2, '_', /extract)
              group = group(1)
              
              ;; Now get the function_name
              ii = where(strcompress(params.ext.name, /remove_all) eq 'extra_'+group+'_func')
              
              if ii(0) ge 0 then begin
                 if strmid(strcompress(params.ext(ii(0)).value, /remove_all), 0, 1) eq "'" then begin
                    pro_name = strsplit(params.ext(ii).value, "'", /extract)
                 endif else pro_name = strcompress(params.ext(ii(0)).value)
                 resolve_routine, pro_name, /compile_full_file
                 call_procedure, pro_name(0)
              endif
           endif
        endif

        if state.ext4_select ne -1L then begin
           ;; Run the extra function #4 if it is enabled 
           widget_control, state.ext4_select, get_value=dum
           ;; Check which group/function is associated with this tab
           if dum eq 0 then begin
              dum2 = widget_info(state.ext4_select, /uname)
              group = strsplit(dum2, '_', /extract)
              group = group(1)
              
              ;; Now get the function_name
              ii = where(strcompress(params.ext.name, /remove_all) eq 'extra_'+group+'_func')
              
              if ii(0) ge 0 then begin
                 if strmid(strcompress(params.ext(ii(0)).value, /remove_all), 0, 1) eq "'" then begin
                    pro_name = strsplit(params.ext(ii).value, "'", /extract)
                 endif else pro_name = strcompress(params.ext(ii(0)).value)
                 resolve_routine, pro_name, /compile_full_file
                 call_procedure, pro_name(0)
              endif
           endif
        endif
        
        ;; Read the model density and temperature structure
        ;; WARNING!!!! The temperature structure will not be appropriate
        ;; unless thermal Monte-Carlo is also run
        dum = findfile('dust_density.inp', count=dcount)
        dum = findfile('dust_temperature.dat', count=tcount)
        if tcount eq 0 then begin
           if dcount eq 0 then begin
              ;!!!!!!!! ERRROR!!! something went wrong!!!!
              dum = dialog_message('Neither the dust_density.inp nor the dust_temperature.dat was found. '+$
                                   'Somehow the dust density file was not created, although it should have been!', $
                                   /error, title='ERROR')
              goto, skip_update
           endif else begin
              struct = read_data(/ddens)
           endelse
        endif else begin
           if dcount eq 0 then begin
             ;!!!!!!!! ERRROR!!! something went wrong!!!!
              dum = dialog_message('Something weird is going on here... The dust_temperature.dat was found but the '+$
                                   'dust_density.inp file was not! We have the dust temperatures but which density distribution it belongs to is unknown...', $
                                   /error, title='ERROR')
              goto, skip_update
           endif else begin
              struct = read_data(/ddens, /dtemp)
           endelse
        endelse
        ;; Read the dust opacity file name
        widget_control, state.ppdisk_dustopac_inp, get_value=dustopac_fname
        ;; Read the wavelength at which the optical depth should be calculated
        widget_control, state.plot_nd_tauwav_inp, get_value=tauwav
        ;; Calculate the optical depth
        tau    = get_tau_1dust(dens=struct, wav=float(tauwav), dustopac_fname=strcompress(dustopac_fname, /remove_all),$
                               /imirt)
        ;; Get the coordinate system type
        ;; widget_control, state.plot_nd_coord_sys_button, get_value=plot_crd_sys
        ;;stop


        ;; Set the disk structure x-, y-range sliders maximum value to the 
        ;;   number of r, theta grid points
        widget_control, state.disk_struct_xslide, set_slider_max=struct.grid.nr-1
       
        if plot_dim eq 1 then begin
           if plot_crd eq 0 then begin
              ;; Update the y-slider
              widget_control, state.disk_struct_xslide, set_slider_max=struct.grid.nr-1
           endif 
           
           if plot_crd eq 1 then begin
              ;; Use the same slider behaviour in both axes
              widget_control, state.disk_struct_xslide, set_slider_max=struct.grid.ntheta-1
           endif
        endif
        
        if plot_dim eq 2 then begin
           if plot_crd_sys eq 0 then begin
              ;; Update the y-slider
              yrange_2d = [0., max(!pi/2.-struct.grid.theta)]
              widget_control, state.disk_struct_yslide, set_slider_max=struct.grid.ntheta-1
           endif 
           
           if plot_crd_sys eq 1 then begin
              ;; Use the same slider behaviour in both axes
              yrange_2d = xrange_2d
              widget_control, state.disk_struct_yslide, set_slider_max=struct.grid.nr-1
           endif
        endif
        
        ;; Now update the screen
        wset, state.disk_struct_drawer
        if plot_dim eq 1 then replot_1d_disk_struct
        if plot_dim eq 2 then replot_2d_disk_struct

skip_update:

     end
;; ------------------------------------
;; Create a postscript file from the disk structure screen
;; ------------------------------------

     'plot_struct_ps_button' : begin
        
        ;; Pick the file name in which we want to plot the disk structure
        ps_filename = DIALOG_PICKFILE(file='disk_struct.eps', /write,$
	 				/overwrite_prompt)
        if strlen(ps_filename) ne 0 then begin
           if strmid(ps_filename, 0, 1, /reverse_offset) ne '/' then begin
              
              set_plot, 'ps'
              device, file=ps_filename, xsize=28, ysize=20, /landscape, /color, bits=8, xoffset=0, yoffset=28, $
                      /encapsulated
              
              if plot_crd_sys eq 0 then begin
                 ;; Update the y-slider
                 yrange_2d = [0., max(!pi/2.-struct.grid.theta)]
                 widget_control, state.disk_struct_yslide, set_slider_max=struct.grid.ntheta-1
              endif 
              
              if plot_crd_sys eq 1 then begin
                 ;; Use the same slider behaviour in both axes
                 yrange_2d = xrange_2d
                 widget_control, state.disk_struct_yslide, set_slider_max=struct.grid.nr-1
              endif
              
              if plot_dim eq 1 then replot_1d_disk_struct, charsize=1.6, xthick=2, ythick=2, charthick=2
              if plot_dim eq 2 then replot_2d_disk_struct, charsize=1.6, xthick=2, ythick=2, charthick=2
              
              device, /close
              set_plot, 'x'
           endif
        endif
        spawn, 'gv '+ps_filename+' &'
     end

     'rad_source_tab_base' : begin
     end

;; ------------------------------------
;; Enable/Disable central star as a radiation source
;; ------------------------------------
     'star_source_select' : begin
        if xrevent.value eq 'Enable' then WIDGET_CONTROL, state.star_source_inp_base, sensitive=1
        if xrevent.value eq 'Disable' then WIDGET_CONTROL, state.star_source_inp_base, sensitive=0
     end
     
     'star_source_temp_inp' : begin
     end
     
     'star_source_rad_inp' : begin
     end

     'star_source_mass_inp' : begin
     end

     'star_source_spec_select' : begin
     end
;; ------------------------------------
;; Enable/Disable continuous starlike radiation source
;; ------------------------------------
     'cstar_source_select' : begin
        if xrevent.value eq 'Enable' then begin
           WIDGET_CONTROL, state.cstar_source_inp_base, sensitive=1
           dum = dialog_message('This feature is under developement (i.e. temporarily disabled)!', /error, $
                             title='ERROR')
        endif
        if xrevent.value eq 'Disable' then WIDGET_CONTROL, state.cstar_source_inp_base, sensitive=0
     end

     'cstar_source_rin_inp' : begin
     end

     'cstar_source_rout_inp' : begin
     end

     'cstar_source_tin_inp' : begin
     end

     'cstar_source_powex_inp' : begin
     end

;; ------------------------------------
;; Enable/Disable external radiation field
;; ------------------------------------
     'ext_source_select' : begin
        if xrevent.value eq 'Enable' then WIDGET_CONTROL, state.ext_source_inp_base, sensitive=1
        if xrevent.value eq 'Disable' then WIDGET_CONTROL, state.ext_source_inp_base, sensitive=0
     end
     
     'ext_source_temp_inp' : begin
     end

;; ------------------------------------
;; Enable/Disable circumstellar disk density distribution
;; ------------------------------------
      'ppdisk_select' : begin
         if xrevent.value eq 'Enable' then WIDGET_CONTROL, state.ppdisk_inp_base, sensitive=1
         if xrevent.value eq 'Disable' then WIDGET_CONTROL, state.ppdisk_inp_base, sensitive=0
      end
      
      'dens_struct_tab_base' : begin
      end
      
      'ppdisk_rin_inp' : begin
      end
      
      'ppdisk_rout_inp' : begin
      end
      
      'ppdisk_sig0_inp' : begin
      end
      
      'ppdisk_plsig1_inp' : begin
      end
      
      'ppdisk_hrdisk_inp' : begin
      end
      
      'ppdisk_plh_inp' : begin
      end
      
      'ppdisk_dustopac_inp' : begin
      end
;; ------------------------------------
;; Plot the dust opacity of the disk
;; ------------------------------------
      'ppdisk_plotopac_button' : begin
         ;; Switch to 'plot_kappa' mode for the postscript plotting
         plot_sed_kappa = 1
         widget_control, state.ppdisk_dustopac_inp, get_value=txt
         txt = strsplit(txt, "'", /extract)
         openr, 1, txt(0)
         readf, 1, iformat
         readf, 1, nlam
         opac = {wave:dblarr(nlam), cabs:dblarr(nlam), csca:dblarr(nlam)}

         if iformat eq 1 then begin
            dum = [0d0, 0d0]
            for i=0, nlam-1 do begin
               readf, 1, dum
               opac.wave(i) = dum(0)
               opac.cabs(i) = dum(1)
            endfor
         endif
         if iformat ge 2 then begin
            dum = [0d0, 0d0, 0d0]
            for i=0, nlam-1 do begin
               readf, 1, dum
               opac.wave(i) = dum(0)
               opac.csca(i) = dum(1)
               opac.csca(i) = dum(2)
            endfor
         endif
         close, 1
         
         ;; Update the x- and y-range
         xrange_sed = [min(opac.wave), max(opac.wave)]
         yrange_sed = [min([opac.cabs]), max([opac.cabs])]

         ;; Now update the screen
         wset, state.spec_imag_drawer
         plot_kappa
      end
;; ------------------------------------
;; Enable/Disable circumstellar envelope density distribution
;; ------------------------------------
     'envelope_select' : begin
        if xrevent.value eq 'Enable' then WIDGET_CONTROL, state.envelope_inp_base, sensitive=1
        if xrevent.value eq 'Disable' then WIDGET_CONTROL, state.envelope_inp_base, sensitive=0
     end

     'envelope_rin_inp' : begin
     end

     'envelope_rout_inp' : begin
     end

     'envelope_rho0_inp' : begin
     end
     
     'envelope_plrho_inp' : begin
     end

     'envelope_cavrad_inp' : begin
     end

     'envelope_cavrfact_inp' : begin
     end

     'envelope_dustopac_inp' : begin
     end
;; ------------------------------------
;; Plot the dust opacity of the envelope
;; ------------------------------------
     'envelope_plotopac_button' : begin
        ;; Switch to 'plot_kappa' mode for the postscript plotting
        plot_sed_kappa = 1
        widget_control, state.envelope_dustopac_inp, get_value=txt
        txt = strsplit(txt, "'", /extract)
        openr, 1, txt(0)
        readf, 1, iformat
        readf, 1, nlam
        opac = {wave:dblarr(nlam), cabs:dblarr(nlam), csca:dblarr(nlam)}
        
        if iformat eq 1 then begin
           dum = [0d0, 0d0]
           for i=0, nlam-1 do begin
              readf, 1, dum
              opac.wave(i) = dum(0)
               opac.cabs(i) = dum(1)
            endfor
         endif
         if iformat ge 2 then begin
            dum = [0d0, 0d0, 0d0]
            for i=0, nlam-1 do begin
               readf, 1, dum
               opac.wave(i) = dum(0)
               opac.csca(i) = dum(1)
               opac.csca(i) = dum(2)
            endfor
         endif
         close, 1
         
         ;; Update the x- and y-range
         xrange_sed = [min(opac.wave), max(opac.wave)]
         yrange_sed = [min([opac.cabs]), max([opac.cabs])]

         ;; Now update the screen
         wset, state.spec_imag_drawer
         plot_kappa
      end
;; ------------------------------------
;; Enable/Disable background density 
;;  (I think it is supposed to be suggested that
;;   this should always be turned on as in RADMC2D) 
;; ------------------------------------
     'bg_select' : begin
        if xrevent.value eq 'Enable' then WIDGET_CONTROL, state.bg_inp_base, sensitive=1
        if xrevent.value eq 'Disable' then WIDGET_CONTROL, state.bg_inp_base, sensitive=0
     end
     
     'bg_rho_inp' : begin
     end   

     'bg_dustopac_inp' : begin
     end   

;; ------------------------------------
;; Enable/Disable extra function #1
;; ------------------------------------

     'ext1_select' : begin
        if xrevent.value eq 'Enable' then begin
           widget_control, state.ext1_inp.par_row_base, sensitive=1
           name = widget_info(state.ext1_select, /uname)
           ii = where(strcompress(params.ext.name, /remove_all) eq name(0))
           if ii(0) ge 0 then params.ext(ii(0)).value = '1'
           id = widget_info(state.ext1_inp.par_row_base, find_by_uname=name)
           widget_control, id, set_value='1'
        endif
           
        if xrevent.value eq 'Disable' then begin
           widget_control, state.ext1_inp.par_row_base, sensitive=0
           name = widget_info(state.ext1_select, /uname)
           ii = where(strcompress(params.ext.name, /remove_all) eq name(0))
           if ii(0) ge 0 then params.ext(ii(0)).value = '0'
           id = widget_info(state.ext1_inp.par_row_base, find_by_uname=name)
           widget_control, id, set_value='0'
        endif
      end

;; ------------------------------------
;; Enable/Disable extra function #2
;; ------------------------------------

     'ext2_select' : begin
        if xrevent.value eq 'Enable' then begin
           widget_control, state.ext2_inp.par_row_base, sensitive=1
           name = widget_info(state.ext2_select, /uname)
           ii = where(strcompress(params.ext.name, /remove_all) eq name(0))
           if ii(0) ge 0 then params.ext(ii(0)).value = '1'
           id = widget_info(state.ext2_inp.par_row_base, find_by_uname=name)
           widget_control, id, set_value='1'
        endif
           
        if xrevent.value eq 'Disable' then begin
           widget_control, state.ext2_inp.par_row_base, sensitive=0
           name = widget_info(state.ext2_select, /uname)
           ii = where(strcompress(params.ext.name, /remove_all) eq name(0))
           if ii(0) ge 0 then params.ext(ii(0)).value = '0'
           id = widget_info(state.ext2_inp.par_row_base, find_by_uname=name)
           widget_control, id, set_value='0'
        endif
      end

;; ------------------------------------
;; Enable/Disable extra function #3
;; ------------------------------------

     'ext3_select' : begin
        if xrevent.value eq 'Enable' then begin
           widget_control, state.ext3_inp.par_row_base, sensitive=1
           name = widget_info(state.ext3_select, /uname)
           ii = where(strcompress(params.ext.name, /remove_all) eq name(0))
           if ii(0) ge 0 then params.ext(ii(0)).value = '1'
           id = widget_info(state.ext3_inp.par_row_base, find_by_uname=name)
           widget_control, id, set_value='1'
        endif
           
        if xrevent.value eq 'Disable' then begin
           widget_control, state.ext3_inp.par_row_base, sensitive=0
           name = widget_info(state.ext3_select, /uname)
           ii = where(strcompress(params.ext.name, /remove_all) eq name(0))
           if ii(0) ge 0 then params.ext(ii(0)).value = '0'
           id = widget_info(state.ext3_inp.par_row_base, find_by_uname=name)
           widget_control, id, set_value='0'
        endif
      end

;; ------------------------------------
;; Enable/Disable extra function #4
;; ------------------------------------

     'ext4_select' : begin
        if xrevent.value eq 'Enable' then begin
           widget_control, state.ext3_inp.par_row_base, sensitive=1
           name = widget_info(state.ext4_select, /uname)
           ii = where(strcompress(params.ext.name, /remove_all) eq name(0))
           if ii(0) ge 0 then params.ext(ii(0)).value = '1'
           id = widget_info(state.ext4_inp.par_row_base, find_by_uname=name)
           widget_control, id, set_value='1'
        endif
           
        if xrevent.value eq 'Disable' then begin
           widget_control, state.ext4_inp.par_row_base, sensitive=0
           name = widget_info(state.ext3_select, /uname)
           ii = where(strcompress(params.ext.name, /remove_all) eq name(0))
           if ii(0) ge 0 then params.ext(ii(0)).value = '0'
           id = widget_info(state.ext4_inp.par_row_base, find_by_uname=name)
           widget_control, id, set_value='0'
        endif
      end



;; ------------------------------------
;; Plot the dust opacity of the background
;; ------------------------------------
     'bg_plotopac_button' : begin
        ;; Switch to 'plot_kappa' mode for the postscript plotting
        plot_sed_kappa = 1
        widget_control, state.bg_dustopac_inp, get_value=txt
        txt = strsplit(txt, "'", /extract)
        openr, 1, txt(0)
        readf, 1, iformat
        readf, 1, nlam
        opac = {wave:dblarr(nlam), cabs:dblarr(nlam), csca:dblarr(nlam)}
        
        if iformat eq 1 then begin
           dum = [0d0, 0d0]
           for i=0, nlam-1 do begin
              readf, 1, dum
              opac.wave(i) = dum(0)
               opac.cabs(i) = dum(1)
            endfor
         endif
         if iformat ge 2 then begin
            dum = [0d0, 0d0, 0d0]
            for i=0, nlam-1 do begin
               readf, 1, dum
               opac.wave(i) = dum(0)
               opac.csca(i) = dum(1)
               opac.csca(i) = dum(2)
            endfor
         endif
         close, 1
         
         ;; Update the x- and y-range
         xrange_sed = [min(opac.wave), max(opac.wave)]
         yrange_sed = [min([opac.cabs]), max([opac.cabs])]

         ;; Now update the screen
         wset, state.spec_imag_drawer
         plot_kappa
      end
               
     'grid_tab_base' : begin
     end
     
     'grid_rin_inp' : begin
     end

     'grid_rout_inp' : begin
     end

     'grid_nr_inp' : begin
     end

     'grid_theta_max_inp' : begin
     end

     'grid_nr_theta_inp' : begin
     end

     'inp_wav_grid_base' : begin
     end

     'run_nr_phot_inp' : begin
     end

;; ------------------------------------
;; Run the thermal Monte-Carlo simulation
;; ------------------------------------
     'run_therm_mc' : begin
        ;; Read out all input widgets
        read_all_widget_values, state, params
        ;; Read out all input widgets in the extra tabs

        dum = uniq(params.ext.group)
        next = n_elements(dum)
        if params.ext(0).group_id eq -1 then next = 0
        if next gt 0 then read_ext_tabs, state.ext1_inp, params, group=params.ext(dum(0)).group
        if next gt 1 then read_ext_tabs, state.ext2_inp, params, group=params.ext(dum(1)).group
        if next gt 2 then read_ext_tabs, state.ext3_inp, params, group=params.ext(dum(2)).group
        if next gt 3 then read_ext_tabs, state.ext4_inp, params, group=params.ext(dum(3)).group
        
        ;; Write the problem_params.pro file
        write_params_radmc3d, params, write_file_name='problem_params.pro'
        ;; Make the model density structure
        spawn, 'idl do_setup'


        ;; Run the extra function #1 if it is enabled 
        if state.ext1_select ne -1L then begin
           widget_control, state.ext1_select, get_value=dum
           ;; Check which group/function is associated with this tab
           if dum eq 0 then begin
              dum2 = widget_info(state.ext1_select, /uname)
              group = strsplit(dum2, '_', /extract)
              group = group(1)
              
              ;; Now get the function_name
              ii = where(strcompress(params.ext.name, /remove_all) eq 'extra_'+group+'_func')
              
              if ii(0) ge 0 then begin
                 if strmid(strcompress(params.ext(ii(0)).value, /remove_all), 0, 1) eq "'" then begin
                    pro_name = strsplit(params.ext(ii).value, "'", /extract)
                 endif else pro_name = strcompress(params.ext(ii(0)).value)
                 resolve_routine, pro_name, /compile_full_file
                 call_procedure, pro_name(0)
              endif
           endif
        endif
        
        if state.ext2_select ne -1L then begin
           ;; Run the extra function #2 if it is enabled 
           widget_control, state.ext2_select, get_value=dum
           ;; Check which group/function is associated with this tab
           if dum eq 0 then begin
              dum2 = widget_info(state.ext2_select, /uname)
              group = strsplit(dum2, '_', /extract)
              group = group(1)
              
              ;; Now get the function_name
              ii = where(strcompress(params.ext.name, /remove_all) eq 'extra_'+group+'_func')
              
              if ii(0) ge 0 then begin
                 if strmid(strcompress(params.ext(ii(0)).value, /remove_all), 0, 1) eq "'" then begin
                    pro_name = strsplit(params.ext(ii).value, "'", /extract)
                 endif else pro_name = strcompress(params.ext(ii(0)).value)
                 resolve_routine, pro_name, /compile_full_file
                 call_procedure, pro_name(0)
              endif
           endif
        endif
        
        if state.ext3_select ne -1L then begin
           ;; Run the extra function #3 if it is enabled 
           widget_control, state.ext3_select, get_value=dum
           ;; Check which group/function is associated with this tab
           if dum eq 0 then begin
              dum2 = widget_info(state.ext3_select, /uname)
              group = strsplit(dum2, '_', /extract)
              group = group(1)
              
              ;; Now get the function_name
              ii = where(strcompress(params.ext.name, /remove_all) eq 'extra_'+group+'_func')
              
              if ii(0) ge 0 then begin
                 if strmid(strcompress(params.ext(ii(0)).value, /remove_all), 0, 1) eq "'" then begin
                    pro_name = strsplit(params.ext(ii).value, "'", /extract)
                 endif else pro_name = strcompress(params.ext(ii(0)).value)
                 resolve_routine, pro_name, /compile_full_file
                 call_procedure, pro_name(0)
              endif
           endif
        endif

        if state.ext4_select ne -1L then begin
           ;; Run the extra function #4 if it is enabled 
           widget_control, state.ext4_select, get_value=dum
           ;; Check which group/function is associated with this tab
           if dum eq 0 then begin
              dum2 = widget_info(state.ext4_select, /uname)
              group = strsplit(dum2, '_', /extract)
              group = group(1)
              
              ;; Now get the function_name
              ii = where(strcompress(params.ext.name, /remove_all) eq 'extra_'+group+'_func')
              
              if ii(0) ge 0 then begin
                 if strmid(strcompress(params.ext(ii(0)).value, /remove_all), 0, 1) eq "'" then begin
                    pro_name = strsplit(params.ext(ii).value, "'", /extract)
                 endif else pro_name = strcompress(params.ext(ii(0)).value)
                 resolve_routine, pro_name, /compile_full_file
                 call_procedure, pro_name(0)
              endif
           endif
        endif

        ;; Read the model density structure (ONLY DENSITY is read just
        ;; to update the maximum values of the sliders)
        dum = findfile('dust_density.inp', count=dcount)
        if dcount eq 0 then begin
        ;!!!!!!!! ERRROR!!! something went wrong!!!!
           dum = dialog_message('Somehow the dust density file was not created, although it should have been!', $
                                /error, title='ERROR')
        endif else begin
           struct = read_data(/ddens)
        endelse
        ;; Set the disk structure x-, y-range sliders maximum value to the 
        ;;   number of r, theta grid points
        widget_control, state.disk_struct_xslide, set_slider_max=struct.grid.nr-1
        widget_control, state.disk_struct_yslide, set_slider_max=struct.grid.ntheta-1
        ;; Read the dust opacity file name
        widget_control, state.ppdisk_dustopac_inp, get_value=dustopac_fname
        ;; Read the wavelength at which the optical depth should be calculated
        widget_control, state.plot_nd_tauwav_inp, get_value=tauwav
        ;; Calculate the optical depth
        tau    = get_tau_1dust(dens=struct, wav=float(tauwav), dustopac_fname=strcompress(dustopac_fname, /remove_all),$
	                               /imirt)
        ;; Run the code
        spawn, 'nice radmc3d mctherm'

        ;; Check if temperature structe was created
        dum = findfile('dust_temperature.dat', count=tcount)
        if tcount eq 0 then begin
           ;!!!!!!!! ERRROR!!! something went wrong!!!!
           dum = dialog_message('Somehow RADMC3D did not finish properly and the dust_temperature.dat file was '+$
                                'not created.', /error, title='ERROR')
        endif else begin
           ;; Re-read the model density and temperature structure
           struct = read_data(/ddens, /dtemp)
                        
           ;; Now update the disk structure drawer
           wset, state.disk_struct_drawer
           if plot_dim eq 1 then replot_1d_disk_struct
           if plot_dim eq 2 then replot_2d_disk_struct
        endelse
     end
;; ------------------------------------
;; Run a single-wavelength Monte-Carlo simulation for scattering
;;  (Currently disabled)
;; ------------------------------------
     'run_scat_mc' : begin
        dum = dialog_message('This feature is under developement (i.e. temporarily disabled)!', /error, $
                             title='ERROR')
     end

      'run_incl_inp' : begin
     end
;; ------------------------------------
;; Calculate the spectrum
;; ------------------------------------
     'run_make_spectrum' : begin
        ;; Get the inclination angle at which the spectrum should be calculated
        widget_control, state.run_incl_inp, get_value=txt
        command = 'radmc3d spectrum globallam incl ' + strcompress(txt, /remove_all)
        ;; Run RADMC3D
        spawn, command
        ;; Read the spectrum
        spec = readspectrum()
        ;; Now update the screen
        wset, state.spec_imag_drawer
        replot_sed
     end
;; ------------------------------------
;; Calculate an image
;;  TODO: I do not know at the moment if it is worth to
;;        develope a brand new stuff here or I should just
;;        link what Kees has already written
;; ------------------------------------
     'run_make_image' : begin
        dum = dialog_message('This feature is under developement (i.e. temporarily disabled)!', /error, $
                             title='ERROR')
     end
  endcase

end
