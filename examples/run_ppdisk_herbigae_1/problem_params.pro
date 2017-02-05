;==========================================================================
;
;                  PARAMETERS FOR THE 2-D DISK MODELS
;
;==========================================================================

;--------------------------------------------------------------------------
;                       Main switch of the simulation
;--------------------------------------------------------------------------
imakedisk = 1           ; =1  Make the disk structure with or without
                        ;     vertical structure iteration. 
                        ; =-1 Use existing dustdens_struct.inp, add abundances
show = 1

;--------------------------------------------------------------------------
;                         The central star
;--------------------------------------------------------------------------
istar  = 10 
case istar of
   0: begin
      rstar  = 2.0d0*RS         ; Radius of the star
      mstar  = 2.5d0*MS         ; Mass of the star
      tstar  = 10500d0          ; Temperature of the star
      kurucz = 1                ; Use Kurucz model (1) or blackbody (0)
   end
   10: begin
      rstar  = 2.0d0*RS         ; Radius of the star
      mstar  = 2.5d0*MS         ; Mass of the star
      tstar  = 10500d0          ; Temperature of the star
      kurucz = 0                ; Use Kurucz model (1) or blackbody (0)
   end
   1: begin
      rstar  = 2.0d0*RS         ; Radius of the star
      mstar  = 1.0d0*MS         ; Mass of the star
      tstar  = 4000d0           ; Temperature of the star
      kurucz = 1                ; Use Kurucz model (1) or blackbody (0)
   end
   11: begin
      rstar  = 2.0d0*RS         ; Radius of the star
      mstar  = 1.0d0*MS         ; Mass of the star
      tstar  = 4000d0           ; Temperature of the star
      kurucz = 0                ; Use Kurucz model (1) or blackbody (0)
   end
   2: begin
      rstar  = 1.0d0*RS         ; Radius of the star
      mstar  = 0.1d0*MS         ; Mass of the star
      tstar  = 3000d0           ; Temperature of the star
      kurucz = 0                ; Use Kurucz model (1) or blackbody (0)
   end
   else: stop
endcase
ifinstar = 1                   ;  Include finite size of the star in calc.

;--------------------------------------------------------------------------
;                 The opacity table and frequency grid
;            (NOTE: See end of this file for mixed opacities)
;--------------------------------------------------------------------------
abunc    = 0.25           ; Mix carbon
mixnames = ['dustkappa_amorph_mix.inp',$
            'dustkappa_cryst_mix.inp']
mixspecs = [['dustkappa_carbon.inp','dustkappa_silicate.inp'],$
            ['dustkappa_carbon.inp','dustkappa_forsterite.inp']]
mixabun  = [[abunc,(1.0-abunc)],$
            [abunc,(1.0-abunc)]]
fresmd = 30             ; Frequency resolution/grid mode
opacs = $              ; Dust opacity input file
  ['amorph_mix','cryst_mix']
pll    = [-1,-1]        ; Power law index at long wavelengths
scat   = 1              ; Include scattering (0=no scattering)?

;--------------------------------------------------------------------------
;                      Some general parameters
;--------------------------------------------------------------------------
gastodust = 100.        ; Gas-to-dust ratio
rhofloor  = 1d-26       ; Minimum density (gas density)
run       = 0           ; Do a run within IDL? If 0, then only setup.

;--------------------------------------------------------------------------
;                         The spatial grid 
;--------------------------------------------------------------------------
nr     = 130            ; number of points in radius
nt     = 50             ; Number of points in theta for main theta grid
ntex   = 10             ; Extra points in theta 
rin    = 0.0*AU         ; Inner radius of radial grid (--> tin is computed)
tin    = 1500.d0        ; Inner temperature (--> rin is computed)
rout   = 200.d0*AU      ; Outer radius of radial grid
hrgrid = 0.5            ; Main grid span in Theta (measured negatively from equator)
hrgmax = 0.9*!pi/2.d0   ; Total grid span in Theta (measured negatively from equator)
thintin= 0              ; If 0, use DDN01 model for calc Rin from Tin, else
                        ; use real opacity for doing so.

;--------------------------------------------------------------------------
;      Grid refinement and smoothing of solution near inner edge
;--------------------------------------------------------------------------
rrefine = {nlevr:3,$    ; higher resolution at level #
           nspanr:3,$   ; number of steps to resolve
           nstepr:3}    ; number of sub-steps to create
drsm   = 0.03           ; smoothing radius

;--------------------------------------------------------------------------
;                    The Sigma(R) setup parameters
;--------------------------------------------------------------------------
rdisk  = 150 * AU       ; Pivot of power law definitions (often used as rdisk)
sig0   = 0.0            ; Sigma at rdisk (either sig0 or mdisk)
mdmstr = 1d-4           ; Mass of the disk (gas+dust) in units of Mstar
mdisk  = mdmstr * mstar ; Mass of the disk (gas+dust) (either sig0 or mdisk)
plsig1 = -1.0d0		; Powerlaw for Sigma(R) for R<Rdisk
plsig2 = -12.0d0	; Powerlaw for Sigma(R) for R>Rdisk
bgdens = 1.d-40         ; Background density in g/cm^3 (can be tau>1 !!!)

;--------------------------------------------------------------------------
;          Settings for radial mixing model (if multi species used)
;--------------------------------------------------------------------------
schm   = 1.0/3.0        ; Schmidt number
ab_r0  = 1.0*AU         ; Abundance function pivot radius
ab_ab0 = [1.0]          ; Array of nspec-1 elem, for abund0 of species 2,3,4...
ab_pl  = [-3.*schm/2.]  ; Array of nspec-1 elem, for abun pl of species 2,3,4...
ab_min = [1d-6]         ; Array of nspec-1 elem, for min abun of species 2,3,4...

;--------------------------------------------------------------------------
;           Settings for H(R) (initial value, minimal value)
;       (will be washed-out if vertical structure is computed)
;--------------------------------------------------------------------------
hrdisk = 0.15           ; Vertical pressure scale height in units of radius
hrmin  = 0.01           ; Lower limit of vertical pressure scale height
plh    = (1./7.)        ; Powerlaw for H(R) in initial setup

;--------------------------------------------------------------------------
;                  Artificial puffed-up inner rim
;       (will be washed-out if vertical structure is computed)
;--------------------------------------------------------------------------
rpfrin = 0.0            ; Start of the puffed-up inner rim: rpuff=rpfrin*rin
hrpuff = 0.0            ; Puffing-up H/R

;--------------------------------------------------------------------------
;           Parameters for RADMC : The radiative transfer part
;--------------------------------------------------------------------------
nphot  = 100000L         ; Nr of photons for RADMC
;npdiff = 15             ; Minimum nr of photons per cell --> diffusion alg
;errtol = 1d-6           ; Error tolerance for diffusion alg

;;--------------------------------------------------------------------------
;;           Parameters for RADMC : The vertical structure part
;;--------------------------------------------------------------------------
;nvstr  = 0              ; Nr of iterations of structure
;vserrt = 0.01           ; Error tolerance for vertical structure
;ivstrt = 1              ; Use dust species 1 to represent T_gas
;dostr  = [1,1]          ; Which species participate in vertical structure?

