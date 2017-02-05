;
; Some natural constants
;
AU  = 1.49598d13     ; Astronomical Unit       [cm]
pc  = 3.08572d18     ; Parsec                  [cm]
MS  = 1.98892d33     ; Solar mass              [g]
TS  = 5.78d3         ; Solar temperature       [K]
LS  = 3.8525d33      ; Solar luminosity        [erg/s]
RS  = 6.96d10        ; Solar radius            [cm]

;
;>>>>>>>>>>>>> Declaration of the input variables (Do not modify this line!)<<<<<<<<
;
;
; Monte Carlo parameters
;
nphot = 100000 ; Number of photons in the Monte-Carlo simulation
;
; Star parameters
;
star_enable = 1 ; Enable a central star as a radiation source
mstar = 1.0*MS ; Mass of the star
rstar = 1.0*RS ; Radius of the star
tstar = 1.0*TS ; Temperature of the star
px = 0. ; Position of the star in the grid (x)
py = 0. ; Position of the star in the grid (y)
pz = 0. ; Position of the star in the grid (z)

;
; Continuous starlike energy source parameters
;
cstar_enable = 0 ; Enable a central star as a radiation source
cstar_rin = 1.0 *AU ; Inner radius
cstar_rout = 1.0 *AU ; Outer radius
cstar_trin = 1000.0 ; Temperature at the inner radius
cstar_powex = -0.75 ; Radial temperature distribution power exponent
;
; External (e.g. interstellar)  energy source parameters
;
ext_enable = 0 ; Enable a central star as a radiation source
ext_temp = 30d0 ; Effective temperature of the radiation field (one single bb)
;
; Disk parameters
;
ppdisk_enable = 1 ; 1 - the disk is added to the model, 0 - no disk is added to the model
ppdisk_rin = 2.01*AU ; Inner radius of the disk
ppdisk_rout = 200*AU ; Outer radius of the disk
ppdisk_mdisk = 0.001d0*MS ; Mass of the disk
ppdisk_sig0 = 1d-4 ; Surface density at the outer radius
ppdisk_plsig1 = -1d0 ; Power exponent of the radial surface density distribution
ppdisk_hrdisk = 0.15 ; Ratio of the pressure scale height to the radius at rout
ppdisk_plh = 2d0/7d0 ; Flare index
ppdisk_dustopac_file = 'dustkappa_silicate.inp' ; File name containing the dust opacity for the disk
;
; Envelope parameters
;
envelope_enable = 1 ; 1 - the envelope is added to the model, 0 - no envelope is added to the model
envelope_rin = 5*AU ; Inner radius of the envelope
envelope_rout = 450*AU ; Outer radius of the envelope
envelope_mass = 0.001d0*MS ; Mass of the envelope
envelope_rho0 = 0.5d-21 ; Mass of the envelope
envelope_plrho = -1.5 ; Radial density distribution power exponent
envelope_cavrad = 30.0 ; Radius of the polar cavity in the envelope (in degrees!)
envelope_cavrfact = 0.0 ; The density in the cavity is reduced by this factor compared to the 'no cavity' case
envelope_dustopac_file = 'dustkappa_silicate.inp' ; File name containing the dust opacity for the envelope

;
; My extra function - example 
;   Cut a gap in the disk reducing the density by an arbitrary factor
;
extra_gap_tab_name = 'My gap function' ; Name of the tab in the GUI
extra_gap_enable = 1 ; 0 - Enable, 1 - Disable this function
extra_gap_func = 'my_gap_function' ; Name of the user defined function to be executed (after the problem_setup.pro)
extra_gap_rin = 5.0*AU ; Inner radius of the gap
extra_gap_rout = 10.*AU ; Outer radius of the gap
extra_gap_drfact = 1d-4 ; Density reduction factor in the gap (with respect to the original density)

;
; My extra function - example #2
;   Cut a gap in the disk reducing the density by an arbitrary factor
;
extra_gap2_tab_name = 'My gap function 2' ; Name of the tab in the GUI
extra_gap2_enable = 1 ; 0 - Enable, 1 - Disable this function
extra_gap2_func = 'my_gap2_function' ; Name of the user defined function to be executed (after the problem_setup.pro)
extra_gap2_rin = 20.0*AU ; Inner radius of the gap
extra_gap2_rout = 30.*AU ; Outer radius of the gap
extra_gap2_drfact = 1d-4 ; Density reduction factor in the gap (with respect to the original density)


;
; Background parameters
;
bg_enable = 1 ; 1 - constant background density is added to the model, 0 - no background density is added
bg_rho = 1d-22 ; Background density
bg_dustopac_file = 'dustkappa_silicate.inp' ; File name containing the dust opacity for the background
;
; Wavelength grid parameters
;
lambda1 = 0.1d0
lambda2 = 5.0d0
lambda3 = 40.0d0
lambda4 = 10000.0d0
n12 = 20
n23 = 20
n34 = 20
;
; Spatial grid parameters (spherical)
;
grid_rin = 2*AU ; Inner radius of the grid
grid_rout = 500*AU ; Outer radius of the grid
grid_nr = 100L ; Number of radial grid points
grid_tmax = 0.9 ; Upper boundary of theta (meridional angular coordinate)
grid_nt = 60L ; Number of theta grid points One-sided only. So the "real" value is twice this.
grid_np = 1L ; Number of phi grid points

;
;>>>>>>>>>>>>> End of input variable declaration (Do not modify this line!)<<<<<
;

