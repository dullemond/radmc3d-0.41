pro my_gap2_function
@problem_params.pro
;;=============================================================================
;; Procedure to cut a gap in the disk density distribution
;;=============================================================================
struct = read_data(/ddens)
rho_new = struct.rho
;;
;; Create a gap in the disk
;;
for ir=0, struct.grid.nr-1 do begin
   if struct.grid.r[ir] ge extra_gap2_rin and struct.grid.r[ir] le extra_gap2_rout then $
      rho_new[ir,*] = struct.rho[ir,*] * extra_gap2_drfact
endfor


;;
;; If the background density is enabled then set the density to the
;; background value wherever the reduced density of the disk falls
;; below the background value
;;
if bg_enable eq 1 then begin
   for ir=0, struct.grid.nr-1 do begin
      ii = where(rho_new[ir,*] lt bg_rho)
      if ii(0) ge 0 then rho_new[ir,ii] = bg_rho
   endfor
endif

;;
;; Write out the new dust density distribution
;;
openw,1,'dust_density.inp'
printf,1,1                                ; Format number
print, grid_nr*grid_nt*grid_np
printf,1,grid_nr*grid_nt*grid_np          ; Nr of cells
printf,1,1                                ; Nr of dust species
for ip=0,grid_np-1 do begin
   for it=0,grid_nt-1 do begin
      for ir=0,grid_nr-1 do begin
         ;printf,1,rho[ip,it,ir]
;
; TODO : I need to change the order of indices in rho
;      to be [phi, theta, r] as originally used by Kees
;
         printf,1,rho_new[ir,it]
      endfor
   endfor
endfor
close,1

end
