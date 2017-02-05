;=============================================================================
;             IMAGE TO FITS OR FITS TO IMAGE ROUTINES FOR RADMC-3D
;
; WARNING: These routines require the astrolib routines readfits() and
;          writefits(). Web site: http://idlastro.gsfc.nasa.gov/
;       
;=============================================================================
@readradmc
@readfits
@writefits

;-----------------------------------------------------------------------------
;                   ROUTINE FOR WRITING FITS IMAGE
;-----------------------------------------------------------------------------
pro writefitsimage,image,filename,freqhz,pixdeg_x,pixdeg_y,radeg,decdeg,$
                   arcsec=arcsec,mas=mas
nx = n_elements(image[*,0])
ny = n_elements(image[0,*])
;
; Compute the conversion factor from erg/cm^2/s/Hz/ster to Jy/pixel
;
pixsurf_ster = pixdeg_x*pixdeg_y/(!RADEG^2)
factor = 1d+23 * pixsurf_ster
;
; Make FITS header information, reverse order
;
header_in = strarr(1)
header_in[0] = 'END     '
;
; ...Rest frequency
;
str = "RESTFREQ=  "+string(freqhz,format='(E19.12)')
header_in = [str,header_in]
;
; ...Zero point of coordinate system
;
str = "CRPIX2  ="+string((ny+1)/2,format='(E21.13)')
header_in = [str,header_in]
;
; ...Pixel scale
;
if keyword_set(arcsec) then begin
   str = "CUNIT2  = 'arcsec  '"
   header_in = [str,header_in]
   str = "CDELT2  ="+string(3.6d3*pixdeg_y,format='(E21.13)')
   header_in = [str,header_in]
endif else begin
   if keyword_set(mas) then begin
      str = "CUNIT2  = 'mas     '"
      header_in = [str,header_in]
      str = "CDELT2  ="+string(3.6d6*pixdeg_y,format='(E21.13)')
      header_in = [str,header_in]
   endif else begin
      str = "CUNIT2  = 'deg     '"
      header_in = [str,header_in]
      str = "CDELT2  ="+string(pixdeg_y,format='(E21.13)')
      header_in = [str,header_in]
   endelse
endelse
;
; ...CRVAL2: value of y-axis at critical pixel
;
if keyword_set(decdeg) then begin
   str = "CRVAL2  =  "+string(decdeg,format='(E19.12)')
   header_in = [str,header_in]
   str = "CTYPE2  = 'DEC--SIN'"
   header_in = [str,header_in]
endif
;
; ...Zero point of coordinate system
;
str = "CRPIX1  ="+string((nx+1)/2,format='(E21.13)')
header_in = [str,header_in]
;
; ...Pixel scale
;
if keyword_set(arcsec) then begin
   str = "CUNIT1  = 'arcsec  '"
   header_in = [str,header_in]
   str = "CDELT1  ="+string(3.6d3*pixdeg_x,format='(E21.13)')
   header_in = [str,header_in]
endif else begin
   if keyword_set(mas) then begin
      str = "CUNIT1  = 'mas     '"
      header_in = [str,header_in]
      str = "CDELT1  ="+string(3.6d6*pixdeg_x,format='(E21.13)')
      header_in = [str,header_in]
   endif else begin
      str = "CUNIT1  = 'deg     '"
      header_in = [str,header_in]
      str = "CDELT1  ="+string(pixdeg_x,format='(E21.13)')
      header_in = [str,header_in]
   endelse
endelse
;
; ...CRVAL1: value of x-axis at critical pixel
;
if keyword_set(radeg) then begin
   str = "CRVAL1  =  "+string(radeg,format='(E19.12)')
   header_in = [str,header_in]
   str = "CTYPE1  = 'RA---SIN'"
   header_in = [str,header_in]
   str = "LONPOLE =   1.800000000000E+02"
   header_in = [str,header_in]
   str = "EPOCH   =   2.000000000000E+03"
   header_in = [str,header_in]
endif
;
; ...Unit of intensity
;
str = "BUNIT   = 'JY/PIXEL'"
header_in = [str,header_in]
;
; ...BZERO
;
str = "BZERO   =   0.000000000000E+00"
header_in = [str,header_in]
;
; ...BSCALE
;
str = "BSCALE   =   1.000000000000E+00"
header_in = [str,header_in]
;
; ...Type of data
;
str = "BTYPE   = 'Intensity'"
header_in = [str,header_in]
;
; Make a FITS file
;
imjypix = factor * image
;fits_write,'./image.fits',imjypix,header_in
fits_write,filename,imjypix,header_in
;
end


;-----------------------------------------------------------------------------
;                              IMAGE TO FITS
;
; RADMC-3D Images are text files (the standard is image.out). This routine
; converts such a file to FITS standard.
;
; ARGUMENTS:
;  imagename               Name of the image file that RADMC-3D produces.
;                          The standard file name RADMC-3D produces is:
;                          'image.out'.
;  fitsname                Name of the output fits file you want to produce.
;  dpc                     Distance of observer to object in units of parsec.
; 
;-----------------------------------------------------------------------------
pro radmcimage_to_fits,imagename,fitsname,dpc,arcsec=arcsec,mas=mas
  cc  = 2.9979245800000d10      ; Light speed             [cm/s]
  pc  = 3.08572d18              ; Parsec                  [cm]
  im  = readimage(file=imagename)
  pixdeg_x = 180.d0*(im.sizepix_x/(dpc*pc))/!dpi
  pixdeg_y = 180.d0*(im.sizepix_y/(dpc*pc))/!dpi
  writefitsimage,im.image,fitsname,1d4*cc/im.lambda,pixdeg_x,pixdeg_y,$
                 arcsec=arcsec,mas=mas
end


