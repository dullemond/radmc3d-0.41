;; **********************************************************************************************
;; Function to create a tab in the dust densit distribution for an
;; arbitrary external module
;;
;; NOTE : The name of the parameter, whose value can be modified
;; in an imput widget, is stored both in uname and in uvalue
;; **********************************************************************************************
function create_ext_tab2, parent=parent, params=params, group=group
;;
;; Set up some defaults/rules
;;
ninp_max_left = 5 ;; Maximum number of input field widgets in the left column
ninp_max_right = 2 ;; Maximum number of input field widgets in the right column
pname_length_max = 18 ;; Maximum length of a parameter name in terms of characters
;;
;; Count the number of parameters in this group
;;
jj = where(params.ext.group eq group(0) and strcompress(params.ext.name, /remove_all) ne 'extra_'+group+'_tab_name')
npar = n_elements(jj) ;; Number of parameters in the group
nshow_par = min([npar, 7]) ;; Number of showed parameters directly on the screen
;;
;; Now build the tab
;;

bspace = " "
par_row_base = WIDGET_BASE(parent, /row, frame=0)
par_left_column_base = WIDGET_BASE(par_row_base, /column, frame=0)
par_right_column_base = WIDGET_BASE(par_row_base, /column, frame=0)

;; ------------------------------------------------------------------------------------------------
;; Par 1
;; ------------------------------------------------------------------------------------------------
if nshow_par gt 0 then begin
   ipar = 0
   dum = strtrim(params.ext(jj(ipar)).name, 2)
   pname_length = strlen(dum)

;; Create a parameter name with a fixed string length

   tagname = dum
   if pname_length lt pname_length_max then begin
      dum2 = strarr(pname_length_max-pname_length)
      dum2(*) = bspace
      pname = strjoin([dum, dum2], /single)
   endif else begin
      if pname_length eq pname_length_max then begin
         pname = dum
      endif else begin
         pname = strmid(dum, 0, pname_length_max)
      endelse
   endelse
;; Create the input widget

;; Check if the parameter is has an extra_groupname_enable-type name
;; and if so make the input field non-editable
   
   dum = strmid(strcompress(pname,/remove_all),5,6,/reverse)
   
   if dum eq 'enable' then begin
      par1_inp = CW_FIELD(par_left_column_base, title=pname, uvalue=tagname, uname=tagname, $
                          value='', xsize=8, /row, /noedit)
   endif else begin
      par1_inp = CW_FIELD(par_left_column_base, title=pname, uvalue=tagname, uname=tagname, $
                          value='', xsize=8, /row)
   endelse
endif else par1_inp = -1L

;; ------------------------------------------------------------------------------------------------
;; Par 2
;; ------------------------------------------------------------------------------------------------
if nshow_par gt 1 then begin
   ipar = 1
   dum = strtrim(params.ext(jj(ipar)).name, 2)
   pname_length = strlen(dum)

;; Create a parameter name with a fixed string length

   tagname = dum
   if pname_length lt pname_length_max then begin
      dum2 = strarr(pname_length_max-pname_length)
      dum2(*) = bspace
      pname = strjoin([dum, dum2], /single)
   endif else begin
      if pname_length eq pname_length_max then begin
         pname = dum
      endif else begin
         pname = strmid(dum, 0, pname_length_max)
      endelse
   endelse
;; Create the input widget

;; Check if the parameter is has an extra_groupname_enable-type name
;; and if so make the input field non-editable
   
   dum = strmid(strcompress(pname,/remove_all),5,6,/reverse)
   
   if dum eq 'enable' then begin
      par2_inp = CW_FIELD(par_left_column_base, title=pname, uvalue=tagname, uname=tagname, $
                          value='', xsize=8, /row, /noedit)
   endif else begin
      par2_inp = CW_FIELD(par_left_column_base, title=pname, uvalue=tagname, uname=tagname, $
                          value='', xsize=8, /row)
   endelse
   
endif else par2_inp = -1L


;; ------------------------------------------------------------------------------------------------
;; Par 3
;; ------------------------------------------------------------------------------------------------
if nshow_par gt 2 then begin
   ipar = 2
   dum = strtrim(params.ext(jj(ipar)).name, 2)
   pname_length = strlen(dum)

;; Create a parameter name with a fixed string length

   tagname = dum
   if pname_length lt pname_length_max then begin
      dum2 = strarr(pname_length_max-pname_length)
      dum2(*) = bspace
      pname = strjoin([dum, dum2], /single)
   endif else begin
      if pname_length eq pname_length_max then begin
         pname = dum
      endif else begin
         pname = strmid(dum, 0, pname_length_max)
      endelse
   endelse
;; Create the input widget

;; Check if the parameter is has an extra_groupname_enable-type name
;; and if so make the input field non-editable
   
   dum = strmid(strcompress(pname,/remove_all),5,6,/reverse)
   
   if dum eq 'enable' then begin
      par3_inp = CW_FIELD(par_left_column_base, title=pname, uvalue=tagname, uname=tagname, $
                          value='', xsize=8, /row, /noedit)
   endif else begin
      par3_inp = CW_FIELD(par_left_column_base, title=pname, uvalue=tagname, uname=tagname, $
                          value='', xsize=8, /row)
   endelse

   
endif else par3_inp = -1L

;; ------------------------------------------------------------------------------------------------
;; Par 4
;; ------------------------------------------------------------------------------------------------
if nshow_par gt 3 then begin
   ipar = 3
   dum = strtrim(params.ext(jj(ipar)).name, 2)
   pname_length = strlen(dum)

;; Create a parameter name with a fixed string length

   tagname = dum
   if pname_length lt pname_length_max then begin
      dum2 = strarr(pname_length_max-pname_length)
      dum2(*) = bspace
      pname = strjoin([dum, dum2], /single)
   endif else begin
      if pname_length eq pname_length_max then begin
         pname = dum
      endif else begin
         pname = strmid(dum, 0, pname_length_max)
      endelse
   endelse
;; Create the input widget

;; Check if the parameter is has an extra_groupname_enable-type name
;; and if so make the input field non-editable
   
   dum = strmid(strcompress(pname,/remove_all),5,6,/reverse)
   
   if dum eq 'enable' then begin
      par4_inp = CW_FIELD(par_left_column_base, title=pname, uvalue=tagname, uname=tagname, $
                          value='', xsize=8, /row, /noedit)
   endif else begin
      par4_inp = CW_FIELD(par_left_column_base, title=pname, uvalue=tagname, uname=tagname, $
                          value='', xsize=8, /row)
   endelse


endif else par4_inp = -1L

;; ------------------------------------------------------------------------------------------------
;; Par 5
;; ------------------------------------------------------------------------------------------------
if nshow_par gt 4 then begin
   ipar = 4
   dum = strtrim(params.ext(jj(ipar)).name, 2)
   pname_length = strlen(dum)

;; Create a parameter name with a fixed string length

   tagname = dum
   if pname_length lt pname_length_max then begin
      dum2 = strarr(pname_length_max-pname_length)
      dum2(*) = bspace
      pname = strjoin([dum, dum2], /single)
   endif else begin
      if pname_length eq pname_length_max then begin
         pname = dum
      endif else begin
         pname = strmid(dum, 0, pname_length_max)
      endelse
   endelse
;; Create the input widget

;; Check if the parameter is has an extra_groupname_enable-type name
;; and if so make the input field non-editable
   
   dum = strmid(strcompress(pname,/remove_all),5,6,/reverse)
   
   if dum eq 'enable' then begin
      par5_inp = CW_FIELD(par_left_column_base, title=pname, uvalue=tagname, uname=tagname, $
                          value='', xsize=8, /row, /noedit)
   endif else begin
      par5_inp = CW_FIELD(par_left_column_base, title=pname, uvalue=tagname, uname=tagname, $
                          value='', xsize=8, /row)
   endelse


   
endif else par5_inp = -1L

;; ------------------------------------------------------------------------------------------------
;; Par 6
;; ------------------------------------------------------------------------------------------------
if nshow_par gt 5 then begin
   ipar = 5
   dum = strtrim(params.ext(jj(ipar)).name, 2)
   pname_length = strlen(dum)

;; Create a parameter name with a fixed string length

   tagname = dum
   if pname_length lt pname_length_max then begin
      dum2 = strarr(pname_length_max-pname_length)
      dum2(*) = bspace
      pname = strjoin([dum, dum2], /single)
   endif else begin
      if pname_length eq pname_length_max then begin
         pname = dum
      endif else begin
         pname = strmid(dum, 0, pname_length_max)
      endelse
   endelse
;; Create the input widget

;; Check if the parameter is has an extra_groupname_enable-type name
;; and if so make the input field non-editable
   
   dum = strmid(strcompress(pname,/remove_all),5,6,/reverse)
   
   if dum eq 'enable' then begin
      par6_inp = CW_FIELD(par_right_column_base, title=pname, uvalue=tagname, uname=tagname, $
                          value='', xsize=8, /row, /noedit)
   endif else begin
      par6_inp = CW_FIELD(par_right_column_base, title=pname, uvalue=tagname, uname=tagname, $
                          value='', xsize=8, /row)
   endelse

   
endif else par6_inp = -1L

;; ------------------------------------------------------------------------------------------------
;; Par 7
;; ------------------------------------------------------------------------------------------------
if nshow_par gt 6 then begin
   ipar = 6
   dum = strtrim(params.ext(jj(ipar)).name, 2)
   pname_length = strlen(dum)

;; Create a parameter name with a fixed string length

   tagname = dum
   if pname_length lt pname_length_max then begin
      dum2 = strarr(pname_length_max-pname_length)
      dum2(*) = bspace
      pname = strjoin([dum, dum2], /single)
   endif else begin
      if pname_length eq pname_length_max then begin
         pname = dum
      endif else begin
         pname = strmid(dum, 0, pname_length_max)
      endelse
   endelse
;; Create the input widget

;; Check if the parameter is has an extra_groupname_enable-type name
;; and if so make the input field non-editable
   
   dum = strmid(strcompress(pname,/remove_all),5,6,/reverse)
   
   if dum eq 'enable' then begin
      par7_inp = CW_FIELD(par_right_column_base, title=pname, uvalue=tagname, uname=tagname, $
                          value='', xsize=8, /row, /noedit)
   endif else begin
      par7_inp = CW_FIELD(par_right_column_base, title=pname, uvalue=tagname, uname=tagname, $
                          value='', xsize=8, /row)
   endelse


endif else par7_inp = -1L

;; ------------------------------------------------------------------------------------------------
;; Par 8
;; ------------------------------------------------------------------------------------------------
if nshow_par gt 7 then begin
   ipar = 7
   dum = strtrim(params.ext(jj(ipar)).name, 2)
   pname_length = strlen(dum)

;; Create a parameter name with a fixed string length

   tagname = dum
   if pname_length lt pname_length_max then begin
      dum2 = strarr(pname_length_max-pname_length)
      dum2(*) = bspace
      pname = strjoin([dum, dum2], /single)
   endif else begin
      if pname_length eq pname_length_max then begin
         pname = dum
      endif else begin
         pname = strmid(dum, 0, pname_length_max)
      endelse
   endelse
;; Create the input widget

;; Check if the parameter is has an extra_groupname_enable-type name
;; and if so make the input field non-editable
   
   dum = strmid(strcompress(pname,/remove_all),5,6,/reverse)
   
   if dum eq 'enable' then begin
      par8_inp = CW_FIELD(par_right_column_base, title=pname, uvalue=tagname, uname=tagname, $
                          value='', xsize=8, /row, /noedit)
   endif else begin
      par8_inp = CW_FIELD(par_right_column_base, title=pname, uvalue=tagname, uname=tagname, $
                          value='', xsize=8, /row)
   endelse

   
endif else par8_inp = -1L

;; ------------------------------------------------------------------------------------------------
;; Par 9
;; ------------------------------------------------------------------------------------------------
if nshow_par gt 8 then begin
   ipar = 8
   dum = strtrim(params.ext(jj(ipar)).name, 2)
   pname_length = strlen(dum)

;; Create a parameter name with a fixed string length

   tagname = dum
   if pname_length lt pname_length_max then begin
      dum2 = strarr(pname_length_max-pname_length)
      dum2(*) = bspace
      pname = strjoin([dum, dum2], /single)
   endif else begin
      if pname_length eq pname_length_max then begin
         pname = dum
      endif else begin
         pname = strmid(dum, 0, pname_length_max)
      endelse
   endelse
;; Create the input widget

;; Check if the parameter is has an extra_groupname_enable-type name
;; and if so make the input field non-editable
   
   dum = strmid(strcompress(pname,/remove_all),5,6,/reverse)
   
   if dum eq 'enable' then begin
      par9_inp = CW_FIELD(par_right_column_base, title=pname, uvalue=tagname, uname=tagname, $
                          value='', xsize=8, /row, /noedit)
   endif else begin
      par9_inp = CW_FIELD(par_right_column_base, title=pname, uvalue=tagname, uname=tagname, $
                          value='', xsize=8, /row)
   endelse
   
endif else par9_inp = -1L

;; ------------------------------------------------------------------------------------------------
;; Par 10
;; ------------------------------------------------------------------------------------------------
if nshow_par gt 9 then begin
   ipar = 9
   dum = strtrim(params.ext(jj(ipar)).name, 2)
   pname_length = strlen(dum)

;; Create a parameter name with a fixed string length

   tagname = dum
   if pname_length lt pname_length_max then begin
      dum2 = strarr(pname_length_max-pname_length)
      dum2(*) = bspace
      pname = strjoin([dum, dum2], /single)
   endif else begin
      if pname_length eq pname_length_max then begin
         pname = dum
      endif else begin
         pname = strmid(dum, 0, pname_length_max)
      endelse
   endelse
;; Create the input widget

;; Check if the parameter is has an extra_groupname_enable-type name
;; and if so make the input field non-editable
   
   dum = strmid(strcompress(pname,/remove_all),5,6,/reverse)
   
   if dum eq 'enable' then begin
      par10_inp = CW_FIELD(par_right_column_base, title=pname, uvalue=tagname, uname=tagname, $
                          value='', xsize=8, /row, /noedit)
   endif else begin
      par10_inp = CW_FIELD(par_right_column_base, title=pname, uvalue=tagname, uname=tagname, $
                          value='', xsize=8, /row)
   endelse
   
endif else par10_inp = -1L


;; ------------------------------------------------------------------------------------------------
;; Create a structure to return
;; ------------------------------------------------------------------------------------------------

ext = {par_row_base:par_row_base, par_left_column_base:par_left_column_base, $
       par_right_column_base:par_right_column_base, $
       par1_inp:par1_inp, par2_inp:par2_inp, par3_inp:par3_inp, par4_inp:par4_inp, par5_inp:par5_inp, $
       par6_inp:par6_inp, par7_inp:par7_inp, par8_inp:par8_inp, par9_inp:par9_inp, par10_inp:par10_inp}


;;
;; Return the whole structure
;;

skip_to_this:

return, ext
end

;; **********************************************************************************************
;; A procedure to fill up all the extra input widgets we created in the density
;; setup field
;;
;; How it works : 
;; Go trhough all the text widgets. If a widget has been created and
;; is shown in the GUI (i.e. ID ne -1L) then get the uvalue which
;; corresponds to the parameter name the value of which is in the text
;; field. Then find this field in the params structure and modify the
;; value. 
;; **********************************************************************************************
pro fillup_ext_tabs, ext, params, group=group

  if ext.par1_inp ne -1 then begin
     widget_control, ext.par1_inp, get_uvalue=pname
     ii = where(strtrim(params.ext.name, 2) eq pname)
     if ii(0) ge 0 then widget_control, ext.par1_inp, set_value=params.ext(ii(0)).value
  endif

  if ext.par2_inp ne -1 then begin 
     widget_control, ext.par2_inp, get_uvalue=pname
     ii = where(strtrim(params.ext.name, 2) eq pname)
     if ii(0) ge 0 then widget_control, ext.par2_inp, set_value=params.ext(ii(0)).value
  endif
  
  if ext.par3_inp ne -1 then begin
     widget_control, ext.par3_inp, get_uvalue=pname
     ii = where(strtrim(params.ext.name, 2) eq pname)
     if ii(0) ge 0 then widget_control, ext.par3_inp, set_value=params.ext(ii(0)).value
  endif
  if ext.par4_inp ne -1 then begin
     widget_control, ext.par4_inp, get_uvalue=pname
     ii = where(strtrim(params.ext.name, 2) eq pname)
     if ii(0) ge 0 then widget_control, ext.par4_inp, set_value=params.ext(ii(0)).value
  endif
  
  if ext.par5_inp ne -1 then begin
     widget_control, ext.par5_inp, get_uvalue=pname
     ii = where(strtrim(params.ext.name, 2) eq pname)
     if ii(0) ge 0 then widget_control, ext.par5_inp, set_value=params.ext(ii(0)).value
  endif
  
  if ext.par6_inp ne -1 then begin
     widget_control, ext.par6_inp, get_uvalue=pname
     ii = where(strtrim(params.ext.name, 2) eq pname)
     if ii(0) ge 0 then widget_control, ext.par6_inp, set_value=params.ext(ii(0)).value
  endif
  
  if ext.par7_inp ne -1 then begin
     widget_control, ext.par7_inp, get_uvalue=pname
     ii = where(strtrim(params.ext.name, 2) eq pname)
     if ii(0) ge 0 then widget_control, ext.par7_inp, set_value=params.ext(ii(0)).value
  endif
  
  if ext.par8_inp ne -1 then begin
     widget_control, ext.par8_inp, get_uvalue=pname
     ii = where(strtrim(params.ext.name, 2) eq pname)
     if ii(0) ge 0 then widget_control, ext.par8_inp, set_value=params.ext(ii(0)).value
  endif
  
  if ext.par9_inp ne -1 then begin
     widget_control, ext.par9_inp, get_uvalue=pname
     ii = where(strtrim(params.ext.name, 2) eq pname)
     if ii(0) ge 0 then widget_control, ext.par9_inp, set_value=params.ext(ii(0)).value
  endif
  
  if ext.par10_inp ne -1 then begin
     widget_control, ext.par10_inp, get_uvalue=pname
     ii = where(strtrim(params.ext.name, 2) eq pname)
     if ii(0) ge 0 then widget_control, ext.par10_inp, set_value=params.ext(ii(0)).value
  endif

end

;; **********************************************************************************************
;; A procedure to fill up all the extra input widgets we created in the density
;; setup field
;; **********************************************************************************************
pro read_ext_tabs, ext, params, group=group

  if ext.par1_inp ne -1 then begin
     widget_control, ext.par1_inp, get_uvalue=pname
     ii = where(strtrim(params.ext.name, 2) eq pname)
     if ii(0) ge 0 then widget_control, ext.par1_inp, get_value=dum 
     params.ext(ii(0)).value = strtrim(dum, 2)
  endif

  if ext.par2_inp ne -1 then begin 
     widget_control, ext.par2_inp, get_uvalue=pname
     ii = where(strtrim(params.ext.name, 2) eq pname)
     if ii(0) ge 0 then widget_control, ext.par2_inp, get_value=dum 
     params.ext(ii(0)).value = strtrim(dum, 2)
  endif
  
  if ext.par3_inp ne -1 then begin
     widget_control, ext.par3_inp, get_uvalue=pname
     ii = where(strtrim(params.ext.name, 2) eq pname)
     if ii(0) ge 0 then widget_control, ext.par3_inp, get_value=dum 
     params.ext(ii(0)).value = strtrim(dum, 2)
  endif
  if ext.par4_inp ne -1 then begin
     widget_control, ext.par4_inp, get_uvalue=pname
     ii = where(strtrim(params.ext.name, 2) eq pname)
     if ii(0) ge 0 then widget_control, ext.par4_inp, get_value=dum 
     params.ext(ii(0)).value = strtrim(dum, 2) 
  endif
  
  if ext.par5_inp ne -1 then begin
     widget_control, ext.par5_inp, get_uvalue=pname
     ii = where(strtrim(params.ext.name, 2) eq pname)
     if ii(0) ge 0 then widget_control, ext.par5_inp, get_value=dum 
     params.ext(ii(0)).value = strtrim(dum, 2) 
  endif
  
  if ext.par6_inp ne -1 then begin
     widget_control, ext.par6_inp, get_uvalue=pname
     ii = where(strtrim(params.ext.name, 2) eq pname)
     if ii(0) ge 0 then widget_control, ext.par6_inp, get_value=dum 
     params.ext(ii(0)).value = strtrim(dum, 2)
  endif
  
  if ext.par7_inp ne -1 then begin
     widget_control, ext.par7_inp, get_uvalue=pname
     ii = where(strtrim(params.ext.name, 2) eq pname)
     if ii(0) ge 0 then widget_control, ext.par7_inp, get_value=dum 
     params.ext(ii(0)).value = strtrim(dum, 2)
  endif
  
  if ext.par8_inp ne -1 then begin
     widget_control, ext.par8_inp, get_uvalue=pname
     ii = where(strtrim(params.ext.name, 2) eq pname)
     if ii(0) ge 0 then widget_control, ext.par8_inp, get_value=dum 
     params.ext(ii(0)).value = strtrim(dum, 2)
  endif
  
  if ext.par9_inp ne -1 then begin
     widget_control, ext.par9_inp, get_uvalue=pname
     ii = where(strtrim(params.ext.name, 2) eq pname)
     if ii(0) ge 0 then widget_control, ext.par9_inp, get_value=dum 
     params.ext(ii(0)).value = strtrim(dum, 2) 
  endif
  
  if ext.par10_inp ne -1 then begin
     widget_control, ext.par10_inp, get_uvalue=pname
     ii = where(strtrim(params.ext.name, 2) eq pname)
     if ii(0) ge 0 then widget_control, ext.par10_inp, get_value=dum 
     params.ext(ii(0)).value = strtrim(dum, 2) 
  endif


end
