;************************************************************************************************
; Functin to read the input parameters of RADMC3D from problem_params.pro 
;    Attila Juhasz, Heidelberg, 03.01.2010
; 
; Compatibility rules in the structure of problem_params.pro
;
; - The parameter section should begin with a line :
;  '; XXXXXXXXXXXXXX Declaration of the input variables XXXXXXXXXXXXXX' 
;  where X marks an arbitrary character
;
; - The parameter section should end with a line :
;  '; XXXXXXXXXXXXXX End of input variable declaration XXXXXXXXXXXXXX' 
;  where X marks an arbitrary character
;
; - Every line outside the parameter section and lines in the parameter
;   section beginning with semicolon are read as one single string and
;   will not be modified 
;
; - Each line containing a parameter should look like :
;     parameter_name = 1.09d10 * AU ; Optional comment to describe what the variable is
;   Such lines will be split into three parts : parameter name,
;   parameter value and comment. These parts are separated by '=' and
;   ';' signs. If there is no comment (i.e. no ';' sign in the line)
;   the string containing the line is split into two parts, parameter
;   name and its value. 
;
; - The line number of each line is stored to make it possible to
;   re-write the read file exactly in the same way as it was 
;
;  Modifications :
;               
;     11.01.2010.
;     - Addition of the 'extra' recognition. If a parameter name
;       begins with 'extra_' then the line will be stored in a different 
;       part of the final structure. 
;     08.04.2010
;     - By creating the 'extra' tabs the CW_FIELD for
;       'extra_xxxxx_enable' will be non-editable
;************************************************************************************************
function read_params_radmc3d, read_file_name=read_file_name

if not keyword_set(read_file_name) then read_file_name = 'problem_params.pro'
;
; Read the whole file
;
openr, 1, read_file_name

txt = ''
while eof(1) eq 0 do begin
   dum = ''
   readf, 1, dum
   txt = [txt, dum]
endwhile
txt = txt(1:n_elements(txt)-1)
close, 1

;
; Count how many parameters we have and
;   how many blocks they are split into
;
npar = 0 ; Number of parameters
nexpar = 0 ; Number of extra parameters (parameter name begins with 'extra_')
param_sec_iline = [9999,9999] ; Line numbers where the parameter section begins and ends

for i=0, n_elements(txt)-1 do begin
;
; Store the line number where the parameter section begins
;
   if strpos(strtrim(txt(i), 2), 'Declaration of the input variables') ge 0 then $
      param_sec_iline(0) = i
;
; Store the line number where the parameter section ends
;
   if strpos(strtrim(txt(i), 2), 'End of input variable declaration') ge 0 then $
      param_sec_iline(1) = i
;
; If we are in the parameter section then count the number of parameters
;
   if i gt param_sec_iline(0) and i lt param_sec_iline(1) then begin
      if strmid(strtrim(txt(i), 2), 0, 1) ne ';' and strlen(strtrim(txt(i),2)) gt 3 then begin
         dum = strmid(strcompress(txt(i), /remove_all), 0, 6)
         if dum eq 'extra_' or dum eq 'Extra_' or dum eq 'EXTRA_' then begin
            nexpar = nexpar + 1
         endif else begin
            npar = npar+1
         endelse
      endif
   endif

endfor

print, npar, ' Parameters have been found'
print, nexpar, ' Extra parameters have been found'

;
; Create a structure in which the content of the problem_params.pro
; file will be stored
;

nl_cont = n_elements(txt)-npar-nexpar
if nexpar gt 0 then begin
   file_content = {par:  replicate({name    : '',$
                                    value   : '',$
                                    comment : '',$
                                    iline   : 0   }, npar),$
                   ext:  replicate({name    : '',$
                                    value   : '',$
                                    group   : '',$
                                    group_id: 0, $
                                    comment : '',$
                                    iline   : 0   }, nexpar),$
                   cont: replicate({text    : '',$
                                    iline   : 0   }, nl_cont)}
endif else begin
   file_content = {par:  replicate({name    : '',$
                                    value   : '',$
                                    comment : '',$
                                    iline   : 0   }, npar),$
                   ext:  replicate({name    : '',$
                                    value   : '',$
                                    group   : '',$
                                    group_id: -1, $
                                    comment : '',$
                                    iline   : -1   }, 1),$
                   cont: replicate({text    : '',$
                                    iline   : 0   }, nl_cont)}
end
   

;
; Put the content of the read_file_name (txt) into the
;  structure above
;

pid = 0 ; index of file_content.par 
eid = 0 ; index of file_content.ext
cid = 0 ; index of file_content.cont
for i=0, n_elements(txt)-1 do begin

   if i gt param_sec_iline(0) and i lt param_sec_iline(1) then begin
;
; We are in the parameter section
;  Now separate name, value and comments of the input variables
;

;
;  (Check if the line is not commented, i.e. does not begin with ';')
;
      if strmid(strtrim(txt(i), 2), 0, 1) ne ';' and strlen(strtrim(txt(i),2)) gt 3 then begin
         dummy = strsplit(txt(i), /extract)

;
; Now check if the parameter name begins with 'extra_' meaning an
; extra parameter or if not it is a normal mandatory parameter
;
         dum = strmid(strcompress(txt(i), /remove_all), 0, 6)
         if dum eq 'extra_' or dum eq 'Extra_' or dum eq 'EXTRA_' then begin

           
;
; We have an 'extra_' parameter
;

;
; Check if we have a comment in the line (e.g. explaining what the
; parameter means)
;
            ii = strpos(txt(i), ';')
;
; If we have some comment
;
            if ii(0) ge 0 then begin
               dummy = strsplit(txt(i), ';', /extract)
               file_content.ext(eid).comment = dummy(1)
               
               dummy2 = strsplit(dummy(0), '=', /extract)
               file_content.ext(eid).name  = dummy2(0)
               file_content.ext(eid).value = dummy2(1)
               file_content.ext(eid).iline = i
               
               dummy3 = strsplit(dummy2(0), '_', /extract)
               if eid eq 0 then begin
                  file_content.ext(eid).group = dummy3(1)
                  file_content.ext(eid).group_id = 0
               endif else begin
                  file_content.ext(eid).group = dummy3(1)
                  if dummy3(1) eq file_content.ext(eid-1).group then $
                     file_content.ext(eid).group_id = file_content.ext(eid-1).group_id 
                  if dummy3(1) ne file_content.ext(eid-1).group then $
                     file_content.ext(eid).group_id = file_content.ext(eid-1).group_id + 1 
                 
               endelse
               eid = eid+1
            
;
; If we do not have any comment
;
            endif else begin
               dummy = strsplit(txt(i), '=', /extract)
               file_content.ext(eid).name  = dummy(0)
               file_content.ext(eid).value = dummy(1)
               file_content.ext(eid).iline = i

               dummy2 = strsplit(dummy2(0), '_', /extract)
               file_content.ext(eid).group = dummy2(1)
               eid = eid+1
            endelse

            
         endif else begin

;
; We have a 'normal' mandatory parameter
;            

;
; Check if we have a comment in the line (e.g. explaining what the
; parameter means)
;
            ii = strpos(txt(i), ';')
;
; If we have some comment
;
            if ii(0) ge 0 then begin
               dummy = strsplit(txt(i), ';', /extract)
               file_content.par(pid).comment = dummy(1)
               
               dummy2 = strsplit(dummy(0), '=', /extract)
               file_content.par(pid).name  = dummy2(0)
               file_content.par(pid).value = dummy2(1)
               file_content.par(pid).iline = i
               pid = pid+1
            
;
; If we do not have any comment
;
            endif else begin
               dummy = strsplit(txt(i), '=', /extract)
               file_content.par(pid).name  = dummy(0)
               file_content.par(pid).value = dummy(1)
               file_content.par(pid).iline = i
               pid = pid+1
            endelse
         
         endelse
      endif else begin
;
; If we have a commented line
;
         file_content.cont(cid).text  = txt(i)
         file_content.cont(cid).iline = i
         cid = cid+1
         
      endelse
   endif else begin
;
; If we are out of the parameter section
;
      file_content.cont(cid).text  = txt(i)
      file_content.cont(cid).iline = i
      cid = cid+1
   endelse
      
endfor

return, file_content
end
;************************************************************************************************
; Procedure to write out the parameter file of RADMC3D
;    Attila Juhasz, Heidelberg, 03.01.2010
;
; How it works: 
;   
; - The input structure 'file_content' (the output of
;   read_params_radmc3d()) is written into 'write_file_name'
; - Lines which are
;
;
;  Modifications :
;               
;     11.01.2010.
;     - Addition of the 'extra' recognition. If a parameter name
;       begins with 'extra_' then the line will be stored in a different 
;       part of the final structure. 
;
;************************************************************************************************
pro write_params_radmc3d, file_content, write_file_name=write_file_name

;
; If no write_file_name is specified then write out the parameters
; into problem_params.pro
;
if not keyword_set(write_file_name) then write_file_name = 'problem_params.pro'

openw, 1, write_file_name
;
; Go line by line and write out the content of the file and modify
; only the parameter values
;
npar  = n_elements(file_content.par)
ncont = n_elements(file_content.cont)
next  = n_elements(file_content.ext)
if file_content.ext(0).group_id eq -1 then next = 0
for i=0, npar + ncont + next -1  do begin

;
; The whole txt file was stored in the structure file_content line by
; line. Now search for the ith line to write it out.
;
   ii = where(file_content.par.iline eq i)

;
; If the ith line contains a parameter
;
   if ii(0) ge 0 then begin
      if strlen(file_content.par(ii(0)).comment) gt 2 then begin
         printf, 1, strtrim(file_content.par(ii(0)).name, 2),' = ',$
                    strtrim(file_content.par(ii(0)).value, 2),$
                 ' ; ', strtrim(file_content.par(ii(0)).comment, 2), format='(A,A,A,A,A)'
      endif else begin
         printf, 1, strtrim(file_content.par(ii(0)).name, 2),' = ', $
                 strtrim(file_content.par(ii(0)).value, 2), format='(A,A,A)'
      endelse
;
; If the ith line does not contain any parameter only floating text
;
   endif else begin
      jj = where(file_content.ext.iline eq i) 
      if jj(0) ge 0 then begin
         if strlen(file_content.ext(jj(0)).comment) gt 2 then begin
            printf, 1, strtrim(file_content.ext(jj(0)).name, 2),' = ',$
                    strtrim(file_content.ext(jj(0)).value, 2),$
                    ' ; ', strtrim(file_content.ext(jj(0)).comment, 2), format='(A,A,A,A,A)'
         endif else begin
            printf, 1, strtrim(file_content.ext(jj(0)).name, 2),' = ', $
                    strtrim(file_content.ext(jj(0)).value, 2), format='(A,A,A)'
         endelse
      endif else begin
         kk = where(file_content.cont.iline eq i)
         if kk(0) ge 0 then begin
            printf, 1, file_content.cont(kk(0)).text, format='(A)'
         endif else begin
;
; If the line number was not found in the structure (ERROR!)
;
            print, ' ERROR! No line number has been found ', i
            stop
         endelse
      endelse
   endelse
endfor
close, 1
end

