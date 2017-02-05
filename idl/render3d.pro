;============================================================================
;              3-D VOLUME AND ISO-SURFACE RENDERING ROUTINE
;                   BASED ON THE D_VECTRACK.PRO ROUTINE
;                         IN THE IDL DISTRIBUTION
;============================================================================


;----------------------------------------------------------------------------
;                          A DRAW ROUTINE 
;----------------------------------------------------------------------------
pro demo_draw, oWindow, oView, debug=debug
;
;Flush and print any accumulated math errors
;
void = check_math(/print)
;
;Silently accumulate any subsequent math errors, unless we are debugging.
;
orig_except = !except
!except = ([0, 2])[keyword_set(debug)]
;
;Draw.
;
oWindow->Draw, oView
;
;Silently (unless we are debugging) flush any accumulated math errors.
;
void = check_math(PRINT=keyword_set(debug))
;
;Restore original math error behavior.
;
!except = orig_except
end



;----------------------------------------------------------------------------
;     MAKE THE VECTOR FIELD TO BE DISPLAYED IN ONE OF THE 2-D SURFACES
;----------------------------------------------------------------------------
PRO d_vectrackMkVectors, u, v, w, fVerts, iConn, X=x, Y=y, Z=z, SCALE=scale,$
               RANDOM=random, NVECTORS=nvectors, STEPSIZE=stepsize,$
               bMAG=bMag, VC=vc
  ;;
  ;; Ensure volumes match in number of elements.
  ;;
  nSamples = N_ELEMENTS(u)
  IF (nSamples NE N_ELEMENTS(v) OR nSamples NE N_ELEMENTS(w)) THEN BEGIN
     MESSAGE,'Number of elements in u, v, and w must match.'
     RETURN
  ENDIF
  sz = SIZE(u)
  ;;
  ;; Handle keywords, set defaults.
  ;;
  IF (N_ELEMENTS(scale) EQ 0) THEN scale = 1.0
  IF (N_ELEMENTS(random) EQ 0) THEN random=0
  IF (N_ELEMENTS(nvectors) EQ 0) THEN nvectors = 100
  IF (N_ELEMENTS(stepsize) eq 0) then stepsize = 1
  IF (N_ELEMENTS(bMag) EQ 0) THEN bMag = BYTSCL(SQRT(u^2+v^2+w^2))
  ;;
  ;; Get plane information.
  ;;
  doPlane = 0
  IF (N_ELEMENTS(x) GT 0) THEN BEGIN
     IF (random EQ 0) THEN BEGIN
        nRow = sz[2] / stepsize
        nCol = sz[3] / stepsize
     ENDIF
     doPlane = 1
  ENDIF
  IF (N_ELEMENTS(y) GT 0) THEN BEGIN
     IF (doPlane) THEN MESSAGE,'X, Y, and Z keywords are mutually exclusive'
     IF (random EQ 0) THEN BEGIN
        nRow = sz[3] / stepsize
        nCol = sz[1] / stepsize
     ENDIF
     doPlane = 2
  ENDIF
  IF (N_ELEMENTS(z) GT 0) THEN BEGIN
     IF (doPlane) THEN MESSAGE,'X, Y, and Z keywords are mutually exclusive'
     IF (random EQ 0) THEN BEGIN
        nRow = sz[2] / stepsize
        nCol = sz[1] / stepsize
     ENDIF
     doPlane = 3
  ENDIF
  IF (doPlane EQ 0) THEN MESSAGE, 'Must specify a plane.'
  ;;
  ;; Grab max, min values in vector volumes.
  ;;
  maxU = MAX(u, MIN=minU)
  maxV = MAX(v, MIN=minV)
  maxW = MAX(w, MIN=minW)
  ;;
  ;; Compute the magnitude.
  ;;
  mag = SQRT((maxU-minU)^2 + (maxV-minV)^2 + (maxW-minW)^2)
  fNorm = scale / mag
  ;;
  ;; Compute radomly spaced vectors.
  ;;
  IF (random) THEN BEGIN
     fVerts = FLTARR(3, 2*nvectors)
     iConn = LONARR(3*nvectors)
     vc = BYTARR(2*nvectors)
     ;;
     CASE doPlane OF
        1: BEGIN                ; X=x
           randomX = REPLICATE(x, nvectors)
           seed = x
           randomY = RANDOMU(seed, nvectors) * (sz[2]-1)
           randomZ = RANDOMU(seed, nvectors) * (sz[3]-1)
        END
        2: BEGIN                ; Y=y
           seed = y
           randomX = RANDOMU(seed, nvectors) * (sz[1]-1)
           randomY = REPLICATE(y, nvectors)
           randomZ = RANDOMU(seed, nvectors) * (sz[3]-1)
        END
        3: BEGIN                ; Z=z
           seed = z
           randomX = RANDOMU(seed, nvectors) * (sz[1]-1)
           randomY = RANDOMU(seed, nvectors) * (sz[2]-1)
           randomZ = REPLICATE(z, nvectors)
        END
     ENDCASE
     ;;
     inds = LINDGEN(nvectors)
     x0 = randomx[inds]
     y0 = randomy[inds]
     z0 = randomz[inds]
     ;;
     v0 = transpose([[x0],[y0],[z0]])
     v1 = transpose([[x0+u[x0,y0,z0]*fNorm],$
                     [y0+v[x0,y0,z0]*fNorm],$
                     [z0+w[x0,y0,z0]*fNorm]])
     fVerts[*,inds*2] = v0[*,inds]
     fVerts[*,inds*2+1] = v1[*,inds]
     iConn[inds*3] = 2
     iConn[inds*3+1] = inds*2
     iConn[inds*3+2] = inds*2+1
     ;;
     vc[inds*2] = bMag[x0,y0,z0]
     vc[inds*2+1] = bMag[x0,y0,z0]
     ;;
  ENDIF ELSE BEGIN
     ;;
     ;; Compute evenly sampled vectors.
     ;;
     nV = nRow*nCol
     ;;
     fVerts = FLTARR(3, 2*nV)
     iConn = LONARR(3*nV)
     vc = BYTARR(2*nV)
     ;;
     CASE doPlane OF
        1: BEGIN                ; X=x
           x0 = REPLICATE(x,nV)
           y0 = REFORM((REPLICATE(1,nCol) # (LINDGEN(nRow)*stepsize)),nV)
           z0 = REFORM(((LINDGEN(nCol)*stepsize) # REPLICATE(1,nRow)),nV)
        END
        2: BEGIN                ; Y=y
           y0 = REPLICATE(y,nRow*nCol)
           z0 = REFORM((REPLICATE(1,nCol) # (LINDGEN(nRow)*stepsize)),nV)
           x0 = REFORM(((LINDGEN(nCol)*stepsize) # REPLICATE(1,nRow)),nV)
        END
        3: BEGIN                ; Z=z
           z0 = REPLICATE(z,nRow*nCol)
           y0 = REFORM((REPLICATE(1,nCol) # (LINDGEN(nRow)*stepsize)),nV)
           x0 = REFORM(((LINDGEN(nCol)*stepsize) # REPLICATE(1,nRow)),nV)
        END
     ENDCASE
     ;;
     inds = LINDGEN(nV)
     v0 = transpose([[x0],[y0],[z0]])
     v1 = transpose([[x0+u[x0,y0,z0]*fNorm],$
                     [y0+v[x0,y0,z0]*fNorm],$
                     [z0+w[x0,y0,z0]*fNorm]])
     fVerts[*,inds*2] = v0[*,inds]
     fVerts[*,inds*2+1] = v1[*,inds]
     iConn[inds*3] = 2
     iConn[inds*3+1] = inds*2
     iConn[inds*3+2] = inds*2+1
     ;;
     vc[inds*2] = bMag[x0,y0,z0]
     vc[inds*2+1] = bMag[x0,y0,z0]
  ENDELSE
;;
END




;----------------------------------------------------------------------------
;                INTEGRATE V TO GET STREAM RIBBONS
;                       (STILL ACTIVE??)
;----------------------------------------------------------------------------
PRO d_vectrackRibbonTrace,vdata,start,auxdata,STEPS=steps,FRAC=frac, $
                          WIDTH=width, UP=up,COLOR=color,  $
                          OUTVERTS = outverts, OUTCONN = outconn,  $
                          VERT_COLORS = vertcolors
  ;;
  ;;  Set defaults
  ;;
  if (N_ELEMENTS(steps) eq 0) then steps = 100
  if (N_ELEMENTS(frac) eq 0) then frac = 1.0
  if (N_ELEMENTS(up) eq 0) then up = [0.0,0.0,1.0]
  if (N_ELEMENTS(width) eq 0) then width = .5
  ;;
  ;;  Trace the vector field using IDL's particle trace routine
  ;;
  PARTICLE_TRACE,vdata,start,overts,oconn,onormals, $
        MAX_ITERATIONS=steps, MAX_STEPSIZE=frac,INTEGRATION=0 $
        ,ANISOTROPY=[1,1,1], SEED_NORMALS=up
  ;;
  ;;  Now make the ribbon around the streamline and add color
  ;;
    if((N_ELEMENTS(oconn) gt 0) and (SIZE(overts, /N_DIMENSIONS) eq 2))  $
        then begin
        STREAMLINE,overts,oconn,onormals*width,outverts,outconn
        cdata = INTERPOLATE(auxdata,outverts[0,*],outverts[1,*],outverts[2,*])
        cdata = REFORM(cdata,N_ELEMENTS(outverts)/3)
        vertcolors = BYTSCL(cdata)
    end
END




;----------------------------------------------------------------------------
;                 INTEGRATE V TO GET STREAMLINES
;                       (STILL ACTIVE??)
;----------------------------------------------------------------------------
PRO d_vectrackStreamlineTrace,vdata,start,auxdata,STEPS=steps,FRAC=frac, $
                              OUTVERTS = outverts, OUTCONN = outconn,  $
                              VERT_COLORS = vertcolors
  ;;
  ;;  How many integration steps (???)
  ;;
  if (N_ELEMENTS(steps) eq 0) then steps = 100
  if (N_ELEMENTS(frac) eq 0) then frac = 1.0
  ;;
  ;;  Use IDL's particle_trace routine to integrate the vector field
  ;;
  PARTICLE_TRACE,vdata,start,outverts,outconn, $
        MAX_ITERATIONS=steps, MAX_STEPSIZE=frac,INTEGRATION=0 $
        ,ANISOTROPY=[1,1,1]
  ;;
  ;;  Put colors on the streamlines
  ;;
  if((N_ELEMENTS(outconn) gt 0) and $
               (SIZE(outverts, /N_DIMENSIONS) eq 2)) then begin
     cdata = INTERPOLATE(auxdata,outverts[0,*],outverts[1,*],outverts[2,*])
     cdata = REFORM(cdata,N_ELEMENTS(outverts)/3)
     vertcolors = BYTSCL(cdata)
  end
END





;----------------------------------------------------------------------------
;                      UPDATE THE ISO SURFACES
;----------------------------------------------------------------------------
PRO d_vectrackIsoSurfUpdate,sState,bUpdate
  ;;
  ;;  Check if we must update
  ;;
  IF (sState.bIsoShow AND ((bUpdate EQ 5) OR (bUpdate EQ 4))) THEN BEGIN
     ;;
     ;;  Get the data
     ;;
     sState.oVols[sState.iIsoVol]->GetProperty,DATA0=vdata,$
        /NO_COPY
     ;;
     ;;  Create the 3-D iso surface in the form of a list of vertices 
     ;;  and polygons
     ;;
     SHADE_VOLUME,vdata,sState.fIsoLevel,fVerts,iConn
     ;;
     ;;  Check if any vertex is created
     ;;
     IF (N_ELEMENTS(fVerts) LE 0) THEN BEGIN
        fVerts=[0,0,0]
        iConn=0
     END
     ;;
     ;;  Pass this list of vertices and polygons to the polygon object
     ;;
     if n_elements(fVerts) gt 3 then begin
        sState.oIsoPolygon->SetProperty,DATA=fVerts,POLYGONS=iConn
     endif else begin
        sState.oIsoPolygon->SetProperty,HIDE=1
     endelse
     ;;
     ;;  Give data back to the volume object
     ;;
     sState.oVols[sState.iIsoVol]->SetProperty,DATA0=vdata,$
        /NO_COPY
  ENDIF
END





;----------------------------------------------------------------------------
;                    UPDATE THE CURRENT 2-D PLANES
;----------------------------------------------------------------------------
PRO d_vectrackPlanesUpdate,sState,bUpdate
  ;;
  ;;  Get some of the information of the data that should be 
  ;;  displayed in the planes
  ;;
  if sState.incl_velo ne 0 then begin
     WIDGET_CONTROL, sState.wEvenField, GET_VALUE=stepsize
     WIDGET_CONTROL, sState.wRandomField, GET_VALUE=nvectors
     WIDGET_CONTROL, sState.wScaleField, GET_VALUE=scale
  endif else begin
     stepsize=1
     nvectors=2
     scale=2
  endelse
  ;;
  ;;  Do 2-d x-plane
  ;; 
  ;;  ...Check which kind of update to do
  ;;
  IF (sState.bShow[0] AND ((bUpdate EQ 1) OR (bUpdate EQ 4))) THEN BEGIN
     ;;
     ;;  Get position of x-slider
     ;;
     WIDGET_CONTROL, sState.wXSlider, GET_VALUE=x
     ;;
     ;;  Check which thing to put onto plane (vector or image)
     ;;
     IF (sState.bImage[0]) THEN BEGIN
        ;;
        ;;  ...Image
        ;;
        sState.oVols[sState.iImgVol]->GetProperty,DATA0=vdata,$
           /NO_COPY,RGB_TABLE0=ctab
        d_vectrackSampleOrthoPlane,vdata,0,x,sState.oSlices[0], $
                                   sState.oImages[0],ctab,sState.iAlpha
        sState.oVols[sState.iImgVol]->SetProperty,DATA0=vdata,$
           /NO_COPY
     END ELSE BEGIN
        ;;
        ;;  ...Vector field
        ;;
        sState.oVols[2]->GetProperty,DATA0=bMag,RGB_TABLE0=pal,/NO_COPY
        d_vectrackMkVectors,sState.u,sState.v,sState.w, $
               fVerts,iConn,X=x,SCALE=scale, $
               RANDOM=sState.bRandom, NVECTORS=nvectors,$
               STEPSIZE=stepsize,BMAG=bMag,Vc=vc
        sState.oVols[2]->SetProperty,DATA0=bMag,/NO_COPY
     END
  ENDIF
  ;;
  ;;  Do 2-d y-plane
  ;; 
  ;;  ...Check which kind of update to do
  ;;
  IF (sState.bShow[1] AND ((bUpdate EQ 2) OR (bUpdate EQ 4))) THEN BEGIN
     ;;
     ;;  Get position of y-slider
     ;;
     WIDGET_CONTROL, sState.wYSlider, GET_VALUE=y
     ;;
     ;;  Check which thing to put onto plane (vector or image)
     ;;
     IF (sState.bImage[1]) THEN BEGIN
        ;;
        ;;  ...Image
        ;;
        sState.oVols[sState.iImgVol]->GetProperty,DATA0=vdata,$
           /NO_COPY,RGB_TABLE0=ctab
        d_vectrackSampleOrthoPlane,vdata,1,y,sState.oSlices[1], $
                   sState.oImages[1],ctab,sState.iAlpha
        sState.oVols[sState.iImgVol]->SetProperty,DATA0=vdata,$
           /NO_COPY
     END ELSE BEGIN
        ;;
        ;;  ...Vector field
        ;;
        sState.oVols[2]->GetProperty,DATA0=bMag,RGB_TABLE0=pal,/NO_COPY
        d_vectrackMkVectors,sState.u,sState.v,sState.w, $
                 fVerts,iConn,Y=y,SCALE=scale, $
                 RANDOM=sState.bRandom, NVECTORS=nvectors,$
                 STEPSIZE=stepsize,BMAG=bMag,Vc=vc
        sState.oVols[2]->SetProperty,DATA0=bMag,/NO_COPY
     END
  ENDIF
  ;;
  ;;  Do 2-d z-plane
  ;; 
  ;;  ...Check which kind of update to do
  ;;
  IF (sState.bShow[2] AND ((bUpdate EQ 3) OR (bUpdate EQ 4))) THEN BEGIN
     ;;
     ;;  Get position of z-slider
     ;;
     WIDGET_CONTROL, sState.wZSlider, GET_VALUE=z
     ;;
     ;;  Check which thing to put onto plane (vector or image)
     ;;
     IF (sState.bImage[2]) THEN BEGIN
        ;;
        ;;  ...Image
        ;;
        sState.oVols[sState.iImgVol]->GetProperty,DATA0=vdata,$
           /NO_COPY,RGB_TABLE0=ctab
        d_vectrackSampleOrthoPlane,vdata,2,z,sState.oSlices[2], $
                                   sState.oImages[2],ctab,sState.iAlpha
        sState.oVols[sState.iImgVol]->SetProperty,DATA0=vdata,$
           /NO_COPY
     END ELSE BEGIN
        ;;
        ;;  ...Vector field
        ;;
        sState.oVols[2]->GetProperty,DATA0=bMag,RGB_TABLE0=pal,/NO_COPY
        d_vectrackMkVectors,sState.u,sState.v,sState.w, $
                 fVerts,iConn,Z=z,SCALE=scale,  $
                 RANDOM=sState.bRandom, NVECTORS=nvectors,$
                 STEPSIZE=stepsize,BMAG=bMag,Vc=vc
        sState.oVols[2]->SetProperty,DATA0=bMag,/NO_COPY
     END
  ENDIF
END





;----------------------------------------------------------------------------
;       ROUTINE TO DO A TIME STEP: CHANGE THE DATA TO NEW TIME SLICE
;              NOTE: ONLY IF MULTIPLE TIME SLICES ARE SPECIFIED
;----------------------------------------------------------------------------
FUNCTION d_vectrackChangeData,sState,iTimeStep
  ;;
  ;;  Check if more than one time step are available
  ;;
  IF (PTR_VALID(sState.pSteps)) THEN BEGIN
     ;;
     ;; Is there anything to do?
     ;;
     IF (iTimeStep GE sState.nSteps) THEN RETURN,0
     IF (iTimeStep EQ sState.curStep) THEN RETURN,0
     ;;
     ;; Update UVW to the new time
     ;;
     sState.u = *((*sState.pSteps).fU[iTimeStep])
     sState.v = *((*sState.pSteps).fV[iTimeStep])
     sState.w = *((*sState.pSteps).fW[iTimeStep])
     ;;
     ;; Update volume data
     ;;
     sState.oVols[0]->SetProperty,$
            DATA0=*((*sState.pSteps).bP[iTimeStep])
     sState.oVols[1]->SetProperty,$
            DATA0=*((*sState.pSteps).bT[iTimeStep])
     sState.oVols[2]->SetProperty,$
            DATA0=*((*sState.pSteps).bM[iTimeStep])
     ;;
     ;; And the current timestep
     ;;
     sState.curStep = iTimeStep
     ;;
     ;; Recompute streamlines [STILL ACTIVE??]
     ;;
     oList = sState.oStreamModel->Get(/ALL)
     sz=size(oList)
     IF (sz[2] EQ 11) THEN BEGIN
        FOR i=0,N_ELEMENTS(oList)-1 DO BEGIN
           oList[i]->GetProperty,UVALUE=xyz
           d_vectrackBuildRibbon,sState,oLIst[i],xyz
        END
     END
     ;;
     ;; Recompute the iso-surface and planes
     ;;
     d_vectrackIsoSurfUpdate,sState,5
     d_vectrackPlanesUpdate,sState,4
     ;;
     ;; Done... Tell that new time slice was installed
     ;;
     RETURN,1
  END
  ;;
  ;; Nothing happened
  ;;
  RETURN,0
END




;----------------------------------------------------------------------------
;       UPDATE THE COLOR BAR TO MATCH THE CURRENT VREND SELECTION
;----------------------------------------------------------------------------
PRO d_vectrackColorBarUpdate,sState
  ;;
  ;;  Get data of the color bar
  ;;
  oList = sState.oCBTop->Get(/ALL)
  ;;
  ;;  Find the min and max range
  ;;
  fMin = sState.fScale[0,sState.iVrendVol]
  fMax = sState.fScale[1,sState.iVrendVol]
  ;;
  ;;  Set the position and orientation of the color bar in 3-D space
  ;;
  sState.oCBTop->SetProperty,TRANSFORM=IDENTITY(4)
  dx = fMax-fMin
  sState.oCBTop->Translate,-(fMax+fMin)*0.5,0.0,0.0
  sState.oCBTop->Scale,1.0/dx,1.0/260.0,1.0
  sState.oCBTop->Translate,0.0,-0.725,0.0
  ;;
  ;;  The Image is (0)
  ;;
  sState.oVols[sState.iVrendVol]->GetProperty,RGB_TABLE0=pal
  pal = TRANSPOSE(pal)
  rgb = REFORM(pal[*,INDGEN(256*16) MOD 256],3,256,16)
  oList[0]->SetProperty,DATA=rgb,DIMENSIONS=[dx,16],LOCATION=[fMin,0.]
  ;;
  ;;  The Axis is (1)
  ;;
  sTitle = sState.datanames[sState.iVrendVol]
  ;;
  ;;  Now install the stuff
  ;;
  oList[1]->GetProperty,TICKTEXT=oText,TITLE=oTitle
  oText->SetProperty,CHAR_DIMENSIONS=[0,0]
  oTitle->SetProperty,CHAR_DIMENSIONS=[0,0],STRING=sTitle
  oList[1]->SetProperty,RANGE=[fMin,fMax]
  ;;
END




;----------------------------------------------------------------------------
;    ROUTINE TO READ AND SETUP THE VOLUME OBJECT PALETTES AND OPACITIES
;----------------------------------------------------------------------------
PRO d_vectrackReadVoxelPalettes,vobj,coltabnr=coltabnr,opacmax=opacmax
  COMMON COLORS, R_orig, G_orig, B_orig, R_curr, G_curr, B_curr  
  ;;
  ;;  Default
  ;;
  if n_elements(coltabnr) eq 0 then coltabnr=-3
  if n_elements(opacmax) eq 0 then opacmax=127
  ;;
  ;;  Declare the color table and opacity table arrays
  ;;
  vColors = BYTARR(256,3,/NOZERO)
  vOpac = BYTARR(256,/NOZERO)
  ;;
  ;;  Read the color table from a file
  ;;    [NEEDS TO BE MODIFIED]
  ;;
  ;q=findfile('colortable.idl',count=count)
  q=file_search('colortable.idl',count=count)
  if count ge 1 and coltabnr eq -1 then begin
     ;;
     ;; Read color table from file
     ;;
     print,'Reading color table from file colortable.idl'
     OPENR, lun, /GET_LUN, 'colortable.idl'
     READU, lun,  vColors
     CLOSE, lun
     FREE_LUN, lun
  endif else begin
     if coltabnr ge 0 then begin
        ;;
        ;; Use the loadct to get a standard color table
        ;;
        loadct,coltabnr
        vcolors[*,0] = R_curr
        vcolors[*,1] = G_curr
        vcolors[*,2] = B_curr
     endif else begin
        ;;
        ;; Create a color table similar to the original d_vectrack.pro 
        ;;
        x=dindgen(256)
        if coltabnr ge -2 then begin
           width = 30.
           minc  = 127
        endif else begin
           width = 45.
           minc  = 50
        endelse
        xr = 128.+1.2*width
        xg = 128.
        xb = 128.-1.2*width
        vColors[*,0] = bytscl(1./(1./((400*exp(-((x-xr)/width)^2)^0.9))^4+(1./256.)^4)^0.25)
        vColors[floor(xr):*,0] = vColors[floor(xr):*,0] > minc
        vColors[*,1] = bytscl(1./(1./((400*exp(-((x-xg)/width)^2)^0.9))^4+(1./256.)^4)^0.25)
        vColors[*,2] = bytscl(1./(1./((400*exp(-((x-xb)/width)^2)^0.9))^4+(1./256.)^4)^0.25)
        vColors[0:floor(xb),2] = vColors[0:floor(xb),2] > minc
     endelse
  endelse
  ;;
  ;;  Read the opacity table from a file 
  ;;
  ;q=findfile('opacty.idl',count=count)
  q=file_search('opacty.idl',count=count)
  if count ge 1 then begin
     ;;
     ;; Read opacity table from file
     ;;
     print,'Reading opacity table from file opacty.idl'
     OPENR, lun, /GET_LUN, 'opacity.idl'
     READU, lun,  vOpac
     CLOSE, lun
     FREE_LUN, lun
  endif else begin
     ;;
     ;; Create opacity file
     ;;
     vOpac[*] = bytscl(dindgen(256)) < opacmax
  endelse
  ;;
  ;;  Install colors
  ;;
  vobj->SetProperty,RGB_TABLE0=vColors,OPACITY_TABLE0=vOpac
  ;;
END





;----------------------------------------------------------------------------
;                UPDATE TEXTURE MAPPED IMAGE DISPLAY OF SLICES
;----------------------------------------------------------------------------
PRO d_vectrackSampleOrthoPlane,data,axis,slice,oPoly,oImage,ctab,alpha
    sz=size(data)-1
    ;;
    ;;  Check for which of the 2-D planes to do this
    ;;
    CASE axis OF
       0: BEGIN
          img = data[slice,*,*]
          img = REFORM(img,sz[2]+1,sz[3]+1,/OVERWRITE)
          fTxCoords = [[0,0],[1,0],[1,1],[0,1]]
          verts=[[slice,0,0],[slice,sz[2],0], $
                 [slice,sz[2],sz[3]],[slice,0,sz[3]]]
       END
       1: BEGIN
          img = data[*,slice,*]
          img = REFORM(img,sz[1]+1,sz[3]+1,/OVERWRITE)
          fTxCoords = [[0,0],[1,0],[1,1],[0,1]]
          verts=[[0,slice,0],[sz[1],slice,0], $
                 [sz[1],slice,sz[3]],[0,slice,sz[3]]]
       END
       2: BEGIN
          img = data[*,*,slice]
          img = REFORM(img,sz[1]+1,sz[2]+1,/OVERWRITE)
          fTxCoords = [[0,0],[1,0],[1,1],[0,1]]
          verts=[[0,0,slice],[sz[1],0,slice], $
                 [sz[1],sz[2],slice],[0,sz[2],slice]]
       END
    END
    ;;
    ;; Convert to 3xNxM or 4xNxM, i.e. convert to colorfull image
    ;;
    sz=size(img)
    IF ((alpha[0] EQ 0) AND (alpha[1] EQ 255)) THEN BEGIN
        rgbtab=TRANSPOSE(ctab)
        rgb=rgbtab[*,img]
        rgb=REFORM(rgb,3,sz[1],sz[2],/OVERWRITE)
    END ELSE BEGIN
        rgbtab=bytarr(4,256)
        rgbtab[0:2,*]=TRANSPOSE(ctab)
        rgbtab[3,*] = 0
        rgbtab[3,alpha[0]:alpha[1]] = 255
        rgb=rgbtab[*,img]
        rgb=REFORM(rgb,4,sz[1],sz[2],/OVERWRITE)
    END
    ;;
    ;;  Install the new image into the image object
    ;;  (which is linked to the quadratic polygon)
    ;;
    oImage->SetProperty,DATA=rgb
    ;;
    ;;  Install the new coordinates of the quadratic textured polygon
    ;;
    oPoly->SetProperty,DATA=verts,TEXTURE_COORD=fTxCoords
    ;;
END






;----------------------------------------------------------------------------
;               CONVERT A MOUSE POINT INTO A STREAMLINE 
;                          [INACTIVE???]
;----------------------------------------------------------------------------
FUNCTION d_vectrackDoStream,sEvent,sState,new_flag
  pick = sState.oWindow->PickData(sState.oView,$
                                  sState.oVols[0], $
                                  [sEvent.x,sEvent.y],dataxyz)
  IF (pick NE 1) THEN RETURN,0

  sState.oVols[2]->GetProperty,xcoord_conv=xc,ycoord_conv=yc,$
                                           zcoord_conv=zc

  IF (new_flag) THEN BEGIN
     WIDGET_CONTROL, sState.wRibbons, GET_VALUE = bRibbons
     IF (bRibbons) THEN BEGIN
        sState.oStreamline = OBJ_NEW('IDLgrPolygon', $
                   SHADING=1,STYLE=2,UVALUE=dataxyz,$
                   xcoord_conv=xc,ycoord_conv=yc,zcoord_conv=zc)
     END ELSE BEGIN
        sState.oStreamline = OBJ_NEW('IDLgrPolyline', $
                UVALUE=dataxyz,$
                xcoord_conv=xc,ycoord_conv=yc,zcoord_conv=zc)
     END
     sState.oStreamModel->Add,sState.oStreamline
  END
  
  d_vectrackBuildRibbon,sState,sState.oStreamline,dataxyz
  
  RETURN,1
END




;--------------------------------------------------------------------------
;                COMPUTE A SINGLE RIBBON / STREAMLINE
;                         [INACTIVE???]
;--------------------------------------------------------------------------
PRO d_vectrackBuildRibbon,sState,oObj,dataxyz

    sState.oVols[sState.iVrendVol]->GetProperty,DATA0=auxdata, $
        RGB_TABLE0=auxpal, /NO_COPY
    auxpal = TRANSPOSE(auxpal)
    grAuxpal = OBJ_NEW('IDLgrPalette',  $
                       auxpal[0, *], auxpal[1, *], auxpal[2, *])
    iStep = 100
    ;fFrac = 1.0/88.0  ; 1.0/max velocity
    fFrac = .5 ;step size

    WIDGET_CONTROL, sState.wRibbons, GET_VALUE = bRibbons
    IF (bRibbons) THEN BEGIN
        d_vectrackRibbonTrace, sState.vdata, dataxyz, auxdata,FRAC=fFrac, $
            STEP=iStep,VERT_COLORS=vertcolors, WIDTH = sState.fWidth, $
            OUTVERTS = outverts, OUTCONN = outconn
    END ELSE BEGIN
        d_vectrackStreamlineTrace, sState.vdata, dataxyz, auxdata,FRAC=fFrac, $
            STEP=iStep,VERT_COLORS=vertcolors, $
            OUTVERTS = outverts, OUTCONN = outconn
    END

    IF ((N_ELEMENTS(outconn) GT 1) AND (SIZE(outverts, /N_DIMENSIONS) EQ 2)) $
        THEN BEGIN
        WIDGET_CONTROL, sState.wRibbons, GET_VALUE = bRibbons
        IF (bRibbons) and OBJ_VALID(oObj) THEN BEGIN 
            ss = N_ELEMENTS(color)*0.5
            oObj->SetProperty, REJECT=0, STYLE=2,SHADING=1,PALETTE=grAuxpal, $
                DATA = outverts, POLYGONS = outconn, VERT_COLORS = vertcolors
        END ELSE BEGIN
            oObj->SetProperty,HIDE=0,DATA=outverts,POLYLINES = outconn, $
                VERT_COLORS=vertcolors, PALETTE = grAuxpal
        END
    END ELSE BEGIN
        oObj->SetProperty,HIDE=1
    END

    sState.oVols[sState.iVrendVol]->SetProperty,DATA0=auxdata,/NO_COPY

END






;-----------------------------------------------------------------
;               CLEANUP EVERYTHING AND DESTROY OBJECTS
;-----------------------------------------------------------------
pro d_vectrackCleanup, wTopBase
  ;;
  ;;  Get the data
  ;;
  WIDGET_CONTROL, wTopBase, GET_UVALUE=sState, /NO_COPY
  ;;
  ;; Destroy the objects.
  ;;
  OBJ_DESTROY, sState.oHolder
  ;;
  ;; Destroy the time step data (if present)
  ;;
  IF (PTR_VALID(sState.pSteps)) THEN BEGIN
     PTR_FREE,(*sState.pSteps).fU
     PTR_FREE,(*sState.pSteps).fV
     PTR_FREE,(*sState.pSteps).fW
     
     PTR_FREE,(*sState.pSteps).bP
     PTR_FREE,(*sState.pSteps).bT
     PTR_FREE,(*sState.pSteps).bM

     PTR_FREE,sState.pSteps
  END
  ;;
  ;;  Restore the color table that was there before this program
  ;;
  TVLCT, sState.colorTable
  ;;
  ;;  If this program is part of a bigger program....
  ;;
  if WIDGET_INFO(sState.groupBase, /VALID_ID) then $
     WIDGET_CONTROL, sState.groupBase, /MAP
  ;;
end   ;  of d_vectrackCleanup






;----------------------------------------------------------------------------
;                    MAIN EVENT HANDLER
;----------------------------------------------------------------------------
PRO d_vectrackEvent, sEvent
  ;;
  ;;  Get the data
  ;;
  WIDGET_CONTROL, sEvent.top, GET_UVALUE=sState, /NO_COPY
  WIDGET_CONTROL, sEvent.top, SET_UVALUE=sState, /NO_COPY
  ;;
  ;;  If event is a kill request...
  ;;
  if (TAG_NAMES(sEvent, /STRUCTURE_NAME) EQ  $
      'WIDGET_KILL_REQUEST') then begin
     WIDGET_CONTROL, sEvent.top, /DESTROY
     RETURN
  endif
  ;;
  ;;  ?? Get the event again ??
  ;;
  WIDGET_CONTROL, sEvent.id, GET_UVALUE=uval
  ;;
  ;;  By default, no updating is needed
  ;;
  bUpdate = 0
  bRedraw = 0
  graphout = 0
  ;;
  ;;  ?? Get the data again ??
  ;;
  WIDGET_CONTROL, sEvent.top, GET_UVALUE=sState, /NO_COPY
  ;;
  ;;
  ;;  ............................................................
  ;;                   Handle the various events
  ;;  ............................................................
  ;;
  ;;if uval ne 'POSTSCRIPT' then sState.volrend=0
  CASE uval OF
     ;;
     ;;  If we have data at multiple times, then install a new time slice
     ;;
     'RIBBON':
     'TIME_STEP': BEGIN
        IF (d_vectrackChangeData(sState,sEvent.value) GT 0) THEN BEGIN
           bRedraw = 1
           bUpdate = 4
       END
     END
     ;;
     ;;  The selector for which quantity to use for volume rendering
     ;; 
     'VREND_VOLSEL': BEGIN
        sState.iVrendVol = sEvent.index
        d_vectrackColorBarUpdate,sState
        bRedraw = 1
     END
     ;;
     ;;  Button to tell the code to do the volume rendering  
     ;;
;     'VREND': BEGIN
;        sState.oVols[sState.iVrendVol]->SetProperty,HIDE=0
;        WIDGET_CONTROL,sEvent.top,/HOURGLASS
;        demo_draw, sState.oWindow, sState.oView, debug=sState.debug
;        sState.oVols[sState.iVrendVol]->SetProperty,HIDE=1
;        sState.volrend=1
;     END
     ;;
     ;;  [INACTIVE?]  Clear the streams
     ;;
     'CLEAR_STREAMS' : BEGIN
        oList = sState.oStreamModel->Get(/ALL)
        sz=size(oList)
        IF (sz[2] EQ 11) THEN BEGIN
           OBJ_DESTROY,oList
           bRedraw = 1
        END
     END
     ;;
     ;;  Select which quantity to use for the 2-d plane images
     ;;
     'IMG_VOLSEL': BEGIN
        sState.iImgVol = sEvent.index
        bUpdate = 4
        bRedraw = 1
     END
     ;;
     ;;  The widgets for setting the transparency of the 2-d surfaces
     ;;
     'ALPHA_LEVEL': BEGIN
        WIDGET_CONTROL, sState.wAlpha[0], GET_VALUE=v1
        WIDGET_CONTROL, sState.wAlpha[1], GET_VALUE=v2
        IF (v1 GE v2) THEN BEGIN
           sState.iAlpha = [v2,v1]
        END ELSE BEGIN
           sState.iAlpha = [v1,v2]
        END
        bUpdate = 4
        bRedraw = 1
     END
     ;;
     ;;  The widget for setting the maximum opacity of 3-D volume rendering
     ;;
     'OPACITY_MAX': BEGIN
        WIDGET_CONTROL, sState.wOpacityMax, GET_VALUE=v1
        sstate.opacmax = v1
        if sState.nr_data ge 1 then d_vectrackReadVoxelPalettes,sState.oVols[0],$
                  coltabnr=sState.loadct,opacmax=sstate.opacmax
        if sState.nr_data ge 2 then d_vectrackReadVoxelPalettes,sState.oVols[1],$
                  coltabnr=sState.loadct,opacmax=sstate.opacmax
        if sState.nr_data ge 3 then d_vectrackReadVoxelPalettes,sState.oVols[2],$
                  coltabnr=sState.loadct,opacmax=sstate.opacmax
        if sState.incl_velo ne 0 then d_vectrackReadVoxelPalettes,$
                  sState.oVols[sState.nr_data_tot-1],$
                  coltabnr=sState.loadct,opacmax=sstate.opacmax
        bRedraw = 1
     END
     ;;
     ;;  The widget for setting the color table run-time
     ;;
     'COLOR_TABLE': BEGIN
        WIDGET_CONTROL, sState.wColorTabSel, GET_VALUE=v1
        sstate.loadct=v1
        if sState.nr_data ge 1 then d_vectrackReadVoxelPalettes,sState.oVols[0],$
                  coltabnr=sState.loadct,opacmax=sstate.opacmax
        if sState.nr_data ge 2 then d_vectrackReadVoxelPalettes,sState.oVols[1],$
                  coltabnr=sState.loadct,opacmax=sstate.opacmax
        if sState.nr_data ge 3 then d_vectrackReadVoxelPalettes,sState.oVols[2],$
                  coltabnr=sState.loadct,opacmax=sstate.opacmax
        if sState.incl_velo ne 0 then d_vectrackReadVoxelPalettes,$
                  sState.oVols[sState.nr_data_tot-1],$
                  coltabnr=sState.loadct,opacmax=sstate.opacmax
        d_vectrackColorBarUpdate,sState
        bUpdate = 4
        bRedraw = 1
     END
     ;;
     ;;  Select which quantity to use for the iso-surface rendering
     ;;
     'ISO_VOLSEL': BEGIN
        sState.iIsoVol = sEvent.index
        Val  = sState.fIsoLevel
        eps  = (val/256.)
        fVal = (1.-eps)*sState.fScale[0,sState.iIsoVol] + $
               eps*sState.fScale[1,sState.iIsoVol]
        WIDGET_CONTROL, sState.wIsotext, SET_VALUE=$
               STRING(fVal,format=$
               sState.dformat[sState.iIsoVol])
        IF sState.bIsoShow or sState.bVolShow THEN BEGIN
           bUpdate = 5
           bRedraw = 1
        END
     END
     ;;
     ;;  Select which level to use for the iso-surface
     ;;
     'ISO_LEVEL': BEGIN
        WIDGET_CONTROL, sState.wIsoLevel, GET_VALUE=Val
        eps  = (val/256.)
        fVal = (1.-eps)*sState.fScale[0,sState.iIsoVol] + $
               eps*sState.fScale[1,sState.iIsoVol]
        WIDGET_CONTROL, sState.wIsotext, SET_VALUE=$
               STRING(fVal,format=$
               sState.dformat[sState.iIsoVol])
        sState.fIsoLevel = Val
        IF (sState.bIsoShow) THEN BEGIN
           bUpdate = 5
           bRedraw = 1
        END
     END
     ;;
     ;;  Activate/deactivate the iso-surface rendering
     ;;
     'ISO_SHOW': BEGIN
        sState.bIsoShow = 1-sState.bIsoShow
        IF (sState.bIsoShow) THEN BEGIN
           bUpdate = 5
           bRedraw = 1
        END ELSE BEGIN
           sState.oIsopolygon->SetProperty,HIDE=1
           bRedraw = 1
        END
     END
     ;;
     ;;  Activate/deactivate the volume rendering
     ;;
     'VOL_SHOW': BEGIN
        sState.bVolShow = 1-sState.bVolShow
        IF (sState.bVolShow) THEN BEGIN
           bUpdate = 5
           bRedraw = 1
           sState.oVols[sState.iVrendVol]->SetProperty,HIDE=0
        END ELSE BEGIN
           sState.oVols[sState.iVrendVol]->SetProperty,HIDE=1
           bRedraw = 1
        END
     END
     ;;
     ;;  Select what to do in the 2D x-plane (vector, image, nothing)
     ;;
     'X_STYLE': BEGIN
        CASE sEvent.index OF
           2: BEGIN             ; Vector.
              sState.bShow[0] = 1
              sState.bImage[0] = 0
              bRedraw = 1
              bUpdate = 1
           END
           1: BEGIN             ; Image.
              sState.bShow[0] = 1
              sState.bImage[0] = 1
              bRedraw = 1
              bUpdate = 1
           END
           0: BEGIN             ; Hide.
              sState.bShow[0] = 0
              sState.bImage[0] = 0
                        ;;;; sState.oXPolyline->SetProperty, HIDE=1
              sState.oSlices[0]->SetProperty, HIDE=1
              bRedraw = 1
           END
        ENDCASE
     END
     ;;
     ;;  Select what to do in the 2D y-plane (vector, image, nothing)
     ;;
     'Y_STYLE': BEGIN
        CASE sEvent.index OF
           2: BEGIN             ; Vector.
              sState.bShow[1] = 1
              sState.bImage[1] = 0
              bRedraw = 1
              bUpdate = 2
           END
           1: BEGIN             ; Image.
              sState.bShow[1] = 1
              sState.bImage[1] = 1
              bRedraw = 1
              bUpdate = 2
           END
           0: BEGIN             ; Hide.
              sState.bShow[1] = 0
              sState.bImage[1] = 0
                   ;;;;     sState.oYPolyline->SetProperty, HIDE=1
              sState.oSlices[1]->SetProperty, HIDE=1
              bRedraw = 1
           END
        ENDCASE
     END
     ;;
     ;;  Select what to do in the 2D z-plane (vector, image, nothing)
     ;;
     'Z_STYLE': BEGIN
        CASE sEvent.index OF
           2: BEGIN             ; Vector.
              sState.bShow[2] = 1
              sState.bImage[2] = 0
              bRedraw = 1
              bUpdate = 3
           END
           1: BEGIN             ; Image.
              sState.bShow[2] = 1
              sState.bImage[2] = 1
              bRedraw = 1
              bUpdate = 3
           END
           0: BEGIN             ; Hide.
              sState.bShow[2] = 0
              sState.bImage[2] = 0
                     ;;;;      sState.oZPolyline->SetProperty, HIDE=1
              sState.oSlices[2]->SetProperty, HIDE=1
              bRedraw = 1
           END
        ENDCASE
     END
     ;;
     ;;  The value of the slider for the position of the 2-d x-plane
     ;;
     'X_VALUE': BEGIN
        IF (sState.bShow[0]) THEN BEGIN
           bUpdate = 1
           bRedraw = 1
        ENDIF
     END
     ;;
     ;;  The value of the slider for the position of the 2-d y-plane
     ;;
     'Y_VALUE': BEGIN
        IF (sState.bShow[1]) THEN BEGIN
           bUpdate = 2
           bRedraw = 1
        ENDIF
     END
     ;;
     ;;  The value of the slider for the position of the 2-d z-plane
     ;;
     'Z_VALUE': BEGIN
        IF (sState.bShow[2]) THEN BEGIN
           bUpdate = 3
           bRedraw = 1
        ENDIF
     END
     ;;
     ;;  For the vector field: even ordering of the vectors
     ;;
     'SAMPLE_EVEN': BEGIN
        sState.bRandom = 0
        WIDGET_CONTROL, sState.wRandomField, SENSITIVE=0
        bUpdate = 4
        bRedraw = 1
     END
     ;;
     ;;  For the vector field: random ordering of the vectors
     ;;
     'SAMPLE_RANDOM': BEGIN
        sState.bRandom = 1
        WIDGET_CONTROL, sState.wRandomField, SENSITIVE=1
        bUpdate = 4
        bRedraw = 1
     END
     ;;
     ;;  Select how many arrows for even ordering
     ;;
     'N_EVEN': BEGIN
        bUpdate = 4
        bRedraw = 1
     END
     ;;
     ;;  Select how many arrows for random ordering
     ;;
     'N_RANDOM': BEGIN
        bUpdate = 4
        bRedraw = 1
     END
     ;;
     ;;  How long do the arrows have to be?
     ;;
     'SCALE': BEGIN
        bUpdate = 4
        bRedraw = 1
     END
     ;;
     ;;  An event of the draw panel: trackball events or ????
     ;;
     'DRAW': BEGIN
        ;;
        ;; Expose.  (???????)
        ;;
        IF (sEvent.type EQ 4) THEN BEGIN
           demo_draw, sState.oWindow, sState.oView, debug=sState.debug
           WIDGET_CONTROL, sEvent.top, SET_UVALUE=sState, /NO_COPY
           RETURN
        ENDIF
        ;;
        ;; Handle trackball updates.
        ;;
        bHaveTransform = sState.oTrack->Update( sEvent, TRANSFORM=qmat )
        IF (bHaveTransform NE 0) THEN BEGIN
           sState.oGroup->GetProperty, TRANSFORM=t
           sState.oGroup->SetProperty, TRANSFORM=t#qmat
           bRedraw = 1
        ENDIF
        ;;
        ;; Handle other events: PICKING, quality changes, etc.
        ;; Button press.
        ;;
        IF (sEvent.type EQ 0) THEN BEGIN
           IF (sEvent.press EQ 4) THEN BEGIN ; Right mouse.
              IF (d_vectrackDoStream(sEvent,sState,1)) THEN BEGIN
                 bRedraw = 1
                 sState.btndown = 4b
                 sState.oWindow->SetProperty, QUALITY=sState.dragq
                 WIDGET_CONTROL, sState.wDraw, /DRAW_MOTION
              END
           END ELSE IF (sEvent.press EQ 1) THEN BEGIN ; other mouse button.
              sState.btndown = 1b
              sState.oWindow->SetProperty, QUALITY=sState.dragq
              WIDGET_CONTROL, sState.wDraw, /DRAW_MOTION
              bRedraw = 1
           ENDIF ELSE BEGIN     ; middle mouse
              oList = sState.oStreamModel->Get(/ALL)
              sz=size(oList)
              IF (sz[2] EQ 11) THEN BEGIN
                 OBJ_DESTROY,oList
                 bRedraw = 1
              END
           END
        ENDIF
        ;;
        ;; Button motion.
        ;;
        IF (sEvent.type EQ 2) THEN BEGIN
           IF (sState.btndown EQ 4b) THEN BEGIN ; Right mouse button.
              status = d_vectrackDoStream(sEvent,sState,0)
              bRedraw = 1
           ENDIF
        ENDIF
        ;;
        ;; Button release.
        ;;
        IF (sEvent.type EQ 1) THEN BEGIN
           IF (sState.btndown EQ 1b) THEN BEGIN
              sState.oWindow->SetProperty, QUALITY=2
              bRedraw = 1
           END ELSE IF (sState.btndown EQ 4b) THEN BEGIN
              status = d_vectrackDoStream(sEvent,sState,0)
              sState.oWindow->SetProperty, QUALITY=2
              bRedraw = 1
              ;; Build a New Polyline...
           ENDIF
           sState.btndown = 0b
           WIDGET_CONTROL, sState.wDraw, DRAW_MOTION=0
        ENDIF
     END
     ;;
     ;;  Print object to idl.ps
     ;;
     'POSTSCRIPT': BEGIN
        bRedraw = 1
        graphout = 1
     END
     ;;
     ;;  The quit button
     ;;
     'QUIT' : BEGIN
        WIDGET_CONTROL, sEvent.top, SET_UVALUE=sState, /NO_COPY
        WIDGET_CONTROL, sEvent.top, /DESTROY
        RETURN
     end
       ;;;;      'ABOUT' : BEGIN
       ;;;;
       ;;;;            ONLINE_HELP, 'd_vectrack', $
       ;;;;               book=demo_filepath("idldemo.adp", $
       ;;;;                       SUBDIR=['examples','demo','demohelp']), $
       ;;;;                       /FULL_PATH
       ;;;;      end   ; of ABOUT
  ENDCASE
  ;;
  ;;
  ;;  ............................................................
  ;;                    Now do the updates
  ;;  ............................................................
  ;;
  ;;  Update the current display
  ;;
  IF (bUpdate NE 0) THEN BEGIN
     ;;
     ;;  Start by hiding everything
     ;;
     hid = (1-sState.bShow)
     sState.oSlices[0]->SetProperty,HIDE=hid[0]+(1-sState.bImage[0])
     sState.oSlices[1]->SetProperty,HIDE=hid[1]+(1-sState.bImage[1])
     sState.oSlices[2]->SetProperty,HIDE=hid[2]+(1-sState.bImage[2])
     sState.oIsopolygon->SetProperty,HIDE=(1-sState.bIsoShow)
     ;;
     ;;  Update 
     ;;
     d_vectrackPlanesUpdate,sState,bUpdate
     d_vectrackIsoSurfUpdate,sState,bUpdate
  ENDIF
  ;;
  ;; Redraw the graphics
  ;;
  IF (bRedraw) THEN BEGIN
     case graphout of
        0: demo_draw, sState.oWindow, sState.oView, debug=sState.debug
        1: begin
           sState.oClip->Draw, sState.oView, VECTOR=0, POSTSCRIPT=1,$
                                  FILENAME="idl.ps"
        end
     endcase
  ENDIF
  ;;
  ;; Restore state.
  ;;
  WIDGET_CONTROL, sEvent.top, SET_UVALUE=sState, /NO_COPY
  ;;
END




;------------------------------------------------------------------------
;                           MAIN PROGRAM:
;                     RENDER 3-D HYDRO DATA SET
; 
; This is a program that renders 3D data volumes of 1 up to 3 different
; data volumes simultaneously using 2-D cross-sections, iso-surfaces
; and volume rendering. The program is based on the IDL example of the
; visualization of a thunderstorm (d_vectrack.pro), but adapted so that
; it has more general use.
;
; ARGUMENTS:
;   d1         The first 3-D data set (for example: density)
;   d2         The second 3-D data set, which is optional (example: T)
;   d3         The third 3-D data set, which is optional (example: P)
;   velo       A fourth 3-D vector field, which is optional (velocity)
;   dname1     The name of the data set 1
;   dname2     The name of the data set 2
;   dname3     The name of the data set 3
;   xyztitle   The title on the xyz axis
;   xyzrange   The range of the xyz axis
;   iaspect    =1 -> Make sure x and y pixes are equal size
;              =2 -> Make sure x,y and z pixels are equal size
;   loadct     Color table used. If not defined, then internal is used.
;              If set to -1, then it is read from file (if present).
;              If >=0 then use loadct color tables
;   simple     =1, then fewer options
;------------------------------------------------------------------------
PRO render3d, d1,d2,d3,$
                dname1=dname1,$
                dname2=dname2,$
                dname3=dname3,$
                xtitle=xtitle,$
                ytitle=ytitle,$
                ztitle=ztitle,$
                xrange=xrange,$
                yrange=yrange,$
                zrange=zrange,$
                velo=velo,$
                iaspect=iaspect,$
                loadct=loadct,$
                simple=simple,$
                WIDTH = fWidth,  $
                GROUP=group, $               ; IN: (opt) group identifier
                DEBUG=debug, $               ; IN: (opt)
                APPTLB = appTLB              ; OUT: (opt) TLB of this application
  ;;
  ;; For ribbon-style stream lines  (if active)
  ;;
  IF (N_ELEMENTS(fWidth) EQ 0) THEN fWidth = 1.0
  ;;
  ;; Defaults
  ;;
  if not keyword_set(xtitle) then xtitle='X'
  if not keyword_set(ytitle) then ytitle='Y'
  if not keyword_set(ztitle) then ztitle='Z'
  if n_elements(loadct) eq 0 then loadct=-3
  opacmax = 32
  ;;
  ;; Check which kind of data is given
  ;;
  sz= size(velo)
  if sz[0] eq 4 then begin
     if(sz[4] eq 3) then begin
        u = velo[*,*,*,0]
        v = velo[*,*,*,0]
        w = velo[*,*,*,0]
     endif else begin
        u = REFORM(velo[0,*,*,*])
        v = REFORM(velo[0,*,*,*])
        w = REFORM(velo[0,*,*,*])
     endelse
     m = sqrt(u^2+v^2+w^2)
     incl_velo = 1
  endif else begin
     m = 0.d0
     u = 0.d0
     v = 0.d0
     w = 0.d0
     incl_velo = 0
  endelse
  ;;
  ;; Check how many datasets are given
  ;;
  if n_elements(d1) eq 0 then begin
     print,'At least one scalar data block must be specified'
     stop
  endif
  nr_data = 1
  if n_elements(d2) gt 0 then nr_data = 2
  if n_elements(d3) gt 0 then nr_data = 3
  nr_data_tot = nr_data
  if incl_velo ne 0 then nr_data_tot = nr_data_tot + 1
  ;;
  ;; Check validity of data sets
  ;;
  if nr_data ge 1 then begin
     if max(d1)-min(d1) eq 0.d0 then begin
        print,'ERROR: Data set 1 has 0 range'
     endif
  endif
  if nr_data ge 2 then begin
     if max(d2)-min(d2) eq 0.d0 then begin
        print,'ERROR: Data set 2 has 0 range'
     endif
  endif
  if nr_data ge 3 then begin
     if max(d3)-min(d3) eq 0.d0 then begin
        print,'ERROR: Data set 3 has 0 range'
     endif
  endif
  ;;
  ;; Set the names of the datasets
  ;;
  datanames = ['data1','data2','data3']
  if keyword_set(dname1) then datanames[0]=dname1
  if keyword_set(dname2) then datanames[1]=dname2
  if keyword_set(dname3) then datanames[2]=dname3
  datanames = datanames[0:nr_data-1]
  if incl_velo ne 0 then datanames=[datanames,'|v|']
  ;;
  ;;  The scales (always spanning the entire data range)
  ;;
  fScale = dblarr(2,nr_data+1)
  if nr_data ge 1 then fscale[*,0] = [min(d1),max(d1)]
  if nr_data ge 2 then fscale[*,1] = [min(d2),max(d2)]
  if nr_data ge 3 then fscale[*,2] = [min(d3),max(d3)]
  if incl_velo ne 0 then fscale[*,nr_data_tot-1] = [min(m),max(m)]
  ;;
  ;;  Dependent on the range, we choose the format of the numbers
  ;;
  dformat = strarr(nr_data_tot)
  for i=0,nr_data_tot-1 do begin
     if max(abs(fscale[*,i])) gt 1d3 or $
        min(abs(fscale[*,i])) lt 1d-3 then begin
        dformat[i] = '(E13.6)'
     endif else begin
        dformat[i] = '(F9.4)'
     endelse
  endfor
  ;;
  ;; Create the byte scaled versions of the data sets
  ;; which are necessary for the cross-cut images
  ;;
  ;; NOTE: Here I hard-set that the color scale covers the entire 
  ;;       range of the data. if a sub-range is required, then 
  ;;       the code must be modified here
  ;;
  if incl_velo ne 0 then mb = bytscl(m)
  if nr_data ge 1 then db1 = bytscl(d1)
  if nr_data ge 2 then db2 = bytscl(d2)
  if nr_data ge 3 then db3 = bytscl(d3)
  ;;
  ;; For now (in this version) we have only 1 time step of the data
  ;;
  nSteps = 0
  pData = PTR_NEW()
  ;;
  ;; Check the validity of the group identifier, in case this application
  ;; is part of a bigger application
  ;;
  ngroup = N_ELEMENTS(group)
  if (ngroup NE 0) then begin
     check = WIDGET_INFO(group, /VALID_ID)
     if (check NE 1) then begin
        print,'Error, the group identifier is not valid'
        print, 'Return to the main application'
        RETURN
     endif
     groupBase = group
  endif else groupBase = 0L
  ;;
  ;;  Get the screen size.
  ;;
  Device, GET_SCREEN_SIZE = screenSize
  ;;
  ;;  Set up dimensions of the drawing (viewing) area.
  ;;
  xdim = screenSize[0]*0.6
  ydim = xdim*0.8
  ;;
  ;;  Get the current color vectors to restore
  ;;  when this application is exited.
  ;;
  TVLCT, savedR, savedG, savedB, /GET
  ;;
  ;;  Build color table from color vectors
  ;;
  colorTable = [[savedR],[savedG],[savedB]]
  ;;
  ;; Get the data size (take the p variable as an example)
  ;;
  sz = SIZE(d1)
  ;;
  ;;  ..............................................
  ;;               Create widgets.
  ;;  ..............................................
  ;;
  if (N_ELEMENTS(group) EQ 0) then begin
     wTopBase = WIDGET_BASE(/COLUMN, $
            TITLE="Data Rendering Tool in 3-D", $
            XPAD=0, YPAD=0, $
            /TLB_KILL_REQUEST_EVENTS, $
            TLB_FRAME_ATTR=1, MBAR=barBase, $
            UNAME='D_VECTRACK:tlb')
  endif else begin
     wTopBase = WIDGET_BASE(/COLUMN, $
            TITLE="Data Rendering Tool in 3-D", $
            XPAD=0, YPAD=0, $
            /TLB_KILL_REQUEST_EVENTS, $
            GROUP_LEADER=group, $
            TLB_FRAME_ATTR=1, MBAR=barBase, $
            UNAME='D_VECTRACK:tlb')
  endelse
  ;;
  ;;  Create the menu bar. It contains the file/quit,
  ;;  edit/ shade-style, help/about.
  ;;
  wFileButton = WIDGET_BUTTON(barBase, VALUE='File', /MENU)
  wPostscriptButton = WIDGET_BUTTON(wFileButton, $
                VALUE='Postscript', UVALUE='POSTSCRIPT')
  wQuitButton = WIDGET_BUTTON(wFileButton, $
                VALUE='Quit', UVALUE='QUIT')
  ;;
  ;;  Create the menu bar item help that contains the about button
  ;;
  ;;;; wHelpButton = WIDGET_BUTTON(barBase, VALUE='About', /HELP, /MENU)
  ;;;; wAboutButton = WIDGET_BUTTON(wHelpButton, $
  ;;;;              VALUE='About Thunderstorm Visualization', UVALUE='ABOUT')
  ;;
  ;;  Now arrange the basic geometry of the GUI
  ;;
  wTopRowBase = WIDGET_BASE(wTopBase,/ROW,/FRAME)
  wGuiBase = WIDGET_BASE(wTopRowBase, /COLUMN)
  wRowBase = WIDGET_BASE(wGuiBase, /ROW )
  wFrameBase = WIDGET_BASE(wRowBase, /COLUMN, /FRAME)
  wLabel = WIDGET_LABEL(wFrameBase,VALUE='Cross-Section Planes:')
  wGuiBase2 = WIDGET_BASE(wFrameBase,/COLUMN)
  ;;
  ;;  Now the droplists for what kind of thing to plot in each direction
  ;;  (i.e. a vector field, an image or nothing)
  ;;
  if incl_velo ne 0 then begin
     stringarray = ['<Off>','Image','Vector']
  endif else begin
     stringarray = ['<Off>','Image']
  endelse
  wXDrop = WIDGET_DROPLIST(wGuiBase2, VALUE=stringarray,$
                           TITLE='X:', UVALUE='X_STYLE', $
                           UNAME='D_VECTRACK:xdroplist')
  wYDrop = WIDGET_DROPLIST(wGuiBase2, VALUE=stringarray,$
                           TITLE='Y:', UVALUE='Y_STYLE', $
                           UNAME='D_VECTRACK:ydroplist')
  wZDrop = WIDGET_DROPLIST(wGuiBase2, VALUE=stringarray,$
                           TITLE='Z:', UVALUE='Z_STYLE', $
                           UNAME='D_VECTRACK:zdroplist')
  ;;
  ;;  Now the sliders for the location of each of the 2d planes
  ;;
  wXSlider = WIDGET_SLIDER(wFrameBase, MAXIMUM=sz[1]-1, $
                           TITLE='X Plane', UVALUE='X_VALUE')
  wYSlider = WIDGET_SLIDER(wFrameBase, MAXIMUM=sz[2]-1, $
                           value=sz[2]/2, $
                           TITLE='Y Plane', UVALUE='Y_VALUE')
  wZSlider = WIDGET_SLIDER(wFrameBase, MAXIMUM=sz[3]-1, $
                           TITLE='Z Plane', UVALUE='Z_VALUE')
  ;;
  ;;  Determine which operating system we have
  ;;
  frame = !version.os_family ne 'unix' ; Workaround problem 7763.
  ;;
  ;;  A new base 
  ;;
  wFrameBase = WIDGET_BASE(wRowBase, /COLUMN, /FRAME)
  ;;
  ;;  The selectors which thing to render/plot
  ;;
  idum = nr_data
  if incl_velo ne 0 then idum=idum+1
  if idum gt 1 then begin
     stringarray = datanames
     wVVolSel = WIDGET_DROPLIST(wFrameBase,VALUE=stringarray,$
                FRAME=frame, TITLE='Vrend vol',UVAL='VREND_VOLSEL')
     wIVolSel = WIDGET_DROPLIST(wFrameBase,VALUE=stringarray,$
                FRAME=frame,TITLE='Img vol',UVAL='IMG_VOLSEL')
     wSVolSel = WIDGET_DROPLIST(wFrameBase,VALUE=stringarray,$
                FRAME=frame, TITLE='Iso vol',UVAL='ISO_VOLSEL')
  endif
  ;;
  ;;  The button to activate/deactivate the iso-surface
  ;;
  wGuiBase2 = WIDGET_BASE(wFrameBase,column=2,/NONEXCLUSIVE)
  wIsotoggle = WIDGET_BUTTON(wGuiBase2,VALUE='Iso',UVALUE='ISO_SHOW')
  wVoltoggle = WIDGET_BUTTON(wGuiBase2,VALUE='Vol',UVALUE='VOL_SHOW')
  ;;
  ;; Volume rendering button
  ;;
  ;;wVRender = WIDGET_BUTTON(wFrameBase, VALUE='Vol Render', $
  ;;              UVALUE='VREND', UNAME='D_VECTRACK:volrendr')
  ;;
  ;;  The slider determining the level of the volume iso-surface
  ;;
  wIsoText = WIDGET_LABEL(wFrameBase,VALUE=$
          strcompress(string(0.5*(fscale[0,0]+fscale[1,0]),$
          format=dformat[0])))
  wIsoLevel = WIDGET_SLIDER(wFrameBase,/SUPPRESS_VALUE, $
                TITLE='Level',UVAL='ISO_LEVEL',MAXIMUM=255,VALUE=128,$
                UNAME='D_VECTRACK:iso_level')
  ;;
  ;;  Sliders for determining the transparency of the volume rendering (??)
  ;;
  wAlphamin = WIDGET_SLIDER(wFrameBase,/SUPPRESS_VALUE,UVAL='ALPHA_LEVEL', $
                MAXIMUM=255,VALUE=0, TITLE='2-D Plane Transp Min.')
  wAlphamax = WIDGET_SLIDER(wFrameBase,/SUPPRESS_VALUE,UVAL='ALPHA_LEVEL', $
                MAXIMUM=255,VALUE=255, TITLE='2-D Plane Transp Max.')
  ;;
  ;;  Opacity slider
  ;;
  if not keyword_set(simple) then begin
     wOpacityMax = WIDGET_SLIDER(wFrameBase,/SUPPRESS_VALUE,UVAL='OPACITY_MAX', $
                MAXIMUM=255,VALUE=opacmax, TITLE='Opacity Max')
  endif else begin
     wOpacityMax = 0
  endelse
  ;;
  ;;  Color table selector
  ;;
  if not keyword_set(simple) then begin
     wColorTabSel = WIDGET_SLIDER(wFrameBase,/SUPPRESS_VALUE,UVAL='COLOR_TABLE', $
                MINIMUM=-3,MAXIMUM=40,VALUE=loadct, TITLE='Color Table')
  endif else begin
     wColorTabSel = 0
  endelse
  
  ;;
  ;;  The stuff for the streamlines (disabled???)
  ;;
  ;;;; wClearStreams = WIDGET_BUTTON(wFrameBase,VALUE='Clear Streamlines', $
  ;;;;               UVALUE='CLEAR_STREAMS',UNAME='D_VECTRACK:clearstreams')
  ;;;; wRibbons = CW_BGROUP(wFrameBase, ['Ribbons'], UVALUE='RIBBON', $
  ;;;;               /NONEXCLUSIVE)
  ;;;; WIDGET_CONTROL, wRibbons, SET_VALUE=1
  ;;
  ;;  Now the stuff for the vectors
  ;;
  if incl_velo ne 0 then begin
     wFrameBase = WIDGET_BASE(wGuiBase, /COLUMN, /FRAME)
     wLabel = WIDGET_LABEL(wFrameBase,VALUE='Vector Sampling:')
     wGuiBase2 = WIDGET_BASE(wFrameBase,/ROW,/EXCLUSIVE)
     wEButton = WIDGET_BUTTON(wGuiBase2, VALUE='Even',$
            UVALUE='SAMPLE_EVEN',/NO_RELEASE, $
            UNAME='D_VECTRACK:even')
     wRButton = WIDGET_BUTTON(wGuiBase2, VALUE='Random',$
            UVALUE='SAMPLE_RANDOM',/NO_RELEASE, $
            UNAME='D_VECTRACK:random')
     wEvenField = CW_FIELD(wFrameBase, /INTEGER, $
            TITLE='Sample Every Nth, N=', XSIZE=2, $
            UVALUE='N_EVEN',VALUE=2,/RETURN_EVENTS)
     wRandomField = CW_FIELD(wFrameBase, /INTEGER, $
            TITLE='Total Number of Samples=', XSIZE=5, $
            UVALUE='N_RANDOM',VALUE=500,/RETURN_EVENTS)
     wScaleField = CW_FIELD(wGuiBase, /FLOAT, TITLE='Vector Length = ', $
            UVALUE='SCALE', VALUE=10.0, XSIZE=8, /RETURN_EVENTS)
  endif else begin
     wEvenField = 0
     wRandomField = 0
     wScaleField = 0
  endelse
  ;;
  ;;  If there exists more than 1 time snapshot, then here is the 
  ;;  time slider
  ;;
  IF (nSteps GT 0) THEN BEGIN
     wTimeStep = WIDGET_SLIDER(wGuiBase,UVAL='TIME_STEP', $
            MAXIMUM=nSteps-1,VALUE=0)
  END
  ;;
  ;; Now set up the main drawing panel
  ;;
  wDraw = WIDGET_DRAW(wTopRowBase, XSIZE=xdim, YSIZE=ydim, UVALUE='DRAW', $
                        RETAIN=0, /EXPOSE_EVENTS, /BUTTON_EVENTS, $
                        GRAPHICS_LEVEL=2, UNAME='D_VECTRACK:draw')
  ;;
  ;;  Realize the base widget.
  ;;
  WIDGET_CONTROL, wTopBase, /REALIZE
  ;;
  ;;  Returns the top level base in the appTLB keyword
  ;;
  appTLB = wTopBase
  ;;
  ;;  Get the window id of the drawable.
  ;;
  WIDGET_CONTROL, wDraw, GET_VALUE=oWindow
  ;;
  ;; Set default widget items.
  ;;
  WIDGET_CONTROL, wXDrop, SET_DROPLIST_SELECT=2
  WIDGET_CONTROL, wYDrop, SET_DROPLIST_SELECT=1
  WIDGET_CONTROL, wZDrop, SET_DROPLIST_SELECT=2
  if incl_velo ne 0 then begin
     WIDGET_CONTROL, wEButton,SET_BUTTON=1
     WIDGET_CONTROL, wRButton,SET_BUTTON=0
     WIDGET_CONTROL, wRandomField, SENSITIVE=0
  endif
  ;;
  ;;  ..............................................
  ;;            Create object graphics
  ;;  ..............................................
  ;;
  ;;  Compute viewplane rect based on aspect ratio.
  ;;
  aspect = FLOAT(xdim) / FLOAT(ydim)
  myview = [-0.5, -0.5, 1, 1] * 1.5
  IF (aspect > 1) THEN BEGIN
     myview[0] = myview[0] - ((aspect-1.0)*myview[2])/2.0
     myview[2] = myview[2] * aspect
  ENDIF ELSE BEGIN
     myview[1] = myview[1] - (((1.0/aspect)-1.0)*myview[3])/2.0
     myview[3] = myview[3] * aspect
  ENDELSE
  ;;
  ;;  Drop the view down a bit to make room for the colorbar
  ;;
  myview[1] = myview[1] - 0.15
  ;;
  ;;  Create (for postscript dump) a clipboard object
  ;;
  oClip = OBJ_NEW('IDLgrClipboard')
  ;;
  ;;  Create view.
  ;;
  oView = OBJ_NEW('IDLgrView', PROJECTION=2,$
                    VIEWPLANE_RECT=myview,COLOR=[50,50,70])
  ;;
  ;;  Create model.
  ;;
  oTop = OBJ_NEW('IDLgrModel')
  oTop->Scale,0.8,0.8,0.8
  oGroup = OBJ_NEW('IDLgrModel')
  oTop->Add, oGroup
  ;;
  ;;  Compute data bounds.
  ;;
  sz = SIZE(d1)
  xMax = sz[1] - 1
  yMax = sz[2] - 1
  zMax = sz[3] - 1
  ;;
  ;; Compute coordinate conversion to normalize.
  ;; Note: x and y are here forced to same scale. z has own scale
  ;;
  if not keyword_set(iaspect) then iaspect=0
  case iaspect of 
     0: begin
        maxDim = MAX([xMax, yMax])
        xs = [-0.5 * xMax/xMax,1.0/xMax]
        ys = [-0.5 * yMax/yMax,1.0/yMax]
        zoom = 0.7
        zs = [-0.5*zMax/zmax, 1.0/zmax]*zoom
     end
     1: begin
        maxDim = MAX([xMax, yMax])
        xs = [-0.5 * xMax/maxDim,1.0/maxDim]
        ys = [-0.5 * yMax/maxDim,1.0/maxDim]
        zoom = 0.7
        zs = [-0.5*zMax/zmax, 1.0/zmax]*zoom
     end
     2: begin
        maxDim = MAX([xMax, yMax, zMax])
        xs = [-0.5 * xMax/maxDim,1.0/maxDim]
        ys = [-0.5 * yMax/maxDim,1.0/maxDim]
        zoom = 0.7
        zs = [-0.5*zMax/maxDim, 1.0/maxDim]*zoom
     end
  endcase
  ;;
  ;;  Create the X-axis object
  ;;
  if keyword_set(xrange) then begin
     range=xrange
     cocon=[-0.5-(range[0])/(range[1]-range[0]),1.d0/(range[1]-range[0])]
  endif else begin
     range=[0,xMax]
     cocon=[xs[0],xs[1]]
  endelse
  oXTitle = OBJ_NEW('IDLgrText', xtitle)
  oAxis = OBJ_NEW('IDLgrAxis', 0, COLOR=[255,255,255],RANGE=range,$
             TITLE=oXTitle,TICKLEN=4*(ymax/100.),/EXACT,extend=0, $
             XCOORD_CONV=cocon, $
             YCOORD_CONV=ys, $
             ZCOORD_CONV=zs)
  oGroup->Add, oAxis
  ;;
  ;;  Create the Y-axis object
  ;;
  if keyword_set(yrange) then begin
     range=yrange
     cocon=[-0.5-(range[0])/(range[1]-range[0]),1.d0/(range[1]-range[0])]
  endif else begin
     range=[0,yMax]
     cocon=[ys[0],ys[1]]
  endelse
  oYTitle = OBJ_NEW('IDLgrText', ytitle)
  oAxis = OBJ_NEW('IDLgrAxis', 1, COLOR=[255,255,255],RANGE=range,$
             TITLE=oYTitle,TICKLEN=4*(xmax/100.),/EXACT,extend=0, $
             XCOORD_CONV=xs, $
             YCOORD_CONV=cocon, $
             ZCOORD_CONV=zs)
  oGroup->Add, oAxis
  ;;
  ;;  Create the Z-axis object
  ;;
  if keyword_set(zrange) then begin
     range=zrange
     cocon=[-0.5-(range[0])/(range[1]-range[0]),1.d0/(range[1]-range[0])]
  endif else begin
     range=[0,zMax]
     cocon=[zs[0],zs[1]]
  endelse
  oZTitle = OBJ_NEW('IDLgrText', ztitle)
  oAxis = OBJ_NEW('IDLgrAxis', 2, COLOR=[255,255,255],RANGE=range,$
             /EXACT, extend=0, $
             TICKLEN=4*(xmax/100.),TITLE=oZTitle,$
             XCOORD_CONV=xs, $
             YCOORD_CONV=[.5 * yMax/yMax,1.0/yMax], $
             ZCOORD_CONV=cocon)
  oGroup->Add, oAxis
  ;;
  ;; wireframe box
  ;;
  oBox = OBJ_NEW('IDLgrPolyline', $
        [[0,0,0],[xMax,0,0],[0,yMax,0],[xMax,yMax,0], $
         [0,0,zMax],[xMax,0,zMax],[0,yMax,zMax],$
         [xMax,yMax,zMax]], $
        COLOR=[200,200,200], $
        POLYLINE=[5,0,1,3,2,0,5,4,5,7,6,4,2,0,4,2,1,5,2,2,6,2,3,7],$
        XCOORD_CONV=xs, YCOORD_CONV=ys, ZCOORD_CONV=zs)
  oGroup->Add, oBox
  ;;
  ;;  slice image-form objects (texture mapped quadratic polygons)
  ;;
  oXImage = OBJ_NEW('IDLgrImage',dist(5),HIDE=1)
  oYImage = OBJ_NEW('IDLgrImage',dist(5),HIDE=1)
  oZImage = OBJ_NEW('IDLgrImage',dist(5),HIDE=1)
  oXPolygon = OBJ_NEW('IDLgrPolygon', COLOR=[255,255,255], HIDE=1,$
              TEXTURE_MAP=oXImage, TEXTURE_INTERP=1, $
                  XCOORD_CONV=xs, YCOORD_CONV=ys, ZCOORD_CONV=zs)
  oYPolygon = OBJ_NEW('IDLgrPolygon', COLOR=[255,255,255], HIDE=1,$
              TEXTURE_MAP=oYImage, TEXTURE_INTERP=1, $
                  XCOORD_CONV=xs, YCOORD_CONV=ys, ZCOORD_CONV=zs)
  oZPolygon = OBJ_NEW('IDLgrPolygon', COLOR=[255,255,255], HIDE=1,$
              TEXTURE_MAP=oZImage, TEXTURE_INTERP=1, $
                  XCOORD_CONV=xs, YCOORD_CONV=ys, ZCOORD_CONV=zs)
  oGroup->Add, oXPolygon
  oGroup->Add, oYPolygon
  oGroup->Add, oZPolygon
  ;;
  ;;  Iso Surface objects
  ;;
  oIsoPolygon = OBJ_NEW('IDLgrPolygon', COLOR=[127,127,127], HIDE=1,$
            SHADING=1,XCOORD_CONV=xs, YCOORD_CONV=ys, ZCOORD_CONV=zs)
  oGroup->Add, oIsoPolygon
  ;;
  ;; Stream line objects holder
  ;;  ????? DISABLED ???? 
  ;;
  oStreamModel = OBJ_NEW('IDLgrModel')
  oGroup->Add, oStreamModel
  ;;
  ;;  Volume objects MUST be added LAST
  ;;
  if nr_data ge 1 then begin
     od1vol = OBJ_NEW('IDLgrVolume', DATA0 = db1, $
            /ZERO_OPACITY_SKIP, /ZBUFFER, HIDE=1, $
            XCOORD_CONV=xs, YCOORD_CONV=ys, ZCOORD_CONV=zs)
  endif
  if nr_data ge 2 then begin
     od2vol = OBJ_NEW('IDLgrVolume', DATA0 = db2, $
            /ZERO_OPACITY_SKIP, /ZBUFFER, HIDE=1, $
            XCOORD_CONV=xs, YCOORD_CONV=ys, ZCOORD_CONV=zs)
  endif
  if nr_data ge 3 then begin
     od3vol = OBJ_NEW('IDLgrVolume', DATA0 = db3, $
            /ZERO_OPACITY_SKIP, /ZBUFFER, HIDE=1, $
            XCOORD_CONV=xs, YCOORD_CONV=ys, ZCOORD_CONV=zs)
  endif
  if incl_velo ne 0 then begin
     oMvol = OBJ_NEW('IDLgrVolume', DATA0 = mb, $
            /ZERO_OPACITY_SKIP, /ZBUFFER, HIDE=1, $
            XCOORD_CONV=xs, YCOORD_CONV=ys, ZCOORD_CONV=zs)
  endif
  ;;
  ;; add the voxel volumes
  ;;
  if nr_data ge 1 then oGroup->Add, od1vol
  if nr_data ge 2 then oGroup->Add, od2vol
  if nr_data ge 3 then oGroup->Add, od3vol
  if incl_velo ne 0 then oGroup->Add, oMvol
  ;;
  ;; setup their palettes/opacity tables
  ;;
  if nr_data ge 1 then d_vectrackReadVoxelPalettes,od1vol,coltabnr=loadct
  if nr_data ge 2 then d_vectrackReadVoxelPalettes,od2vol,coltabnr=loadct
  if nr_data ge 3 then d_vectrackReadVoxelPalettes,od3vol,coltabnr=loadct
  if incl_velo ne 0 then d_vectrackReadVoxelPalettes,oMvol,coltabnr=loadct
  ;;
  ;; For 0-MAX volumes, replace the opacity table
  ;;    ***** TO BE MODIFIED *****
  ;;
  if incl_velo ne 0 then begin
     oMvol->SetProperty,OPACITY_TABLE0=((findgen(256)/16)^2.0)/4.0
  endif
  ;;
  ;; Grab the Volume color table for use by the image objects
  ;;
  od1vol->GetProperty,RGB_TABLE0=ctab
  oPal = OBJ_NEW('IDLgrPalette',ctab[*,0],ctab[*,1],ctab[*,2])
  oXImage->SetProperty,PALETTE=oPal
  oYImage->SetProperty,PALETTE=oPal
  oZImage->SetProperty,PALETTE=oPal
  ;;
  ;; Create some lights.
  ;;
  oLight = OBJ_NEW('IDLgrLight', LOCATION=[2,2,2], TYPE=1, INTENSITY=0.8)
  oTop->Add, oLight
  oLight = OBJ_NEW('IDLgrLight', TYPE=0, INTENSITY=0.5)
  oTop->Add, oLight
  ;;
  ;; Create the color bar annotation (lower left corner)
  ;;
  oCBTop = OBJ_NEW('IDLgrModel')
  rgb = bytarr(3,256,16)
  rgb[0,*,*] = indgen(256*16) MOD 256
  rgb[1,*,*] = indgen(256*16) MOD 256
  rgb[2,*,*] = indgen(256*16) MOD 256
  oImage = OBJ_NEW('IDLgrImage',rgb,DIMENSIONS=[256,16])
  oCBTitle = OBJ_NEW('IDLgrText','Grayscale')
  oAxis = OBJ_NEW('IDLgrAxis',range=[0,256],/EXACT,COLOR=[255,255,255], $
             TICKLEN=15,MAJOR=5,TITLE=oCBTitle)
  oCBTop->Add,oImage
  oCBTop->Add,oAxis
  oView->Add,oCBTop
  ;;
  ;; Place the model in the view.
  ;;
  oView->Add, oTop
  ;;
  ;; Rotate to standard view for first draw.
  ;;
  oGroup->Rotate, [1,0,0], -90
  oGroup->Rotate, [0,1,0], 30
  oGroup->Rotate, [1,0,0], 30
  ;;
  ;; Create a trackball.
  ;;
  oTrack = OBJ_NEW('Trackball', [xdim/2, ydim/2.], xdim/2.)
  ;;
  ;; Create a holder object for easy destruction.
  ;;
  oHolder = OBJ_NEW('IDL_Container')
  oHolder->Add, oView
  oHolder->Add, oTrack
  oHolder->Add, oXTitle
  oHolder->Add, oYTitle
  oHolder->Add, oZTitle
     ;;;;    oHolder->Add, oTickText
  oHolder->Add, oCBTitle
  oHolder->Add, oPal
  oHolder->Add, oXImage
  oHolder->Add, oYImage
  oHolder->Add, oZImage
  oHolder->Add, oClip
  ;;
  ;;  ..............................................
  ;;   Now put all data in state and start widgets
  ;;  ..............................................
  ;;
  if incl_velo ne 0 then begin
     vdims = SIZE(u, /DIMENSIONS)
     vdata = FLTARR(3, vdims[0], vdims[1], vdims[2])
     vdata[0, *, *, *] = FLOAT(u)
     vdata[1, *, *, *] = FLOAT(v)
     vdata[2, *, *, *] = FLOAT(w)
  endif else begin
     vdata = 0.d0
  endelse
  ;;
  ;; Bundle all volume rendering data
  ;; 
  if nr_data ge 1 then oVols   = [od1vol]
  if nr_data ge 2 then oVols   = [oVols,od2vol]
  if nr_data ge 3 then oVols   = [oVols,od3vol]
  if incl_velo ne 0 then begin
     oVols = [oVols,oMvol]
  endif
  ;;
  ;; Save state.
  ;;
  sState = {nSteps: nSteps,        $
        curStep: 0,            $
        pSteps: pData,         $
        btndown: 0b,           $
        dragq: 1,          $
        u:u,                   $
        v:v,                   $
        w:w,                   $
        vdata:vdata,           $
        bShow: [0b,0b,0b],     $
        bImage: [0b,0b,0b],    $
        bRandom: 0b,           $
        oHolder: oHolder,      $
        oTrack:oTrack,         $
        wXSlider:wXSlider,     $
        wYSlider:wYSlider,     $
        wZSlider:wZSlider,     $
        wEvenField:wEvenField, $
        wRandomField:wRandomField, $
        wScaleField:wScaleField, $
        wIsoLevel: wIsoLevel,  $
        wIsoText: wIsoText,    $
        wAlpha: [wAlphamin,wAlphamax],    $
        wDraw: wDraw,          $
        wLabel: wLabel,        $
;;;;        wRibbons: wRibbons, $
        oWindow: oWindow,      $
        oClip:oClip,           $
        oGroup: oGroup,        $
        oCBTop: oCBTop,        $
        iAlpha: [0,255],       $
        iImgVol: 0,            $
        iVrendVol: 0,          $
        oVols: oVols, $
        oSlices: [oXPolygon,oYPolygon,oZPolygon], $
        oIsopolygon: oIsopolygon, $
        oStreamModel: oStreamModel, $
        oStreamline: OBJ_NEW(), $
        oImages: [oXImage,oYImage,oZImage], $
        iIsoVol: 0,            $
        fIsoLevel: 128.0,      $
        fWidth: fWidth, $
        bIsoShow: 0,           $
        bVolShow: 0,           $
        fScale:fScale, $
        ColorTable: colorTable, $ ; Color table to restore at exit
        debug: keyword_set(debug), $
        groupBase: groupBase, $               ; Base of Group Leader
            ;;
            ;;  Here comes stuff by C.P. Dullemond
            ;;
            incl_velo: incl_velo, $
            datanames: datanames, $
            nr_data: nr_data, $
            nr_data_tot: nr_data_tot, $
            dformat: dformat, $
            wOpacityMax:wOpacityMax, $
            opacmax: opacmax, $
            loadct: loadct, $
            wColorTabSel:wColorTabSel, $
            oView: oView      $
        }
  ;;
  ;; Initialize with an interesting view
  ;;
  sState.bShow[1]=1
  sState.bImage[1]=1
  sState.bIsoShow=1
  widget_control, wIsotoggle, /set_button
  d_vectrackColorBarUpdate,sState
  d_vectrackIsoSurfUpdate,sState,5
  d_vectrackPlanesUpdate,sState,4
  sState.oSlices[1]->SetProperty,HIDE=0
  sState.oIsopolygon->SetProperty,HIDE=0
  ;;
  ;; Put sstate to the main widget
  ;;
  WIDGET_CONTROL, wTopBase, SET_UVALUE=sState, /NO_COPY
  ;;
  ;; Now blow life into the whole thing
  ;;
  XMANAGER, 'RENDER3D', wTopBase, $
        EVENT_HANDLER='d_vectrackEvent', $
        /NO_BLOCK, $
        CLEANUP='d_vectrackCleanup'
  ;;
END





