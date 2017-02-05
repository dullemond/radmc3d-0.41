;==============================================================================
;
; read a FLASH HDF5 file into IDL.  
;
; This routine reads in FLASH data in the HDF 5 format.  
;
; This version is for IDL >= 5.6, as it uses the built in HDF5
; support.
;
; Arguments:  filename -- the name of the file to read in
;
;             var_name -- a 4-character string specifying the variable
;                         name to read in.  If this keyword is not
;                         present, all variables are read in.
;
;             double   -- read the variables in in double precision
;
;==============================================================================

pro read_amr_hdf5, filename, $
                   VAR_NAME=var_name, $
                   DOUBLE=double, $
                   TREE=tree, $
                   DATA=unk, $
                   PARAMETERS=params, $
                   STORED_VARS=unk_names, $
                   NUM_PARTICLES=numParticles, $
                   PARTICLES=particles, $
                   INT_PROP_NAMES=IntPropNames, $
                   REAL_PROP_NAMES=RealPropNames, $
                   GEOMETRY=geometry


; clean up any pre-existing arrays
if n_elements(unk) GT 0 then undefine, unk
if n_elements(particles) gt 0 then undefine, particles

; check to see which variable to read; if var is not defined, read them all
if n_elements(var_name) EQ 0 then var_name = 'none'
print, 'var_name = ', var_name

;itype = determine_file_type(filename)
itype = 1
;if(itype EQ -1) then begin
;  message,'Error: Problem determining filetype of file ' + filename
;endif

if(itype EQ 1) then begin


;------------------------------------------------------------------------------
; open up the file for read now and read in the header information
;------------------------------------------------------------------------------

print, 'read_amr_hdf5: try to open ',filename

file_identifier = H5F_OPEN(filename)

;------------------------------------------------------------------------------
; grab the list of variables, and make sure that our requested
; variable (if any) exists in the dataset
;------------------------------------------------------------------------------

dataset = H5D_OPEN(file_identifier, "unknown names")

unk_names =  H5D_READ(dataset)
unk_names = reform(temporary(unk_names))

H5D_CLOSE, dataset



; we can also now set the number of variables
nvar = (size(unk_names))[1]

if (var_name NE 'none') then begin
    var = (where(unk_names EQ var_name))[0]

    if (var EQ -1) then begin
        print, 'ERROR: requested variable not found in dataset'
    endif

    print, 'reading in only the variable ', var_name
endif else begin
    var = -1
endelse



; read in the simulation paramters
dataset = H5D_OPEN(file_identifier, "simulation parameters")
datatype = H5D_GET_TYPE(dataset)

idl_type = H5T_IDLTYPE(datatype, STRUCTURE=sim_params)

sim_params = H5D_READ(dataset)

H5D_CLOSE, dataset
H5T_CLOSE, datatype



; get the dimensionality
dataset = H5D_OPEN(file_identifier, "coordinates")
dataspace = H5D_GET_SPACE(dataset)
dims = H5S_GET_SIMPLE_EXTENT_DIMS(dataspace)

ndim = dims[0]

H5D_CLOSE,dataset
H5S_CLOSE, dataspace



; figure out if we are dealing with corners
corners = 0

if (strpos(filename, 'crn_') GT 0) then begin
    corners = 1
endif

nfaces = 2*ndim
nchild = 2^ndim


if (n_elements(geometry) eq 0) then begin
; get the geometry
    group = H5G_OPEN(file_identifier, "/")
    attribute = H5A_OPEN_NAME(group, "geometry name")
    
    geometry = H5A_READ(attribute)
        
    H5A_CLOSE, attribute
    H5G_CLOSE, group
endif 

; ---- setup the structures to pass the data ----------------------------------
if (n_elements(double) EQ 0) then double = 0

params = {totBlocks:sim_params.total_blocks, $
          corners:corners, $
          ndim:ndim, $
          nvar:nvar, $
          nxb:sim_params.nxb, $
          nyb:sim_params.nyb, $
          nzb:sim_params.nzb, $
          ntopx:1, $
          ntopy:1, $
          ntopz:1, $
          time:sim_params.time, $
          dt:sim_params.timestep, $
          redshift:1.d0, $
          geometry:"unknown"}

params.geometry = strcompress(geometry, /REMOVE_ALL)

; for some reason, the redshift field is not stored in the plotfiles
if ( (where(tag_names(sim_params) EQ "REDSHIFT"))[0] NE -1) then begin
        params.redshift = sim_params.redshift
endif


;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
; We now have the dimension, total number of blocks, and the number of
; zones per block in each direction.  We can use this to read the tree
; information, coordinate information, and unknowns
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


;----------------------------------------------------------------------------
; read in the coordinate infomation
;----------------------------------------------------------------------------

dataset = H5D_OPEN(file_identifier, "coordinates")
coord = H5D_READ(dataset)
H5D_CLOSE, dataset

dataset = H5D_OPEN(file_identifier, "block size")
size = H5D_READ(dataset)
H5D_CLOSE, dataset

dataset = H5D_OPEN(file_identifier, "bounding box")
bnd_box = H5D_READ(dataset)
H5D_CLOSE, dataset

help, coord, size, bnd_box

if size(coord,/type) EQ 4 then dbl_tree = 0 $
else dbl_tree = 1

if (NOT dbl_tree) then begin
    tree = {lrefine:0l, $
            nodeType:0l, $
            processorNumber:0l, $
            gid:lonarr(nfaces+1+nchild), $
            coord:fltarr(ndim), $
            size:fltarr(ndim), $
            bndBox:fltarr(2,ndim)}
endif else begin
    tree = {lrefine:0l, $
            nodeType:0l, $
            processorNumber:0l, $
            gid:lonarr(nfaces+1+nchild), $
            coord:dblarr(ndim), $
            size:dblarr(ndim), $
            bndBox:dblarr(2,ndim)}
endelse

tree = replicate(tree,sim_params.total_blocks)

; save the coordinate information into the tree structure
for i = 0, ndim-1 do begin
    tree[*].coord[i] = reform(coord[i,*])
    tree[*].size[i] = reform(size[i,*])
    tree[*].bndBox[0,i] = reform(bnd_box[0,i,*])
    tree[*].bndBox[1,i] = reform(bnd_box[1,i,*])
endfor

undefine, coord
undefine, size
undefine, bnd_box

;----------------------------------------------------------------------------
; read in the tree information
;----------------------------------------------------------------------------
ngid = long(2^ndim + 1 + 2*ndim)

dataset = H5D_OPEN(file_identifier, "refine level")
refine_level = H5D_READ(dataset)
H5D_CLOSE, dataset

dataset = H5D_OPEN(file_identifier, "node type")
node_type = H5D_READ(dataset)
H5D_CLOSE, dataset

dataset = H5D_OPEN(file_identifier, "gid")
gid = H5D_READ(dataset)
H5D_CLOSE, dataset

file_contents = h5_parse(filename)

if ( (where(tag_names(file_contents) EQ "PROCESSOR_NUMBER"))[0] NE -1) then begin
    dataset = H5D_OPEN(file_identifier, "processor number")
    proc_number = H5D_READ(dataset)
    H5D_CLOSE, dataset
endif else begin
    proc_number = 0
endelse

tree[*].lrefine = refine_level
tree[*].nodeType = node_type
tree[*].processorNumber = proc_number

for i = 0, ngid-1 do begin
    tree[*].gid[i] = reform(gid[i,*])
endfor

undefine, refine_level
undefine, node_type
undefine, gid
undefine, proc_number


; read in the unknowns

if (NOT double) then begin

    if var LT 0 then begin
        
        unk = fltarr(nvar,params.nxb,params.nyb,params.nzb,params.totBlocks)

        for i = 0, nvar-1 do begin

            dataset = H5D_OPEN(file_identifier, unk_names[i])
            var = H5D_READ(dataset)
            H5D_CLOSE, dataset

            help, var
            unk[i,*,*,*,*] = float(var)

        endfor

        undefine, var
	var = -1

    endif else begin
        unk = fltarr(params.nxb,params.nyb,params.nzb,params.totBlocks)

        dataset = H5D_OPEN(file_identifier, unk_names[var])
        var = H5D_READ(dataset)
        H5D_CLOSE, dataset

        help, var

        unk = reform(temporary(float(var)), $
                     1,params.nxb,params.nyb,params.nzb,params.totBlocks)

    endelse

endif else begin

    if var LT 0 then begin

        unk = dblarr(nvar,params.nxb,params.nyb,params.nzb,params.totBlocks)

        for i = 0, nvar-1 do begin

            dataset = H5D_OPEN(file_identifier, unk_names[i])
            var = H5D_READ(dataset)
            H5D_CLOSE, dataset

            help, var
            unk[i,*,*,*,*] = var

        endfor

        undefine, var
	var = -1

    endif else begin
        unk = dblarr(params.nxb,params.nyb,params.nzb,params.totBlocks)

        dataset = H5D_OPEN(file_identifier, unk_names[var])
        var = H5D_READ(dataset)
        H5D_CLOSE, dataset

        help, var

        unk = reform(temporary(unk), $
                     1,params.nxb,params.nyb,params.nzb,params.totBlocks)


    endelse

endelse
; sometimes, a simulation exists that has only a single block.  IDL
; has this annoying habit or automagically dropping a dimension from
; an array if it's index is 1.  A lot of other routines want this
; index to be here, so restore it.
if (params.totBlocks EQ 1) then begin
    if (var EQ -1) then begin
        unk = reform(unk, nvar,params.nxb,params.nyb,params.nzb,1)
    endif else begin
        unk = reform(unk, 1,params.nxb,params.nyb,params.nzb,1)
    endelse
endif

; switch the block index to the front
unk = transpose(temporary(unk),[0,4,1,2,3])


; get the number of particles
numParticles = 0l

; sometimes we don't store the particle information -- check first
file_contents = h5_parse(filename)
if ( (where(tag_names(file_contents) EQ $
            "PARTICLE_TRACERS"))[0] NE -1) then begin

    dataset = H5D_OPEN(file_identifier, "particle tracers")
    datatype = H5D_GET_TYPE(dataset)

    idl_type = H5T_IDLTYPE(datatype, STRUCTURE=particles)

    particles = H5D_READ(dataset)

    H5D_CLOSE, dataset
    H5T_CLOSE, datatype

    numParticles = (size(particles))[1]
endif else begin
    numParticles = 0l
    particles = 1
endelse


H5F_CLOSE, file_identifier

endif ; itype = HDF5 

if(itype EQ 2) then begin
 
file_identifier = NCDF_OPEN(filename)

;------------------------------------------------------------------------------
; grab the list of variables, and make sure that our requested
; variable (if any) exists in the dataset
;------------------------------------------------------------------------------

gdata = ncdf_inquire(file_identifier)
if (strpos(filename, 'plt') GT 0) then begin
	unk_vars = gdata.nvars - 6
	unk_names = strarr(1,unk_vars)
	for i=6, gdata.nvars-1 DO BEGIN
		unkname = ncdf_varinq(file_identifier, i)
		unk_names[i-6] = unkname.name
		endfor
	endif else begin
	unk_vars = gdata.nvars - 8
        unk_names = strarr(1,unk_vars)
	for i=8, gdata.nvars-1 DO BEGIN
                unkname = ncdf_varinq(file_identifier, i)
                unk_names[i-8] = unkname.name
                endfor
endelse

unk_names = reform(unk_names)


; we can also now set the number of variables
nvar = (size(unk_names))[1]

if (var_name NE 'none') then begin
    var = (where(unk_names EQ var_name))[0]

    if (var EQ -1) then begin
        print, 'ERROR: requested variable not found in dataset'
    endif

    print, 'reading in only the variable ', var_name
endif else begin
    var = -1
endelse



; read in the simulation paramters
ncdf_attget, file_identifier, /GLOBAL, "total_blocks",total_blocks
ncdf_attget, file_identifier, /GLOBAL, "time", time
ncdf_attget, file_identifier, /GLOBAL, "nsteps", nsteps 
ncdf_attget, file_identifier, /GLOBAL, "timestep", timestep
ncdf_attget, file_identifier, /GLOBAL, "redshift", redshift 
ncdf_attget, file_identifier, /GLOBAL, "nxb" , nxb
ncdf_attget, file_identifier, /GLOBAL, "nyb" , nyb
ncdf_attget, file_identifier, /GLOBAL, "nzb" , nzb

sim_params = create_struct('total_blocks' , total_blocks, $
		            'nsteps', nsteps, 'redshift', redshift, $		
		            'time', time, 'timestep', timestep, $
			    'nxb', nxb, 'nyb', nyb, 'nzb', nzb)


; get the dimensionality
ncdf_diminq, file_identifier, 5, nndim, ndim

; figure out if we are dealing with corners
corners = 0

if (strpos(filename, 'crn_') GT 0) then begin
    corners = 1
endif

nfaces = 2*ndim
nchild = 2^ndim



; get the geometry
if (n_elements(geometry) eq 0) then begin
    ncdf_attget,file_identifier,/GLOBAL,"geometry_name",geometry
endif 

; ---- setup the structures to pass the data ----------------------------------
if (n_elements(double) EQ 0) then double = 0


if (NOT double) then begin
    tree = {lrefine:0l, $
            nodeType:0l, $
            gid:lonarr(nfaces+1+nchild), $
            coord:fltarr(ndim), $
            size:fltarr(ndim), $
            bndBox:fltarr(2,ndim)}
endif else begin
    tree = {lrefine:0l, $
            nodeType:0l, $
            gid:lonarr(nfaces+1+nchild), $
            coord:dblarr(ndim), $
            size:dblarr(ndim), $
            bndBox:dblarr(2,ndim)}
endelse

tree = replicate(tree,sim_params.total_blocks)


if (NOT double) then begin
    params = {totBlocks:sim_params.total_blocks, $
              corners:corners, $
              ndim:ndim, $
              nvar:nvar, $
              nxb:sim_params.nxb, $
              nyb:sim_params.nyb, $
              nzb:sim_params.nzb, $
              ntopx:1, $
              ntopy:1, $
              ntopz:1, $
              time:float(sim_params.time), $
              dt:float(sim_params.timestep), $
 redshift:1.0, $
              geometry:"unknown"}
endif else begin
    params = {totBlocks:sim_params.total_blocks, $
              corners:corners, $
              ndim:ndim, $
              nvar:nvar, $
              nxb:sim_params.nxb, $
              nyb:sim_params.nyb, $
              nzb:sim_params.nzb, $
              ntopx:1, $
              ntopy:1, $
              ntopz:1, $
              time:sim_params.time, $
              dt:sim_params.timestep, $
              redshift:1.d0, $
              geometry:"unknown"}
endelse

params.geometry = strcompress(geometry, /REMOVE_ALL)

; for some reason, the redshift field is not stored in the plotfiles
if ( (where(tag_names(sim_params) EQ "REDSHIFT"))[0] NE -1) then begin

    if (NOT double) then begin
        params.redshift = float(sim_params.redshift)
    endif else begin
        params.redshift = sim_params.redshift
    endelse

endif


;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
; We now have the dimension, total number of blocks, and the number of
; zones per block in each direction.  We can use this to read the tree
; information, coordinate information, and unknowns
;xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

;----------------------------------------------------------------------------
; read in the tree information
;----------------------------------------------------------------------------
ngid = long(2^ndim + 1 + 2*ndim)

ncdf_varget, file_identifier, "lrefine", refine_level
ncdf_varget, file_identifier, "nodetype", node_type
ncdf_varget, file_identifier, "gid", gid

tree[*].lrefine = refine_level
tree[*].nodeType = node_type

for i = 0, ngid-1 do begin
    tree[*].gid[i] = reform(gid[i,*])
endfor

undefine, refine_level
undefine, node_type
undefine, gid



;----------------------------------------------------------------------------
; read in the coordinate infomation
;----------------------------------------------------------------------------

ncdf_varget, file_identifier, "coordinates" , coord
ncdf_varget, file_identifier, "blocksize", size
ncdf_varget, file_identifier, "bndbox", bnd_box

; save the coordinate information into the tree structure
for i = 0, ndim-1 do begin
    tree[*].coord[i] = reform(coord[i,*])
    tree[*].size[i] = reform(size[i,*])
    tree[*].bndBox[0,i] = reform(bnd_box[0,i,*])
    tree[*].bndBox[1,i] = reform(bnd_box[1,i,*])
endfor

undefine, coord
undefine, size
undefine, bnd_box




; read in the unknowns

if (NOT double) then begin

    if var LT 0 then begin

        unk = fltarr(nvar,params.nxb,params.nyb,params.nzb,params.totBlocks)

        for i = 0, nvar-1 do begin

            ncdf_varget, file_identifier,unk_names[i],var		
            help, var
            unk[i,*,*,*,*] = float(var)

 endfor

        undefine, var

    endif else begin
        unk = fltarr(params.nxb,params.nyb,params.nzb,params.totBlocks)

 	ncdf_varget,file_identifier,unk_names[var],var

        help, var

        unk = reform(temporary(float(var)), $
                     1,params.nxb,params.nyb,params.nzb,params.totBlocks)

    endelse

endif else begin

    if var LT 0 then begin

        unk = dblarr(nvar,params.nxb,params.nyb,params.nzb,params.totBlocks)

        for i = 0, nvar-1 do begin
		
            ncdf_varget, file_identifier, unk_names[i], var		

            help, var
            unk[i,*,*,*,*] = var

        endfor

        undefine, var

    endif else begin
        unk = dblarr(params.nxb,params.nyb,params.nzb,params.totBlocks)

 	ncdf_varget,file_identifier,unk_names[var],var

        help, var

        unk = reform(temporary(unk), $
                     1,params.nxb,params.nyb,params.nzb,params.totBlocks)

    endelse

endelse

; sometimes, a simulation exists that has only a single block.  IDL
; has this annoying habit or automagically dropping a dimension from
; an array if it's index is 1.  A lot of other routines want this
; index to be here, so restore it.
if (params.totBlocks EQ 1) then begin
    if (var EQ -1) then begin
        unk = reform(unk, nvar,params.nxb,params.nyb,params.nzb,1)
    endif else begin
        unk = reform(unk, 1,params.nxb,params.nyb,params.nzb,1)
    endelse
endif

; switch the block index to the front
unk = transpose(temporary(unk),[0,4,1,2,3])

; for now set numParticles and particles to default values in pnetcdf case

    numParticles = 0l
    particles = 1

NCDF_CLOSE, file_identifier

endif ; itype = NCDF

;------------------------------------------------------------------------------
; compute the number of top level blocks in each direction
;------------------------------------------------------------------------------
top_blocks = where(tree[*].lrefine EQ 1)

ntopx = 1
ntopy = 1
ntopz = 1


; for 1d problems
case ndim of
    1: begin
        ntopx = (size(where(tree[top_blocks].bndBox[0,0] EQ $
                            min(tree[*].bndBox[0,0]))))[1]
    end

;for 2d problems
    2: begin

; find the number of level 1 blocks whose lower y coord is the
; bottom of the domain
        ntopx1 = (size(where(tree[top_blocks].bndBox[0,1] EQ $
                             min(tree[*].bndBox[0,1]))))[1]

; now find the number of top level blocks at the upper y coord
        ntopx2 = (size(where(tree[top_blocks].bndBox[1,1] EQ $
                             max(tree[*].bndBox[1,1]))))[1]

; take the max of these, so we can deal with L shaped domains
        ntopx = ntopx1 > ntopx2

; find the number of level 1 blocks whose min x coord is the
        ntopy = (size(where(tree[top_blocks].bndBox[0,0] EQ $
                            min(tree[*].bndBox[0,0]))))[1]
    end


;for 3d problems
    3: begin

; find the number of level 1 blocks whose lower y coord is the minimum
; y value and whose lower z coord is the minimum z value
        ntopx = (size(where(tree[top_blocks].bndBox[0,1] EQ $
                            min(tree[*].bndBox[0,1]) AND $
                            tree[top_blocks].bndBox[0,2] EQ $
                            min(tree[*].bndBox[0,2]))))[1]

; find the number of level 1 blocks whose lower x coord is the minimum
; x value and whose lower z coord is the minimum z value
        ntopy = (size(where(tree[top_blocks].bndBox[0,0] EQ $
                            min(tree[*].bndBox[0,0]) AND $
                            tree[top_blocks].bndBox[0,2] EQ $
                            min(tree[*].bndBox[0,2]))))[1]

; find the number of level 1 blocks whose lower x coord is the minimum
; x value and whose lower y coord is the minimum y value
        ntopz = (size(where(tree[top_blocks].bndBox[0,0] EQ $
                            min(tree[*].bndBox[0,0]) AND $
                            tree[top_blocks].bndBox[0,1] EQ $
                            min(tree[*].bndBox[0,1]))))[1]

    end
endcase


params.ntopx = ntopx
params.ntopy = ntopy
params.ntopz = ntopz



end

