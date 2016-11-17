PRO cut_cube,do_cut,suff
; cut_cube,0,'' to see the slices z
; cut_cube,1,'' to cut the cube with BAO
; cut_cube,1,'0' to cut the cube 0 with noBAO
; for i=0,9 do  cut_cube,1,strtrim(i,2)

Nz_tot = 700 
dlref = 3200
n=70
cell=8

dir='/sps/lsst/data/rcecile/'
;if (strlen(suff) eq 0) then dir = dir+"TJP_BAO/" else  dir = dir+"TJP_noBAO/"
;if (strlen(suff) eq 0) then dir = dir+"Planck_BAO/" else  dir = dir+"Planck_noBAO/"
if (strlen(suff) eq 0) then dir = dir+"toy_cube/" 
print,dir

if (do_cut eq 1) then begin
   file=dir+"simu"+suff+"_r.fits"
   print,file
   c=mrdfits(file,0,hh)
   help,c
   sxaddpar,hh,'DREF',dlref
   print,"dlref read, =", dlref
endif

dlmin = dlref - Nz_tot/2 * cell
dlmax= Nz_tot * cell + dlmin

Nz = Nz_tot / n
print,'Nz = ',Nz

h    = indgen(n)* Nz * cell + dlmin + Nz/2*cell
hmin = h - Nz / 2 * cell
hmax = h + Nz / 2 * cell
hdbl = 1.d*h
zref=dblarr(n)
href=zref 

for i=0,n-1 do zref[i] = zfrlos(hdbl[i],7)
href = hdbl
 
print,'============================ INFO PAVES : ============================'

for i=0,n-1 do print,'slice ',i,i*Nz,hdbl[i],zref[i],href[i]


if do_cut eq 1 then begin
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; DO CUT ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 
   print,hh
   
   for i=0,n-1 do begin
     print,'slice ',i,i*Nz,hdbl[i],zref[i],href[i]
     print,i*Nz , i*Nz + Nz-1
     help,c
     my_c = c[i*Nz : i*Nz + Nz-1, *, *]
     WRITEFITS,dir+"simu70"+suff+"_Slice"+strtrim(i,2)+"_r.fits", my_c,hh
     sxaddpar,hh,'NZ',Nz
     sxaddpar,hh,'ZREF',zref[i]
     sxaddpar,hh,'KZREF',1.*Nz/2.
     sxaddpar,hh,'HREF',href[i]
     sxaddpar,hh,'DREF',hdbl[i]
     modfits,dir+"simu70"+suff+"_Slice"+strtrim(i,2)+"_r.fits",0,hh    
   endfor
      
endif


END
