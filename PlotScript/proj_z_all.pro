PRO proj_z_all,z,tot,nd,nErr,nz,snx,ngrids

suff='cube_'
err=['','_G0.03','_errP','_errPBDT9','_errPBDT8']
print,err

dirn='/sps/lsst/data/rcecile/Planck_noBAO_grids/'
ngrids=10
dir='/sps/lsst/data/rcecile/Planck_BAO_grids/'

D=['0.5','0.9','1.5']
snx=['160','300','225']
nz=[140,125,150]

nD = n_elements(D)
nErr = n_elements(err)

z = dblarr(nD,max(nz))
for id = 0,nD-1 do begin
   grid=dir+'grid_'+suff+strtrim(snx[id],2)+'_z'+D[id]+'.fits'
   print,grid
   c=mrdfits(grid,6,h)
   if (n_elements(c) gt 1) then begin
     ; print,h
      z[id,0:nz[id]-1]=c[0,0,*]
      ngal  = sxpar(h,'NGRID')
      print,'grid_'+suff+strtrim(snx[id],2)+'_z'+D[id],ngal
   endif
   c=0
endfor

tot = dblarr(ngrids+1,nD,nErr,max(nz))
help,tot

for in=0, ngrids do begin
   for id = 0,nD-1 do begin
      for ie = 0,nErr-1 do begin
         if (in eq 0) then mydir=dir else mydir = dirn
         if (in eq 0) then mysuff=suff else mysuff = suff + strtrim(in-1,2)+'_'

         grid=mydir+'grid_'+mysuff+strtrim(snx[id],2)+'_z'+D[id]+err[ie]+'.fits'
         print,grid
         c=mrdfits(grid,1,h)
         if (n_elements(c) gt 1) then begin
            for iz=0,nz[id]-1 do tot[in,id,ie,iz]=total(c[*,*,iz])
            if (ie gt 0) then begin
               h = headfits(grid,ext=2)
               ngal  = sxpar(h,'NGRID')
               print,'grid_'+strtrim(snx[id],2)+'_z'+D[id]+err[ie],ngal
            endif
         endif
         c=0
      endfor
   endfor

endfor



END
