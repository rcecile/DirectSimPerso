
PRO proj_z,z,tot,nd,nErr,nz,snx

suff='_G2sfmean'
err=['','_G0.03','_errP','_errPpodds']
print,err

;dir='/sps/lsst/data/rcecile/TJP_noBAO_grids/'
;suff='G2_0_'
dir='/sps/lsst/data/rcecile/TJP_BAO_grids/'
suff='G2_'

D=['0.7','1.4']
snx=['450', '875']
nz=[200,200]

D=['0.9','1.3','1.8']
snx=['640','900','500']
nz=[125,75,75]


nD = n_elements(D)
nErr = n_elements(err)

z = dblarr(nD,max(nz))
for id = 0,nD-1 do begin
   grid=dir+'grid_'+suff+strtrim(snx[id],2)+'_z'+D[id]+'.fits'
   print,grid
   c=mrdfits(grid,3,h)
   if (n_elements(c) gt 1) then begin
     ; print,h
      z[id,0:nz[id]-1]=c[0,0,*]
      ngal  = sxpar(h,'NGRID')
      print,'grid_'+suff+strtrim(snx[id],2)+'_z'+D[id],ngal
   endif
endfor

tot = dblarr(nD,nErr,max(nz))
help,tot

for id = 0,nD-1 do begin
   for ie = 0,nErr-1 do begin
      grid=dir+'grid_'+suff+strtrim(snx[id],2)+'_z'+D[id]+err[ie]+'.fits'
      print,grid
      c=mrdfits(grid,1,h)
      if (n_elements(c) gt 1) then begin
         for iz=0,nz[id]-1 do tot[id,ie,iz]=total(c[*,*,iz])
         if (ie gt 0) then begin
            h = headfits(grid,ext=2)
            ngal  = sxpar(h,'NGRID')
            print,'grid_'+strtrim(snx[id],2)+'_z'+D[id]+err[ie],ngal
         endif
      endif
   endfor

endfor



END
