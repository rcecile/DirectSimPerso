PRO test_zp

dir='/sps/lsst/data/rcecile/TJP_BAO/'
nslice =70
binz=0.01
zmax=3.
x=findgen(zmax/binz+1)*binz
nz=n_elements(x)
new = dindgen(nz)
for i = 0,nslice-1 do begin
   file='file'+strtrim(i,2)+'.fits'
   print,file
   m=mrdfits(dir+file,1,hh,col=['ZP','PODDS'])
   ok = where(m.(1) gt 0.894,nok)
   print,'N GAL ', nok
   h=histogram((m.(0))[ok],bin=binz,max=zmax,min=0)
   new = new +  h
endfor
old = dindgen(nz)
for i = 0,10 do begin
   file='Cat_G2_Slice'+strtrim(i,2)+'.fits'
   print,file
   m=mrdfits(dir+file,1,hh,col=['ZP','ZS'])
   ok = where(m.(1) le 0.31137100,nok)
   print,'N GAL ', nok
   old = old +  histogram((m.(0))[ok],bin=binz,max=zmax,min=0)
endfor

END
