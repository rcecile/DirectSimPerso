PRO prepa_hist_zs_zp

dir='/sps/lsst/data/rcecile/TJP_BAO/'

binz=0.01
zmax=3.
ncell = fix(zmax/binz+1)
z=findgen(ncell)*binz

hp = dblarr(ncell,ncell)
hd = hp

diff_bin=0.001
diff_max=3.
hpdiff = dblarr(fix(diff_max/diff_bin*2+1))
hddiff = hpdiff

; better to do this part in qlogin mode
for i = 0,69 do begin
   file='cat_G2_Slice'+strtrim(i,2)+'_errP.fits'
   print,file
   m=mrdfits(dir+file,1,h,col=['ZS','ZP'])
   h = hist_2d(m.(0),m.(1),bin1=binz,bin2=binz,max1=zmax,max2=zmax,min1=0,min2=0)
   hp = hp + h

   h = histogram( (m.(0)-m.(1))/(1.+m.(0)) ,bin=diff_bin,min=-1*diff_max,max=diff_max)
   hpdiff = hpdiff + h

   file='cat_G2_Slice'+strtrim(i,2)+'_errPpodds.fits'
   print,file
   m=mrdfits(dir+file,1,h,col=['ZS','ZP'])
   h = hist_2d(m.(0),m.(1),bin1=binz,bin2=binz,max1=zmax,max2=zmax,min1=0,min2=0)
   hd = hd + h

   h = histogram( (m.(0)-m.(1))/(1.+m.(0)) ,bin=diff_bin,min=-1*diff_max,max=diff_max)
   hddiff = hddiff + h

   print,'Save up to slice ',i,total(hp),total(hd)
   save,hp,hd,hpdiff,hddiff,z,binz,zmax,i,file=dir+'histo_zs_zp_0.01_diag.sav'
endfor

END
