PRO prepa_hist_zs_zp

dir='/sps/lsst/data/rcecile/TJP_BAO/'
nslice =70

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

pinterquart = dblarr(nslice)
pbias = dblarr(nslice)
poutlier = dblarr(nslice)
dinterquart = dblarr(nslice)
dbias = dblarr(nslice)
doutlier = dblarr(nslice)
zstat = dblarr(nslice)

; better to do this part in qlogin mode
for i = 50,nslice-1 do begin
   file='cat_G2_Slice'+strtrim(i,2)+'_errP.fits'
   print,file
   m=mrdfits(dir+file,1,h,col=['ZS','ZP'])
;   h = hist_2d(m.(0),m.(1),bin1=binz,bin2=binz,max1=zmax,max2=zmax,min1=0,min2=0)
;   hp = hp + h

;   h = histogram( (m.(0)-m.(1))/(1.+m.(0)) ,bin=diff_bin,min=-1*diff_max,max=diff_max)
;   hpdiff = hpdiff + h
   zstat[i] = median(m.(0))
   t = (m.(1)-m.(0))/(1.+m.(0))
   RSTAT,t, Med, Hinge1, Hinge2
   pbias[i] = abs(Med)
   pinterquart[i] = Hinge2-Hinge1
   out = where(abs(t) gt 3.*pinterquart[i],nout)
   print,minmax(t),' et ',3.*pinterquart[i]
   poutlier[i] = double(nout)/double(n_elements(m.(0)))*100.
   print,'ICI ',nout,double(nout),double(n_elements(m.(0))),poutlier[i]

   file='cat_G2_Slice'+strtrim(i,2)+'_errPpodds.fits'
   print,file
   m=mrdfits(dir+file,1,h,col=['ZS','ZP'])
;   h = hist_2d(m.(0),m.(1),bin1=binz,bin2=binz,max1=zmax,max2=zmax,min1=0,min2=0)
;   hd = hd + h

;   h = histogram( (m.(0)-m.(1))/(1.+m.(0)) ,bin=diff_bin,min=-1*diff_max,max=diff_max)
;   hddiff = hddiff + h

   t = (m.(1)-m.(0))/(1.+m.(0))
   RSTAT,t, Med, Hinge1, Hinge2
   dbias[i] = abs(Med)
   dinterquart[i] = Hinge2-Hinge1
   out = where(abs(t) gt 3.*dinterquart[i],nout)
   doutlier[i] = double(nout)/double(n_elements(m.(0)))*100.

;   print,'Save up to slice ',i,total(hp),total(hd)
;   save,hp,hd,hpdiff,hddiff,z,binz,zmax,i,file=dir+'histo_zs_zp_0.01_diag.sav'
   print,'stat ',i,zstat[i],'       bias ',pbias[i],dbias[i],'    sigma ',pinterquart[i],dinterquart[i],'   outlier ',poutlier[i],doutlier[i]
   print,'----------------------------------------------------------------------------------------------------------------------'
   save,zstat,pinterquart,pbias,poutlier,dinterquart,dbias,doutlier,file=dir+'histo_zs_zp_requirement.sav'
endfor

END
