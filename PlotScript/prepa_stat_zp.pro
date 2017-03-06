PRO prepa_stat_zp,isuff

dir='/sps/lsst/data/rcecile/Planck_BAO/'

nslice = 100
suff = ['_errP','_errPBDT9','_errPBDT8']
pinterquart = dblarr(nslice)
pbias = dblarr(nslice)
poutlier = dblarr(nslice)
z_max = 3.

z = findgen(nslice+1)*z_max/nslice


; better to do this part in qlogin mode
for i = 3,nslice-2 do begin
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   file='cat_zOrdZP_Slice'+strtrim(i,2)+suff[isuff]+'.fits'
   print,file
   m=mrdfits(dir+file,1,h,col=['ZS','ZP'])
   t = (m.(1)-m.(0))/(1.+m.(0))
   
   RSTAT,t, Med, Hinge1, Hinge2
   pbias[i] = Med
   pinterquart[i] = Hinge2-Hinge1
   out = where(abs(t) gt 0.15,nout)
   poutlier[i] = double(nout)/double(n_elements(t))*100.

   print,'As a function of z_s --------------------------------------------------'
   print,'-------------------------------------------------------------------------------------------------'
   print,'stat ',i,'       bias ',pbias[i],'    sigma ',pinterquart[i],'   outlier ',poutlier[i]
   print,'-------------------------------------------------------------------------------------------------'
   save,z,pinterquart,pbias,poutlier,file=dir+'stat_zp_1pzs'+suff[isuff]+'.sav'
endfor

END
