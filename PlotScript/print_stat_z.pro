PRO print_stat_z

surf = !pi / !pi * 180. / !pi * 180. *60. * 60. ; 1/4 du ciel en arc-min2
suff = ['_err0.03','_errP','_errPBDT9','_errPBDT8']
dir='/sps/lsst/data/rcecile/Planck_BAO2/'
theta = 60. /180. * !pi

;;;;;;;; TABLE 2
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
print,'!!!!!!!!!!! TABLE 2 !!!!!!!!!!!'

;restore,dir+'histo_zs_zp'+suff[0]+'.sav'
restore,dir+'histo_zs_statZp'+suff[0]+'.sav'

n = n_elements(z)
nz = dblarr(n)
n_mpc = nz

d = dblarr(n,n)
for i=0,n-1 do d[i,*] = (z[i] - z)/ (1.+z[i])
sidx1 = where( abs(d) lt 0.02)
sidx3 = where( abs(d) lt 0.06)
sidx5 = where( abs(d) lt 0.05)
sidx9 = where( abs(d) lt 0.15)

myformat = '(F7.1, F7.1, F7.1, F7.1)'

print,'[%] within 0.02, 0.06, 0.05, 0.15 for z_p EN FONCTION DE ZS ou ZP'
for i=0,n_elements(suff)-1 do begin
   ;restore,dir+'histo_zs_zp'+suff[i]+'.sav'
   restore,dir+'histo_zs_statZp'+suff[i]+'.sav'

   if (i eq 0) then all_hist_ref = all_hist
   tot=total(all_hist)
   tot_ref=total(all_hist_ref)
   print,' '
   print,suff[i]
   n = [total(all_hist[sidx1]),total(all_hist[sidx3]),total(all_hist[sidx5]),total(all_hist[sidx9])]
   sn_zs = string(n/tot*100.,format=myformat)
   print,'fct de Z_S ',sn_zs

endfor

read,xx
;;;;;;;; TABLE 1

print,'!!!!!!!!!!! TABLE 1 !!!!!!!!!!!'
nz = ['z_s','z_p','z_{BDT9}','z_{BDT8}']
zmin = [0.2,0.5,1.0,1.5,0.2]
zmax = [0.5,1.0,1.5,2.5,2.5]
n_z = n_elements(zmin)
n = dblarr(n_z)

myformat = '(A10, F7.1,F7.1, F7.1, F7.1, F17.1)'

tit= '['+strtrim(zmin,2)+'-'+strtrim(zmax,2)+'['
print,tit
for i=0,n_elements(suff)-1 do begin

   restore,dir+'histo_zs_statZp'+suff[i]+'.sav'

   for iz=0,n_z-1 do begin
      ok = where (z ge zmin[iz] and z lt zmax[iz])
      if (i eq 0) then  n[iz] = total(all_hist[ok,*]) else n[iz] = total(all_hist[*,ok])
   endfor
   res = string(nz[i],n/surf,format=myformat)
   print,res
endfor



END
