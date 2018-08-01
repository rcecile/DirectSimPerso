PRO print_test_lf
nsuff=['All','LF_Dahlen','LF_DahlenTEST_type','LF_DahlenAll','LF_DahlenAll_ini']
n_prod = n_elements(nsuff)
!p.charsize=1.75

loadct,39
restore,'temp_ascii_new.sav'        ; contient dir
dir="/sps/lsst/data/rcecile/Planck_BAO2/"

h=dblarr(n_prod)
for i=10,60,10 do begin
   for ip = 0,n_prod-1 do begin
      f="/sps/lsst/data/rcecile/Planck_BAO2/cat_"+nsuff[ip]+"_Slice"+strtrim(i,2)+".fits"
      hh=headfits(f,ext=1,/silent)
      h[ip] = 1.d * sxpar(hh,'NAXIS2')
   endfor
   z0 = 1.d * sxpar(hh,'ZCAT_MIN')+1.d
   z1 = 1.d * sxpar(hh,'ZCAT_MAX')+1.d
   print,'Slice ',i,z0,'-',z1,' : ',h[0],h[1],h[2],h[3],h[4]
endfor


END
