PRO cata_concat,simu
; cata_concat,''
; cata_concat,'0'
; for i=0,9 do cata_concat,strtrim(i,2)

dir='/sps/lsst/data/rcecile/'
if (strlen(simu) eq 0) then dir = dir+"Planck_BAO/" else  dir = dir+"Planck_noBAO/"
if (strlen(simu) eq 0) then cas = "" else  cas = "_"+simu
print,dir
nsuff = ['_errP_bin','_errPBDT9_bin','_errPBDT8_bin'];,'_err0.03']



loadct,39
!p.charsize=2
!p.thick=3
;recopie une tranche sauf les catastrophiques (> 3 sigma de G 0.03), regroupes
i_min=3
i_max=99
     
for isuff=1,2 do begin
   suff = nsuff[isuff]
   nok_tot= lonarr(i_max-i_min+1)
; do get the reasons of the cases where we have no podds value
   n10_tot = nok_tot
   n20_tot = nok_tot
   n50_tot = nok_tot
   n60_tot =  nok_tot
   n_tot = nok_tot
   z_tot = dblarr(i_max-i_min+1)
   
   ifirst=0
   for i=i_min,i_max do begin
      
      name=dir+"cat_zOrd"+cas+"_Slice"+strtrim(i,2)+suff+"_cata.fits"
      check = FILE_TEST(name)
      print,name
      if (check eq 0) then continue
      print,'  ok '
      m=mrdfits(name,1,h) 
      ok = where(m.(6) lt 10,nok)
      if (ifirst eq 0) then mtot = m[ok] else mtot=[mtot,m[ok]]
      nok_tot[i-i_min] = nok
      print,name,n_elements(m.(6)),nok
      n_tot[i-i_min] = n_elements(m.(6))
      z_tot[i-i_min] = mean(m.(3))
      
      if (strcmp(suff,'_errPpodds') eq 1) then begin
         pb10 = where(m.(6) eq 10,n10)
         pb20 = where(m.(6) eq 20,n20)
         pb50 = where(m.(6) eq 50,n50)
         pb60 = where(m.(6) eq 60,n60)
         n10_tot[i-i_min] = n10
         n20_tot[i-i_min] = n20
         n50_tot[i-i_min] = n50
         n60_tot[i-i_min] = n60
         print,name,'                  et  >= 10 : ',n10,n20,n50,n60
      endif
      ifirst=1
      
   endfor
   
   plot,z_tot,n_tot,/xs,/ys,th=3
   oplot,z_tot,nok_tot,th=3,col=123
   oplot,z_tot,n10_tot+n20_tot+n50_tot+n60_tot,th=3,col=80
   
   nametot=dir+"cat_gold"+cas+"_AllSlice"+suff+"_cataLT10.fits"
   MWRFITS,mtot, nametot, h
   
;   save,z_tot,nok_tot,n_tot,n10_tot,n20_tot,n50_tot,n60_tot,file=dir+'/StatCata'+cas+'_'+suff+'.sav'
   
endfor

END
