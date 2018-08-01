
PRO stat_basique

dir='/sps/lsst/data/rcecile/Planck_BAO2/'
!p.charsize=2
!p.thick=3

h0=dblarr(70)
h1=h0
h3=h0
h4=h0
h0z=h0
h1z=h0
h3z=h0
h4z=h0

for i=0,69 do begin
   hh = headfits(dir+'cat_All_Slice'+strtrim(i,2)+'_ZONLY.fits',ext=1)
   h0[i] = 1.d * sxpar(hh,'NAXIS2')
   hh = headfits(dir+'cat_ELSB_Slice'+strtrim(i,2)+'_ZONLY.fits',ext=1)
   h1[i] = 1.d * sxpar(hh,'NAXIS2')
   hh = headfits(dir+'cat_LFdahlen_mag24_13all_Slice'+strtrim(i,2)+'_ZONLY.fits',ext=1)
   h3[i] = 1.d * sxpar(hh,'NAXIS2')
   hh = headfits(dir+'cat_LFdahlen_mag24_13_Slice'+strtrim(i,2)+'_ZONLY.fits',ext=1)
   h4[i] = 1.d * sxpar(hh,'NAXIS2')

   hh = headfits(dir+'cat_All_Slice'+strtrim(i,2)+'.fits',ext=1)
   h0z[i] = 1.d * sxpar(hh,'NAXIS2')
   hh = headfits(dir+'cat_ELSB_Slice'+strtrim(i,2)+'.fits',ext=1)
   h1z[i] = 1.d * sxpar(hh,'NAXIS2')
   hh = headfits(dir+'cat_LFdahlen_mag24_13all_Slice'+strtrim(i,2)+'.fits',ext=1)
   h3z[i] = 1.d * sxpar(hh,'NAXIS2')
   hh = headfits(dir+'cat_LFdahlen_mag24_13_Slice'+strtrim(i,2)+'.fits',ext=1)
   h4z[i] = 1.d * sxpar(hh,'NAXIS2')

   print,i
endfor

   plot,h0,/xs,/ys,yra=[1e5,max([max(h),max(h1)])*1.1],xtit="tranche",ytit="Ngal avant et apres golden cut",/yl,/nodata
   oplot,x,h0,col=123
   oplot,x,h1,col=23
   oplot,x,h3,col=123,li=2
   oplot,x,h4,col=23,li=2
   oplot,x,h0z,col=123
   oplot,x,h1z,col=23
   oplot,x,h3z,col=123,li=2
   oplot,x,h4z,col=23,li=2
   legend,['old all=true','old all=false' ,'new all=true','new all=false'],col=[123,23,123,23],li=[0,0,2,2],/fill,/left,/top


;write_jpeg,'/sps/lsst/dev/rcecile/BAO_InOut/check_ngal_zonly.jpg' ,tvrd(true=3),true=3
stop
END
