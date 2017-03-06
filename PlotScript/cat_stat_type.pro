PRO cat_stat_type,doplot
;cat_stat_type,0 : to compute the histogram (long), must be done first
;cat_stat_type,1 : to do the plot (fast)
loadct,39
!p.charsize=2.
!p.thick=4
lcol  = [0,80]

if (doplot eq 0) then begin

   dir='/sps/lsst/data/rcecile/Planck_BAO/'
   nSlice=100
   
   n_type=3
   ntype=dblarr(nSlice,n_type)
   z=dblarr(nSlice)

  
   for i=3,nSlice-1 do begin
      name=dir+"cat_zOrd_Slice"+strtrim(i,2)+".fits"
      m=mrdfits(name,1,h,col=['TYPE'])
      zmin = sxpar(h,'ZCAT_MIN')
      zmax = sxpar(h,'ZCAT_MAX')
      z[i] = (zmin+zmax)/2.
      for it=0,n_type-1 do begin
         ok = where(fix(m.(0)) eq it+1,nok)
         ntype[i,it] = double(nok)
      endfor
      PRINT,NAME,ntype[i,*]
      m=0
   endfor

   mytot=dblarr(n_type)
   for i=0,n_type-1 do mytot[i] = total(ntype[*,i])
   myformat='(E8.2)'
   myout = string(mytot,format=myformat)

   nall=dblarr(nSlice)
   for i=0,nSlice-1 do nall[i] = total(ntype[i,*])

   save,z,ntype,n_type,nall,myout,file='/sps/lsst/data/rcecile/Planck_BAO/cat_stat_type.save'
   stop
    
endif


if (doplot eq 2) then begin
; current plotting device.
   mydevice = !D.NAME
   SET_PLOT, 'PS'
   DEVICE, FILENAME='/sps/lsst/dev/rcecile/Fig/cat_stat_type.eps', /PORTRAIT,/COLOR,XSIZE=8.8,YSIZE=8.8,FONT_SIZE=4

endif

if (doplot ge 1) then restore,'/sps/lsst/data/rcecile/Planck_BAO/cat_stat_type.save'
restore,'temp_ascii.sav'
all = read_ascii('/sps/lsst/data/rcecile/Planck_BAO_grids/SelFunc_gold__ZONLY.txt' , template =  TEMP_SEL_FUNC)
gold = read_ascii('/sps/lsst/data/rcecile/Planck_BAO_grids/SelFunc_gold__specz_nofz.txt' , template =  TEMP_SEL_FUNC)


nonly = dblarr(n_elements(z))
nall = nonly
for i=0,n_elements(z)-2 do begin
   ok = where(all.(0) ge z[i] and all.(0) lt z[i+1])
   nonly[i] = total( (all.(1))[ok] )
   nall[i] = total( (gold.(1)*all.(1))[ok] ) 
endfor

;h=dloscom(zall)
;ang=1600.*8. / h / sqrt(2.)
;correct = where(ang lt 120./180*!pi)
;nonlyc= nonly
;nallc = nall
;nonlyc[correct] = nonlyc[correct] / ang[correct]*120./180*!pi
;nallc[correct] = nallc[correct] / ang[correct]*120./180*!pi

myformat='(E8.1)'
ok = where(z ge 0.2 and z le 2.45)
myout = string([total(nonly[ok]),total(nall[ok])],format=myformat)

!p.multi=[0,1,2]
plot,z,alog10(smooth(nonly,5)),/xs,/ys,ytit='log10(Ngal)',yra=[5,8.9],xmar=[7,1],ymar=[0,2],xtickname=replicate(" ",10),xra=[0.2,2.45];,xticklen=1
oplot,z,alog10(smooth(nall,5)),col=lcol[1]
;oplot,z,alog10(nonlyc),li=2
;oplot,z,alog10(nallc),col=lcol[1],li=2
what=['No magnitude cut',$
      'Golden sample']+', Ngal ='+myout
legend,what,col=[lcol[0],lcol[1]],line=0,box=1,/fill,/left,/bottom,charsize=1.75

plot,z,alog10([1,100]),xtit='Redshift',/xs,/ys,ytit='log10(Ngal)',psym=10,yra=[2,7.9],xmar=[7,1],ymar=[4,0],/nodata,xra=[0.2,2.45];,xticklen=1
for i=0,n_type-1 do oplot,z,alog10(ntype[*,i]),col=lcol[1],lin=2-i
what=['Early type    ',$
      'Late type     ',$
      'StarBurst type' ]

legend,what,col=lcol[1],line=[2-indgen(n_type)],box=1,/fill,/left,/bottom,charsize=1.75

if (doplot eq 2) then begin
   DEVICE, /CLOSE
   SET_PLOT, mydevice
endif

stop

end

