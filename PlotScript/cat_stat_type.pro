
PRO cat_stat_type,doplot,icase
loadct,39
!p.charsize=1.5
!p.thick=4
lcol  = [0,80]


if (doplot eq 1) then begin
; current plotting device.
   mydevice = !D.NAME
   SET_PLOT, 'PS'
   DEVICE, FILENAME='/sps/lsst/users/rcecile/Fig/cat_stat_type.eps', /PORTRAIT,/COLOR,XSIZE=8.8,YSIZE=8.8,FONT_SIZE=4

endif

restore,'temp_ascii_new.sav'
;restore,'/sps/lsst/data/rcecile/Planck_BAO2/cat_stat_All1_z.save' ;;;
;avant

if (icase eq 0) then begin
   restore,'/sps/lsst/data/rcecile/Planck_BAO2/cat_stat_type_zOld.save'
   all = read_ascii('/sps/lsst/data/rcecile/Planck_BAO2_grids/SelFunc_All_ZONLY.txt' , template =  TEMP_SEL_FUNC)
   gold = read_ascii('/sps/lsst/data/rcecile/Planck_BAO2_grids/SelFunc_All_specz_nofz.txt' , template =  TEMP_SEL_FUNC)
endif
if (icase eq 1) then begin
   restore,'/sps/lsst/data/rcecile/Planck_BAO2/cat_stat_typelfZuccaAllFalse.save'
   all = read_ascii('/sps/lsst/data/rcecile/Planck_BAO2_grids/SelFunc_Zucca_ZONLY.txt' , template =  TEMP_SEL_FUNC)
   gold = read_ascii('/sps/lsst/data/rcecile/Planck_BAO2_grids/SelFunc_Zucca_specz_nofz.txt' , template =  TEMP_SEL_FUNC)
endif

nonly = dblarr(n_elements(z))
nall = nonly
for i=0,n_elements(z)-2 do begin
   ok = where(all.(0) ge z[i] and all.(0) lt z[i+1])
   nonly[i] = total( (all.(1))[ok] )
   nall[i] = total( (gold.(1)*all.(1))[ok] ) 
endfor

myformat='(E8.1)'
ok = where(z ge 0.2 and z le 2.45)
myout = string([total(nonly[ok]),total(nall[ok])],format=myformat)
print,total(nonly[ok]),total(nall[ok])

;window,icase

!p.multi=[0,1,2]
plot,z,alog10(smooth(nonly,5)),/xs,/ys,ytit='log10(Ngal)',yra=[5,9.25],xmar=[7,1],ymar=[0,2],xtickname=replicate(" ",10),xra=[0.2,2.45],psym=10;,xticklen=1
oplot,z,alog10(smooth(nall,5)),col=lcol[1],psym=10
what=['No magnitude cut',$
      'Golden sample']+', Ngal ='+myout
legend,what,col=[lcol[0],lcol[1]],line=0,box=1,/fill,/left,/bottom,charsize=1.75

what=['Early type    ',$
      'Late type     ',$
      'StarBurst type' ]

tot=z123[*,0]
for i=0,n_elements(z)-2 do tot[i] = z123[i,0]+z123[i,1]+z123[i,2]
plot,z,alog10([1,100]),xtit='Redshift',/xs,/ys,ytit='Fraction [%]',psym=10,yra=[0,99],xmar=[7,1],ymar=[4,0],/nodata,xra=[0.2,2.45];,xticklen=1
for i=0,2 do oplot,z,z123[*,i]/tot*100.,col=lcol[1],psym=-1*(i*2+2),th=2,symsi=2

legend,what,col=lcol[1],psym=indgen(3)*2+2,box=1,/fill,/right,/center,charsize=1.75,th=2;,symsi=2

;window,icase+10
;!p.multi=0
;plot,z,alog10([1,100]),xtit='Redshift',/xs,/ys,ytit='log10(Ngal)',yra=[1,1.1*max(alog10(z123[*,*]))],xmar=[7,1],ymar=[4,0],/nodata,xra=[0.2,2.45];,psym=10;,xticklen=1
;for i=0,2 do oplot,z,alog10(z123[*,i]),col=lcol[1],lin=2-i;,psym=10

;legend,what,col=lcol[1],line=[2-indgen(3)],box=1,/fill,/left,/bottom,charsize=1.75

if (doplot eq 1) then begin
   DEVICE, /CLOSE
   SET_PLOT, mydevice
endif

stop
if (icase eq 2) then    write_jpeg,'/sps/lsst/users/rcecile/BAO_InOut/check_typeRamos.jpg' ,tvrd(true=3),true=3
if (icase eq 1) then    write_jpeg,'/sps/lsst/users/rcecile/BAO_InOut/check_typeDahlenAll.jpg' ,tvrd(true=3),true=3
if (icase eq 0) then    write_jpeg,'/sps/lsst/users/rcecile/BAO_InOut/check_typeOld.jpg' ,tvrd(true=3),true=3
if (icase eq 4) then    write_jpeg,'/sps/lsst/users/rcecile/BAO_InOut/check_typeZucca.jpg' ,tvrd(true=3),true=3

end

