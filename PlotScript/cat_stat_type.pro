PRO cat_stat_type2,doplot
;cat_stat_type,0 : to compute the histogram (long), must be done first
;cat_stat_type,1 : to do the plot (fast)
loadct,39
!p.charsize=2.
!p.thick=4
lcol  = [0,80]

if (doplot eq 0) then begin

   dir='/sps/lsst/data/rcecile/Planck_BAO2/'
   nSlice=100
   
   n_type2=3
   ntype2=dblarr(nSlice,n_type2)
   z2=dblarr(nSlice)

  
   for i=3,nSlice-1 do begin
      name=dir+"cat_zType_Slice"+strtrim(i,2)+".fits"
      m=mrdfits(name,1,h,col=['TYPE'])
      zmin = sxpar(h,'ZCAT_MIN')
      zmax = sxpar(h,'ZCAT_MAX')
      z2[i] = (zmin+zmax)/2.
      for it=0,n_type2-1 do begin
         ok = where(fix(m.(0)) eq it+1,nok)
         ntype2[i,it] = double(nok)
      endfor
      PRINT,NAME,ntype2[i,0],ntype2[i,1],ntype2[i,2]
      m=0
   endfor

   mytot2=dblarr(n_type2)
   for i=0,n_type2-1 do mytot2[i] = total(ntype2[*,i])
   myformat='(E8.2)'
   myout2 = string(mytot2,format=myformat)

   nall2=dblarr(nSlice)
   for i=0,nSlice-1 do nall2[i] = total(ntype2[i,*])

   save,z2,ntype2,n_type2,nall2,myout2,file='/sps/lsst/data/rcecile/Planck_BAO2/cat_stat_ztype.save'
   stop
    
endif


if (doplot ge 1) then restore,'/sps/lsst/data/rcecile/Planck_BAO/cat_stat_type.save'
if (doplot ge 1) then restore,'/sps/lsst/data/rcecile/Planck_BAO2/cat_stat_ztype.save'
restore,'temp_ascii.sav'
all = read_ascii('/sps/lsst/data/rcecile/Planck_BAO_grids/SelFunc_gold__ZONLY.txt' , template =  TEMP_SEL_FUNC)
gold = read_ascii('/sps/lsst/data/rcecile/Planck_BAO_grids/SelFunc_gold__specz_nofz.txt' , template =  TEMP_SEL_FUNC)

all2 = read_ascii('/sps/lsst/data/rcecile/Planck_BAO2_grids/SelFunc_type__ZONLY.txt' , template =  TEMP_SEL_FUNC)
gold2 = read_ascii('/sps/lsst/data/rcecile/Planck_BAO2_grids/SelFunc_type__specz_nofz.txt' , template =  TEMP_SEL_FUNC)


nonly = dblarr(n_elements(z))
nall = nonly
for i=0,n_elements(z)-2 do begin
   ok = where(all.(0) ge z[i] and all.(0) lt z[i+1])
   nonly[i] = total( (all.(1))[ok] )
   nall[i] = total( (gold.(1)*all.(1))[ok] ) 
endfor
nonly2 = dblarr(n_elements(z2))
nall2 = nonly2
for i=0,n_elements(z2)-2 do begin
   ok = where(all2.(0) ge z2[i] and all2.(0) lt z2[i+1])
   nonly2[i] = total( (all2.(1))[ok] )
   nall2[i] = total( (gold2.(1)*all2.(1))[ok] ) 
endfor

myformat='(E8.1)'
ok = where(z ge 0.2 and z le 2.45)
myout = string([total(nonly[ok]),total(nall[ok])],format=myformat)
myout2 = string([total(nonly2[ok]),total(nall2[ok])],format=myformat)

!p.multi=[0,1,2]


;window,0
plot,z,alog10(smooth(nonly,5)),/xs,/ys,ytit='log10(Ngal)',yra=[5,8.9],xmar=[7,1],ymar=[0,2],xtickname=replicate(" ",10),xra=[0.2,2.45];,xticklen=1
oplot,z,alog10(smooth(nall,5)),col=lcol[1]
oplot,z,alog10(smooth(nonly2,5)),col=150,li=2
oplot,z,alog10(smooth(nall2,5)),col=lcol[1]+150,li=2

what=['No magnitude cut',$
      'Golden sample']+', Ngal ='+myout
legend,what,col=[lcol[0],lcol[1]],line=0,box=1,/fill,/left,/bottom,charsize=1.75
what=['No magnitude cut',$
      'Golden sample']+', Ngal ='+myout2
legend,what,col=[lcol[0],lcol[1]]+150,line=2,box=1,/fill,/right,/center,charsize=1.75

plot,z,alog10([1,150]),xtit='Redshift',/xs,/ys,ytit='log10(Ngal)',psym=10,yra=[2,7.9],xmar=[7,1],ymar=[4,0],/nodata,xra=[0.2,2.45];,xticklen=1
for i=0,n_type-1 do oplot,z,alog10(ntype[*,i]),col=lcol[1],lin=2-i
for i=0,n_type-1 do oplot,z,alog10(ntype2[*,i]),col=lcol[1]+150,lin=2-i

what=['Early type    ',$
      'Late type     ',$
      'StarBurst type' ]
legend,what,col=lcol[1]+150,line=[2-indgen(n_type)],box=1,/fill,/left,/bottom,charsize=1.75

;write_jpeg,'/sps/lsst/dev/rcecile/BAO_InOut/check_ngal.jpg' ,tvrd(true=3),true=3

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

surf = !pi / !pi * 180. / !pi * 180. *60. * 60.; 1/4 du ciel en arc-min2 + changement de bining

!p.multi=0

plot,z, (smooth(nall,5))/surf /0.03*0.05,/xs,/ys,/yl,yra=[0.04,2],xra=[0.2,2.45],xtit='redshift',ytit='Ngal/arc-min2/bin en z 0.05',th=5
oplot,z,(smooth(nall2,5))/surf/0.03*0.05,col=150,th=5
for i=0,n_type-1 do oplot,z,(ntype[*,i])/surf/0.03*0.05,lin=5-i*2
for i=0,n_type-1 do oplot,z,(ntype2[*,i])/surf/0.03*0.05,col=150,lin=5-i*2

what=['Early type    ',$
      'Late type     ',$
      'StarBurst type' ]
legend,what,col=0,line=[5-indgen(n_type)*2],box=1,/fill,/right,/top,charsize=1.75

write_jpeg,'/sps/lsst/dev/rcecile/BAO_InOut/check_ngal_FH.jpg' ,tvrd(true=3),true=3

stop
end
