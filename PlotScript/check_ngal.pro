PRO check_ngal


restore,'temp_ascii_new.sav'
!p.charsize=1.5
!p.multi=0
!p.thick=3

dirg='/sps/lsst/data/rcecile/Planck_BAO_grids/'
dirg2='/sps/lsst/data/rcecile/Planck_BAO2_grids/'

fs=dirg+'SelFunc_gold__specz_nofz.txt'
fs2=dirg2+'SelFunc_gold__specz_nofz.txt'
fs3=dirg2+'SelFunc_type__specz_nofz.txt'

s = read_ascii(fs, template =  TEMP_SEL_FUNC)
s2 = read_ascii(fs2, template =  TEMP_SEL_FUNC)
s3 = read_ascii(fs3, template =  TEMP_SEL_FUNC)

!p.multi=[0,1,2]
plot,s.(0),smooth(s.(1),5),/xs,/ys,ytit='Selection fct',xra=[0.2,2.45],yra=[0.0003,1],/yl
oplot,s2.(0),smooth(s2.(1),5),col=123
oplot,s3.(0),smooth(s3.(1),5),col=234,li=2
legend,['old','new','new+debug'],col=[0,123,234],li=[0,0,2],/fill,/left,/bottom
plot,s.(0),smooth(s2.(1)/s.(1),5),xra=[0.2,2.45],yra=[0.9,2],ytit='SF new/old'
oplot,s.(0),smooth(s2.(1)/s.(1),5),col=123
oplot,s.(0),smooth(s3.(1)/s.(1),5),col=234

;write_jpeg,'/sps/lsst/dev/rcecile/BAO_InOut/check_ngal_sf_debug.jpg' ,tvrd(true=3),true=3

sr = read_ascii(dirg+'SelFunc_gold__ZONLY.txt', template =  TEMP_SEL_FUNC)
sr2 = read_ascii(dirg2+'SelFunc_gold__ZONLY.txt', template =  TEMP_SEL_FUNC)
sr3 = read_ascii(dirg2+'SelFunc_type__ZONLY.txt', template =  TEMP_SEL_FUNC)
h=sr.(1) *s.(1)
h2=sr2.(1) *s2.(1)
h3=sr3.(1) *s3.(1)

!p.multi=0
plot,s.(0),h,/xs,/ys,/yl,xra=[0.2,2.45],yra=[1e4,.5e7]
oplot,s.(0),h2,col=123
oplot,s.(0),h3,col=234,li=2
legend,['old','new','new+debug'],col=[0,123,234],li=[0,0,2],/fill,/right,/top


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
surf = !pi / !pi * 180. / !pi * 180. *60. * 60. ; 1/4 du ciel en arc-min2
z=dindgen(51)/50.*2.5
ok = where(s.(0) ge 0.2 and s.(0) lt 2.5)

plot,s.(0),h/surf /0.002*0.05,/xs,/ys,/yl,yra=[0.04,2],xra=[0.2,2.45],xtit='redshift',ytit='Ngal/arc-min2/bin en z 0.05'
oplot,s.(0),h2/surf /0.002*0.05,col=123
oplot,s.(0),h3/surf /0.002*0.05,col=234,li=2
;write_jpeg,'/sps/lsst/dev/rcecile/BAO_InOut/check_ngal_FH.jpg' ,tvrd(true=3),true=3
legend,['old','new','new+debug'],col=[0,123,234],li=[0,0,2],/fill,/right,/top

zmin=[0.36,0.72,1.16]
zmax=[0.68,1.19,2.12]

for i=0,n_elements(zmin)-1 do begin
   ok = where(s.(0) ge zmin[i] and s.(0) le zmax[i])
   t=total(h[ok])/surf
   t2=total(h2[ok])/surf
   print,'z range =', zmin[i],zmax[i],', Ngal/arc-min2 = (avant, apres, diff[%]',t,t2, (t2-t)/t*100.
endfor
print,' '
print,' '


zmin=[0.2,0.5,1.,1.5,0.2]
zmax=[0.5,1.,1.5,2.5,2.5]

for i=0,n_elements(zmin)-1 do begin
   ok = where(s.(0) ge zmin[i] and s.(0) le zmax[i])
   t=total(h[ok])/surf
   t2=total(h2[ok])/surf
   print,'z range =', zmin[i],zmax[i],', Ngal/arc-min2 = (avant, apres, diff[%]',t,t2, (t2-t)/t*100.
endfor

print,'************ AVANT *********************'
for i=10,29,4 do begin
;   m=mrdfits("/sps/lsst/data/rcecile/Planck_BAO2/cat_zType_Slice"+strtrim(i,2)+".fits",1,range=[0,1e6])
   m=mrdfits("/sps/lsst/data/rcecile/Planck_BAO/cat_zOrd_Slice"+strtrim(i,2)+".fits",1,range=[0,1e6])
   ok=where(m.(4) eq 2.06)
   ok2=where(m.(4) eq 2.15)
   ok3=where(m.(4) eq 2.34)
   print,'z range = ' ,minmax(m.(3)), '   n type 6 = ', n_elements(ok), '   n type 15 = ', n_elements(ok2), '   n type 34 = ', n_elements(ok3)
endfor

print,'************ TYPE *********************'
for i=10,29,4 do begin
;   m=mrdfits("/sps/lsst/data/rcecile/Planck_BAO2/cat_zType_Slice"+strtrim(i,2)+".fits",1,range=[0,1e6])
   m=mrdfits("/sps/lsst/data/rcecile/Planck_BAO2/cat_zOrd_Slice"+strtrim(i,2)+".fits",1,range=[0,1e6])
   
   ok=where(m.(4) eq 2.05)
   ok2=where(m.(4) eq 2.06)
   ok3=where(m.(4) eq 2.34)
   print,'z range = ' ,minmax(m.(3)), '   n type 6 = ', n_elements(ok), '   n type 15 = ', n_elements(ok2), '   n type 34 = ', n_elements(ok3)
endfor

print,'************ TYPE sans 0.5 *********************'
for i=10,29,4 do begin
   m=mrdfits("/sps/lsst/data/rcecile/Planck_BAO2/cat_zType_Slice"+strtrim(i,2)+".fits",1,range=[0,1e6])
   ok=where(m.(4) eq 2.05)
   ok2=where(m.(4) eq 2.06)
   ok3=where(m.(4) eq 2.34)
   print,'z range = ' ,minmax(m.(3)), '   n type 6 = ', n_elements(ok), '   n type 15 = ', n_elements(ok2), '   n type 34 = ', n_elements(ok3)
endfor

stop

;m=mrdfits("/sps/lsst/data/rcecile/Planck_BAO2/cat_goldALL_Slice10.fits",1,col=["ZS","TYPE"],range=[0,1e7])
;z=dindgen(11)/10.*.5 + 0.25

;m=mrdfits("/sps/lsst/data/rcecile/Planck_BAO2/cat_goldALL_Slice60.fits",1,col=["ZS","TYPE"],range=[10e7,11e7])
;z=dindgen(101)/100.*.5 + 1.75

m=mrdfits("/sps/lsst/data/rcecile/Planck_BAO2/cat_zType_Slice25.fits",1,col=["ZS","TYPE"])
z=dindgen(101)/100.*.5 + 1.

mz = (m.(0))
mt = (m.(1))

a1=0
b1=4                            ;   //early (1)
a2=5
b2=34                           ;   //late  (2)
a3=35
b3=50                           ;   //starburst (3)

type=fltarr(51)
for i=a1,b1 do type[i] = 1.+ i/100.
for i=a2,b2 do type[i] = 2.+ i/100.
for i=a3,b3 do type[i] = 3.+ i/100.

n = dblarr(n_elements(z),51)

for in=0,50 do begin
   for i=0,n_elements(z)-2 do begin
      n[i,in] = 1.d*n_elements(where(mz ge z[i] and mz lt z[i+1] and mt eq type[in]))
   endfor
endfor

n1=a1
for n2=a1+1,b1 do begin
   ok = where(n[*,n2] gt 1,nok)
   if (nok lt 10) then continue
   plot,z,n[ok,n1]/n[ok,n2],/xs,/ys,psym=10
   PRINT,n1,n2,MEAN(n[ok,n1]/n[ok,n2])
endfor
n1=a2
for n2=a2+1,b2 do begin
   ok = where(n[*,n2] gt 1,nok)
   if (nok lt 10) then continue
   plot,z,n[ok,n1]/n[ok,n2],/xs,/ys,psym=10
   PRINT,n1,n2,MEAN(n[ok,n1]/n[ok,n2])
endfor
n1=a3
for n2=a3+1,b3 do begin
   ok = where(n[*,n2] gt 1,nok)
   if (nok lt 10) then print,'rien en type=',n2
   if (nok lt 10) then continue
   plot,z,n[ok,n1]/n[ok,n2],/xs,/ys,psym=10
   PRINT,n1,n2,MEAN(n[ok,n1]/n[ok,n2])
endfor

;write_jpeg,'/sps/lsst/dev/rcecile/BAO_InOut/check_ngal_n6n34_slice60.jpg' ,tvrd(true=3),true=3
;  int a1=0,  b1=5;    //early (1)
;  int a2=5,  b2=35;   //late  (2)
;  int a3=35, b3=51;   //starburst (3)

!p.multi=[0,1,3]
m=mrdfits("/sps/lsst/data/rcecile/Planck_BAO2/cat_type_Slice12.fits",1,range=[0,1e6])
m=mrdfits("/sps/lsst/data/rcecile/Planck_BAO2/cat_zType_Slice3.fits",1,range=[0,1e6])
h=histogram(m.(4),min=0,max=3.6,bin=0.01)
xh = findgen(n_elements(h))/(n_elements(h)-1)*3.6
plot,xh,h,/xs,/ys,psym=10,xra=[0.98,1.06]
plot,xh,h,/xs,/ys,psym=10,xra=[2.03,2.36]
plot,xh,h,/xs,/ys,psym=10,xra=[3.33,3.53]
stop

END
