loadct,39
!p.charsize=2.
!p.thick=3
lsuff=['lfZuccaAllFalse']
isuff=0
suff = lsuff[isuff]

dir='/sps/lsst/data/rcecile/Planck_BAO2/'
is=8

name=dir+"cat_"+suff+"_Slice"+strtrim(is,2)+".fits"
m=mrdfits(name,1,h)

ok1 = where(fix(m.(4)) eq 1,n1)
ok2 = where(fix(m.(4)) eq 2,n2)
ok3 = where(fix(m.(4)) eq 3,n3)
print,n1,n2,n3,n1+n2+n3
print,mean(m.(3)),100.*n1/(n1+n2+n3),100.*n2/(n1+n2+n3),100.*n3/(n1+n2+n3)


h = histogram((m.(5)),min=-24,max=-16,bin=0.1)
h1 = histogram((m.(5))[ok1],min=-24,max=-16,bin=0.1)
h2 = histogram((m.(5))[ok2],min=-24,max=-16,bin=0.1)
h3 = histogram((m.(5))[ok3],min=-24,max=-16,bin=0.1)
x=findgen(n_elements(h1))*0.1-24.
mycol=[250,200,80]
window,10
plot ,x,h,/xs,/ys,yra=[1e2,2e6],/yl,xra=[-24,-16.5]
oplot ,x,h1,col=mycol[0]
oplot,x,h2,col=mycol[1]
oplot,x,h3,col=mycol[2]
   legend,['All','Early = type 1 Zucca','Late = type 2 Zucca','StarBurst = type 3+4 Zucca'],line=0,col=[0,mycol],box=1,/fill,/right,/bottom,charsize=1.5

;   write_jpeg,'/sps/lsst/users/rcecile/BAO_InOut/check_MagAbs.jpg' ,tvrd(true=3),true=3

z=mrdfits("test_Zucca.fits",1)
d=mrdfits("test_Dahlen.fits",1)
r1=mrdfits("test_Ramos1.fits",1)
r2=mrdfits("test_Ramos2.fits",1)
plot,z.(0),z.(8),/xs,/ys,th=4
oplot,z.(0),z.(5),psym=1
oplot,z.(0),z.(6),psym=4
oplot,z.(0),z.(7),psym=-2
oplot,d.(0),d.(8),col=123,th=4
oplot,z.(0),d.(5),psym=1,col=123
oplot,z.(0),d.(6),psym=4,col=123
oplot,z.(0),zd(7),psym=-2,col=123

surf = !pi / !pi * 180. / !pi * 180. *60. * 60. ; 1/4 du ciel en arc-min2
zmin = [0.2,0.5,1.0,1.5,0.2]
zmax = [0.5,1.0,1.5,2.5,2.5]
n_z = n_elements(zmin)
.r
for iz=0,n_z-1 do begin
   ok = where (z.(0) ge zmin[iz] and z.(0) lt zmax[iz],nok)
   print, zmin[iz],zmax[iz],'  ',total((z.(8))[ok]),total((d.(8))[ok]),total((r1.(8))[ok]),total((r2.(8))[ok])
endfor
end
                                   Zucca           Dahlen
     0.200000     0.500000         14.231554       15.619235
     0.500000      1.00000         30.156607       20.769728
      1.00000      1.50000         12.389560       6.6521939
      1.50000      2.50000         3.6481854       2.2359551

     0.200000      2.50000         60.425907       45.277113
