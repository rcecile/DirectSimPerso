PRO plot_cube_grids,doplot

if (doplot) then begin
; current plotting device.
   mydevice = !D.NAME
   SET_PLOT, 'PS'
   DEVICE, FILENAME='/sps/lsst/dev/rcecile/Fig/plot_cube_2grids.eps',/PORTRAIT,/COLOR,XSIZE=8.8,YSIZE=9.9,FONT_SIZE=4
endif 

open_ang = !pi/3.

loadct,39
!p.charsize=2.
!x.ticklen=0
!y.ticklen=0
cell_z=8.
Nz=700.
Nx=1600.
Nx=1250.
h0=3200.
xc=findgen(Nx)*cell_z - Nx/2*cell_z
y0=fltarr(Nx)

hmin = h0-Nz/2.*cell_z
hmax = h0+Nz/2.*cell_z
xmax = Nx * cell_z /2.
xmin = -1*xmax
theta=findgen(1000)/1000.*2*open_ang -open_ang
mnt = min(theta)
mxt = max(theta)

na=1000

plot,[0,0],[0,0],/xs,/ys,th=3,xra=[0,5300],yra=[0,6100],/nodata,xtit='Distance [comoving Mpc]',ytit='Distance [comoving Mpc]',xmar=[8,1]

nslice=70
h= h0 + (findgen(nslice+1)-nslice/2)*Nz*cell_z/nslice
ang=Nx*cell_z / h / sqrt(2.)
for i=0,nslice-1 do oplot,[xmin,xmax],[h[i],h[i]],th=2,col=190
for i=0,nslice-1 do begin 
   x=findgen(na)/na*ang[i] -ang[i]/2.
;   oplot,x,y0+h[i],th=2,col=190
;   oplot,h[i]*sin(x),h[i]*cos(x),th=2,col=190
endfor
ang=Nx*cell_z / h 
for i=0,nslice-1 do oplot,[xmin,xmax],[h[i],h[i]],th=2,col=190
for i=0,nslice-1 do begin 
   x=findgen(na)/na*ang[i] -ang[i]/2.
;   oplot,x,y0+h[i],th=2,col=190
;   oplot,h[i]*sin(x),h[i]*cos(x),th=2,col=190,li=2
endfor
oplot,[0,cos(mnt+!pi/2.)*hmax*2],[0,sin(mnt+!pi/2.)*hmax*2],th=3
oplot,[0,cos(mxt+!pi/2.)*hmax*2],[0.,sin(mxt+!pi/2.)*hmax*2],th=3

zc=[0.7,1.4]
nzc=[200,200]

zc=[0.9,1.6]
nzc = [125.,150.] ;Nz du cube produit par cat_grid

zc=[0.9,1.3,1.8,1.8]
nzc = [125.,75.,65.,75.] ;Nz du cube produit par cat_grid
cell_g=[8, 8, 8, 16]

zc=[0.5,0.9,1.3,1.8]
nzc = [140., 125.,75.,65.] ;Nz du cube produit par cat_grid
cell_g=[8, 8, 8, 8]

ng = n_elements(zc)
hg=dloscom(zc)
h1=(dloscom(zc)-nzc/2.*cell_g)
h2=(dloscom(zc)+nzc/2.*cell_g)
for i=0,ng-1 do print,'redshift min/max = ',zfrlos(h1[i],10),zfrlos(h2[i],10)

Nx_parfait = h1* 2*open_ang/8.;/sqrt(2)
print,'Nx_parfait =',Nx_parfait

nxc=[350,640,900,1024]
nxc=[350,640,900,1200]
print,'Nx used =',nxc
print,'Hgrids = ',hg

coeff=[1.,sqrt(2.)]
mycol=[50,95,210,250]

for ia=0,1 do begin
   mynxc = nxc * coeff[ia]
   ang1=myNxc*cell_g / h1 / sqrt(2.)
   ang2=myNxc*cell_g / h2 / sqrt(2.)
   for i=0,ng-1 do begin
      x1=findgen(na)/na*ang1[i] -ang1[i]/2.
      oplot,cos(x1+!pi/2.)*h1[i],sin(x1+!pi/2.)*h1[i],col=mycol[i],th=5,li=ia*2
     
      x2=findgen(na)/na*ang2[i] -ang2[i]/2.
      oplot,cos(x2+!pi/2.)*h2[i],sin(x2+!pi/2.)*h2[i],col=mycol[i],th=5,li=ia*2
     
      oplot,[cos(x1[0]+!pi/2.)*h1[i],cos(x2[0]+!pi/2.)*h2[i]],[sin(x1[0]+!pi/2.)*h1[i],sin(x2[0]+!pi/2.)*h2[i]],col=mycol[i],th=5,li=ia*2
      oplot,[cos(x1[na-1]+!pi/2.)*h1[i],cos(x2[na-1]+!pi/2.)*h2[i]],[sin(x1[na-1]+!pi/2.)*h1[i],sin(x2[na-1]+!pi/2.)*h2[i]],col=mycol[i],th=5,li=ia*2

   endfor
endfor
pos = intarr(ng)+300
xyouts,pos,hg,'z='+string(zc,format='(f4.2)')
;xyouts,1800,700,'Comoving radial distance',ori=28

print,'Volume (Gpc^3) = ',1.d*cell_g^3*nxc*nxc*nzc/1.e9
if (doplot) then begin
   DEVICE, /CLOSE
   SET_PLOT, mydevice
endif

stop



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
window,1,xs=600,ys=600
!p.charsize=2
!x.ticklen=0
!y.ticklen=0
phi=findgen(1000)/1000.*2.*!pi
r = h2*open_ang
plot,cos(phi)*r[ng-1],sin(phi)*r[ng-1],th=5,/nodata,xtit='Transverse distance [comoving Mpc]',ytit='Transverse distance [comoving Mpc]',yma=[6,2]
for i=ng-1,0,-1 do oplot,cos(phi)*r[i],sin(phi)*r[i],th=5,col=(i+1)*80,li=1
r = h1*open_ang
for i=ng-1,0,-1 do oplot,cos(phi)*r[i],sin(phi)*r[i],th=5,col=(i+1)*80
nxc *= cell_z
for i=0,ng-1 do begin
   oplot,[-nxc[i]/2,nxc[i]/2],[-nxc[i]/2,-nxc[i]/2],col=(i+1)*80,th=5
   oplot,[-nxc[i]/2,nxc[i]/2],[nxc[i]/2,nxc[i]/2],col=(i+1)*80,th=5
   oplot,[-nxc[i]/2,-nxc[i]/2],[-nxc[i]/2,nxc[i]/2],col=(i+1)*80,th=5
   oplot,[nxc[i]/2,nxc[i]/2],[-nxc[i]/2,nxc[i]/2],col=(i+1)*80,th=5 
endfor

print,nxc^2 /4./!pi/h1^2
print,4.*!pi*h1^2
print,nxc^2


stop

END

