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

plot,[0,0],[0,0],/xs,/ys,th=3,xra=[0,5300],yra=[0,6100],/nodata,xtit='Distance along the x Euclidian axis [comoving Mpc]',ytit='Distance along the z Euclidian axis [comoving Mpc]',xmar=[8,1],ymar=[3.5,.5]

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


;zc=[0.5,0.9,1.3,1.8]
;nzc = [140., 125.,75.,65.] ;Nz du cube produit par cat_grid
;cell_g=[8, 8, 8, 8]

;zc=[0.5,0.9,1.5]
;nzc = [125., 125.,175.] ;Nz du cube produit par cat_grid
zc=[0.5,0.9,1.3]
nzc = [125., 125.,125.] ;Nz du cube produit par cat_grid
cell_g=[8, 8, 8, 8]

ng = n_elements(zc)
hg=dloscom(zc)
h1=(dloscom(zc)-nzc/2.*cell_g)
h2=(dloscom(zc)+nzc/2.*cell_g)
for i=0,ng-1 do print,'redshift min/max = ',zfrlos(h1[i],10),zfrlos(h2[i],10)

Nx_parfait = h1* 2*open_ang/8.
Nx_parfait = h1*open_ang / cell_z
print,'Nx_parfait =',Nx_parfait

nxc=[160,300,225]
nxc=[120,225,320]
nxc=[120,225,300]
print,'Nx used =',nxc
print,'Hgrids = ',hg

mycol=[80,205,240]
xx=findgen(na)/na*open_ang*2 -open_ang

ang1=Nxc/2.*cell_g / h1       
ang2=Nxc/2.*cell_g / h2  
ang_rot = 40. *!pi/180.
cosT = cos(ang_rot)
sinT = sin(ang_rot)
      
for i=0,ng-1 do begin

   x = findgen(nxc[i]/2.)*cell_g[i]
   xtot = findgen(nxc[i])*cell_g[i] - nxc[i]/2.*cell_g[i]
   y1 = x*0. + h1[i]
   y2 = x*0. + h2[i]
   y1tot = xtot*0. + h1[i]
   y2tot = xtot*0. + h2[i]

   y = findgen(nzc[i])*cell_g[i] -nzc[i]/2.*cell_g[i] + dloscom(zc[i])
   x1 = y*0. -nxc[i]/2.*cell_g[i]
   x2 = y*0. +nxc[i]/2.*cell_g[i]

   oplot,x,y1,th=5,col=mycol[i]
   oplot,x,y2,th=5,col=mycol[i]
   oplot,x2,y,th=5,col=mycol[i]

; trace 2ieme grille tournée de 40 degrés
   oplot,xtot*cosT+y1tot*sinT,-xtot*sinT + y1tot*cosT,th=5,col=mycol[i]
   oplot,xtot*cosT+y2tot*sinT,-xtot*sinT + y2tot*cosT,th=5,col=mycol[i]
   oplot,x1*cosT+y*sinT,-x1*sinT + y*cosT,th=5,col=mycol[i]
   oplot,x2*cosT+y*sinT,-x2*sinT + y*cosT,th=5,col=mycol[i]
   
   oplot,cos(xx+!pi/2.)*dloscom(zc[i]),sin(xx+!pi/2.)*dloscom(zc[i]),li=2,col=mycol[i],th=2
   xyouts,x2[0]*cosT+y[0]*sinT,-x2[0]*sinT + y[0]*cosT-150,'redshift='+string(zc[i],format='(f4.2)'),col=mycol[i],charth=2
   xyouts,x2[0]*cosT+y[0]*sinT,-x2[0]*sinT + y[0]*cosT-150,'redshift='+string(zc[i],format='(f4.2)')

endfor

oplot,cos(xx+!pi/2.)*dloscom(0.2),sin(xx+!pi/2.)*dloscom(0.2),li=2,th=2

oplot,cos(xx+!pi/2.)*dloscom(2.45),sin(xx+!pi/2.)*dloscom(2.45),li=2,th=2
;oplot,cos(xx+!pi/2.)*dloscom(.5),sin(xx+!pi/2.)*dloscom(.5),li=2
;oplot,cos(xx+!pi/2.)*dloscom(.75),sin(xx+!pi/2.)*dloscom(.75),li=2


xyouts,700,300,'redshift=0.2'
xyouts,4000,4200,'redshift=2.45'
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

print,dloscom(0.2)
print,dloscom(2.45)

stop

END

