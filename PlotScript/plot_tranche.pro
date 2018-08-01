PRO plot_tranche,n,i
; pour tracer la tranche i dans une grille de n tranches
; plot_tranche,128,0 par ex

t="/sps/lsst/data/benaadad/Planck_BAO/cat_z_08_1.txt"

;TEMP_pourDTFE = ASCII_TEMPLATE(t)
;save,TEMP_pourDTFE,file='temp_merieme.sav'

restore,'temp_merieme.sav'
   
g = read_ascii(t, template =  TEMP_pourDTFE)     
x = g.(0)
y = g.(1)
z = g.(2)

dz = (max(z)-min(z))/n
my_zmin = min(z) + i * dz 
my_zmax = min(z) + (i+1) * dz 
ok = where(z ge my_zmin and z lt my_zmax,nok)
x = x[ok]
y = y[ok]
print,nok,minmax(z),dz,minmax(z[ok])


window,i,XSIZE=600,YSIZE=600
plot,x,y,/xs,/ys,psym=1

stop

END
 
