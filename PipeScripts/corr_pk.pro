PRO corr_pk,iz
; bias : corrige le z ref
; sigma : damp le spectre
; outlier : shot noise addtionnel


!p.charsize=2
!p.thick=2
!p.symsize=4
h= 0.679
h3 = h*h*h

loadct,39
restore,'temp_ascii.sav'        ; contient dir
dir="/sps/lsst/data/rcecile/TJP_BAO_PS/"
z=[0.9, 1.3,1.8,1.8]
namez = ['0.9','1.3','1.8','1.8']
nx= ['_640', '_900', '_1024', '_500']
gsuff = ['', '', ' thin grid', ' thick grid']
cell=[8,8,8,16]
kmax = [12, 12, 12, 6]
 ;       z      biasx1000       sigmax100      %outlier
;     0.900000      0.79657436       2.4108059       2.5204028
;      1.30000     -0.30511629       2.6356436       2.1756421
;      1.80000      -7.0835856       3.8389667       3.0180539
;      1.80000      -5.8074459       3.2689693       4.4343087
bias1000 = [0.79657436,-0.30511629,-7.0835856,-5.8074459]
sigma100 = [2.4108059,2.6356436,3.8389667,3.2689693]
outlier = [2.5204028,2.1756421,3.0180539,4.4343087]
   
; read data
t  =  dir + 'PS_G2_2fitErrk'+nx[iz]+'_z'+namez[iz]+mysuff+'_wngal.txt'
pcorr = read_ascii(t, template =  TEMP_2FIT)     
xs = pcorr.(0)
   
; computed biased reference spectrum
zb = bias1000[iz]/1000.*(1+z[iz])+z[iz]
command="/sps/lsst/dev/rcecile/BAOProgs/SimLSS/exe/cmvginit3d -o /sps/lsst/data/rcecile/TJP_BAO_PS/simu_psBIASED -u 0.7,0.3,0.05,0.7,-1.,0.,9.06581172724113E-05 -8 0.8154,-8. -G 2 -s 1000,0.001,0.5 -x 10,10 -y 10,10 -z10,10 -Z "+ strtrim(zb,2)+" -O 2,2"
print,command
spawn,command

ftheo  = dir + 'simu_psBIASED_ntpk.txt'
ptheo = read_ascii(ftheo, template =  TEMP_POW_SPEC_TH)     
xt=ptheo.(0)
ok = intarr(n_elements(xs))
for i=0,n_elements(xs)-1 do ok[i] = where( abs(xs[i]-xt) eq min(abs(xs[i]-xt)))

; 

;igma
;grille de random avec 1 gal/cell, cell 1Mpc. calculer spectre avec damping selon kz (tous les sigma/pas de 0.01  de 1 à 5). Verifier (in)dépendance / taille de la grille

;cell size
;grille avec random de 1/cell de 1Mpc. Calculer spectre avec grille de cell =8 ou 16 Mpc et 1Mpc.

;outlier

;somme (outliers / selfunc) * <selfunc>_deltaZ = f_out (ca doit pas être bon)
;SN = SN*sqrt(f_out) ????????

;prendre le Pk_theo

;calculer Pk_degr = Pk_theo(zb) * Pk_cell / Pk_1Mpc * Pk_sigma/Pk_sigma0 * (1. + SN_out)
;puis osc = (Pk_obs-SNobs) /  Pk_degr



END
