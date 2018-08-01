PRO print_psfit,ierr

h= 0.679
h3 = h*h*h

restore,'temp_ascii_new.sav'        ; contient dir
dir="/sps/lsst/data/rcecile/Planck_BAO2_PS/"

z=[0.5,0.9, 1.5]
namez = ['0.5','0.9','1.3']
nx= ['_120','_225', '_300']
nz = n_elements(namez)

nsuff=['','_err0.03','_errP','_errPBDT9','_errPBDT8']
lerr=['spectroZ','Gaussian 0.03','photoZ','photoZ BDT 90%','photoZ BDT 80%']
myformat ='(F7.1)'
kmin=[4,3,2]
kmax=[20,20,20,20,20,    20,15,15,15,15,    20,12,12,12,12]
nk05 = ['_k0.0' +strtrim(kmin[0],2)+'_0.'+strtrim(kmax[0],2),'_k0.0' +strtrim(kmin[0],2)+'_0.'+strtrim(kmax[1],2),'_k0.0' +strtrim(kmin[0],2)+'_0.'+strtrim(kmax[2],2),'_k0.0' +strtrim(kmin[0],2)+'_0.'+strtrim(kmax[3],2),'_k0.0' +strtrim(kmin[0],2)+'_0.'+strtrim(kmax[4],2)]
nk09 = ['_k0.0' +strtrim(kmin[1],2)+'_0.'+strtrim(kmax[5],2),'_k0.0' +strtrim(kmin[1],2)+'_0.'+strtrim(kmax[6],2),'_k0.0' +strtrim(kmin[1],2)+'_0.'+strtrim(kmax[7],2),'_k0.0' +strtrim(kmin[1],2)+'_0.'+strtrim(kmax[8],2),'_k0.0' +strtrim(kmin[1],2)+'_0.'+strtrim(kmax[9],2)]
nk15 = ['_k0.0' +strtrim(kmin[2],2)+'_0.'+strtrim(kmax[10],2),'_k0.0' +strtrim(kmin[2],2)+'_0.'+strtrim(kmax[11],2),'_k0.0' +strtrim(kmin[2],2)+'_0.'+strtrim(kmax[12],2),'_k0.0' +strtrim(kmin[2],2)+'_0.'+strtrim(kmax[13],2),'_k0.0' +strtrim(kmin[2],2)+'_0.'+strtrim(kmax[14],2)]
nk=[nk05,nk09,nk15]

mysuff = nsuff[ierr]
print,' '
print,lerr[ierr]
kmax = intarr(n_elements(namez), n_elements(nsuff))
kmax[*,0] = [0,1,3]
kmax[*,1] = [1,2,5]
kmax[*,2] = [1,2,5]
kmax[*,3] = [2,2,5]
kmax[*,4] = [2,2,5]


for iz=0,2 do begin
   print, '     iz = ', namez[iz],' for kmax =',kmax[iz,ierr]
   k = kmax[iz,ierr]
   ffit = dir + 'fit_lfZuccaAllFalse'+nk[iz*3+ierr]+nx[iz]+'_z'+namez[iz]+mysuff+'_result.txt'
   pfitw = read_ascii(ffit, template =  TEMP_POW_SPEC_FIT_RES)     
   
   print,'kmin/max= ',nk[iz*3+ierr],' s = ',(pfitw.(0))[0],' +',(pfitw.(1))[0],' - ',(pfitw.(2))[0],'], chi = ',(pfitw.(4))[0],'], amp = ',(pfitw.(3))[0]
   print,' ' 
endfor

end
