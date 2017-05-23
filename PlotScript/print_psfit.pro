PRO print_psfit,ierr

h= 0.679
h3 = h*h*h

restore,'temp_ascii.sav'        ; contient dir
dir="/sps/lsst/data/rcecile/Planck_BAO_PS/"

z=[0.5,0.9, 1.5]
namez = ['0.5','0.9','1.5']
nx= ['_120','_225', '_320']
nz = n_elements(namez)

nsuff=['','_err0.03','_errP','_errPBDT9','_errPBDT8']
lerr=['spectroZ','Gaussian 0.03','photoZ','photoZ BDT 90%','photoZ BDT 80%']
skmax = ['_k0.20','_k0.18','_k0.16','_k0.14','_k0.12','_k0.10']
myformat ='(F7.1)'

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
   ffit = dir + 'fit_Ref'+skmax[k]+nx[iz]+'_z'+namez[iz]+mysuff+'_result.txt'
   pfitw = read_ascii(ffit, template =  TEMP_POW_SPEC_FIT_SA)     
   
   print,'kmax= ',skmax[k],' s = ',(pfitw.(1))[0],' +',(pfitw.(3))[0],' - ',(pfitw.(4))[0],'], chi = ',(pfitw.(6))[0],'], amp = ',(pfitw.(5))[0]
   print,' ' 
endfor


end
