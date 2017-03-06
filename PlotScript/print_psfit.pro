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
skmax = ['','_k0.20','_k0.18','_k0.16','_k0.14','_k0.12','_k0.10','_k0.08']
myformat ='(F7.1)'

mysuff = nsuff[ierr]
print,' '
print,lerr[ierr]

for iz=0,2 do begin
   print, '     iz = ', namez[iz]
   for k= 0,n_elements(skmax)-1 do begin
      ffit = dir + 'fit_newErr'+skmax[k]+nx[iz]+'_z'+namez[iz]+mysuff+'_result.txt'
      pfit = read_ascii(ffit, template =  TEMP_POW_SPEC_FIT)     
      
      kpt = string(pfit.(1)*1000.,format=myformat)+' +/-'+string(pfit.(3)*1000.,format=myformat)
      spt = string(2.*!pi/pfit.(1),format=myformat)+' +/-'+string(2.*!pi*pfit.(3)/pfit.(1)/pfit.(1),format=myformat)

      chi2fit = dir + 'fit_newErr'+skmax[k]+nx[iz]+'_z'+namez[iz]+mysuff+'_chisq.txt'
      c2fit = read_ascii(chi2fit, template =  TEMP_POW_SPEC_FITCHI2)     
                                ;  print, '         kmax = ', skmax[k]+ '   k_a ='+kpt + '      chi2 = ' + strtrim(min(c2fit.(1)),2)  
      print, '         s = ', skmax[k]+ '   s ='+spt + '      chi2 = ' + string(min(c2fit.(1)),format=myformat)  
      
   endfor
endfor


end
