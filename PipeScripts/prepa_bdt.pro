PRO prepa_bdt,first

n=70

ngal_per_slice = 10000ll

dir='/sps/lsst/data/rcecile/TJP_BAO/'
dir2='/sps/lsst/data/rcecile/TJP_BAO_BDT/'

nkeep = 0
nbdt = 0l
for i=first,n-1 do begin
     name = dir + "cat_G2_Slice"+strtrim(i,2)+".fits"
     m=mrdfits(name,1,hm)
     ngal = sxpar(hm,'NAXIS2')
     nslice = (ngal + nkeep) / ngal_per_slice
     print,'slice ',i,name,nslice
     for j=0,nslice-1 do begin    
        if (nbdt gt 36999) then stop
        name2 = dir2 + "cat_Slice"+strtrim(nbdt,2)+".fits"

        if (j eq 0 and nkeep gt 0) then mym = [mkeep, m[0:ngal_per_slice-nkeep-1]]  $
        else mym = m[j*ngal_per_slice -nkeep : (j+1)*ngal_per_slice-1 -nkeep]

        if (j mod 100 eq 0) then $
           if (j eq 0 and nkeep gt 0) then $
              print,j,name2,n_elements(mym),0,nkeep,(j+1)*ngal_per_slice-1 -nkeep $
           else $
              print,j,name2,n_elements(mym),j*ngal_per_slice -nkeep,(j+1)*ngal_per_slice-1 -nkeep 
        h=hm
        if (nbdt eq 36999) then mwrfits,mym,name2, h
        sxaddpar,h,'ZCAT_MIN',min(mym.(3))
        sxaddpar,h,'ZCAT_MAX',max(mym.(3))
        if (nbdt eq 36999) then modfits,name2,0,h,exten_no=1 ;Update header of ext 1, let data unchaged 

        nbdt ++
     endfor
     mkeep = m[j*ngal_per_slice  -nkeep : *]
     nkeep = n_elements(mkeep)
      
endfor


END