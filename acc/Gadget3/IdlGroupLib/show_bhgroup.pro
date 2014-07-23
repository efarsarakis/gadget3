
BaseDir = "/afs/rzg/bc-b/deboras/S1_csfm"   ; The output-directory of the simulation
SnapBase= "snap_S1"                  ; The basename of the snapshot files

BoxSize  = 10000.0
Num = 11                          ; The number of the snapshot we look at


;;; The object-file of the compile C-library for accessing the group catalogue

ObjectFile = "idlgrouplib.so"
Num=long(Num)
exts='000'
exts=exts+strcompress(string(num),/remove_all)
exts=strmid(exts,strlen(exts)-3,3)
Outputdir  = Basedir + "/"





;;;;;;;;;; First, we get the number of groups

Ngroups = CALL_EXTERNAL(ObjectFile, $
                       'get_total_number_of_groups', /UL_VALUE, $
                        OutputDir, $
                        Num)

print, "Number of groups in the catalogue: ", Ngroups

;;;;;;;;;; Now we load the group catalogue


GroupLen = lonarr(Ngroups)
GroupFileNr = lonarr(Ngroups)
GroupNr = lonarr(Ngroups)
GroupTypeLen = lonarr(6, Ngroups)
GroupTypeMass = dblarr(6, Ngroups)
GroupCenter = fltarr(3, Ngroups)
GroupSfr = fltarr(Ngroups)

Ngroups = CALL_EXTERNAL(ObjectFile, $
                       'get_group_catalogue', /UL_VALUE, $
                        OutputDir, $
                        Num, $
                        GroupLen, $
                        GroupFileNr, $
                        GroupNr, $
                        GroupTypeLen, $
                        GroupTypeMass, $
                        GroupCenter, $
			GroupSfr)

;;;;;;;;;; We determine the total particle number, and the number of files

  
Files = 0L

NumPart = CALL_EXTERNAL(ObjectFile, $
                       'get_total_particle_count', /UL_VALUE, $
                        OutputDir, $
                        Num, $
			SnapBase, $
                        Files)

;;;;;;;;;; We determine the gas particle number

  
NumGas = CALL_EXTERNAL(ObjectFile, $
                       'get_gas_particle_count', /UL_VALUE, $
                        OutputDir, $
                        Num, $
			SnapBase, $
                        Files)


;;;;;;;;;; We determine the dm particle number

  
NumDm = CALL_EXTERNAL(ObjectFile, $
                       'get_dm_particle_count', /UL_VALUE, $
                        OutputDir, $
                        Num, $
			SnapBase, $
                        Files)

;;;;;;;;;; We determine the particle number with variable mass

  
NumMass = CALL_EXTERNAL(ObjectFile, $
                       'get_mass_particle_count', /UL_VALUE, $
                        OutputDir, $
                        Num, $
			SnapBase, $
                        Files)

;;;;;;;;;; We determine the star particle number 

  
NumStars = CALL_EXTERNAL(ObjectFile, $
                       'get_stars_particle_count', /UL_VALUE, $
                        OutputDir, $
                        Num, $
			SnapBase, $
                        Files)

;;;;;;;;;; Now we load the particle data

Pos = fltarr(3, NumPart)
Vel = fltarr(3, NumPart)
Mass = fltarr(NumMass)
U = fltarr(NumGas) 
Rho =fltarr(NumGas)
ID = ulonarr(NumPart)
N_e = fltarr(NumGas)
N_H = fltarr(NumGas)
Hsml = fltarr(NumGas)
SFR = fltarr(NumGas)
stellarage = fltarr(NumStars)
Metallicity = fltarr(NumGas+NumStars)
Type = lonarr(NumPart)
RankID = lonarr(NumPart)
RankIDgas = lonarr(NumGas)
RankIDmass = lonarr(NumMass)
RankIDstars = lon64arr(NumStars)

result  = CALL_EXTERNAL(ObjectFile, $
                       'get_particle_data', /UL_VALUE, $
                        OutputDir, $
                        Num, $
			SnapBase, $
                        Files, $
			NumPart, $
                        NumGas, $
                        NumDm, $
                        NumMass, $
                        NumStars, $
                        Pos, $
                        Vel, $
			ID, $
                        Mass, $
                        U, $
                        Rho, $
                        N_e, $
                        N_H, $
                        Hsml, $
                        SFR, $
                        stellarage, $
                        Metallicity, $
                        Type, $
                        RankID, $
                        RankIDgas, $
                        RankIDmass, $
			RankIDstars)


;;;;;; Now we are all set to read out individual groups (repeatedly if desired)

;;; In this example, we show the 3rd, 4th, 5th and 6th most massive groups.

window, xsize=1000, ysize=1000
!p.multi=[0,2,2]

grstart = 0

for rep=grstart, grstart+3 do begin


  N= rep  ; group number

  finr=  GroupFileNr(N) ; determines in which file this group is stored
  grnr = GroupNr(N)     ; gives the group number within this file
  Len = GroupLen(N)     ; gives the group length
  Lengas = GroupTypeLen(0,N)
  Lendm = GroupTypeLen(1,N) + GroupTypeLen(2,N) + GroupTypeLen(3,N)
  Lenstars = GroupTypeLen(4,N)
  Indices = Lonarr(Len)
  Indicesgas= lonarr(Lengas)
  Indicesstars= lonarr(Lenstars)

  print, "GroupNr=", N, " length=", Len, " gas length=", Lengas, " stars lenght=",  Lenstars, " dm lenght=",Lendm

  result = CALL_EXTERNAL(ObjectFile, $
                       'get_group_indices', /UL_VALUE, $
                        OutputDir, $
                        Num, $
                        SnapBase, $
                        grnr, $
                        finr, $
			Len, $
                        Lengas, $ 
                        Lendm, $ 
                        Lenstars, $ 
			NumPart, $ 
                        NumGas, $
                        NumDm, $ 
                        NumStars, $ 
			ID, $
			RankID, $
                        RankIDgas, $
                        RankIDstars, $ 
                        Indices, $
                        Indicesgas, $ 
                        Indicesstars) 

   print, "Center: ", GroupCenter(*,N)

   xx = Pos(0,Indices) - GroupCenter(0,N)
   yy = Pos(1,Indices) - GroupCenter(1,N)


   ind = where(xx lt -Boxsize/2)
   if ind(0) ne -1 then xx(ind) = xx(ind) + BoxSize
   ind = where(xx gt BoxSize/2)
   if ind(0) ne -1 then xx(ind) = xx(ind) - BoxSize
   ind = where(yy lt -BoxSize/2)
   if ind(0) ne -1 then yy(ind) = yy(ind) + BoxSize
   ind = where(yy gt BoxSize/2)
   if ind(0) ne -1 then yy(ind) = yy(ind) - BoxSize

   len = max([abs(xx),abs(yy)])

   Plot, xx, yy, psym=3 , xrange = [-len,len], yrange=[-len,len], $
        xstyle =1, ystyle=1

   
   ind0 = where(type(indices) eq 0)
   ind1 = where(type(indices) eq 1)
   ind4 = where(type(indices) eq 4)
   ind5 = where(type(indices) eq 5)

   N0 = n_elements(ind0)
   N1 = n_elements(ind1)
   N4 = n_elements(ind4)
   N5 = n_elements(ind5)

   if ind0(0) eq -1 then N0 = 0
   if ind4(0) eq -1 then N4 = 0
   if ind5(0) eq -1 then N5 = 0
   
   print, "N0=", N0, "  N1=", N1, "  N4=", N4, "   N5=", N5   

   if N0 gt 0 then begin
       oplot, xx(ind0), yy(ind0), psym=3, color=255*256L*256L
   endif


   if N4 gt 0 then begin
       oplot, xx(ind4), yy(ind4), color=255*256L, psym=4
   endif


   if N5 gt 0 then begin

       oplot, xx(ind5), yy(ind5), psym=4, color=255, thick=3.0, symsize=1.5
   endif

endfor

end
