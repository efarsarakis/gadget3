
Num = 63 

Base= "/gpfs/mpa/sbonoli/Milli_BH/Basic/"


FLAG_Group_VelDisp = 0  ; Set this to one if the SO-properties computed by SUBFIND for
                        ; the FOF halos contain velocity dispersions

;--------------------------------------------------

if num ge 1000 then begin
   exts='0000'
   exts=exts+strcompress(string(Num),/remove_all)
   exts=strmid(exts,strlen(exts)-4,4)
endif else begin
   exts='000'
   exts=exts+strcompress(string(Num),/remove_all)
   exts=strmid(exts,strlen(exts)-3,3)
endelse


skip = 0L
skip_sub = 0L
fnr = 0L

repeat begin

    f = Base + "/groups_" + exts +"/subhalo_tab_"+exts +"."+strcompress(string(fnr),/remove_all)

    Ngroups = 0L
    TotNgroups = 0L
    Nids = 0L
    TotNids = 0LL
    NTask = 0L
    Nsubgroups = 0L
    TotNsubgroups = 0L

    openr,1,f
    readu,1, Ngroups, TotNgroups, Nids, TotNids, NTask, Nsubgroups, TotNsubgroups

    if fnr eq 0 then begin
        GroupLen = lonarr(TotNgroups)
        GroupOffset = lonarr(TotNgroups)
        GroupMass = fltarr(TotNgroups)
        GroupPos = fltarr(3, TotNgroups)
        Group_M_Mean200 = fltarr(TotNgroups)
        Group_R_Mean200 = fltarr(TotNgroups)
        Group_M_Crit200 = fltarr(TotNgroups)
        Group_R_Crit200 = fltarr(TotNgroups)
        Group_M_TopHat200 = fltarr(TotNgroups)
        Group_R_TopHat200 = fltarr(TotNgroups)
        Group_VelDisp_Mean200 = fltarr(TotNgroups)
        Group_VelDisp_Crit200 = fltarr(TotNgroups)
        Group_VelDisp_TopHat200 = fltarr(TotNgroups)
        GroupContaminationCount = lonarr(TotNgroups)
        GroupContaminationMass = fltarr(TotNgroups)
        GroupNsubs = lonarr(TotNgroups)
        GroupFirstSub = lonarr(TotNgroups)

        SubhaloLen = lonarr(TotNsubgroups)
        SubhaloOffset = lonarr(TotNsubgroups)
        SubhaloParent = lonarr(TotNsubgroups)
        SubhaloMass = fltarr(TotNsubgroups)
        SubhaloPos = fltarr(3, TotNsubgroups)
        SubhaloVel = fltarr(3, TotNsubgroups)
        SubhaloCM = fltarr(3, TotNsubgroups)
        SubhaloSpin = fltarr(3, TotNsubgroups)
        SubhaloVelDisp = fltarr(TotNsubgroups)
        SubhaloVmax = fltarr(TotNsubgroups)
        SubhaloVmaxRad = fltarr(TotNsubgroups)
        SubhaloHalfmassRad = fltarr(TotNsubgroups)
        SubhaloIDMostbound = lonarr(TotNsubgroups)
        SubhaloGrNr = lonarr(TotNsubgroups)
        SubhaloMassTab = fltarr(6, TotNsubgroups)
    endif

    if Ngroups gt 0 then begin

        locLen = lonarr(Ngroups)
        locOffset = lonarr(Ngroups)
        locMass = fltarr(Ngroups)
        locPos = fltarr(3, Ngroups)
        loc_M_Mean200 = fltarr(Ngroups)
        loc_R_Mean200 = fltarr(Ngroups)
        loc_M_Crit200 = fltarr(Ngroups)
        loc_R_Crit200 = fltarr(Ngroups)
        loc_M_TopHat200 = fltarr(Ngroups)
        loc_R_TopHat200 = fltarr(Ngroups)
        loc_VelDisp_Mean200 = fltarr(Ngroups)
        loc_VelDisp_Crit200 = fltarr(Ngroups)
        loc_VelDisp_TopHat200 = fltarr(Ngroups)
        locContaminationCount = lonarr(Ngroups)
        locContaminationMass = fltarr(Ngroups)
        locNsubs = lonarr(Ngroups)
        locFirstSub = lonarr(Ngroups)

        readu,1, loclen
        readu,1, locOffset
        readu,1, locMass
        readu,1, locPos
        readu,1, loc_M_Mean200
        readu,1, loc_R_Mean200
        readu,1, loc_M_Crit200
        readu,1, loc_R_Crit200
        readu,1, loc_M_TopHat200
        readu,1, loc_R_TopHat200
        if FLAG_Group_VelDisp ne 0 then begin
          readu,1, loc_VelDisp_Mean200
          readu,1, loc_VelDisp_Crit200
          readu,1, loc_VelDisp_TopHat200
        endif
        readu,1, locContaminationCount
        readu,1, locContaminationMass
        readu,1, locNsubs
        readu,1, locFirstSub

        GroupLen(skip:skip+Ngroups-1) = locLen(*)
        GroupOffset(skip:skip+Ngroups-1) = locOffset(*)
        GroupMass(skip:skip+Ngroups-1) = locMass(*)
        GroupPos(*, skip:skip+Ngroups-1) = locPos(*,*)
        Group_M_Mean200(skip:skip+Ngroups-1) = loc_M_Mean200(*)
        Group_R_Mean200(skip:skip+Ngroups-1) = loc_R_Mean200(*)
        Group_M_Crit200(skip:skip+Ngroups-1) = loc_M_Crit200(*)
        Group_R_Crit200(skip:skip+Ngroups-1) = loc_R_Crit200(*)
        Group_M_TopHat200(skip:skip+Ngroups-1) = loc_M_TopHat200(*)
        Group_R_TopHat200(skip:skip+Ngroups-1) = loc_R_TopHat200(*)
        Group_VelDisp_Mean200(skip:skip+Ngroups-1) = loc_VelDisp_Mean200(*)
        Group_VelDisp_Crit200(skip:skip+Ngroups-1) = loc_VelDisp_Crit200(*)
        Group_VelDisp_TopHat200(skip:skip+Ngroups-1) = loc_VelDisp_TopHat200(*)
        GroupContaminationCount(skip:skip+Ngroups-1) = locContaminationCount(*)
        GroupContaminationMass(skip:skip+Ngroups-1) = locContaminationMass(*)
        GroupNsubs(skip:skip+Ngroups-1) = locNsubs(*)
        GroupFirstsub(skip:skip+Ngroups-1) = locFirstSub(*)

        skip+= Ngroups
    endif


    if Nsubgroups gt 0 then begin

        locLen = lonarr(Nsubgroups)
        locOffset = lonarr(Nsubgroups)
        locParent = lonarr(Nsubgroups)
        locMass = fltarr(Nsubgroups)
        locPos = fltarr(3, Nsubgroups)
        locVel = fltarr(3, Nsubgroups)
        locCM = fltarr(3, Nsubgroups)
        locSpin= fltarr(3, Nsubgroups)
        locVelDisp = fltarr(Nsubgroups)
        locVmax = fltarr(Nsubgroups)
        locVmaxRad = fltarr(Nsubgroups)
        locHalfMassRad = fltarr(Nsubgroups)
        locIDMostBound = lonarr(Nsubgroups)
        locGrNr = lonarr(Nsubgroups)
        locSubhaloMassTab = fltarr(6, Nsubgroups)

        readu,1, loclen
        readu,1, locOffset
        readu,1, locParent
        readu,1, locMass
        readu,1, locPos
        readu,1, locVel
        readu,1, locCM
        readu,1, locSpin
        readu,1, locVelDisp
        readu,1, locVmax
        readu,1, locVmaxRad
        readu,1, locHalfMassRad
        readu,1, locIDMostBound
        readu,1, locGrNr
        readu,1, locSubhaloMassTab

        SubhaloLen(skip_sub:skip_sub+Nsubgroups-1) = locLen(*)
        SubhaloOffset(skip_sub:skip_sub+Nsubgroups-1) = locOffset(*)
        SubhaloParent(skip_sub:skip_sub+Nsubgroups-1) = locParent(*)
        SubhaloMass(skip_sub:skip_sub+Nsubgroups-1) = locMass(*)
        SubhaloPos(*, skip_sub:skip_sub+Nsubgroups-1) = locPos(*,*)
        SubhaloVel(*, skip_sub:skip_sub+Nsubgroups-1) = locVel(*,*)
        SubhaloCM(*, skip_sub:skip_sub+Nsubgroups-1) = locCM(*,*)
        SubhaloSpin(*, skip_sub:skip_sub+Nsubgroups-1) = locSpin(*,*)
        SubhaloVeldisp(skip_sub:skip_sub+Nsubgroups-1) = locVeldisp(*)
        SubhaloVmax(skip_sub:skip_sub+Nsubgroups-1) = locVmax(*)
        SubhaloVmaxRad(skip_sub:skip_sub+Nsubgroups-1) = locVmaxRad(*)
        SubhaloHalfmassRad(skip_sub:skip_sub+Nsubgroups-1) = locHalfmassRad(*)
        SubhaloIDMostBound(skip_sub:skip_sub+Nsubgroups-1) = locIDMostBound(*)
        SubhaloGrNr(skip_sub:skip_sub+Nsubgroups-1) = locGrNr(*)
        SubhaloMassTab(*, skip_sub:skip_sub+Nsubgroups-1) =  locSubhaloMassTab(*,*)

        skip_sub+= Nsubgroups
    endif

    close, 1

    fnr++

endrep until fnr eq NTask

print
print, "TotNgroups   =", TotNgroups
print, "TotNsubgroups=", TotNsubgroups
print
print, "Largest group of length ", GroupLen(0)," has", GroupNsubs(0)," substructures"
print





;;;; At this point, we have read-in the subhalo catalogue
;;;; We can now read for an abitrary subhalo the detailed particle data from the
;;;; "posvel" files.



SubhaloNr = 5 ; this is the number we want

Off = SubhaloOffset(SubhaloNr)
Len = SubhaloLen(SubhaloNr)



Pos = fltarr(3, Len)
Vel = fltarr(3, Len)
Type = lonarr(Len)
Mass = fltarr(Len)
FormationTime = fltarr(Len) ;; for stars




left = Len
skip = Off
done = 0L

fnr = 0L

repeat begin

    f = Base + "/groups_" + exts +"/subhalo_posvel_"+exts +"."+strcompress(string(fnr),/remove_all)

    Ngroups = 0L
    TotNgroups = 0L
    Nids = 0L
    TotNids = 0LL
    NTask = 0L
    fileoff= 0L	
    scalefactor = 0.0D


    openr,1,f
    readu,1, Ngroups, TotNgroups, Nids, TotNids, NTask, fileoff, scalefactor

    if skip lt Nids then begin

	if (Nids - skip) ge left then begin
            loclen = left
        endif else begin
	    loclen = Nids - skip
	endelse

        print, "loclen=", loclen
	
	locPos = fltarr(3, locLen)
	locVel = fltarr(3, locLen)
	locType = bytarr(locLen)
	locMass = fltarr(locLen)
	locFormationTime = fltarr(locLen)

        point_lun, 1, 36 + skip*12L
	readu, 1, locPos
        point_lun, 1, 36 + Nids*12L+ skip*12L
	readu, 1, locVel
        point_lun, 1, 36 + 2L*Nids*12L+ skip
	readu, 1, locType
        point_lun, 1, 36 + 2L*Nids*12L+ Nids + skip*4L
	readu, 1, locMass
        point_lun, 1, 36 + 2L*Nids*12L+ Nids*5L + skip*4L
	readu, 1, locFormationTime

	Pos(*, done:done+locLen-1) = locPos(*,*)	
	Vel(*, done:done+locLen-1) = locVel(*,*)	
	Type(done:done+locLen-1) = locType(*)	
	Mass(done:done+locLen-1) = locMass(*)	
	FormationTime(done:done+locLen-1) = locFormationTime(*)	

        done += locLen
	left -= locLen
	skip = 0

    endif else begin	

	skip -= Nids

    endelse

    close,1
   
    fnr++
        
endrep until (fnr eq NTask) or (left eq 0)




; get group center

cx = subhalopos(0, subhalonr)
cy = subhalopos(1, subhalonr)
cz = subhalopos(2, subhalonr)


plot, pos(0, *) - cx, pos(1,*) - cy, psym=3

ind = where(type eq 4)

if ind(0) ne -1 then begin

  oplot, pos(0,ind) -cx, pos(1,ind) - cy, psym=4, color=255*256L+255
	
endif

print, total(mass), subhalomass(subhalonr)



;; check whether the group has a black hole

ind = where(Type eq 5)

if ind(0) ne -1 then begin
   
   oplot, pos(0,ind) - cx, pos(1, ind)- cy, psym=5, color=255, thick=4.0

endif




end












