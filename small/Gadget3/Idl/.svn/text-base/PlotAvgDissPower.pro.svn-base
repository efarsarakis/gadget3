
window, xsize=1200, ysize=1200


count = 0L
E_inj_sum = 0
E_dis_sum = 0

for Num = 0,199 do begin

    path= "//hits/tap/springel/SubsonicTurbulence/arepo/box_64_cont/"


    Ngrid = 64
    box = 1068.0
    knyquist = 2*!PI/box * Ngrid/2
    
    
    exts='000'
    exts=exts+strcompress(string(Num),/remove_all)
    exts=strmid(exts,strlen(exts)-3,3)
    
    fname = path+"/powerspec_inj_" +exts+".txt"


    openr, 1, fname
    Time = 0.0D 
    Grid = 0L
    Bins = 0L
    readf, 1, Time
    readf, 1, Grid
    readf, 1, Bins
    vel_disp = dblarr(3)
    readf, 1, vel_disp
    da= fltarr(4, bins)
    readf, 1, da
    close,1


    K = da(0,*)
    ModePow = da(1,*)
    ModeCount = da(2,*)
    SumPower =  da(3,*)


    print, num, total(ModeCount * ModePow, /double), total(vel_disp)


; we will do a band averaging of the finely binned points, 
; subject to two conditions:
; We want enough modes per bin in order to reduce the variance in a bin,
; and simultaneously, for large k, we don't want the bins become too narrow.
;
; The first condition is set by "MinModeCount",
; the second by "TargetBinNummer", which is used to compute a minimum 
; logarithmic bin-size.


    MinModeCount = 200
    TargetBinNummer = 60

    MinDlogK = (alog10(max(K)) - alog10(min(K)))/TargetbinNummer


    istart=0
    ind=[istart]

    k_list     = [0]
    power_list = [0]
    count_list = [0]

    repeat begin
        count = total(ModeCount(ind))
        deltak =  (alog10(max(K(ind))) - alog10(min(K(ind))))

        if (deltak ge mindlogk) and (count ge MinModeCount) then begin
            d2 = total(SumPower(ind))/total(ModeCount(ind))
            kk = total(K(ind)*ModeCount(ind))/total(ModeCount(ind))

            
            k_list = [k_list, kk]
            power_list = [power_list, d2]
            count_list = [count_list, total(ModeCount(ind))]
            istart = istart + 1
            ind = [istart]
        endif else begin
            istart = istart + 1
            ind = [ind, istart]
        endelse
    endrep until istart ge Bins

    k_list     = k_list(1:*)
    power_list = power_list(1:*)
    count_list = count_list(1:*)



    E_list = k_list^3 * power_list ; convert to energy-spectrum, E(k) = k^2 * P(k)


    E_inj_sum += E_list




    ;;;;;;;;;;;;;;;;;;;;;;;;




    fname = path+"/powerspec_disp_" +exts+".txt"


    openr, 1, fname
    Time = 0.0D 
    Grid = 0L
    Bins = 0L
    readf, 1, Time
    readf, 1, Grid
    readf, 1, Bins
    vel_disp = dblarr(3)
    readf, 1, vel_disp
    da= fltarr(4, bins)
    readf, 1, da
    close,1


    K = da(0,*)
    ModePow = da(1,*)
    ModeCount = da(2,*)
    SumPower =  da(3,*)


    print, num, total(ModeCount * ModePow, /double), total(vel_disp)


; we will do a band averaging of the finely binned points, 
; subject to two conditions:
; We want enough modes per bin in order to reduce the variance in a bin,
; and simultaneously, for large k, we don't want the bins become too narrow.
;
; The first condition is set by "MinModeCount",
; the second by "TargetBinNummer", which is used to compute a minimum 
; logarithmic bin-size.


    MinModeCount = 200
    TargetBinNummer = 60

    MinDlogK = (alog10(max(K)) - alog10(min(K)))/TargetbinNummer


    istart=0
    ind=[istart]

    k_list     = [0]
    power_list = [0]
    count_list = [0]

    repeat begin
        count = total(ModeCount(ind))
        deltak =  (alog10(max(K(ind))) - alog10(min(K(ind))))

        if (deltak ge mindlogk) and (count ge MinModeCount) then begin
            d2 = total(SumPower(ind))/total(ModeCount(ind))
            kk = total(K(ind)*ModeCount(ind))/total(ModeCount(ind))

            
            k_list = [k_list, kk]
            power_list = [power_list, d2]
            count_list = [count_list, total(ModeCount(ind))]
            istart = istart + 1
            ind = [istart]
        endif else begin
            istart = istart + 1
            ind = [ind, istart]
        endelse
    endrep until istart ge Bins

    k_list     = k_list(1:*)
    power_list = power_list(1:*)
    count_list = count_list(1:*)



    E_list = k_list^3 * power_list ; convert to energy-spectrum, E(k) = k^2 * P(k)


    E_dis_sum += E_list



    count++
endfor



    E_inj = E_inj_sum / count
    E_dis = E_dis_sum / count


    plot, K_list,  E_inj, /xlog, /ylog , xtitle = "k [ h/kpc ]", ytitle = "E(k)", charsize=2.0

    oplot, K_list,  E_inj, psym=4, color=255

    oplot, K_list, E_inj(4) * (k_list/K_list(4))^(-2.0/3), linestyle=1


    oplot, knyquist*[1,1],[1.0e-20,1.0e30], linestyle=1



    oplot, K_list,  E_dis, color=255*256L^2

    oplot, K_list,  E_dis, psym=4, color=255*256L^2



end
