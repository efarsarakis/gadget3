
window, xsize=1200, ysize=1200

MinModeCount = 4
TargetBinNummer = 30

Ngrid = 64
box = 1068.0
knyquist = 2*!PI/box * Ngrid/2


count = 0L
Pow_vel_sum = 0
Pow_rnd_sum = 0
Pow_dis1_sum = 0
Pow_dis2_sum = 0

for Num = 0,40 do begin

    path= "//hits/tap/springel/SubsonicTurbulence/sph/box_64_cont/"

   
    
    exts='000'
    exts=exts+strcompress(string(Num),/remove_all)
    exts=strmid(exts,strlen(exts)-3,3)
    
    fname = path+"/powerspec_vel_" +exts+".txt"
 ;  fname = path+"/powerspec_vorticity_" +exts+".txt"


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


    Pow_vel_sum += power_list



    fname = path+"/powerspec_random_" +exts+".txt"


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


    Pow_rnd_sum += power_list






    fname = path+"/powerspec_dis1_" +exts+".txt"


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


    Pow_dis1_sum += power_list




    fname = path+"/powerspec_dis2_" +exts+".txt"


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


    Pow_dis2_sum += power_list







    count++
endfor



    Pow_vel = Pow_vel_sum / count
    Pow_rnd = Pow_rnd_sum / count

    Pow_dis1 = Pow_dis1_sum / count
    Pow_dis2 = Pow_dis2_sum / count
    Pow_dis = Pow_dis1 - Pow_dis2
    
    E_k =  K_list^3 * Pow_vel

    Diss_k = K_list^3 * Pow_dis




    plot, K_list, E_k, /xlog, /ylog , xtitle = "k [ h/kpc ]", ytitle = "E(k)", charsize=2.0, xrange=[0.5*2*!PI/Box,2*knyquist], xstyle=1

    oplot, K_list, E_k, psym=4, color=255
    
;    oplot, K_list, E_k / (Pow_rnd/ampl), linestyle=2

    oplot, knyquist*[1,1],[1.0e-20,1.0e30], linestyle=1

    oplot, K_list, E_K(4) * (k_list/K_list(4))^(-2.0/3), linestyle=1


    oplot, K_list, Diss_k*50, color=255
    
    

end
