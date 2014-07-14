
Num = 004

exts='000'
exts=exts+strcompress(string(Num),/remove_all)
exts=strmid(exts,strlen(exts)-3,3)

fname = "./powerspec_" +exts+".txt"

openr, 1, fname
Time = 0.0D 
Bins = 0L
Mass = 0.0D
Npart = 0LL
readf, 1, Time
readf, 1, Bins
readf, 1, Mass
readf, 1, Npart
da1= fltarr(10, bins)
readf, 1, da1
readf, 1, Time
readf, 1, Bins
readf, 1, Mass
readf, 1, Npart
da2= fltarr(10, bins)
readf, 1, da2
close,1


K_A = da1(0,*)
Delta2_A = da1(1,*)
Shot_A = da1(2,*)
ModePow_A = da1(3,*)
ModeCount_A = da1(4,*)
Delta2Uncorrected_A = da1(5,*)
ModePowUncorrected_A = da1(6,*)
Specshape_A =  da1(7,*)
SumPower_A =  da1(8,*)
ConvFac_A =  da1(9,*)

K_B = da2(0,*)
Delta2_B = da2(1,*)
Shot_B = da2(2,*)
ModePow_B = da2(3,*)
ModeCount_B = da2(4,*)
Delta2Uncorrected_B = da2(5,*)
ModePowUncorrected_B = da2(6,*)
Specshape_B =  da2(7,*)
SumPower_B =  da2(8,*)
ConvFac_B =  da2(9,*)


; we will do a band averaging of the finely binned points, 
; subject to two conditions:
; We want enough modes per bin in order to reduce the variance in a bin,
; and simultaneously, for large k, we don't want the bins become too narrow.
;
; The first condition is set by "MinModeCount",
; the second by "TargetBinNummer", which is used to compute a minimum 
; logarithmic bin-size.


MinModeCount = 20
TargetBinNummer = 60

MinDlogK = (alog10(max(K_A)) - alog10(min(K_A)))/TargetbinNummer


istart=0
ind=[istart]
k_list_A = [0]
Delta2_list_A = [0]
count_list_A = [0]
repeat begin
    count = total(modecount_a(ind))
    deltak =  (alog10(max(K_A(ind))) - alog10(min(K_A(ind))))

    if (deltak ge mindlogk) and (count ge MinModeCount) then begin
       d2 = total(SumPower_A(ind))/total(ModeCount_A(ind))
       b = fix(total(double(ind)*ModeCount_A(ind))/total(ModeCount_A(ind)))
       kk = K_A(b)
       d2 = ConvFac_A(b)*d2*Specshape_A(b)
       k_list_A = [k_list_A, kk]
       Delta2_list_A = [Delta2_list_A, d2]
       count_list_A = [count_list_A, total(ModeCount_A(ind))]
       istart = istart + 1
       ind = [istart]
    endif else begin
       istart = istart + 1
       ind = [ind, istart]
    endelse
endrep until istart ge Bins
K_list_A = k_list_A(1:*)
Delta2_list_A = delta2_list_A(1:*)
Count_list_A = count_list_A(1:*)



istart=0
ind=[istart]
k_list_B = [0]
Delta2_list_B = [0]
count_list_B = [0]
repeat begin
   count = total(modecount_B(ind))
   deltak =  (alog10(max(K_B(ind))) - alog10(min(K_B(ind))))

   if (deltak ge mindlogk) and (count ge MinModeCount) then begin
      d2 = total(SumPower_B(ind))/total(ModeCount_B(ind))
      b = fix(total(double(ind)*ModeCount_B(ind))/total(ModeCount_B(ind)))
      kk = K_B(b)
      d2 = ConvFac_B(b) * d2 * Specshape_B(b)
      k_list_B = [k_list_B, kk]
      Delta2_list_B = [Delta2_list_B, d2]
      count_list_B = [count_list_B, total(ModeCount_B(ind))]
      istart = istart + 1
      ind = [istart]
   endif else begin
      istart = istart + 1
      ind = [ind, istart]
   endelse
endrep until istart ge Bins
K_list_B = k_list_B(1:*)
Delta2_list_B = delta2_list_B(1:*)
Count_list_B = count_list_B(1:*)



; we discard the last bin of the top-level mesh, and then switch to the folded one




plot, [K_list_A, K_list_B], [Delta2_list_A, Delta2_list_B], /xlog, /ylog , $
   xtitle = "k [ h/kpc ]", ytitle = "Delta^2(k)", psym=4

oplot, K_list_B, Delta2_list_B, color=255, psym=4  
oplot, K_list_A, Delta2_list_A, color=255*256L, psym=4 

oplot, K_list_B, Delta2_list_B, color=255
oplot, K_list_A, Delta2_list_A, color=255*256L

oplot, K_B, Shot_B, color=255+256L*255

ind=where(Delta2_B gt 0)
;oplot, K_B(ind), Delta2_B(ind)

ind=where(Delta2_A gt 0)
;oplot, K_A(ind), Delta2_A(ind)


end
