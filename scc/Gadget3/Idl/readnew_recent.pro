; example:
; 
; Reading the header information:
; readnew,'snap_000',head,'HEAD'
; print,head.npart
;
; Reading the Poisitions:
; readnew,'snap_000',x,'POS '
; print,'X-Range:',min(x(0,*)),max(x(0,*))
; print,'Y-Range:',min(x(1,*)),max(x(1,*))
;
; Reading the Masses
; readnew,'snap_000',m,'MASS'
; ATTENTION! GADGET allows to save space by storing the masses
;            in the massarr if all particle of a certain species
;            have the same mass !!!!
;            Therefore its recomendet to have a look at the header first !
; readnew,'snapo_000',head,'HEAD',type='HEAD'
; if head.massarr(0) EQ 0 then print,m(0) else print,head.massarr(0)


FUNCTION IS_DEF, x
aux = SIZE(x)
RETURN, aux(N_ELEMENTS(aux)-2) NE 0
END

PRO readnew,name,x,label,type=type,debug=debug,quiet=quiet,parttype=parttype,$
            status=status,ndim=ndim,is_present=is_present,double=double,ignorema=ignorema
@setconst
   status=-1
   IF NOT(IS_DEF(name)) THEN BEGIN
      print,'Usage: readnew,name,x,label,type=type,debug=debug'
      print,'       name  : filename'
      print,'       x     : variable containing the result'
      print,'       label : names the block identifier ("HEAD","POS ",...)'
      print,'       type  : names the data type ("FLOAT","FLOAT3","LONG","HEAD")'
      print,'       debug : give some debug information'
      return
   END

   blabel=BYTE(label+"    ")
   label=STRING(blabel(0:3))
   IF NOT(IS_DEF(is_present)) THEN is_present=[1,0,0,0,0,0]
   IF NOT(IS_DEF(TYPE)) THEN BEGIN
      IF strcmp(label,"HEAD")  THEN  type="HEAD"
      IF strcmp(label,"ID  ") OR strcmp(label,"MASS") OR strcmp(label,"IDU ") THEN BEGIN
         IF strcmp(label,"ID  ") OR strcmp(label,"IDU ") THEN type="LONG" ELSE type="FLOAT"
         is_present=[1,1,1,1,1,1]
      END
      IF strcmp(label,"POT ") or strcmp(label,"TSTP") THEN BEGIN
         type="FLOAT"
         is_present=[1,1,1,1,1,1]
      END
      IF strcmp(label,"TSTP")  THEN BEGIN
         type="FLOAT"
         is_present=[1,1,1,1,1,1]
      END
      IF strcmp(label,"POS ") OR strcmp(label,"VEL ") OR strcmp(label,"ACCE") THEN BEGIN
         type="FLOATN"
         ndim=3
         is_present=[1,1,1,1,1,1]
      END
      IF strcmp(label,"BFLD") OR strcmp(label,"BFSM") OR $ 
         strcmp(label,"ROTB") OR strcmp(label,"SRTB") THEN BEGIN
         type="FLOATN"
         ndim=3
         is_present=[1,0,0,0,0,0]
      END
      IF strcmp(label,"Zs  ") THEN BEGIN
         type="FLOATN"
         IF NOT(IS_DEF(ndim)) THEN ndim=8
         is_present=[1,0,0,0,1,0]
      END   
      IF strcmp(label,"Z   ")  THEN BEGIN
         type="FLOATN"
         IF NOT(IS_DEF(ndim)) THEN ndim=1
         is_present=[1,0,0,0,1,0]
      END   
      IF strcmp(label,"AGE ") OR strcmp(label,"iM  ") OR $
         strcmp(label,"SLg ") OR strcmp(label,"HSMS") OR $
         strcmp(label,"SUB ") THEN BEGIN
         type="FLOAT"
         is_present=[0,0,0,0,1,0]
      END
      IF strcmp(label,"HSBP") THEN BEGIN
         type="FLOAT"
         is_present=[0,0,1,1,0,0]    
      END

      IF strcmp(label,"TIPS") THEN BEGIN
         type="FLOATN"
         ndim=9
         is_present=[1,1,1,1,1,1]
      END

      IF strcmp(label,"DIPS") THEN BEGIN
         type="FLOATN"
         ndim=36
         is_present=[1,1,1,1,1,1]
      END

      IF strcmp(label,"CACO") THEN BEGIN
         type="FLOATN"
         ndim=1
         is_present=[1,1,1,1,1,1]
      END

      IF strcmp(label,"FLDE") THEN BEGIN
         type="FLOATN"
         ndim=1
         is_present=[1,1,1,1,1,1]
      END

      IF strcmp(label,"STDE") THEN BEGIN
         type="FLOATN"
         ndim=1
         is_present=[1,1,1,1,1,1]
      END

      IF strcmp(label,"PSDE") THEN BEGIN
         type="FLOATN"
         ndim=1
         is_present=[1,1,1,1,1,1]
      END

      IF strcmp(label,"ANRA") THEN BEGIN
         type="FLOATN"
         ndim=3
         is_present=[1,1,1,1,1,1]
      END

      IF strcmp(label,"LACA") THEN BEGIN
         type="FLOATN"
         ndim=20
         is_present=[1,1,1,1,1,1]
      END

      IF strcmp(label,"SHIN") THEN BEGIN
         type="FLOATN"
         ndim=3
         is_present=[1,1,1,1,1,1]
      END

      IF strcmp(label,"SHOR") THEN BEGIN
         type="FLOATN"
         ndim=9
         is_present=[1,1,1,1,1,1]
      END

      IF strcmp(label,"INDE") THEN BEGIN
         type="FLOATN"
         ndim=1
         is_present=[1,1,1,1,1,1]
      END

      IF NOT(IS_DEF(TYPE)) THEN type="FLOAT"
   END
   IF NOT(IS_DEF(ndim)) THEN ndim=1

   b4=bytarr(4)
   bl=0L

   get_lun,myfile
   myname=name
   numfiles=1

   ss=SIZE(FINDFILE(myname))
   IF ss(0) EQ 0 THEN BEGIN
      myname=name+'.0'
      ss=SIZE(FINDFILE(myname))
      IF ss(0) EQ 0 THEN BEGIN
         print,'Cant find file ',name,' or ',myname,' ...'
         free_lun,myfile
         return
      END ELSE BEGIN
         npart=lonarr(6)
         readnew,myname,h,'HEAD',debug=debug,quiet=quiet,ignorema=ignorema,double=double
         npart=npart+h.npart
         WHILE ss(0) NE 0 DO BEGIN
            myname=name+'.'+strcompress(string(numfiles),/remove_all)
            ss=SIZE(FINDFILE(myname))
            IF ss(0) NE 0 THEN BEGIN
               numfiles=numfiles+1
               readnew,myname,h,'HEAD',debug=debug,quiet=quiet,ignorema=ignorema,double=double
               npart=npart+h.npart
            END
         END
         IF IS_DEF(debug) THEN print,'Found',numfiles,' files ...'
         IF IS_DEF(debug) THEN print,'TotNP',npart
      END
   END

   IF numfiles GT 1 THEN BEGIN
      IF strcmp(type,"HEAD") THEN BEGIN
         h.npart(*)=npart(*)
         x=h
      END ELSE BEGIN
         IF IS_DEF(parttype) THEN BEGIN
            nalloc=is_present(parttype)*npart(parttype)
         END ELSE BEGIN
; Keep explicite sum, as TOTAL(npart) would be a real !!!!
            nalloc=0L
            FOR j=0,5 DO nalloc=nalloc+is_present(j)*npart(j)
         END

         IF IS_DEF(double) THEN BEGIN
            IF strcmp(type,"FLOAT") AND (nalloc GT 0) THEN x=dblarr(nalloc) 
            IF strcmp(type,"FLOATN") AND (nalloc GT 0) THEN x=dblarr(ndim,nalloc)
            IF strcmp(type,"LONG") AND (nalloc GT 0) THEN x=ulon64arr(nalloc)
         END ELSE BEGIN
            IF strcmp(type,"FLOAT") AND (nalloc GT 0) THEN x=fltarr(nalloc) 
            IF strcmp(type,"FLOATN") AND (nalloc GT 0) THEN x=fltarr(ndim,nalloc)
            IF strcmp(type,"LONG") AND (nalloc GT 0) THEN x=ulonarr(nalloc)
         END

         nstart=0L
         FOR i=0L,numfiles-1 DO BEGIN
            myname=name+'.'+strcompress(string(i),/remove_all)
            readnew,myname,htmp,'HEAD',debug=debug,quiet=quiet,status=status,ignorema=ignorema,double=double


;    IF IS_DEF(ignorema) THEN BEGIN
;print,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
;print,"!!!!!!!!!!!!!!!!! Patching !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
;print,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
;            htmp.massarr(*)=0
;print,"Setting massarr members to zero"
;print,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
;print,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
;    END

            gof=1
            IF IS_DEF(parttype) THEN BEGIN
               IF htmp.npart(parttype) EQ 0 THEN gof=0
            END
            IF gof EQ 1 THEN BEGIN
               rfile=1
               IF strcmp(label,"MASS") THEN BEGIN
                  IF IS_DEF(parttype) THEN BEGIN
                     IF (htmp.npart(parttype) EQ 0) THEN BEGIN
                        rfile=0
                     END ELSE BEGIN
                        IF (htmp.massarr(parttype) GT 0) THEN BEGIN
                           rfile=0
                           status=1
                           IF IS_DEF(double) THEN BEGIN
                              xpart=dblarr(htmp.npart(parttype))
                           END ELSE BEGIN
                              xpart=fltarr(htmp.npart(parttype))
                           END
                           xpart(*)=htmp.massarr(parttype)
                        END
                     END
                  END ELSE BEGIN
                     nmassread=0L
                     nmasstable=0L
                     FOR im=0,5 DO BEGIN
                        IF htmp.npart(im) GT 0 THEN BEGIN
                           nmasstable=nmasstable+htmp.npart(im)
                           IF htmp.massarr(im) EQ 0 THEN nmassread=nmassread+1
                        END
                     END
                     IF nmassread EQ 0 THEN BEGIN
                        rfile=0
                        If IS_DEF(double) THEN BEGIN
                            IF nmasstable GT 0 THEN xpart=dblarr(nmasstable) ELSE xpart=0.0d0
                        END ELSE BEGIN
                            IF nmasstable GT 0 THEN xpart=fltarr(nmasstable) ELSE xpart=0.0
                        END
                        istart=0L
                        FOR im=0,5 DO BEGIN
                           IF htmp.npart(im) GT 0 THEN BEGIN
                              xpart(istart:istart+htmp.npart(im)-1) = htmp.massarr(im) 
                              istart=istart+htmp.npart(im)
                           END
                        END
                     END
                  END
               END
               IF rfile EQ 1 THEN $
                   readnew,myname,xpart,label,type=type,debug=debug,$
                           quiet=quiet,parttype=parttype,status=status,$
                           ndim=ndim,is_present=is_present,ignorema=ignorema,$
                           double=double
               IF (status GT -1) THEN BEGIN
                  ss=SIZE(xpart)
                  IF strcmp(type,"FLOATN") THEN BEGIN
                     IF ss(0) EQ 1 THEN nmax=1L ELSE nmax=ss(2)
                     FOR ip=0L,nmax-1 DO BEGIN
                        FOR k=0,ndim-1 DO BEGIN
                           x(k,nstart)=xpart(k,ip)
                        END
                        nstart=nstart+1
                     END
                     IF IS_DEF(debug) THEN $
                         print,'Reading data',nstart,nstart+ss(2)-1,' ...'
                  END ELSE BEGIN
                     x(nstart:nstart+ss(1)-1)=xpart(*)
                     nstart=nstart+ss(1)
                     IF IS_DEF(debug) THEN $
                         print,'Reading data',nstart,nstart+ss(1)-1,' ...'
                  END
               END
            END ELSE BEGIN
               IF IS_DEF(debug) THEN $
                   print,'Skipping file, does not contain particle of type ',parttype
            END
         END
      END
   END ELSE BEGIN

      IF IS_DEF(debug) THEN print,'Testing File ',myname
      openr,myfile,myname
      readu,myfile,bl
      close,myfile

      IF bl EQ 8 THEN BEGIN
         openr,myfile,name,/f77 
         IF IS_DEF(debug) THEN print,'Open file in normal mode ..,'
      END ELSE BEGIN
         openr,myfile,name,/f77,/SWAP_ENDIAN
         IF IS_DEF(debug) THEN print,'Open file with SWAP_ENDIAN ..,'
      END

      found=0

      WHILE(NOT(EOF(myfile)) AND found EQ 0) DO BEGIN
         readu,myfile,b4,bl
         thislabel=STRING(b4)
         IF IS_DEF(debug) THEN print,thislabel," <=> ",label,bl
         IF strcmp(thislabel,'HEAD') NE 0 THEN BEGIN
            IF bl NE 264 THEN BEGIN
               bl=264L
               print,'Patching length for HEAD !!!'
            END 
         END

         IF strcmp(thislabel,label) EQ 0 THEN BEGIN
            point_lun,-myfile,pos
            point_lun,myfile,ulong64(pos)+ulong64(bl)
         END ELSE BEGIN
            found=1
         END
      ENDWHILE

      IF found EQ 0 THEN BEGIN
         status=-1
         IF NOT(IS_DEF(quiet)) THEN $
           print,'File ',name,' does not contain a block labled with "',label,'" !'

         xH=0.76

         IF strcmp(label,"TEMP") THEN BEGIN
            IF NOT(IS_DEF(quiet)) THEN $ 
              print,"WARNING: calculting temperature, assuming no starformation !!!"
              readnew,name,u,'U   ',quiet=quiet,parttype=parttype,status=status,ignorema=ignorema,double=double      
              IF status NE -1 THEN BEGIN
                 g1 = 5.0 / 3.0 - 1.

                 readnew,name,xt,'DI  ',quiet=quiet,parttype=parttype,ignorema=ignorema,double=double
                 IF NOT(IS_DEF(xt)) THEN BEGIN
                    readnew,name,xNe,'NE  ',quiet=quiet,parttype=parttype,ignorema=ignorema,double=double
                    IF NOT(IS_DEF(xNe)) THEN BEGIN
                       xNe=1.0
                       nHe_fak=3.
                    END ELSE BEGIN
                       nHe_fak=1.
                    END
                 END ELSE BEGIN
                      xNe=1.0
                      nHe_fak=3.
                 END
                 yhelium = (1. - xH) / (4. * xH)
                 mu = (1. + 4. * yhelium) / (1. + nHe_fak * yhelium + xNe)

                 x = FLOAT(g1 / kcgs * mp * u * e_unit / m_unit * mu)
                status=1
            END
         END
         IF strcmp(label,"MHOT") THEN BEGIN
            readnew,name,h,'HEAD',ignorema=ignorema,double=double
            IF h.npart(0) GT 0 THEN BEGIN
               IF h.massarr(0) EQ 0.0 THEN BEGIN
                  readnew,name,m,'MASS',parttype=parttype,ignorema=ignorema,double=double
                  gmas=m(0:h.npart(0)-1)
               END ELSE BEGIN
                  If IS_DEF(double) THEN BEGIN
                     gmas=dblarr(h.npart(0))
                  END ELSE BEGIN
                     gmas=fltarr(h.npart(0))
                  END
                  gmas(*)=h.massarr(0) 
               END
               readnew,name,mhi,'MHI ',quiet=quiet,parttype=parttype,ignorema=ignorema,double=double
               IF IS_DEF(mhi) THEN BEGIN
                  IF NOT(IS_DEF(quiet)) THEN print,"Calculting MHOT from MHI ..."
                  xfr=mhi/gmas/xH
               END ELSE BEGIN
                  xfr=0.
               END
               x = (1 - xfr) * gmas
               status=1
            END
         END
         IF strcmp(label,"RHOT") THEN BEGIN
            readnew,name,rho,'RHO ',parttype=parttype,ignorema=ignorema,double=double
            readnew,name,mhi,'MHI ',quiet=quiet,parttype=parttype,ignorema=ignorema,double=double
            IF NOT(IS_DEF(mhi)) THEN BEGIN
               IF NOT(IS_DEF(quiet)) THEN print,"Returning normal density ..."
               x=rho
            END ELSE BEGIN
               IF NOT(IS_DEF(quiet)) THEN print,"Calculting RHOT from MHI ..."
               readnew,name,h,'HEAD',ignorema=ignorema,double=double
               IF h.massarr(0) EQ 0.0 THEN BEGIN
                  readnew,name,m,'MASS',parttype=parttype,ignorema=ignorema,double=double
                  gmas=m(0:h.npart(0)-1)
               END ELSE BEGIN
                  gmas=h.massarr(0) 
               END
               xfr=mhi/gmas/xH
               x = (1 - xfr) * rho
            END
            status=1
         END
         IF strcmp(label,"MCLD") THEN BEGIN
            readnew,name,mhi,'MHI ',quiet=quiet,parttype=parttype,ignorema=ignorema,double=double
            IF IS_DEF(mhi) THEN BEGIN
               IF NOT(IS_DEF(quiet)) THEN print,"Calculting MCLD from MHI ..."
               readnew,name,h,'HEAD',ignorema=ignorema,double=double
               IF h.massarr(0) EQ 0.0 THEN BEGIN
                  readnew,name,m,'MASS',parttype=parttype,ignorema=ignorema,double=double
                  gmas=m(0:h.npart(0)-1)
               END ELSE BEGIN
                  gmas=h.massarr(0) 
               END
               xfr=mhi/gmas/xH
               x = xfr * gmas
               status=1
            END
         END
         IF strcmp(label,"RCLD") THEN BEGIN
            readnew,name,rho,'RHO ',parttype=parttype,ignorema=ignorema,double=double
            readnew,name,mhi,'MHI ',quiet=quiet,parttype=parttype,ignorema=ignorema,double=double
            IF NOT(IS_DEF(mhi)) THEN BEGIN
               IF NOT(IS_DEF(quiet)) THEN print,"Returning normal density ..."
               x=rho
            END ELSE BEGIN
               IF NOT(IS_DEF(quiet)) THEN print,"Calculting RHOT from MHI ..."
               readnew,name,h,'HEAD',ignorema=ignorema,double=double
               IF h.massarr(0) EQ 0.0 THEN BEGIN
                  readnew,name,m,'MASS',parttype=parttype,ignorema=ignorema,double=double
                  gmas=m(0:h.npart(0)-1)
               END ELSE BEGIN
                  gmas=h.massarr(0) 
               END
               xfr=mhi/gmas/xH
               x = xfr * rho
            END
            status=1
         END
         IF strcmp(label,"MASS") THEN BEGIN
            readnew,name,h,'HEAD',ignorema=ignorema,double=double
            ntotm=TOTAL(h.npart,/int)
            If IS_DEF(double) THEN BEGIN
               x=dblarr(ntotm)
            END ELSE BEGIN
               x=fltarr(ntotm)
            END
            st=0L
            str=0L
            FOR i=0,5 DO BEGIN
               IF h.npart(i) GT 0 THEN BEGIN
;                  print,'filling ',st,st+h.npart(i)-1
                  IF h.massarr(i) GT 0 THEN BEGIN
                     x(st:st+h.npart(i)-1)=h.massarr(i)
                     st=st+h.npart(i)
                  END ELSE BEGIN
                     print,"Something is wrong ! No MASS block and no massarr entry !"
                  END
               END
            END
            IF IS_DEF(parttype) THEN BEGIN
               IF h.npart(parttype) GT 0 THEN valid=1 ELSE valid=0
            END ELSE valid=1
            IF valid EQ 1 THEN BEGIN
               IF IS_DEF(parttype) THEN BEGIN
                  IF ntotm NE h.npart(parttype) THEN BEGIN
                     pstart=0L
                     FOR i=0,parttype-1 DO pstart=pstart+h.npart(i)*is_present(i)
                     IF NOT(IS_DEF(quiet)) THEN $
                        print,'Restricting return value to',pstart,pstart+h.npart(parttype)-1
                     IF strcmp(type,"FLOATN") THEN $
                        x=x(*,pstart:pstart+h.npart(parttype)-1) $
                     ELSE $
                        x=x(pstart:pstart+h.npart(parttype)-1) 
                  END
               END
            END
            status=1
         END
      END ELSE BEGIN
         IF strcmp(type,"HEAD")   THEN BEGIN
            npart=lonarr(6)	
            massarr=dblarr(6)
            time=0.0D
            redshift=0.0D
            flag_sfr=0L
            flag_feedback=0L
            partTotal=lonarr(6)
            flag_cooling=0L
            num_files=0L
            BoxSize=0.0D
            Omega0=0.0D
            OmegaLambda=0.0D
            HubbleParam=0.0D

            bytesleft=256-6*4 - 6*8 - 8 - 8 - 2*4 - 6*4 - 2*4 - 4*8
            la=bytarr(bytesleft)
            readu,myfile,npart,massarr,time,redshift,flag_sfr,flag_feedback,partTotal, $
                         flag_cooling,num_files,BoxSize,Omega0,OmegaLambda,HubbleParam,la
            x = { head , npart:npart,$
                         massarr:massarr,$
                         time:time,$
                         redshift:redshift,$
                         flag_sfr:flag_sfr,$
                         flag_feedback:flag_feedback,$
		         partTotal:partTotal,$
		         flag_cooling:flag_cooling,$
                         num_files:num_files,$
                         BoxSize:BoxSize,$
                         Omega0:Omega0,$
                         OmegaLambda:OmegaLambda,$
                         HubbleParam:HubbleParam,$
                         la:la}
                IF IS_DEF(ignorema) THEN BEGIN
print,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
print,"!!!!!!!!!!!!!!!!! Patching !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
print,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
            x.massarr(*)=0
		  IF NOT(IS_DEF(quiet)) THEN print,"Setting massarr members to zero"
print,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
print,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
                END

         END ELSE BEGIN
            status=1
            IF strcmp(label,"MASS") THEN BEGIN
               readnew,name,h,'HEAD',ignorema=ignorema,double=double
               ntotm=TOTAL(h.npart,/int)
               If IS_DEF(double) THEN BEGIN
                  ninput=(bl-8)/8
                  x=dblarr(ntotm)
                  xr=dblarr(ninput)
               END ELSE BEGIN
                  ninput=(bl-8)/4
                  x=fltarr(ntotm)
                  xr=fltarr(ninput)
               END
               readu,myfile,xr
               st=0L
               str=0L
               FOR i=0,5 DO BEGIN
                  IF h.npart(i) GT 0 THEN BEGIN
;                     print,'filling ',st,st+h.npart(i)-1
                     IF h.massarr(i) GT 0 THEN BEGIN
                        x(st:st+h.npart(i)-1)=h.massarr(i)
                        st=st+h.npart(i)
                     END ELSE BEGIN
                        x(st:st+h.npart(i)-1)=xr(str:str+h.npart(i)-1)
                        st=st+h.npart(i)
                        str=str+h.npart(i)
                     END
                  END
               END
               IF IS_DEF(parttype) THEN BEGIN
                  IF h.npart(parttype) GT 0 THEN valid=1 ELSE valid=0
               END ELSE valid=1
            END ELSE BEGIN
               IF IS_DEF(parttype) THEN BEGIN
                  readnew,name,h,'HEAD',ignorema=ignorema,double=double
                  IF h.npart(parttype) GT 0 THEN valid=1 ELSE valid=0
               END ELSE valid=1
	       IF IS_DEF(double) THEN BEGIN
                  IF strcmp(type,"FLOATN") THEN ninput=(bl-8)/8/ndim $
                                           ELSE ninput=(bl-8)/8
               END ELSE BEGIN
                  IF strcmp(type,"FLOATN") THEN ninput=(bl-8)/4/ndim $
                                           ELSE ninput=(bl-8)/4
               END
               IF valid EQ 1 THEN BEGIN
                  If IS_DEF(double) THEN BEGIN
                     IF strcmp(type,"FLOAT")  THEN x=dblarr(ninput) 
                     IF strcmp(type,"FLOATN") THEN x=dblarr(ndim,ninput)
                     IF strcmp(type,"LONG")   THEN x=ulon64arr(ninput)
                  END ELSE BEGIN 
                     IF strcmp(type,"FLOAT")  THEN x=fltarr(ninput) 
                     IF strcmp(type,"FLOATN") THEN x=fltarr(ndim,ninput)
                     IF strcmp(type,"LONG")   THEN x=ulonarr(ninput) 
                  END
               END
               IF valid EQ 1 THEN readu,myfile,x
            END
            IF valid EQ 1 THEN BEGIN
               IF IS_DEF(parttype) THEN BEGIN
                  IF ninput NE h.npart(parttype) THEN BEGIN
                     pstart=0L
                     FOR i=0,parttype-1 DO pstart=pstart+h.npart(i)*is_present(i)
                     IF NOT(IS_DEF(quiet)) THEN $
                        print,'Restricting return value to',pstart,pstart+h.npart(parttype)-1
                     IF strcmp(type,"FLOATN") THEN $
                        x=x(*,pstart:pstart+h.npart(parttype)-1) $
                     ELSE $
                        x=x(pstart:pstart+h.npart(parttype)-1) 
                  END
               END
            END ELSE BEGIN
               status=-1
               IF IS_DEF(parttype) THEN $
                  print,'File contains no "',type,'" for particle oy type',parttype,' !' $
               ELSE $
                  print,'Unknown type "',type,'" !'
            END
         END
      END

      close,myfile
   END
   free_lun,myfile
 END
   

pro extract_gas,snr
   IF NOT(IS_DEF(snr)) THEN snr='045'

   FOR i=0,31 DO BEGIN
      ext='.'+STRING(i,format='(i2)')
      IF i LT 10 THEN ext='.'+STRING(i,format='(i1)')

      iname='snapdir_'+snr+'/snap_'+snr+ext
      name='snap_'+snr+ext

      readnew,iname,h,'HEAD'
      h.npart(1:5)=0
      h.partTotal(1:5)=0
      write_head,name,h

      readnew,iname,x,'POS',part=0
      add_block,name,x,'POS'

      readnew,iname,x,'VEL',part=0
      add_block,name,x,'VEL'

      readnew,iname,x,'MASS',part=0
      add_block,name,x,'MASS'

      readnew,iname,x,'RHO'
      add_block,name,x,'RHO'

      readnew,iname,x,'HSML'
      add_block,name,x,'HSML'

      readnew,iname,x,'TEMP'
      add_block,name,x,'TEMP'

      readnew,iname,x,'MHI '
      add_block,name,x,'MHI '

   END

end

pro extract_all
;   add_temp
   snr=['031','032','034','035','036','038','039','041','042','044']
;   snr=['045']
   ss=size(snr)

  for i=0,ss(1)-1 do begin
     extract_gas,snr(i)
  end
end
