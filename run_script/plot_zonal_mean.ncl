;-------------------------------------------------------------------------
; Plot zonal averages of CondiDiag output.
; One QoI on each page. First panel shows the QoI at checkpoint CLDMIC01.
; Other panels show differences between two adjacent checkpoints.
; Before running this script, make sure you have remapped the E3SM output 
; to lat-lon grid, e.g., using the script remap_to_latlon.bash in the 
; same directory.
;-------------------------------------------------------------------------
; Hui Wan (2021-12-11)
;-------------------------------------------------------------------------

begin

  daytag = "05"
  itime  = 0

  data_path = "/compyfs/wanh895/scidac_scratch/macmic_subcycles_compy_F2010_ne30pg2_r05_oECv3_CondiDiag_for_MichaelBrunke/run/"
  data_file = "macmic_subcycles_compy_F2010_ne30pg2_r05_oECv3_CondiDiag_for_MichaelBrunke.eam.h0.0001-01-"+daytag+"-00000.latlon.nc"

  QoIs      = (/"CLOUD","FREQL","AWNC","CCN3","UP2","VP2","WP2","WP3","num_a1","num_a2","num_a3","num_a4"/)
  chkpts    = (/"CLDMIC01","CLDMIC02","CLDMIC03","CLDMIC04","CLDMIC05","CLDMIC06"/)
  cndprefix = "cnd01" 

  nqois   = dimsizes(QoIs)
  nchkpts = dimsizes(chkpts)

  plot = new(nchkpts,"graphic")
  File = addfile(data_path + data_file,"r")
  wks  = gsn_open_wks("pdf","day"+daytag)


  do iqoi = 0,nqois-1
     do ichkpt = 0,nchkpts-1

       res = True
       res@gsnDraw       = False
       res@gsnFrame      = False
       res@cnFillOn      = True 
       res@cnLinesOn     = False 
       res@cnFillPalette = "BlAqGrYeOrReVi200"
       res@tiYAxisString = "Norminal pressure (hPa)"
       res@trYMinF       = 600.
       res@tmYRMode      = "Automatic"

        varname = cndprefix+"_"+QoIs(iqoi)+"_"+chkpts(ichkpt)
        zonmean = dim_avg_Wrap(File->$varname$(itime,:,:,:))

        if (ichkpt.eq.0) then
            res@gsnLeftString = QoIs(iqoi)+" at "+chkpts(ichkpt)
            plot(ichkpt) = gsn_csm_pres_hgt(wks,zonmean,res)
        else
            res@gsnLeftString = "~F5~D~F~"+QoIs(iqoi)+", "+chkpts(ichkpt)+" - "+chkpts(ichkpt-1)
            diff = zonmean
            diff = zonmean - save
            symMinMaxPlt(diff,12,False,res)
            plot(ichkpt) = gsn_csm_pres_hgt(wks,diff,res)
        end if

        save = zonmean
        delete(res)
        print("plotting done for "+varname)
     end do

     resp = True
     resp@gsnMaximize = True
     gsn_panel(wks, plot,(/3,2/),resp)

  end do

 
end
