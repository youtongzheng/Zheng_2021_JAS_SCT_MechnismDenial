PRO RUNMLM_BW97_EXP1_BASIC_FXDLHF_20201027

  Path = 'D:\Results\Onetime\20201027_revisit_BW97\highres\'

  ;;0. Set consts and coefficients
  ;0.1. consts
  g         = 9.80665d ; gravity const.
  BOLTZMAN  = 1.380658e-23;                 Boltzman constant (J/K)
  AVOGADRO  = .602214199e24;                Avogadro constant (1/mol)
  MD        = 28.9644e-3;                   molar mass dry air (kg/mol)
  MV        = 18.0153e-3;                   molae mass water vapor (kg/mol)
  r_v       = (AVOGADRO)*(BOLTZMAN)/(MV);   gas constant for water vapor (J/(kg-K))
  r_d       = (AVOGADRO)*(BOLTZMAN)/(MD);   gas constant for dry air (J/(kg-K))
  cp        = 7./2d*(r_d);                  specific heat of air (J/(kg-K))
  L        = 2.5008e6;                      latent heat of vaporization (J/kg)
  rho_ref       = 1.d ;                     density of air (kg/m3)
  rho_lw        = 1000.;                    density of liquid water (kg/m3)
  psi = 0.61 ;                              for transferring from theta to theta_v (theta_v = theta + psi*q_v - q_l)

  ;0.2. empirical coefficients
  mu = 0.93; from BW97
  epsilon = 0.1; from BW97
  sigma = 0.35; from BW97
  betaa = 0.5; from BW97
  K_C = 2.44e10 ; equation (12) in Wood07
  z_evap = 475. ; equation (13) in  Wood07
  T_ref = epsilon*L/cp

  ;0.3. resolution
  ;0.3.1. vertical grids
  deltaZ = 1.
  z = deltaZ*FINDGEN(3001.)
  nz = N_ELEMENTS(z)

  ;0.3.2. time steps
  deltat = 1200d ;s
  nday = 7.
  nt = (nday*24.)*60.*60./deltat + 1

  ;time (s)
  t_plot = FLTARR(nt)
  ;time (hr)
  t_hr_plot = FINDGEN(nt)/3.

  ;time (local hr)
  t_hr_local_plot = FLTARR(nt)
  ind = WHERE(t_hr_plot MOD 24 EQ 0)
  t_hr_local_plot[ind[0]:ind[1] - 1] = t_hr_plot[ind[0]:ind[1] - 1]
  t_hr_local_plot[ind[1]:ind[2] - 1] = t_hr_plot[ind[1]:ind[2] - 1] - 24.
  t_hr_local_plot[ind[2]:ind[3] - 1] = t_hr_plot[ind[2]:ind[3] - 1] - 24.*2
  t_hr_local_plot[ind[3]:ind[4] - 1] = t_hr_plot[ind[3]:ind[4] - 1] - 24.*3
  t_hr_local_plot[ind[4]:ind[5] - 1] = t_hr_plot[ind[4]:ind[5] - 1] - 24.*4
  t_hr_local_plot[ind[5]:ind[6] - 1] = t_hr_plot[ind[5]:ind[6] - 1] - 24.*5

  ;specified radiative cooling
  F_rad0 = MAKE_ARRAY(nt, /FLOAT, VALUE = 60.)

  ;0.3.3. specify a series of runs
  ;0.3.4. Generate array for holding output
  na = 4

  ;0.3.3.1. profiles
  ;flux variables
  E = FLTARR(na, nt, nz)
  W = FLTARR(na, nt, nz)
  h_flux = FLTARR(na, nt, nz)
  qt_flux = FLTARR(na, nt, nz)
  sv_flux = FLTARR(na,nt, nz)
  b_flux = FLTARR(na,nt, nz)

  F_rad = FLTARR(na,nt, nz)
  F_precip = FLTARR(na,nt, nz)

  ;state variables
  t_profile = FLTARR(na,nt, nz)
  pres_profile = FLTARR(na,nt, nz)
  s_profile = FLTARR(na,nt, nz)
  h_profile = FLTARR(na,nt, nz)
  qv_profile = FLTARR(na,nt, nz)
  ql_profile = FLTARR(na,nt, nz)
  sv_profile = FLTARR(na,nt, nz)
  svl_profile = FLTARR(na,nt, nz)
  theta_profile = FLTARR(na,nt, nz)
  thetae_profile = FLTARR(na,nt, nz)
  thetav_profile = FLTARR(na,nt, nz)
  tv_profile = FLTARR(na,nt, nz)
  rho_profile = FLTARR(na,nt, nz)

  ;0.3.3.2. bulk parameters
  zi_plot = FLTARR(na, nt)
  zb_plot = FLTARR(na, nt)
  SST_plot = FLTARR(na, nt)
  LHF_plot = FLTARR(na, nt)
  SHF_plot = FLTARR(na, nt)
  we_plot = FLTARR(na, nt)
  evap_plot = FLTARR(na, nt)
  CTRC_plot = FLTARR(na, nt)
  BIR_plot = MAKE_ARRAY(na, nt, /FLOAT, VALUE = 999.)
  Bfx_mean_plot = FLTARR(na, nt)
  F_precip_cb_plot = FLTARR(na, nt)
  F_precip_s_plot = FLTARR(na, nt)
  LWP_plot = FLTARR(na, nt)

  svl_M_plot = FLTARR(na, nt)
  A_plot = FLTARR(na, nt)
  delta_inv_svl_plot = FLTARR(na, nt)
  delta_inv_h_plot = FLTARR(na, nt)
  delta_inv_qt_plot = FLTARR(na, nt)

  Ent_CLEB_plot = FLTARR(na, nt)
  Diab_CLEB_plot = FLTARR(na, nt)
  Rad_CLEB_plot = FLTARR(na, nt)
  Prec_CLEB_plot = FLTARR(na, nt)
  Res_CLEB_plot = FLTARR(na, nt)
  Store_CLEB_plot = FLTARR(na, nt)

  ;0.4. Simulation-specific boundary and initial forcing
  ;0.4.1. unvaried forcing
  Nd = 50.e6;     m-3
  V  = 7.1;      surface wind speed, m/s
  ps = 1022.;      surface pressure, hPa
  p_ref = 1000.;  reference  pressure, hPa

  ;0.4.2. initial forcing
  zi_0 = 413.;           boundary layer top (m)
  s_plus0 = 298.*cp ; inversion-layer top static energy
  SST_0 = 285.;          sea surface temp (K)
  s_M_0 = 284.8*cp;        static energy
  qt_M_0 = 7.69/1000.;    total water mixing ratio (g kg-1 -> kg/kg)
  D00  = 4.e-6;     s-1
  qt_plus0 = 5./1000.;  g/kg -> kg/kg
  dSST_dt = 1.5

  ;1. run simulations
  FOR ia = 0, na - 1 DO BEGIN
    if ia eq 3 then s_plus0 = 290.*cp
    
    zi = zi_0
    s_M = s_M_0
    qt_M = qt_M_0
    h_M = s_M + L*qt_M
    D = D00
    qt_plus = qt_plus0

    c0 = FLTARR(nz)
    c1 = FLTARR(nz)
    
    FOR t = 0L, nday*24*60*60L, deltat DO BEGIN
      it = t/deltat

      ;(D #0) Update boundary forcing and other state parameters at reference heights
      SST = SST_0 + t*(dSST_dt/(24.*60.*60.))
      s_plus = s_plus0 + (3.36*(zi - zi_0)/1000.)*cp
      h_plus = s_plus + L*qt_plus
      sv_plus =  h_plus - mu*L*qt_plus
      svl_plus =  sv_plus

      s_M = h_M - L*qt_M;     ML static energy
      T2 = s_M/cp;          near-surface air temperature (K)
      theta_M = T2*(p_ref/ps)^(r_d/cp);       ML potential temperature (K)
      svl_M = h_M - mu*L*qt_M

      rho_ref = 100.*967/(r_d*(SST - 4.5)*(1. + 0.608*qt_M));     surface density
      qt_s = MIXR_SAT(SST - 273.15, 1022)/1000.
      h_s = SST*cp + L*qt_s
      RH_s = MIXR2RH(1000.*qt_M,ps,T2);   surface relative humidity (%)
      zb_ori = ROMPLCL(ps*100.,T2,RH_s/100.);   lifting condensation level (m)
      zb = deltaZ*ROUND(zb_ori)
      ; ---------- END (D #0) ----------

      ; (D #1) Diagnose vertical profiles

      ;---sub-cloud layer---
      ind = WHERE(z LT zb)
      z_ind = z[ind]
      pres_profile[ia, it, ind] = ps/EXP((g*z[ind])/(r_d*theta_M))

      s_profile[ia, it, ind] = s_M
      h_profile[ia, it, ind] = h_M
      qv_profile[ia, it, ind] = qt_M
      ql_profile[ia, it, ind] = 0.
      ;sv_profile[ia, it, ind] = h_M - mu*L*qt_M
      ;svl_profile[ia, it, ind] = h_M - mu*L*qt_M

      t_profile[ia, it, ind] = (s_M - g*z[ind])/cp
      tv_profile[ia, it, ind] =  t_profile[ia, it, ind]*(1. + 0.608*qv_profile[ia, it, ind])

      sv_profile[ia, it, ind] = cp*tv_profile[ia, it, ind]+ g*z[ind]
      svl_profile[ia, it, ind] = sv_profile[ia, it, ind]

      F_rad[ia, it,ind] = 0.; W/m2

      ;---cloud layer---
      ;diagnose temperature at cloud base
      T_cb = (s_M - g*zb)/cp; cloud-base temperature
      p_cb = ps/EXP((g*zb)/(r_d*theta_M)); cloud-base pressure

      ql_cld = PSEUDOADIABT_10_highres(zb, P_cb, T_cb, $
        H_output = IZV,$
        Pres_output = PV,$
        Temp_output = TV)

      ind = WHERE(z GE zb AND z LE zi, count)

      ql_profile[ia, it, ind] = ql_cld[ind - ind[0]]/1000.
      qv_profile[ia, it, ind] = qt_M - ql_profile[ia, it, ind]

      s_profile[ia, it, ind] = h_M - L*qv_profile[ia, it, ind]
      ;s_profile[ia, it, ind] = s_M + (1.-(1. + 0.608)*0.1)*L*ql_cld[ind - ind[0]]/1000.
      t_profile[ia, it, ind] = (s_profile[ia, it, ind] - g*z[ind])/cp
      tv_profile[ia, it, ind] =  t_profile[ia, it, ind]*(1. + 0.608*qv_profile[ia, it, ind] - ql_profile[ia, it, ind] )
      sv_profile[ia, it, ind] = cp*tv_profile[ia, it, ind]+ g*z[ind]
      svl_profile[ia, it, ind] = sv_profile[ia, it, ind] - (1.-(1. + 0.608)*0.1)*L*ql_profile[ia, it, ind]

      pres_profile[ia, it, ind] = PV[ind - ind[0]]

      F_rad[ia, it,ind] = 0.

      ;diagnose LWP
      ql_cld_ind = ql_cld[ind - ind[0]]/1000.
      z_ind = z[ind]
      LWP = 0.5*2.*((zi - zb)/1000.)^2.

      ;--Turton and Nicholl, 1987 parameterization scheme. (Eq.(17))
      LWC_zi = ql_profile[ia, it, ind[count - 1]] ;kg/kg
      pres_zi = pres_profile[ia, it, ind[count - 1]]
      Tv_zi = tv_profile[ia, it, ind[count - 1]]
      rho_zi = pres_zi*100./(287.*Tv_zi)

      LWC_zi = LWC_zi*rho_zi ; kg/kg to kg/m3

      rv_star = ((LWC_zi/Nd)/((4./3.)*3.1415926*rho_lw))^(1./3.)
      rv_star = rv_star*(1.e6)

      IF rv_star GE 10. THEN BEGIN
        g_rvstar = 1.
      ENDIF ELSE BEGIN
        g_rvstar = (rv_star/10.)^3.
      ENDELSE

      F_precip_cb = -3.*g_rvstar*SQRT(LWP)*(1.e-5)
      F_precip[ia, it, ind] = F_precip_cb

      ;modulate subcloud precipitation
      ind = WHERE(z LT zb)
      z_ind = z[ind]

      decrs_rate_precip = 0.86e-5 ; ms-1km-1
      F_precip[ia, it, ind] = F_precip_cb + (zb - z_ind)*decrs_rate_precip/1000.

      F_precip_s = F_precip[ia, it, 0]

      ;---above-cloud layer---
      ind = WHERE(z GT zi)
      s_profile[ia, it,ind] = (s_plus0/cp + (3.36/1000.)*(z[ind] - z[ind[0]]))*cp
      qv_profile[ia, it,ind] = qt_plus ;
      ql_profile[ia, it, ind] = 0.
      t_profile[ia, it, ind] = (s_profile[ia, it, ind] - g*z[ind])/cp
      tv_profile[ia, it, ind] =  t_profile[ia, it, ind]*(1. + 0.608*qv_profile[ia, it, ind] - ql_profile[ia, it, ind])
      sv_profile[ia, it, ind] = cp*tv_profile[ia, it, ind]+ g*z[ind]
      svl_profile[ia, it, ind] = sv_profile[ia, it, ind] - (1.-(1. + 0.608)*0.1)*L*ql_profile[ia, it, ind]

      F_rad[ia, it,ind] = F_rad0[it]
      F_rad_plus = F_rad0[it]

      ; ---------- END (D #1) ----------

      ; (D #2) Diagnose surface fluxes
      Ct = 0.001*(1. + 0.07*V)
      E_0 = Ct*V*(h_s - h_M) + F_rad[ia, it,0]/rho_ref
      W_0 = Ct*V*(qt_s - qt_M) + F_precip[ia, it,0]

      SHF = 1216.*Ct*V*(SST - T2)
      LHF = L*Ct*V*(qt_s - qt_M)

      IF it EQ 0 THEN BEGIN
        h_s0 = h_s
        h_M0 = h_M
        qt_s0 = qt_s
        qt_M0 = qt_M
        SST0 = SST
        T20 = T2
      ENDIF

      IF ia EQ 1 OR ia EQ 3 THEN BEGIN
        E_0 = Ct*V*(h_s0 - h_M0) + F_rad[ia, it,0]/rho_ref
        W_0 = Ct*V*(qt_s0 - qt_M0) + F_precip[ia, it,0]

        SHF = 1216.*Ct*V*(SST0 - T20)
        LHF = L*Ct*V*(qt_s0 - qt_M0)
      ENDIF

      ; ---------- END (D #2) ----------

      ; (D #3) Diagnose entrainment efficiency A based on Eq.13 in BW97
      ;cloud top
      ind = WHERE(z LE zi, count)
      ind2 = ind[count - 1]

      tmp2 = t_profile[ia, it,ind2]
      pres2 = pres_profile[ia, it,ind2]
      qv2 = qv_profile[ia, it,ind2]
      ql2 = qt_M - qv2
      qt2 = qv2 + ql2

      theta2 = tmp2*(p_ref/pres2)^(r_d/cp)
      thetae2 = theta2 + L*qv2/cp
      thetav2 = theta2*(1. + psi*qv2 - ql2)

      ;cloud top + 1
      ind1 = ind2 + 1
      tmp1 = t_profile[ia, it,ind1]
      pres1 = pres2 ; approximation
      qv1 = qv_profile[ia, it,ind1]
      ql1 = 0.
      qt1 = qv1 + ql1

      theta1 = tmp1*(p_ref/pres1)^(r_d/cp)
      thetae1 = theta1 + L*qv1/cp
      thetav1 = theta1*(1. + psi*qv1 - ql1)

      ;calculate the delta_m
      epsilonm = 0.01*FINDGEN(101)

      delta_thetav_lim = epsilonm*(thetae1 - thetae2) $
        - epsilonm*(qt1 - qt2)*(L/cp - psi*thetae2) $
        + (ql1*epsilonm - ql2)*(L/cp - (1. + psi)*thetae2)

      delta_m = 2.*INT_TABULATED(epsilonm, delta_thetav_lim)

      delta_thetav = thetav1 - thetav2

      A = 0.2*(1. + 60.*(1. - delta_m/delta_thetav))
      If ia eq 2 then A = 1.0983624345811129
      
      ; ---------- END (D #3) ----------

      ; (D #4) Diagnose entrainment
      sv_M = h_M - mu*L*qt_M
      sv_0 = sv_M

      alpha = 2.5*A*sv_0/(g*zi*(sv_plus - sv_M))

      d0 = (1. - z/zi)*E_0 + (z/zi)*(F_rad_plus/rho_ref) - F_rad[ia, it, *]/rho_ref
      e0 = (1. - z/zi)*W_0 - F_precip[ia, it, *]

      d1 = -(z/zi)*(h_plus - h_M)
      e1 = -(z/zi)*(qt_plus - qt_M)

      ind = WHERE(z LT zb)
      c0[ind] = (d0[ind] - mu*L*e0[ind])*(g/sv_0)
      c1[ind] = (d1[ind] - mu*L*e1[ind])*(g/sv_0)

      ind = WHERE(z GE zb AND z LE zi)
      c0[ind] = (betaa*d0[ind] - epsilon*L*e0[ind])*(g/sv_0)
      c1[ind] = (betaa*d1[ind] - epsilon*L*e1[ind])*(g/sv_0)

      ind = WHERE(z LE zi)
      c0_ind = c0[ind]
      c1_ind = c1[ind]
      z_ind = z[ind]

      I0 = INT_TABULATED(z_ind, c0_ind)
      I1 = INT_TABULATED(z_ind, c1_ind)

      we = alpha*I0/(1. - alpha*I1)
      ; ---------- END (D #4) ----------

      ; (D #5) Diagnose fluxes
      E_zi = -we*(h_plus - h_M) + F_rad_plus/rho_ref
      W_zi = -we*(qt_plus - qt_M)

      ;buoyancy flux
      b_flux[ia, it, *] = c0 + c1*we
      sv_flux[ia, it, *] = b_flux[ia, it, *]/(g/sv_0)

      h_flux[ia, it, *] = d0 + d1*we
      qt_flux[ia, it, *] = e0 + e1*we

      E[ia, it, *] = h_flux[ia, it, *] + F_rad[ia, it, *]
      W[ia, it, *] = qt_flux[ia, it, *] + F_precip[ia, it, *]

      ind = WHERE(z LE zi)
      bflux_mean = MEAN(sv_flux[ia, it, ind])
      we_check = 2.5*A*bflux_mean/(sv_plus - sv_M)
      ; ---------- END (D #5) ----------

      ; (D #6) Diagnose other useful variables
      Ent_CLEB = rho_ref*we*(svl_plus - svl_M)
      Rad_CLEB = -F_rad0[it]
      Prec_CLEB = -mu*L*F_precip_cb
      Diab_CLEB = Rad_CLEB + Prec_CLEB


      tmp = sv_flux[ia, it, *]
      ind = WHERE(z LT zb, count)
      svl_flux_cb = tmp[ind[count - 1]]
      Res_CLEB = -rho_ref*svl_flux_cb

      ; ---------- END (D #6) ----------

      ; (D #8) Diagnose decoupling
      ind = WHERE(z LE zi)
      sv_flux_temp = sv_flux[ia, it, ind]
      z_temp = z[ind]

      ind_neg = WHERE(sv_flux_temp LT 0., count)
      sv_flux_temp_neg = sv_flux_temp[ind_neg]
      z_temp_neg = z_temp[ind_neg]

      sv_flux_temp[ind_neg] = 0.

      IF count LE 1 THEN BEGIN
        BIR_plot[ia, it] = 0.
      ENDIF ELSE BEGIN
        numerator_BIR =  INT_TABULATED(z_temp_neg, sv_flux_temp_neg)
        denominator_BIR = INT_TABULATED(z_temp, sv_flux_temp)
        BIR_plot[ia, it] = -numerator_BIR/denominator_BIR
      ENDELSE

      ; ---------- END (D #8) ----------


      ;document parameters
      it = t/deltat
      t_plot[it] = t
      zi_plot[ia, it] = zi   ;m
      zb_plot[ia, it] = zb_ori   ;m
      SST_plot[ia, it] = SST   ;m
      LHF_plot[ia, it] = LHF   ;m
      SHF_plot[ia, it] = SHF   ;m
      we_plot[ia, it] = we   ;m/s
      evap_plot[ia, it] = delta_m/delta_thetav   ;
      CTRC_plot[ia, it] = -999.   ;m
      Bfx_mean_plot[ia, it] = bflux_mean   ;Wm-2
      F_precip_cb_plot[ia, it] = rho_ref*L*F_precip_cb   ;Wm-2
      F_precip_s_plot[ia, it] = rho_ref*L*F_precip_s   ;Wm-2
      LWP_plot[ia, it] = LWP
      A_plot[ia, it] = A
      svl_M_plot[ia, it] = svl_M

      Ent_CLEB_plot[ia, it] = Ent_CLEB
      Diab_CLEB_plot[ia, it] = Diab_CLEB
      Res_CLEB_plot[ia, it] = Res_CLEB
      Rad_CLEB_plot[ia, it] = Rad_CLEB
      Prec_CLEB_plot[ia, it] = Prec_CLEB

      delta_inv_svl_plot[ia, it] = svl_plus - svl_M
      delta_inv_h_plot[ia, it] = h_plus - h_M
      delta_inv_qt_plot[ia, it] = qt_plus - qt_M

      z_plot = z
      
      IF BIR_plot[ia, it] GT 0.25 THEN BREAK
      
      ;Advance prognostic variables in time
      zi_new = zi + deltat*(we -D*zi)
      h_M_new = h_m + (deltat/zi)*(-(E_zi - E_0))
      qt_m_new = qt_m + (deltat/zi)*(-(W_zi - W_0))

      zi = zi_new
      h_m = h_m_new
      qt_m = qt_m_new
      
      PRINT, t/3600., A, 1000.*we, sv_flux[ia, it, 0], zi, zb, h_M/cp, 1000.*qt_M, L*F_precip_s, BIR_plot[ia, it]
      PRINT, ' '
    ENDFOR

    ;2. Compute the tendency terms
    tmp = svl_M_plot[ia,*]
    ntmp = N_ELEMENTS(tmp)

    dtmp = (tmp[2: ntmp - 1] - tmp[0:ntmp - 3])/(2.*deltat)

    Store_CLEB_plot[ia,1:ntmp - 2] =  (zi_plot[ia, 1:ntmp - 2] - zb_plot[ia, 1:ntmp - 2])*rho_ref*dtmp
  ENDFOR


  ;3. Output the data
  filename = Path + 'MLMOutput_basic_fxdlhf_20210330_highres.nc'

  ; Create a new NetCDF file with the filename inquire.nc:
  id = NCDF_CREATE(filename, /CLOBBER)
  ; Fill the file with default values:
  NCDF_CONTROL, id, /FILL

  ; Define dimensions
  xid = NCDF_DIMDEF(id, 'x', na)
  yid = NCDF_DIMDEF(id, 'y', nt)
  zid = NCDF_DIMDEF(id, 'z', nz)

  ; Define variables:
  tid = NCDF_VARDEF(id, 't', [yid], /FLOAT)
  thrid = NCDF_VARDEF(id, 't_hr', [yid], /FLOAT)
  thrlocalid = NCDF_VARDEF(id, 't_hr_local', [yid], /FLOAT)
  frad0id = NCDF_VARDEF(id, 'F_rad0', [yid], /FLOAT)

  zzid = NCDF_VARDEF(id, 'z', [zid], /FLOAT)

  ziid = NCDF_VARDEF(id, 'zi', [xid, yid], /FLOAT)
  zbid = NCDF_VARDEF(id, 'zb', [xid, yid], /FLOAT)
  sstid = NCDF_VARDEF(id, 'sst', [xid, yid], /FLOAT)
  LHFid = NCDF_VARDEF(id, 'LHF', [xid, yid], /FLOAT)
  SHFid = NCDF_VARDEF(id, 'SHF', [xid, yid], /FLOAT)
  weid = NCDF_VARDEF(id, 'we', [xid, yid], /FLOAT)
  evapid = NCDF_VARDEF(id, 'evap', [xid, yid], /FLOAT)
  Bfxmeanid = NCDF_VARDEF(id, 'Bfx_mean', [xid, yid], /FLOAT)
  Fprecipcbid = NCDF_VARDEF(id, 'F_precip_cb', [xid, yid], /FLOAT)
  Fprecipsid = NCDF_VARDEF(id, 'F_precip_s', [xid, yid], /FLOAT)
  LWPid = NCDF_VARDEF(id, 'LWP', [xid, yid], /FLOAT)
  BIRid = NCDF_VARDEF(id, 'BIR', [xid, yid], /FLOAT)
  Aid = NCDF_VARDEF(id, 'A', [xid, yid], /FLOAT)
  deltainvsvlid = NCDF_VARDEF(id, 'delta_inv_svl', [xid, yid], /FLOAT)
  deltainvhid = NCDF_VARDEF(id, 'delta_inv_h', [xid, yid], /FLOAT)
  deltainvqtid = NCDF_VARDEF(id, 'delta_inv_qt', [xid, yid], /FLOAT)
  EntCLEBid = NCDF_VARDEF(id, 'Ent_CLEB', [xid, yid], /FLOAT)
  DiabCLEBid = NCDF_VARDEF(id, 'Diab_CLEB', [xid, yid], /FLOAT)
  RadCLEBid = NCDF_VARDEF(id, 'Rad_CLEB', [xid, yid], /FLOAT)
  PrecCLEBid = NCDF_VARDEF(id, 'Prec_CLEB', [xid, yid], /FLOAT)
  ResCLEBid = NCDF_VARDEF(id, 'Res_CLEB', [xid, yid], /FLOAT)
  StoreCLEBid = NCDF_VARDEF(id, 'Store_CLEB', [xid, yid], /FLOAT)

  Eid = NCDF_VARDEF(id, 'E', [xid,yid,zid], /FLOAT)
  Wid = NCDF_VARDEF(id, 'W', [xid,yid,zid], /FLOAT)
  hfluxid = NCDF_VARDEF(id, 'h_flux', [xid,yid,zid], /FLOAT)
  qtfluxid = NCDF_VARDEF(id, 'qt_flux', [xid,yid,zid], /FLOAT)
  svfluxid = NCDF_VARDEF(id, 'sv_flux', [xid,yid,zid], /FLOAT)
  bfluxid = NCDF_VARDEF(id, 'b_flux', [xid,yid,zid], /FLOAT)
  Fradid = NCDF_VARDEF(id, 'F_rad', [xid,yid,zid], /FLOAT)
  Fprecipid = NCDF_VARDEF(id, 'F_precip', [xid,yid,zid], /FLOAT)
  tempid = NCDF_VARDEF(id, 'temp', [xid,yid,zid], /FLOAT)
  tempvid = NCDF_VARDEF(id, 'tempv', [xid,yid,zid], /FLOAT)
  presid = NCDF_VARDEF(id, 'pres', [xid,yid,zid], /FLOAT)
  sid = NCDF_VARDEF(id, 's', [xid,yid,zid], /FLOAT)
  hid = NCDF_VARDEF(id, 'h', [xid,yid,zid], /FLOAT)
  qvid = NCDF_VARDEF(id, 'qv', [xid,yid,zid], /FLOAT)
  qlid = NCDF_VARDEF(id, 'ql', [xid,yid,zid], /FLOAT)
  svid = NCDF_VARDEF(id, 'sv', [xid,yid,zid], /FLOAT)
  svlid = NCDF_VARDEF(id, 'svl', [xid,yid,zid], /FLOAT)

  ; Put file in data mode:
  NCDF_CONTROL, id, /ENDEF

  ; Input data:
  NCDF_VARPUT, id, zzid, z_plot
  NCDF_VARPUT, id, tid, t_plot
  NCDF_VARPUT, id, thrid, t_hr_plot
  NCDF_VARPUT, id, thrlocalid, t_hr_local_plot
  NCDF_VARPUT, id, frad0id, F_rad0

  FOR j=0,nt - 1 DO NCDF_VARPUT, id, ziid, $
    REFORM(zi_plot[*,j]), OFFSET=[0,j]

  FOR j=0,nt - 1 DO NCDF_VARPUT, id, zbid, $
    REFORM(zb_plot[*,j]), OFFSET=[0,j]

  FOR j=0,nt - 1 DO NCDF_VARPUT, id, sstid, $
    REFORM(sst_plot[*,j]), OFFSET=[0,j]

  FOR j=0,nt - 1 DO NCDF_VARPUT, id, LHFid, $
    REFORM(LHF_plot[*,j]), OFFSET=[0,j]

  FOR j=0,nt - 1 DO NCDF_VARPUT, id, SHFid, $
    REFORM(SHF_plot[*,j]), OFFSET=[0,j]

  FOR j=0,nt - 1 DO NCDF_VARPUT, id, weid, $
    REFORM(we_plot[*,j]), OFFSET=[0,j]

  FOR j=0,nt - 1 DO NCDF_VARPUT, id, evapid, $
    REFORM(evap_plot[*,j]), OFFSET=[0,j]

  FOR j=0,nt - 1 DO NCDF_VARPUT, id, Bfxmeanid, $
    REFORM(Bfx_mean_plot[*,j]), OFFSET=[0,j]

  FOR j=0,nt - 1 DO NCDF_VARPUT, id, Fprecipcbid, $
    REFORM(F_precip_cb_plot[*,j]), OFFSET=[0,j]

  FOR j=0,nt - 1 DO NCDF_VARPUT, id, Fprecipsid, $
    REFORM(F_precip_s_plot[*,j]), OFFSET=[0,j]

  FOR j=0,nt - 1 DO NCDF_VARPUT, id, LWPid, $
    REFORM(LWP_plot[*,j]), OFFSET=[0,j]

  FOR j=0,nt - 1 DO NCDF_VARPUT, id, BIRid, $
    REFORM(BIR_plot[*,j]), OFFSET=[0,j]

  FOR j=0,nt - 1 DO NCDF_VARPUT, id, Aid, $
    REFORM(A_plot[*,j]), OFFSET=[0,j]

  FOR j=0,nt - 1 DO NCDF_VARPUT, id, deltainvsvlid, $
    REFORM(delta_inv_svl_plot[*,j]), OFFSET=[0,j]

  FOR j=0,nt - 1 DO NCDF_VARPUT, id, deltainvhid, $
    REFORM(delta_inv_h_plot[*,j]), OFFSET=[0,j]

  FOR j=0,nt - 1 DO NCDF_VARPUT, id, deltainvqtid, $
    REFORM(delta_inv_qt_plot[*,j]), OFFSET=[0,j]

  FOR j=0,nt - 1 DO NCDF_VARPUT, id, EntCLEBid, $
    REFORM(Ent_CLEB_plot[*,j]), OFFSET=[0,j]

  FOR j=0,nt - 1 DO NCDF_VARPUT, id, DiabCLEBid, $
    REFORM(Diab_CLEB_plot[*,j]), OFFSET=[0,j]

  FOR j=0,nt - 1 DO NCDF_VARPUT, id, RadCLEBid, $
    REFORM(Rad_CLEB_plot[*,j]), OFFSET=[0,j]

  FOR j=0,nt - 1 DO NCDF_VARPUT, id, PrecCLEBid, $
    REFORM(Prec_CLEB_plot[*,j]), OFFSET=[0,j]

  FOR j=0,nt - 1 DO NCDF_VARPUT, id, ResCLEBid, $
    REFORM(Res_CLEB_plot[*,j]), OFFSET=[0,j]

  FOR j=0,nt - 1 DO NCDF_VARPUT, id, StoreCLEBid, $
    REFORM(Store_CLEB_plot[*,j]), OFFSET=[0,j]

  FOR k = 0, nz - 1 DO BEGIN
    FOR j=0,nt - 1 DO NCDF_VARPUT, id, Eid, $
      REFORM(E[*,j,k]), OFFSET=[0,j,k]
  ENDFOR

  FOR k = 0, nz - 1 DO BEGIN
    FOR j=0,nt - 1 DO NCDF_VARPUT, id, Wid, $
      REFORM(W[*,j,k]), OFFSET=[0,j,k]
  ENDFOR

  FOR k = 0, nz - 1 DO BEGIN
    FOR j=0,nt - 1 DO NCDF_VARPUT, id, hfluxid, $
      REFORM(h_flux[*,j,k]), OFFSET=[0,j,k]
  ENDFOR

  FOR k = 0, nz - 1 DO BEGIN
    FOR j=0,nt - 1 DO NCDF_VARPUT, id, qtfluxid, $
      REFORM(qt_flux[*,j,k]), OFFSET=[0,j,k]
  ENDFOR

  FOR k = 0, nz - 1 DO BEGIN
    FOR j=0,nt - 1 DO NCDF_VARPUT, id, svfluxid, $
      REFORM(sv_flux[*,j,k]), OFFSET=[0,j,k]
  ENDFOR

  FOR k = 0, nz - 1 DO BEGIN
    FOR j=0,nt - 1 DO NCDF_VARPUT, id, bfluxid, $
      REFORM(b_flux[*,j,k]), OFFSET=[0,j,k]
  ENDFOR

  FOR k = 0, nz - 1 DO BEGIN
    FOR j=0,nt - 1 DO NCDF_VARPUT, id, Fradid, $
      REFORM(F_rad[*,j,k]), OFFSET=[0,j,k]
  ENDFOR

  FOR k = 0, nz - 1 DO BEGIN
    FOR j=0,nt - 1 DO NCDF_VARPUT, id, Fprecipid, $
      REFORM(F_precip[*,j,k]), OFFSET=[0,j,k]
  ENDFOR

  FOR k = 0, nz - 1 DO BEGIN
    FOR j=0,nt - 1 DO NCDF_VARPUT, id, tempid, $
      REFORM(t_profile[*,j,k]), OFFSET=[0,j,k]
  ENDFOR


  FOR k = 0, nz - 1 DO BEGIN
    FOR j=0,nt - 1 DO NCDF_VARPUT, id, tempvid, $
      REFORM(tv_profile[*,j,k]), OFFSET=[0,j,k]
  ENDFOR


  FOR k = 0, nz - 1 DO BEGIN
    FOR j=0,nt - 1 DO NCDF_VARPUT, id, presid, $
      REFORM(pres_profile[*,j,k]), OFFSET=[0,j,k]
  ENDFOR


  FOR k = 0, nz - 1 DO BEGIN
    FOR j=0,nt - 1 DO NCDF_VARPUT, id, sid, $
      REFORM(s_profile[*,j,k]), OFFSET=[0,j,k]
  ENDFOR


  FOR k = 0, nz - 1 DO BEGIN
    FOR j=0,nt - 1 DO NCDF_VARPUT, id, hid, $
      REFORM(h_profile[*,j,k]), OFFSET=[0,j,k]
  ENDFOR


  FOR k = 0, nz - 1 DO BEGIN
    FOR j=0,nt - 1 DO NCDF_VARPUT, id, qvid, $
      REFORM(qv_profile[*,j,k]), OFFSET=[0,j,k]
  ENDFOR


  FOR k = 0, nz - 1 DO BEGIN
    FOR j=0,nt - 1 DO NCDF_VARPUT, id, qlid, $
      REFORM(ql_profile[*,j,k]), OFFSET=[0,j,k]
  ENDFOR


  FOR k = 0, nz - 1 DO BEGIN
    FOR j=0,nt - 1 DO NCDF_VARPUT, id, svid, $
      REFORM(sv_profile[*,j,k]), OFFSET=[0,j,k]
  ENDFOR


  FOR k = 0, nz - 1 DO BEGIN
    FOR j=0,nt - 1 DO NCDF_VARPUT, id, svlid, $
      REFORM(svl_profile[*,j,k]), OFFSET=[0,j,k]
  ENDFOR
  NCDF_CLOSE, id ; Close the NetCDF file.

  PRINT,''

END

FUNCTION MIXINGRATIO, P, T
  ES = 6.112*EXP(17.67*(T - 273.)/(T - 29.5))
  WW = 0.622*ES/(P - ES)
  RETURN, WW
END

FUNCTION  PSEUDOADIABT_10_highres, Z0, $     ;Cloud base height, unit:m
  P0, $     ;Cloud base pressure, unit:hPa
  T0, $     ;Cloud base temperature, unit: K
  Ice = ice, $     ;The phase of cloud droplet
  H_output = IZV, $     ;Height profile, unit: m
  Pres_output = PV, $     ;Pressure profile, unit: hPa
  Temp_output = TV, $     ;Temperature profile, unit: C
  gm3_output = gm3


  ;Define constants
  drylap = 0.0098
  IR = 287.04
  ICP = 1004.67
  Epsilon = 0.622

  ideltaz = 1.

  ;  ;Define arrays
  ;  EE = FLTARR(1000)
  ;  IZV = FLTARR(1000)
  ;  PV = FLTARR(1000)
  ;  TV = FLTARR(1000)
  ;  WV = FLTARR(1000)
  ;  BWLOAD = FLTARR(1000)
  ;  gm3 = FLTARR(1000)

  EE = DBLARR(1001)
  IZV = DBLARR(1001)
  PV = DBLARR(1001)
  TV = DBLARR(1001)
  WV = DBLARR(1001)
  BWLOAD = DBLARR(1001)
  gm3 = DBLARR(1001)

  ;Start calculatoin
  FL = 0.
  IF KEYWORD_SET(ice) THEN FL = 333718.

  W0 = MIXINGRATIO(P0, T0); unit:kg/kg

  T = T0
  P = P0
  W = W0
  IZ = Z0
  ii = 0
  FOR i = 0, 1000 DO BEGIN
    IF IZ MOD 1 EQ 0 OR i EQ 0 THEN BEGIN
      IZV[ii] = IZ
      PV[ii] = P
      TV[ii] = T
      WV[ii] = W
      BWLOAD[ii] = W0 - WV[ii]
      BWLOAD[ii] = 1000.*BWLOAD[ii] ; from kg/kg to g/kg
      gm3[ii] = BWLOAD[ii]*(PV[ii]*100./(287.*TV[ii])); from g/kg to gm-3
      ii = ii + 1
    ENDIF
    L = 3133541. - 2317.*T + FL
    A = 1. + L*W/(IR*T)
    B = 1. + L*L*Epsilon*W/(IR*ICP*T*T)
    plap = drylap*A/B

    W = W + ideltaz*ICP/L*(plap - drylap/(1. + 0.6*W))
    EE = (W0 - W)/W0
    T = T - plap*ideltaz
    P = P*EXP(-ideltaz*9.8/(IR*T*(1+0.6*W)))
    IZ = IZ + ideltaz
  ENDFOR

  TV = TV - 273.15

  RETURN, BWLOAD

END