      program wsheet20
c-------------------------------------------------------------------------------
c---- Program WSHEET20.F
c---- Integral Equation Solution for Wide-band Thin Sheet EM modeling in whole
c---- space.
c---- Number of Sheets are two and the dip and strike of each sheet can be
c---- arbitrary.
c---- Both of the sheets and the whole space can be set polarizable with
c---- Cole - Cole parameter.
c---- LU Decomposition Version
c---- Author : Yoonho Song, Korea Institute of Geoscience and Mineral Resources
c-------------------------------------------------------------------------------
c---- input parameter file = 'WSHEET20.PAR'
c----      no_sheet : The number of sheet (1 or 2)
c----      na(i),nb(i) : No. of cell division of i-th sheet
c----      x_top(i),y_top(i),z_top(i) : x,y,z-coordinates of the center of
c----                  of the i-th sheet after rotation.
c----                  After rotation depth to the top of the sheet remains
c----                  unchanged
c----      aa(i),bb(i),alp,bet : strike length and depth extent (or x-length
c----                  when dip=0), strike (deg.) with respect to x-axis, and
c----                  dip angle (deg.)
c----                  Coordinates rotation is with respect to the center
c----                  position of the sheets
c----      thick(i) : thickness of i-th sheet (meter)
c----      res_l(i), res_h(i), res_tau(i), res_alp(i) : Cole-Cole parameter
c----                  resistivity(low),resistivity(high), time constant, and
c----                  alpha of i-th sheet
c----      eps_l(i), eps_h(i), eps_tau(i), eps_alp(i) : Cole-Cole parameter
c----                  of electric permittivity of the i-th sheet
c----      resn_l, resn_h, resn_tau, resn_alp : Cole-Cole parameter of the
c----                  resistivity of whole space
c----      epsn_l, epsn_h, epsn_tau, epsn_alp : Cole-Cole parameter of the
c----                  whole space
c----      freq : frequency (Hz)
c----      ns_type, ns_dir, angle : source type (0 for plane wave, 1 for
c----                  electric dipole, 2 for magnetic dipole)
c----                  source direction (1 for x, 2 for y, 3 for z if plane
c----                  wave 1 for TM 2 for TE)
c----                  angle is incidence for plane wave
c----      n_tx : No. of Tx
c----      xs_s,ys_s,zs_s : x,y,z-coordinates of the initial source position 
c----      xs_i,ys_i,zs_i : increment of source position
c----      n_rx : No. of Rx
c----      xr_s,yr_s,zr_s : x,y,z-coordinates of the initial receiver position 
c----      xr_i,yr_i,zr_i : increment of receiver position
c----      it_s : control parameter for Green's tensor integral and LUD
c----                  0 means first run and write the LUD results
c----                  1 means first run and do not write the LUD results
c----                  2 means using former results when frequency,
c----                    properties and geometries of the sheets are same 
c----      nprt(i) : control parameter for writing the output 
c----                  nprt(1) : write (0) the incident and scattering fields
c----                            in the sheets or not (1)
c----                  nprt(2) : write (0) electric field output or not (1)
c----                  nprt(3) : write (0) magnetic field output or not (1)
c----                  nprt(4) : secondary field and total field (0)
c----                            normalized secondary and total field by
c----                            primary field (1)
c-------------------------------------------------------------------------------
c---- ouput files :
c----    'INCIDENT.OUT' and 'SCATTER.OUT' when nprt(1) is zero
c----    'E-FIELD.DAT' when nprt(2) is zero
c----    'H-FIELD.DAT' when nprt(3) is zero
c----    'DECOMP.LUD' for storing the LUD results
c-------------------------------------------------------------------------------
c---- Required Subprogram : 'SUB_LUD.F'
c----    contains LU Decompsition, Romberg Integral and Bi-qubic spline
c-------------------------------------------------------------------------------
c---- July 18th, 1997
c---- Latest Revision  Sep. 16, 1997 for making the conductance of the sheets
c----       be the products of the thichnesses and the conductivity differences
c----       between host and sheets
c---- Nov. 7, 1997 for replacing the conductivity parameters including Cole-
c----       Cole model by those of resistivity
c---- March 11, 1998 for changing the coordinate and dip convention
c----       Same as HFSHEET which is the layered earth version
c---- Version 2.0, March 2002
c----       for using dynamic memory allocation to increase cell division
c-------------------------------------------------------------------------------
c
      dimension res_l(2), res_h(2), res_tau(2), res_alp(2)
      dimension eps_l(2), eps_h(2), eps_tau(2), eps_alp(2), thick(2)
      complex admit, imped, i_no, c_null, tau, cole, i_om,wave_s
      complex sig,eps,sig_n,eps_n,res,res_n
      complex f_prim(2,3), f_sec(2,3), ikk
      integer nprt(4)
      complex, allocatable :: sys(:,:), d_psi(:,:), d_phi(:,:), pot(:)
	real, allocatable :: x_loc(:), y_loc(:), z_loc(:)
	integer, allocatable :: indx(:)
      common /geo/ alpha(2), beta(2), aa(2), bb(2)
      common /locate/ x_top(2), y_top(2), z_top(2)
      common /prop/ admit, imped, wave_s, tau(2)
      common /parm/ omega, pi
      common /txrx/ ns_type, ns_dir, angle
      common /const/ i_no, c_null, conv 
      common /cell/ na(2), nb(2), max_n
      common /simp1/ ikk,i_d
c
      pi = 4.*atan(1.)
      eps0 = 8.85419e-12
      amyu0 = 4.e-7*pi
      conv = pi/180.
      i_no = cmplx (0.0, 1.0)
      c_null = cmplx (0.0, 0.0)
c
c---- initialization of parameters
c
      do i = 1, 2
         na(i) = 0
         nb(i) = 0
         aa(i) = 0.0
         bb(i) = 0.0
         alpha(i) = 0.0
         beta(i) = 0.0
         x_top(i) = 0.0
         y_top(i) = 0.0
         z_top(i) = 0.0
         tau(i) = c_null
         thick(i) = 0.0
      enddo
c
      open (2,file='wsheet20.par',status='old')
      read (2,*)
      read (2,*) no_sheet
      read (2,*)
      do i = 1, no_sheet
         read (2,*) na(i), nb(i)
      enddo
      read (2,*)
      do i = 1, no_sheet
         read (2,*) x_top(i), y_top(i), z_top(i)
      enddo
      read (2,*)
      do i = 1, no_sheet
         read (2,*) aa(i), bb(i), alp, bet, thick(i)
c
c---- Matching the convention of dip
c
         if (bet .ne. 0.) bet = 180. - bet
         alpha(i) = alp * conv
         beta(i) = bet * conv
      enddo
      read (2,*)
      do i = 1, no_sheet
         read (2,*) res_l(i), res_h(i), res_tau(i), res_alp(i)
      enddo
      read (2,*)
      do i = 1, no_sheet
         read (2,*) eps_l(i), eps_h(i), eps_tau(i), eps_alp(i)
      enddo
      read (2,*)
      read (2,*) resn_l, resn_h, resn_tau, resn_alp
      read (2,*)
      read (2,*) epsn_l, epsn_h, epsn_tau, epsn_alp
      read (2,*)
      read (2,*) freq
      read (2,*)
      read (2,*) ns_type, ns_dir, angle
      read (2,*)
      read (2,*) n_tx
      read (2,*)
      read (2,*) xs_s, ys_s, zs_s
      read (2,*) xs_i, ys_i, zs_i
      read (2,*)
      read (2,*) n_rx
      read (2,*)
      read (2,*) xr_s, yr_s, zr_s
      read (2,*) xr_i, yr_i, zr_i
      read (2,*)
      read (2,*) it_s
      read (2,*)
      read (2,*) (nprt(ip),ip=1,4)
      close (2)
c
      n_col = 2*(na(1)*nb(1) + na(2)*nb(2))
      ns=no_sheet
	max_na = na(1)
	max_nb = nb(1)
	do i = 1, ns
	   if (na(ns) .gt. na(1)) max_na = na(ns)
	   if (nb(ns) .gt. nb(1)) max_nb = nb(ns)
	enddo
	max_n = max (max_na, max_nb)
      angle = angle * conv
      omega = 2.*pi*freq
      i_om = i_no*omega
c
c---- Allocate and initialize arrays
c
      allocate (sys(n_col,n_col), d_psi(n_col,2), d_phi(n_col,2), 
     2          pot(n_col), x_loc(n_col), y_loc(n_col), z_loc(n_col),
	3          indx(n_col), stat=istat)
	if (istat .ne. 0) then
	    write(*,*)'>> Array allocation error !!'
	    stop
	endif
      do ij = 1, n_col
         indx(ij) = 0
         x_loc(ij) = 0.0
         y_loc(ij) = 0.0
         z_loc(ij) = 0.0
         pot(ij) = c_null
         do kl = 1, n_col
            sys(ij,kl) = c_null
         enddo
         do ic = 1, 2
            d_psi(ij,ic) = c_null
            d_phi(ij,ic) = c_null
         enddo
      enddo

c
c---- assign Cole-Cole relaxation to conductivity and dielectric constant
c
      res_n = cole(resn_l,resn_h,resn_tau,resn_alp)
      eps_n = cole(epsn_l,epsn_h,epsn_tau,epsn_alp) * eps0
	sig_n = 0.0
      if (cabs(res_n). gt. 0.) sig_n = 1./res_n
      do i = 1, ns
         res = cole(res_l(i),res_h(i),res_tau(i),res_alp(i))
         eps = cole(eps_l(i),eps_h(i),eps_tau(i),eps_alp(i)) * eps0
         sig = 1./res
         sig = sig - sig_n
         eps = eps - eps_n
         tau(i) = (sig + i_om*eps) * thick(i)
      enddo
c     
      admit = sig_n + i_om*eps_n
      imped = i_om*amyu0
      wave_s = -admit*imped
      ikk = i_no * csqrt(wave_s)
c
      if (nprt(1) .eq. 0) then
          open (3,file='incident.out')!,status='new',err=999)
          open (4,file='scatter.out')!,status='new',err=999)
      endif
      if (nprt(2) .eq. 0) then   
          open (8,file='e-field.dat')!,status='new',err=999)
          write(8,10)
      endif
      if (nprt(3) .eq. 0) then   
          open (9,file='h-field.dat')!,status='new',err=999)
          write(9,10)
      endif
c
c---- determine the x,y,z-coordinates of the center locations of the 
c---- cells in each sheets
c
      do is = 1, ns
         call ts_grid (is,n_col,x_loc,y_loc,z_loc)
      enddo
c
c---- Evaluation of Green's tensor integral over the sheets and perform
c---- Singular Value Decomposition
c
      if (it_s .ne. 2) then
          write(*,*) '>>> GREEN FUNCTION $ POTENTIAL EVALUATION'
          do is = 1, ns
             write(*,*) '     For Sheet No. S', is
             call green_0 (is,ns,n_col,sys)
          enddo
          write(*,*) '>>> LU DECOMPOSITION '
          call ludcmp (sys,n_col,n_col,indx,d)
          if (it_s .eq. 0) then
              open (1,file='decomp.lud')
              do ii = 1,n_col
                 write(1,*) ii, indx(ii)
              enddo
              do i = 1, n_col
                 do j = 1, n_col
                    write(1,*) i, j, sys(i,j)
                 enddo
              enddo
              close(1)
          endif
        else
          open (1,file='decomp.lud',status='old',err=888)
          do ii = 1, n_col
             read(1,*) id,indx(ii)
          enddo
          do i = 1, n_col
             do j = 1, n_col
                read(1,*)id, jd, sys(i,j)
                if ((id.ne.i) .or. (jd.ne.j)) then
                    write(*,*) '>>> INCORRECT LU-DECOMPOSITION DATA !!!'
                    stop
                endif
            enddo
          enddo
          close(1)
      endif
c
c---- Computing the incident electric fields in the sheets for each Tx location
c
      do it = 1, n_tx
         xs = xs_s + (it-1.)*xs_i
         ys = ys_s + (it-1.)*ys_i
         zs = zs_s + (it-1.)*zs_i
         do is = 1, ns
            call incident (is,nprt(1),n_col,xs,ys,zs,
     2                     x_loc,y_loc, z_loc, pot)
         enddo
c
c---- Solve the matrix equation and computing the scattering currents 
c---- in the sheets
c
         call lubksb(sys,n_col,n_col,indx,pot)
         do is = 1, ns
            call j_sheet (is, d_psi, d_phi, pot, n_col, nprt(1))
         enddo
c
c---- Computing the primary and secondary EM fields at the receiver location
c---- and write the results
c
         if (n_tx .eq. 1) then
            do ir = 1, n_rx
               write(*,*)
     2         '>>> Primary and Secondary Field computation for Rx', ir
               xr = xr_s + (ir-1.)*xr_i
               yr = yr_s + (ir-1.)*yr_i
               zr = zr_s + (ir-1.)*zr_i
               call primary(xs,ys,zs,xr,yr,zr,f_prim)
               call second_0 (ns, n_col,xr,yr,zr,d_psi,d_phi,f_sec)
               call out_data (n_tx,nprt,xs,ys,zs,xr,yr,zr,f_prim,f_sec)
            enddo
          else
            write(*,*)
     2      '>>> Primary and Secondary Field computation for Tx', it
            xr = xr_s + (it-1.)*xr_i
            yr = yr_s + (it-1.)*yr_i
            zr = zr_s + (it-1.)*zr_i
            call primary(xs,ys,zs,xr,yr,zr,f_prim)
            call second_0 (ns, n_col,xr,yr,zr,d_psi,d_phi,f_sec)
            call out_data (n_tx,nprt,xs,ys,zs,xr,yr,zr,f_prim,f_sec)
         endif
      enddo
      if (nprt(1) .eq. 0) then
          close(3)
          close(4)
      endif
      if (nprt(2) .eq. 0) then
          close(8)
      endif
      if (nprt(3) .eq. 0) then
          close(9)
      endif
c
      deallocate (sys, d_psi, d_phi, pot, x_loc, y_loc, z_loc, indx)
10    format (3x,'X',7x,'Y',7x,'Z', 7x,'Secondary - X',
     2        10x,'Total - X',12x,'Secondary - Y',10x,'Total - Y',
     3        12x,'Secondary - Z',10x,'Total - Z'/)
      stop
999   write (*,*) '>>> OUTPUT FILE ALREADY EXISTS !!!'
      stop
888   write (*,*) ' >>> DECOMPOSED DATA DO NOT EXISTS !!!'
      end
c
c-------------------------------------------------------------------------------
c
      subroutine out_data (n_tx,nprt,xs,ys,zs,xr,yr,zr,f_prim,f_sec)
      complex f_prim(2,3), f_sec(2,3), total(3), rel(3), c_norm
      integer nprt(4)
      common /txrx/ ns_type, ns_dir, angle
c
c---- For fixed source
c
      if (n_tx .eq. 1) then
           xm = xr
           ym = yr
           zm = zr
c
c---- For moving source and receiver, writes the coordinates of middle point
c
        else
           xm = 0.5 * (xs + xr)
           ym = 0.5 * (ys + yr)
           zm = 0.5 * (zs + zr)
      endif
c
      do ip = 2,3
         ip1 = ip - 1
         ifr = ip + 6
         if (nprt(ip) .eq. 0) then
             do i = 1,3
                rel(i) = f_sec(ip1,i)
                total(i) = f_prim(ip1,i) + rel(i) 
             enddo
             if (nprt(4) .ne. 0) then
                 do i = 1,3
                    c_norm = f_prim(ip1,i)
                    if (nprt(4) .eq. 2) c_norm = f_prim(ip1,ns_dir)
                    rel(i) = rel(i) / c_norm
                    total(i) = total(i) / c_norm
                 enddo
             endif 
             write(ifr,20) xm,ym,zm,(rel(i),total(i),i=1,3)
         endif
      enddo    
c
20    format (3f8.1,12(1x,1p1e11.4))
      return
      end
c
c-------------------------------------------------------------------------------
c---- subroutine calculates the scattering currents in the sheets using the
c---- numerical differentiation of potentials
c-------------------------------------------------------------------------------
c
      subroutine j_sheet (is,d_psi, d_phi, pot, n_col, icon)
      complex pot(n_col), scat_a, scat_b
      complex, allocatable :: psi(:,:), phi(:,:)
      complex admit, imped, i_no, c_null, tau, wave_s
      complex d_psi(n_col,2), d_phi(n_col,2)
      complex d_psi_a, d_psi_b, d_phi_a, d_phi_b
      dimension da(0:1,0:1), db(0:1,0:1)
      common /geo/ alpha(2), beta(2), aa(2), bb(2)
      common /prop/ admit, imped, wave_s, tau(2)
      common /parm/ omega, pi
      common /const/ i_no, c_null, conv 
      common /cell/ na(2), nb(2), max_n
c
      nab1 = na(1)*nb(1)
      nap = na(is)
      nbp = nb(is)
      del_a = aa(is) / nap
      del_b = bb(is) / nbp
c
      mn = max_n + 1
      allocate (psi(0:mn,0:mn),phi(0:mn,0:mn), stat=istat)
	if (istat .ne. 0) then
	    write(*,*) '>> Array allocation in J_SHEET error!!'
	    stop
	endif
      do k = 0, mn
         do l = 0, mn
            psi(k,l) = c_null
            phi(k,l) = c_null
         enddo
      enddo
c
c---- sign for numerical derivative 
c
      da(0,0) = 1.
      da(0,1) = 1.
      da(1,0) = -1.
      da(1,1) = -1.
      db(0,0) = 1.
      db(0,1) = -1.
      db(1,0) = 1.
      db(1,1) = -1.
c
c---- rearrange for divergence-free(psi) and curl-free(phi) potentials
c
      iskip = 2*(is-1)*nab1
      do k = 1, nap-1
         kk = iskip + (k-1)*(nbp-1)
         do l = 1, nbp-1
            kl = kk + l
            psi(k,l) = pot(kl)
         enddo
      enddo
      kskip = kl 
      do k = 0, nap-1
         kk = kskip + k*(nbp+1)
         do l = 0, nbp
            kl = kk + l + 1
            phi(k,l) = pot(kl)
         enddo
      enddo
      klast = kl
      do l = 0, nbp-2
         kl = klast+l+1
         phi(nap,l) = pot(kl)
      enddo 
c
c---- computation of scattering currents via numerical derivation of potentials
c---- 1 for second indice of d_psi and d_phi means derivative with respect to a
c---- 2 means with respect to b
c
      do i = 1, nap
         ii = (i-1)*nbp + iskip/2
         do j = 1, nbp
            ij = ii + j
            d_psi_a = c_null
            d_psi_b = c_null
            d_phi_a = c_null
            d_phi_b = c_null
            do m = 0, 1
               k = i - m
               do n = 0, 1
                  l = j - n
                  d_psi_a = d_psi_a + da(m,n)*psi(k,l)
                  d_psi_b = d_psi_b + db(m,n)*psi(k,l)
                  d_phi_a = d_phi_a + da(m,n)*phi(k,l)
                  d_phi_b = d_phi_b + db(m,n)*phi(k,l)
               enddo
            enddo
            d_psi(ij,1) = 0.5 * d_psi_a / del_a
            d_psi(ij,2) = 0.5 * d_psi_b / del_b
            d_phi(ij,1) = 0.5 * d_phi_a / del_a
            d_phi(ij,2) = 0.5 * d_phi_b / del_b
            scat_a =  d_psi(ij,2) + wave_s*d_phi(ij,1)
            scat_b = -d_psi(ij,1) + wave_s*d_phi(ij,2)
            if (icon .eq. 0) then
                write(4,12) i,j,scat_a,scat_b
            endif
         enddo
      enddo
c
      deallocate (psi,phi)
12    format(2i3,4(1x,g17.10))
      return
      end  
c
c-------------------------------------------------------------------------------
c---- Subroutine computing the secondary EM fields at the receiver located in
c---- the whole space.
c---- via numerical integration of electric and magnetic Green's function
c---- over the sheet and multiplying corresponding scattering currents.
c---- All the coordinates are transformed to those of the sheet.
c---- Note that the surface integral is to be done with respect to source
c---- coordinates, while the differentiation according to each EM field
c---- direction is with respect to receiver coordinates, which are sign
c---- reversal from each other.
c-------------------------------------------------------------------------------
c
      subroutine second_0 (ns,n_col,xr,yr,zr, d_psi,d_phi, f_sec)
      parameter (maxp = 9, maxpp = maxp*(maxp+15))
      complex f_sec(2,3), d_psi(n_col,2), d_phi(n_col,2)
      complex admit, imped, i_no, c_null, tau 
      complex dsaa, dfaa, dfab, dfbb, wave_s, pot, dhba
      complex dfca, dfcb, dfaa_i,dfaa_f,dfbb_i,dfbb_f
      complex e_aa, e_ba, e_ca, e_ab, e_bb, e_cb, e_a, e_b, e_c
      complex h_ba, h_ca, h_ab, h_cb, h_a, h_b, h_c
      complex deriv, ikk, ikr, dsca, dscb, d_pot
      complex scat_a, scat_b
      real r_quad(maxp,maxp), i_quad(maxp,maxp), wk(maxpp)
      real r_quad1(maxp,maxp), i_quad1(maxp,maxp)
      dimension d_i(maxp),d_j(maxp)
      common /geo/ alpha(2), beta(2), aa(2), bb(2)
      common /locate/ x_top(2), y_top(2), z_top(2)
      common /prop/ admit, imped, wave_s, tau(2)
      common /parm/ omega, pi
      common /const/ i_no, c_null, conv 
      common /cell/ na(2), nb(2), max_n
      common /simp1/ ikk, i_d
      common /simp2/ xx, yy, zz
      data np /9/
      external deriv
c
      dp1 = np - 1.0
      nab_1 = na(1)*nb(1)
      pi4 = 0.25 / pi
      do ic = 1, 2
         do id = 1, 3
            f_sec(ic,id) = c_null
         enddo
      enddo
      do is = 1, ns
         ijs = (is-1)*nab_1
c
c---- Rotation of coordinates of the receiver with respect to the center of
c---- is-th sheet to the strike and dip direction
c
         a_cos = ddcos(alpha(is))
         a_sin = sin(alpha(is))
         b_cos = ddcos(beta(is))
         b_sin = sin(beta(is))
         a_init = - 0.5*aa(is)
         b_init = - 0.5*bb(is)
         z_rise = b_init*b_sin
c   
         x_t = xr - x_top(is)
         y_t = yr - y_top(is)
         z_t = zr - z_top(is) + z_rise
         ap = x_t*a_cos + y_t*a_sin
         bp = -x_t*a_sin*b_cos + y_t*a_cos*b_cos + z_t*b_sin
         cp = x_t*a_sin*b_sin - y_t*a_cos*b_sin + z_t*b_cos
c
         del_a = aa(is) / na(is)
         del_b = bb(is) / nb(is)
         do i = 1, na(is)
            ii = (i-1)*nb(is) + ijs
            a_ii = a_init + (i-1.)*del_a
            a_ff = a_ii + del_a
            a_pi = ap - a_ii
            a_pf = ap - a_ff
            do j = 1, nb(is)
               ij = ii + j
               b_ii = b_init + (j-1.)*del_b
               b_ff = b_ii + del_b
               b_pi = bp - b_ii
               b_pf = bp - b_ff
               r_apbp = sqrt (a_pf**2. + b_pf**2. + cp**2.)
               r_ambp = sqrt (a_pi**2. + b_pf**2. + cp**2.)
               r_apbm = sqrt (a_pf**2. + b_pi**2. + cp**2.)
               r_ambm = sqrt (a_pi**2. + b_pi**2. + cp**2.)
               dfab = cexp(-ikk*r_apbp) / r_apbp 
     2              - cexp(-ikk*r_apbm) / r_apbm
     3              - cexp(-ikk*r_ambp) / r_ambp
     4              + cexp(-ikk*r_ambm) / r_ambm
               i_d = 1
               zz = cp
               yy = bp
               xx = a_pi
               call qromb(deriv,b_ii,b_ff,dfaa_i)
               xx = a_pf
               call qromb(deriv,b_ii,b_ff,dfaa_f)
               dfaa = a_pf*dfaa_f - a_pi*dfaa_i
               dfca = dfaa_f - dfaa_i 
               yy = ap
               xx = b_pi
               call qromb(deriv,a_ii,a_ff,dfbb_i)
               xx = b_pf
               call qromb(deriv,a_ii,a_ff,dfbb_f)
               dfbb = b_pf*dfbb_f - b_pi*dfbb_i
               dfcb = dfbb_f - dfbb_i 
               i_d = 0
               yy = bp
               xx = a_pi
               call qromb(deriv,b_ii,b_ff,dfaa_i)
               xx = a_pf
               call qromb(deriv,b_ii,b_ff,dfaa_f)
               dsca = -dfaa_f + dfaa_i
               yy = ap
               xx = b_pi
               call qromb(deriv,a_ii,a_ff,dfbb_i)
               xx = b_pf
               call qromb(deriv,a_ii,a_ff,dfbb_f)
               dscb = -dfbb_f + dfbb_i
               dsaa = c_null
               dhba = c_null  
               do im = 1, np
                  d_i(im) = a_ii + (im-1.)*del_a/dp1
                  do jm = 1, np
                     d_j(jm) = b_ii + (jm-1.)*del_b/dp1    
                     r_ij = sqrt((ap-d_i(im))**2. + 
     2                            (bp-d_j(jm))**2.+cp**2.)
                     ikr = ikk * r_ij
                     pot = cexp(-ikr)/r_ij
                     d_pot = (1.+ikr)*pot/r_ij**2. 
                     r_quad(im,jm) = real(pot)
                     i_quad(im,jm) = aimag(pot)
                     r_quad1(im,jm) = real(d_pot)
                     i_quad1(im,jm) = aimag(d_pot)
                  enddo
               enddo
	         a_ff = d_i(np) 
	         b_ff = d_j(np) 
               call DBCQDU (r_quad,maxp,d_i,np,d_j,np,a_ii,a_ff,
     2                      b_ii,b_ff,drr,WK,IER)
               if (ier. ne. 0) write(*,*) 'ier1 =', ier
               call DBCQDU (i_quad,maxp,d_i,np,d_j,np,a_ii,a_ff,
     2                      b_ii,b_ff,dii,WK,IER)
               if (ier. ne. 0) write(*,*) 'ier2 =', ier
               dsaa = cmplx (drr, dii)
               call DBCQDU (r_quad1,maxp,d_i,np,d_j,np,a_ii,a_ff,
     2                      b_ii,b_ff,drr,WK,IER)
               if (ier. ne. 0) write(*,*) 'ier3 =', ier
               call DBCQDU (i_quad1,maxp,d_i,np,d_j,np,a_ii,a_ff,
     2                      b_ii,b_ff,dii,WK,IER)
               if (ier. ne. 0) write(*,*) 'ier4 =', ier
               dhba = cmplx (drr, dii)

c
c---- Secondary field computation in sheet axis, a,b,c
c
               e_aa = dsaa*d_psi(ij,2)
     2               + (wave_s*dsaa + dfaa)*d_phi(ij,1)
               e_ba = dfab*d_phi(ij,1)
               e_ca = cp*dfca*d_phi(ij,1)
               e_ab = dfab*d_phi(ij,2)
               e_bb = -dsaa*d_psi(ij,1)
     2               + (wave_s*dsaa + dfbb)*d_phi(ij,2)
               e_cb = cp*dfcb*d_phi(ij,2)
c
               scat_a =  d_psi(ij,2) + wave_s*d_phi(ij,1)
               scat_b = -d_psi(ij,1) + wave_s*d_phi(ij,2)
               h_ba = -cp*dhba*scat_a
               h_ca = dscb*scat_a
               h_ab = -cp*dhba*scat_b
               h_cb = dsca*scat_b
c
               e_a = -imped*pi4*(e_aa+e_ab)
               e_b = -imped*pi4*(e_ba+e_bb)
               e_c = -imped*pi4*(e_ca+e_cb)
               h_a = - pi4*h_ab
               h_b = pi4*h_ba
               h_c = pi4*(-h_ca+h_cb)
c
c---- Rotation of coordinates for receiver axis, x,y,z
c
               f_sec(1,1) = f_sec(1,1) + ( e_a*a_cos - e_b*a_sin*b_cos 
     2                                 + e_c*a_sin*b_sin )
               f_sec(1,2) = f_sec(1,2) + ( e_a*a_sin + e_b*a_cos*b_cos
     2                                 - e_c*a_cos*b_sin )
               f_sec(1,3) = f_sec(1,3) + ( e_b*b_sin + e_c*b_cos )
               f_sec(2,1) = f_sec(2,1) + ( h_a*a_cos - h_b*a_sin*b_cos
     2                                 + h_c*a_sin*b_sin )
               f_sec(2,2) = f_sec(2,2) + ( h_a*a_sin + h_b*a_cos*b_cos
     2                                 - h_c*a_cos*b_sin )
               f_sec(2,3) = f_sec(2,3) + ( h_b*b_sin + h_c*b_cos )
c
            enddo
         enddo
      enddo
c
      return
      end    
c
c-------------------------------------------------------------------------------
c---- Evaluate the integral of primary Green's functions which is for the whole
c---- space and construct the system matrix which is composed of divergence-
c---- free and curl-free potentials at each nodal points.
c---- When only one sheet exists or strikes and dips of two sheets are same, 
c---- evaluation fo Green's tensor integral is simple
c---- because of symmetric properties of the whole space Green's function and
c---- the shape of the cell as following. 
c----                  S0_aa = S0_bb, S0_ab = S0_ba = 0
c----                  F0_aa, F0_bb, F0_ab = F0_ba
c---- However, in general, all the eight components of Green's tensor should
c---- stored separately. 
c---- np : No. of subdivision of cell for bi-spline quadrature integral
c-------------------------------------------------------------------------------
c  
      subroutine green_0 (is,ns,n_col,sys)
      parameter (maxp=9, maxpp=maxp*(maxp+15) )
      complex sys(n_col,n_col)
      complex, allocatable :: saa(:,:),sab(:,:),sba(:,:),sbb(:,:)
      complex, allocatable :: faa(:,:),fab(:,:),fba(:,:),fbb(:,:)
      complex admit, imped, i_no, ikk, c_null
      complex dsaa, dfaa, dfab, dfbb, wave_s, tau, ctau
      complex dsxa, dfxa, dfya, dfza, dfyb, dfzb
      complex paa, pab, pba, pbb, qaa, qab, qba, qbb
      complex ssaa, sfaa, sfbb,dfaa_i,dfaa_f,dfbb_i,dfbb_f
      complex deriv, pot
      real r_quad(maxp,maxp), i_quad(maxp,maxp), wk(maxpp)
      dimension d_i(maxp),d_j(maxp)
      dimension da(0:1,0:1), db(0:1,0:1)
      common /geo/ alpha(2), beta(2), aa(2), bb(2)
      common /locate/ x_top(2), y_top(2), z_top(2)
      common /prop/ admit, imped, wave_s, tau(2)
      common /parm/ omega, pi
      common /const/ i_no, c_null, conv 
      common /cell/ na(2), nb(2), max_n
      common /simp1/ ikk, i_d
      common /simp2/ xx, yy, zz
      data np /9/
      external deriv
c
      mn = max_n + 1
c
      allocate (saa(0:mn,0:mn),sab(0:mn,0:mn),sba(0:mn,0:mn),
     2          sbb(0:mn,0:mn),faa(0:mn,0:mn),fab(0:mn,0:mn),
	3          fba(0:mn,0:mn),fbb(0:mn,0:mn), stat=istat)
	if (istat .ne. 0) then
	    write(*,*) '>> Array allocation in GREEN_0 error!!'
	    stop
	endif
c
      do k = 0, mn
         do l = 0, mn
            saa(k,l) = c_null
            faa(k,l) = c_null
            fab(k,l) = c_null
            fbb(k,l) = c_null
         enddo
      enddo
c
      do ip = 1, maxp
         do jp = 1, maxp
            r_quad(ip,jp) = 0.0
            i_quad(ip,jp) = 0.0
         enddo
      enddo
      dp1 = np - 1.0
      pi4 = 0.25 / pi
      ins = 3 - is
      nap = na(is)
      nbp = nb(is)
      nas = na(ins)
      nbs = nb(ins)
      del_a = aa(is) / nap 
      del_b = bb(is) / nbp 
      nabp = nap * nbp
      nabp_1 = (nap-1) * (nbp-1)
      nabs = nas * nbs
      nabs_1 = (nas-1) * (nbs-1)
      a_init = - 0.5*aa(is)
      b_init = - 0.5*bb(is)
      as_init = - 0.5*aa(ins)
      bs_init = - 0.5*bb(ins) 
      i_d = 1
c
c---- sign for numerical derivative 
c
      da(0,0) = 1.
      da(0,1) = 1.
      da(1,0) = -1.
      da(1,1) = -1.
      db(0,0) = 1.
      db(0,1) = -1.
      db(1,0) = 1.
      db(1,1) = -1.
c
c---- For singular integral when r_ij = r_kl
c
      call singular (del_a, del_b, ssaa, sfaa, sfbb)
c
c-------------------------------------------------------------------------------
c---- Evaluation of Green's tensor integral for the same sheet
c---- surface integral and derivatives are with respect to k,l coordinates
c-------------------------------------------------------------------------------
c
      ctau = tau(is)
      zz = 0.0
      do i = 1, nap
         ii = (i-1)*nbp + 2*(is-1)*nabs
         a_i = a_init + (i-0.5)*del_a
         do j = 1, nbp
            ij = ii + j
            b_j = b_init + (j-0.5)*del_b
            do k = 1, nap
               a_ii = a_init + (k-1.)*del_a
               a_ff = a_ii + del_a
               a_di = a_i - a_ii
               a_df = a_i - a_ff
               do l = 1, nbp
                  b_ii = b_init + (l-1.)*del_b
                  b_ff = b_ii + del_b
                  b_di = b_j - b_ii
                  b_df = b_j - b_ff
c
c---- For singular cell
c
                  if ((i.eq.k) .and. (j.eq.l)) then
                        saa(k,l) = ssaa 
                        faa(k,l) = sfaa
                        fbb(k,l) = sfbb
                        fab(k,l) = c_null
c
c---- Green's function with r = sqrt (r_ij - r_kl)
c
                      else
                        r_apbp = sqrt (a_df**2. + b_df**2.)
                        r_apbm = sqrt (a_df**2. + b_di**2.)
                        r_ambp = sqrt (a_di**2. + b_df**2.)
                        r_ambm = sqrt (a_di**2. + b_di**2.)
                        dfab = cexp(-ikk*r_apbp) / r_apbp 
     2                       - cexp(-ikk*r_apbm) / r_apbm
     3                       - cexp(-ikk*r_ambp) / r_ambp
     4                       + cexp(-ikk*r_ambm) / r_ambm
                        yy = a_i
                        xx = b_di
                        call qromb(deriv,a_ii,a_ff,dfbb_i)
                        xx = b_df
                        call qromb(deriv,a_ii,a_ff,dfbb_f)
                        dfbb = b_df*dfbb_f - b_di*dfbb_i 
                        yy = b_j
                        xx = a_di
                        call qromb(deriv,b_ii,b_ff,dfaa_i)
                        xx = a_df
                        call qromb(deriv,b_ii,b_ff,dfaa_f)
                        dfaa = a_df*dfaa_f - a_di*dfaa_i 
                        do im = 1, np
                           d_i(im) = a_ii + (im-1.)*del_a/dp1
                           do jm = 1, np
                              d_j(jm) = b_ii + (jm-1.)*del_b/dp1    
                              r_ij = sqrt((a_i-d_i(im))**2.  
     2                              + (b_j-d_j(jm))**2.)
                              pot = cexp(-ikk*r_ij)/r_ij
                              r_quad(im,jm) = real(pot)
                              i_quad(im,jm) = aimag(pot)
                           enddo
                        enddo
	                  a_ff = d_i(np)
	                  b_ff = d_j(np)
                        call DBCQDU (r_quad,maxp,d_i,np,d_j,np,
     2                               a_ii,a_ff,b_ii,b_ff,drr,WK,IER)
                        if (ier. ne. 0) write(*,*) 'ier_p1 =', ier
                        call DBCQDU (i_quad,maxp,d_i,np,d_j,np,
     2                               a_ii,a_ff,b_ii,b_ff,dii,WK,IER)
                        if (ier. ne. 0) write(*,*) 'ier_p2 =', ier
                        dsaa = cmplx (drr, dii)
                        saa(k,l) = pi4 * dsaa
                        faa(k,l) = pi4 * dfaa
                        fbb(k,l) = pi4 * dfbb
                        fab(k,l) = pi4 * dfab
                  endif
               enddo
            enddo
c
c---- rearrange the Green's tensor integral to approximate the numerical
c---- derivatives of diverence-free (P) and curl-free (Q) potentials
c---- and construct system matrix
c
            do k = 0, nap
               kp = (k-1) * (nbp-1)
               kq = k * (nbp+1)
               do l = 0, nbp
                  kl = k*l
                  kl_p = kp + l + 2*(is-1)*nabs
                  kl_q = kq + l + 1 + nabp_1 + 2*(is-1)*nabs
c
c---- The last terms on the indice indicate the shift of element locations
c---- for the Green's tensor due to same sheet
c
c                     2*Na1*Nb1     +     2*Na2*Nb2      = N_COL
c          (Na1-1)(Nb1-1) | (Na1+1)(Nb1+1)-2 
c                |------------------|------------------|  
c                |        |         |                  |
c                |   PA   |  QA     |                  |
c                |        |         |                  |
c                |--------|---------|                  |
c                |      S1|<- S1    |      S1 <- S2    | 2*Na1*Nb1
c                |   PB   |  QB     |                  |
c           ij   |        |         |                  |
c                |------------------|------------------|  
c                |                  |                  |
c                |                  |                  |
c                |      S2 <- S1    |      S2 <- S2    |
c                |                  |                  | 2*Na2*Nb2
c                |                  |                  |
c                |                  |                  |
c                |------------------|------------------|  
c                    kl_p      kl_q
c
                  paa = c_null
                  pbb = c_null
                  qaa = c_null
                  qab = c_null
                  qba = c_null
                  qbb = c_null
                  do m = 0, 1
                     mk = m+k
                     do n = 0, 1
                        nl = n+l
                        paa = paa + db(m,n)*imped*saa(mk,nl)
                        pbb = pbb + da(m,n)*imped*saa(mk,nl)
                        qaa = qaa + da(m,n)*imped
     2                            *(wave_s*saa(mk,nl)+faa(mk,nl))
                        qab = qab + db(m,n)*imped*fab(mk,nl)
                        qba = qba + da(m,n)*imped*fab(mk,nl)
                        qbb = qbb + db(m,n)*imped
     2                            *(wave_s*saa(mk,nl)+fbb(mk,nl))
                        if ((mk.eq.i) .and. (nl.eq.j)) then
                              paa = paa + db(m,n) / ctau
                              pbb = pbb + da(m,n) / ctau
                              qaa = qaa + da(m,n)*wave_s/ctau
                              qbb = qbb + db(m,n)*wave_s/ctau
                        endif
                     enddo
                  enddo
             
                  if (kl.gt.0) then
                     if ((k.lt.nap) .and. (l.lt.nbp)) then
c---- PA_ij
                         sys(ij,kl_p) = 0.5*paa/del_b
c---- PB_ij
                         sys(ij+nabp,kl_p) = -0.5*pbb/del_a
                     endif
                  endif
                  if (k.lt.nap) then
c---- QA_ij
                      sys(ij,kl_q) = 0.5*(qaa/del_a + qab/del_b)
c---- QB_ij
                      sys(ij+nabp,kl_q) = 0.5*
     2                                   (qba/del_a + qbb/del_b)
                    elseif (l.lt.(nbp-1)) then
c---- QA_ij
                      sys(ij,kl_q) = 0.5*(qaa/del_a + qab/del_b)
c---- QB_ij
                      sys(ij+nabp,kl_q) = 0.5*
     2                                   (qba/del_a + qbb/del_b)
                  endif
               enddo
            enddo
         enddo
      enddo
      if (ns .eq. 1) return
c
c-------------------------------------------------------------------------------
c---- Green's tensor intgral over the other sheet
c---- Rotation of coordinates so as for the axis of the other sheet to be
c---- reference coordinates
c-------------------------------------------------------------------------------
c
c---- make the origin and the coordinates be those of the other sheet
c
      del_as = aa(ins) / nas
      del_bs = bb(ins) / nbs
      alp = alpha(is) - alpha(ins)
      bet = beta(is) - beta(ins)
      a_cos = ddcos(alp)
      a_sin = sin(alp)
      b_cos = ddcos(bet)
      b_sin = sin(bet)
      a_cos1 = ddcos(alpha(ins))
      a_sin1 = sin(alpha(ins))
      b_cos1 = ddcos(beta(ins))
      b_sin1 = sin(beta(ins))
      x_p1 = x_top(is) - x_top(ins) 
      y_p1 = y_top(is) - y_top(ins)
      z_p1 = z_top(is) - z_top(ins)
      x_p = x_p1*a_cos1 + y_p1*a_sin1
      y_p = -x_p1*a_sin1*b_cos1 + y_p1*a_cos1*b_cos1 + z_p1*b_sin1
      z_p = x_p1*a_sin1*b_sin1 - y_p1*a_cos1*b_sin1 + z_p1*b_cos1
      do k = 0, mn
         do l = 0, mn
            saa(k,l) = c_null
            sab(k,l) = c_null
            sba(k,l) = c_null
            sbb(k,l) = c_null
            faa(k,l) = c_null
            fab(k,l) = c_null
            fba(k,l) = c_null
            fbb(k,l) = c_null
         enddo
      enddo
c
      do i = 1, nap
         ii = (i-1)*nbp + 2*(is-1)*nabs
         a_temp = a_init + (i-0.5)*del_a
         do j = 1, nbp
            ij = ii + j
            b_temp = b_init + (j-0.5)*del_b
            xs = x_p + a_temp*a_cos - b_temp*a_sin*b_cos
            ys = y_p + a_temp*a_sin + b_temp*a_cos*b_cos
            zs = z_p + b_temp*b_sin 
            do ks = 1, nas
               a_ii = as_init + (ks-1.)*del_as
               a_ff = a_ii + del_as
               as_i = xs - a_ii
               as_f = xs - a_ff
               do ls = 1, nbs
                  b_ii = bs_init + (ls-1.)*del_bs
                  b_ff = b_ii + del_bs
                  bs_i = ys - b_ii
                  bs_f = ys - b_ff
                  r_apbp = sqrt (as_f**2. + bs_f**2. + zs**2.)
                  r_ambp = sqrt (as_i**2. + bs_f**2. + zs**2.)
                  r_apbm = sqrt (as_f**2. + bs_i**2. + zs**2.)
                  r_ambm = sqrt (as_i**2. + bs_i**2. + zs**2.)
                  dfya = cexp(-ikk*r_apbp) / r_apbp 
     2                 - cexp(-ikk*r_apbm) / r_apbm
     3                 - cexp(-ikk*r_ambp) / r_ambp
     4                 + cexp(-ikk*r_ambm) / r_ambm
                  zz = zs
                  yy = ys
                  xx = as_i
                  call qromb(deriv,b_ii,b_ff,dfaa_i)
                  xx = as_f
                  call qromb(deriv,b_ii,b_ff,dfaa_f)
                  dfxa = as_f*dfaa_f - as_i*dfaa_i 
                  dfza = zs*(dfaa_f - dfaa_i)
                  yy = xs
                  xx = bs_i
                  call qromb(deriv,a_ii,a_ff,dfbb_i)
                  xx = bs_f
                  call qromb(deriv,a_ii,a_ff,dfbb_f)
                  dfyb = bs_f*dfbb_f - bs_i*dfbb_i
                  dfzb = zs*(dfbb_f - dfbb_i)
                  dsxa = c_null
                  do im = 1, np
                     d_i(im) = a_ii + (im-1.)*del_as/dp1
                     do jm = 1, np
                        d_j(jm) = b_ii + (jm-1.)*del_bs/dp1    
                        r_ij = sqrt((xs-d_i(im))**2.  
     2                             + (ys-d_j(jm))**2.+zs**2.)
                        pot = cexp(-ikk*r_ij)/r_ij
                        r_quad(im,jm) = real(pot)
                        i_quad(im,jm) = aimag(pot)
                     enddo
                  enddo
	            a_ff = d_i(np)
	            b_ff = d_j(np)
                  call DBCQDU (r_quad,maxp,d_i,np,d_j,np,
     2                         a_ii,a_ff,b_ii,b_ff,drr,WK,IER)
                  if (ier. ne. 0) write(*,*) 'ier_rs =', ier
                  call DBCQDU (i_quad,maxp,d_i,np,d_j,np,
     2                         a_ii,a_ff,b_ii,b_ff,dii,WK,IER)
                  if (ier. ne. 0) write(*,*) 'ier_is=', ier
                  dsxa = cmplx (drr,dii)
c
                  saa(ks,ls) = pi4 * a_cos * dsxa 
                  sab(ks,ls) = pi4 * a_sin * dsxa 
                  sba(ks,ls) = -pi4 * a_sin * b_cos * dsxa 
                  sbb(ks,ls) = pi4 * a_cos * b_cos * dsxa
                  faa(ks,ls) = pi4 * (a_cos*dfxa + a_sin*dfya)
                  fab(ks,ls) = pi4 * (a_cos*dfya + a_sin*dfyb) 
                  fba(ks,ls) = pi4 * (-a_sin*b_cos*dfxa
     2                               + a_cos*b_cos*dfya 
     3                               + b_sin*dfza)
                  fbb(ks,ls) = pi4 * (-a_sin*b_cos*dfya
     2                               + a_cos*b_cos*dfyb
     3                               + b_sin*dfzb)
               enddo
            enddo
c
c---- rearrange the Green's tensor integral to approximate the numerical
c---- derivatives of diverence-free (P) and curl-free (Q) potentials
c---- and compose system matrix.
c---- There is no singular cell in this case
c
            do k = 0, nas
               kp = (k-1) * (nbs-1)
               kq = k * (nbs+1)
               do l = 0, nbs
                  kl = k*l
                  kl_p = kp + l + 2*(ins-1)*nabp
                  kl_q = kq + l + 1 + nabs_1 + 2*(ins-1)*nabp
                  paa = c_null
                  pab = c_null
                  pba = c_null
                  pbb = c_null
                  qaa = c_null
                  qab = c_null
                  qba = c_null
                  qbb = c_null
                  do m = 0, 1
                     do n = 0, 1
                        mk = m+k
                        nl = n+l
                        paa = paa + db(m,n)*imped*saa(mk,nl)
                        pab = pab + da(m,n)*imped*sab(mk,nl)
                        pba = pba + db(m,n)*imped*sba(mk,nl)
                        pbb = pbb + da(m,n)*imped*sbb(mk,nl)
                        qaa = qaa + da(m,n)*imped
     2                            *(wave_s*saa(mk,nl)+faa(mk,nl))
                        qab = qab + db(m,n)*imped
     2                            *(wave_s*sab(mk,nl)+fab(mk,nl))
                        qba = qba + da(m,n)*imped
     2                            *(wave_s*sba(mk,nl)+fba(mk,nl))
                        qbb = qbb + db(m,n)*imped
     2                            *(wave_s*sbb(mk,nl)+fbb(mk,nl))
                     enddo
                  enddo
                  if (kl.gt.0) then
                        if ((k.lt.nas) .and. (l.lt.nbs)) then
c---- PA_ij
                            sys(ij,kl_p) = 0.5*
     2                                 (paa/del_bs - pab/del_as)
c---- PB_ij
                            sys(ij+nabp,kl_p) = 0.5*
     2                                 (pba/del_bs - pbb/del_as)
                        endif
                  endif
                  if (k.lt.nas) then
c---- QA_ij
                         sys(ij,kl_q) = 0.5*
     2                                 (qaa/del_as + qab/del_bs)
c---- QB_ij
                         sys(ij+nabp,kl_q) = 0.5*
     2                                 (qba/del_as + qbb/del_bs)
                     elseif (l.lt.(nbs-1)) then
c---- QA_ij
                         sys(ij,kl_q) = 0.5*
     2                                 (qaa/del_as + qab/del_bs)
c---- QB_ij
                         sys(ij+nabp,kl_q) = 0.5*
     2                                 (qba/del_as + qbb/del_bs)
                  endif
               enddo
            enddo
         enddo
      enddo
c
      deallocate (saa,sab,sba,sbb,faa,fab,fba,fbb)

      return
      end
c
c-------------------------------------------------------------------------------
c    
      complex function deriv(ab)
      complex ikk, ikr, pot
      common /simp1/ ikk, i_d
      common /simp2/ xx, yy, zz
c
      rr = sqrt (xx**2. + (yy-ab)**2. + zz**2.)
      ikr = ikk * rr
      pot = cexp(-ikr) / rr
c
      if (i_d. eq. 0) then
            deriv = pot
         else
            deriv = pot * (1.+ikr)/rr**2.
      endif
c
      return
      end        
c
c-------------------------------------------------------------------------------
c---- Evaluate the singular Green's function in which the observing and source
c---- points are indentical
c---- Analytic forms are availble for square when del_a and del_b are equal.
c---- Numerical integration for asymmetric parts
c-------------------------------------------------------------------------------
c
      subroutine singular (del_a, del_b, saa, faa, fbb)
      parameter (maxp=20, maxpp=maxp*(maxp+15))
      complex saa, faa, fbb, dfaa_i, dfaa_f, dfbb, ctemp
      complex i_no, ikrho, e_ikrho, c_null, ikk
      complex deriv, ikr, pot
      real r_quad(maxp,maxp), i_quad(maxp,maxp), wk(maxpp)
      dimension d_i(maxp),d_j(maxp)
      common /parm/ omega, pi
      common /const/ i_no, c_null, conv 
      common /simp1/ ikk, i_d
      common /simp2/ xx, yy, zz
      external deriv
      data np /15/
c
      del = min (del_a, del_b)
      del2 = del / 2.
      dif = abs(del_a - del_b) / 2.
      f_l = del2 + dif
      saa = c_null
      zz = 0.0
c
c---- Analytic computation of singular integral of Green's function
c 
      rho = del / sqrt(pi)
      ikrho = ikk*rho
      e_ikrho = cexp(-ikrho)
      faa = -0.25 * (1.+ikrho) * e_ikrho / rho
      fbb = faa
c
c---- Series expansion at small wave propagation constant for stability
c
      if (cabs(ikk) .gt. 1.e-3) then
             saa = 0.5 * (1.0 - e_ikrho) / ikk
          elseif (cabs(ikk) .eq. 0.0) then
             saa = 0.5 * rho
          else
             saa = cmplx(1.0, 0.0)
             do i = 2, 10
                saa = saa + (-ikrho)**(i-1.) / facto(i)
             enddo
             saa = 0.5 * rho * saa
      endif
c
c---- Numerical integration over asymmetric part of singular cell
c 
      if (dif .ne. 0.0) then
            i_d = 1
            yy = 0.0
            xx = del2
            call qromb(deriv,0.0,del2,dfaa_f)
            xx = f_l
            call qromb(deriv,0.0,del2,dfaa_i)
            faa = faa + (del2*dfaa_f - f_l*dfaa_i) / pi
            xx = del2
            call qromb(deriv,del2,f_l,dfbb)
            fbb = fbb - del2*dfbb/pi
            do im = 1, np
               d_i(im) = (im-1)*del2/(np-1)
               do jm = 1, np
                  d_j(jm) = del2 + (jm-1)*dif/(np-1)    
                  r_ij = sqrt(d_i(im)**2. + d_j(jm)**2.)
                  ikr = ikk * r_ij
                  pot = cexp(-ikr)/r_ij
                  r_quad(im,jm) = real(pot)
                  i_quad(im,jm) = aimag(pot)
               enddo
            enddo
	      d_i(np) = del2
	      d_j(np) = f_l
            call DBCQDU (r_quad,maxp,d_i,np,d_j,np,
     2                   0.0,del2,del2,f_l,drr,WK,IER)
            if (ier. ne. 0) write(*,*) 'ier_sr =', ier
            call DBCQDU (i_quad,maxp,d_i,np,d_j,np,
     2                   0.0,del2,del2,f_l,dii,WK,IER)
            if (ier. ne. 0) write(*,*) 'ier_si =', ier
            saa = saa + cmplx (drr, dii)/pi
            if (del_b. gt. del_a) then
                    ctemp = fbb
                    fbb = faa
                    faa = ctemp
            endif
      endif    
c
      return
      end 
c
c-------------------------------------------------------------------------------
c
      real function facto (nmax)
c
      fact = 1.0
c 
      do i = 1, nmax
         fact = fact * float(i)
      enddo
      facto = fact
c
      return
      end
c
c-------------------------------------------------------------------------------
c
      subroutine primary (xs,ys,zs,xr,yr,zr,f_prim)
      complex admit, imped, wavep, i_no, ikr, e_ikr, c_null, tau
      complex Ex, Ey, Ez, c_ratio, deriv_0, deriv_1, deriv_2
      complex Hx, Hy, Hz, f_prim(2,3), wave_s
      common /prop/ admit, imped, wave_s, tau(2)
      common /const/ i_no, c_null, conv 
      common /parm/ omega, pi
      common /txrx/ ns_type, ns_dir, angle
c
      pi4 = 0.25 / pi
      wavep = csqrt(wave_s)
c
      x_dist = xr - xs
      y_dist = yr - ys
      z_dist = zr - zs
      rr = sqrt(x_dist**2. + y_dist**2. + z_dist**2.)
      ikr = i_no*wavep*rr
      e_ikr = cexp(-ikr)
      dx = x_dist / rr
      dy = y_dist / rr
      dz = z_dist / rr
      deriv_0 = e_ikr / rr
      deriv_1 = e_ikr * (ikr + 1.) / rr**2.
      deriv_2 = e_ikr * (ikr*ikr + 3.*ikr + 3) / rr**3.
c
c---- electric dipole source
c
      if (ns_type .eq. 1) then
           c_ratio = imped / admit
           if (ns_dir .eq. 1) then
                Ex = - admit*deriv_0  
     2               + (-deriv_1/rr + deriv_2*dx**2.) / imped
                Ey = deriv_2*dx*dy / imped
                Ez = deriv_2*dz*dx / imped
                Hx = c_null
                Hy = -deriv_1*dz
                Hz = deriv_1*dy
              elseif (ns_dir .eq. 2) then
                Ex = deriv_2*dx*dy / imped  
                Ey = -admit*deriv_0 
     2               + (-deriv_1/rr + deriv_2*dy**2.) / imped
                Ez = deriv_2*dy*dz / imped
                Hx = deriv_1*dz 
                Hy = c_null
                Hz = -deriv_1*dx
              else
                Ex = deriv_2*dz*dx / imped
                Ey = deriv_2*dx*dy / imped
                Ez = -admit*deriv_0 
     2               + (-deriv_1/rr + deriv_2*dz**2.) / imped
                Hx = -deriv_1*dy
                Hy = deriv_1*dx
                Hz = c_null
           endif
           f_prim(1,1) = Ex * c_ratio * pi4 
           f_prim(1,2) = Ey * c_ratio * pi4 
           f_prim(1,3) = Ez * c_ratio * pi4
           f_prim(2,1) = Hx * pi4
           f_prim(2,2) = Hy * pi4 
           f_prim(2,3) = Hz * pi4 
c
c---- magnetic dipole source
c
         else
           if (ns_dir .eq. 1) then
                Ex = c_null
                Ey = deriv_1*dz
                Ez = -deriv_1*dy
                Hx = wave_s*deriv_0 -deriv_1/rr + deriv_2*dx**2.
                Hy = deriv_2*dx*dy
                Hz = deriv_2*dz*dx              
              elseif (ns_dir .eq. 2) then
                Ex = -deriv_1*dz
                Ey = c_null
                Ez = deriv_1*dx
                Hx = deriv_2*dx*dy
                Hy = wave_s*deriv_0 -deriv_1/rr + deriv_2*dy**2.
                Hz = deriv_2*dy*dz              
              else
                Ex = deriv_1*dy
                Ey = -deriv_1*dx
                Ez = c_null
                Hx = deriv_2*dz*dx
                Hy = deriv_2*dy*dz              
                Hz = wave_s*deriv_0 -deriv_1/rr + deriv_2*dz**2.
           endif
           f_prim(1,1) = pi4 * imped * Ex 
           f_prim(1,2) = pi4 * imped * Ey 
           f_prim(1,3) = pi4 * imped * Ez 
           f_prim(2,1) = pi4 * Hx
           f_prim(2,2) = pi4 * Hy
           f_prim(2,3) = pi4 * Hz
      endif
c
      return
      end
c
c-------------------------------------------------------------------------------
c
      subroutine incident (is,nprt,n_col,xs,ys,zs,
     2                     x_loc,y_loc,z_loc,pot)
      dimension x_loc(n_col), y_loc(n_col), z_loc(n_col)
      complex admit, imped, wavep, i_no, ikr, e_ikr, c_null
      complex prim_a, prim_b, tau, wave_s, pot(n_col)
      complex Ex, Ey, Ez, c_ratio, deriv_0, deriv_1, deriv_2
      common /geo/ alpha(2), beta(2), aa(2), bb(2)
      common /locate/ x_top(2), y_top(2), z_top(2)
      common /prop/ admit, imped, wave_s, tau(2)
      common /parm/ omega, pi
      common /txrx/ ns_type, ns_dir, angle
      common /cell/ na(2), nb(2), max_n
      common /const/ i_no, c_null, conv 
c
      pi4 = 0.25 / pi
      dt_cos = ddcos(angle)
      dt_sin = sin(angle)
      wavep = csqrt(wave_s)
      naa = na(is)
      nbb = nb(is)
      nab = naa*nbb
      iskip = (is-1)*2*na(1)*nb(1)
      da_cos = ddcos(alpha(is))
      da_sin = sin(alpha(is))
      db_cos = ddcos(beta(is))
      db_sin = sin(beta(is))
c
      do i = 1, naa
         ii = (i-1) * nbb + iskip/2
         ijs = (i-1) * nbb + iskip
         do j = 1, nbb
            ij = ii + j
            ij_a = ijs + j
            ij_b = ij_a + nab
            x_dist = x_loc(ij) - xs
            y_dist = y_loc(ij) - ys
            z_dist = z_loc(ij) - zs
            rr = sqrt(x_dist**2. + y_dist**2. + z_dist**2.)
            ikr = i_no*wavep*rr
            e_ikr = cexp(-ikr)
            dx = x_dist / rr
            dy = y_dist / rr
            dz = z_dist / rr
            deriv_0 = e_ikr / rr
            deriv_1 = e_ikr * (ikr + 1.) / rr**2.
            deriv_2 = e_ikr * (ikr*ikr + 3.*ikr + 3) / rr**3.
c
c---- plane wave source
c
            if (ns_type .eq. 0) then
c---- TM
                if (ns_dir .eq. 1) then
                    Ex = dt_cos * e_ikr
                    Ey = c_null
                    Ez = -dt_sin * e_ikr
                 else
c---- TE
                    Ex = c_null
                    Ey = e_ikr
                    Ez = c_null
                endif
c
c---- electric dipole source
c
              elseif (ns_type .eq. 1) then
                c_ratio = imped / admit
                if (ns_dir .eq. 1) then
                    Ex = - admit*deriv_0 
     2                   + (-deriv_1/rr + deriv_2*dx**2.) / imped
                    Ey = deriv_2*dx*dy / imped
                    Ez = deriv_2*dz*dx
                 elseif (ns_dir .eq. 2) then
                    Ex = deriv_2*dx*dy / imped  
                    Ey = - admit*deriv_0 
     2                   + (-deriv_1/rr + deriv_2*dy**2.) / imped
                    Ez = deriv_2*dy*dz / imped
                 else
                    Ex = deriv_2*dz*dx / imped
                    Ey = deriv_2*dx*dy / imped
                    Ez = - admit*deriv_0 
     2                   + (-deriv_1/rr + deriv_2*dz**2.) / imped
                endif
                Ex = pi4 * Ex * c_ratio  
                Ey = pi4 * Ey * c_ratio  
                Ez = pi4 * Ez * c_ratio  
c
c---- magnetic dipole source
c
              else
                if (ns_dir .eq. 1) then
                    Ex = c_null
                    Ey = deriv_1*dz
                    Ez = -deriv_1*dy              
                 elseif (ns_dir .eq. 2) then
                    Ex = -deriv_1*dz
                    Ey = c_null
                    Ez = deriv_1*dx
                 else
                    Ex = deriv_1*dy
                    Ey = -deriv_1*dx
                    Ez = c_null
                endif
                Ex = pi4 * imped * Ex
                Ey = pi4 * imped * Ey
                Ez = pi4 * imped * Ez
            endif
c
            prim_a = Ex*da_cos + Ey*da_sin
            prim_b = -Ex*da_sin*db_cos + Ey*da_cos*db_cos
     2                   + Ez*db_sin
            pot(ij_a) = prim_a
            pot(ij_b) = prim_b
         enddo
      enddo
c
c---- write the incident electric fields depending on option nprt
c
      if (nprt .eq. 0) then
          do i = 1, naa
             ii = (i-1) * nbb + iskip
             do j = 1, nbb
                ij_a = ii + j
                ij_b = ij_a + nab
                write(3,11) i,j,pot(ij_a),pot(ij_b)
             enddo
          enddo
      endif
c
11    format(2i3,4(1x,g17.10))
      return
      end
c
c-----------------------------------------------------------------------------
c---- subroutine calculates the xyz coordinates of the centers each cells
c---- in the sheet by means of coordinate rotation
c-----------------------------------------------------------------------------
c
      subroutine ts_grid (is, n_col, x_loc, y_loc, z_loc)
      dimension x_loc(n_col), y_loc(n_col), z_loc(n_col)
      common /geo/ alpha(2), beta(2), aa(2), bb(2)
      common /locate/ x_top(2), y_top(2), z_top(2)
      common /cell/ na(2), nb(2), max_n
c
      naa = na(is)
      nbb = nb(is)
      iskip = (is-1)*na(1)*nb(1)
c
      del_a = aa(is) / naa
      del_b = bb(is) / nbb
      x_t = x_top(is)
      y_t = y_top(is)
      z_t = z_top(is)
      da_cos = ddcos(alpha(is))
      da_sin = sin(alpha(is))
      db_cos = ddcos(beta(is))
      db_sin = sin(beta(is))
      a_init = - 0.5*aa(is)
      b_init = - 0.5*bb(is)
      z_rise = b_init*db_sin
c   
      do i = 1, naa
         ii = (i-1)*nbb + iskip
         a_temp = a_init + (i-0.5)*del_a
         do j = 1, nbb
            ij = ii + j
            b_temp = b_init + (j-0.5)*del_b
            x_loc (ij) = x_t + a_temp*da_cos - b_temp*da_sin*db_cos
            y_loc (ij) = y_t + a_temp*da_sin + b_temp*da_cos*db_cos
            z_loc (ij) = z_t + b_temp*db_sin - z_rise
         enddo
      enddo
c
      return
      end
c
c-------------------------------------------------------------------------------
c
      complex function cole(z_l,z_h,tau,alpha)   
      complex i_om
      common /parm/ omega, pi
c
      i_om = cmplx(0.0, 1.0) * omega
      if ((z_l.lt.z_h) .or. (alpha.gt.1.0)) then
            write(*,*) '>> INVALID COLE-COLE PARAMETERS !!'
            stop
      endif
      if (tau .lt. 1.e-12) then
            cole = cmplx(1.0, 0.0) * z_l
          else
            cole = z_h + (z_l - z_h) / (1.+(i_om*tau)**alpha)
      endif
c
      return
      end
c
c-----------------------------------------------------------------------------
c
      real function ddcos(arg)
      common /parm/ omega, pi
c           
      pi2 = pi / 2.0
      diff = abs(pi2 - arg)
      if (diff. lt. 1.e-2) then
             ddcos = 0.0
          else
             ddcos = cos(arg)
      endif
      return
      end
 
 
