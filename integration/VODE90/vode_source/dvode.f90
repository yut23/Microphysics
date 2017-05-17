module dvode_module

  use dvode_type_module, only: dvode_t 
  use dvode_output_module, only: xerrwd
  use bl_types, only: dp_t
  
  implicit none

  real(dp_t), parameter :: ZERO = 0.0D0
  real(dp_t), parameter :: ONE  = 1.0D0
  real(dp_t), parameter :: HALF = 0.5D0
  real(dp_t), parameter :: TWO  = 2.0D0
  real(dp_t), parameter :: FOUR = 4.0D0
  real(dp_t), parameter :: SIX  = 6.0D0
  real(dp_t), parameter :: HUN  = 100.0D0
  real(dp_t), parameter :: THOU = 1000.0D0

  public :: dvode
  
contains

  function dumach() result(dum)
    ! ***BEGIN PROLOGUE  DUMACH
    ! ***PURPOSE  Compute the unit roundoff of the machine.
    ! ***CATEGORY  R1
    ! ***TYPE      DOUBLE PRECISION (RUMACH-S, DUMACH-D)
    ! ***KEYWORDS  MACHINE CONSTANTS
    ! ***AUTHOR  Hindmarsh, Alan C., (LLNL)
    ! ***DESCRIPTION
    !  *Usage:
    !         DOUBLE PRECISION  A, DUMACH
    !         A = DUMACH()
    ! 
    !  *Function Return Values:
    !      A : the unit roundoff of the machine.
    ! 
    !  *Description:
    !      The unit roundoff is defined as the smallest positive machine
    !      number u such that  1.0 + u .ne. 1.0.  This is computed by DUMACH
    !      in a machine-independent manner.
    ! 
    ! ***REFERENCES  (NONE)
    ! ***ROUTINES CALLED  (NONE)
    ! ***REVISION HISTORY  (YYMMDD)
    !    930216  DATE WRITTEN
    !    930818  Added SLATEC-format prologue.  (FNF)
    ! ***END PROLOGUE  DUMACH
    ! 
    ! *Internal Notes:
    ! -----------------------------------------------------------------------
    !  Subroutines/functions called by DUMACH.. None
    ! -----------------------------------------------------------------------
    ! **End
    ! 
    real(dp_t) :: U, dum
    dum = EPSILON(U)
  end function dumach
  
  subroutine dewset(N, ITOL, RTOL, ATOL, YCUR, EWT)
    ! ***BEGIN PROLOGUE  DEWSET
    ! ***SUBSIDIARY
    ! ***PURPOSE  Set error weight vector.
    ! ***TYPE      DOUBLE PRECISION (SEWSET-S, DEWSET-D)
    ! ***AUTHOR  Hindmarsh, Alan C., (LLNL)
    ! ***DESCRIPTION
    ! 
    !   This subroutine sets the error weight vector EWT according to
    !       EWT(i) = RTOL(i)*ABS(YCUR(i)) + ATOL(i),  i = 1,...,N,
    !   with the subscript on RTOL and/or ATOL possibly replaced by 1 above,
    !   depending on the value of ITOL.
    ! 
    ! ***SEE ALSO  DLSODE
    ! ***ROUTINES CALLED  (NONE)
    ! ***REVISION HISTORY  (YYMMDD)
    !    791129  DATE WRITTEN
    !    890501  Modified prologue to SLATEC/LDOC format.  (FNF)
    !    890503  Minor cosmetic changes.  (FNF)
    !    930809  Renamed to allow single/double precision versions. (ACH)
    ! ***END PROLOGUE  DEWSET
    ! **End

    integer    :: I, N, ITOL
    real(dp_t) :: RTOL(:), ATOL(:), YCUR(N), EWT(N)

    GO TO (10, 20, 30, 40), ITOL
10  CONTINUE
    do I = 1,N
       EWT(I) = RTOL(1)*ABS(YCUR(I)) + ATOL(1)
    end do
    RETURN
20  CONTINUE
    do I = 1,N
       EWT(I) = RTOL(1)*ABS(YCUR(I)) + ATOL(I)
    end do
    RETURN
30  CONTINUE
    do I = 1,N
       EWT(I) = RTOL(I)*ABS(YCUR(I)) + ATOL(1)
    end do
    RETURN
40  CONTINUE
    do I = 1,N
       EWT(I) = RTOL(I)*ABS(YCUR(I)) + ATOL(I)
    end do
    RETURN
  end subroutine dewset

  function dvnorm(N, V, W) result(dvn)
    ! ***BEGIN PROLOGUE  DVNORM
    ! ***SUBSIDIARY
    ! ***PURPOSE  Weighted root-mean-square vector norm.
    ! ***TYPE      DOUBLE PRECISION (SVNORM-S, DVNORM-D)
    ! ***AUTHOR  Hindmarsh, Alan C., (LLNL)
    ! ***DESCRIPTION
    ! 
    !   This function routine computes the weighted root-mean-square norm
    !   of the vector of length N contained in the array V, with weights
    !   contained in the array W of length N:
    !     DVNORM = SQRT( (1/N) * SUM( V(i)*W(i) )**2 )
    ! 
    ! ***SEE ALSO  DLSODE
    ! ***ROUTINES CALLED  (NONE)
    ! ***REVISION HISTORY  (YYMMDD)
    !    791129  DATE WRITTEN
    !    890501  Modified prologue to SLATEC/LDOC format.  (FNF)
    !    890503  Minor cosmetic changes.  (FNF)
    !    930809  Renamed to allow single/double precision versions. (ACH)
    ! ***END PROLOGUE  DVNORM
    ! **End
    integer    :: N, I
    real(dp_t) :: V(N), W(N), SUM, dvn

    SUM = 0.0D0
    do I = 1,N
       SUM = SUM + (V(I)*W(I))**2
    end do
    dvn = SQRT(SUM/N)
    RETURN
  end function dvnorm
  
  subroutine dvhin(N, T0, Y0, YDOT, F, RPAR, IPAR, TOUT, UROUND, &
       EWT, ITOL, ATOL, Y, TEMP, H0, NITER, IER)
    ! -----------------------------------------------------------------------
    !  Call sequence input -- N, T0, Y0, YDOT, F, RPAR, IPAR, TOUT, UROUND,
    !                         EWT, ITOL, ATOL, Y, TEMP
    !  Call sequence output -- H0, NITER, IER
    !  COMMON block variables accessed -- None
    ! 
    !  Subroutines called by DVHIN:  F
    !  Function routines called by DVHI: DVNORM
    ! -----------------------------------------------------------------------
    !  This routine computes the step size, H0, to be attempted on the
    !  first step, when the user has not supplied a value for this.
    ! 
    !  First we check that TOUT - T0 differs significantly from zero.  Then
    !  an iteration is done to approximate the initial second derivative
    !  and this is used to define h from w.r.m.s.norm(h**2 * yddot / 2) = 1.
    !  A bias factor of 1/2 is applied to the resulting h.
    !  The sign of H0 is inferred from the initial values of TOUT and T0.
    ! 
    !  Communication with DVHIN is done with the following variables:
    ! 
    !  N      = Size of ODE system, input.
    !  T0     = Initial value of independent variable, input.
    !  Y0     = Vector of initial conditions, input.
    !  YDOT   = Vector of initial first derivatives, input.
    !  F      = Name of subroutine for right-hand side f(t,y), input.
    !  RPAR, IPAR = Dummy names for user's real and integer work arrays.
    !  TOUT   = First output value of independent variable
    !  UROUND = Machine unit roundoff
    !  EWT, ITOL, ATOL = Error weights and tolerance parameters
    !                    as described in the driver routine, input.
    !  Y, TEMP = Work arrays of length N.
    !  H0     = Step size to be attempted, output.
    !  NITER  = Number of iterations (and of f evaluations) to compute H0,
    !           output.
    !  IER    = The error flag, returned with the value
    !           IER = 0  if no trouble occurred, or
    !           IER = -1 if TOUT and T0 are considered too close to proceed.
    ! -----------------------------------------------------------------------
    EXTERNAL F
    real(dp_t) :: T0, Y0(:), YDOT(:), RPAR(:), TOUT, UROUND
    real(dp_t) :: EWT(:), ATOL(:), Y(:), TEMP(:), H0
    integer    :: N, IPAR(:), ITOL, NITER, IER

    real(dp_t) :: AFI, ATOLI, DELYI, H, HG, HLB, HNEW, HRAT
    real(dp_t) :: HUB, T1, TDIST, TROUND, YDDNRM
    integer    :: I, ITER

    real(dp_t), parameter :: PT1 = 0.1D0

    NITER = 0
    write(*,*) 'TOUT = ', TOUT
    write(*,*) 'T0 = ', T0
    TDIST = ABS(TOUT - T0)
    write(*,*) 'TDIST = ', TDIST
    TROUND = UROUND*MAX(ABS(T0),ABS(TOUT))
    write(*,*) 'TDIST = ', TDIST
    write(*,*) 'TWO = ', TWO
    write(*,*) 'TROUND = ', TROUND
    IF (TDIST .LT. TWO*TROUND) GO TO 100

    ! Set a lower bound on h based on the roundoff level in T0 and TOUT. ---
    HLB = HUN*TROUND
    ! Set an upper bound on h based on TOUT-T0 and the initial Y and YDOT. -
    HUB = PT1*TDIST
    ATOLI = ATOL(1)
    do I = 1, N
       IF (ITOL .EQ. 2 .OR. ITOL .EQ. 4) ATOLI = ATOL(I)
       DELYI = PT1*ABS(Y0(I)) + ATOLI
       AFI = ABS(YDOT(I))
       IF (AFI*HUB .GT. DELYI) HUB = DELYI/AFI
    end do

    ! Set initial guess for h as geometric mean of upper and lower bounds. -
    ITER = 0
    HG = SQRT(HLB*HUB)
    ! If the bounds have crossed, exit with the mean value. ----------------
    IF (HUB .LT. HLB) THEN
       H0 = HG
       GO TO 90
    ENDIF

    ! Looping point for iteration. -----------------------------------------
50  CONTINUE
    ! Estimate the second derivative as a difference quotient in f. --------
    H = SIGN (HG, TOUT - T0)
    T1 = T0 + H
    do I = 1, N
       Y(I) = Y0(I) + H*YDOT(I)
    end do
    CALL F (N, T1, Y, TEMP, RPAR, IPAR)
    do I = 1, N
       TEMP(I) = (TEMP(I) - YDOT(I))/H
    end do
    YDDNRM = DVNORM (N, TEMP, EWT)
    ! Get the corresponding new value of h. --------------------------------
    IF (YDDNRM*HUB*HUB .GT. TWO) THEN
       HNEW = SQRT(TWO/YDDNRM)
    ELSE
       HNEW = SQRT(HG*HUB)
    ENDIF
    ITER = ITER + 1
    ! -----------------------------------------------------------------------
    !  Test the stopping conditions.
    !  Stop if the new and previous h values differ by a factor of .lt. 2.
    !  Stop if four iterations have been done.  Also, stop with previous h
    !  if HNEW/HG .gt. 2 after first iteration, as this probably means that
    !  the second derivative value is bad because of cancellation error.
    ! -----------------------------------------------------------------------
    IF (ITER .GE. 4) GO TO 80
    HRAT = HNEW/HG
    IF ( (HRAT .GT. HALF) .AND. (HRAT .LT. TWO) ) GO TO 80
    IF ( (ITER .GE. 2) .AND. (HNEW .GT. TWO*HG) ) THEN
       HNEW = HG
       GO TO 80
    ENDIF
    HG = HNEW
    GO TO 50

    ! Iteration done.  Apply bounds, bias factor, and sign.  Then exit. ----
80  H0 = HNEW*HALF
    IF (H0 .LT. HLB) H0 = HLB
    IF (H0 .GT. HUB) H0 = HUB
90  H0 = SIGN(H0, TOUT - T0)
    NITER = ITER
    IER = 0
    RETURN
    ! Error return for TOUT - T0 too small. --------------------------------
100 IER = -1
    RETURN
  end subroutine dvhin
  
  subroutine dvindy(T, K, YH, LDYH, DKY, IFLAG, dvode_state)
    ! -----------------------------------------------------------------------
    !  Call sequence input -- T, K, YH, LDYH
    !  Call sequence output -- DKY, IFLAG
    !  COMMON block variables accessed:
    !      /DVOD01/ --  H, TN, UROUND, L, N, NQ
    !      /DVOD02/ --  HU
    ! 
    !  Subroutines called by DVINDY: DSCAL, XERRWD
    !  Function routines called by DVINDY: None
    ! -----------------------------------------------------------------------
    !  DVINDY computes interpolated values of the K-th derivative of the
    !  dependent variable vector y, and stores it in DKY.  This routine
    !  is called within the package with K = 0 and T = TOUT, but may
    !  also be called by the user for any K up to the current order.
    !  (See detailed instructions in the usage documentation.)
    ! -----------------------------------------------------------------------
    !  The computed values in DKY are gotten by interpolation using the
    !  Nordsieck history array YH.  This array corresponds uniquely to a
    !  vector-valued polynomial of degree NQCUR or less, and DKY is set
    !  to the K-th derivative of this polynomial at T.
    !  The formula for DKY is:
    !               q
    !   DKY(i)  =  sum  c(j,K) * (T - TN)**(j-K) * H**(-j) * YH(i,j+1)
    !              j=K
    !  where  c(j,K) = j*(j-1)*...*(j-K+1), q = NQCUR, TN = TCUR, H = HCUR.
    !  The quantities  NQ = NQCUR, L = NQ+1, N, TN, and H are
    !  communicated by COMMON.  The above sum is done in reverse order.
    !  IFLAG is returned negative if either K or T is out of bounds.
    ! 
    !  Discussion above and comments in driver explain all variables.
    ! -----------------------------------------------------------------------
    !
    
    type(dvode_t) :: dvode_state
    real(dp_t) :: T, YH(LDYH, dvode_state % LMAX), DKY(:)
    integer    :: K, LDYH, IFLAG


    real(dp_t) :: C, R, S, TFUZZ, TN1, TP
    integer    :: I, IC, J, JB, JB2, JJ, JJ1, JP1
    character (len=80) :: MSG

    IFLAG = 0
    IF (K .LT. 0 .OR. K .GT. dvode_state % NQ) GO TO 80
    TFUZZ = HUN * dvode_state % UROUND * (dvode_state % TN + dvode_state % HU)
    TP = dvode_state % TN - dvode_state % HU - TFUZZ
    TN1 = dvode_state % TN + TFUZZ
    IF ((T-TP)*(T-TN1) .GT. ZERO) GO TO 90

    S = (T - dvode_state % TN)/dvode_state % H
    IC = 1
    IF (K .EQ. 0) GO TO 15
    JJ1 = dvode_state % L - K
    do JJ = JJ1, dvode_state % NQ
       IC = IC*JJ
    end do
15  C = REAL(IC)
    do I = 1, dvode_state % N
       DKY(I) = C*YH(I,dvode_state % L)
    end do
    IF (K .EQ. dvode_state % NQ) GO TO 55
    JB2 = dvode_state % NQ - K
    do JB = 1, JB2
       J = dvode_state % NQ - JB
       JP1 = J + 1
       IC = 1
       IF (K .EQ. 0) GO TO 35
       JJ1 = JP1 - K
       do JJ = JJ1, J
          IC = IC*JJ
       end do
35     C = REAL(IC)
       do I = 1, dvode_state % N
          DKY(I) = C*YH(I,JP1) + S*DKY(I)
       end do
    end do
    IF (K .EQ. 0) RETURN
55  R = dvode_state % H**(-K)
    CALL DSCAL (dvode_state % N, R, DKY, 1)
    RETURN

80  MSG = 'DVINDY-- K (=I1) illegal      '
    CALL XERRWD (MSG, 30, 51, 1, 1, K, 0, 0, ZERO, ZERO)
    IFLAG = -1
    RETURN
90  MSG = 'DVINDY-- T (=R1) illegal      '
    CALL XERRWD (MSG, 30, 52, 1, 0, 0, 0, 1, T, ZERO)
    MSG='      T not in interval TCUR - HU (= R1) to TCUR (=R2)      '
    CALL XERRWD (MSG, 60, 52, 1, 0, 0, 0, 2, TP, dvode_state % TN)
    IFLAG = -2
    RETURN
  end subroutine dvindy
  
  subroutine dvode(F, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK, &
       ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, MF, &
       RPAR, IPAR)
    
!    use rpar_indices, only: n_rpar_comps, n_ipar_comps

    external F, JAC

    integer    :: NEQ, ITOL, ITASK, ISTATE, IOPT, LRW, LIW, MF
    integer    :: IWORK(LIW)
!    integer    :: IPAR(n_ipar_comps)
    integer    :: IPAR(:)    
    real(dp_t) :: T, TOUT
    real(dp_t) :: Y(NEQ), RTOL(NEQ), ATOL(NEQ), RWORK(LRW)
!    real(dp_t) :: RPAR(n_rpar_comps)
    real(dp_t) :: RPAR(:)

    logical    :: IHIT
    real(dp_t) :: ATOLI, BIG, EWTI, H0, HMAX, HMX
    real(dp_t) :: RH, RTOLI, SIZE, TCRIT, TNEXT, TOLSF, TP
    integer    :: I, IER, IFLAG, IMXER, JCO, KGO, LENIW, LENJ, LENP, LENRW
    integer    :: LENWM, LF0, MBAND, MFA, ML, MU, NITER
    integer    :: NSLAST
    integer, dimension(2) :: MORD = [12, 5]
    character (len=80) :: MSG

    type(dvode_t) :: dvode_state

    ! Parameter declarations
    integer, parameter :: MXSTP0 = 500
    integer, parameter :: MXHNL0 = 10
    real(dp_t), parameter :: PT2 = 0.2D0
    real(dp_t), parameter :: HUN = 100.0D0
    
    ! -----------------------------------------------------------------------
    !  Block A.
    !  This code block is executed on every call.
    !  It tests ISTATE and ITASK for legality and branches appropriately.
    !  If ISTATE .gt. 1 but the flag INIT shows that initialization has
    !  not yet been done, an error return occurs.
    !  If ISTATE = 1 and TOUT = T, return immediately.
    ! -----------------------------------------------------------------------

    IF (ISTATE .LT. 1 .OR. ISTATE .GT. 3) GO TO 601
    IF (ITASK .LT. 1 .OR. ITASK .GT. 5) GO TO 602
    IF (ISTATE .EQ. 1) GO TO 10
    IF (dvode_state % INIT .NE. 1) GO TO 603
    IF (ISTATE .EQ. 2) GO TO 200
    GO TO 20
10  dvode_state % INIT = 0
    IF (TOUT .EQ. T) RETURN
    
    ! -----------------------------------------------------------------------
    !  Block B.
    !  The next code block is executed for the initial call (ISTATE = 1),
    !  or for a continuation call with parameter changes (ISTATE = 3).
    !  It contains checking of all input and various initializations.
    !
    !  First check legality of the non-optional input NEQ, ITOL, IOPT,
    !  MF, ML, and MU.
    ! -----------------------------------------------------------------------
    
20  IF (NEQ .LE. 0) GO TO 604
    IF (ISTATE .EQ. 1) GO TO 25
    IF (NEQ .GT. dvode_state % N) GO TO 605
25  dvode_state % N = NEQ
    IF (ITOL .LT. 1 .OR. ITOL .GT. 4) GO TO 606
    IF (IOPT .LT. 0 .OR. IOPT .GT. 1) GO TO 607
    dvode_state % JSV = SIGN(1,MF)
    MFA = ABS(MF)
    dvode_state % METH = MFA/10
    dvode_state % MITER = MFA - 10*dvode_state % METH
    IF (dvode_state % METH .LT. 1 .OR. dvode_state % METH .GT. 2) GO TO 608
    IF (dvode_state % MITER .LT. 0 .OR. dvode_state % MITER .GT. 5) GO TO 608
    IF (dvode_state % MITER .LE. 3) GO TO 30
    ML = IWORK(1)
    MU = IWORK(2)
    IF (ML .LT. 0 .OR. ML .GE. dvode_state % N) GO TO 609
    IF (MU .LT. 0 .OR. MU .GE. dvode_state % N) GO TO 610
30  CONTINUE

    ! Next process and check the optional input. ---------------------------
      
    IF (IOPT .EQ. 1) GO TO 40
    dvode_state % MAXORD = MORD(dvode_state % METH)
    dvode_state % MXSTEP = MXSTP0
    dvode_state % MXHNIL = MXHNL0
    IF (ISTATE .EQ. 1) H0 = ZERO
    dvode_state % HMXI = ZERO
    dvode_state % HMIN = ZERO
    GO TO 60
40  dvode_state % MAXORD = IWORK(5)
    IF (dvode_state % MAXORD .LT. 0) GO TO 611
    IF (dvode_state % MAXORD .EQ. 0) dvode_state % MAXORD = 100
    dvode_state % MAXORD = MIN(dvode_state % MAXORD,MORD(dvode_state % METH))
    dvode_state % MXSTEP = IWORK(6)
    IF (dvode_state % MXSTEP .LT. 0) GO TO 612
    IF (dvode_state % MXSTEP .EQ. 0) dvode_state % MXSTEP = MXSTP0
    dvode_state % MXHNIL = IWORK(7)
    IF (dvode_state % MXHNIL .LT. 0) GO TO 613
    !      EDIT 07/16/2016 -- see comments above about MXHNIL
    !      IF (MXHNIL .EQ. 0) MXHNIL = MXHNL0
    IF (ISTATE .NE. 1) GO TO 50
    H0 = RWORK(5)
    IF ((TOUT - T)*H0 .LT. ZERO) GO TO 614
50  HMAX = RWORK(6)
    IF (HMAX .LT. ZERO) GO TO 615
    dvode_state % HMXI = ZERO
    IF (HMAX .GT. ZERO) dvode_state % HMXI = ONE/HMAX
    dvode_state % HMIN = RWORK(7)
    IF (dvode_state % HMIN .LT. ZERO) GO TO 616
    
    ! -----------------------------------------------------------------------
    !  Set work array pointers and check lengths LRW and LIW.
    !  Pointers to segments of RWORK and IWORK are named by prefixing L to
    !  the name of the segment.  E.g., the segment YH starts at RWORK(LYH).
    !  Segments of RWORK (in order) are denoted  YH, WM, EWT, SAVF, ACOR.
    !  Within WM, LOCJS is the location of the saved Jacobian (JSV .gt. 0).
    ! -----------------------------------------------------------------------
    
60  dvode_state % LYH = 21
    IF (ISTATE .EQ. 1) dvode_state % NYH = dvode_state % N
    dvode_state % LWM = dvode_state % LYH + (dvode_state % MAXORD + 1)*dvode_state % NYH
    JCO = MAX(0,dvode_state % JSV)
    IF (dvode_state % MITER .EQ. 0) LENWM = 0
    IF (dvode_state % MITER .EQ. 1 .OR. dvode_state % MITER .EQ. 2) THEN
       LENWM = 2 + (1 + JCO)*dvode_state % N*dvode_state % N
       dvode_state % LOCJS = dvode_state % N*dvode_state % N + 3
    ENDIF
    IF (dvode_state % MITER .EQ. 3) LENWM = 2 + dvode_state % N
    IF (dvode_state % MITER .EQ. 4 .OR. dvode_state % MITER .EQ. 5) THEN
       MBAND = ML + MU + 1
       LENP = (MBAND + ML)*dvode_state % N
       LENJ = MBAND*dvode_state % N
       LENWM = 2 + LENP + JCO*LENJ
       dvode_state % LOCJS = LENP + 3
    ENDIF
    dvode_state % LEWT = dvode_state % LWM + LENWM
    dvode_state % LSAVF = dvode_state % LEWT + dvode_state % N
    dvode_state % LACOR = dvode_state % LSAVF + dvode_state % N
    LENRW = dvode_state % LACOR + dvode_state % N - 1
    IWORK(17) = LENRW
    dvode_state % LIWM = 1
    LENIW = 30 + dvode_state % N
    IF (dvode_state % MITER .EQ. 0 .OR. dvode_state % MITER .EQ. 3) LENIW = 30
    IWORK(18) = LENIW
    IF (LENRW .GT. LRW) GO TO 617
    IF (LENIW .GT. LIW) GO TO 618
    ! Check RTOL and ATOL for legality. ------------------------------------
    RTOLI = RTOL(1)
    ATOLI = ATOL(1)
    do I = 1,dvode_state % N
       IF (ITOL .GE. 3) RTOLI = RTOL(I)
       IF (ITOL .EQ. 2 .OR. ITOL .EQ. 4) ATOLI = ATOL(I)
       IF (RTOLI .LT. ZERO) GO TO 619
       IF (ATOLI .LT. ZERO) GO TO 620
    end do
    IF (ISTATE .EQ. 1) GO TO 100
    ! If ISTATE = 3, set flag to signal parameter changes to DVSTEP. -------
    dvode_state % JSTART = -1
    IF (dvode_state % NQ .LE. dvode_state % MAXORD) GO TO 90
    ! MAXORD was reduced below NQ.  Copy YH(*,MAXORD+2) into SAVF. ---------
    CALL DCOPY (dvode_state % N, RWORK(dvode_state % LWM), 1, RWORK(dvode_state % LSAVF), 1)
    ! Reload WM(1) = RWORK(LWM), since LWM may have changed. ---------------
90  IF (dvode_state % MITER .GT. 0) RWORK(dvode_state % LWM) = SQRT(dvode_state % UROUND)
    GO TO 200

    ! -----------------------------------------------------------------------
    !  Block C.
    !  The next block is for the initial call only (ISTATE = 1).
    !  It contains all remaining initializations, the initial call to F,
    !  and the calculation of the initial step size.
    !  The error weights in EWT are inverted after being loaded.
    ! -----------------------------------------------------------------------
       
100 dvode_state % UROUND = DUMACH()
    dvode_state % TN = T
    IF (ITASK .NE. 4 .AND. ITASK .NE. 5) GO TO 110
    TCRIT = RWORK(1)
    IF ((TCRIT - TOUT)*(TOUT - T) .LT. ZERO) GO TO 625
    if (H0 .NE. ZERO .AND. (T + H0 - TCRIT)*H0 .GT. ZERO) then
       H0 = TCRIT - T
    end if
110 dvode_state % JSTART = 0
    IF (dvode_state % MITER .GT. 0) RWORK(dvode_state % LWM) = SQRT(dvode_state % UROUND)
    dvode_state % CCMXJ = PT2
    dvode_state % MSBJ = 50
    dvode_state % NHNIL = 0
    dvode_state % NST = 0
    dvode_state % NJE = 0
    dvode_state % NNI = 0
    dvode_state % NCFN = 0
    dvode_state % NETF = 0
    dvode_state % NLU = 0
    dvode_state % NSLJ = 0
    NSLAST = 0
    dvode_state % HU = ZERO
    dvode_state % NQU = 0
    ! Initial call to F.  (LF0 points to YH(*,2).) -------------------------
    LF0 = dvode_state % LYH + dvode_state % NYH
    CALL F (dvode_state % N, T, Y, RWORK(LF0), RPAR, IPAR)
    dvode_state % NFE = 1
    ! Load the initial value vector in YH. ---------------------------------
    CALL DCOPY (dvode_state % N, Y, 1, RWORK(dvode_state % LYH), 1)
    ! Load and invert the EWT array.  (H is temporarily set to 1.0.) -------
    dvode_state % NQ = 1
    dvode_state % H = ONE
    CALL DEWSET (dvode_state % N, ITOL, RTOL, ATOL, RWORK(dvode_state % LYH), RWORK(dvode_state % LEWT))
    do I = 1,dvode_state % N
       IF (RWORK(I+dvode_state % LEWT-1) .LE. ZERO) GO TO 621
       RWORK(I+dvode_state % LEWT-1) = ONE/RWORK(I+dvode_state % LEWT-1)
    end do
    IF (H0 .NE. ZERO) GO TO 180
    ! Call DVHIN to set initial step size H0 to be attempted. --------------
    CALL DVHIN (dvode_state % N, T, &
         RWORK(dvode_state % LYH:dvode_state % LYH + dvode_state % N - 1), &
         RWORK(LF0:LF0 + dvode_state % N - 1), F, RPAR, IPAR, TOUT, &
         dvode_state % UROUND, &
         RWORK(dvode_state % LEWT:dvode_state % LEWT + dvode_state % N - 1), &
         ITOL, ATOL, Y, &
         RWORK(dvode_state % LACOR:dvode_state % LACOR + dvode_state % N - 1), &
         H0, NITER, IER)
    dvode_state % NFE = dvode_state % NFE + NITER
    IF (IER .NE. 0) GO TO 622
    ! Adjust H0 if necessary to meet HMAX bound. ---------------------------
180 RH = ABS(H0)*dvode_state % HMXI
    IF (RH .GT. ONE) H0 = H0/RH
    ! Load H with H0 and scale YH(*,2) by H0. ------------------------------
    dvode_state % H = H0
    CALL DSCAL (dvode_state % N, H0, RWORK(LF0), 1)
    GO TO 270
    
    ! -----------------------------------------------------------------------
    !  Block D.
    !  The next code block is for continuation calls only (ISTATE = 2 or 3)
    !  and is to check stop conditions before taking a step.
    ! -----------------------------------------------------------------------
    
200 NSLAST = dvode_state % NST
    dvode_state % KUTH = 0
    GO TO (210, 250, 220, 230, 240), ITASK
210 IF ((dvode_state % TN - TOUT) * dvode_state % H .LT. ZERO) GO TO 250
    CALL DVINDY (TOUT, 0, RWORK(dvode_state % LYH), dvode_state % NYH, Y, IFLAG, dvode_state)
    IF (IFLAG .NE. 0) GO TO 627
    T = TOUT
    GO TO 420
220 TP = dvode_state % TN - dvode_state % HU * (ONE + HUN * dvode_state % UROUND)
    IF ((TP - TOUT) * dvode_state % H .GT. ZERO) GO TO 623
    IF ((dvode_state % TN - TOUT) * dvode_state % H .LT. ZERO) GO TO 250
    GO TO 400
230 TCRIT = RWORK(1)
    IF ((dvode_state % TN - TCRIT) * dvode_state % H .GT. ZERO) GO TO 624
    IF ((TCRIT - TOUT) * dvode_state % H .LT. ZERO) GO TO 625
    IF ((dvode_state % TN - TOUT) * dvode_state % H .LT. ZERO) GO TO 245
    CALL DVINDY (TOUT, 0, RWORK(dvode_state % LYH), dvode_state % NYH, Y, IFLAG, dvode_state)
    IF (IFLAG .NE. 0) GO TO 627
    T = TOUT
    GO TO 420
240 TCRIT = RWORK(1)
    IF ((dvode_state % TN - TCRIT) * dvode_state % H .GT. ZERO) GO TO 624
245 HMX = ABS(dvode_state % TN) + ABS(dvode_state % H)
    IHIT = ABS(dvode_state % TN - TCRIT) .LE. HUN * dvode_state % UROUND * HMX
    IF (IHIT) GO TO 400
    TNEXT = dvode_state % TN + dvode_state % HNEW*(ONE + FOUR * dvode_state % UROUND)
    IF ((TNEXT - TCRIT) * dvode_state % H .LE. ZERO) GO TO 250
    dvode_state % H = (TCRIT - dvode_state % TN)*(ONE - FOUR * dvode_state % UROUND)
    dvode_state % KUTH = 1
    
    ! -----------------------------------------------------------------------
    !  Block E.
    !  The next block is normally executed for all calls and contains
    !  the call to the one-step core integrator DVSTEP.
    !
    !  This is a looping point for the integration steps.
    !
    !  First check for too many steps being taken, update EWT (if not at
    !  start of problem), check for too much accuracy being requested, and
    !  check for H below the roundoff level in T.
    ! -----------------------------------------------------------------------
    
250 CONTINUE
    IF ((dvode_state % NST-NSLAST) .GE. dvode_state % MXSTEP) GO TO 500
    CALL DEWSET (dvode_state % N, ITOL, RTOL, ATOL, RWORK(dvode_state % LYH), RWORK(dvode_state % LEWT))
    do I = 1,dvode_state % N
       IF (RWORK(I+dvode_state % LEWT-1) .LE. ZERO) GO TO 510
       RWORK(I+dvode_state % LEWT-1) = ONE/RWORK(I+dvode_state % LEWT-1)
    end do
270 TOLSF = dvode_state % UROUND * DVNORM (dvode_state % N, RWORK(dvode_state % LYH), RWORK(dvode_state % LEWT))
    IF (TOLSF .LE. ONE) GO TO 280
    TOLSF = TOLSF*TWO
    IF (dvode_state % NST .EQ. 0) GO TO 626
    GO TO 520
280 IF ((dvode_state % TN + dvode_state % H) .NE. dvode_state % TN) GO TO 290
    dvode_state % NHNIL = dvode_state % NHNIL + 1
    IF (dvode_state % NHNIL .GT. dvode_state % MXHNIL) GO TO 290
    MSG = 'DVODE--  Warning: internal T (=R1) and H (=R2) are'
    CALL XERRWD (MSG, 50, 101, 1, 0, 0, 0, 0, ZERO, ZERO)
    MSG='      such that in the machine, T + H = T on the next step  '
    CALL XERRWD (MSG, 60, 101, 1, 0, 0, 0, 0, ZERO, ZERO)
    MSG = '      (H = step size). solver will continue anyway'
    CALL XERRWD (MSG, 50, 101, 1, 0, 0, 0, 2, dvode_state % TN, dvode_state % H)
    IF (dvode_state % NHNIL .LT. dvode_state % MXHNIL) GO TO 290
    MSG = 'DVODE--  Above warning has been issued I1 times.  '
    CALL XERRWD (MSG, 50, 102, 1, 0, 0, 0, 0, ZERO, ZERO)
    MSG = '      it will not be issued again for this problem'
    CALL XERRWD (MSG, 50, 102, 1, 1, dvode_state % MXHNIL, 0, 0, ZERO, ZERO)
290 CONTINUE
    
    ! -----------------------------------------------------------------------
    !  CALL DVSTEP (Y, YH, NYH, YH, EWT, SAVF, VSAV, ACOR,
    !               WM, IWM, F, JAC, F, DVNLSD, RPAR, IPAR)
    ! -----------------------------------------------------------------------
    
    CALL DVSTEP(Y, RWORK(dvode_state % LYH), dvode_state % NYH, &
         RWORK(dvode_state % LYH:dvode_state % LYH + dvode_state % NYH * dvode_state % LMAX - 1), &
         RWORK(dvode_state % LEWT:dvode_state % LEWT + NEQ - 1), &
         RWORK(dvode_state % LSAVF:dvode_state % LSAVF + NEQ - 1), &
         Y, &
         RWORK(dvode_state % LACOR:dvode_state % LACOR + NEQ - 1), &
         RWORK(dvode_state % LWM:dvode_state % LEWT - 1), &
         IWORK, F, JAC, F, DVNLSD, RPAR, IPAR, dvode_state)
    KGO = 1 - dvode_state % KFLAG
    ! Branch on KFLAG.  Note: In this version, KFLAG can not be set to -3.
    !  KFLAG .eq. 0,   -1,  -2
    GO TO (300, 530, 540), KGO
    
    ! -----------------------------------------------------------------------
    !  Block F.
    !  The following block handles the case of a successful return from the
    !  core integrator (KFLAG = 0).  Test for stop conditions.
    ! -----------------------------------------------------------------------
    
300 dvode_state % INIT = 1
    dvode_state % KUTH = 0
    GO TO (310, 400, 330, 340, 350), ITASK
    ! ITASK = 1.  If TOUT has been reached, interpolate. -------------------
310 IF ((dvode_state % TN - TOUT) * dvode_state % H .LT. ZERO) GO TO 250
    CALL DVINDY (TOUT, 0, RWORK(dvode_state % LYH), dvode_state % NYH, Y, IFLAG, dvode_state)
    T = TOUT
    GO TO 420
    ! ITASK = 3.  Jump to exit if TOUT was reached. ------------------------
330 IF ((dvode_state % TN - TOUT) * dvode_state % H .GE. ZERO) GO TO 400
    GO TO 250
    ! ITASK = 4.  See if TOUT or TCRIT was reached.  Adjust H if necessary.
340 IF ((dvode_state % TN - TOUT) * dvode_state % H .LT. ZERO) GO TO 345
    CALL DVINDY (TOUT, 0, RWORK(dvode_state % LYH), dvode_state % NYH, Y, IFLAG, dvode_state)
    T = TOUT
    GO TO 420
345 HMX = ABS(dvode_state % TN) + ABS(dvode_state % H)
    IHIT = ABS(dvode_state % TN - TCRIT) .LE. HUN * dvode_state % UROUND * HMX
    IF (IHIT) GO TO 400
    TNEXT = dvode_state % TN + dvode_state % HNEW*(ONE + FOUR * dvode_state % UROUND)
    IF ((TNEXT - TCRIT) * dvode_state % H .LE. ZERO) GO TO 250
    dvode_state % H = (TCRIT - dvode_state % TN)*(ONE - FOUR * dvode_state % UROUND)
    dvode_state % KUTH = 1
    GO TO 250
    ! ITASK = 5.  See if TCRIT was reached and jump to exit. ---------------
350 HMX = ABS(dvode_state % TN) + ABS(dvode_state % H)
    IHIT = ABS(dvode_state % TN - TCRIT) .LE. HUN * dvode_state % UROUND * HMX
    
    ! -----------------------------------------------------------------------
    !  Block G.
    !  The following block handles all successful returns from DVODE.
    !  If ITASK .ne. 1, Y is loaded from YH and T is set accordingly.
    !  ISTATE is set to 2, and the optional output is loaded into the work
    !  arrays before returning.
    ! -----------------------------------------------------------------------
    
400 CONTINUE
    CALL DCOPY (dvode_state % N, RWORK(dvode_state % LYH), 1, Y, 1)
    T = dvode_state % TN
    IF (ITASK .NE. 4 .AND. ITASK .NE. 5) GO TO 420
    IF (IHIT) T = TCRIT
420 ISTATE = 2
    RWORK(11) = dvode_state % HU
    RWORK(12) = dvode_state % HNEW
    RWORK(13) = dvode_state % TN
    IWORK(11) = dvode_state % NST
    IWORK(12) = dvode_state % NFE
    IWORK(13) = dvode_state % NJE
    IWORK(14) = dvode_state % NQU
    IWORK(15) = dvode_state % NEWQ
    IWORK(19) = dvode_state % NLU
    IWORK(20) = dvode_state % NNI
    IWORK(21) = dvode_state % NCFN
    IWORK(22) = dvode_state % NETF

    return
    
    ! -----------------------------------------------------------------------
    !  Block H.
    !  The following block handles all unsuccessful returns other than
    !  those for illegal input.  First the error message routine is called.
    !  if there was an error test or convergence test failure, IMXER is set.
    !  Then Y is loaded from YH, and T is set to TN.
    !  The optional output is loaded into the work arrays before returning.
    ! -----------------------------------------------------------------------
    
    ! The maximum number of steps was taken before reaching TOUT. ----------
500 MSG = 'DVODE--  At current T (=R1), MXSTEP (=I1) steps   '
    CALL XERRWD (MSG, 50, 201, 1, 0, 0, 0, 0, ZERO, ZERO)
    MSG = '      taken on this call before reaching TOUT     '
    CALL XERRWD (MSG, 50, 201, 1, 1, dvode_state % MXSTEP, 0, 1, dvode_state % TN, ZERO)
    ISTATE = -1
    GO TO 580
    ! EWT(i) .le. 0.0 for some i (not at start of problem). ----------------
510 EWTI = RWORK(dvode_state % LEWT+I-1)
    MSG = 'DVODE--  At T (=R1), EWT(I1) has become R2 .le. 0.'
    CALL XERRWD (MSG, 50, 202, 1, 1, I, 0, 2, dvode_state % TN, EWTI)
    ISTATE = -6
    GO TO 580
    ! Too much accuracy requested for machine precision. -------------------
520 MSG = 'DVODE--  At T (=R1), too much accuracy requested  '
    CALL XERRWD (MSG, 50, 203, 1, 0, 0, 0, 0, ZERO, ZERO)
    MSG = '      for precision of machine:   see TOLSF (=R2) '
    CALL XERRWD (MSG, 50, 203, 1, 0, 0, 0, 2, dvode_state % TN, TOLSF)
    RWORK(14) = TOLSF
    ISTATE = -2
    GO TO 580
    ! KFLAG = -1.  Error test failed repeatedly or with ABS(H) = HMIN. -----
530 MSG = 'DVODE--  At T(=R1) and step size H(=R2), the error'
    CALL XERRWD (MSG, 50, 204, 1, 0, 0, 0, 0, ZERO, ZERO)
    MSG = '      test failed repeatedly or with abs(H) = HMIN'
    CALL XERRWD (MSG, 50, 204, 1, 0, 0, 0, 2, dvode_state % TN, dvode_state % H)
    ISTATE = -4
    GO TO 560
    ! KFLAG = -2.  Convergence failed repeatedly or with ABS(H) = HMIN. ----
540 MSG = 'DVODE--  At T (=R1) and step size H (=R2), the    '
    CALL XERRWD (MSG, 50, 205, 1, 0, 0, 0, 0, ZERO, ZERO)
    MSG = '      corrector convergence failed repeatedly     '
    CALL XERRWD (MSG, 50, 205, 1, 0, 0, 0, 0, ZERO, ZERO)
    MSG = '      or with abs(H) = HMIN   '
    CALL XERRWD (MSG, 30, 205, 1, 0, 0, 0, 2, dvode_state % TN, dvode_state % H)
    ISTATE = -5
    ! Compute IMXER if relevant. -------------------------------------------
560 BIG = ZERO
    IMXER = 1
    do I = 1,dvode_state % N
       SIZE = ABS(RWORK(I+dvode_state % LACOR-1)*RWORK(I+dvode_state % LEWT-1))
       IF (BIG .GE. SIZE) exit
       BIG = SIZE
       IMXER = I
    end do
    IWORK(16) = IMXER
    ! Set Y vector, T, and optional output. --------------------------------
580 CONTINUE
    CALL DCOPY (dvode_state % N, RWORK(dvode_state % LYH), 1, Y, 1)
    T = dvode_state % TN
    RWORK(11) = dvode_state % HU
    RWORK(12) = dvode_state % H
    RWORK(13) = dvode_state % TN
    IWORK(11) = dvode_state % NST
    IWORK(12) = dvode_state % NFE
    IWORK(13) = dvode_state % NJE
    IWORK(14) = dvode_state % NQU
    IWORK(15) = dvode_state % NQ
    IWORK(19) = dvode_state % NLU
    IWORK(20) = dvode_state % NNI
    IWORK(21) = dvode_state % NCFN
    IWORK(22) = dvode_state % NETF

    return
    
    ! -----------------------------------------------------------------------
    !  Block I.
    !  The following block handles all error returns due to illegal input
    !  (ISTATE = -3), as detected before calling the core integrator.
    !  First the error message routine is called.   If the illegal input
    !  is a negative ISTATE, the run is aborted (apparent infinite loop).
    ! -----------------------------------------------------------------------
    
601 MSG = 'DVODE--  ISTATE (=I1) illegal '
    CALL XERRWD (MSG, 30, 1, 1, 1, ISTATE, 0, 0, ZERO, ZERO)
    IF (ISTATE .LT. 0) GO TO 800
    GO TO 700
602 MSG = 'DVODE--  ITASK (=I1) illegal  '
    CALL XERRWD (MSG, 30, 2, 1, 1, ITASK, 0, 0, ZERO, ZERO)
    GO TO 700
603 MSG='DVODE--  ISTATE (=I1) .gt. 1 but DVODE not initialized      '
    CALL XERRWD (MSG, 60, 3, 1, 1, ISTATE, 0, 0, ZERO, ZERO)
    GO TO 700
604 MSG = 'DVODE--  NEQ (=I1) .lt. 1     '
    CALL XERRWD (MSG, 30, 4, 1, 1, NEQ, 0, 0, ZERO, ZERO)
    GO TO 700
605 MSG = 'DVODE--  ISTATE = 3 and NEQ increased (I1 to I2)  '
    CALL XERRWD (MSG, 50, 5, 1, 2, dvode_state % N, NEQ, 0, ZERO, ZERO)
    GO TO 700
606 MSG = 'DVODE--  ITOL (=I1) illegal   '
    CALL XERRWD (MSG, 30, 6, 1, 1, ITOL, 0, 0, ZERO, ZERO)
    GO TO 700
607 MSG = 'DVODE--  IOPT (=I1) illegal   '
    CALL XERRWD (MSG, 30, 7, 1, 1, IOPT, 0, 0, ZERO, ZERO)
    GO TO 700
608 MSG = 'DVODE--  MF (=I1) illegal     '
    CALL XERRWD (MSG, 30, 8, 1, 1, MF, 0, 0, ZERO, ZERO)
    GO TO 700
609 MSG = 'DVODE--  ML (=I1) illegal:  .lt.0 or .ge.NEQ (=I2)'
    CALL XERRWD (MSG, 50, 9, 1, 2, ML, NEQ, 0, ZERO, ZERO)
    GO TO 700
610 MSG = 'DVODE--  MU (=I1) illegal:  .lt.0 or .ge.NEQ (=I2)'
    CALL XERRWD (MSG, 50, 10, 1, 2, MU, NEQ, 0, ZERO, ZERO)
    GO TO 700
611 MSG = 'DVODE--  MAXORD (=I1) .lt. 0  '
    CALL XERRWD (MSG, 30, 11, 1, 1, dvode_state % MAXORD, 0, 0, ZERO, ZERO)
    GO TO 700
612 MSG = 'DVODE--  MXSTEP (=I1) .lt. 0  '
    CALL XERRWD (MSG, 30, 12, 1, 1, dvode_state % MXSTEP, 0, 0, ZERO, ZERO)
    GO TO 700
613 MSG = 'DVODE--  MXHNIL (=I1) .lt. 0  '
    CALL XERRWD (MSG, 30, 13, 1, 1, dvode_state % MXHNIL, 0, 0, ZERO, ZERO)
    GO TO 700
614 MSG = 'DVODE--  TOUT (=R1) behind T (=R2)      '
    CALL XERRWD (MSG, 40, 14, 1, 0, 0, 0, 2, TOUT, T)
    MSG = '      integration direction is given by H0 (=R1)  '
    CALL XERRWD (MSG, 50, 14, 1, 0, 0, 0, 1, H0, ZERO)
    GO TO 700
615 MSG = 'DVODE--  HMAX (=R1) .lt. 0.0  '
    CALL XERRWD (MSG, 30, 15, 1, 0, 0, 0, 1, HMAX, ZERO)
    GO TO 700
616 MSG = 'DVODE--  HMIN (=R1) .lt. 0.0  '
    CALL XERRWD (MSG, 30, 16, 1, 0, 0, 0, 1, dvode_state % HMIN, ZERO)
    GO TO 700
617 CONTINUE
    MSG='DVODE--  RWORK length needed, LENRW (=I1), exceeds LRW (=I2)'
    CALL XERRWD (MSG, 60, 17, 1, 2, LENRW, LRW, 0, ZERO, ZERO)
    GO TO 700
618 CONTINUE
    MSG='DVODE--  IWORK length needed, LENIW (=I1), exceeds LIW (=I2)'
    CALL XERRWD (MSG, 60, 18, 1, 2, LENIW, LIW, 0, ZERO, ZERO)
    GO TO 700
619 MSG = 'DVODE--  RTOL(I1) is R1 .lt. 0.0        '
    CALL XERRWD (MSG, 40, 19, 1, 1, I, 0, 1, RTOLI, ZERO)
    GO TO 700
620 MSG = 'DVODE--  ATOL(I1) is R1 .lt. 0.0        '
    CALL XERRWD (MSG, 40, 20, 1, 1, I, 0, 1, ATOLI, ZERO)
    GO TO 700
621 EWTI = RWORK(dvode_state % LEWT+I-1)
    MSG = 'DVODE--  EWT(I1) is R1 .le. 0.0         '
    CALL XERRWD (MSG, 40, 21, 1, 1, I, 0, 1, EWTI, ZERO)
    GO TO 700
622 CONTINUE
    MSG='DVODE--  TOUT (=R1) too close to T(=R2) to start integration'
    CALL XERRWD (MSG, 60, 22, 1, 0, 0, 0, 2, TOUT, T)
    GO TO 700
623 CONTINUE
    MSG='DVODE--  ITASK = I1 and TOUT (=R1) behind TCUR - HU (= R2)  '
    CALL XERRWD (MSG, 60, 23, 1, 1, ITASK, 0, 2, TOUT, TP)
    GO TO 700
624 CONTINUE
    MSG='DVODE--  ITASK = 4 or 5 and TCRIT (=R1) behind TCUR (=R2)   '
    CALL XERRWD (MSG, 60, 24, 1, 0, 0, 0, 2, TCRIT, dvode_state % TN)
    GO TO 700
625 CONTINUE
    MSG='DVODE--  ITASK = 4 or 5 and TCRIT (=R1) behind TOUT (=R2)   '
    CALL XERRWD (MSG, 60, 25, 1, 0, 0, 0, 2, TCRIT, TOUT)
    GO TO 700
626 MSG = 'DVODE--  At start of problem, too much accuracy   '
    CALL XERRWD (MSG, 50, 26, 1, 0, 0, 0, 0, ZERO, ZERO)
    MSG='      requested for precision of machine:   see TOLSF (=R1) '
    CALL XERRWD (MSG, 60, 26, 1, 0, 0, 0, 1, TOLSF, ZERO)
    RWORK(14) = TOLSF
    GO TO 700
627 MSG='DVODE--  Trouble from DVINDY.  ITASK = I1, TOUT = R1.       '
    CALL XERRWD (MSG, 60, 27, 1, 1, ITASK, 0, 1, TOUT, ZERO)
    
700 CONTINUE
    ISTATE = -3

    return
    
800 MSG = 'DVODE--  Run aborted:  apparent infinite loop     '
    CALL XERRWD (MSG, 50, 303, 2, 0, 0, 0, 0, ZERO, ZERO)

    return
  end subroutine dvode

  subroutine dvsol(WM, IWM, X, IERSL, dvode_state)
    ! -----------------------------------------------------------------------
    !  Call sequence input -- WM, IWM, X
    !  Call sequence output -- X, IERSL
    !  COMMON block variables accessed:
    !      /DVOD01/ -- H, RL1, MITER, N
    ! 
    !  Subroutines called by DVSOL: DGESL, DGBSL
    !  Function routines called by DVSOL: None
    ! -----------------------------------------------------------------------
    !  This routine manages the solution of the linear system arising from
    !  a chord iteration.  It is called if MITER .ne. 0.
    !  If MITER is 1 or 2, it calls DGESL to accomplish this.
    !  If MITER = 3 it updates the coefficient H*RL1 in the diagonal
    !  matrix, and then computes the solution.
    !  If MITER is 4 or 5, it calls DGBSL.
    !  Communication with DVSOL uses the following variables:
    !  WM    = Real work space containing the inverse diagonal matrix if
    !          MITER = 3 and the LU decomposition of the matrix otherwise.
    !          Storage of matrix elements starts at WM(3).
    !          WM also contains the following matrix-related data:
    !          WM(1) = SQRT(UROUND) (not used here),
    !          WM(2) = HRL1, the previous value of H*RL1, used if MITER = 3.
    !  IWM   = Integer work space containing pivot information, starting at
    !          IWM(31), if MITER is 1, 2, 4, or 5.  IWM also contains band
    !          parameters ML = IWM(1) and MU = IWM(2) if MITER is 4 or 5.
    !  X     = The right-hand side vector on input, and the solution vector
    !          on output, of length N.
    !  IERSL = Output flag.  IERSL = 0 if no trouble occurred.
    !          IERSL = 1 if a singular matrix arose with MITER = 3.
    ! -----------------------------------------------------------------------
    ! 

    real(dp_t) :: WM(:), X(:)
    integer    :: IWM(:), IERSL
    type(dvode_t) :: dvode_state

    integer    :: I, MEBAND, ML, MU
    real(dp_t) :: DI, HRL1, PHRL1, R

    IERSL = 0
    GO TO (100, 100, 300, 400, 400), dvode_state % MITER
100 CALL DGESL (WM(3), dvode_state % N, dvode_state % N, IWM(31), X, 0)
    RETURN

300 PHRL1 = WM(2)
    HRL1 = dvode_state % H*dvode_state % RL1
    WM(2) = HRL1
    IF (HRL1 .EQ. PHRL1) GO TO 330
    R = HRL1/PHRL1
    do I = 1,dvode_state % N
       DI = ONE - R*(ONE - ONE/WM(I+2))
       IF (ABS(DI) .EQ. ZERO) GO TO 390
       WM(I+2) = ONE/DI
    end do

330 do I = 1,dvode_state % N
       X(I) = WM(I+2)*X(I)
    end do
    RETURN
390 IERSL = 1
    RETURN

400 ML = IWM(1)
    MU = IWM(2)
    MEBAND = 2*ML + MU + 1
    CALL DGBSL (WM(3), MEBAND, dvode_state % N, ML, MU, IWM(31), X, 0)
    RETURN
  end subroutine dvsol
  
  subroutine dacopy(NROW, NCOL, A, NROWA, B, NROWB)
    ! -----------------------------------------------------------------------
    !  Call sequence input -- NROW, NCOL, A, NROWA, NROWB
    !  Call sequence output -- B
    !  COMMON block variables accessed -- None
    ! 
    !  Subroutines called by DACOPY: DCOPY
    !  Function routines called by DACOPY: None
    ! -----------------------------------------------------------------------
    !  This routine copies one rectangular array, A, to another, B,
    !  where A and B may have different row dimensions, NROWA and NROWB.
    !  The data copied consists of NROW rows and NCOL columns.
    ! -----------------------------------------------------------------------
    integer    :: NROW, NCOL, NROWA, NROWB
    real(dp_t) :: A(NROWA,NCOL), B(NROWB,NCOL)
    integer    :: IC

    do IC = 1,NCOL
       CALL DCOPY (NROW, A(1,IC), 1, B(1,IC), 1)
    end do
    RETURN
  end subroutine dacopy

  subroutine dvjac(Y, YH, LDYH, EWT, FTEM, SAVF, WM, IWM, F, JAC, &
       IERPJ, RPAR, IPAR, dvode_state)
    ! -----------------------------------------------------------------------
    !  Call sequence input -- Y, YH, LDYH, EWT, FTEM, SAVF, WM, IWM,
    !                         F, JAC, RPAR, IPAR
    !  Call sequence output -- WM, IWM, IERPJ
    !  COMMON block variables accessed:
    !      /DVOD01/  CCMXJ, DRC, H, RL1, TN, UROUND, ICF, JCUR, LOCJS,
    !                MITER, MSBJ, N, NSLJ
    !      /DVOD02/  NFE, NST, NJE, NLU
    ! 
    !  Subroutines called by DVJAC: F, JAC, DACOPY, DCOPY, DGBFA, DGEFA,
    !                               DSCAL
    !  Function routines called by DVJAC: DVNORM
    ! -----------------------------------------------------------------------
    !  DVJAC is called by DVNLSD to compute and process the matrix
    !  P = I - h*rl1*J , where J is an approximation to the Jacobian.
    !  Here J is computed by the user-supplied routine JAC if
    !  MITER = 1 or 4, or by finite differencing if MITER = 2, 3, or 5.
    !  If MITER = 3, a diagonal approximation to J is used.
    !  If JSV = -1, J is computed from scratch in all cases.
    !  If JSV = 1 and MITER = 1, 2, 4, or 5, and if the saved value of J is
    !  considered acceptable, then P is constructed from the saved J.
    !  J is stored in wm and replaced by P.  If MITER .ne. 3, P is then
    !  subjected to LU decomposition in preparation for later solution
    !  of linear systems with P as coefficient matrix. This is done
    !  by DGEFA if MITER = 1 or 2, and by DGBFA if MITER = 4 or 5.
    ! 
    !  Communication with DVJAC is done with the following variables.  (For
    !  more details, please see the comments in the driver subroutine.)
    !  Y          = Vector containing predicted values on entry.
    !  YH         = The Nordsieck array, an LDYH by LMAX array, input.
    !  LDYH       = A constant .ge. N, the first dimension of YH, input.
    !  EWT        = An error weight vector of length N.
    !  SAVF       = Array containing f evaluated at predicted y, input.
    !  WM         = Real work space for matrices.  In the output, it containS
    !               the inverse diagonal matrix if MITER = 3 and the LU
    !               decomposition of P if MITER is 1, 2 , 4, or 5.
    !               Storage of matrix elements starts at WM(3).
    !               Storage of the saved Jacobian starts at WM(LOCJS).
    !               WM also contains the following matrix-related data:
    !               WM(1) = SQRT(UROUND), used in numerical Jacobian step.
    !               WM(2) = H*RL1, saved for later use if MITER = 3.
    !  IWM        = Integer work space containing pivot information,
    !               starting at IWM(31), if MITER is 1, 2, 4, or 5.
    !               IWM also contains band parameters ML = IWM(1) and
    !               MU = IWM(2) if MITER is 4 or 5.
    !  F          = Dummy name for the user supplied subroutine for f.
    !  JAC        = Dummy name for the user supplied Jacobian subroutine.
    !  RPAR, IPAR = Dummy names for user's real and integer work arrays.
    !  RL1        = 1/EL(2) (input).
    !  IERPJ      = Output error flag,  = 0 if no trouble, 1 if the P
    !               matrix is found to be singular.
    !  JCUR       = Output flag to indicate whether the Jacobian matrix
    !               (or approximation) is now current.
    !               JCUR = 0 means J is not current.
    !               JCUR = 1 means J is current.
    ! -----------------------------------------------------------------------
    ! 

    EXTERNAL F, JAC
    type(dvode_t) :: dvode_state
    real(dp_t) :: Y(:), YH(LDYH, dvode_state % LMAX), EWT(:), FTEM(:)
    real(dp_t) :: SAVF(:), WM(:), RPAR(:)
    integer    :: LDYH, IWM(:), IERPJ, IPAR(:)


    real(dp_t) :: CON, DI, FAC, HRL1, R, R0, SRUR, YI, YJ, YJJ
    integer    :: I, I1, I2, IER, II, J, J1, JJ, JOK, LENP, MBA, MBAND
    integer    :: MEB1, MEBAND, ML, ML3, MU, NP1

    ! Parameter declarations
    real(dp_t), parameter :: PT1 = 0.1D0

    IERPJ = 0
    HRL1 = dvode_state % H*dvode_state % RL1
    ! See whether J should be evaluated (JOK = -1) or not (JOK = 1). -------
    JOK = dvode_state % JSV
    IF (dvode_state % JSV .EQ. 1) THEN
       IF (dvode_state % NST .EQ. 0 .OR. dvode_state % NST .GT. dvode_state % NSLJ+dvode_state % MSBJ) JOK = -1
       IF (dvode_state % ICF .EQ. 1 .AND. dvode_state % DRC .LT. dvode_state % CCMXJ) JOK = -1
       IF (dvode_state % ICF .EQ. 2) JOK = -1
    ENDIF
    ! End of setting JOK. --------------------------------------------------

    IF (JOK .EQ. -1 .AND. dvode_state % MITER .EQ. 1) THEN
       ! If JOK = -1 and MITER = 1, call JAC to evaluate Jacobian. ------------
       dvode_state % NJE = dvode_state % NJE + 1
       dvode_state % NSLJ = dvode_state % NST
       dvode_state % JCUR = 1
       LENP = dvode_state % N * dvode_state % N
       do I = 1,LENP
          WM(I+2) = ZERO
       end do
       CALL JAC (dvode_state % N, dvode_state % TN, Y, 0, 0, WM(3), dvode_state % N, RPAR, IPAR)
       IF (dvode_state % JSV .EQ. 1) CALL DCOPY (LENP, WM(3), 1, WM(dvode_state % LOCJS), 1)
    ENDIF

    IF (JOK .EQ. -1 .AND. dvode_state % MITER .EQ. 2) THEN
       ! If MITER = 2, make N calls to F to approximate the Jacobian. ---------
       dvode_state % NJE = dvode_state % NJE + 1
       dvode_state % NSLJ = dvode_state % NST
       dvode_state % JCUR = 1
       FAC = DVNORM (dvode_state % N, SAVF, EWT)
       R0 = THOU*ABS(dvode_state % H) * dvode_state % UROUND * REAL(dvode_state % N)*FAC
       IF (R0 .EQ. ZERO) R0 = ONE
       SRUR = WM(1)
       J1 = 2
       do J = 1,dvode_state % N
          YJ = Y(J)
          R = MAX(SRUR*ABS(YJ),R0/EWT(J))
          Y(J) = Y(J) + R
          FAC = ONE/R
          CALL F (dvode_state % N, dvode_state % TN, Y, FTEM, RPAR, IPAR)
          do I = 1,dvode_state % N
             WM(I+J1) = (FTEM(I) - SAVF(I))*FAC
          end do
          Y(J) = YJ
          J1 = J1 + dvode_state % N
       end do
       dvode_state % NFE = dvode_state % NFE + dvode_state % N
       LENP = dvode_state % N * dvode_state % N
       IF (dvode_state % JSV .EQ. 1) CALL DCOPY (LENP, WM(3), 1, WM(dvode_state % LOCJS), 1)
    ENDIF

    IF (JOK .EQ. 1 .AND. (dvode_state % MITER .EQ. 1 .OR. dvode_state % MITER .EQ. 2)) THEN
       dvode_state % JCUR = 0
       LENP = dvode_state % N * dvode_state % N
       CALL DCOPY (LENP, WM(dvode_state % LOCJS), 1, WM(3), 1)
    ENDIF

    IF (dvode_state % MITER .EQ. 1 .OR. dvode_state % MITER .EQ. 2) THEN
       ! Multiply Jacobian by scalar, add identity, and do LU decomposition. --
       CON = -HRL1
       CALL DSCAL (LENP, CON, WM(3), 1)
       J = 3
       NP1 = dvode_state % N + 1
       do I = 1,dvode_state % N
          WM(J) = WM(J) + ONE
          J = J + NP1
       end do
       dvode_state % NLU = dvode_state % NLU + 1
       CALL DGEFA (WM(3), dvode_state % N, dvode_state % N, IWM(31), IER)
       IF (IER .NE. 0) IERPJ = 1
       RETURN
    ENDIF
    ! End of code block for MITER = 1 or 2. --------------------------------

    IF (dvode_state % MITER .EQ. 3) THEN
       ! If MITER = 3, construct a diagonal approximation to J and P. ---------
       dvode_state % NJE = dvode_state % NJE + 1
       dvode_state % JCUR = 1
       WM(2) = HRL1
       R = dvode_state % RL1*PT1
       do I = 1,dvode_state % N
          Y(I) = Y(I) + R*(dvode_state % H*SAVF(I) - YH(I,2))
       end do
       CALL F (dvode_state % N, dvode_state % TN, Y, WM(3), RPAR, IPAR)
       dvode_state % NFE = dvode_state % NFE + 1
       do I = 1,dvode_state % N
          R0 = dvode_state % H*SAVF(I) - YH(I,2)
          DI = PT1*R0 - dvode_state % H*(WM(I+2) - SAVF(I))
          WM(I+2) = ONE
          IF (ABS(R0) .LT. dvode_state % UROUND/EWT(I)) cycle
          IF (ABS(DI) .EQ. ZERO) GO TO 330
          WM(I+2) = PT1*R0/DI
       end do
       RETURN
330    IERPJ = 1
       RETURN
    ENDIF
    ! End of code block for MITER = 3. -------------------------------------

    ! Set constants for MITER = 4 or 5. ------------------------------------
    ML = IWM(1)
    MU = IWM(2)
    ML3 = ML + 3
    MBAND = ML + MU + 1
    MEBAND = MBAND + ML
    LENP = MEBAND * dvode_state % N

    if (JOK .EQ. -1 .AND. dvode_state % MITER .EQ. 4) then
       ! If JOK = -1 and MITER = 4, call JAC to evaluate Jacobian. ------------
       dvode_state % NJE = dvode_state % NJE + 1
       dvode_state % NSLJ = dvode_state % NST
       dvode_state % JCUR = 1
       do I = 1,LENP
          WM(I+2) = ZERO
       end do
       CALL JAC (dvode_state % N, dvode_state % TN, Y, ML, MU, WM(ML3), MEBAND, RPAR, IPAR)
       if (dvode_state % JSV .EQ. 1) then
          CALL DACOPY(MBAND, dvode_state % N, &
               WM(ML3:ML3 + MEBAND * dvode_state % N - 1), &
               MEBAND, &
               WM(dvode_state % LOCJS:dvode_state % LOCJS + MBAND * dvode_state % N - 1), &
               MBAND)
       end if

    else if (JOK .EQ. -1 .AND. dvode_state % MITER .EQ. 5) then
       ! If MITER = 5, make ML+MU+1 calls to F to approximate the Jacobian. ---
       dvode_state % NJE = dvode_state % NJE + 1
       dvode_state % NSLJ = dvode_state % NST
       dvode_state % JCUR = 1
       MBA = MIN(MBAND,dvode_state % N)
       MEB1 = MEBAND - 1
       SRUR = WM(1)
       FAC = DVNORM (dvode_state % N, SAVF, EWT)
       R0 = THOU*ABS(dvode_state % H) * dvode_state % UROUND * REAL(dvode_state % N)*FAC
       IF (R0 .EQ. ZERO) R0 = ONE
       do J = 1,MBA
          do I = J,dvode_state % N,MBAND
             YI = Y(I)
             R = MAX(SRUR*ABS(YI),R0/EWT(I))
             Y(I) = Y(I) + R
          end do
          CALL F (dvode_state % N, dvode_state % TN, Y, FTEM, RPAR, IPAR)
          do JJ = J,dvode_state % N,MBAND
             Y(JJ) = YH(JJ,1)
             YJJ = Y(JJ)
             R = MAX(SRUR*ABS(YJJ),R0/EWT(JJ))
             FAC = ONE/R
             I1 = MAX(JJ-MU,1)
             I2 = MIN(JJ+ML,dvode_state % N)
             II = JJ*MEB1 - ML + 2
             do I = I1,I2
                WM(II+I) = (FTEM(I) - SAVF(I))*FAC
             end do
          end do
       end do
       dvode_state % NFE = dvode_state % NFE + MBA
       if (dvode_state % JSV .EQ. 1) then
          CALL DACOPY(MBAND, dvode_state % N, &
               WM(ML3:ML3 + MEBAND * dvode_state % N - 1), &
               MEBAND, &
               WM(dvode_state % LOCJS:dvode_state % LOCJS + MBAND * dvode_state % N - 1), &
               MBAND)
       end if
    end if

    IF (JOK .EQ. 1) THEN
       dvode_state % JCUR = 0
       CALL DACOPY(MBAND, dvode_state % N, &
            WM(dvode_state % LOCJS:dvode_state % LOCJS + MBAND * dvode_state % N - 1), &
            MBAND, &
            WM(ML3:ML3 + MEBAND * dvode_state % N - 1), &
            MEBAND)
    ENDIF

    ! Multiply Jacobian by scalar, add identity, and do LU decomposition.
    CON = -HRL1
    CALL DSCAL (LENP, CON, WM(3), 1 )
    II = MBAND + 2
    do I = 1,dvode_state % N
       WM(II) = WM(II) + ONE
       II = II + MEBAND
    end do
    dvode_state % NLU = dvode_state % NLU + 1
    CALL DGBFA (WM(3), MEBAND, dvode_state % N, ML, MU, IWM(31), IER)
    if (IER .NE. 0) then
       IERPJ = 1
    end if
    RETURN
    ! End of code block for MITER = 4 or 5. --------------------------------
  end subroutine dvjac
  
  subroutine dvnlsd(Y, YH, LDYH, VSAV, SAVF, EWT, ACOR, IWM, WM, &
       F, JAC, PDUM, NFLAG, RPAR, IPAR, dvode_state)
    ! -----------------------------------------------------------------------
    !  Call sequence input -- Y, YH, LDYH, SAVF, EWT, ACOR, IWM, WM,
    !                         F, JAC, NFLAG, RPAR, IPAR
    !  Call sequence output -- YH, ACOR, WM, IWM, NFLAG
    !  COMMON block variables accessed:
    !      /DVOD01/ ACNRM, CRATE, DRC, H, RC, RL1, TQ(5), TN, ICF,
    !                 JCUR, METH, MITER, N, NSLP
    !      /DVOD02/ HU, NCFN, NETF, NFE, NJE, NLU, NNI, NQU, NST
    ! 
    !  Subroutines called by DVNLSD: F, DAXPY, DCOPY, DSCAL, DVJAC, DVSOL
    !  Function routines called by DVNLSD: DVNORM
    ! -----------------------------------------------------------------------
    !  Subroutine DVNLSD is a nonlinear system solver, which uses functional
    !  iteration or a chord (modified Newton) method.  For the chord method
    !  direct linear algebraic system solvers are used.  Subroutine DVNLSD
    !  then handles the corrector phase of this integration package.
    ! 
    !  Communication with DVNLSD is done with the following variables. (For
    !  more details, please see the comments in the driver subroutine.)
    ! 
    !  Y          = The dependent variable, a vector of length N, input.
    !  YH         = The Nordsieck (Taylor) array, LDYH by LMAX, input
    !               and output.  On input, it contains predicted values.
    !  LDYH       = A constant .ge. N, the first dimension of YH, input.
    !  VSAV       = Unused work array.
    !  SAVF       = A work array of length N.
    !  EWT        = An error weight vector of length N, input.
    !  ACOR       = A work array of length N, used for the accumulated
    !               corrections to the predicted y vector.
    !  WM,IWM     = Real and integer work arrays associated with matrix
    !               operations in chord iteration (MITER .ne. 0).
    !  F          = Dummy name for user supplied routine for f.
    !  JAC        = Dummy name for user supplied Jacobian routine.
    !  PDUM       = Unused dummy subroutine name.  Included for uniformity
    !               over collection of integrators.
    !  NFLAG      = Input/output flag, with values and meanings as follows:
    !               INPUT
    !                   0 first call for this time step.
    !                  -1 convergence failure in previous call to DVNLSD.
    !                  -2 error test failure in DVSTEP.
    !               OUTPUT
    !                   0 successful completion of nonlinear solver.
    !                  -1 convergence failure or singular matrix.
    !                  -2 unrecoverable error in matrix preprocessing
    !                     (cannot occur here).
    !                  -3 unrecoverable error in solution (cannot occur
    !                     here).
    !  RPAR, IPAR = Dummy names for user's real and integer work arrays.
    ! 
    !  IPUP       = Own variable flag with values and meanings as follows:
    !               0,            do not update the Newton matrix.
    !               MITER .ne. 0, update Newton matrix, because it is the
    !                             initial step, order was changed, the error
    !                             test failed, or an update is indicated by
    !                             the scalar RC or step counter NST.
    ! 
    !  For more details, see comments in driver subroutine.
    ! -----------------------------------------------------------------------
    !
    EXTERNAL F, JAC, PDUM
    type(dvode_t) :: dvode_state
    real(dp_t) :: Y(dvode_state % N), YH(LDYH, dvode_state % LMAX), VSAV(:), SAVF(:)
    real(dp_t) :: EWT(:), ACOR(:), WM(:), RPAR(:)
    integer    :: LDYH, IWM(:), NFLAG, IPAR(:)
    
    real(dp_t) :: CSCALE, DCON, DEL, DELP
    integer    :: I, IERPJ, IERSL, M

    ! Parameter declarations
    real(dp_t), parameter :: CCMAX = 0.3D0
    real(dp_t), parameter :: CRDOWN = 0.3D0
    real(dp_t), parameter :: RDIV  = 2.0D0
    integer, parameter :: MAXCOR = 3
    integer, parameter :: MSBP = 20

    ! -----------------------------------------------------------------------
    !  On the first step, on a change of method order, or after a
    !  nonlinear convergence failure with NFLAG = -2, set IPUP = MITER
    !  to force a Jacobian update when MITER .ne. 0.
    ! -----------------------------------------------------------------------
    IF (dvode_state % JSTART .EQ. 0) dvode_state % NSLP = 0
    IF (NFLAG .EQ. 0) dvode_state % ICF = 0
    IF (NFLAG .EQ. -2) dvode_state % IPUP = dvode_state % MITER
    IF ( (dvode_state % JSTART .EQ. 0) .OR. (dvode_state % JSTART .EQ. -1) ) dvode_state % IPUP = dvode_state % MITER
    ! If this is functional iteration, set CRATE .eq. 1 and drop to 220
    IF (dvode_state % MITER .EQ. 0) THEN
       dvode_state % CRATE = ONE
       GO TO 220
    ENDIF
    ! -----------------------------------------------------------------------
    !  RC is the ratio of new to old values of the coefficient H/EL(2)=h/l1.
    !  When RC differs from 1 by more than CCMAX, IPUP is set to MITER
    !  to force DVJAC to be called, if a Jacobian is involved.
    !  In any case, DVJAC is called at least every MSBP steps.
    ! -----------------------------------------------------------------------
    dvode_state % DRC = ABS(dvode_state % RC-ONE)
    IF (dvode_state % DRC .GT. CCMAX .OR. dvode_state % NST .GE. dvode_state % NSLP+MSBP) dvode_state % IPUP = dvode_state % MITER
    ! -----------------------------------------------------------------------
    !  Up to MAXCOR corrector iterations are taken.  A convergence test is
    !  made on the r.m.s. norm of each correction, weighted by the error
    !  weight vector EWT.  The sum of the corrections is accumulated in the
    !  vector ACOR(i).  The YH array is not altered in the corrector loop.
    ! -----------------------------------------------------------------------
220 M = 0
    DELP = ZERO
    !    CALL DCOPY (dvode_state % N, YH(1,1), 1, Y, 1)
    CALL DCOPY (dvode_state % N, YH(1:LDYH, 1:dvode_state % LMAX), &
         1, Y(1:dvode_state % N), 1)
    CALL F (dvode_state % N, dvode_state % TN, Y, SAVF, RPAR, IPAR)
    dvode_state % NFE = dvode_state % NFE + 1
    IF (dvode_state % IPUP .LE. 0) GO TO 250
    ! -----------------------------------------------------------------------
    !  If indicated, the matrix P = I - h*rl1*J is reevaluated and
    !  preprocessed before starting the corrector iteration.  IPUP is set
    !  to 0 as an indicator that this has been done.
    ! -----------------------------------------------------------------------
    CALL DVJAC (Y, YH, LDYH, EWT, ACOR, SAVF, WM, IWM, F, JAC, IERPJ, &
         RPAR, IPAR, dvode_state)
    dvode_state % IPUP = 0
    dvode_state % RC = ONE
    dvode_state % DRC = ZERO
    dvode_state % CRATE = ONE
    dvode_state % NSLP = dvode_state % NST
    ! If matrix is singular, take error return to force cut in step size. --
    IF (IERPJ .NE. 0) GO TO 430
250 do I = 1,dvode_state % N
       ACOR(I) = ZERO
    end do
    ! This is a looping point for the corrector iteration. -----------------
270 IF (dvode_state % MITER .NE. 0) GO TO 350
    ! -----------------------------------------------------------------------
    !  In the case of functional iteration, update Y directly from
    !  the result of the last function evaluation.
    ! -----------------------------------------------------------------------
    do I = 1,dvode_state % N
       SAVF(I) = dvode_state % RL1*(dvode_state % H*SAVF(I) - YH(I,2))
    end do
    do I = 1,dvode_state % N
       Y(I) = SAVF(I) - ACOR(I)
    end do
    DEL = DVNORM (dvode_state % N, Y, EWT)
    do I = 1,dvode_state % N
       Y(I) = YH(I,1) + SAVF(I)
    end do
    CALL DCOPY(dvode_state % N, SAVF, 1, ACOR, 1)
    GO TO 400
    ! -----------------------------------------------------------------------
    !  In the case of the chord method, compute the corrector error,
    !  and solve the linear system with that as right-hand side and
    !  P as coefficient matrix.  The correction is scaled by the factor
    !  2/(1+RC) to account for changes in h*rl1 since the last DVJAC call.
    ! -----------------------------------------------------------------------
350 do I = 1,dvode_state % N
       Y(I) = (dvode_state % RL1*dvode_state % H)*SAVF(I) - (dvode_state % RL1*YH(I,2) + ACOR(I))
    end do
    CALL DVSOL (WM, IWM, Y, IERSL, dvode_state)
    dvode_state % NNI = dvode_state % NNI + 1
    IF (IERSL .GT. 0) GO TO 410
    IF (dvode_state % METH .EQ. 2 .AND. dvode_state % RC .NE. ONE) THEN
       CSCALE = TWO/(ONE + dvode_state % RC)
       CALL DSCAL (dvode_state % N, CSCALE, Y, 1)
    ENDIF
    DEL = DVNORM (dvode_state % N, Y, EWT)
    CALL DAXPY (dvode_state % N, ONE, Y, 1, ACOR, 1)
    do I = 1,dvode_state % N
       Y(I) = YH(I,1) + ACOR(I)
    end do
    ! -----------------------------------------------------------------------
    !  Test for convergence.  If M .gt. 0, an estimate of the convergence
    !  rate constant is stored in CRATE, and this is used in the test.
    ! -----------------------------------------------------------------------
400 IF (M .NE. 0) dvode_state % CRATE = MAX(CRDOWN*dvode_state % CRATE,DEL/DELP)
    DCON = DEL*MIN(ONE,dvode_state % CRATE)/dvode_state % TQ(4)
    IF (DCON .LE. ONE) GO TO 450
    M = M + 1
    IF (M .EQ. MAXCOR) GO TO 410
    IF (M .GE. 2 .AND. DEL .GT. RDIV*DELP) GO TO 410
    DELP = DEL
    CALL F (dvode_state % N, dvode_state % TN, Y, SAVF, RPAR, IPAR)
    dvode_state % NFE = dvode_state % NFE + 1
    GO TO 270
    
410 IF (dvode_state % MITER .EQ. 0 .OR. dvode_state % JCUR .EQ. 1) GO TO 430
    dvode_state % ICF = 1
    dvode_state % IPUP = dvode_state % MITER
    GO TO 220
    
430 CONTINUE
    NFLAG = -1
    dvode_state % ICF = 2
    dvode_state % IPUP = dvode_state % MITER
    RETURN

    ! Return for successful step. ------------------------------------------
450 NFLAG = 0
    dvode_state % JCUR = 0
    dvode_state % ICF = 0
    IF (M .EQ. 0) dvode_state % ACNRM = DEL
    IF (M .GT. 0) dvode_state % ACNRM = DVNORM (dvode_state % N, ACOR, EWT)
    RETURN
  end subroutine dvnlsd
  
  subroutine dvjust(YH, LDYH, IORD, dvode_state)
    ! -----------------------------------------------------------------------
    !  Call sequence input -- YH, LDYH, IORD
    !  Call sequence output -- YH
    !  COMMON block input -- NQ, METH, LMAX, HSCAL, TAU(13), N
    !  COMMON block variables accessed:
    !      /DVOD01/ -- HSCAL, TAU(13), LMAX, METH, N, NQ,
    ! 
    !  Subroutines called by DVJUST: DAXPY
    !  Function routines called by DVJUST: None
    ! -----------------------------------------------------------------------
    !  This subroutine adjusts the YH array on reduction of order,
    !  and also when the order is increased for the stiff option (METH = 2).
    !  Communication with DVJUST uses the following:
    !  IORD  = An integer flag used when METH = 2 to indicate an order
    !          increase (IORD = +1) or an order decrease (IORD = -1).
    !  HSCAL = Step size H used in scaling of Nordsieck array YH.
    !          (If IORD = +1, DVJUST assumes that HSCAL = TAU(1).)
    !  See References 1 and 2 for details.
    ! -----------------------------------------------------------------------
    !
    type(dvode_t) :: dvode_state
    real(dp_t) :: YH(LDYH, dvode_state % LMAX)
    integer    :: LDYH, IORD

    real(dp_t) :: ALPH0, ALPH1, HSUM, PROD, T1, XI,XIOLD
    integer    :: I, IBACK, J, JP1, LP1, NQM1, NQM2, NQP1

    IF ((dvode_state % NQ .EQ. 2) .AND. (IORD .NE. 1)) RETURN
    NQM1 = dvode_state % NQ - 1
    NQM2 = dvode_state % NQ - 2
    GO TO (100, 200), dvode_state % METH
    ! -----------------------------------------------------------------------
    !  Nonstiff option...
    !  Check to see if the order is being increased or decreased.
    ! -----------------------------------------------------------------------
100 CONTINUE
    IF (IORD .EQ. 1) GO TO 180
    ! Order decrease. ------------------------------------------------------
    do J = 1, dvode_state % LMAX
       dvode_state % EL(J) = ZERO
    end do
    dvode_state % EL(2) = ONE
    HSUM = ZERO
    do J = 1, NQM2
       ! Construct coefficients of x*(x+xi(1))*...*(x+xi(j)). -----------------
       HSUM = HSUM + dvode_state % TAU(J)
       XI = HSUM/dvode_state % HSCAL
       JP1 = J + 1
       do IBACK = 1, JP1
          I = (J + 3) - IBACK
          dvode_state % EL(I) = dvode_state % EL(I)*XI + dvode_state % EL(I-1)
       end do
    end do
    ! Construct coefficients of integrated polynomial. ---------------------
    do J = 2, NQM1
       dvode_state % EL(J+1) = REAL(dvode_state % NQ) * dvode_state % EL(J)/REAL(J)
    end do
    ! Subtract correction terms from YH array. -----------------------------
    do J = 3, dvode_state % NQ
       do I = 1, dvode_state % N
          YH(I,J) = YH(I,J) - YH(I,dvode_state % L) * dvode_state % EL(J)
       end do
    end do
    RETURN
    ! Order increase. ------------------------------------------------------
    ! Zero out next column in YH array. ------------------------------------
180 CONTINUE
    LP1 = dvode_state % L + 1
    do I = 1, dvode_state % N
       YH(I,LP1) = ZERO
    end do
    RETURN
    ! -----------------------------------------------------------------------
    !  Stiff option...
    !  Check to see if the order is being increased or decreased.
    ! -----------------------------------------------------------------------
200 CONTINUE
    IF (IORD .EQ. 1) GO TO 300
    ! Order decrease. ------------------------------------------------------
    do J = 1, dvode_state % LMAX
       dvode_state % EL(J) = ZERO
    end do
    dvode_state % EL(3) = ONE
    HSUM = ZERO
    do J = 1,NQM2
       ! Construct coefficients of x*x*(x+xi(1))*...*(x+xi(j)). ---------------
       HSUM = HSUM + dvode_state % TAU(J)
       XI = HSUM/dvode_state % HSCAL
       JP1 = J + 1
       do IBACK = 1, JP1
          I = (J + 4) - IBACK
          dvode_state % EL(I) = dvode_state % EL(I) * XI + dvode_state % EL(I-1)
       end do
    end do
    ! Subtract correction terms from YH array. -----------------------------
    do J = 3,dvode_state % NQ
       do I = 1, dvode_state % N
          YH(I,J) = YH(I,J) - YH(I,dvode_state % L) * dvode_state % EL(J)
       end do
    end do
    RETURN
    ! Order increase. ------------------------------------------------------
300 do J = 1, dvode_state % LMAX
       dvode_state % EL(J) = ZERO
    end do
    dvode_state % EL(3) = ONE
    ALPH0 = -ONE
    ALPH1 = ONE
    PROD = ONE
    XIOLD = ONE
    HSUM = dvode_state % HSCAL
    IF (dvode_state % NQ .EQ. 1) GO TO 340
    do J = 1, NQM1
       ! Construct coefficients of x*x*(x+xi(1))*...*(x+xi(j)). ---------------
       JP1 = J + 1
       HSUM = HSUM + dvode_state % TAU(JP1)
       XI = HSUM/dvode_state % HSCAL
       PROD = PROD*XI
       ALPH0 = ALPH0 - ONE/REAL(JP1)
       ALPH1 = ALPH1 + ONE/XI
       do IBACK = 1, JP1
          I = (J + 4) - IBACK
          dvode_state % EL(I) = dvode_state % EL(I) * XIOLD + dvode_state % EL(I-1)
       end do
       XIOLD = XI
    end do
340 CONTINUE
    T1 = (-ALPH0 - ALPH1)/PROD
    ! Load column L + 1 in YH array. ---------------------------------------
    LP1 = dvode_state % L + 1
    do I = 1, dvode_state % N
       YH(I,LP1) = T1*YH(I,dvode_state % LMAX)
    end do
    ! Add correction terms to YH array. ------------------------------------
    NQP1 = dvode_state % NQ + 1
    do J = 3, NQP1
       CALL DAXPY (dvode_state % N, dvode_state % EL(J), YH(1,LP1), 1, YH(1,J), 1 )
    end do
    RETURN
  end subroutine dvjust

  subroutine dvset(dvode_state)
    ! -----------------------------------------------------------------------
    !  Call sequence communication: None
    !  COMMON block variables accessed:
    !      /DVOD01/ -- EL(13), H, TAU(13), TQ(5), L(= NQ + 1),
    !                  METH, NQ, NQWAIT
    ! 
    !  Subroutines called by DVSET: None
    !  Function routines called by DVSET: None
    ! -----------------------------------------------------------------------
    !  DVSET is called by DVSTEP and sets coefficients for use there.
    ! 
    !  For each order NQ, the coefficients in EL are calculated by use of
    !   the generating polynomial lambda(x), with coefficients EL(i).
    !       lambda(x) = EL(1) + EL(2)*x + ... + EL(NQ+1)*(x**NQ).
    !  For the backward differentiation formulas,
    !                                      NQ-1
    !       lambda(x) = (1 + x/xi*(NQ)) * product (1 + x/xi(i) ) .
    !                                      i = 1
    !  For the Adams formulas,
    !                               NQ-1
    !       (d/dx) lambda(x) = c * product (1 + x/xi(i) ) ,
    !                               i = 1
    !       lambda(-1) = 0,    lambda(0) = 1,
    !  where c is a normalization constant.
    !  In both cases, xi(i) is defined by
    !       H*xi(i) = t sub n  -  t sub (n-i)
    !               = H + TAU(1) + TAU(2) + ... TAU(i-1).
    ! 
    ! 
    !  In addition to variables described previously, communication
    !  with DVSET uses the following:
    !    TAU    = A vector of length 13 containing the past NQ values
    !             of H.
    !    EL     = A vector of length 13 in which vset stores the
    !             coefficients for the corrector formula.
    !    TQ     = A vector of length 5 in which vset stores constants
    !             used for the convergence test, the error test, and the
    !             selection of H at a new order.
    !    METH   = The basic method indicator.
    !    NQ     = The current order.
    !    L      = NQ + 1, the length of the vector stored in EL, and
    !             the number of columns of the YH array being used.
    !    NQWAIT = A counter controlling the frequency of order changes.
    !             An order change is about to be considered if NQWAIT = 1.
    ! -----------------------------------------------------------------------
    !

    type(dvode_t), intent(inout) :: dvode_state
    real(dp_t) :: EM(13)
    real(dp_t) :: AHATN0, ALPH0, CNQM1, CSUM, ELP
    real(dp_t) :: EM0, FLOTI, FLOTL, FLOTNQ, HSUM, RXI, RXIS, S
    real(dp_t) :: T1, T2, T3, T4, T5, T6, XI
    integer    :: I, IBACK, J, JP1, NQM1, NQM2

    ! Parameter declaration
    real(dp_t), parameter :: CORTES = 0.1D0

    FLOTL = REAL(dvode_state % L)
    NQM1 = dvode_state % NQ - 1
    NQM2 = dvode_state % NQ - 2
    GO TO (100, 200), dvode_state % METH

    !  Set coefficients for Adams methods. ----------------------------------
100 IF (dvode_state % NQ .NE. 1) GO TO 110
    dvode_state % EL(1) = ONE
    dvode_state % EL(2) = ONE
    dvode_state % TQ(1) = ONE
    dvode_state % TQ(2) = TWO
    dvode_state % TQ(3) = SIX*dvode_state % TQ(2)
    dvode_state % TQ(5) = ONE
    GO TO 300
110 HSUM = dvode_state % H
    EM(1) = ONE
    FLOTNQ = FLOTL - ONE
    do I = 2, dvode_state % L
       EM(I) = ZERO
    end do
    do J = 1, NQM1
       IF ((J .NE. NQM1) .OR. (dvode_state % NQWAIT .NE. 1)) GO TO 130
       S = ONE
       CSUM = ZERO
       do I = 1, NQM1
          CSUM = CSUM + S*EM(I)/REAL(I+1)
          S = -S
       end do
       dvode_state % TQ(1) = EM(NQM1)/(FLOTNQ*CSUM)
130    RXI = dvode_state % H/HSUM
       do IBACK = 1, J
          I = (J + 2) - IBACK
          EM(I) = EM(I) + EM(I-1)*RXI
       end do
       HSUM = HSUM + dvode_state % TAU(J)
    end do
    ! Compute integral from -1 to 0 of polynomial and of x times it. -------
    S = ONE
    EM0 = ZERO
    CSUM = ZERO
    do I = 1, dvode_state % NQ
       FLOTI = REAL(I)
       EM0 = EM0 + S*EM(I)/FLOTI
       CSUM = CSUM + S*EM(I)/(FLOTI+ONE)
       S = -S
    end do
    ! In EL, form coefficients of normalized integrated polynomial. --------
    S = ONE/EM0
    dvode_state % EL(1) = ONE
    do I = 1, dvode_state % NQ
       dvode_state % EL(I+1) = S*EM(I)/REAL(I)
    end do
    XI = HSUM/dvode_state % H
    dvode_state % TQ(2) = XI*EM0/CSUM
    dvode_state % TQ(5) = XI/dvode_state % EL(dvode_state % L)
    IF (dvode_state % NQWAIT .NE. 1) GO TO 300
    ! For higher order control constant, multiply polynomial by 1+x/xi(q). -
    RXI = ONE/XI
    do IBACK = 1, dvode_state % NQ
       I = (dvode_state % L + 1) - IBACK
       EM(I) = EM(I) + EM(I-1)*RXI
    end do
    ! Compute integral of polynomial. --------------------------------------
    S = ONE
    CSUM = ZERO
    do I = 1, dvode_state % L
       CSUM = CSUM + S*EM(I)/REAL(I+1)
       S = -S
    end do
    dvode_state % TQ(3) = FLOTL*EM0/CSUM
    GO TO 300

    ! Set coefficients for BDF methods. ------------------------------------
200 do I = 3, dvode_state % L
       dvode_state % EL(I) = ZERO
    end do
    dvode_state % EL(1) = ONE
    dvode_state % EL(2) = ONE
    ALPH0 = -ONE
    AHATN0 = -ONE
    HSUM = dvode_state % H
    RXI = ONE
    RXIS = ONE
    IF (dvode_state % NQ .EQ. 1) GO TO 240
    do J = 1, NQM2
       ! In EL, construct coefficients of (1+x/xi(1))*...*(1+x/xi(j+1)). ------
       HSUM = HSUM + dvode_state % TAU(J)
       RXI = dvode_state % H/HSUM
       JP1 = J + 1
       ALPH0 = ALPH0 - ONE/REAL(JP1)
       do IBACK = 1, JP1
          I = (J + 3) - IBACK
          dvode_state % EL(I) = dvode_state % EL(I) + dvode_state % EL(I-1)*RXI
       end do
    end do
    ALPH0 = ALPH0 - ONE/REAL(dvode_state % NQ)
    RXIS = -dvode_state % EL(2) - ALPH0
    HSUM = HSUM + dvode_state % TAU(NQM1)
    RXI = dvode_state % H/HSUM
    AHATN0 = -dvode_state % EL(2) - RXI
    do IBACK = 1, dvode_state % NQ
       I = (dvode_state % NQ + 2) - IBACK
       dvode_state % EL(I) = dvode_state % EL(I) + dvode_state % EL(I-1)*RXIS
    end do
240 T1 = ONE - AHATN0 + ALPH0
    T2 = ONE + REAL(dvode_state % NQ)*T1
    dvode_state % TQ(2) = ABS(ALPH0*T2/T1)
    dvode_state % TQ(5) = ABS(T2/(dvode_state % EL(dvode_state % L)*RXI/RXIS))
    IF (dvode_state % NQWAIT .NE. 1) GO TO 300
    CNQM1 = RXIS/dvode_state % EL(dvode_state % L)
    T3 = ALPH0 + ONE/REAL(dvode_state % NQ)
    T4 = AHATN0 + RXI
    ELP = T3/(ONE - T4 + T3)
    dvode_state % TQ(1) = ABS(ELP/CNQM1)
    HSUM = HSUM + dvode_state % TAU(dvode_state % NQ)
    RXI = dvode_state % H/HSUM
    T5 = ALPH0 - ONE/REAL(dvode_state % NQ+1)
    T6 = AHATN0 - RXI
    ELP = T2/(ONE - T6 + T5)
    dvode_state % TQ(3) = ABS(ELP*RXI*(FLOTL + ONE)*T5)
300 dvode_state % TQ(4) = CORTES*dvode_state % TQ(2)
    RETURN
  end subroutine dvset
  
  subroutine dvstep(Y, YH, LDYH, YH1, EWT, SAVF, VSAV, ACOR, &
       WM, IWM, F, JAC, PSOL, VNLS, RPAR, IPAR, dvode_state)
    ! -----------------------------------------------------------------------
    !  Call sequence input -- Y, YH, LDYH, YH1, EWT, SAVF, VSAV,
    !                         ACOR, WM, IWM, F, JAC, PSOL, VNLS, RPAR, IPAR
    !  Call sequence output -- YH, ACOR, WM, IWM
    !  COMMON block variables accessed:
    !      /DVOD01/  ACNRM, EL(13), H, HMIN, HMXI, HNEW, HSCAL, RC, TAU(13),
    !                TQ(5), TN, JCUR, JSTART, KFLAG, KUTH,
    !                L, LMAX, MAXORD, N, NEWQ, NQ, NQWAIT
    !      /DVOD02/  HU, NCFN, NETF, NFE, NQU, NST
    ! 
    !  Subroutines called by DVSTEP: F, DAXPY, DCOPY, DSCAL,
    !                                DVJUST, VNLS, DVSET
    !  Function routines called by DVSTEP: DVNORM
    ! -----------------------------------------------------------------------
    !  DVSTEP performs one step of the integration of an initial value
    !  problem for a system of ordinary differential equations.
    !  DVSTEP calls subroutine VNLS for the solution of the nonlinear system
    !  arising in the time step.  Thus it is independent of the problem
    !  Jacobian structure and the type of nonlinear system solution method.
    !  DVSTEP returns a completion flag KFLAG (in COMMON).
    !  A return with KFLAG = -1 or -2 means either ABS(H) = HMIN or 10
    !  consecutive failures occurred.  On a return with KFLAG negative,
    !  the values of TN and the YH array are as of the beginning of the last
    !  step, and H is the last step size attempted.
    ! 
    !  Communication with DVSTEP is done with the following variables:
    ! 
    !  Y      = An array of length N used for the dependent variable vector.
    !  YH     = An LDYH by LMAX array containing the dependent variables
    !           and their approximate scaled derivatives, where
    !           LMAX = MAXORD + 1.  YH(i,j+1) contains the approximate
    !           j-th derivative of y(i), scaled by H**j/factorial(j)
    !           (j = 0,1,...,NQ).  On entry for the first step, the first
    !           two columns of YH must be set from the initial values.
    !  LDYH   = A constant integer .ge. N, the first dimension of YH.
    !           N is the number of ODEs in the system.
    !  YH1    = A one-dimensional array occupying the same space as YH.
    !  EWT    = An array of length N containing multiplicative weights
    !           for local error measurements.  Local errors in y(i) are
    !           compared to 1.0/EWT(i) in various error tests.
    !  SAVF   = An array of working storage, of length N.
    !           also used for input of YH(*,MAXORD+2) when JSTART = -1
    !           and MAXORD .lt. the current order NQ.
    !  VSAV   = A work array of length N passed to subroutine VNLS.
    !  ACOR   = A work array of length N, used for the accumulated
    !           corrections.  On a successful return, ACOR(i) contains
    !           the estimated one-step local error in y(i).
    !  WM,IWM = Real and integer work arrays associated with matrix
    !           operations in VNLS.
    !  F      = Dummy name for the user supplied subroutine for f.
    !  JAC    = Dummy name for the user supplied Jacobian subroutine.
    !  PSOL   = Dummy name for the subroutine passed to VNLS, for
    !           possible use there.
    !  VNLS   = Dummy name for the nonlinear system solving subroutine,
    !           whose real name is dependent on the method used.
    !  RPAR, IPAR = Dummy names for user's real and integer work arrays.
    ! -----------------------------------------------------------------------
    EXTERNAL F, JAC, PSOL, VNLS
    type(dvode_t) :: dvode_state
    real(dp_t) :: Y(dvode_state % N), YH(LDYH, dvode_state % LMAX)
    real(dp_t) :: YH1(:), EWT(:), SAVF(:)
    real(dp_t) :: VSAV(:), ACOR(:), WM(:), RPAR(:)
    integer    :: LDYH, IWM(:), IPAR(:)
      
    real(dp_t) :: CNQUOT, DDN, DSM, DUP, TOLD
    real(dp_t) :: ETAQ, ETAQM1, ETAQP1, FLOTL, R
    integer    :: I, I1, I2, IBACK, J, JB, NCF, NFLAG

    ! Parameter declarations
    integer, parameter :: KFC = -3
    integer, parameter :: KFH = -7
    integer, parameter :: MXNCF = 10
    real(dp_t), parameter :: ADDON = 1.0D-6
    real(dp_t), parameter :: BIAS1 = 6.0D0
    real(dp_t), parameter :: BIAS2 = 6.0D0
    real(dp_t), parameter :: BIAS3 = 10.0D0
    real(dp_t), parameter :: ETACF = 0.25D0
    real(dp_t), parameter :: ETAMIN = 0.1D0
    real(dp_t), parameter :: ETAMXF = 0.2D0
    real(dp_t), parameter :: ETAMX1 = 1.0D4
    real(dp_t), parameter :: ETAMX2 = 10.0D0
    real(dp_t), parameter :: ETAMX3 = 10.0D0
    real(dp_t), parameter :: ONEPSM = 1.00001D0
    real(dp_t), parameter :: THRESH = 1.5D0

    ETAQ   = ONE
    ETAQM1 = ONE

    dvode_state % KFLAG = 0
    TOLD = dvode_state % TN
    NCF = 0
    dvode_state % JCUR = 0
    NFLAG = 0
    IF (dvode_state % JSTART .GT. 0) GO TO 20
    IF (dvode_state % JSTART .EQ. -1) GO TO 100
    ! -----------------------------------------------------------------------
    !  On the first call, the order is set to 1, and other variables are
    !  initialized.  ETAMAX is the maximum ratio by which H can be increased
    !  in a single step.  It is normally 10, but is larger during the
    !  first step to compensate for the small initial H.  If a failure
    !  occurs (in corrector convergence or error test), ETAMAX is set to 1
    !  for the next increase.
    ! -----------------------------------------------------------------------
    dvode_state % LMAX = dvode_state % MAXORD + 1
    dvode_state % NQ = 1
    dvode_state % L = 2
    dvode_state % NQNYH = dvode_state % NQ*LDYH
    dvode_state % TAU(1) = dvode_state % H
    dvode_state % PRL1 = ONE
    dvode_state % RC = ZERO
    dvode_state % ETAMAX = ETAMX1
    dvode_state % NQWAIT = 2
    dvode_state % HSCAL = dvode_state % H
    GO TO 200
    ! -----------------------------------------------------------------------
    !  Take preliminary actions on a normal continuation step (JSTART.GT.0).
    !  If the driver changed H, then ETA must be reset and NEWH set to 1.
    !  If a change of order was dictated on the previous step, then
    !  it is done here and appropriate adjustments in the history are made.
    !  On an order decrease, the history array is adjusted by DVJUST.
    !  On an order increase, the history array is augmented by a column.
    !  On a change of step size H, the history array YH is rescaled.
    ! -----------------------------------------------------------------------
20  CONTINUE
    IF (dvode_state % KUTH .EQ. 1) THEN
       dvode_state % ETA = MIN(dvode_state % ETA,dvode_state % H/dvode_state % HSCAL)
       dvode_state % NEWH = 1
    ENDIF
50  IF (dvode_state % NEWH .EQ. 0) GO TO 200
    IF (dvode_state % NEWQ .EQ. dvode_state % NQ) GO TO 150
    IF (dvode_state % NEWQ .LT. dvode_state % NQ) THEN
       CALL DVJUST (YH, LDYH, -1, dvode_state)
       dvode_state % NQ = dvode_state % NEWQ
       dvode_state % L = dvode_state % NQ + 1
       dvode_state % NQWAIT = dvode_state % L
       GO TO 150
    ENDIF
    IF (dvode_state % NEWQ .GT. dvode_state % NQ) THEN
       CALL DVJUST (YH, LDYH, 1, dvode_state)
       dvode_state % NQ = dvode_state % NEWQ
       dvode_state % L = dvode_state % NQ + 1
       dvode_state % NQWAIT = dvode_state % L
       GO TO 150
    ENDIF
    ! -----------------------------------------------------------------------
    !  The following block handles preliminaries needed when JSTART = -1.
    !  If N was reduced, zero out part of YH to avoid undefined references.
    !  If MAXORD was reduced to a value less than the tentative order NEWQ,
    !  then NQ is set to MAXORD, and a new H ratio ETA is chosen.
    !  Otherwise, we take the same preliminary actions as for JSTART .gt. 0.
    !  In any case, NQWAIT is reset to L = NQ + 1 to prevent further
    !  changes in order for that many steps.
    !  The new H ratio ETA is limited by the input H if KUTH = 1,
    !  by HMIN if KUTH = 0, and by HMXI in any case.
    !  Finally, the history array YH is rescaled.
    ! -----------------------------------------------------------------------
100 CONTINUE
    dvode_state % LMAX = dvode_state % MAXORD + 1
    IF (dvode_state % N .EQ. LDYH) GO TO 120
    I1 = 1 + (dvode_state % NEWQ + 1)*LDYH
    I2 = (dvode_state % MAXORD + 1)*LDYH
    IF (I1 .GT. I2) GO TO 120
    do I = I1, I2
       YH1(I) = ZERO
    end do
120 IF (dvode_state % NEWQ .LE. dvode_state % MAXORD) GO TO 140
    FLOTL = REAL(dvode_state % LMAX)
    IF (dvode_state % MAXORD .LT. dvode_state % NQ-1) THEN
       DDN = DVNORM (dvode_state % N, SAVF, EWT)/dvode_state % TQ(1)
       dvode_state % ETA = ONE/((BIAS1*DDN)**(ONE/FLOTL) + ADDON)
    ENDIF
    IF (dvode_state % MAXORD .EQ. dvode_state % NQ .AND. dvode_state % NEWQ .EQ. dvode_state % NQ+1) dvode_state % ETA = ETAQ
    IF (dvode_state % MAXORD .EQ. dvode_state % NQ-1 .AND. dvode_state % NEWQ .EQ. dvode_state % NQ+1) THEN
       dvode_state % ETA = ETAQM1
       CALL DVJUST (YH, LDYH, -1, dvode_state)
    ENDIF
    IF (dvode_state % MAXORD .EQ. dvode_state % NQ-1 .AND. dvode_state % NEWQ .EQ. dvode_state % NQ) THEN
       DDN = DVNORM (dvode_state % N, SAVF, EWT)/dvode_state % TQ(1)
       dvode_state % ETA = ONE/((BIAS1*DDN)**(ONE/FLOTL) + ADDON)
       CALL DVJUST (YH, LDYH, -1, dvode_state)
    ENDIF
    dvode_state % ETA = MIN(dvode_state % ETA,ONE)
    dvode_state % NQ = dvode_state % MAXORD
    dvode_state % L = dvode_state % LMAX
140 IF (dvode_state % KUTH .EQ. 1) dvode_state % ETA = MIN(dvode_state % ETA,ABS(dvode_state % H/dvode_state % HSCAL))
    IF (dvode_state % KUTH .EQ. 0) dvode_state % ETA = MAX(dvode_state % ETA,dvode_state % HMIN/ABS(dvode_state % HSCAL))
    dvode_state % ETA = dvode_state % ETA/MAX(ONE,ABS(dvode_state % HSCAL)*dvode_state % HMXI*dvode_state % ETA)
    dvode_state % NEWH = 1
    dvode_state % NQWAIT = dvode_state % L
    IF (dvode_state % NEWQ .LE. dvode_state % MAXORD) GO TO 50
    ! Rescale the history array for a change in H by a factor of ETA. ------
150 R = ONE
    do J = 2, dvode_state % L
       R = R * dvode_state % ETA
       CALL DSCAL (dvode_state % N, R, YH(1,J), 1 )
    end do
    dvode_state % H = dvode_state % HSCAL * dvode_state % ETA
    dvode_state % HSCAL = dvode_state % H
    dvode_state % RC = dvode_state % RC * dvode_state % ETA
    dvode_state % NQNYH = dvode_state % NQ*LDYH
    ! -----------------------------------------------------------------------
    !  This section computes the predicted values by effectively
    !  multiplying the YH array by the Pascal triangle matrix.
    !  DVSET is called to calculate all integration coefficients.
    !  RC is the ratio of new to old values of the coefficient H/EL(2)=h/l1.
    ! -----------------------------------------------------------------------
200 dvode_state % TN = dvode_state % TN + dvode_state % H
    I1 = dvode_state % NQNYH + 1
    do JB = 1, dvode_state % NQ
       I1 = I1 - LDYH
       do I = I1, dvode_state % NQNYH
          YH1(I) = YH1(I) + YH1(I+LDYH)
       end do
    end do
    CALL DVSET(dvode_state)
    dvode_state % RL1 = ONE/dvode_state % EL(2)
    dvode_state % RC = dvode_state % RC * (dvode_state % RL1/dvode_state % PRL1)
    dvode_state % PRL1 = dvode_state % RL1
    ! 
    !  Call the nonlinear system solver. ------------------------------------
    ! 
    CALL VNLS (Y, YH, LDYH, VSAV, SAVF, EWT, ACOR, IWM, WM, &
         F, JAC, PSOL, NFLAG, RPAR, IPAR, dvode_state)

    IF (NFLAG .EQ. 0) GO TO 450
    ! -----------------------------------------------------------------------
    !  The VNLS routine failed to achieve convergence (NFLAG .NE. 0).
    !  The YH array is retracted to its values before prediction.
    !  The step size H is reduced and the step is retried, if possible.
    !  Otherwise, an error exit is taken.
    ! -----------------------------------------------------------------------
    NCF = NCF + 1
    dvode_state % NCFN = dvode_state % NCFN + 1
    dvode_state % ETAMAX = ONE
    dvode_state % TN = TOLD
    I1 = dvode_state % NQNYH + 1
    do JB = 1, dvode_state % NQ
       I1 = I1 - LDYH
       do I = I1, dvode_state % NQNYH
          YH1(I) = YH1(I) - YH1(I+LDYH)
       end do
    end do
    IF (NFLAG .LT. -1) GO TO 680
    IF (ABS(dvode_state % H) .LE. dvode_state % HMIN*ONEPSM) GO TO 670
    IF (NCF .EQ. MXNCF) GO TO 670
    dvode_state % ETA = ETACF
    dvode_state % ETA = MAX(dvode_state % ETA,dvode_state % HMIN/ABS(dvode_state % H))
    NFLAG = -1
    GO TO 150
    ! -----------------------------------------------------------------------
    !  The corrector has converged (NFLAG = 0).  The local error test is
    !  made and control passes to statement 500 if it fails.
    ! -----------------------------------------------------------------------
450 CONTINUE
    DSM = dvode_state % ACNRM/dvode_state % TQ(2)
    IF (DSM .GT. ONE) GO TO 500
    ! -----------------------------------------------------------------------
    !  After a successful step, update the YH and TAU arrays and decrement
    !  NQWAIT.  If NQWAIT is then 1 and NQ .lt. MAXORD, then ACOR is saved
    !  for use in a possible order increase on the next step.
    !  If ETAMAX = 1 (a failure occurred this step), keep NQWAIT .ge. 2.
    ! -----------------------------------------------------------------------
    dvode_state % KFLAG = 0
    dvode_state % NST = dvode_state % NST + 1
    dvode_state % HU = dvode_state % H
    dvode_state % NQU = dvode_state % NQ
    do IBACK = 1, dvode_state % NQ
       I = dvode_state % L - IBACK
       dvode_state % TAU(I+1) = dvode_state % TAU(I)
    end do
    dvode_state % TAU(1) = dvode_state % H
    do J = 1, dvode_state % L
       CALL DAXPY (dvode_state % N, dvode_state % EL(J), ACOR, 1, YH(1,J), 1 )
    end do
    dvode_state % NQWAIT = dvode_state % NQWAIT - 1
    IF ((dvode_state % L .EQ. dvode_state % LMAX) .OR. (dvode_state % NQWAIT .NE. 1)) GO TO 490
    CALL DCOPY (dvode_state % N, ACOR, 1, YH(1,dvode_state % LMAX), 1 )
    dvode_state % CONP = dvode_state % TQ(5)
490 IF (dvode_state % ETAMAX .NE. ONE) GO TO 560
    IF (dvode_state % NQWAIT .LT. 2) dvode_state % NQWAIT = 2
    dvode_state % NEWQ = dvode_state % NQ
    dvode_state % NEWH = 0
    dvode_state % ETA = ONE
    dvode_state % HNEW = dvode_state % H
    GO TO 690
    ! -----------------------------------------------------------------------
    !  The error test failed.  KFLAG keeps track of multiple failures.
    !  Restore TN and the YH array to their previous values, and prepare
    !  to try the step again.  Compute the optimum step size for the
    !  same order.  After repeated failures, H is forced to decrease
    !  more rapidly.
    ! -----------------------------------------------------------------------
500 dvode_state % KFLAG = dvode_state % KFLAG - 1
    dvode_state % NETF = dvode_state % NETF + 1
    NFLAG = -2
    dvode_state % TN = TOLD
    I1 = dvode_state % NQNYH + 1
    do JB = 1, dvode_state % NQ
       I1 = I1 - LDYH
       do I = I1, dvode_state % NQNYH
          YH1(I) = YH1(I) - YH1(I+LDYH)
       end do
    end do
    IF (ABS(dvode_state % H) .LE. dvode_state % HMIN*ONEPSM) GO TO 660
    dvode_state % ETAMAX = ONE
    IF (dvode_state % KFLAG .LE. KFC) GO TO 530
    ! Compute ratio of new H to current H at the current order. ------------
    FLOTL = REAL(dvode_state % L)
    dvode_state % ETA = ONE/((BIAS2*DSM)**(ONE/FLOTL) + ADDON)
    dvode_state % ETA = MAX(dvode_state % ETA,dvode_state % HMIN/ABS(dvode_state % H),ETAMIN)
    IF ((dvode_state % KFLAG .LE. -2) .AND. (dvode_state % ETA .GT. ETAMXF)) dvode_state % ETA = ETAMXF
    GO TO 150
    ! -----------------------------------------------------------------------
    !  Control reaches this section if 3 or more consecutive failures
    !  have occurred.  It is assumed that the elements of the YH array
    !  have accumulated errors of the wrong order.  The order is reduced
    !  by one, if possible.  Then H is reduced by a factor of 0.1 and
    !  the step is retried.  After a total of 7 consecutive failures,
    !  an exit is taken with KFLAG = -1.
    ! -----------------------------------------------------------------------
530 IF (dvode_state % KFLAG .EQ. KFH) GO TO 660
    IF (dvode_state % NQ .EQ. 1) GO TO 540
    dvode_state % ETA = MAX(ETAMIN,dvode_state % HMIN/ABS(dvode_state % H))
    CALL DVJUST (YH, LDYH, -1, dvode_state)
    dvode_state % L = dvode_state % NQ
    dvode_state % NQ = dvode_state % NQ - 1
    dvode_state % NQWAIT = dvode_state % L
    GO TO 150
540 dvode_state % ETA = MAX(ETAMIN,dvode_state % HMIN/ABS(dvode_state % H))
    dvode_state % H = dvode_state % H * dvode_state % ETA
    dvode_state % HSCAL = dvode_state % H
    dvode_state % TAU(1) = dvode_state % H
    CALL F (dvode_state % N, dvode_state % TN, Y, SAVF, RPAR, IPAR)
    dvode_state % NFE = dvode_state % NFE + 1
    do I = 1, dvode_state % N
       YH(I,2) = dvode_state % H*SAVF(I)
    end do
    dvode_state % NQWAIT = 10
    GO TO 200
    ! -----------------------------------------------------------------------
    !  If NQWAIT = 0, an increase or decrease in order by one is considered.
    !  Factors ETAQ, ETAQM1, ETAQP1 are computed by which H could
    !  be multiplied at order q, q-1, or q+1, respectively.
    !  The largest of these is determined, and the new order and
    !  step size set accordingly.
    !  A change of H or NQ is made only if H increases by at least a
    !  factor of THRESH.  If an order change is considered and rejected,
    !  then NQWAIT is set to 2 (reconsider it after 2 steps).
    ! -----------------------------------------------------------------------
    !  Compute ratio of new H to current H at the current order. ------------
560 FLOTL = REAL(dvode_state % L)
    ETAQ = ONE/((BIAS2*DSM)**(ONE/FLOTL) + ADDON)
    IF (dvode_state % NQWAIT .NE. 0) GO TO 600
    dvode_state % NQWAIT = 2
    ETAQM1 = ZERO
    IF (dvode_state % NQ .EQ. 1) GO TO 570
    ! Compute ratio of new H to current H at the current order less one. ---
    DDN = DVNORM (dvode_state % N, YH(1,dvode_state % L), EWT)/dvode_state % TQ(1)
    ETAQM1 = ONE/((BIAS1*DDN)**(ONE/(FLOTL - ONE)) + ADDON)
570 ETAQP1 = ZERO
    IF (dvode_state % L .EQ. dvode_state % LMAX) GO TO 580
    ! Compute ratio of new H to current H at current order plus one. -------
    CNQUOT = (dvode_state % TQ(5)/dvode_state % CONP)*(dvode_state % H/dvode_state % TAU(2))**dvode_state % L
    do I = 1, dvode_state % N
       SAVF(I) = ACOR(I) - CNQUOT*YH(I,dvode_state % LMAX)
    end do
    DUP = DVNORM (dvode_state % N, SAVF, EWT)/dvode_state % TQ(3)
    ETAQP1 = ONE/((BIAS3*DUP)**(ONE/(FLOTL + ONE)) + ADDON)
580 IF (ETAQ .GE. ETAQP1) GO TO 590
    IF (ETAQP1 .GT. ETAQM1) GO TO 620
    GO TO 610
590 IF (ETAQ .LT. ETAQM1) GO TO 610
600 dvode_state % ETA = ETAQ
    dvode_state % NEWQ = dvode_state % NQ
    GO TO 630
610 dvode_state % ETA = ETAQM1
    dvode_state % NEWQ = dvode_state % NQ - 1
    GO TO 630
620 dvode_state % ETA = ETAQP1
    dvode_state % NEWQ = dvode_state % NQ + 1
    CALL DCOPY (dvode_state % N, ACOR, 1, YH(1,dvode_state % LMAX), 1)
    ! Test tentative new H against THRESH, ETAMAX, and HMXI, then exit. ----
630 IF (dvode_state % ETA .LT. THRESH .OR. dvode_state % ETAMAX .EQ. ONE) GO TO 640
    dvode_state % ETA = MIN(dvode_state % ETA,dvode_state % ETAMAX)
    dvode_state % ETA = dvode_state % ETA/MAX(ONE,ABS(dvode_state % H)*dvode_state % HMXI * dvode_state % ETA)
    dvode_state % NEWH = 1
    dvode_state % HNEW = dvode_state % H * dvode_state % ETA
    GO TO 690
640 dvode_state % NEWQ = dvode_state % NQ
    dvode_state % NEWH = 0
    dvode_state % ETA = ONE
    dvode_state % HNEW = dvode_state % H
    GO TO 690
    ! -----------------------------------------------------------------------
    !  All returns are made through this section.
    !  On a successful return, ETAMAX is reset and ACOR is scaled.
    ! -----------------------------------------------------------------------
660 dvode_state % KFLAG = -1
    GO TO 720
670 dvode_state % KFLAG = -2
    GO TO 720
680 IF (NFLAG .EQ. -2) dvode_state % KFLAG = -3
    IF (NFLAG .EQ. -3) dvode_state % KFLAG = -4
    GO TO 720
690 dvode_state % ETAMAX = ETAMX3
    IF (dvode_state % NST .LE. 10) dvode_state % ETAMAX = ETAMX2
    R = ONE/dvode_state % TQ(2)
    CALL DSCAL (dvode_state % N, R, ACOR, 1)
720 dvode_state % JSTART = 1
    RETURN
  end subroutine dvstep
      
end module dvode_module
