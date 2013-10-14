Module TempMod
   
   USE NWTC_Library
   USE SubDyn_Types
   USE SubDyn_Output

   IMPLICIT NONE

   PRIVATE
   
   PUBLIC :: SubDyn_UpdateStates                   ! Loose coupling routine for solving for constraint states, integrating
                                                 
   PUBLIC :: SubDyn_CalcOutput                     ! Routine for computing outputs

   PUBLIC :: SubDyn_CalcContStateDeriv              ! Tight coupling routine for computing derivatives of continuous states
   
   CONTAINS
   
    !----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE SubDyn_UpdateStates( t, n, Inputs, InputTimes, p, x, xd, z, OtherState, ErrStat, ErrMsg )
! Loose coupling routine for solving for constraint states, integrating continuous states, and updating discrete states
! Constraint states are solved for input time, t; Continuous and discrete states are updated for t + Interval
! A guess for the contstraint states at t + Interval is also calculated.
!..................................................................................................................................

      REAL(DbKi),                         INTENT(IN   ) :: t               ! Current simulation time in seconds
      INTEGER(IntKi),                     INTENT(IN   ) :: n               ! Current step of the simulation: t = n*Interval
      TYPE(SD_InputType),                 INTENT(INOUT) :: Inputs(:)       ! Inputs at Times
      REAL(DbKi),                         INTENT(IN   ) :: InputTimes(:)   ! Times in seconds associated with Inputs
      TYPE(SD_ParameterType),             INTENT(IN   ) :: p               ! Parameters
      TYPE(SD_ContinuousStateType),       INTENT(INOUT) :: x               ! Input: Continuous states at t;
                                                                           !   Output: Continuous states at t + Interval
      TYPE(SD_DiscreteStateType),         INTENT(INOUT) :: xd              ! Input: Discrete states at t;
                                                                           !   Output: Discrete states at t + Interval
      TYPE(SD_ConstraintStateType),       INTENT(INOUT) :: z               ! Input: Initial guess of constraint states at t;
                                                                           !   Output: Constraint states at t
      TYPE(SD_OtherStateType),            INTENT(INOUT) :: OtherState      ! Other/optimization states
      INTEGER(IntKi),                     INTENT(  OUT) :: ErrStat         ! Error status of the operation
      CHARACTER(*),                       INTENT(  OUT) :: ErrMsg          ! Error message if ErrStat /= ErrID_None

         ! Local variables

      TYPE(SD_ContinuousStateType)                 :: dxdt        ! Continuous state derivatives at t
      TYPE(SD_InputType)                           :: u           ! Instantaneous inputs
      INTEGER(IntKi)                               :: ErrStat2    ! Error status of the operation (occurs after initial error)
      CHARACTER(LEN(ErrMsg))                       :: ErrMsg2     ! Error message if ErrStat2 /= ErrID_None

      
         ! Initialize variables

      ErrStat   = ErrID_None           ! no error has occurred
      ErrMsg    = ""
            
      
         ! Get the inputs, based on the array of values sent by the glue code:
         
    CALL SD_CopyInput( Inputs(1), u, MESH_NEWCOPY, ErrStat, ErrMsg )          ! bjj: this will need to be changed when the routine is implemented

    
      IF ( ErrStat >= AbortErrLev ) THEN
         RETURN
      END IF
      
   

         ! Get first time derivatives of continuous states (dxdt): NOT SURE THIS IS NEEDED SINCE it is BEING CALLED WITHIN THE INTEGRATOR

      CALL SubDyn_CalcContStateDeriv( t, u, p, x, xd, z, OtherState, dxdt, ErrStat, ErrMsg )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL SD_DestroyContState( dxdt, ErrStat2, ErrMsg2)
         ErrMsg = TRIM(ErrMsg)//' '//TRIM(ErrMsg2)
         RETURN
      END IF

     !x%qm=x%qm + dxdt%qm*p%SDDELTAt
     !x%qmdot=x%qmdot+dxdt%qmdot*p%SDDELTAt
            ! Integrate (update) continuous states (x) here:
        !LET US USE INTEGRATOR
        
      !x = function of dxdt and x
      IF (p%IntMethod .eq. 1) THEN
 
         CALL SD_RK4( t, n, Inputs, InputTimes, p, x, xd, z, OtherState, ErrStat, ErrMsg )

      ELSEIF (p%IntMethod .eq. 2) THEN

         CALL SD_AB4( t, n, Inputs, InputTimes, p, x, xd, z, OtherState, ErrStat, ErrMsg )

      ELSE 
          
         CALL SD_ABM4( t, n, Inputs, InputTimes, p, x, xd, z, OtherState, ErrStat, ErrMsg )

      END IF

    !Assign the acceleration to the x variable since it will be used for output file purposes for SSqmdd01-99, and dxdt will disappear
    x%qmdotdot=dxdt%qmdot
  
    ! Destroy dxdt because it is not necessary for the rest of the subroutine
    
      CALL SD_DestroyContState( dxdt, ErrStat, ErrMsg)

      IF ( ErrStat >= AbortErrLev ) RETURN


END SUBROUTINE SubDyn_UpdateStates


!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE SubDyn_CalcOutput( t, u, p, x, xd, z, OtherState, y, ErrStat, ErrMsg )
! Routine for computing outputs, used in both loose and tight coupling.
!..................................................................................................................................

      REAL(DbKi),                   INTENT(IN   )  :: t           ! Current simulation time in seconds
      TYPE(SD_InputType),           INTENT(IN   )  :: u           ! Inputs at t
      TYPE(SD_ParameterType),       INTENT(IN   )  :: p           ! Parameters
      TYPE(SD_ContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at t
      TYPE(SD_DiscreteStateType),   INTENT(IN   )  :: xd          ! Discrete states at t
      TYPE(SD_ConstraintStateType), INTENT(IN   )  :: z           ! Constraint states at t
      TYPE(SD_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
      TYPE(SD_OutputType),          INTENT(INOUT)  :: y           ! Outputs computed at t (Input only so that mesh con-
                                                                       !   nectivity information does not have to be recalculated)
      INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      !locals
      INTEGER(IntKi)                               :: L1,L2,L3,L4,L5      ! partial Lengths of state and input arrays
      INTEGER(IntKi)                               :: I,J,K,K2,L      ! Counters
      INTEGER(IntKi)                               :: UnJck      ! JCkOutputFile
      REAL(DbKi), DIMENSION(12)                    :: Junk      ! temporary storage for output stuff
      REAL(ReKi)                                   :: AllOuts(0:p%MaxOutPts+p%OutAllInt*p%OutAllDims)
      INTEGER(IntKi), SAVE                         :: Decimat=0
      REAL(ReKi)                                   :: rotations(3)
      REAL(ReKi)                                   :: u_TP(6), udot_TP(6), udotdot_TP(6)
      REAL(ReKi)                                   :: UFL(p%DOFL)
      REAL(ReKi)                                   :: U_L(p%NNodes_L*6), Udot_L(p%NNodes_L*6)
      REAL(ReKi)                                   :: Y1(6)
      INTEGER(IntKi)                               :: startDOF
      REAL(ReKi)                                   :: DCM(3,3)
      ! Initialize ErrStat
      ErrStat = ErrID_None
      ErrMsg  = ""

      !x%qm should be allocated already
          
      L1=p%qmL /2   !Length of array qm (half length of x)
           
      L2=p%TPdofL*2+1     !start index for array U_TP_dotdot (Subarray of u)
      L3=p%TPdofL*3+1       !start index for array FL (Subarray of u)
   
      L4=p%URbarL+p%DOFL+1    !start index for subarray UR_dot
      L5=2*p%URbarL+p%DOFL+1    !start index for subarray UL_dot
      
      
      
        
      !u_TP       = u%UFL(1:p%TPdofL)
      
         ! Compute the small rotation angles given the input direction cosine matrix
      rotations  = GetSmllRotAngs(u%TPMesh%Orientation(:,:,1), ErrStat, Errmsg)
      u_TP       = (/u%TPMesh%TranslationDisp(:,1), rotations/)
   
      !udot_TP    = u%UFL(p%TPdofL+1:L2-1)
      udot_TP    = (/u%TPMesh%TranslationVel(:,1), u%TPMesh%RotationVel(:,1)/)
      
      !udotdot_TP = u%UFL(L2:L3-1)
      udotdot_TP = (/u%TPMesh%TranslationAcc(:,1), u%TPMesh%RotationAcc(:,1)/)
      
      
      
      
      
      !Y2= Structure node displacements and velocities  
      !OtherState%Y2(1:p%URbarL)        = matmul(p%D2_11,u%UFL(1:p%TPdofL)) 
      !OtherState%Y2(p%URbarL+1:L4-1)   = matmul(p%C2_21,x%qm)  
      !write(*, *)  ' UL(9)=',OtherState%Y2(9), 'qm(1)=', x%qm(1) !DEBUG!debug
      !OtherState%Y2(p%URbarL+1:L4-1)   =OtherState%Y2(p%URbarL+1:L4-1)   + matmul(p%D2_21,u%UFL(1:p%TPdofL)) 
     ! OtherState%Y2(L4:L5-1)        = matmul(p%D2_32,u%UFL(p%TPdofL+1:L2-1)) 
     ! OtherState%Y2(L5:p%Y2L)       = matmul(p%C2_42,x%qmdot)  + matmul(p%D2_42,u%UFL(p%TPdofL+1:L2-1)) 
      
          
         ! TODO: The URbarL qu      
      OtherState%Y2(1:p%URbarL)        = matmul( p%D2_11, u_TP )
      
      
      OtherState%Y2(p%URbarL+1:L4-1) = matmul( p%C2_21, x%qm )  + matmul( p%D2_21, u_TP ) 
      U_L =     matmul( p%C2_21, x%qm )  + matmul( p%D2_21, u_TP ) 
      
      
      OtherState%Y2(L4:L5-1)        = matmul( p%D2_32, udot_TP ) 
      
      
      OtherState%Y2(L5:p%Y2L) = matmul( p%C2_42, x%qmdot )  + matmul( p%D2_42, udot_TP ) 
      Udot_L = matmul( p%C2_42, x%qmdot )  + matmul( p%D2_42, udot_TP ) 
      
      !UFL =  u%UFL(L3:p%uL)
      
         ! Place the outputs onto the Y2 output mesh 
      DO I = 1, p%NNodes_L       
             
            ! starting index in the master arrays for the current node    
         startDOF = (I-1)*6 + 1
         
             ! Construct UFL array from the Force and Moment fields of the input mesh
         UFL ( startDOF   : startDOF + 2 ) = u%LMesh%Force (:,I)
         UFL ( startDOF+3 : startDOF + 5 ) = u%LMesh%Moment(:,I)
         
            ! Construct the direction cosine matrix given the output angles
         CALL SmllRotTrans( 'U_L input angles', U_L(startDOF + 3), U_L(startDOF + 4), U_L(startDOF + 5), DCM, '', ErrStat, ErrMsg )
         !TODO: Add error handling
         
            ! Y2 = Interior node displacements and velocities  for use as inputs to HydroDyn
         y%Y2mesh%TranslationDisp (:,I)     = U_L    ( startDOF     : startDOF + 2 )
         y%Y2mesh%Orientation     (:,:,I)   = DCM
         y%Y2mesh%TranslationVel  (:,I)     = Udot_L ( startDOF     : startDOF + 2 )
         y%Y2mesh%RotationVel     (:,I)     = Udot_L ( startDOF + 3 : startDOF + 5 )
         
      END DO
      
      
      !Y1= TP reaction Forces, i.e. force that the jacket exerts onto the TP and above  
      
      Y1 = -( matmul(p%C1_11, x%qm)       + matmul(p%C1_12,x%qmdot)    + matmul(p%D1_11, u_TP) + &
              matmul(p%D1_13, udotdot_TP) + matmul(p%D1_14, UFL)       + p%FY )
      
      !Y1 = -( matmul(p%C1_11,x%qm) + matmul(p%C1_12,x%qmdot)  + matmul(p%D1_11,u%UFL(1:p%TPdofL)) + &
      !        matmul(p%D1_13,u%UFL(L2:L3-1))+ matmul(p%D1_14,u%UFL(L3:p%uL)) + p%FY )
      
      y%Y1Mesh%Force (:,1) = Y1(1:3)
      y%Y1Mesh%Moment(:,1) = Y1(4:6)

      
      
      !Now I need to calculate the Output File Stuff
  
   !Calculate accelerations even though they may not be necessary, in the future we may put 
   ! a condition on the type of input to speed up, if forces are requested we need them
    ! OtherState%Udotdot(1:p%URbarL)= matmul(p%Dbar_13,u%UFL(L2:L3-1))
      
    ! OtherState%Udotdot(p%URbarL+1:p%URbarL+p%DOFL)= matmul(p%Cbar_21,x%qm) + matmul(p%Cbar_22,x%qmdot) + & 
    !                                         matmul(p%Dbar_23,u%UFL(L2:L3-1)) +  matmul(p%Dbar_24,u%UFL(L3:p%uL)) + &
    !                                               p%Fbar_21
OtherState%Udotdot(1:p%URbarL) = matmul( p%Dbar_13, udotdot_TP ) ! This is Ubardotdot_R which is not being used at this point. GJH 5/23/13
   
   
   OtherState%Udotdot(p%URbarL+1:p%URbarL+p%DOFL) =                               & 
                  matmul( p%Cbar_21, x%qm       ) + matmul( p%Cbar_22, x%qmdot ) + & 
                  matmul( p%Dbar_23, udotdot_TP ) + matmul( p%Dbar_24, UFL     ) + &
                          p%Fbar_21
                                  !_____________________________________!
                                ! CALCULATE OUTPUT TO BE WRITTEN TO FILE !
                                  !_____________________________________!
                                
         ! OutSwtch determines whether or not to actually output results via the WriteOutput array
         ! 1 = SubDyn will generate an output file of its own.  2 = the caller will handle the outputs, but
         ! SubDyn needs to provide them.  3 = Both 1 and 2, 0 = No one needs the SubDyn outputs provided
         ! via the WriteOutput array.
      
      IF ((Decimat .EQ. p%OutDec) .OR. (Decimat .EQ. 0))  THEN
       Decimat=1  !reset counter
       IF ( p%OutSwtch > 0 ) THEN
         
            ! Map calculated results into the AllOuts Array + perform averaging and all necessary extra calculations
         CALL SDOut_MapOutputs(t, u,p,x,y, OtherState, AllOuts, ErrStat, ErrMsg)
         
         ! Put the output data in the WriteOutput array
         DO I = 1,p%NumOuts+p%OutAllInt*p%OutAllDims
            y%WriteOutput(I) = p%OutParam(I)%SignM * AllOuts( p%OutParam(I)%Indx )
         END DO
         
         ! Generate output into the output file
            
         IF ( p%OutSwtch == 1 .OR. p%OutSwtch == 3 ) THEN
            CALL SDOut_WriteOutputs( p%UnJckF, t, y, p, ErrStat, ErrMsg )         
         END IF
         
        ENDIF           
    ELSE      
      Decimat=Decimat+1
    ENDIF
  
END SUBROUTINE SubDyn_CalcOutput
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE SubDyn_CalcContStateDeriv( t, u, p, x, xd, z, OtherState, dxdt, ErrStat, ErrMsg )
! Tight coupling routine for computing derivatives of continuous states
!..................................................................................................................................

      REAL(DbKi),                        INTENT(IN   )  :: t           ! Current simulation time in seconds
      TYPE(SD_InputType),           INTENT(IN   )  :: u           ! Inputs at t
      TYPE(SD_ParameterType),       INTENT(IN   )  :: p           ! Parameters
      TYPE(SD_ContinuousStateType), INTENT(INOUT)  :: x           ! Continuous states at t
      TYPE(SD_DiscreteStateType),   INTENT(IN   )  :: xd          ! Discrete states at t
      TYPE(SD_ConstraintStateType), INTENT(IN   )  :: z           ! Constraint states at t
      TYPE(SD_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
      TYPE(SD_ContinuousStateType), INTENT(  OUT)  :: dxdt        ! Continuous state derivatives at t
      INTEGER(IntKi),                    INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                      INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None
    
      INTEGER(IntKi)                               ::L1,L2,L3,L4      ! partial Lengths of state and input arrays
      REAL(ReKi)                                   :: udotdot_TP(6)
      REAL(ReKi)                                   :: UFL(p%DOFL)
      INTEGER(IntKi)                               :: I, startDOF
         ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = ""

      !x%qm should be allocated already
          
      L1=p%qmL /2   !Length of array qm (half length of x)
      L4=p%uL       !Length of array u
      
      L2=p%TPdofL*2+1     !start index for array U_TP_dotdot (Subarray of u)
      L3=p%TPdofL*3+1       !start index for array FL (Subarray of u)
         ! Compute the first time derivatives of the continuous states here:
        
      dxdt%DummyContState = 0
      !How is it possible that we have to check this all the time?
!bjj: INTENT(OUT) automatically deallocates the array on entry.      
      ALLOCATE(dxdt%qm(p%qmL),STAT=ErrStat)  
      ALLOCATE(dxdt%qmdot(p%qmL),STAT=ErrStat)  
      IF ( ErrStat/= 0 ) THEN
         ErrStat = ErrID_Fatal
         ErrMsg  = 'Error allocating states derivatives in SubDyn_CalcContStateDeriv'
         RETURN
      END IF
      !X=Ax + Bu + Fx
      dxdt%qm= x%qmdot

      !udotdot_TP = u%UFL(L2:L3-1)
      udotdot_TP = (/u%TPMesh%TranslationAcc(:,1), u%TPMesh%RotationAcc(:,1)/)
      
      !UFL = u%UFL(L3:L4)
      DO I = 1, p%NNodes_L
         startDOF = (I-1)*6
         UFL( startDOF + 1 : startDOF + 3 ) = u%LMesh%Force (:,I)
         UFL( startDOF + 4 : startDOF + 6 ) = u%LMesh%Moment(:,I)
      END DO
      
      dxdt%qmdot= matmul(p%A_21,x%qm) + matmul(p%A_22,x%qmdot)+ matmul(p%B_23,udotdot_TP)  + matmul(p%B_24,UFL) + p%FX

     ! dxdt%qmdot= matmul(p%A_21,x%qm) + matmul(p%A_22,x%qmdot)+ matmul(p%B_23,u%UFL(L2:L3-1))  + matmul(p%B_24,u%UFL(L3:L4)) + p%FX


END SUBROUTINE SubDyn_CalcContStateDeriv

    !----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE SD_RK4( t, n, u, utimes, p, x, xd, z, OtherState, ErrStat, ErrMsg )
!
! This subroutine implements the fourth-order Runge-Kutta Method (RK4) for numerically integrating ordinary differential equations:
!
!   Let f(t, x) = xdot denote the time (t) derivative of the continuous states (x). 
!   Define constants k1, k2, k3, and k4 as 
!        k1 = dt * f(t        , x_t        )
!        k2 = dt * f(t + dt/2 , x_t + k1/2 )
!        k3 = dt * f(t + dt/2 , x_t + k2/2 ), and
!        k4 = dt * f(t + dt   , x_t + k3   ).
!   Then the continuous states at t = t + dt are
!        x_(t+dt) = x_t + k1/6 + k2/3 + k3/3 + k4/6 + O(dt^5)
!
! For details, see:
! Press, W. H.; Flannery, B. P.; Teukolsky, S. A.; and Vetterling, W. T. "Runge-Kutta Method" and "Adaptive Step Size Control for 
!   Runge-Kutta." §16.1 and 16.2 in Numerical Recipes in FORTRAN: The Art of Scientific Computing, 2nd ed. Cambridge, England: 
!   Cambridge University Press, pp. 704-716, 1992.
!
!..................................................................................................................................

      REAL(DbKi),                     INTENT(IN   )  :: t           ! Current simulation time in seconds
      INTEGER(IntKi),                 INTENT(IN   )  :: n           ! time step number
      TYPE(SD_InputType),             INTENT(INOUT)  :: u(:)        ! Inputs at t
      REAL(DbKi),                     INTENT(IN   )  :: utimes(:)   ! times of input
      TYPE(SD_ParameterType),         INTENT(IN   )  :: p           ! Parameters
      TYPE(SD_ContinuousStateType),   INTENT(INOUT)  :: x           ! Continuous states at t on input at t + dt on output
      TYPE(SD_DiscreteStateType),     INTENT(IN   )  :: xd          ! Discrete states at t
      TYPE(SD_ConstraintStateType),   INTENT(IN   )  :: z           ! Constraint states at t (possibly a guess)
      TYPE(SD_OtherStateType),        INTENT(INOUT)  :: OtherState  ! Other/optimization states
      INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      ! local variables
         
      TYPE(SD_ContinuousStateType)                 :: xdot        ! time derivatives of continuous states      
      TYPE(SD_ContinuousStateType)                 :: k1          ! RK4 constant; see above
      TYPE(SD_ContinuousStateType)                 :: k2          ! RK4 constant; see above 
      TYPE(SD_ContinuousStateType)                 :: k3          ! RK4 constant; see above 
      TYPE(SD_ContinuousStateType)                 :: k4          ! RK4 constant; see above 
      TYPE(SD_ContinuousStateType)                 :: x_tmp       ! Holds temporary modification to x
      TYPE(SD_InputType)                           :: u_interp    ! interpolated value of inputs 

      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = "" 

      ! Initialize interim vars
      k1=x
      x_tmp=x
      k2=x
      k3=x
      k4=x
      
      CALL MeshCopy ( SrcMesh  = u(1)%TPMesh         &
                    , DestMesh = u_interp%TPMesh     &
                    , CtrlCode = MESH_NEWCOPY        &
                    , ErrStat  = ErrStat             &
                    , ErrMess  = ErrMsg               )

      CALL MeshCopy ( SrcMesh  = u(1)%LMesh          &
                    , DestMesh = u_interp%LMesh      &
                    , CtrlCode = MESH_NEWCOPY        &
                    , ErrStat  = ErrStat             &
                    , ErrMess  = ErrMsg               )
      
     ! ALLOCATE(u_interp%UFL(p%uL), STAT=ErrStat)   !need to 
     ! IF ( ErrStat/= ErrID_None ) THEN
     !    ErrStat = ErrID_Fatal
     !    ErrMsg  = 'Error allocating input u_interp%UFL in  SD_RK4'
     !    RETURN
     ! END IF
     !u_interp%UFL=0
     
      ! interpolate u to find u_interp = u(t)
      CALL SD_Input_ExtrapInterp( u, utimes, u_interp, t, ErrStat, ErrMsg )

      ! find xdot at t
      CALL SubDyn_CalcContStateDeriv( t, u_interp, p, x, xd, z, OtherState, xdot, ErrStat, ErrMsg )

      k1%qm    = p%SDDeltaT * xdot%qm
      k1%qmdot = p%SDDeltaT * xdot%qmdot
  
      x_tmp%qm    = x%qm    + 0.5 * k1%qm
      x_tmp%qmdot = x%qmdot + 0.5 * k1%qmdot

      ! interpolate u to find u_interp = u(t + dt/2)
      CALL SD_Input_ExtrapInterp(u, utimes, u_interp, t+0.5*p%SDDeltaT, ErrStat, ErrMsg)

      ! find xdot at t + dt/2
      CALL SubDyn_CalcContStateDeriv( t + 0.5*p%SDDeltaT, u_interp, p, x_tmp, xd, z, OtherState, xdot, ErrStat, ErrMsg )

      k2%qm    = p%SDDeltaT * xdot%qm
      k2%qmdot = p%SDDeltaT * xdot%qmdot

      x_tmp%qm    = x%qm    + 0.5 * k2%qm
      x_tmp%qmdot = x%qmdot + 0.5 * k2%qmdot

      ! find xdot at t + dt/2
      CALL SubDyn_CalcContStateDeriv( t + 0.5*p%SDDeltaT, u_interp, p, x_tmp, xd, z, OtherState, xdot, ErrStat, ErrMsg )
     
      k3%qm    = p%SDDeltaT * xdot%qm
      k3%qmdot = p%SDDeltaT * xdot%qmdot

      x_tmp%qm    = x%qm    + k3%qm
      x_tmp%qmdot = x%qmdot + k3%qmdot

      ! interpolate u to find u_interp = u(t + dt)
      CALL SD_Input_ExtrapInterp(u, utimes, u_interp, t + p%SDDeltaT, ErrStat, ErrMsg)

      ! find xdot at t + dt
      CALL SubDyn_CalcContStateDeriv( t + p%SDDeltaT, u_interp, p, x_tmp, xd, z, OtherState, xdot, ErrStat, ErrMsg )

      k4%qm    = p%SDDeltaT * xdot%qm
      k4%qmdot = p%SDDeltaT * xdot%qmdot

      x%qm    = x%qm    +  ( k1%qm    + 2. * k2%qm    + 2. * k3%qm    + k4%qm    ) / 6.      
      x%qmdot = x%qmdot +  ( k1%qmdot + 2. * k2%qmdot + 2. * k3%qmdot + k4%qmdot ) / 6.      

      CALL SD_DestroyInput( u_interp, ErrStat, ErrMsg )
      
END SUBROUTINE SD_RK4
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE SD_AB4( t, n, u, utimes, p, x, xd, z, OtherState, ErrStat, ErrMsg )
!
! This subroutine implements the fourth-order Adams-Bashforth Method (RK4) for numerically integrating ordinary differential 
! equations:
!
!   Let f(t, x) = xdot denote the time (t) derivative of the continuous states (x). 
!
!   x(t+dt) = x(t)  + (dt / 24.) * ( 55.*f(t,x) - 59.*f(t-dt,x) + 37.*f(t-2.*dt,x) - 9.*f(t-3.*dt,x) )
!
!  See, e.g.,
!  http://en.wikipedia.org/wiki/Linear_multistep_method
!
!  or
!
!  K. E. Atkinson, "An Introduction to Numerical Analysis", 1989, John Wiley & Sons, Inc, Second Edition.
!
!..................................................................................................................................

      REAL(DbKi),                     INTENT(IN   )  :: t           ! Current simulation time in seconds
      INTEGER(IntKi),                 INTENT(IN   )  :: n           ! time step number
      TYPE(SD_InputType),             INTENT(INOUT)  :: u(:)        ! Inputs at t
      REAL(DbKi),                     INTENT(IN   )  :: utimes(:)   ! times of input
      TYPE(SD_ParameterType),         INTENT(IN   )  :: p           ! Parameters
      TYPE(SD_ContinuousStateType),   INTENT(INOUT)  :: x           ! Continuous states at t on input at t + dt on output
      TYPE(SD_DiscreteStateType),     INTENT(IN   )  :: xd          ! Discrete states at t
      TYPE(SD_ConstraintStateType),   INTENT(IN   )  :: z           ! Constraint states at t (possibly a guess)
      TYPE(SD_OtherStateType),        INTENT(INOUT)  :: OtherState  ! Other/optimization states
      INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None


      ! local variables
      TYPE(SD_ContinuousStateType) :: xdot       ! Continuous state derivs at t
      TYPE(SD_InputType)           :: u_interp
         

      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = "" 

      ! need xdot at t
      CALL SD_Input_ExtrapInterp(u, utimes, u_interp, t, ErrStat, ErrMsg)
      CALL SubDyn_CalcContStateDeriv( t, u_interp, p, x, xd, z, OtherState, xdot, ErrStat, ErrMsg )

      if (n .le. 2) then

         OtherState%n = n

         OtherState%xdot ( 3 - n ) = xdot

         CALL SD_RK4(t, n, u, utimes, p, x, xd, z, OtherState, ErrStat, ErrMsg )

      else

         if (OtherState%n .lt. n) then

            OtherState%n = n
            OtherState%xdot(4)    = OtherState%xdot(3)
            OtherState%xdot(3)    = OtherState%xdot(2)
            OtherState%xdot(2)    = OtherState%xdot(1)

         elseif (OtherState%n .gt. n) then
 
            ErrStat = ErrID_Fatal
            ErrMsg = ' Backing up in time is not supported with a multistep method '
            RETURN

         endif

         OtherState%xdot ( 1 )     = xdot  ! make sure this is most up to date

         x%qm    = x%qm    + (p%SDDeltaT / 24.) * ( 55.*OtherState%xdot(1)%qm - 59.*OtherState%xdot(2)%qm    + 37.*OtherState%xdot(3)%qm  &
                                       - 9. * OtherState%xdot(4)%qm )

         x%qmdot = x%qmdot + (p%SDDeltaT / 24.) * ( 55.*OtherState%xdot(1)%qmdot - 59.*OtherState%xdot(2)%qmdot  &
                                          + 37.*OtherState%xdot(3)%qmdot  - 9.*OtherState%xdot(4)%qmdot )

      endif


END SUBROUTINE SD_AB4
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE SD_ABM4( t, n, u, utimes, p, x, xd, z, OtherState, ErrStat, ErrMsg )
!
! This subroutine implements the fourth-order Adams-Bashforth-Moulton Method (RK4) for numerically integrating ordinary 
! differential equations:
!
!   Let f(t, x) = xdot denote the time (t) derivative of the continuous states (x). 
!
!   Adams-Bashforth Predictor:
!   x^p(t+dt) = x(t)  + (dt / 24.) * ( 55.*f(t,x) - 59.*f(t-dt,x) + 37.*f(t-2.*dt,x) - 9.*f(t-3.*dt,x) )
!
!   Adams-Moulton Corrector:
!   x(t+dt) = x(t)  + (dt / 24.) * ( 9.*f(t+dt,x^p) + 19.*f(t,x) - 5.*f(t-dt,x) + 1.*f(t-2.*dt,x) )
!
!  See, e.g.,
!  http://en.wikipedia.org/wiki/Linear_multistep_method
!
!  or
!
!  K. E. Atkinson, "An Introduction to Numerical Analysis", 1989, John Wiley & Sons, Inc, Second Edition.
!
!..................................................................................................................................

      REAL(DbKi),                     INTENT(IN   )  :: t           ! Current simulation time in seconds
      INTEGER(IntKi),                 INTENT(IN   )  :: n           ! time step number
      TYPE(SD_InputType),             INTENT(INOUT)  :: u(:)        ! Inputs at t
      REAL(DbKi),                     INTENT(IN   )  :: utimes(:)   ! times of input
      TYPE(SD_ParameterType),         INTENT(IN   )  :: p           ! Parameters
      TYPE(SD_ContinuousStateType),   INTENT(INOUT)  :: x           ! Continuous states at t on input at t + dt on output
      TYPE(SD_DiscreteStateType),     INTENT(IN   )  :: xd          ! Discrete states at t
      TYPE(SD_ConstraintStateType),   INTENT(IN   )  :: z           ! Constraint states at t (possibly a guess)
      TYPE(SD_OtherStateType),        INTENT(INOUT)  :: OtherState  ! Other/optimization states
      INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      ! local variables

      TYPE(SD_InputType)            :: u_interp        ! Continuous states at t
      TYPE(SD_ContinuousStateType)  :: x_pred          ! Continuous states at t
      TYPE(SD_ContinuousStateType)  :: xdot_pred       ! Continuous states at t

      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = "" 

      CALL SD_CopyContState(x, x_pred, 0, ErrStat, ErrMsg)

      CALL SD_AB4( t, n, u, utimes, p, x_pred, xd, z, OtherState, ErrStat, ErrMsg )

      if (n .gt. 2) then

         CALL SD_Input_ExtrapInterp(u, utimes, u_interp, t + p%SDDeltaT, ErrStat, ErrMsg)

         CALL SubDyn_CalcContStateDeriv(t + p%SDDeltaT, u_interp, p, x_pred, xd, z, OtherState, xdot_pred, ErrStat, ErrMsg )

         x%qm    = x%qm    + (p%SDDeltaT / 24.) * ( 9. * xdot_pred%qm +  19. * OtherState%xdot(1)%qm - 5. * OtherState%xdot(2)%qm &
                                          + 1. * OtherState%xdot(3)%qm )
   
         x%qmdot = x%qmdot + (p%SDDeltaT / 24.) * ( 9. * xdot_pred%qmdot + 19. * OtherState%xdot(1)%qmdot - 5. * OtherState%xdot(2)%qmdot &
                                          + 1. * OtherState%xdot(3)%qmdot )

      else

         x%qm    = x_pred%qm
         x%qmdot = x_pred%qmdot

      endif

END SUBROUTINE SD_ABM4

!----------------------------------------------------------------------------------------------------------------------------------
!SUBROUTINE SD_Input_ExtrapInterp( u, t, u_out, t_out, ErrStat, ErrMsg )
!!
!! This subroutine calculates a extrapolated (or interpolated) input u_out at time t_out, from previous/future time
!! values of u (which has values associated with times in t).  Order of the interpolation is given by the size of u
!!  
!!  expressions below based on either
!!
!!  f(t) = a 
!!  f(t) = a + b * t, or
!!  f(t) = a + b * t + c * t**2
!!
!!  where a, b and c are determined as the solution to
!!  f(t1) = u1, f(t2) = u2, f(t3) = u3  (as appropriate)
!!
!!..................................................................................................................................
!
!      TYPE(SD_InputType),      INTENT(IN   )  :: u(:)      ! Inputs at t1 > t2 > t3
!      REAL(DbKi),                INTENT(IN   )  :: t(:)      ! Times associated with the inputs 
!      TYPE(SD_InputType),      INTENT(  OUT)  :: u_out     ! Inputs at t1 > t2 > t3; initialize
!      REAL(DbKi),                INTENT(IN   )  :: t_out     ! time to be extrap/interp'd to                                   
!      INTEGER(IntKi),            INTENT(  OUT)  :: ErrStat   ! Error status of the operation
!      CHARACTER(*),              INTENT(  OUT)  :: ErrMsg    ! Error message if ErrStat /= ErrID_None                          
!
!      ! local variables
!
!      INTEGER(IntKi)        :: order     ! order of polynomial fit (max 2)
!   
!      REAL(ReKi)            :: a(SIZE(u(1)%UFL))        ! constant for extrapolaton/interpolation
!      REAL(ReKi)            :: b(SIZE(u(1)%UFL))        ! constant for extrabpolation/interpolation
!      REAL(ReKi)            :: c(SIZE(u(1)%UFL))     ! constant for extrabpolation/interpolation
!      
!      ! Initialize ErrStat
!
!      ErrStat = ErrID_None
!      ErrMsg  = "" 
!
!      if ( size(t) .ne. size(u)) then
!         ErrStat = ErrID_Fatal
!         ErrMsg = ' Error in SD_Input_ExtrapInterp: size(t) must equal size(u) '
!         RETURN
!      endif
!
!      if (size(u) .gt. 3) then
!         ErrStat = ErrID_Fatal
!         ErrMsg  = ' Error in SD_Input_ExtrapInterp: size(u) must be less than 4 '
!         RETURN
!      endif
!
!      !!How is it possible that we have to check this all the time?
!      u_out=u(1) !Initialize
!      !IF ( ErrStat/= 0 ) THEN
!      !   ErrStat = ErrID_Fatal
!      !   ErrMsg  = 'Error allocating states derivatives in SubDyn_CalcContStateDeriv'
!      !   RETURN
!      !END IF
!      order = SIZE(u) - 1
!
!      if (order .eq. 0) then
!
!         u_out%UFL = u(1)%UFL
!
!      elseif (order .eq. 1) then
!
!         IF ( EqualRealNos( t(1), t(2) ) ) THEN
!           ErrStat = ErrID_Fatal
!           ErrMsg  = ' Error in SD_Input_ExtrapInterp: t(1) must not equal t(2) to avoid a division-by-zero error.'
!           RETURN
!         END IF
!       
!         a = -(( t(2)*u(1)%UFL - t(1)*u(2)%UFL ) / ( t(1) - t(2) ) )
!         b = -((-u(1)%UFL + u(2)%UFL)/(t(1) - t(2)))
!
!         u_out%UFL = a + b * t_out
! 
!      elseif (order .eq. 2) then
!
!         IF ( EqualRealNos( t(1), t(2) ) ) THEN
!           ErrStat = ErrID_Fatal
!           ErrMsg  = ' Error in SD_Input_ExtrapInterp: t(1) must not equal t(2) to avoid a division-by-zero error.'
!           RETURN
!         END IF
!         IF ( EqualRealNos( t(2), t(3) ) ) THEN
!           ErrStat = ErrID_Fatal
!           ErrMsg  = ' Error in SD_Input_ExtrapInterp: t(2) must not equal t(3) to avoid a division-by-zero error.'
!           RETURN
!         END IF
!         IF ( EqualRealNos( t(1), t(3) ) ) THEN
!           ErrStat = ErrID_Fatal
!           ErrMsg  = ' Error in SD_Input_ExtrapInterp: t(1) must not equal t(3) to avoid a division-by-zero error.'
!           RETURN
!         END IF
! 
!         a =  (t(1)*t(3)*(-t(1) + t(3))*u(2)%UFL + t(2)**2*(t(3)*u(1)%UFL - t(1)*u(3)%UFL) &
!              + t(2)*(-(t(3)**2*u(1)%UFL) + t(1)**2*u(3)%UFL))/ & 
!              ((t(1) - t(2))*(t(1) - t(3))*(t(2) - t(3)))
!         b =  (t(3)**2*(u(1)%UFL - u(2)%UFL) + t(1)**2*(u(2)%UFL - u(3)%UFL) + t(2)**2*(-u(1)%UFL &
!              + u(3)%UFL))/((t(1) - t(2))*(t(1) - t(3))*(t(2) - t(3)))
!         c =  (t(3)*(-u(1)%UFL + u(2)%UFL) + t(2)*(u(1)%UFL - u(3)%UFL) + t(1)*(-u(2)%UFL + u(3)%UFL)) &
!              /((t(1) - t(2))*(t(1) - t(3))*(t(2) - t(3)))
!
!         u_out%UFL = a + b * t_out + c * t_out**2
!
!      endif
!
!END SUBROUTINE SD_Input_ExtrapInterp
!
End Module TempMod
