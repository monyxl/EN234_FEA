!
!    ABAQUS format UEL subroutine
!
!    This file is compatible with both EN234_FEA and ABAQUS/Standard
!
!    The example implements a standard fully integrated 3D linear elastic continuum element
!
!    The file also contains the following subrouines:
!          abq_UEL_2D_CH_integrationpoints           - defines integration points for 2D continuum elements
!          abq_UEL_2D_CH_shapefunctions              - defines shape functions for 2D continuum elements
!          abq_UEL_1D_CH_integrationpoints(n_points, n_nodes, xi, w)  = defines integration points for 1D line integral
!=========================== ABAQUS format user element subroutine ===================

      SUBROUTINE UEL_CH(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     1     PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     2     KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,NPREDF,
     3     LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,PERIOD)
    !
      INCLUDE 'ABA_PARAM.INC'
    !
    !
      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),
     1   SVARS(*),ENERGY(8),COORDS(MCRD,NNODE),U(NDOFEL),
     2   DU(MLVARX,*),V(NDOFEL),A(NDOFEL),TIME(2),PARAMS(*),
     3   JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),
     4   PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)

    !
    !       Variables that must be computed in this routine
    !       RHS(i)                     Right hand side vector.  In EN234_FEA the dimensions are always RHS(MLVARX,1)
    !       AMATRX(i,j)                Stiffness matrix d RHS(i)/ d DU(j)
    !       SVARS(1:NSVARS)            Element state variables.  Must be updated in this routine
    !       ENERGY(1:8)
    !                                  Energy(1) Kinetic Energy
    !                                  Energy(2) Elastic Strain Energy
    !                                  Energy(3) Creep Dissipation
    !                                  Energy(4) Plastic Dissipation
    !                                  Energy(5) Viscous Dissipation
    !                                  Energy(6) Artificial strain energy
    !                                  Energy(7) Electrostatic energy
    !                                  Energy(8) Incremental work done by loads applied to the element
    !       PNEWDT                     Allows user to control ABAQUS time increments.
    !                                  If PNEWDT<1 then time step is abandoned and computation is restarted with
    !                                  a time increment equal to PNEWDT*DTIME
    !                                  If PNEWDT>1 ABAQUS may increase the time increment by a factor PNEWDT
    !
    !       Variables provided for information
    !       NDOFEL                     Total # DOF for the element
    !       NRHS                       Dimension variable
    !       NSVARS                     Total # element state variables
    !       PROPS(1:NPROPS)            User-specified properties of the element
    !       NPROPS                     No. properties
    !       JPROPS(1:NJPROPS)          Integer valued user specified properties for the element
    !       NJPROPS                    No. integer valued properties
    !       COORDS(i,N)                ith coordinate of Nth node on element
    !       MCRD                       Maximum of (# coords,minimum of (3,#DOF)) on any node
    !       U                          Vector of DOF at the end of the increment
    !       DU                         Vector of DOF increments
    !       V                          Vector of velocities (defined only for implicit dynamics)
    !       A                          Vector of accelerations (defined only for implicit dynamics)
    !       JTYPE                      Integer identifying element type (the number n in the Un specification in the input file)
    !       TIME(1:2)                  TIME(1)   Current value of step time
    !                                  TIME(2)   Total time
    !       DTIME                      Time increment
    !       KSTEP                      Current step number (always 1 in EN234_FEA)
    !       KINC                       Increment number
    !       JELEM                      User assigned element number in ABAQUS (internally assigned in EN234_FEA)
    !       PARAMS(1:3)                Time increment parameters alpha, beta, gamma for implicit dynamics
    !       NDLOAD                     Number of user-defined distributed loads defined for this element
    !       JDLTYP(1:NDLOAD)           Integers n defining distributed load types defined as Un or (if negative) UnNU in input file
    !       ADLMAG(1:NDLOAD)           Distributed load magnitudes
    !       DDLMAG(1:NDLOAD)           Increment in distributed load magnitudes
    !       PREDEF(1:2,1:NPREDF,1:NNODE)   Predefined fields.
    !       PREDEF(1,...)              Value of predefined field
    !       PREDEF(2,...)              Increment in predefined field
    !       PREDEF(1:2,1,k)            Value of temperature/temperature increment at kth node
    !       PREDEF(1:2,2:NPREDF,k)     Value of user defined field/field increment at kth node (not used in EN234FEA)
    !       NPREDF                     Number of predefined fields (1 for en234FEA)
    !       LFLAGS                     Control variable
    !       LFLAGS(1)                  Defines procedure type
    !       LFLAGS(2)                  0 => small displacement analysis  1 => Large displacement (NLGEOM option)
    !       LFLAGS(3)                   1 => Subroutine must return both RHS and AMATRX (always true in EN234FEA)
    !                                   2 => Subroutine must return stiffness AMATRX = -dF/du
    !                                   3 => Subroutine must return daming matrix AMATRX = -dF/dudot
    !                                   4 => Subroutine must return mass matrix AMATRX = -dF/duddot
    !                                   5 => Define the RHS only
    !                                   6 => Define the mass matrix for the initial acceleration calculation
    !                                   100 => Define perturbation quantities for output
    !       LFLAGS(4)                   0 => General step   1 => linear perturbation step
    !       LFLAGS(5)                   0 => current approximation to solution based on Newton correction; 1 => based on extrapolation
    !       MLVARX                      Dimension variable (equal to NDOFEL in EN234FEA)
    !       PERIOD                      Time period of the current step
    !
    ! Local Variables
      integer      :: i,j,n_points,kint, nfacenodes, ipoin, ksize
      integer      :: n_points_2, n_points_tot,NNODE_2, NNODE_tot
      integer      :: face_node_list(3)                       ! List of nodes on an element face
    !
      double precision  ::  xi(2,9)                           ! Area integration points for disp
      double precision  ::  xi2(2,9)                          ! Area integration points for concentration
      double precision  ::  w(9)                              ! Area integration weights
      double precision  ::  w2(9)                             ! Area integration weights
      double precision  ::  N(9)                              ! 2D shape functions for disp
      double precision  ::  N2(9)                             ! 2D shape functions for concentration
      double precision  ::  dNdxi(9,2)                        ! 2D shape function derivatives for disp
      double precision  ::  dN2dxi(9,2)                       ! 2D shape function derivatives for concentration
      double precision  ::  dNdx(9,2)                         ! Spatial derivatives for disp
      double precision  ::  dN2dx(9,2)                      ! Spatial derivatives for concentration
      double precision  ::  dxdxi(2,2)                        ! Derivative of spatial coords wrt normalized coords
      double precision  ::  dxdxi2(2,2)                       ! Derivative of spatial coords wrt normalized coords

    !   Variables below are for computing integrals over element faces
      double precision  ::  face_coords(2,3)                  ! Coords of nodes on an element face
      double precision  ::  xi1(6)                            ! 1D integration points
      double precision  ::  w1(6)                             ! Integration weights
      double precision  ::  N1(3)                             ! 1D shape functions
      double precision  ::  dN1dxi(3)                         ! 1D shape function derivatives
      double precision  ::  norm(2)                           ! Normal to an element face
      double precision  ::  dxdxi1(2)                         ! Derivative of 1D spatial coord wrt normalized areal coord
    !
      double precision  ::  strain(4)                         ! Strain vector contains [e11, e22, e33, 2e12]
      double precision  ::  stress(4)                         ! Stress vector contains [s11, s22, s33, s12]
      double precision  ::  DEL(4,4)                          ! D matrix inside element: DEL
      double precision  ::  B(9,24)                           ! strain = B*(dof_total)
      double precision  ::  ktemp(22,22)                      ! Temporary stiffness (for incompatible mode elements)
      double precision  ::  rhs_temp(22)                      ! Temporary RHS vector (for incompatible mode elements)
      double precision  ::  kuu(18,18)                        ! Upper diagonal stiffness
      double precision  ::  kaa(4,4),kaainv(4,4)              ! Lower diagonal stiffness
      double precision  ::  kau(4,18)                         ! Lower quadrant of stiffness
      double precision  ::  kua(18,4)                         ! Upper quadrant of stiffness
      double precision  ::  alpha(4)                          ! Internal DOF for incompatible mode element
      double precision  ::  Utemp(12)                         ! Internal DOF for incompatible mode element
      double precision  ::  dxidx(2,2),  det                  ! Jacobian inverse and determinant
      double precision  ::  dxidx2(2,2), det2                 ! Jacobian inverse and determinant
      double precision  ::  E, xnu, D44, D11, D12             ! Material properties
      double precision  ::  Omega, Wgibbs, Kappa, Diff, Theta      ! Material properties
      double precision  ::  D(9,9)                            ! stress=D*(strain)  (NOTE FACTOR OF 2 in shear strain)
      double precision  ::  q(9)                              ! q vector
      double precision  ::  F, dFdc                           ! Free energy and its derivatives
      double precision  ::  sol(9), dsol(9)                   ! generalized disp and disp inc
      double precision  ::  c                                 ! concentration
      double precision  ::  diag(4)                           ! diagonal vec [1 1 1 0]
      

    !
    !     Example ABAQUS UEL implementing 2D linear elastic elements
    !     Includes option for incompatible mode elements
    !     El props are:

    !     PROPS(1)         Young's modulus
    !     PROPS(2)         Poisson's ratio


      if (NNODE == 3) n_points = 1                ! Linear triangle
      if (NNODE == 4) n_points = 4                ! Linear rectangle
      if (NNODE == 6) n_points = 4                ! Quadratic triangle
      if (NNODE == 8) n_points = 9                ! Serendipity rectangle
      if (NNODE == 9) n_points = 9                ! Quadratic rect
      
      NNODE_2 = 4                                 ! int pts for concentration
      n_points = 4                              ! # OF int pts for concentration
      NNODE_tot = NNODE + NNODE_2
      
    ! Write your code for a 2D element below
      
      call abq_UEL_2D_CH_integrationpoints(n_points, NNODE, xi, w)
      call abq_UEL_2D_CH_integrationpoints(n_points, NNODE_2, xi2, w2)
      
      !if (MLVARX<2*NNODE) then
      !  write(6,*) ' Error in abaqus UEL '
      !  write(6,*) ' Variable MLVARX must exceed 2*NNODE'
      !  write(6,*) ' MLVARX = ',MLVARX,' NNODE = ',NNODE
      !  stop
      !endif
      
      RHS(1:MLVARX,1) = 0.d0
      AMATRX(1:NDOFEL,1:NDOFEL) = 0.d0
      
      DEL = 0.d0
      D = 0.d0
      
      E = PROPS(1)
      xnu = PROPS(2)
      Omega = PROPS(3)
      Wgibbs = PROPS(4)
      Kappa = PROPS(5)
      Diff = PROPS(6)
      Theta = PROPS(7)
      
      D44 = 0.5D0*E/(1+xnu)
      D11 = (1.D0-xnu)*E/( (1.D0+xnu)*(1.D0-2.D0*xnu) )
      D12 = xnu*E/( (1.D0+xnu)*(1.D0-2.D0*xnu) )
      
      DEL(1:3,1:3) = D12
      DEL(1,1) = D11
      DEL(2,2) = D11
      DEL(4,4) = D44
      
      D(1:2,1:2) = DEL(1:2,1:2)
      D(3,3) = DEL(4,4)
      D(1,5) = -Omega*(DEL(1,1)+DEL(1,2)+DEL(1,3))/3.d0
      D(2,5) = -Omega*(DEL(2,1)+DEL(2,2)+DEL(2,3))/3.d0
      D(4,1) = -Omega*(DEL(1,1)+DEL(1,2)+DEL(1,3))/3.d0
      D(4,2) = -Omega*(DEL(2,1)+DEL(2,2)+DEL(2,3))/3.d0
      D(4,4) = 1.d0
      D(4,5) = 0.d0 ! WILL CALCULATE LATER
      D(5,5) = 1.d0/DTIME
      D(6,8) = -Kappa
      D(7,9) = -Kappa
      D(8,6) = Theta*Diff
      D(9,7) = Theta*Diff
      
      
      ENERGY(1:8) = 0.d0
      det = 0.d0
      det2 = 0.d0
      
      !     --  Loop over integration points
      do kint = 1, n_points
      !     --  Loop over disp integration points
        call abq_UEL_2D_CH_shapefunctions(xi(1:2,kint),NNODE,N,dNdxi)
        dxdxi = matmul(coords(1:2,1:NNODE),dNdxi(1:NNODE,1:2))
        call abq_inverse_LU_CH(dxdxi,dxidx,2)
        det = dxdxi(1,1)*dxdxi(2,2)-dxdxi(1,2)*dxdxi(2,1)
        IF (det==0.d0) THEN
          write(6,*) ' Error in disp'
          write(6,*) ' A 2x2 matrix has a zero determinant'
          stop
        endif
        dNdx(1:NNODE,1:2) = matmul(dNdxi(1:NNODE,1:2),dxidx)
      !     --  Loop over concentration integration points
        call abq_UEL_2D_CH_shapefunctions(xi2(1:2,kint),NNODE_2,N2,
     1                                                       dN2dxi)
        dxdxi2 = matmul(coords(1:2,1:NNODE_2),dN2dxi(1:NNODE_2,1:2))
        call abq_inverse_LU_CH(dxdxi2,dxidx2,2)
        det2 = dxdxi2(1,1)*dxdxi2(2,2)-dxdxi2(1,2)*dxdxi2(2,1)
        IF (det2==0.d0) THEN
          write(6,*) ' Error in concentration'
          write(6,*) ' A 2x2 matrix has a zero determinant'
          stop
        endif
        dN2dx(1:NNODE_2,1:2) = matmul(dN2dxi(1:NNODE_2,1:2),dxidx2)
        !dN2dx(1:NNODE_2,1:2) = matmul(dN2dxi(1:NNODE_2,1:2),dxidx)
        
        B = 0.d0
        !B(1,1:2*NNODE-1:2) = dNdx(1:NNODE,1)
        !B(2,2:2*NNODE:2) = dNdx(1:NNODE,2)
        !B(4,1:2*NNODE-1:2) = dNdx(1:NNODE,2)
        !B(4,2:2*NNODE:2) = dNdx(1:NNODE,1)
        
        ! First line
        B(1,1:4*NNODE_2-3:4) = dNdx(1:NNODE_2,1)
        B(1,4*NNODE_2+1:2*NNODE_tot-1:2) = dNdx(NNODE_2+1:NNODE,1)
        !Second line
        B(2,2:4*NNODE_2-2:4) = dNdx(1:NNODE_2,2)
        B(2,4*NNODE_2+2:2*NNODE_tot:2) = dNdx(NNODE_2+1:NNODE,2)
        !Third line
        B(3,1:4*NNODE_2-3:4) = dNdx(1:NNODE_2,2)
        B(3,4*NNODE_2+1:2*NNODE_tot-1:2) = dNdx(NNODE_2+1:NNODE,2)
        B(3,2:4*NNODE_2-2:4) = dNdx(1:NNODE_2,1)
        B(3,4*NNODE_2+2:2*NNODE_tot:2) = dNdx(NNODE_2+1:NNODE,1)
        ! Fourth line
        B(4,3:4*NNODE_2-1:4) = N2(1:NNODE_2)
        ! Fifth line
        B(5,4:4*NNODE_2:4) = N2(1:NNODE_2)
        ! Sixth line
        B(6,3:4*NNODE_2-1:4) = dN2dx(1:NNODE_2,1)
        ! Seventh line
        B(7,3:4*NNODE_2-1:4) = dN2dx(1:NNODE_2,2)
        ! Eighth line
        B(8,4:4*NNODE_2:4) = dN2dx(1:NNODE_2,1)
        ! Nineth line
        B(9,4:4*NNODE_2:4) = dN2dx(1:NNODE_2,2)
        
        sol = 0.d0
        dsol = 0.d0
        sol(1:9) = matmul(B(1:9,1:24),U(1:24))        ! SOL(5) is concentration
        dsol(1:9) = matmul(B(1:9,1:24),DU(1:24,1))    ! dSOL(5) is concentration inc
        
        c = sol(5)
        
        ! CALCULATE strain and stress
        strain = 0.d0
        stress = 0.d0
        diag = 0.d0
        diag(1:3) = 1.d0
        strain(1:2) = sol(1:2)
        strain(4) = sol(3)
        stress(1:4) = matmul(DEL(1:4,1:4)
     1        ,(strain(1:4)-diag(1:4)*Omega*c/3.d0))
        
        ! CALCULATE dFdc
    !   dFdc = Wgibbs*(12.d0*(c+dsol(5))**2.d0
    !1         -12.d0*(c+dsol(5))+2.d0)
        dFdc = Wgibbs*(12.d0*(c)**2.d0
     1         -12.d0*(c)+2.d0)
        ! CALCULATE F(C+DC)
    !   F = 2.d0*Wgibbs*(c+dsol(5))*(c+dsol(5)-1.d0)
    !1      *(2.d0*(c+dsol(5))-1.d0)
        
        F = 2.d0*Wgibbs*(c)*(c-1.d0)
     1      *(2.d0*(c)-1.d0)
        
        ! CALCULATE q vector
        q = 0.d0
        q(1:2) = stress(1:2)
        q(3) = stress(4)
        q(4) = sol(4)-F-Omega*
     1              sum(stress(1:3))/3.d0
        q(5) = dsol(5)/DTIME
        !q(6) = -Kappa*(sol(8)+dsol(8))
        !q(7) = -Kappa*(sol(9)+dsol(9))
        !q(8) = Diff*(sol(6)+Theta*dsol(6))
        !q(9) = Diff*(sol(7)+Theta*dsol(7))
        q(6) = -Kappa*(sol(8))
        q(7) = -Kappa*(sol(9))
        q(8) = Diff*(sol(6)+(Theta-1.d0)*dsol(6))
        q(9) = Diff*(sol(7)+(Theta-1.d0)*dsol(7))
        ! CALCULATE  D(4,5)
        D(4,5) = -dFdc + (Omega**2.d0)*(sum(DEL(1:3,1:3)))/9.d0
        
        
        ! CALCULATE  RHS
        RHS(1:2*NNODE_tot,1) = RHS(1:2*NNODE_tot,1)
     1   - matmul(transpose(B(1:9,1:2*NNODE_tot)),q(1:9))*
     2                                          w(kint)*det
        
        AMATRX(1:2*NNODE_tot,1:2*NNODE_tot) = 
     1  AMATRX(1:2*NNODE_tot,1:2*NNODE_tot)
     2    + matmul(transpose(B(1:9,1:2*NNODE_tot)),
     3      matmul(D(1:9,1:9),B(1:9,1:2*NNODE_tot)))*w(kint)*det
        
        ENERGY(2) = ENERGY(2)
     1   + 0.5D0*dot_product(stress,strain)*w(kint)*det           ! Store the elastic strain energy
        
        !stress(4)=1.d0
        if (NSVARS>=n_points*4) then   ! Store stress at each integration point (if space was allocated to do so)
            SVARS(4*kint-3:4*kint) = stress(1:4)
        endif
        
      end do
      
      PNEWDT = 1.d0          ! This leaves the timestep unchanged (ABAQUS will use its own algorithm to determine DTIME)
    !
    !   Apply distributed loads
    !
    !   Distributed loads are specified in the input file using the Un option in the input file.
    !   n specifies the face number, following the ABAQUS convention
    !
    !
!       do j = 1,NDLOAD
!          
!        call abq_facenodes_2D_CH(NNODE,iabs(JDLTYP(j,1)),
!     1                                     face_node_list,nfacenodes)
!        
!        do i = 1,nfacenodes
!            face_coords(1:2,i) = coords(1:2,face_node_list(i))
!        end do
!     
!        if (nfacenodes == 3) n_points = 3
!        if (nfacenodes == 6) n_points = 4
!        if (nfacenodes == 4) n_points = 4
!        if (nfacenodes == 8) n_points = 9
!        call abq_UEL_1D_CH_integrationpoints(n_points, nfacenodes, xi1, w)
!        
!        do kint = 1,n_points
!            dxdxi1 = matmul(face_coords(1:2,1:nfacenodes),
!     1                           dN1dxi(1:nfacenodes))
!            !norm(1)=(dxdxi2(2,1)*dxdxi2(3,2))-(dxdxi2(2,2)*dxdxi2(3,1))
!            !norm(2)=(dxdxi2(1,1)*dxdxi2(3,2))-(dxdxi2(1,2)*dxdxi2(3,1))
!            !norm(3)=(dxdxi2(1,1)*dxdxi2(2,2))-(dxdxi2(1,2)*dxdxi2(2,1))
!     
!            do i = 1,nfacenodes
!                ipoin = 2*face_node_list(i)-1
!                RHS(ipoin:ipoin+1,1) = RHS(ipoin:ipoin+1,1)
!     1              - N1(1:nfacenodes)*adlmag(j,1)*dxdxi1(1:2)**w(kint)
!            end do
!        end do
!      end do
      
      return
      
      END SUBROUTINE UEL_CH



      subroutine abq_UEL_2D_CH_integrationpoints(n_points, n_nodes, 
     1                                                       xi, w)

      implicit none
      integer, intent(in) :: n_points
      integer, intent(in) :: n_nodes

      double precision, intent(out) :: xi(2,*)
      double precision, intent(out) :: w(*)

      integer :: i,j,k,n

      double precision :: cn,w1,w2,w11,w12,w22

    !         Defines integration points and weights for 2D continuum elements

      if ( n_points==1 ) then
        if ( n_nodes==4 .or. n_nodes==9 ) then    !     ---   4 or 9 noded quad
            xi(1, 1) = 0.D0
            xi(2, 1) = 0.D0
            w(1) = 4.D0
        else if ( n_nodes==3 .or. n_nodes==6 ) then !     ---   3 or 6 noded triangle
            xi(1, 1) = 1.D0/3.D0
            xi(2, 1) = 1.D0/3.D0
            w(1) = 1.D0/2.D0
        end if
      else if ( n_points==3 ) then
        xi(1, 1) = 0.5D0
        xi(2, 1) = 0.5D0
        w(1) = 1.D0/6.D0
        xi(1, 2) = 0.D0
        xi(2, 2) = 0.5D0
        w(2) = w(1)
        xi(1, 3) = 0.5D0
        xi(2, 3) = 0.D0
        w(3) = w(1)
      else if ( n_points==4 ) then
        if ( n_nodes==4 .or. n_nodes==8 .or. n_nodes==9 ) then
            !     2X2 GAUSS INTEGRATION POINTS FOR QUADRILATERAL
            !     43
            !     12
            cn = 0.5773502691896260D0
            xi(1, 1) = -cn
            xi(1, 2) = cn
            xi(1, 3) = cn
            xi(1, 4) = -cn
            xi(2, 1) = -cn
            xi(2, 2) = -cn
            xi(2, 3) = cn
            xi(2, 4) = cn
            w(1) = 1.D0
            w(2) = 1.D0
            w(3) = 1.D0
            w(4) = 1.D0
        else if ( n_nodes==3 .or. n_nodes==6 ) then
            !     xi integration points for triangle
            xi(1, 1) = 1.D0/3.D0
            xi(2, 1) = xi(1, 1)
            w(1) = -27.D0/96.D0
            xi(1, 2) = 0.6D0
            xi(2, 2) = 0.2D0
            w(2) = 25.D0/96.D0
            xi(1, 3) = 0.2D0
            xi(2, 3) = 0.6D0
            w(3) = w(2)
            xi(1, 4) = 0.2D0
            xi(2, 4) = 0.2D0
            w(4) = w(2)
        end if

      else if ( n_points==7 ) then
        ! Quintic integration for triangle
        xi(1,1) = 1.d0/3.d0
        xi(2,1) = xi(1,1)
        w(1) = 0.1125d0
        xi(1,2) = 0.0597158717d0
        xi(2,2) = 0.4701420641d0
        w(2) = 0.0661970763d0
        xi(1,3) = xi(2,2)
        xi(2,3) = xi(1,2)
        w(3) = w(2)
        xi(1,4) = xi(2,2)
        xi(2,4) = xi(2,2)
        w(4) = w(2)
        xi(1,5) = 0.7974269853d0
        xi(2,5) = 0.1012865073d0
        w(5) = 0.0629695902d0
        xi(1,6) = xi(2,5)
        xi(2,6) = xi(1,5)
        w(6) = w(5)
        xi(1,7) = xi(2,5)
        xi(2,7) = xi(2,5)
        w(7) = w(5)
      else if ( n_points==9 ) then
        !     3X3 GAUSS INTEGRATION POINTS
        !     789
        !     456
        !     123
        cn = 0.7745966692414830D0
        xi(1, 1) = -cn
        xi(1, 2) = 0.D0
        xi(1, 3) = cn
        xi(1, 4) = -cn
        xi(1, 5) = 0.D0
        xi(1, 6) = cn
        xi(1, 7) = -cn
        xi(1, 8) = 0.D0
        xi(1, 9) = cn
        xi(2, 1) = -cn
        xi(2, 2) = -cn
        xi(2, 3) = -cn
        xi(2, 4) = 0.D0
        xi(2, 5) = 0.D0
        xi(2, 6) = 0.D0
        xi(2, 7) = cn
        xi(2, 8) = cn
        xi(2, 9) = cn
        w1 = 0.5555555555555560D0
        w2 = 0.8888888888888890D0
        w11 = w1*w1
        w12 = w1*w2
        w22 = w2*w2
        w(1) = w11
        w(2) = w12
        w(3) = w11
        w(4) = w12
        w(5) = w22
        w(6) = w12
        w(7) = w11
        w(8) = w12
        w(9) = w11
      end if

      return

      end subroutine abq_UEL_2D_CH_integrationpoints




      subroutine abq_UEL_2D_CH_shapefunctions(xi,n_nodes,f,df)

      implicit none
      integer, intent(in) :: n_nodes

      double precision, intent(in) :: xi(2)
      double precision, intent(out) :: f(*)
      double precision, intent(out) :: df(9,2)
      double precision g1, g2, g3, dg1, dg2, dg3
      double precision h1, h2, h3, dh1, dh2, dh3
      double precision z,dzdp, dzdq

            if ( n_nodes==3 ) then        !     SHAPE FUNCTIONS FOR 3 NODED TRIANGLE
                f(1) = xi(1)
                f(2) = xi(2)
                f(3) = 1.D0 - xi(1) - xi(2)
                df(1, 1) = 1.D0
                df(1, 2) = 0.D0
                df(2, 1) = 0.D0
                df(2, 2) = 1.D0
                df(3, 1) = -1.D0
                df(3, 2) = -1.D0
            else if ( n_nodes==4 ) then
                !     SHAPE FUNCTIONS FOR 4 NODED QUADRILATERAL
                !     43
                !     12
                g1 = 0.5D0*(1.D0 - xi(1))
                g2 = 0.5D0*(1.D0 + xi(1))
                h1 = 0.5D0*(1.D0 - xi(2))
                h2 = 0.5D0*(1.D0 + xi(2))
                f(1) = g1*h1
                f(2) = g2*h1
                f(3) = g2*h2
                f(4) = g1*h2
                dg1 = -0.5D0
                dg2 = 0.5D0
                dh1 = -0.5D0
                dh2 = 0.5D0
                df(1, 1) = dg1*h1
                df(2, 1) = dg2*h1
                df(3, 1) = dg2*h2
                df(4, 1) = dg1*h2
                df(1, 2) = g1*dh1
                df(2, 2) = g2*dh1
                df(3, 2) = g2*dh2
                df(4, 2) = g1*dh2

            else if ( n_nodes==6 ) then

                !     SHAPE FUNCTIONS FOR 6 NODED TRIANGLE
                !          3

                !       6      5

                !     1    4     2

                !     P = L1
                !     Q = L2
                !     Z = 1 - P - Q = L3

                z = 1.D0 - xi(1) - xi(2)
                f(1) = (2.D0*xi(1) - 1.D0)*xi(1)
                f(2) = (2.D0*xi(2) - 1.D0)*xi(2)
                f(3) = (2.D0*z - 1.D0)*z
                f(4) = 4.D0*xi(1)*xi(2)
                f(5) = 4.D0*xi(2)*z
                f(6) = 4.D0*xi(1)*z
                dzdp = -1.D0
                dzdq = -1.D0
                df(1, 1) = 4.D0*xi(1) - 1.D0
                df(2, 1) = 0.D0
                df(3, 1) = 4.D0*z*dzdp - dzdp
                df(4, 1) = 4.D0*xi(2)
                df(5, 1) = 4.D0*xi(2)*dzdp
                df(6, 1) = 4.D0*z + 4.D0*xi(1)*dzdp
                df(1, 2) = 0.D0
                df(2, 2) = 4.D0*xi(2) - 1.D0
                df(3, 2) = 4.D0*z*dzdq - dzdq
                df(4, 2) = 4.D0*xi(1)
                df(5, 2) = 4.D0*z + 4.D0*xi(2)*dzdq
                df(6, 2) = 4.D0*xi(1)*dzdq

            else if ( n_nodes==8 ) then
                !     SHAPE FUNCTIONS FOR 8 NODED SERENDIPITY ELEMENT
                 f(1) = -0.25*(1.-xi(1))*(1.-xi(2))*(1.+xi(1)+xi(2));
                 f(2) = 0.25*(1.+xi(1))*(1.-xi(2))*(xi(1)-xi(2)-1.);
                 f(3) = 0.25*(1.+xi(1))*(1.+xi(2))*(xi(1)+xi(2)-1.);
                 f(4) = 0.25*(1.-xi(1))*(1.+xi(2))*(xi(2)-xi(1)-1.);
                 f(5) = 0.5*(1.-xi(1)*xi(1))*(1.-xi(2));
                 f(6) = 0.5*(1.+xi(1))*(1.-xi(2)*xi(2));
                 f(7) = 0.5*(1.-xi(1)*xi(1))*(1.+xi(2));
                 f(8) = 0.5*(1.-xi(1))*(1.-xi(2)*xi(2));
                 df(1,1) = 0.25*(1.-xi(2))*(2.*xi(1)+xi(2));
                 df(1,2) = 0.25*(1.-xi(1))*(xi(1)+2.*xi(2));
                 df(2,1) = 0.25*(1.-xi(2))*(2.*xi(1)-xi(2));
                 df(2,2) = 0.25*(1.+xi(1))*(2.*xi(2)-xi(1));
                 df(3,1) = 0.25*(1.+xi(2))*(2.*xi(1)+xi(2));
                 df(3,2) = 0.25*(1.+xi(1))*(2.*xi(2)+xi(1));
                 df(4,1) = 0.25*(1.+xi(2))*(2.*xi(1)-xi(2));
                 df(4,2) = 0.25*(1.-xi(1))*(2.*xi(2)-xi(1));
                 df(5,1) = -xi(1)*(1.-xi(2));
                 df(5,2) = -0.5*(1.-xi(1)*xi(1));
                 df(6,1) = 0.5*(1.-xi(2)*xi(2));
                 df(6,2) = -(1.+xi(1))*xi(2);
                 df(7,1) = -xi(1)*(1.+xi(2));
                 df(7,2) = 0.5*(1.-xi(1)*xi(1));
                 df(8,1) = -0.5*(1.-xi(2)*xi(2));
                 df(8,2) = -(1.-xi(1))*xi(2);
            else if ( n_nodes==9 ) then
                !     SHAPE FUNCTIONS FOR 9 NODED LAGRANGIAN ELEMENT
                !     789
                !     456
                !     123
                g1 = -.5D0*xi(1)*(1.D0 - xi(1))
                g2 = (1.D0 - xi(1))*(1.D0 + xi(1))
                g3 = .5D0*xi(1)*(1.D0 + xi(1))
                h1 = -.5D0*xi(2)*(1.D0 - xi(2))
                h2 = (1.D0 - xi(2))*(1.D0 + xi(2))
                h3 = .5D0*xi(2)*(1.D0 + xi(2))
                dg1 = xi(1) - 0.5d0
                dg2 = -2.d0*xi(1)
                dg3 = xi(1) + 0.5d0
                dh1 = xi(2)-0.5d0
                dh2 = -2.d0*xi(2)
                dh3 = xi(2) + 0.5d0
                f(1) = g1*h1
                f(2) = g2*h1
                f(3) = g3*h1
                f(4) = g1*h2
                f(5) = g2*h2
                f(6) = g3*h2
                f(7) = g1*h3
                f(8) = g2*h3
                f(9) = g3*h3
                df(1,1) = dg1*h1
                df(1,2) = g1*dh1
                df(2,1) = dg2*h1
                df(2,2) = g2*dh1
                df(3,1) = dg3*h1
                df(3,2) = g3*dh1
                df(4,1) = dg1*h2
                df(4,2) = g1*dh2
                df(5,1) = dg2*h2
                df(5,2) = g2*dh2
                df(6,1) = dg3*h2
                df(6,2) = g3*dh2
                df(7,1) = dg1*h3
                df(7,2) = g1*dh3
                df(8,1) = dg2*h3
                df(8,2) = g2*dh3
                df(9,1) = dg3*h3
                df(9,2) = g3*dh3
            end if

      end subroutine abq_UEL_2D_CH_shapefunctions


      subroutine abq_UEL_1D_CH_integrationpoints(n_points, n_nodes, 
     1                                                        xi, w)


      implicit none
      integer, intent(in) :: n_points
      integer, intent(in) :: n_nodes

      double precision, intent(out) :: xi(*)
      double precision, intent(out) :: w(*)

      integer :: i,j,k,n

      double precision x1D(4), w1D(4)



      select case ( n_points )
        case (2)
            xi(1) = .5773502691896257D+00
            xi(2) = -.5773502691896257D+00
            w(1) = .1000000000000000D+01
            w(2) = .1000000000000000D+01
            return
        case (3)
            xi(1) = 0.7745966692414834D+00
            xi(2) = .0000000000000000D+00
            xi(3) = -.7745966692414834D+00
            w(1) = .5555555555555556D+00
            w(2) = .8888888888888888D+00
            w(3) = .5555555555555556D+00
            return
        case (4)
            xi(1) = .8611363115940526D+00
            xi(2) = .3399810435848563D+00
            xi(3) = -.3399810435848563D+00
            xi(4) = -.8611363115940526D+00
            w(1) = .3478548451374538D+00
            w(2) = .6521451548625461D+00
            w(3) = .6521451548625461D+00
            w(4) = .3478548451374538D+00
            return
        case (5)
            xi(1) = .9061798459386640D+00
            xi(2) = .5384693101056831D+00
            xi(3) = .0000000000000000D+00
            xi(4) = -.5384693101056831D+00
            xi(5) = -.9061798459386640D+00
            w(1) = .2369268850561891D+00
            w(2) = .4786286704993665D+00
            w(3) = .5688888888888889D+00
            w(4) = .4786286704993665D+00
            w(5) = .2369268850561891D+00
            return
        case (6)
            xi(1) = .9324695142031521D+00
            xi(2) = .6612093864662645D+00
            xi(3) = .2386191860831969D+00
            xi(4) = -.2386191860831969D+00
            xi(5) = -.6612093864662645D+00
            xi(6) = -.9324695142031521D+00
            w(1) = .1713244923791703D+00
            w(2) = .3607615730481386D+00
            w(3) = .4679139345726910D+00
            w(4) = .4679139345726910D+00
            w(5) = .3607615730481386D+00
            w(6) = .1713244923791703D+00
            return
        case DEFAULT
            write(6,*)'Error in subroutine 
     1                    abq_UEL_1D_CH_integrationpoints'
            write(6,*) ' Invalid number of integration points for 1D'
            write(6,*) ' n_points must be between 1 and 6'
            stop
      end select







      end subroutine abq_UEL_1D_CH_integrationpoints



      subroutine abq_facenodes_2D_CH(nelnodes,face,list,nfacenodes)

      implicit none

      integer, intent (in)      :: nelnodes
      integer, intent (in)      :: face
      integer, intent (out)     :: list(*)
      integer, intent (out)     :: nfacenodes
    !
    !        Subroutine to return list of nodes on an element face for standard 2D solid elements
    !
      integer :: i3(3)
      integer :: i4(4)

      i3(1:3) = [2,3,1]
      i4(1:4) = [2,3,4,1]

      if (nelnodes == 3) then
        nfacenodes = 2
        list(1) = face
        list(2) = i3(face)
      else if (nelnodes == 4) then
        nfacenodes = 2
        list(1) = face
        list(2) = i4(face)
      else if (nelnodes == 6) then
        nfacenodes = 3
        list(1) = face
        list(2) = i3(face)
        list(3) = face+3
      else if (nelnodes == 8) then
        nfacenodes = 3
        list(1) = face
        list(2) = i4(face)
        list(3) = face+4
      else if (nelnodes == 9) then
        nfacenodes = 3
        if (face==1) list(1:3) = (/1,3,2/)
        if (face==2) list(1:3) = (/3,9,6/)
        if (face==3) list(1:3) = (/9,7,8/)
        if (face==4) list(1:3) = (/7,1,4/)
      endif

      end subroutine abq_facenodes_2D_CH

      subroutine abq_inverse_LU_CH(Ain,A_inverse,n)  ! Compute the inverse of an arbitrary matrix by LU decomposition

        implicit none

        integer, intent(in)  :: n

        double precision, intent(in)    :: Ain(n,n)
        double precision, intent(out)   :: A_inverse(n,n)

        double precision :: A(n,n), L(n,n), U(n,n), b(n), d(n), x(n)
        double precision :: coeff
        integer :: i, j, k

        A(1:n,1:n) = Ain(1:n,1:n)
        L=0.d0
        U=0.d0
        b=0.d0

        do k=1, n-1
            do i=k+1,n
                coeff=a(i,k)/a(k,k)
                L(i,k) = coeff
                A(i,k+1:n) = A(i,k+1:n)-coeff*A(k,k+1:n)
            end do
        end do

        forall (i=1:n)  L(i,i) = 1.d0
        forall (j=1:n) U(1:j,j) = A(1:j,j)
        do k=1,n
            b(k)=1.d0
            d(1) = b(1)
            do i=2,n
                d(i)=b(i)
                d(i) = d(i) - dot_product(L(i,1:i-1),d(1:i-1))
            end do
            x(n)=d(n)/U(n,n)
            do i = n-1,1,-1
                x(i) = d(i)
                x(i) = x(i)-dot_product(U(i,i+1:n),x(i+1:n))
                x(i) = x(i)/U(i,i)
            end do
            A_inverse(1:n,k) = x(1:n)
            b(k)=0.d0
        end do

      end subroutine abq_inverse_LU_CH


