!
!    ABAQUS format UEL subroutine
!
!    This file is compatible with both EN234_FEA and ABAQUS/Standard
!
!    The example implements a standard fully integrated 3D linear elastic continuum element
!
!    The file also contains the following subrouines:
!          abq_UEL_3D_integrationpoints           - defines integration ponits for 3D continuum elements
!          abq_UEL_3D_shapefunctions              - defines shape functions for 3D continuum elements
!          abq_UEL_invert3D                       - computes the inverse and determinant of a 3x3 matrix
!          abq_facenodes_3D                       - returns list of nodes on the face of a 3D element
!
!=========================== ABAQUS format user element subroutine ===================

      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
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
    !
    ! Local Variables
      integer      :: i,j,m,l,n_points,kint, nfacenodes, ipoin, ie
      integer      :: face_node_list(8)                       ! List of nodes on an element face
    !
      double precision  ::  xi(3,64)                          ! Volumetric Integration points
      double precision  ::  w(64)                             ! Integration weights
      double precision  ::  N(20)                             ! 3D Shape functions
      double precision  ::  dNdxi(20,3)                       ! 3D Shape function derivatives
      double precision  ::  dxdxi(3,3)                        ! Derivative of position wrt normalized coords
      double precision  ::  dNdx(20,3)                        ! Derivative of shape functions wrt reference coords
      double precision  ::  dNdy(20,3)                        ! Derivative of shape functions wrt deformed coords
    !
    !   Variables below are for computing integrals over element faces
      double precision  ::  face_coords(3,8)                  ! Coords of nodes on an element face
      double precision  ::  xi2(2,9)                          ! Area integration points
      double precision  ::  N2(9)                             ! 2D shape functions
      double precision  ::  dNdxi2(9,2)                       ! 2D shape function derivatives
      double precision  ::  norm(3)                           ! Normal to an element face
      double precision  ::  dxdxi2(3,2)                       ! Derivative of spatial coord wrt normalized areal coord
    !
      double precision  ::  stress(6)                         ! Sigma Stress vector contains [s11, s22, s33, s12, s13, s23]
      double precision  ::  stressmat(3,3)                    ! Sigma Stress matrix
      double precision  ::  Cstressmat(3,3)                   ! Cauchy Stress matrix
      double precision  ::  Cstress(6)                        ! Cauchy Stress vector
      double precision  ::  F(3,3)                            ! Deformation gradient
      double precision  ::  Finv(3,3)                         ! Inverse of deformation gradient
      double precision  ::  C(3,3)                            ! C-G deformation tensor
      double precision  ::  JJ                                ! det(F)
      !double precision  ::  G(6,9)
      double precision  ::  D(6,6)                            ! Material tangent
      double precision  ::  H(6,9)
      double precision  ::  Bstar(9,60)                       ! F = Bstar*(dof_total)
      double precision  ::  q(9,1), q0(3,3)
      !double precision  ::  Pvec(3*NNODE)
      double precision  ::  Y(60,60)
      double precision  ::  Gstif(20,20)
      double precision  ::  HBstar(6,60)
      double precision  ::  DHB(6,60)
      double precision  ::  HDHB(9,60)

      double precision  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant


    !
    !     Example ABAQUS UEL implementing 3D linear elastic elements
    !     El props are:

    !     PROPS(1)         Young's modulus
    !     PROPS(2)         Poisson's ratio

      if (NNODE == 4) n_points = 1               ! Linear tet
      if (NNODE == 10) n_points = 4              ! Quadratic tet
      if (NNODE == 8) n_points = 8               ! Linear Hex
      if (NNODE == 20) n_points = 27             ! Quadratic hex

      call abq_UEL_3D_integrationpoints(n_points, NNODE, xi, w)

      if (MLVARX<3*NNODE) then
        write(6,*) ' Error in abaqus UEL '
        write(6,*) ' Variable MLVARX must exceed 3*NNODE'
        write(6,*) ' MLVARX = ',MLVARX,' NNODE = ',NNODE
        stop
      endif

      RHS(1:MLVARX,1) = 0.d0
      AMATRX(1:NDOFEL,1:NDOFEL) = 0.d0


      ENERGY(1:8) = 0.d0

    !     --  Loop over integration points
      do kint = 1, n_points
        call abq_UEL_3D_shapefunctions(xi(1:3,kint),NNODE,N,dNdxi)
        dxdxi = matmul(coords(1:3,1:NNODE),dNdxi(1:NNODE,1:3))
        call abq_UEL_invert3d(dxdxi,dxidx,determinant)
        dNdx(1:NNODE,1:3) = matmul(dNdxi(1:NNODE,1:3),dxidx)

        Bstar = 0.d0
        Bstar(1,1:3*NNODE-2:3) = dNdx(1:NNODE,1)
        Bstar(2,2:3*NNODE-1:3) = dNdx(1:NNODE,2)
        Bstar(3,3:3*NNODE:3)   = dNdx(1:NNODE,3)
        Bstar(4,1:3*NNODE-2:3) = dNdx(1:NNODE,2)
        Bstar(5,2:3*NNODE-1:3) = dNdx(1:NNODE,1)
        Bstar(6,1:3*NNODE-2:3) = dNdx(1:NNODE,3)
        Bstar(7,3:3*NNODE:3)   = dNdx(1:NNODE,1)
        Bstar(8,2:3*NNODE-1:3) = dNdx(1:NNODE,3)
        Bstar(9,3:3*NNODE:3)   = dNdx(1:NNODE,2)
        
        ! Caculate the deformation gradient
        do i = 1,3
            ie = 3*(NNODE-1)+i
            F(i,1:3) = matmul(U(i:ie:3),dNdx(1:NNODE,1:3))
            F(i,i) = F(i,i) + 1.d0
        end do
        
        ! RIGHT CAUCHY-GREEN TENSOR
        C = matmul(transpose(F),F)
        
        call abq_UEL_invert3d(F,Finv,JJ)
        dNdy(1:NNODE,1:3) = matmul(dNdx(1:NNODE,1:3),Finv)
        
        call secondPK(PROPS(1:NPROPS),NPROPS,F,JJ,stress,D)
        stressmat = 0.d0

        stressmat(1,2) = stress(4)
        stressmat(1,3) = stress(5)
        stressmat(2,3) = stress(6)
        stressmat = stressmat + transpose(stressmat)
        stressmat(1,1) = stress(1)
        stressmat(2,2) = stress(2)
        stressmat(3,3) = stress(3)
        
        q = 0.d0
        q0 = 0.d0
        q0 = matmul(stressmat,transpose(F))
        
        do i=1,3
            q(1,1) = q(1,1) + q0(1,1)
            q(2,1) = q(2,1) + q0(2,1)
            q(3,1) = q(3,1) + q0(3,1)
            q(4,1) = q(4,1) + q0(2,1)
            q(5,1) = q(5,1) + q0(1,2)
            q(6,1) = q(6,1) + q0(3,1)
            q(7,1) = q(7,1) + q0(1,3)
            q(8,1) = q(8,1) + q0(3,2)
            q(9,1) = q(9,1) + q0(2,3)
        end do 
        
        !do i=1,3
        !    q(1,1) = q(1,1) + stressmat(1,i) * F(1,i)
        !    q(2,1) = q(2,1) + stressmat(2,i) * F(1,i)
        !    q(3,1) = q(3,1) + stressmat(3,i) * F(1,i)
        !    q(4,1) = q(4,1) + stressmat(2,i) * F(1,i)
        !    q(5,1) = q(5,1) + stressmat(1,i) * F(2,i)
        !    q(6,1) = q(6,1) + stressmat(3,i) * F(1,i)
        !    q(7,1) = q(7,1) + stressmat(1,i) * F(3,i)
        !    q(8,1) = q(8,1) + stressmat(3,i) * F(2,i)
        !    q(9,1) = q(9,1) + stressmat(2,i) * F(3,i)
        !end do 
            
        H = 0.d0
         
        H(1,1) = F(1,1)
        H(1,5) = F(2,1)
        H(1,7) = F(3,1)
        H(2,2) = F(2,2)
        H(2,4) = F(1,2)
        H(2,9) = F(3,2)
        H(3,3) = F(3,3)
        H(3,6) = F(1,3)
        H(3,8) = F(2,3)
        H(4,1) = F(1,2)
        H(4,2) = F(2,1)
        H(4,4) = F(1,1)
        H(4,5) = F(2,2) 
        H(4,7) = F(3,2)
        H(4,9) = F(3,1)
        H(5,1) = F(1,3)
        H(5,3) = F(3,1)
        H(5,5) = F(2,3)
        H(5,6) = F(1,1)
        H(5,7) = F(3,3)
        H(5,8) = F(2,1)
        H(6,2) = F(2,3)
        H(6,3) = F(3,2)
        H(6,4) = F(1,3)
        H(6,6) = F(1,2)
        H(6,8) = F(2,2)
        H(6,9) = F(3,3)
         
        Gstif = 0.d0
      !   do i = 1,NNODE
      !     	do j = 1,NNODE
      !             do m = 1,3
      !     			do l = 1,3
      !                     Gstif(i,j) = Gstif(i,j) + dNdx(i,m)*
      !1                                stressmat(m,l)*dNdx(j,l)
      !                  end do
      !              end do
      !         end do
      !   end do
        Gstif = matmul(dNdx,matmul(stressmat,transpose(dNdx)))
                
        Y = 0.d0
        do i=1,NNODE
            do j=1,NNODE
                Y(3*i-2,3*j-2) = Gstif(i,j)
                Y(3*i-1,3*j-1) = Gstif(i,j)
                Y(3*i,3*j) = Gstif(i,j)
            end do
        end do
        
        RHS(1:3*NNODE,1) = RHS(1:3*NNODE,1)
     1   - matmul(transpose(Bstar(1:9,1:3*NNODE)),q(1:9,1))*
     2                                          w(kint)*determinant
        
        HBstar(1:6,1:3*NNODE) = matmul(H(1:6,1:9),Bstar(1:9,1:3*NNODE))
        DHB(1:6,1:3*NNODE) = matmul(D(1:6,1:6),HBstar(1:6,1:3*NNODE))
        HDHB(1:9,1:3*NNODE) = matmul(transpose(H(1:6,1:9)),
     1                                   DHB(1:6,1:3*NNODE))
        
  
        AMATRX(1:3*NNODE,1:3*NNODE) = AMATRX(1:3*NNODE,1:3*NNODE)
     1     + matmul(transpose(Bstar(1:9,1:3*NNODE)),
     2     HDHB(1:9,1:3*NNODE))*w(kint)*determinant
     3     + Y(1:3*NNODE,1:3*NNODE)*w(kint)*determinant

        Cstressmat = 0.d0
        Cstressmat = matmul(F,matmul(stressmat,transpose(F)))/JJ
        
        Cstress(1) = Cstressmat(1,1)
        Cstress(2) = Cstressmat(2,2)
        Cstress(3) = Cstressmat(3,3)
        Cstress(4) = Cstressmat(1,2)
        Cstress(5) = Cstressmat(1,3)
        Cstress(6) = Cstressmat(2,3)
        
        
        
        if (NSVARS>=n_points*6) then   ! Store Cauchy stress at each integration point (if space was allocated to do so)
            SVARS(6*kint-5:6*kint) =  Cstress(1:6)/JJ
        endif
      end do


      PNEWDT = 1.d0          ! This leaves the timestep unchanged (ABAQUS will use its own algorithm to determine DTIME)
    !
    !   Apply distributed loads
    !
    !   Distributed loads are specified in the input file using the Un option in the input file.
    !   n specifies the face number, following the ABAQUS convention.
    !
    !   This is coded to apply nominal tractions to the element face (the residual force does not change as the element deforms)
    !
    !
      do j = 1,NDLOAD

        call abq_facenodes_3D(NNODE,iabs(JDLTYP(j,1)),
     1                                     face_node_list,nfacenodes)

        do i = 1,nfacenodes
            face_coords(1:3,i) = coords(1:3,face_node_list(i))
        end do

        if (nfacenodes == 3) n_points = 3
        if (nfacenodes == 6) n_points = 4
        if (nfacenodes == 4) n_points = 4
        if (nfacenodes == 8) n_points = 9

        call abq_UEL_2D_integrationpoints(n_points, nfacenodes, xi2, w)

        do kint = 1,n_points
            call abq_UEL_2D_shapefunctions(xi2(1:2,kint),
     1                        nfacenodes,N2,dNdxi2)
            dxdxi2 = matmul(face_coords(1:3,1:nfacenodes),
     1                           dNdxi2(1:nfacenodes,1:2))
            norm(1)=(dxdxi2(2,1)*dxdxi2(3,2))-(dxdxi2(2,2)*dxdxi2(3,1))
            norm(2)=(dxdxi2(1,1)*dxdxi2(3,2))-(dxdxi2(1,2)*dxdxi2(3,1))
            norm(3)=(dxdxi2(1,1)*dxdxi2(2,2))-(dxdxi2(1,2)*dxdxi2(2,1))

            do i = 1,nfacenodes
                ipoin = 3*face_node_list(i)-2
                RHS(ipoin:ipoin+2,1) = RHS(ipoin:ipoin+2,1)
     1                 - N2(1:nfacenodes)*adlmag(j,1)*norm(1:3)*w(kint)      ! Note determinant is already in normal
            end do
        end do
      end do

      return

      END SUBROUTINE UEL


       subroutine secondPK(ele_prop,n_properties,F,JJ,pkstress,D)

       implicit none

       integer, intent(in)           :: n_properties
       double precision, intent(in)  :: ele_prop(n_properties)
       double precision, intent(in)  :: F(3,3)
       double precision, intent(in)  :: JJ
       double precision, intent(out) :: pkstress(6)
       double precision, intent(out) :: D(6,6)

       double precision :: C(3,3)
       double precision :: Cinv(3,3)
       double precision :: Cvec(6)
       double precision :: Cbar(6)
       double precision :: Cstar(6)
       double precision :: Cstarbar(6)
       double precision :: Cinvvec(6)
       double precision :: eyevec(6)
       double precision :: ss
       double precision :: trC
       double precision :: mu,K,G11,G22,G33,G44                !material properties
       double precision :: G(6,6)
       double precision :: P(6)
       double precision :: Q
       double precision :: omega(6,6)
       double precision :: D1(6,6),D2(6,6),D3(6,6),D4(6,6),D5(6,6)
       

       integer :: i
       


       mu = ele_prop(1)
       K  = ele_prop(2)
       G11 = ele_prop(3)
       G22 = ele_prop(4)
       G33 = ele_prop(5)
       G44 = ele_prop(6)
       
       G=0.d0
       G(1,1) = G11
       G(2,2) = G22
       G(3,3) = G33
       G(4,4) = G44
       G(5,5) = G44
       G(6,6) = G44
       
       
       

       !  This subroutine calculates the Kirchhoff stress tau = J*cauchy_stress (stored as a vector stress(i) = [tau_11, tau_22, tau_33, etc]
       !  and the tangent matrix D[I,J] = [dtau_11/dB_11, dtau_11/dB_22,... 
       !                                   dtau_22/dB_11, dtau_22/dB_22,
       !                                    etc
       
       C = matmul(transpose(F),F)
       Cvec(1) = C(1,1)
       Cvec(2) = C(2,2)
       Cvec(3) = C(3,3)
       Cvec(4) = C(1,2)
       Cvec(5) = C(1,3)
       Cvec(6) = C(2,3)
       
       call abq_UEL_invert3d(C,Cinv,ss)
       ss = JJ**(-2.d0/3.d0)
       
       Cbar = 0.d0
       Cbar = Cvec/ss
       
       Cinvvec = 0.d0
       Cinvvec(1) = Cinv(1,1)
       Cinvvec(2) = Cinv(2,2)
       Cinvvec(3) = Cinv(3,3)
       Cinvvec(4) = Cinv(1,2)
       Cinvvec(5) = Cinv(1,3)
       Cinvvec(6) = Cinv(2,3)
       
       Cstar = 0.d0
       Cstar(1) = C(1,1)
       Cstar(2) = C(2,2)
       Cstar(3) = C(3,3)
       Cstar(4) = C(1,2)*2.d0
       Cstar(5) = C(1,3)*2.d0
       Cstar(6) = C(2,3)*2.d0
       
       Cstarbar = 0.d0
       Cstarbar = Cstar/ss
       
       eyevec(1:3) = 1.d0
       eyevec(4:6) = 0.d0
       
       Q = 0.d0
       Q = dot_product(Cstarbar-eyevec,matmul(G,Cstarbar-eyevec))/4.d0
       
       P=0.d0
       
       P=(matmul(G,Cstarbar-eyevec)-dot_product(Cstar,
     1    matmul(G,Cstarbar-eyevec)*Cinvvec/3.d0))/(2.d0*ss)
       
       pkstress = 0.d0
       pkstress = mu*exp(Q)*P - K*JJ*(JJ-1.d0)*Cinvvec
       
       omega =0.d0

       omega(1,2) = Cinvvec(4)*Cinvvec(4)
       omega(1,3) = Cinvvec(5)*Cinvvec(5)
       omega(1,4) = Cinvvec(1)*Cinvvec(4)
       omega(1,5) = Cinvvec(1)*Cinvvec(5)
       omega(1,6) = Cinvvec(4)*Cinvvec(5)
       omega(2,3) = Cinvvec(6)*Cinvvec(6)
       omega(2,4) = Cinvvec(4)*Cinvvec(2)
       omega(2,5) = Cinvvec(4)*Cinvvec(6)
       omega(2,6) = Cinvvec(2)*Cinvvec(6)
       omega(3,4) = Cinvvec(5)*Cinvvec(6)
       omega(3,5) = Cinvvec(5)*Cinvvec(3)
       omega(4,5) = (Cinvvec(1)*Cinvvec(6)+
     1               Cinvvec(5)*Cinvvec(4))/2.d0
       omega(4,6) = (Cinvvec(4)*Cinvvec(6)+
     1               Cinvvec(5)*Cinvvec(2))/2.d0
       omega(5,6) = (Cinvvec(4)*Cinvvec(3)+
     1               Cinvvec(5)*Cinvvec(6))/2.d0

       omega = omega + transpose(omega)
       
       omega(1,1) = Cinvvec(1)*Cinvvec(1)/2.d0
       omega(2,2) = Cinvvec(2)*Cinvvec(2)/2.d0
       omega(3,3) = Cinvvec(3)*Cinvvec(3)/2.d0
       omega(4,4) = (Cinvvec(1)*Cinvvec(2)+
     1               Cinvvec(4)*Cinvvec(4))/2.d0
       omega(5,5) = (Cinvvec(1)*Cinvvec(3)+
     1               Cinvvec(5)*Cinvvec(5))/2.d0
       omega(6,6) = (Cinvvec(2)*Cinvvec(3)+
     1               Cinvvec(6)*Cinvvec(6))/2.d0
       
       D = 0.d0
       D1 = G-(spread(matmul(G,Cstar),dim=2,ncopies=6)*
     1   spread(Cinvvec,dim=1,ncopies=6)+spread(Cinvvec,dim=2,
     2   ncopies=6)*spread(matmul(G,Cstar),dim=1,ncopies=6))/3.d0
       
       D2 = -ss/3.d0*dot_product(Cstar,
     1   matmul(G,Cstarbar-eyevec))*omega
       
       D3 = 1.d0/9.d0*dot_product(Cstar,matmul(G,Cstar))
     1   *spread(Cinvvec,dim=2,ncopies=6)
     2   *spread(Cinvvec,dim=1,ncopies=6)
       
       D4 = 2.d0*spread(P,dim=2,ncopies=6)*
     1   spread(P-(Cinvvec/3.d0),dim=1,ncopies=6)-1.d0/(3.d0*ss)*
     2   spread(Cinvvec,dim=2,ncopies=6)*
     3   spread(matmul(G,Cstarbar-eyevec),dim=1,ncopies=6)
       
       D5 = K*JJ*((2.d0*JJ-1.d0)*spread(Cinvvec,dim=2,ncopies=6)
     1   *spread(Cinvvec,dim=1,ncopies=6)+2.d0*(JJ-1.d0)*omega)
       
       
       D = mu*exp(Q)*(D1+D2+D3)/(ss**2) + mu*exp(Q)*D4 + D5

       return

      end subroutine secondPK