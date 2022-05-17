!*** s groh, may 19 2010.
!* problem or questions, contact sebstephg@gmail.com
!* create dislocation structure for different caracter, geometry and crystal structures
!*** what are the possibility
!*** FCC
!* EDGE in PAD
!* EDGE in Cylinder
!* Perfect with edge orientation
!* SCREW in Cylinder

!*** BCC
!* EDGE in PAD
!* EDGE in Cylinder
!* Perfect with edge orientation
!* SCREW in Cylinder

!*** screw dislocations are not implemented for the PAD. that's something that needs to be done in the future

program Dislocation

!use varglob
!use displacement

implicit none

integer                                           :: CrystalStructure, FCCMaterial, BCCMaterial, HCPMaterial
integer                                           :: DislocationCaracter, GeometryChoice, CoreChoice
integer                                           :: HCPOrientation
integer                                           :: AnswerPerfect, B
integer                                           :: atomspercel, ncelltot, numatomstot
integer                                           :: i, j, k, l
integer                                           :: NumCellP, numdel, numcell, numhide
integer                                           :: m, iun, isotropic
integer,dimension(4)                              :: N
integer,dimension(:),allocatable                  :: itype
integer,dimension(:,:),allocatable                :: itypenumcell
integer :: crack, deleteatoms
real(kind=8) :: aa,bb
real(kind=8) :: ye, ze, yye, zze
real(kind=8) :: xcrack, ymaxcrack, ymincrack,ycrack

real(kind=8)                                      :: R_2, R_3, R_6, epsilon, rzero
real(kind=8)                                      :: nu, nuFe, nuNb, nuTa
real(kind=8)                                      :: burgers, avalue, acvalue, answerlattice
real(kind=8)                                      :: pivalue, pisurtrois, factor, pioversix
real(kind=8)                                      :: radius, XDislo, YDislo, distance, xx, yy
real(kind=8)                                      :: scaleX, bx, by, nua, nub
real(kind=8)                                      :: ux, uy, theta, x1, x2
real(kind=8)                                      :: c11, c12 ,c44, H, cp44, cp45, cp55, Angle
real(kind=8),dimension(2)                         :: center
real(kind=8),dimension(3)                         :: midle, midle2, dimeunit, r
real(kind=8),dimension(3)                         :: UnitCellDim, sr1, Vburger, px, py
real(kind=8),dimension(3,2)                       :: dimbox
real(kind=8),dimension(:,:), allocatable          :: AtomsCellCoor
real(kind=8),dimension(:,:,:),allocatable         :: AtomsTot
logical,dimension(:,:),allocatable                :: hide

!*** constant
!*** 
rzero=0.0d0
R_3 = sqrt(3.0d0)
R_6 = sqrt(6.0d0)
R_2 = sqrt(2.0d0)
pivalue = 2.0d0*asin(1.0d0) !* given in radians
pisurtrois = pivalue / 3.0d0
pioversix = pivalue / 6.0d0

write(*,*) 'Which crystal structure do you want to consider?'
write(*,*) '!++++++++++++++++++++++++!'
write(*,*) '!++     1 == FCC       ++!'
write(*,*) '!++     2 == BCC       ++!'
write(*,*) '!++     3 == HCP       ++!'
write(*,*) '!++     4 == Diamonds  ++!'
write(*,*) '!++++++++++++++++++++++++!'
read(*,*) CrystalStructure

select case (CrystalStructure)

case(1)
  write(*,*) 'You are about to construct a dislocation in FCC'
  !*** general parameters
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(*,*) '++   what kind of material do you want to consider      ++'
  write(*,*) '++ 1 == Ni                                              ++'
  write(*,*) '++ 2 == Cu                                              ++'
  write(*,*) '++ 3 == Al                                              ++'
  write(*,*) '++ 4 == open value for lattice parameter                ++'
  write(*,*) '=========================================================='
  read(*,*) FCCMaterial
  if (FCCMaterial == 1) then
   avalue = 3.52d0
  else if (FCCMaterial == 2) then
   avalue = 3.62d0
  else if (FCCMaterial == 3) then
   avalue = 4.05d0
  else if (FCCMaterial == 4) then
   write(*,*) 'what is the desired value of lattice parameter?'
  read(*,*) avalue
  endif

  !*** magnitude of the Burgers vector
  burgers = avalue/R_2

  !*** number of atoms per unit cell
  atomspercel = 6
  allocate(AtomsCellCoor(atomspercel,3))
  !*** position of the atoms in the reference unit cell
  !## atoms #1
  AtomsCellCoor(1,1) = avalue/2/R_2 
  AtomsCellCoor(1,2) = 0.0            
  AtomsCellCoor(1,3) = avalue/2/R_6*3
  !## atoms #2
  AtomsCellCoor(2,1) = avalue/2/R_2 
  AtomsCellCoor(2,2) = avalue/R_3     
  AtomsCellCoor(2,3) = avalue*5*R_6/12
  !## atoms #3
  AtomsCellCoor(5,1) = avalue/2/R_2 
  AtomsCellCoor(5,2) = avalue*2/R_3   
  AtomsCellCoor(5,3) = avalue/2/R_6
  !## atoms #4
  AtomsCellCoor(3,1) = 0.0d0        
  AtomsCellCoor(3,2) = 0.0d0          
  AtomsCellCoor(3,3) = 0.0d0
  !## atoms #5
  AtomsCellCoor(4,1) = 0.0d0        
  AtomsCellCoor(4,2) = avalue/R_3     
  AtomsCellCoor(4,3) = avalue/R_6
  !## atoms #6
  AtomsCellCoor(6,1) = 0.0d0        
  AtomsCellCoor(6,2) = avalue*2/R_3   
  AtomsCellCoor(6,3) = avalue*2/R_6
  
  !*** size of the unit cell
  UnitCellDim(1) = avalue/R_2; UnitCellDim(2) = avalue*R_3; UnitCellDim(3) = 3/R_6*avalue
  
  !** dimension of the reference unit cell
  write(*,*) 'dimension reference unit cell :',UnitCellDim
  ! N(1) : must be an odd number
  ! N(2) : must be an even number
  ! N(3) : no constrain
  write(*,*) 'number of unit cells along X, Y and Z'
  read(*,*) N(1),N(2),N(3)

  !*** total number of unit cell
  ncelltot = N(1)*N(2)*N(3)

  allocate(AtomsTot(ncelltot,atomspercel,3))
  !** generate the atomic positions
  NumCellP = 0
  do i=1,N(1)
   do j=1,N(2)
    do k=1,N(3)
     NumCellP = NumCellP + 1
     do l=1,atomspercel
       AtomsTot(NumCellP,l,1) = AtomsCellCoor(l,1) + (i-1)*UnitCellDim(1)
       AtomsTot(NumCellP,l,2) = AtomsCellCoor(l,2) + (j-1)*UnitCellDim(2)
       AtomsTot(NumCellP,l,3) = AtomsCellCoor(l,3) + (k-1)*UnitCellDim(3)
     enddo
    enddo
   enddo
  enddo
  
  write(*,*) ' What kind of dislocation do you want?'    
  write(*,*) '!++++++++++++++++++++++++++++++++++++++++++!'
  write(*,*) '!++       1 == Edge                      ++!'
  write(*,*) '!++       2 == Screw                     ++!'
  write(*,*) '!++++++++++++++++++++++++++++++++++++++++++!'
  read(*,*) DislocationCaracter
  if (DislocationCaracter == 1) then
      write(*,*) 'You are about the build an EDGE dislocation in FCC crystal'
      !*** here we need to define the type of geometry we want to use
      write(*,*) '-------------------------------------------'
      write(*,*) '!** What kind of structure do you want? **!'
      write(*,*) '!** 1 == cylinder                       **!'
      write(*,*) '!** 2 == PAD                            **!'
      write(*,*) '!** 3 == Perfect FCC crystal            **!'
      write(*,*) '-------------------------------------------'
      read(*,*) GeometryChoice 
      
      if (GeometryChoice == 3) then
        write(*,*) 'you are about to build a perfect FCC crystal with edge orientation' 
         write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
         write(*,*) '!++ You are going to create a perfect FCC crystal with orientation : ++!'
         write(*,*) '!++ X -> [110]                                                       ++!'
         write(*,*) '!++ Y -> [111]                                                       ++!'
         write(*,*) '!++ Z -> [112]                                                       ++!' 
         write(*,*) '!++ What kind of shape do you want to build?                         ++!'
         write(*,*) '!++ 1 == cubic                                                       ++!'
         write(*,*) '!++ 2 == cylinder                                                    ++!'
         write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
         read(*,*) AnswerPerfect
         if (AnswerPerfect == 1) then
           numatomstot = N(1)*N(2)*N(3)*atomspercel
           dimbox(:,:) = 0.0d0
           dimbox(1,2) = N(1)*UnitCellDim(1)
           dimbox(2,2) = N(2)*UnitCellDim(2)
           dimbox(3,2) = N(3)*UnitCellDim(3)
           open(1,file='atoms.fcc.perfect.cubic',status='unknown')
           write(1,*) 'Position data for Ni File'
           write(1,*) 
           write(1,*) numatomstot, 'atoms'
           write(1,*) ' 1 atom types'
           write(1,*) dimbox(1,1:2),' xlo xhi'
           write(1,*) dimbox(2,1:2),' ylo yhi'
           write(1,*) dimbox(3,1:2),' zlo zhi'
           write(1,*) 
           write(1,*) 'Atoms'
           write(1,*)
           m = 1
           iun = 1
           NumCellP=0
           do i=1,N(1)
            do j=1,N(2)
             do k=1,N(3)
               NumCellP = NumCellP+1
               do l=1,atomspercel
                write(1,10) m,iun,AtomsTot(NumCellP,l,1:3)
                m = m+1
               enddo
             enddo
            enddo
           enddo 
           close(1)
          
         else if (AnswerPerfect == 2) then
            allocate(hide(ncelltot,6))

            hide(:,:) = .true.
            !** dimension of the simulation cell
            dimbox(:,:) = 0.0d0
            dimbox(1,2) = N(1)*UnitCellDim(1)
            dimbox(2,2) = N(2)*UnitCellDim(2)
            dimbox(3,2) = N(3)*UnitCellDim(3)
            XDislo = 0.5d0*(dimbox(1,2)+dimbox(1,1))
            YDislo = 0.5d0*(dimbox(2,2)+dimbox(2,1))
            !* hide if outside of the cylinder of the radius R
            write(*,*) 'Center of the cylinder = ', XDislo, YDislo
            write(*,*) 'Box size =',dimbox(1,2),dimbox(2,2),dimbox(3,2)
            write(*,*) 'radius of the cylinder'
            read(*,*) radius
            numdel = 0
            NumCellP = 0
            do i=1,N(1)
             do j=1,N(2)
              do k=1,N(3)
                NumCellP = NumCellP + 1
                do l=1,atomspercel
                  !** check if atom is inside of cylinder of radius with axis on XDislo, YDislo
                  xx = AtomsTot(NumCellP,l,1) - XDislo
                  yy = AtomsTot(NumCellP,l,2) - YDislo      
                  distance = sqrt(xx*xx+yy*yy)
                  if (distance .gt. radius) then
                   hide(NumCellP,l) = .false.
                   numdel = numdel + 1
                  endif
                enddo
              enddo
             enddo
            enddo
 
           numatomstot = N(1)*N(2)*N(3)*atomspercel - numdel
           !*** save the perfect crystal structure
           open(1,file='atoms.fcc.perfect.cylinder',status='unknown')
           write(1,*) 'Position data for Mg File'
           write(1,*) 
           write(1,*) numatomstot, 'atoms'
           write(1,*) ' 1 atom types'
           write(1,*) dimbox(1,1:2),' xlo xhi'
           write(1,*) dimbox(2,1:2),' ylo yhi'
           write(1,*) dimbox(3,1:2),' zlo zhi'
           write(1,*) 
           write(1,*) 'Atoms'
           write(1,*)
           m = 1
           iun = 1
           NumCellP=0
           do i=1,N(1)
            do j=1,N(2)
             do k=1,N(3)
               NumCellP = NumCellP+1
               do l=1,6
                if (hide(NumCellP,l)) then
                 write(1,10) m,iun,AtomsTot(NumCellP,l,1:3)
                 m = m+1
                endif
               enddo
             enddo
            enddo
           enddo

           close(1)
           deallocate(hide)
         
         endif
    
      else if (GeometryChoice == 2) then
        write(*,*) 'you are about to build a PAD made of edge dislocation in FCC'
        allocate(hide(ncelltot,6))
        hide(:,:) = .false.
        !** remove half a plane of atoms
        NumCellP = 0
        numdel = 0
        do i=1,N(1)
         do j=1,N(2)
          do k=1,N(3)
            NumCellP = NumCellP +1
            if (i.eq.N(1)) then
             if (j.le.N(2)/2) then
               hide(NumCellP,1:6) = .true.
               numdel = numdel + 6
             endif
            endif
          enddo
         enddo
        enddo

        numatomstot = N(1)*N(2)*N(3)*atomspercel-numdel
        
        !** scaling by +b/2 for atoms with j <= N(2)/2
        !** scaling by -b/2 for atoms with j > N(2)/2
        numcell = 0
        do i=1,N(1)
         do j=1,N(2)
          do k=1,N(3)
           numcell = numcell + 1
           do l=1,atomspercel
            if (j.le.(N(2)/2)) then
             !*** scale all the x coordinates by +b/2
             xx=AtomsTot(numcell,l,1)
             scaleX = 0.5d0*xx*burgers/((N(1)-1)*UnitCellDim(1))
             AtomsTot(numcell,l,1) = AtomsTot(numcell,l,1) + scaleX
            elseif (j .gt. (N(2)/2)) then
             !*** scale all the x coordinates by -b/2
             xx=AtomsTot(numcell,l,1)
             scaleX = -0.5d0*xx*burgers/((N(1)-1)*UnitCellDim(1))
             AtomsTot(numcell,l,1) = AtomsTot(numcell,l,1) + scaleX
            endif
           enddo
          enddo
         enddo
        enddo
        !*** dimension of the simulation cell
        dimbox(:,:) = 0.0d0
        dimbox(1,2) = (N(1)-1)*UnitCellDim(1)+burgers*0.5d0
        dimbox(2,2) = N(2)*UnitCellDim(2)
        dimbox(3,2) = N(3)*UnitCellDim(3)


        !*** Save the atomic position using LAMMPS format
        open(1,file='atoms.fcc.edge.pad',status='unknown')
        write(1,*) 'Position data for Ni File'
        write(1,*) 
        write(1,*) numatomstot, 'atoms'
        write(1,*) ' 1 atom types'
        write(1,*) dimbox(1,1:2),' xlo xhi'
        write(1,*) dimbox(2,1:2),' ylo yhi'
        write(1,*) dimbox(3,1:2),' zlo zhi'
        write(1,*) 
        write(1,*) 'Atoms'
        write(1,*)
        m = 1
        iun = 1
        NumCellP=0
        do i=1,N(1)
         do j=1,N(2)
          do k=1,N(3)
            NumCellP = NumCellP+1
            do l=1,atomspercel
             if (.not.hide(NumCellP,l)) then
              write(1,10) m,iun,AtomsTot(NumCellP,l,1:3)
              m = m+1
             endif
            enddo
          enddo
         enddo
        enddo
        
        deallocate(hide)      
        write(*,*) 'for blocked atoms on the bottom and top surface'
        write(*,*) 'ymin == ',UnitCellDim(2)-0.1d0
        write(*,*) 'ymax == ',(N(2)-1)*UnitCellDim(2) - 0.1d0

  
      else if (GeometryChoice == 1) then
        write(*,*) 'you are about to build single edge dislocation in a cylinder in FCC'
        allocate(hide(ncelltot,atomspercel))
        hide(:,:) = .false.

        NumCellP = 0
        numdel = 0
        do i=1,N(1)
         do j=1,N(2)
          do k=1,N(3)
            NumCellP = NumCellP +1
            do l=1,atomspercel
             if (i.eq.((N(1)-1)/2+1)) then
              if (j.le.N(2)/2) then
                hide(NumCellP,l) = .true.
                numdel = numdel + 1
              endif
             endif
            enddo
          enddo
         enddo
        enddo

        numatomstot = N(1)*N(2)*N(3)*atomspercel-numdel      

        XDislo = UnitCellDim(1)*((N(1)-1)/2+0.5d0) 
        !XDislo = ((N(1)-1)/2)*UnitCellDim(1)
        YDislo = UnitCellDim(2)*(N(2)/2)-0.1
        sr1(1)=XDislo
        sr1(2) = YDislo
        sr1(3)=0.0d0
        factor = burgers/2.0d0/pivalue
        Vburger(1) = burgers;Vburger(2) = 0.0d0;Vburger(3) = 0.0d0
        bx = Vburger(1);by=Vburger(2)
        !*** projection on X
        px(1) = 1.0d0;px(2) = 0.0d0;px(3) = 0.0d0
        !*** projection on Y
        py(1) = 0.0d0;py(2) = 1.0d0;py(3) = 0.0d0
        nua = 1.0d0/4.0d0/(1-0.33d0)/2/pivalue
        nub = (1.0d0-2.0d0*0.33d0)*nua      
      
        !*** apply the displacement field of a dislocation
        NumCellP = 0
        do i=1,N(1)
         do j=1,N(2)
          do k=1,N(3)
           NumCellP = NumCellP + 1
           do l = 1,atomspercel
            if (.not.hide(NumCellP,l)) then
             xx = AtomsTot(NumCellP,l,1) - XDislo
             yy = AtomsTot(NumCellP,l,2) - YDislo

             theta=datan2(-xx,yy)
    
             ux = factor*(theta+xx*yy/2/(1-0.33d0)/(xx*xx+yy*yy))
             uy = -1.0d0*factor*((1.0d0-2.0d0*0.33d0)/4/(1-0.33d0)*log(xx*xx+yy*yy)+(xx*xx-yy*yy)/4.0d0/(1-0.33d0)/(xx*xx+yy*yy))
    
             AtomsTot(NumCellP,l,1) = AtomsTot(NumCellP,l,1) + ux
             AtomsTot(NumCellP,l,2) = AtomsTot(NumCellP,l,2) + uy
            endif
           enddo
          enddo
         enddo
        enddo
        write(*,*) 'Dislo position'
        write(*,*) XDislo, YDislo
        write(*,*) 'radius of the cylinder'
        read(*,*) radius
        
        NumCellP = 0
        do i=1,N(1)
         do j=1,N(2)
          do k=1,N(3)
            NumCellP = NumCellP + 1
            do l=1,atomspercel
             if (.not.hide(NumCellP,l)) then
              !** check if atom is inside of cylinder of radius with axis on XDislo, YDislo
              xx = AtomsTot(NumCellP,l,1) - XDislo
              yy = AtomsTot(NumCellP,l,2) - YDislo      
              distance = sqrt(xx*xx+yy*yy)
              if (distance .gt. radius) then
               hide(NumCellP,l) = .true.
               numdel = numdel + 1
              endif
             endif
            enddo
          enddo
         enddo
        enddo
        
        !*** total number of atoms
        numatomstot = N(1)*N(2)*N(3)*atomspercel-numdel
        
        !*** dimension of the box
        dimbox(:,:)=0.0d0
        dimbox(1,2) = N(1)*UnitCellDim(1)
        dimbox(2,2) = N(2)*UnitCellDim(2)
        dimbox(3,2) = N(3)*UnitCellDim(3)

        open(1,file='atoms.fcc.edge.cylinder',status='unknown')
        write(1,*) 'Position data for Ni File'
        write(1,*) 
        write(1,*) numatomstot, 'atoms'
        write(1,*) ' 1 atom types'
        write(1,*) dimbox(1,1:2),' xlo xhi'
        write(1,*) dimbox(2,1:2),' ylo yhi'
        write(1,*) dimbox(3,1:2),' zlo zhi'
        write(1,*) 
        write(1,*) 'Atoms'
        write(1,*)
        m = 1
        iun = 1
        NumCellP=0
        do i=1,N(1)
         do j=1,N(2)
          do k=1,N(3)
            NumCellP = NumCellP+1
            do l=1,atomspercel
             if (.not.hide(NumCellP,l)) then
              write(1,10) m,iun,AtomsTot(NumCellP,l,1:3)
              m = m+1
             endif
            enddo
          enddo
         enddo
        enddo

       close(1)

      deallocate(hide)
        
      endif
  
  else if (DislocationCaracter == 2) then
    write(*,*) 'you are about the build a SCREW dislocation in FCC crystal'
    !*** here we need to define the type of geometry we want to use
    write(*,*) '-------------------------------------------'
    write(*,*) '!** What kind of structure do you want? **!'
    write(*,*) '!** 1 == cylinder                       **!'
    write(*,*) '!** 2 == PAD                            **!'
    write(*,*) '-------------------------------------------'
    read(*,*) GeometryChoice 
    if (GeometryChoice == 1) then
      !*** define the position of the dislocation line
      midle(:) = 0.0d0
      midle2(:) = 0.0d0
      epsilon = 0.0001D-01
      midle(1) = N(1)/2*UnitCellDim(1) - epsilon
      midle(3) = N(3)/2*UnitCellDim(3) - 20.0d0   !** Z-
      midle2(3) = N(3)/2*UnitCellDim(3) + 20.0d0  !** Z+
      midle(2) = N(2)/2*UnitCellDim(2)           !** Y
   
      center(2) = N(2)/2.0d0*UnitCellDim(2) !** Y
      center(1) = N(3)/2.0d0*UnitCellDim(3) !** Z

      write(*,*) 'dislocation line (X, Y, Z)'
      write(*,*) midle(1),center(2),center(1)

      allocate(hide(ncelltot,atomspercel))
      hide(:,:) = .false.
   
      !*** define the radius of the cylinder
      write(*,*) 'What is the cylinder radius?'
      read(*,*) radius

   
      !*** apply the dislocation displacement from linear elasticity
      !*** cut a cylinder around the dislocation line
      NumCellP = 0
      numhide=0
      do i=1,N(1)
       do j=1,N(2)
        do k=1,N(3)
         NumCellP = NumCellP + 1
         do l=1,6
      
          !* tan**(1)(y/z)
           xx = AtomsTot(NumCellP,l,2)
           yy = AtomsTot(NumCellP,l,3)
        
           !*** screw dislocation displacement
           AtomsTot(NumCellP,l,1) = AtomsTot(NumCellP,l,1) + burgers/4.0d0/pivalue*(datan2((xx-midle(2)),(yy-midle(3))))
           AtomsTot(NumCellP,l,1) = AtomsTot(NumCellP,l,1) + burgers/4.0d0/pivalue*(datan2((xx-midle(2)),(yy-midle2(3))))

           !*** check if atoms inside the cylinder
           xx = xx - center(2)
           yy = yy - center(1)
           distance = sqrt(xx*xx+yy*yy)
           if (distance .gt. radius) then
             numhide = numhide + 1
             hide(NumCellP,l) = .true.
           endif

         enddo
        enddo
       enddo
      enddo
      
      !*** total number of atoms
      numatomstot = N(1)*N(2)*N(3)*atomspercel-numhide
      
      dimbox(:,:)=0.0d0
      dimbox(1,2) = N(1)*UnitCellDim(1)
      dimbox(2,2) = N(2)*UnitCellDim(2)
      dimbox(3,2) = N(3)*UnitCellDim(3)
   
      open(1,file='atoms.fcc.screw.cylinder',status='unknown')
   
      write(1,*) 'Position data for Cu File'
      write(1,*) 
      write(1,*) numatomstot, ' atoms'
      write(1,*) ' 1 atom types'
      write(1,*) dimbox(1,1:2),' xlo xhi'
      write(1,*) dimbox(2,1:2),' ylo yhi'
      write(1,*) dimbox(3,1:2),' zlo zhi'
      write(1,*)
      write(1,*) 'Atoms'
      write(1,*)
      m=1
      iun=1
      NumCellP=0
      do i=1,N(1)
       do j=1,N(2)
        do k=1,N(3)
         NumCellP = NumCellP+1
         do l=1,6
          if (.not.hide(NumCellP,l)) then
            write(1,10) m, iun, AtomsTot(NumCellP,l,1:3)
            m  = m+1
          endif
         enddo
        enddo
       enddo
      enddo
   
     close(1)


      
    else if (GeometryChoice == 2) then
      write(*,*) 'Infinite screw dislocation in a PAD for FCC'
      write(*,*) 'NOT IMPLEMENTED YET'
     !*** deallocate the memory
     deallocate(AtomsCellCoor)
     deallocate(AtomsTot)     
     stop
    endif
    
  
  endif

  !*** deallocate the memory
  deallocate(AtomsCellCoor)
  deallocate(AtomsTot)
  
case(2)
  write(*,*) 'you are about to construct a dislocation in BCC'
  !*** general parameters
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(*,*) '++   what kind of material do you want to consider      ++'
  write(*,*) '++ 1 == Fe                                              ++'
  write(*,*) '++ 2 == Nb                                              ++'
  write(*,*) '++ 3 == Ta                                              ++'
  write(*,*) '++ 4 == Mo                                              ++'
  write(*,*) '++ 5 == W (Zhou)                                        ++'
  write(*,*) '++ 6 == W (Auckland)                                    ++'
  write(*,*) '=========================================================='
  read(*,*) BCCMaterial
  if (BCCMaterial == 1) then
   avalue = 2.855324d0 ! Fe EAM
   nuFe = 0.3571429d0 ! Fe EAM
   nu = nuFe ! Fe EAM
   avalue = 2.850999d0 ! Fe MEAM
   nuFe = 0.324 ! Fe MEAM
   nu = nuFe ! Fe MEAM

  else if (BCCMaterial == 2) then
    avalue = 3.308d0 ! Nb
    nuNb = 0.4047619d0
    nu = nuNb
  else if (BCCMaterial == 3) then
    avalue = 3.302633d0
    nuTa = 0.2276639d0
    nu = nuTa
  else if (BCCMaterial == 4) then
    avalue = 3.1501534d0 ! Mo
    nu = 0.1493893d0
  else if (BCCMaterial == 5) then
    avalue = 3.1648495d0 !W
    nu = 0.1457d0
  else if (BCCMaterial == 6) then
    avalue = 3.1652d0 !W
    nu = 0.161d0
  endif
  
  burgers = avalue*R_3/2
  pivalue = 2.0d0*asin(1.0d0) !* given in radians
  pisurtrois = pivalue / 3.0d0
  !*** atoms for one unit cell
  atomspercel = 6
  allocate(AtomsCellCoor(atomspercel,3))
  dimeunit(1) = avalue*R_3/2.0d0
  dimeunit(2) = avalue*R_2
  dimeunit(3) = avalue*R_6

  !** atom #1
  AtomsCellCoor(1,1)=0.0d0
  AtomsCellCoor(1,2)=0.0d0
  AtomsCellCoor(1,3)=0.0d0
  !** atom #2
  AtomsCellCoor(2,1)=dimeunit(1)*(2.0d0/3.0d0)
  AtomsCellCoor(2,2)=dimeunit(2)*0.0d0
  AtomsCellCoor(2,3)=dimeunit(3)*(1.0d0/3.0d0)
  !** atom #3
  AtomsCellCoor(3,1)=dimeunit(1)*(1.0d0/3.0d0)
  AtomsCellCoor(3,2)=dimeunit(2)*(0.0d0)
  AtomsCellCoor(3,3)=dimeunit(3)*(2.0d0/3.0d0)
  !** atom #4
  AtomsCellCoor(4,1)=dimeunit(1)*(1.0d0/3.0d0)
  AtomsCellCoor(4,2)=dimeunit(2)*(0.5d0)
  AtomsCellCoor(4,3)=dimeunit(3)*(1.0d0/6.0d0)
  !** atom #5
  AtomsCellCoor(5,1)=dimeunit(1)*(0.0d0/3.0d0)
  AtomsCellCoor(5,2)=dimeunit(2)*(0.5d0)
  AtomsCellCoor(5,3)=dimeunit(3)*(1.0d0/2.0d0)
  !** atom #6
  AtomsCellCoor(6,1)=dimeunit(1)*(2.0d0/3.0d0)
  AtomsCellCoor(6,2)=dimeunit(2)*(0.5d0)
  AtomsCellCoor(6,3)=dimeunit(3)*(5.0d0/6.0d0)

  !*** Size of the unit cell
  UnitCellDim(1) = avalue/R_2; UnitCellDim(2) = avalue*R_3; UnitCellDim(3) = 3/R_6*avalue

  UnitCellDim(1)=dimeunit(1) 
  UnitCellDim(2)=dimeunit(2)
  UnitCellDim(3)=dimeunit(3)

  write(*,*) 'Number of unit cells along X (odd*',UnitCellDim(1),') , Y (even*',UnitCellDim(2),') and Z',UnitCellDim(3)
  read(*,*) N(1), N(2), N(3)

  ncelltot = N(1)*N(2)*N(3)

  allocate(AtomsTot(ncelltot,atomspercel,3))

  !** generate the atomic positions
  NumCellP = 0
  do i=1,N(1)
   do j=1,N(2)
    do k=1,N(3)
     NumCellP = NumCellP + 1
     do l=1,atomspercel
       AtomsTot(NumCellP,l,1) = AtomsCellCoor(l,1) + (i-1)*UnitCellDim(1)
       AtomsTot(NumCellP,l,2) = AtomsCellCoor(l,2) + (j-1)*UnitCellDim(2)
       AtomsTot(NumCellP,l,3) = AtomsCellCoor(l,3) + (k-1)*UnitCellDim(3)
     enddo
    enddo
   enddo
  enddo
  write(*,*) ' What kind of dislocation do you want?'    
  write(*,*) '!++++++++++++++++++++++++++++++++++++++++++!'
  write(*,*) '!++       1 == Edge                      ++!'
  write(*,*) '!++       2 == Screw                     ++!'
  write(*,*) '!++++++++++++++++++++++++++++++++++++++++++!'
  read(*,*) DislocationCaracter
  if (DislocationCaracter == 1) then 
      write(*,*) 'You are about the build an EDGE dislocation in BCC crystal'
      !*** here we need to define the type of geometry we want to use
      write(*,*) '-------------------------------------------'
      write(*,*) '!** What kind of structure do you want? **!'
      write(*,*) '!** 1 == cylinder                       **!'
      write(*,*) '!** 2 == PAD                            **!'
      write(*,*) '!** 3 == Perfect FCC crystal            **!'
      write(*,*) '-------------------------------------------'
      read(*,*) GeometryChoice 
      if (GeometryChoice == 3) then
        write(*,*) 'You are about to build a perfect BCC crystal with edge orientation' 
         write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
         write(*,*) '!++ You are going to create a perfect FCC crystal with orientation : ++!'
         write(*,*) '!++                X -> burgers vector                               ++!'
         write(*,*) '!++                Y -> normal to slip plane                         ++!'
         write(*,*) '!++                Z -> line direction                               ++!' 
         write(*,*) '!++   1 == cubic                                                       ++!'
         write(*,*) '!++   2 == cylinder                                                    ++!'
         write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
         read(*,*) AnswerPerfect       
         if (AnswerPerfect == 1) then
         
           !*** dimension of the box
           dimbox(:,:)=0.0d0
           dimbox(1,2) = N(1)*UnitCellDim(1)
           dimbox(2,2) = N(2)*UnitCellDim(2)
           dimbox(3,2) = N(3)*UnitCellDim(3)

           !*** total number of atoms
            numatomstot = N(1)*N(2)*N(3)*atomspercel
 
           open(1,file='atoms.bcc.perfect.cubic',status='unknown')
           write(1,*) 'Position data for Fe File'
           write(1,*)
           write(1,*) numatomstot, 'atoms'
           write(1,*) ' 1 atom types'
           write(1,*) dimbox(1,1:2),' xlo xhi'
           write(1,*) dimbox(2,1:2),' ylo yhi'
           write(1,*) dimbox(3,1:2),' zlo zhi'
           write(1,*)
           write(1,*) 'Atoms'
           write(1,*)
           m = 1
           iun = 1
           NumCellP=0
           do i=1,N(1)
            do j=1,N(2)
             do k=1,N(3)
               NumCellP = NumCellP+1
               do l=1,6
                 write(1,10) m,iun,AtomsTot(NumCellP,l,1:3)
                 m = m+1
               enddo
             enddo
            enddo
           enddo

           close(1)
           
         else if (AnswerPerfect == 2) then
           write(*,*) 'Not implemented yet'
           deallocate(AtomsCellCoor)
           deallocate(AtomsTot)
           stop
        
         endif
          
      else if (GeometryChoice == 2) then
        write(*,*) 'You are about to build a PAD made of edge dislocation in BCC'
        allocate(hide(ncelltot,atomspercel))
        hide(:,:) = .false.
        !** atoms in the units cell with y <= N(2)/2 and x = N(1) are hiden
        NumCellP = 0
        numdel = 0
        do i=1,N(1)
         do j=1,N(2)
          do k=1,N(3)
            NumCellP = NumCellP +1
            if (i.eq.N(1)) then
             if (j.le.N(2)/2) then
               hide(NumCellP,1:6) = .true.
               numdel = numdel + 6
             endif
            endif
          enddo
         enddo
        enddo  
        !*** dimension of the box
        dimbox(:,:)=0.0d0
        dimbox(1,2) = N(1)*UnitCellDim(1)
        dimbox(2,2) = N(2)*UnitCellDim(2)
        dimbox(3,2) = N(3)*UnitCellDim(3)
        
        !*** total number of atoms
        numatomstot = N(1)*N(2)*N(3)*atomspercel-numdel      

        !*** apply extension of b/2 on the atoms with y <= N(2)/2
        !*** apply a contraction of b/2 on the atoms with y > N(2)/2
        numcell=0
        do i=1,N(1)
         do j=1,N(2)
          do k=1,N(3)
           numcell = numcell + 1
           do l=1,atomspercel
            if (.not. hide(numcell,l)) then
             if (j .le. (N(2)/2)) then
              !*** scale all the x coordinates by +b/2
              xx = AtomsTot(numcell,l,1)
              scaleX = 0.5d0*xx*burgers/((N(1)-1)*UnitCellDim(1))
              AtomsTot(numcell,l,1) = AtomsTot(numcell,l,1) + scaleX
             else if (j .gt. (N(2)/2)) then
              !*** scale all the x coordinates by -b/2
              xx = AtomsTot(numcell,l,1)
              scaleX = -0.5d0*xx*burgers/(N(1)*UnitCellDim(1))
              AtomsTot(numcell,l,1) = AtomsTot(numcell,l,1) + scaleX     
             endif
            endif
           enddo
          enddo
         enddo
        enddo
        
        !*** dimension of the box corrected by the deformation due to the finite dimension
        dimbox(:,:)=0.0d0
        dimbox(1,2) = (N(1)-1)*UnitCellDim(1)+burgers*0.5d0
        dimbox(2,2) = N(2)*UnitCellDim(2)
        dimbox(3,2) = N(3)*UnitCellDim(3)

        open(1,file='atoms.bcc.edge.pad',status='unknown')
        write(1,*) 'Position data for Fe File'
        write(1,*) 
        write(1,*) numatomstot, 'atoms'
        write(1,*) ' 1 atom types'
        write(1,*) dimbox(1,1:2),' xlo xhi'
        write(1,*) dimbox(2,1:2),' ylo yhi'
        write(1,*) dimbox(3,1:2),' zlo zhi'
        write(1,*) 
        write(1,*) 'Atoms'
        write(1,*)
        m = 1
        iun = 1
        NumCellP=0
        do i=1,N(1)
         do j=1,N(2)
          do k=1,N(3)
            NumCellP = NumCellP+1
            do l=1,atomspercel
             if (.not.hide(NumCellP,l)) then
              write(1,10) m,iun,AtomsTot(NumCellP,l,1:3)
              m = m+1
             endif
            enddo
          enddo
         enddo
        enddo

        write(*,*) 'for blocked atoms on the bottom and top surface'
        write(*,*) 'ymin == ',UnitCellDim(2)
        write(*,*) 'ymax == ',(N(2)-1)*UnitCellDim(2) - 0.2d0

        close(1)

        deallocate(hide)
        
        
        
      else if (GeometryChoice == 1) then
        write(*,*) 'you are about to build single edge dislocation in a cylinder in BCC'
        allocate(hide(ncelltot,atomspercel))
        hide(:,:) = .false.

        NumCellP = 0
        numdel = 0
        do i=1,N(1)
         do j=1,N(2)
          do k=1,N(3)
            NumCellP = NumCellP +1
            if (i.eq.((N(1)-1)/2+1)) then
             if (j.le.N(2)/2) then
               hide(NumCellP,1:atomspercel) = .true.
               numdel = numdel + atomspercel
             endif
            endif
          enddo
         enddo
        enddo
      
  
      !*** total number of atoms
      numatomstot = N(1)*N(2)*N(3)*atomspercel-numdel 

      XDislo = UnitCellDim(1)*((N(1)-1)/2+0.5d0) 
      XDislo = ((N(1)-1)/2+0.5)*UnitCellDim(1)
      YDislo = UnitCellDim(2)*(N(2)/2)-0.1
      sr1(1)=XDislo
      sr1(2) = YDislo
      sr1(3)=0.0d0
 
      !** calculate the displacement ux and uy
      factor = burgers/2.0d0/pivalue
      Vburger(1) = burgers;Vburger(2) = 0.0d0;Vburger(3) = 0.0d0
      bx = Vburger(1);by=Vburger(2)
      !*** projection on X
      px(1) = 1.0d0;px(2) = 0.0d0;px(3) = 0.0d0
      !*** projection on Y
      py(1) = 0.0d0;py(2) = 1.0d0;py(3) = 0.0d0
 
      nua = 1.0d0/4.0d0/(1-nu)/2/pivalue
      nub = (1-2.0d0*nu)*nua
      NumCellP = 0
      do i=1,N(1)
       do j=1,N(2)
        do k=1,N(3)
         NumCellP = NumCellP + 1
         do l = 1,atomspercel
          if (.not.hide(NumCellP,l)) then
          xx = AtomsTot(NumCellP,l,1) - XDislo
          yy = AtomsTot(NumCellP,l,2) - YDislo
 
          theta=datan2(-xx,yy)
     
          ux = factor*(theta+xx*yy/2/(1-nu)/(xx*xx+yy*yy))
          uy = -1.0d0*factor*((1-2.0d0*nu)/4/(1-nu)*log(xx*xx+yy*yy)+(xx*xx-yy*yy)/4.0d0/(1-nu)/(xx*xx+yy*yy))
     
          AtomsTot(NumCellP,l,1) = AtomsTot(NumCellP,l,1) + ux
          AtomsTot(NumCellP,l,2) = AtomsTot(NumCellP,l,2) + uy
 
          endif
         enddo
        enddo
       enddo
      enddo

      !***
      write(*,*) 'Dislo position'
      write(*,*) XDislo, YDislo
      write(*,*) 'radius of the cylinder'
      read(*,*) radius

      NumCellP = 0
      do i=1,N(1)
       do j=1,N(2)
        do k=1,N(3)
          NumCellP = NumCellP + 1
          do l=1,atomspercel
           if (.not.hide(NumCellP,l)) then
            !** check if atom is inside of cylinder of radius with axis on XDislo, YDislo
            xx = AtomsTot(NumCellP,l,1) - XDislo
            yy = AtomsTot(NumCellP,l,2) - YDislo      
            distance = sqrt(xx*xx+yy*yy)
            if (distance .gt. radius) then
             hide(NumCellP,l) = .true.
             numdel = numdel + 1
            endif
           endif
          enddo
        enddo
       enddo
      enddo
 
      !*** total number of atoms
      numatomstot = N(1)*N(2)*N(3)*atomspercel-numdel
 
      dimbox(:,:) = 0.0d0
      dimbox(1,2) = N(1)*UnitCellDim(1)
      dimbox(2,2) = N(2)*UnitCellDim(2)
      dimbox(3,2) = N(3)*UnitCellDim(3)

      open(1,file='atoms.bcc.edge.cylinder',status='unknown')
      write(1,*) 'Position data for Fe File'
      write(1,*) 
      write(1,*) numatomstot, 'atoms'
      write(1,*) ' 1 atom types'
      write(1,*) dimbox(1,1:2),' xlo xhi'
      write(1,*) dimbox(2,1:2),' ylo yhi'
      write(1,*) dimbox(3,1:2),' zlo zhi'
      write(1,*) 
      write(1,*) 'Atoms'
      write(1,*)
      m = 1
      iun = 1
      NumCellP=0
      do i=1,N(1)
       do j=1,N(2)
        do k=1,N(3)
          NumCellP = NumCellP+1
          do l=1,atomspercel
           if (.not.hide(NumCellP,l)) then
            write(1,10) m,iun,AtomsTot(NumCellP,l,1:3)
            m = m+1
           endif
          enddo
        enddo
       enddo
      enddo
 
      close(1)
      deallocate(hide)

      endif
        
  else if (DislocationCaracter == 2) then
    write(*,*) 'You are about to build a SCREW dislocation'
    !*** here we need to define the type of geometry we want to use
    write(*,*) '-------------------------------------------'
    write(*,*) '!** What kind of structure do you want? **!'
    write(*,*) '!** 1 == cylinder                       **!'
    write(*,*) '!** 2 == PAD                            **!'
    write(*,*) '-------------------------------------------'
    read(*,*) GeometryChoice 
    if (GeometryChoice == 1) then
      !*** define the position of the dislocation line
      midle(:) = 0.0d0
      midle2(:) = 0.0d0
      epsilon = 0.0001D-01
      midle(1) = N(1)/2*UnitCellDim(1) - epsilon
      midle(3) = N(3)/2*UnitCellDim(3) - 0.0d0   !** Z-
      midle2(3) = N(3)/2*UnitCellDim(3) + 0.0d0  !** Z+
      midle(2) = N(2)/2*UnitCellDim(2)           !** Y
   
      center(2) = N(2)/2.0d0*UnitCellDim(2) !** Y
      center(1) = N(3)/2.0d0*UnitCellDim(3) !** Z
   
      write(*,*) 'dislocation line (X, Y, Z)'
      write(*,*) midle(1),center(2),center(1)

      allocate(hide(ncelltot,atomspercel))
      hide(:,:) = .false.
    
      !*** define the radius of the cylinder
      write(*,*) 'What is the cylinder radius?'
      read(*,*) radius    

      NumCellP = 0
      numhide=0
      do i=1,N(1)
       do j=1,N(2)
        do k=1,N(3)
         NumCellP = NumCellP + 1
         do l=1,6
      
          !* tan**(1)(y/z)
           xx = AtomsTot(NumCellP,l,2)
           yy = AtomsTot(NumCellP,l,3)
        
           !*** screw dislocation displacement
           AtomsTot(NumCellP,l,1) = AtomsTot(NumCellP,l,1) + burgers/2.0d0/pivalue*(datan2((xx-midle(2)),(yy-midle(3))))
!           AtomsTot(NumCellP,l,1) = AtomsTot(NumCellP,l,1) + burgers/4.0d0/pivalue*(datan2((xx-midle(2)),(yy-midle2(3))))

           !*** check if atoms inside the cylinder
           xx = xx - center(2)
           yy = yy - center(1)
           distance = sqrt(xx*xx+yy*yy)
           if (distance .gt. radius) then
             numhide = numhide + 1
             hide(NumCellP,l) = .true.
           endif
   
         enddo
        enddo
       enddo
      enddo   
   
      dimbox(:,:)=0.0d0
      dimbox(1,2) = N(1)*UnitCellDim(1)
      dimbox(2,2) = N(2)*UnitCellDim(2)
      dimbox(3,2) = N(3)*UnitCellDim(3)
   
      open(1,file='atoms.bcc.screw.cylinder',status='unknown')
   
      write(1,*) 'Position data for Fe File'
      write(1,*) 
      write(1,*) 6*ncelltot - numhide, ' atoms'
      write(1,*) ' 1 atom types'
      write(1,*) dimbox(1,1:2),' xlo xhi'
      write(1,*) dimbox(2,1:2),' ylo yhi'
      write(1,*) dimbox(3,1:2),' zlo zhi'
      write(1,*)
      write(1,*) 'Atoms'
      write(1,*)
      m=1
      iun=1
      NumCellP=0
      do i=1,N(1)
       do j=1,N(2)
        do k=1,N(3)
         NumCellP = NumCellP+1
         do l=1,6
          if (.not.hide(NumCellP,l)) then
            write(1,10) m, iun, AtomsTot(NumCellP,l,1:3)
            m  = m+1
          endif
         enddo
        enddo
       enddo
      enddo
   
      close(1)
      deallocate(hide)
    
    else if (GeometryChoice == 2) then
      write(*,*) 'you are about to build an infinite screw dislocation in a cylinder in bcc'
      write(*,*) 'NOT IMPLEMENTED YET'
    
    endif
      
  
  endif

  deallocate(AtomsCellCoor)
  deallocate(AtomsTot)
  
  
case(3)
  write(*,*) 'you are about to construct a dislocation in HCP'
  !*** general parameters
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(*,*) '++   what kind of material do you want to consider      ++'
  write(*,*) '++ 1 == Mg (Sun et al, 2006)                            ++'
  write(*,*) '++ 2 == Mg (Liu et al, 1996)                            ++'
  write(*,*) '++ 3 == Mg (experimental value)                         ++'
  write(*,*) '++ 4 == open value                                      ++'
  write(*,*) '=========================================================='
  read(*,*) HCPMaterial
  if (HCPMaterial == 1) then
    !*** for Sun et al potential (2006)
    avalue = 3.184d0
    acvalue = 1.628d0*avalue
  else if (HCPMaterial == 2) then
    !*** for Liu et al potential (1996)
    avalue = 3.206d0
    acvalue = 1.623d0*avalue
  else if (HCPMaterial == 3) then
    !*** experimental value for magnesium
    avalue = 3.203d0
    acvalue = 0.994d0*sqrt(8.0d0/3.0d0)*avalue
  else if (HCPMaterial == 4) then
    write(*,*) 'avalue ??'
    read(*,*) avalue
    acvalue=sqrt(8.0d0/3.0d0)*avalue
  endif
  write(*,*) HCPMaterial
  
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  write(*,*) '++    Which Orientation do you want to consider?        ++'
  write(*,*) '++    1 == Basal                                        ++'
  write(*,*) '++    2 == prismatic                                    ++'
  write(*,*) '++    3 == pyramidal-a                                  ++'
  write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  read(*,*)HCPOrientation 
  write(*,*) HCPOrientation
  if (HCPOrientation == 1) then
    !** basal orientation
    burgers = avalue
    write(*,*) ' What kind of dislocation do you want?'    
    write(*,*) '!++++++++++++++++++++++++++++++++++++++++++!'
    write(*,*) '!++       1 == Edge                      ++!'
    write(*,*) '!++       2 == Screw                     ++!'
    write(*,*) '!++++++++++++++++++++++++++++++++++++++++++!'
    read(*,*) DislocationCaracter
    if (DislocationCaracter == 1) then

    !*** atoms for one unit cell
    r(1) = 2.0d0/3.0d0*avalue*dcos(pisurtrois)+1.0d0/3.0d0*avalue*dcos(pisurtrois)
    r(2) = 0.5*acvalue
    r(3) = 2.0d0/3.0d0*avalue*dsin(pisurtrois)-1.0d0/3.0d0*avalue*dsin(pisurtrois)
    !*** number of atoms per unit cell
    atomspercel = 4
    allocate(AtomsCellCoor(atomspercel, 3))
    !** Atoms #1
    AtomsCellCoor(1,1) = 0.0d0
    AtomsCellCoor(1,2) = 0.0d0
    AtomsCellCoor(1,3) = 0.0d0
    !** Atoms #2
    AtomsCellCoor(2,1) = avalue*sin(pivalue/6.0d0)
    AtomsCellCoor(2,2) = 0.0d0
    AtomsCellCoor(2,3) = avalue*cos(pivalue/6.0d0)
    !** Atoms #3
    AtomsCellCoor(3,1) = AtomsCellCoor(1,1) + r(1)
    AtomsCellCoor(3,2) = AtomsCellCoor(1,2) + r(2)
    AtomsCellCoor(3,3) = AtomsCellCoor(1,3) + r(3)
    !** Atoms #4
    AtomsCellCoor(4,1) = AtomsCellCoor(2,1) + r(1)
    AtomsCellCoor(4,2) = AtomsCellCoor(2,2) + r(2)
    AtomsCellCoor(4,3) = AtomsCellCoor(2,3) + r(3)
    !*** Size of the unit cell
    UnitCellDim(1) = avalue; UnitCellDim(2) = acvalue; UnitCellDim(3) = 2*avalue*cos(pivalue/6.0d0)
    !** dimension of the reference unit cell
    write(*,*) 'dimension reference unit cell :',UnitCellDim
    ! N(1) : must be an odd number
    ! N(2) : must be an even number
    ! N(3) : no constrain
    write(*,*) 'number of unit cells along X, Y and Z'
    read(*,*) N(1),N(2),N(3)

    !*** total number of unit cell
    ncelltot = N(1)*N(2)*N(3)
    allocate(AtomsTot(ncelltot,atomspercel,3))
    !** generate the atomic positions
    NumCellP = 0
    do i = 1, N(1)
      do j = 1, N(2)
        do k = 1, N(3)
          NumCellP = NumCellP + 1
      
          do l=1,atomspercel
            AtomsTot(NumCellP,l,1) = AtomsCellCoor(l,1) + (i-1)*UnitCellDim(1)
            AtomsTot(NumCellP,l,2) = AtomsCellCoor(l,2) + (j-1)*UnitCellDim(2)
            AtomsTot(NumCellP,l,3) = AtomsCellCoor(l,3) + (k-1)*UnitCellDim(3)
          enddo

        enddo
      enddo
    enddo
    

        write(*,*) 'You are about the build an EDGE dislocation in FCC crystal'
        !*** here we need to define the type of geometry we want to use
        write(*,*) '-------------------------------------------'
        write(*,*) '!** What kind of structure do you want? **!'
        write(*,*) '!** 1 == cylinder                       **!'
        write(*,*) '!** 2 == PAD                            **!'
        write(*,*) '!** 3 == Perfect FCC crystal            **!'
        write(*,*) '-------------------------------------------'
        read(*,*) GeometryChoice 
      
        if (GeometryChoice == 3) then
           !*** dimension of the simulation box
           dimbox(:,:) = 0.0d0
           dimbox(1,2) = N(1)*UnitCellDim(1)
           dimbox(2,2) = N(2)*UnitCellDim(2)
           dimbox(3,2) = N(3)*UnitCellDim(3)
           write(*,*) 'you are about to build a perfect hcp structure oriented for edge basal slip'
           write(*,*) '-------------------------------------------'
           write(*,*) '!** What kind of structure do you want? **!'
           write(*,*) '!** 1 == cylinder                       **!'
           write(*,*) '!** 2 == PAD                            **!'
           write(*,*) '-------------------------------------------'         
           read(*,*) AnswerPerfect
           if (AnswerPerfect == 1) then
           !*** cylinder
             allocate(hide(ncelltot, atomspercel))
             hide(:,:)=.false.
             center(1) = N(1)*0.5d0*UnitCellDim(1)
             center(2) = 0.5d0*N(2)*UnitCellDim(3)
             write(*,*) 'center of the cylinder'
             write(*,*) center(1:2)
             write(*,*) 'what is the radius of the cylinder'
             read(*,*) radius
             
             NumCellP=0
             numdel = 0
             do i=1,N(1)
              do j=1,N(2)
               do k=1,N(3)
                NumCellP = NumCellP + 1
                do l=1,atomspercel
                 xx = AtomsTot(NumCellP,l,1) - center(1)
                 yy = AtomsTot(NumCellP,l,2) - center(2)
                 
                 distance = sqrt(xx*xx+yy*yy)
                 
                 if (distance .gt. radius) then
                  hide(NumCellP,l) = .true.
                  numdel = numdel + 1
                 endif
                enddo
               enddo
              enddo
             enddo
             
             numatomstot = N(1)*N(2)*N(3)*atomspercel - numdel
             open(1,file='atoms.hcp.basal.perfect.cylinder',status='unknown')
             write(1,*) 'Position data for Mg File'
             write(1,*) 
             write(1,*) numatomstot, 'atoms'
             write(1,*) ' 1 atom types'
             write(1,*) dimbox(1,1:2),' xlo xhi'
             write(1,*) dimbox(2,1:2),' ylo yhi'
             write(1,*) dimbox(3,1:2),' zlo zhi'
             write(1,*) 
             write(1,*) 'Atoms'
             write(1,*)
             m = 1
             iun = 1
             NumCellP=0
             do i=1,N(1)
              do j=1,N(2)
               do k=1,N(3)
                 NumCellP = NumCellP+1
                 do l=1,atomspercel
                  if (.not.hide(NumCellP,l)) then
                   write(1,10) m,iun,AtomsTot(NumCellP,l,1:3)
                   m = m+1
                  endif
                 enddo
               enddo
              enddo
             enddo     
             
             close(1)             
             
             deallocate(hide)
           
           else if (AnswerPerfect == 2) then

             !*** number of atoms
             numatomstot = N(1)*N(2)*N(3)*atomspercel
             open(1,file='atoms.hcp.basal.perfect.cubic',status='unknown')
             write(1,*) 'Position data for Mg File'
             write(1,*) 
             write(1,*) numatomstot, 'atoms'
             write(1,*) ' 1 atom types'
             write(1,*) dimbox(1,1:2),' xlo xhi'
             write(1,*) dimbox(2,1:2),' ylo yhi'
             write(1,*) dimbox(3,1:2),' zlo zhi'
             write(1,*) 
             write(1,*) 'Atoms'
             write(1,*)
             m = 1
             iun = 1
             NumCellP=0
             do i=1,N(1)
              do j=1,N(2)
               do k=1,N(3)
                 NumCellP = NumCellP+1
                 do l=1,atomspercel
                   write(1,10) m,iun,AtomsTot(NumCellP,l,1:3)
                   m = m+1
                 enddo
               enddo
              enddo
             enddo     
             
             close(1)
                   
           endif
         
        else if (GeometryChoice == 2) then
         write(*,*) 'you are about to build a PA of edge dislocation for basal slip'
         allocate(hide(ncelltot, atomspercel))
         hide(:,:) = .false.
         
         !** atoms in the unit cells with y <= N(2)/2 and x = N(1) are hiden
          NumCellP = 0
          numdel = 0
          do i=1,N(1)
           do j=1,N(2)
            do k=1,N(3)
             NumCellP = NumCellP + 1
             if (i == N(1)) then
              if (j .le. N(2)/2) then
               hide(NumCellP,1:4) = .true.
               numdel = numdel + 4
              endif
             endif
            enddo
           enddo
          enddo 

          !*** dimension of the box
          dimbox(:,:) = 0.0d0
          dimbox(1,2) = N(1)*UnitCellDim(1)
          dimbox(2,2) = N(2)*UnitCellDim(2)
          dimbox(3,2) = N(3)*UnitCellDim(3)
          
          !*** number of atoms
          numatomstot = N(1)*N(2)*N(3)*atomspercel-numdel
          
          !*** apply extension of b/2 on the atoms with y <= N(2)/2
          !*** apply a contraction of b/2 on the atoms with y > N(2)/2
          numcell = 0
          do i=1,N(1)
           do j=1,N(2)
            do k=1,N(3)
             numcell = numcell + 1
             do l=1,4
              if (.not. hide(numcell,l)) then
                if (j .le. (N(2)/2)) then
                 !*** scale all the x coordinates by +b/2
                 xx = AtomsTot(numcell,l,1)
                 scaleX = 0.5d0*xx*burgers/((N(1)-1)*UnitCellDim(1))
                 AtomsTot(numcell,l,1) = AtomsTot(numcell,l,1) + scaleX
                else if (j .gt. (N(2)/2)) then
                 !*** scale all the x coordinates by -b/2
                 xx = AtomsTot(numcell,l,1)
                 scaleX = -0.5d0*xx*burgers/(N(1)*UnitCellDim(1))
                 AtomsTot(numcell,l,1) = AtomsTot(numcell,l,1) + scaleX     
                endif
              endif
             enddo
            enddo
           enddo
          enddo
            
          !*** dimension of the box
          dimbox(:,:) = 0.0d0
          dimbox(1,2) = (N(1)-1)*UnitCellDim(1)+burgers*0.5d0
          dimbox(2,2) = N(2)*UnitCellDim(2)
          dimbox(3,2) = N(3)*UnitCellDim(3)
            
          open(1,file='atoms.hcp.basal.edge.pad',status='unknown')
          write(1,*) 'Position data for Mg File'
          write(1,*) 
          write(1,*) numatomstot, 'atoms'
          write(1,*) ' 1 atom types'
          write(1,*) dimbox(1,1:2),' xlo xhi'
          write(1,*) dimbox(2,1:2),' ylo yhi'
          write(1,*) dimbox(3,1:2),' zlo zhi'
          write(1,*) 
          write(1,*) 'Atoms'
          write(1,*)
          m = 1
          iun = 1
          NumCellP=0
          do i=1,N(1)
           do j=1,N(2)
            do k=1,N(3)
              NumCellP = NumCellP+1
              do l=1,4
               if (.not.hide(NumCellP,l)) then
                write(1,10) m,iun,AtomsTot(NumCellP,l,1:3)
                m = m+1
               endif
              enddo
            enddo
           enddo
          enddo
          write(*,*) 'for blocked atoms on the bottom and top surface'
          write(*,*) 'ymin == ',UnitCellDim(2)-0.1d0
          write(*,*) 'ymax == ',(N(2)-1)*UnitCellDim(2) - 0.2d0

         
         deallocate(hide)
         
        else if (GeometryChoice == 1) then
          write(*,*) 'you are about to build an infinite edge dislocation in a cylinder for basal slip'
          allocate(hide(ncelltot,atomspercel))
          hide(:,:) = .false.

          NumCellP = 0
          numdel = 0
          do i=1,N(1)
           do j=1,N(2)
            do k=1,N(3)
              NumCellP = NumCellP +1
              if (i.eq.((N(1)-1)/2+1)) then
               if (j.le.N(2)/2) then
                 hide(NumCellP,1:atomspercel) = .true.
                 numdel = numdel + atomspercel
               endif
              endif
            enddo
           enddo
          enddo

          numatomstot = N(1)*N(2)*N(3)*atomspercel-numdel

          XDislo = UnitCellDim(1)*((N(1)-1)/2+0.5d0) 
          XDislo = ((N(1)-1)/2+0.5)*UnitCellDim(1)
          YDislo = UnitCellDim(2)*(N(2)/2)-0.1
          sr1(1) = XDislo
          sr1(2) = YDislo
          sr1(3)=0.0d0

         !** calculate the displacement ux and uy
         factor = burgers/2.0d0/pivalue
         Vburger(1) = burgers;Vburger(2) = 0.0d0;Vburger(3) = 0.0d0
         bx = Vburger(1);by=Vburger(2)
         !*** projection on X
         px(1) = 1.0d0;px(2) = 0.0d0;px(3) = 0.0d0
         !*** projection on Y
         py(1) = 0.0d0;py(2) = 1.0d0;py(3) = 0.0d0
 
         nua = 1.0d0/4.0d0/(1-0.33d0)/2/pivalue
         nub = (1-2.0d0*0.33d0)*nua
         NumCellP = 0
         do i=1,N(1)
          do j=1,N(2)
           do k=1,N(3)
            NumCellP = NumCellP + 1
            do l = 1,4
             if (.not.hide(NumCellP,l)) then
             xx = AtomsTot(NumCellP,l,1) - XDislo
             yy = AtomsTot(NumCellP,l,2) - YDislo
 
             theta=datan2(-xx,yy)
     
             ux = factor*(theta+xx*yy/2/(1-0.33d0)/(xx*xx+yy*yy))
             uy = -1.0d0*factor*((1-2.0d0*0.33d0)/4/(1-0.33d0)*log(xx*xx+yy*yy)+(xx*xx-yy*yy)/4.0d0/(1-0.33d0)/(xx*xx+yy*yy))
     
             AtomsTot(NumCellP,l,1) = AtomsTot(NumCellP,l,1) + ux
             AtomsTot(NumCellP,l,2) = AtomsTot(NumCellP,l,2) + uy
 
             endif
            enddo
           enddo
          enddo
         enddo

         !***
         write(*,*) 'Dislo position'
         write(*,*) XDislo, YDislo
         write(*,*) 'radius of the cylinder'
         read(*,*) radius
 
         NumCellP = 0
         do i=1,N(1)
          do j=1,N(2)
           do k=1,N(3)
             NumCellP = NumCellP + 1
             do l=1,4
              if (.not.hide(NumCellP,l)) then
               !** check if atom is inside of cylinder of radius with axis on XDislo, YDislo
               xx = AtomsTot(NumCellP,l,1) - XDislo
               yy = AtomsTot(NumCellP,l,2) - YDislo      
               distance = sqrt(xx*xx+yy*yy)
               if (distance .gt. radius) then
                hide(NumCellP,l) = .true.
                numdel = numdel + 1
               endif
              endif
             enddo
           enddo
          enddo
         enddo
  
         numatomstot = N(1)*N(2)*N(3)*atomspercel-numdel
 
         dimbox(:,:)=0.0d0
         dimbox(1,2) = N(1)*UnitCellDim(1)
         dimbox(2,2) = N(2)*UnitCellDim(2)
         dimbox(3,2) = N(3)*UnitCellDim(3)

         open(1,file='atoms.hcp.basal.edge.cylinder',status='unknown')
         write(1,*) 'Position data for Mg File'
         write(1,*) 
         write(1,*) numatomstot, 'atoms'
         write(1,*) ' 1 atom types'
         write(1,*) dimbox(1,1:2),' xlo xhi'
         write(1,*) dimbox(2,1:2),' ylo yhi'
         write(1,*) dimbox(3,1:2),' zlo zhi'
         write(1,*) 
         write(1,*) 'Atoms'
         write(1,*)
         m = 1
         iun = 1
         NumCellP=0
         do i=1,N(1)
          do j=1,N(2)
           do k=1,N(3)
             NumCellP = NumCellP+1
             do l=1,4
              if (.not.hide(NumCellP,l)) then
               write(1,10) m,iun,AtomsTot(NumCellP,l,1:3)
               m = m+1
              endif
             enddo
           enddo
          enddo
         enddo
 
         close(1)

         deallocate(hide)
                  
        endif
        
        deallocate(AtomsCellCoor)
    else if (DislocationCaracter == 2) then
      write(*,*) 'you are about to build a screw dislocation in hcp oriented for basal slip'
      !*** atoms for one unit cell
      atomspercel = 4
      allocate(AtomsCellCoor(atomspercel, 3))
      !* atom #1
      AtomsCellCoor(1,2) = 0.0d0                             
      AtomsCellCoor(1,3) = 0.0d0          
      AtomsCellCoor(1,1) = 0.0d0
      !* atom #2
      AtomsCellCoor(2,2) = avalue*cos(pioversix)             
      AtomsCellCoor(2,3) = 0.0d0          
      AtomsCellCoor(2,1) = avalue*sin(pioversix)
      !* atom #3
      AtomsCellCoor(3,2) = 1.0d0/3.0d0*avalue*cos(pioversix) 
      AtomsCellCoor(3,3) = 0.5d0*acvalue  
      AtomsCellCoor(3,1) = avalue*sin(pioversix)
      !* atom #4
      AtomsCellCoor(4,2) = 4.0d0/3.0d0*avalue*cos(pioversix) 
      AtomsCellCoor(4,3) = 0.5d0*acvalue  
      AtomsCellCoor(4,1) = 2*avalue*sin(pioversix)

      !*** dimension of the simulation cell
      UnitCellDim(2) = 2*avalue*cos(pioversix); UnitCellDim(3) = acvalue; UnitCellDim(1) = avalue
      write(*,*) '--------------------------------'
      write(*,*) 'dimension of the unit along X (b), Y (d), Z (n)'
      write(*,*) UnitCellDim(1:3)
      ! N(1) : must be an odd number (line of dislocation oriented along N(1) with periodicity
      ! N(2) : must be an even number
      ! N(3) : no constrain
      write(*,*) 'number of unit cells along X(b), Y(n) and Z(d)'
      read(*,*) N(1),N(2),N(3)
      !*** total number of cell
      ncelltot = N(1)*N(2)*N(3)
      !*** location of the dislocation
      midle(1) = N(1)/2*UnitCellDim(1) - epsilon
      midle(2) = N(2)/2*UnitCellDim(2) - 1.0d0
      midle2(2) = N(2)/2*UnitCellDim(2) + 1.0d0
      midle(3) = N(3)/2*UnitCellDim(3) - 0.25*acvalue
      center(1) = midle(3)
      center(2) = N(2)/2*UnitCellDim(2)

      allocate(AtomsTot(ncelltot,atomspercel,3))
      allocate(hide(ncelltot,atomspercel))
      
      write(*,*) '-------------------------------------------'
      write(*,*) '!** What kind of structure do you want? **!'
      write(*,*) '!** 1 == cylinder                       **!'
      write(*,*) '!** 2 == PAD                            **!'
      write(*,*) '-------------------------------------------'
      read(*,*) GeometryChoice
      
      if (GeometryChoice == 1) then
      write(*,*) center(1),center(2)
      write(*,*) 'what is the radius of the cylinder?'
      read(*,*) radius
      
      hide(:,:)=.false.
      numhide = 0
      NumCellP = 0
      do i = 1, N(1)
        do j = 1, N(2)
          do k = 1, N(3)
            NumCellP = NumCellP + 1
      
            do l=1,atomspercel
              AtomsTot(NumCellP,l,1) = AtomsCellCoor(l,1) + (i-1)*UnitCellDim(1)
              AtomsTot(NumCellP,l,2) = AtomsCellCoor(l,2) + (j-1)*UnitCellDim(2)
              AtomsTot(NumCellP,l,3) = AtomsCellCoor(l,3) + (k-1)*UnitCellDim(3)
      
              xx = AtomsTot(NumCellP,l,2)
              yy = AtomsTot(NumCellP,l,3)        

              AtomsTot(NumCellP,l,1) = AtomsTot(NumCellP,l,1) + burgers/4/pivalue*(datan2((yy-midle(3)),(xx-midle(2))))
              AtomsTot(NumCellP,l,1) = AtomsTot(NumCellP,l,1) + burgers/4/pivalue*(datan2((yy-midle(3)),(xx-midle2(2))))

              !* check if atoms inside the cylinder
              xx = xx - center(2)
              yy = yy - center(1)
              distance=sqrt(xx*xx+yy*yy)
              if (distance .gt. radius) then
               numhide = numhide + 1
               hide(NumCellP,l) = .true.
              endif

            enddo
          enddo
        enddo
      enddo

      else if (GeometryChoice == 2) then
        write(*,*) 'screw dislocation for basal slip in a PAD'
        numhide = 0
        NumCellP=0
        do i = 1, N(1)
          do j = 1, N(2)
            do k = 1, N(3)
              NumCellP = NumCellP + 1
      
              do l=1,4
                AtomsTot(NumCellP,l,1) = AtomsCellCoor(l,1) + (i-1)*UnitCellDim(1)
                AtomsTot(NumCellP,l,2) = AtomsCellCoor(l,2) + (j-1)*UnitCellDim(2)
                AtomsTot(NumCellP,l,3) = AtomsCellCoor(l,3) + (k-1)*UnitCellDim(3)

                xx = AtomsTot(NumCellP,l,2)
                yy = AtomsTot(NumCellP,l,3)

                AtomsTot(NumCellP,l,1) = AtomsTot(NumCellP,l,1) + burgers/4.0d0/pivalue*(datan2((yy-midle(3)),(xx-midle(2))))
                AtomsTot(NumCellP,l,1) = AtomsTot(NumCellP,l,1) + burgers/4.0d0/pivalue*(datan2((yy-midle(3)),(xx-midle2(2))))

              enddo


          
              if ((i.eq.1) .and. (j.eq.1) .and. (k.eq.1)) then
               
               write(*,*) 'to set up the boundary conditions on the bottom'
               write(*,*) AtomsTot(NumCellP, 3, 3)
              
              elseif ((i.eq.1) .and. (j.eq.1) .and. (k.eq.N(3))) then
               write(*,*) 'to set up the boundary conditions on the top'
               write(*,*) AtomsTot(NumCellP, 2, 3)       
              endif

            enddo
          enddo
        enddo
      endif

      dimbox(1,1) = 0.0d0; dimbox(2,1) = 0.0d0; dimbox(3,1) = 0.0d0
      dimbox(1,2) = N(1) * UnitCellDim(1)
      dimbox(2,2) = N(2) * UnitCellDim(2)
      dimbox(3,2) = N(3) * UnitCellDim(3) 
      
      numatomstot = N(1)*N(2)*N(3)*atomspercel - numhide  
      if (GeometryChoice == 2   ) then
        open(1,file='atoms.hcp.basal.screw.pad',status='unknown')
      else if (GeometryChoice == 1) then
        open(1,file='atoms.hcp.basal.screw.cylinder',status='unknown')
      endif
      write(1,*) 'Position data for Mg File'
      write(1,*) 
      write(1,*) numatomstot, 'atoms'
      write(1,*) ' 1 atom types'
      write(1,*) dimbox(1,1:2),' xlo xhi'
      write(1,*) dimbox(2,1:2),' ylo yhi'
      write(1,*) dimbox(3,1:2),' zlo zhi'
      if (GeometryChoice == 2) then
        write(1,*) 0.5*burgers, rzero, rzero, 'xy xz yz'
      endif
      write(1,*) 
      write(1,*) 'Atoms'
      write(1,*)
      m = 1
      iun = 1
      NumCellP = 0
      do i = 1, N(1)
        do j = 1, N(2)
          do k = 1, N(3)
            NumCellP = NumCellP + 1
            do l=1,atomspercel
              if (.not.hide(NumCellP,l)) then
                write(1,10) m, iun, AtomsTot(NumCellP, l, 1:3)
                m = m + 1
              endif
            enddo
          enddo
        enddo
      enddo

      close(1)

      deallocate(AtomsCellCoor)
      deallocate(hide)
      deallocate(AtomsTot)      

    endif        
      
  else if (HCPOrientation == 2) then
    !** prismatic orientation
    write(*,*) 'not implemented yet'

  else if (HCPOrientation == 3) then
    !** pyramidal-a orientation
    write(*,*) 'not implemented yet'
  
  endif
  
  
  case(4)

  write(*,*) 'diamond structure'
  write(*,*) 'what is the lattice parameter a (cBN = default value = 1)?'
  read(*,*) answerlattice
  if (answerlattice == 1) then
   avalue = 3.609d0
  else
   avalue =answerlattice
  endif
  
  !*** magnitude of the Burgers vector
  burgers = avalue/R_2
  
  !*** atoms in unit cell
  atomspercel = 12
  allocate(AtomsCellCoor(atomspercel,3))
  allocate(itype(atomspercel))
  !*** position of the atoms in the reference unit cell
!  # atoms from 1 to 6 == B
!  !## atoms #1
!  AtomsCellCoor(1,1) = avalue/2/R_2
!  AtomsCellCoor(1,2) = 0.0
!  AtomsCellCoor(1,3) = avalue/2/R_6*3
!  itype(1) = 1
!  !## atoms #2
!  AtomsCellCoor(2,1) = avalue/2/R_2
!  AtomsCellCoor(2,2) = avalue/R_3
!  AtomsCellCoor(2,3) = avalue*5*R_6/12
!  itype(2) = 1
!  !## atoms #3
!!  AtomsCellCoor(5,1) = avalue/2/R_2
!  AtomsCellCoor(5,2) = avalue*2/R_3
!  AtomsCellCoor(5,3) = avalue/2/R_6
!  itype(5) = 1
!  !## atoms #4
!  AtomsCellCoor(3,1) = 0.0d0
!  AtomsCellCoor(3,2) = 0.0d0
!  AtomsCellCoor(3,3) = 0.0d0
!  itype(3) = 1
!  !## atoms #5
!  AtomsCellCoor(4,1) = 0.0d0
!  AtomsCellCoor(4,2) = avalue/R_3
!  AtomsCellCoor(4,3) = avalue/R_6
!  itype(4) = 1
!  !## atoms #6
!  AtomsCellCoor(6,1) = 0.0d0
!  AtomsCellCoor(6,2) = avalue*2/R_3
!  AtomsCellCoor(6,3) = avalue*2/R_6
!  itype(6) = 1
!! atoms from 7 to 12 = N
!  !## atoms #7
!  AtomsCellCoor(7,1) = AtomsCellCoor(1,1) + 0.25d0*avalue
!  AtomsCellCoor(7,2) = AtomsCellCoor(1,2) + 0.25d0*avalue
!  AtomsCellCoor(7,3) = AtomsCellCoor(1,3) + 0.25d0*avalue
!  itype(7) = 2
!  !## atoms #8
!  AtomsCellCoor(8,1) = AtomsCellCoor(2,1) + 0.25d0*avalue
!  AtomsCellCoor(8,2) = AtomsCellCoor(2,2) + 0.25d0*avalue
!  AtomsCellCoor(8,3) = AtomsCellCoor(2,3) + 0.25d0*avalue
!  itype(8) = 2
!  !## atoms #9
!  AtomsCellCoor(9,1) = AtomsCellCoor(3,1) + 0.25d0*avalue
!  AtomsCellCoor(9,2) = AtomsCellCoor(3,2) + 0.25d0*avalue
!  AtomsCellCoor(9,3) = AtomsCellCoor(3,3) + 0.25d0*avalue
!  itype(9) = 2
!  !## atoms #10
!  AtomsCellCoor(10,1) = AtomsCellCoor(4,1) + 0.25d0*avalue
!  AtomsCellCoor(10,2) = AtomsCellCoor(4,2) + 0.25d0*avalue
!  AtomsCellCoor(10,3) = AtomsCellCoor(4,3) + 0.25d0*avalue
!  itype(10) = 2
!  !## atoms #11
!  AtomsCellCoor(11,1) = AtomsCellCoor(5,1) + 0.25d0*avalue
!  AtomsCellCoor(11,2) = AtomsCellCoor(5,2) + 0.25d0*avalue
!  AtomsCellCoor(11,3) = AtomsCellCoor(5,3) + 0.25d0*avalue
!  itype(11) = 2
!  !## atoms #12
!  AtomsCellCoor(12,1) = AtomsCellCoor(6,1) + 0.25d0*avalue
!  AtomsCellCoor(12,2) = AtomsCellCoor(6,2) + 0.25d0*avalue
!  AtomsCellCoor(12,3) = AtomsCellCoor(6,3) + 0.25d0*avalue
!  itype(12) = 2
!




  !## atoms #1
  AtomsCellCoor(1,1) = avalue/R_2/2.0d0
  AtomsCellCoor(1,2) = avalue*R_3*11.0d0/12.0d0
  AtomsCellCoor(1,3) = 3.0d0*avalue/R_6/6.0d0
  itype(1) = 1
  !## atoms #2
  AtomsCellCoor(2,1) = avalue/R_2/2.0d0
  AtomsCellCoor(2,2) = avalue*R_3*2.0d0/3.0d0
  AtomsCellCoor(2,3) = 3.0d0*avalue/R_6/6.0d0
  itype(2) = 2
  !## atoms #3
  AtomsCellCoor(3,1) = avalue/R_2*0.0d0
  AtomsCellCoor(3,2) = avalue*R_3*11.0d0/12.0d0
  AtomsCellCoor(3,3) = 3.0d0*avalue/R_6*2.0d0/3.0d0
  itype(3) = 1
  !## atoms #4
  AtomsCellCoor(4,1) = avalue/R_2*0.0d0
  AtomsCellCoor(4,2) = avalue*R_3*7.0d0/12.0d0
  AtomsCellCoor(5,3) = 3.0d0*avalue/R_6*0.0d0
  itype(4) = 1  
  !## atoms #5
  AtomsCellCoor(5,1) = avalue/R_2/2.0d0
  AtomsCellCoor(5,2) = avalue*R_3*7.0d0/12.0d0
  AtomsCellCoor(5,3) = 3.0d0*avalue/R_6/2.0d0
  itype(5) = 1  
  !## atoms #6
  AtomsCellCoor(6,1) = avalue/R_2*0.0d0
  AtomsCellCoor(6,2) = avalue*R_3*2.0d0/3.0d0
  AtomsCellCoor(6,3) = 3.0d0*avalue/R_6*2.0d0/3.0d0
  itype(6) = 2
  !## atoms #7
  AtomsCellCoor(7,1) = avalue/R_2*0.0d0
  AtomsCellCoor(7,2) = avalue*R_3/3.0d0
  AtomsCellCoor(7,3) = 3.0d0*avalue/R_6*0.0d0
  itype(7) = 2
  !## atoms #8
  AtomsCellCoor(8,1) = avalue/R_2/2.0d0
  AtomsCellCoor(8,2) = avalue*R_3/3.0d0
  AtomsCellCoor(8,3) = 3.0d0*avalue/R_6/2.0d0
  itype(8) = 2
 !## atoms #9
  AtomsCellCoor(9,1) = avalue/R_2*0.0d0
  AtomsCellCoor(9,2) = avalue*R_3/4.0d0
  AtomsCellCoor(9,3) = 3.0d0*avalue/R_6/3.0d0
  itype(9) = 1
  !## atoms #10
  AtomsCellCoor(10,1) = avalue/R_2/2.0d0
  AtomsCellCoor(10,2) = avalue*R_3/4.0d0
  AtomsCellCoor(10,3) = 3.0d0*avalue/R_6*5.0d0/6.0d0
  itype(10) = 1  
  !## atoms #11
  AtomsCellCoor(11,1) = avalue/R_2*0.0d0
  AtomsCellCoor(11,2) = avalue*R_3*0.0d0
  AtomsCellCoor(11,3) = 3.0d0*avalue/R_6/3.0d0
  itype(11) = 2
  !## atoms #12
  AtomsCellCoor(12,1) = avalue/R_2/2.0d0
  AtomsCellCoor(12,2) = avalue*R_3*0.0d0
  AtomsCellCoor(12,3) = 3.0d0*avalue/R_6*5.0d0/6.0d0
  itype(12) = 2

  !*** size of the unit cell
  UnitCellDim(1) = avalue/R_2; UnitCellDim(2) = avalue*R_3; UnitCellDim(3) = 3/R_6*avalue
  
  !** dimension of the reference unit cell
  write(*,*) 'dimension reference unit cell :',UnitCellDim
  ! N(1) : must be an odd number
  ! N(2) : must be an even number
  ! N(3) : no constrain
  write(*,*) 'number of unit cells along X, Y and Z (oriented for edge dislocation x == b; y == n)'
  read(*,*) N(1),N(2),N(3)

  !*** total number of unit cell
  ncelltot = N(1)*N(2)*N(3)

  allocate(AtomsTot(ncelltot,atomspercel,3))
  allocate(itypenumcell(ncelltot,atomspercel))
  !** generate the atomic positions
  NumCellP = 0
  do i=1,N(1)
   do j=1,N(2)
    do k=1,N(3)
     NumCellP = NumCellP + 1
     do l=1,atomspercel
       AtomsTot(NumCellP,l,1) = AtomsCellCoor(l,1) + (i-1)*UnitCellDim(1)
       AtomsTot(NumCellP,l,2) = AtomsCellCoor(l,2) + (j-1)*UnitCellDim(2)
       AtomsTot(NumCellP,l,3) = AtomsCellCoor(l,3) + (k-1)*UnitCellDim(3)
       itypenumcell(NumCellP,l) = itype(l)
     enddo
    enddo
   enddo
  enddo
  
  write(*,*) ' What kind of dislocation do you want?'    
  write(*,*) '!++++++++++++++++++++++++++++++++++++++++++!'
  write(*,*) '!++       1 == Edge                      ++!'
  write(*,*) '!++       2 == Screw                     ++!'
  write(*,*) '!++++++++++++++++++++++++++++++++++++++++++!'
  read(*,*) DislocationCaracter
  
  if (DislocationCaracter == 1) then
      write(*,*) 'You are about the build an EDGE dislocation in diamond crystal structure'
      !*** here we need to define the type of geometry we want to use
      write(*,*) '----------------------------------------------'
      write(*,*) '!** What kind of structure do you want?    **!'
      write(*,*) '!** 1 == cylinder                          **!'
      write(*,*) '!** 2 == PAD                               **!'
      write(*,*) '!** 3 == Perfect diamond crystal structure **!'
      write(*,*) '----------------------------------------------'
      read(*,*) GeometryChoice 
      if (GeometryChoice == 1) then
        write(*,*) 'you are about to build single edge dislocation in a cylinder in FCC'
        allocate(hide(ncelltot,atomspercel))
        hide(:,:) = .false.

        NumCellP = 0
        numdel = 0
        do i=1,N(1)
         do j=1,N(2)
          do k=1,N(3)
            NumCellP = NumCellP +1
            do l=1,atomspercel
             if (i.eq.((N(1)-1)/2+1)) then
              if (j.le.N(2)/2) then
                hide(NumCellP,l) = .true.
                numdel = numdel + 1
              endif
             endif
            enddo
          enddo
         enddo
        enddo

        numatomstot = N(1)*N(2)*N(3)*atomspercel-numdel      

        XDislo = UnitCellDim(1)*((N(1)-1)/2+0.25d0) 
        !XDislo = ((N(1)-1)/2)*UnitCellDim(1)
        YDislo = UnitCellDim(2)*(N(2)/2)-0.1
        sr1(1)=XDislo
        sr1(2) = YDislo
        sr1(3)=0.0d0
        factor = burgers/2.0d0/pivalue
        Vburger(1) = burgers;Vburger(2) = 0.0d0;Vburger(3) = 0.0d0
        bx = Vburger(1);by=Vburger(2)
        !*** projection on X
        px(1) = 1.0d0;px(2) = 0.0d0;px(3) = 0.0d0
        !*** projection on Y
        py(1) = 0.0d0;py(2) = 1.0d0;py(3) = 0.0d0
        nua = 1.0d0/4.0d0/(1-0.33d0)/2/pivalue
        nub = (1.0d0-2.0d0*0.33d0)*nua      
      
        !*** apply the displacement field of a dislocation
        NumCellP = 0
        do i=1,N(1)
         do j=1,N(2)
          do k=1,N(3)
           NumCellP = NumCellP + 1
           do l = 1,atomspercel
            if (.not.hide(NumCellP,l)) then
             xx = AtomsTot(NumCellP,l,1) - XDislo
             yy = AtomsTot(NumCellP,l,2) - YDislo

             theta=datan2(-xx,yy)
    
             ux = factor*(theta+xx*yy/2/(1-0.33d0)/(xx*xx+yy*yy))
             uy = -1.0d0*factor*((1.0d0-2.0d0*0.33d0)/4/(1-0.33d0)*log(xx*xx+yy*yy)+(xx*xx-yy*yy)/4.0d0/(1-0.33d0)/(xx*xx+yy*yy))
            !ux = 0.0d0
            !uy = 0.0d0 
             AtomsTot(NumCellP,l,1) = AtomsTot(NumCellP,l,1) + ux
             AtomsTot(NumCellP,l,2) = AtomsTot(NumCellP,l,2) + uy
            endif
           enddo
          enddo
         enddo
        enddo
        write(*,*) 'Dislo position'
        write(*,*) XDislo, YDislo
        write(*,*) 'radius of the cylinder'
        read(*,*) radius
        
        NumCellP = 0
        do i=1,N(1)
         do j=1,N(2)
          do k=1,N(3)
            NumCellP = NumCellP + 1
            do l=1,atomspercel
             if (.not.hide(NumCellP,l)) then
              !** check if atom is inside of cylinder of radius with axis on XDislo, YDislo
              xx = AtomsTot(NumCellP,l,1) - XDislo
              yy = AtomsTot(NumCellP,l,2) - YDislo      
              distance = sqrt(xx*xx+yy*yy)
              if (distance .gt. radius) then
               hide(NumCellP,l) = .true.
               numdel = numdel + 1
              endif
             endif
            enddo
          enddo
         enddo
        enddo
        
        !*** total number of atoms
        numatomstot = N(1)*N(2)*N(3)*atomspercel-numdel
        
        !*** dimension of the box
        dimbox(:,:)=0.0d0
        dimbox(1,2) = N(1)*UnitCellDim(1)
        dimbox(2,2) = N(2)*UnitCellDim(2)
        dimbox(3,2) = N(3)*UnitCellDim(3)

        open(1,file='atoms.diamond.edge.cylinder',status='unknown')
        write(1,*) 'Position data for Ni File'
        write(1,*) 
        write(1,*) numatomstot, 'atoms'
        write(1,*) ' 2 atom types'
        write(1,*) dimbox(1,1:2),' xlo xhi'
        write(1,*) dimbox(2,1:2),' ylo yhi'
        write(1,*) dimbox(3,1:2),' zlo zhi'
        write(1,*) 
        write(1,*) 'Atoms'
        write(1,*)
        m = 1
        iun = 1
        NumCellP=0
        do i=1,N(1)
         do j=1,N(2)
          do k=1,N(3)
            NumCellP = NumCellP+1
            do l=1,atomspercel
             if (.not.hide(NumCellP,l)) then
              write(1,10) m,itypenumcell(NumCellP,l),AtomsTot(NumCellP,l,1:3)
              m = m+1
             endif
            enddo
          enddo
         enddo
        enddo

       close(1)

      deallocate(hide)
        
        
        
        
        
        
        
        
        
        
      elseif (GeometryChoice == 2) then
        write(*,*) 'PAD'
        
         allocate(hide(ncelltot,atomspercel))
        hide(:,:) = .false.
        !** remove half a plane of atoms
        NumCellP = 0
        numdel = 0
        do i=1,N(1)
         do j=1,N(2)
          do k=1,N(3)
            NumCellP = NumCellP +1
            if (i.eq.N(1)) then
             if (j.le.N(2)/2) then
               hide(NumCellP,1:atomspercel) = .true.
               numdel = numdel + atomspercel
             endif
            endif
          enddo
         enddo
        enddo

        numatomstot = N(1)*N(2)*N(3)*atomspercel-numdel
        
        !** scaling by +b/2 for atoms with j <= N(2)/2
!SG april 2011, scaling by +b        !** scaling by b/2 for atoms with j <= N(2)/2
! SG April 2011 deleted        !** scaling by -b/2 for atoms with j > N(2)/2
        numcell = 0
        do i=1,N(1)
         do j=1,N(2)
          do k=1,N(3)
           numcell = numcell + 1
           do l=1,atomspercel
            if (j.le.(N(2)/2)) then
             !*** scale all the x coordinates by +b/2
             xx=AtomsTot(numcell,l,1)
             scaleX = xx*burgers/((N(1)-1)*UnitCellDim(1))
             AtomsTot(numcell,l,1) = AtomsTot(numcell,l,1) + scaleX
            elseif (j .gt. (N(2)/2)) then
             !*** scale all the x coordinates by -b/2
         !    xx=AtomsTot(numcell,l,1)
         !    scaleX = -0.5d0*xx*burgers/((N(1)-1)*UnitCellDim(1))
         !    AtomsTot(numcell,l,1) = AtomsTot(numcell,l,1) + scaleX
            endif
           enddo
          enddo
         enddo
        enddo  
        
        !*** dimension of the simulation cell
        dimbox(:,:) = 0.0d0
        dimbox(1,2) = (N(1))*UnitCellDim(1) !+burgers*0.5d0
        dimbox(2,2) = N(2)*UnitCellDim(2)
        dimbox(3,2) = N(3)*UnitCellDim(3)


        !*** Save the atomic position using LAMMPS format
        open(1,file='atoms.diamond.edge.pad',status='unknown')
        write(1,*) 'Position data for BN File'
        write(1,*) 
        write(1,*) numatomstot, 'atoms'
        write(1,*) ' 2 atom types'
        write(1,*) dimbox(1,1:2),' xlo xhi'
        write(1,*) dimbox(2,1:2),' ylo yhi'
        write(1,*) dimbox(3,1:2),' zlo zhi'
        write(1,*) 
        write(1,*) 'Atoms'
        write(1,*)
        m = 1
        NumCellP=0
        do i=1,N(1)
         do j=1,N(2)
          do k=1,N(3)
            NumCellP = NumCellP+1
            do l=1,atomspercel
             if (.not.hide(NumCellP,l)) then
              write(1,10) m,itypenumcell(NumCellP,l),AtomsTot(NumCellP,l,1:3)
              m = m+1
             endif
            enddo
          enddo
         enddo
        enddo
        
        deallocate(hide)      
        write(*,*) 'for blocked atoms on the bottom and top surface'
        write(*,*) 'ymin == ',UnitCellDim(2)-0.1d0
        write(*,*) 'ymax == ',(N(2)-1)*UnitCellDim(2) - 0.1d0
        
        
        
      elseif (GeometryChoice == 3) then
        write(*,*) 'perfect crystal structure'
        
           write(*,*) 'is it a for B (yes = 1)?'
           read(*,*) B
 
           write(*,*) 'do you want a crack geometry? 1 = yes'
           read(*,*) crack
           if (crack == 1) then
            write(*,*) 'half axis = (z, y)?'
            read(*,*) aa,bb
           ye = N(2)/2.0d0*UnitCellDim(2)
           ze = 0.0d0
           endif

             deleteatoms = 0
           if (crack ==1 ) then
             ncelltot = N(1)*N(2)*N(3)
             deleteatoms = 0
             allocate(hide(ncelltot,atomspercel))
             hide(:,:) = .false.
             NumCellP = 0
             do i=1,N(1)
              do j=1,N(2)
               do k=1,N(3)
                NumCellP = NumCellP+1
                do l=1,12
!                 if (k .le. 32) then
!                  if (j == N(2)/2) then
!                    hide(NumCellP,1:12) = .false.
!                    deleteatoms = deleteatoms + 0
 !*** test sharp crack
!                  !  hide(NumCellP,1:12) = .true.
!                  !  deleteatoms = deleteatoms + 12
!                  endif
!                 endif
                yye = AtomsTot(NumCellP,l,2) - ye
                zze = AtomsTot(NumCellP,l,3) - ze
                 if (yye*yye/bb/bb+zze*zze/aa/aa .lt. 1.0d0) then
                   hide(NumCellP,l) = .true.
                   deleteatoms = deleteatoms + 1
                 endif
                enddo
               enddo
              enddo
             enddo
           write(*,*) N(1)*N(2)*N(3)*atomspercel, deleteatoms
           endif ! end crack == 1

           numatomstot = N(1)*N(2)*N(3)*atomspercel-deleteatoms
           dimbox(:,:) = 0.0d0
           dimbox(1,2) = N(1)*UnitCellDim(1)
           dimbox(2,2) = N(2)*UnitCellDim(2)
           dimbox(3,2) = N(3)*UnitCellDim(3)
           open(1,file='atoms.diamond.perfect.cubic',status='unknown')
           write(1,*) 'Position data for BN File'
           write(1,*) 
           write(1,*) numatomstot, 'atoms'
           write(1,*) ' 2 atom types'
           write(1,*) dimbox(1,1:2),' xlo xhi'
           write(1,*) dimbox(2,1:2),' ylo yhi'
           write(1,*) dimbox(3,1:2),' zlo zhi'
           write(1,*) 
           write(1,*) 'Atoms'
           write(1,*)
           m = 1
           iun = 1
           NumCellP=0
           do i=1,N(1)
            do j=1,N(2)
             do k=1,N(3)
               NumCellP = NumCellP+1
               do l=1,atomspercel
                if (B==1) itypenumcell(NumCellP,l) = 1
                 if (crack == 1) then
                  if (.not. hide(NumCellP,l)) then
                   write(1,10) m,itypenumcell(NumCellP,l),AtomsTot(NumCellP,l,1:3)
                   m = m+1
                  endif
                 else
                   write(1,10) m,itypenumcell(NumCellP,l),AtomsTot(NumCellP,l,1:3)
                   m = m+1
                 endif
               enddo
             enddo
            enddo
           enddo 
           close(1) 
 
        deallocate(itypenumcell)
        deallocate(hide)
        
      else
        write(*,*) 'Geometry choice = ',GeometryChoice,' is not defined'
        deallocate(itypenumcell)
        stop
      endif
      
      
  
   elseif (DislocationCaracter == 2) then
     write(*,*) 'You are about the build a SCREW dislocation in diamond crystal structure'
     !*** here we need to define the type of geometry we want to use
     write(*,*) '-------------------------------------------'
     write(*,*) '!** What kind of structure do you want? **!'
     write(*,*) '!** 1 == cylinder                       **!'
     write(*,*) '!** 2 == PAD                            **!'
     write(*,*) '-------------------------------------------'
     read(*,*) GeometryChoice
     if (GeometryChoice == 1) then
      !*** define the position of the dislocation line
      midle(:) = 0.0d0
      midle2(:) = 0.0d0
      epsilon = 0.0001d0
      midle(1) = N(1)/2*UnitCellDim(1) !- epsilon
      midle(3) = N(3)/2*UnitCellDim(3) - 0.0d0   !** Z-
      midle2(3) = N(3)/2*UnitCellDim(3) + 0.0d0  !** Z+
      midle(2) = N(2)/2*UnitCellDim(2) + 0.5d0          !** Y


      write(*,*) 'which core  : C1 = 1; C2 = 2'
      read(*,*) CoreChoice
      
      if (CoreChoice == 1) then
      !** core C1
      x1=AtomsCellCoor(2,2)
      x2=AtomsCellCoor(9,2)
      midle(2)=N(2)/2*UnitCellDim(2)+(x1+x2)*0.5d0
      x1=AtomsCellCoor(5,3)
      x2=AtomsCellCoor(4,3)
      midle(3) = N(3)/2*UnitCellDim(3)+0.5d0*(x1+x2)
      midle2(3) = N(3)/2*UnitCellDim(3)+0.5d0*(x1+x2)
  
 !     center(2) = N(2)/2.0d0*UnitCellDim(2) !** Y
 !     center(1) = N(3)/2.0d0*UnitCellDim(3) !** Z

      center(2) = midle(2)
      center(1) = midle(3)
      else if (CoreChoice == 2) then
      ! core C2
      x1=AtomsCellCoor(8,2)
      x2=AtomsCellCoor(9,2)
      midle(2)=(N(2)/2-1)*UnitCellDim(2)+(x1+x2)*0.5d0
      x1=AtomsCellCoor(8,3)
      x2=AtomsCellCoor(9,3)
      midle(3) = (N(3)/2-1)*UnitCellDim(3)+0.5d0*(x1+x2)
      midle2(3) = (N(3)/2-1)*UnitCellDim(3)+0.5d0*(x1+x2)
  
 !     center(2) = N(2)/2.0d0*UnitCellDim(2) !** Y
 !     center(1) = N(3)/2.0d0*UnitCellDim(3) !** Z

      center(2) = midle(2)
      center(1) = midle(3)     
      
      
      
      
      endif

      write(*,*) 'dislocation line (Y, Z)'
      !write(*,*) midle(1),center(2),center(1)
      write(*,*) center(2), center(1)

      allocate(hide(ncelltot,atomspercel))
      hide(:,:) = .false.

      !*** define the radius of the cylinder
      write(*,*) 'What is the cylinder radius?'
      read(*,*) radius
      
      write(*,*) 'Isotropic solution? (1 == yes)'
      read(*,*) isotropic
      

      !*** apply the dislocation displacement from linear elasticity
      !*** cut a cylinder around the dislocation line
      NumCellP = 0
      numhide=0
      if (isotropic == 1) then
      !* isotropic solution
      do i=1,N(1)
       do j=1,N(2)
        do k=1,N(3)
         NumCellP = NumCellP + 1
         do l=1,atomspercel

          !* tan**(1)(y/z)
           xx = AtomsTot(NumCellP,l,2)
           yy = AtomsTot(NumCellP,l,3)

           !*** screw dislocation displacement
!           AtomsTot(NumCellP,l,1) = AtomsTot(NumCellP,l,1) + burgers/2.0d0/pivalue*(datan2((xx-midle(2)),(yy-midle(3))))
           AtomsTot(NumCellP,l,1) = AtomsTot(NumCellP,l,1) + burgers/2.0d0/pivalue*(datan2((yy-midle(3)),(xx-midle(2))))
!           AtomsTot(NumCellP,l,1) = AtomsTot(NumCellP,l,1) + burgers/4.0d0/pivalue*(datan2((xx-midle(2)),(yy-midle2(3))))

           !*** check if atoms inside the cylinder
           xx = xx - center(2)
           yy = yy - center(1)
           distance = sqrt(xx*xx+yy*yy)
           if (distance .gt. radius) then
             numhide = numhide + 1
             hide(NumCellP,l) = .true.
           endif

         enddo
        enddo
       enddo
      enddo
      
      else if (isotropic == 2) then
      !* anisotropic solution
      c11 = 911d0
      c12 = 166d0
      c44 = 465d0

      H = 2.0d0*c44+c12-c11

      cp44 = c44 - H/3.0d0
      cp45 = -1.0d0*sqrt(2.0d0)*H/6.0d0
      cp55 = c44 - H/6.0d0


      
      
      do i=1,N(1)
       do j=1,N(2)
        do k=1,N(3)
         NumCellP = NumCellP + 1
         do l=1,atomspercel

          !* tan**(1)(y/z)
           xx = AtomsTot(NumCellP,l,2) - midle(3)
           yy = AtomsTot(NumCellP,l,3) - midle(2)

           !*** screw dislocation displacement
!           AtomsTot(NumCellP,l,1) = AtomsTot(NumCellP,l,1) + burgers/2.0d0/pivalue*(datan2((xx-midle(2)),(yy-midle(3))))
!           AtomsTot(NumCellP,l,1) = AtomsTot(NumCellP,l,1) + burgers/2.0d0/pivalue*(datan2((yy-midle(3)),(xx-midle(2))))
            Angle = datan2((sqrt(cp44*cp55-cp45*cp45)*xx),(cp44*yy+cp45*xx))
            AtomsTot(NumCellP,l,1) = AtomsTot(NumCellP,l,1) + burgers/2.0d0/pivalue*Angle
!           AtomsTot(NumCellP,l,1) = AtomsTot(NumCellP,l,1) + burgers/4.0d0/pivalue*(datan2((xx-midle(2)),(yy-midle2(3))))

           !*** check if atoms inside the cylinder
           xx = (xx + midle(3))- center(2)
           yy = (yy + midle(2)) - center(1)
           distance = sqrt(xx*xx+yy*yy)
           if (distance .gt. radius) then
             numhide = numhide + 1
             hide(NumCellP,l) = .true.
           endif

         enddo
        enddo
       enddo
      enddo      
      endif

      !*** total number of atoms
      numatomstot = N(1)*N(2)*N(3)*atomspercel-numhide

      dimbox(:,:)=0.0d0
      dimbox(1,2) = N(1)*UnitCellDim(1)
      dimbox(2,2) = N(2)*UnitCellDim(2)
      dimbox(3,2) = N(3)*UnitCellDim(3)

      open(1,file='atoms.diamond.screw.cylinder',status='unknown')

      write(1,*) 'Position data for BN File'
      write(1,*)
      write(1,*) numatomstot, ' atoms'
      write(1,*) ' 2 atom types'
      write(1,*) dimbox(1,1:2),' xlo xhi'
      write(1,*) dimbox(2,1:2),' ylo yhi'
      write(1,*) dimbox(3,1:2),' zlo zhi'
      write(1,*)
      write(1,*) 'Atoms'
      write(1,*)
      m=1
      iun=1
      NumCellP=0
      do i=1,N(1)
       do j=1,N(2)
        do k=1,N(3)
         NumCellP = NumCellP+1
         do l=1,atomspercel
          if (.not.hide(NumCellP,l)) then
            write(1,10) m, itypenumcell(NumCellP,l), AtomsTot(NumCellP,l,1:3)
            m  = m+1
          endif
         enddo
        enddo
       enddo
      enddo

     close(1)


     else if (GeometryChoice == 2) then
      write(*,*) 'you are about to build SCREW dislocation in diamond structure with PAD'
      write(*,*) 'not implemented yet !!!!'
     else
      write(*,*) 'unknown Geometry'
     endif


     deallocate(itypenumcell)
     
  endif

end select


deallocate(itype)
deallocate(AtomsCellCoor)
deallocate(AtomsTot)

10 format(i8, i8, 3(1x, e12.5))

end program

