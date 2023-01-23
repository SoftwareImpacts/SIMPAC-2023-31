
   
 program StochHyProp
 
 ! Quirijn de Jong van Lier - 2021, 2022, 2023  - qdjvlier@usp.br
 ! Perform CHolesky decomposition and generate stochastic parameter realizations
 ! Input file: [StochHyProp].stc
 ! Output files: [StochHyProp].rep, [StochHyProp]_###.out, [StochHyProp]_###.csv
   
   implicit none


   character*24 commandname, BaseName, ProgName
   character*64 InFile, OutFile, ExtFile, CsvFile, MrgFile, CsvOpen, OutOpen
   character*10 DumChar, Version
   character*80 Tmp1, Tmp2, Form
   character*160 Tmp3, Tmp4
   character*3, st
   character*4, InpFileExt

   Logical Fail, Key, ConsStDev, Extern, Merged, CreateMerged
   
   integer NParm, NRealizations, NSamples
   integer i, j, k, ip, n, mode, start, finish, s, asc
   integer*2 NumChar

   real*8 CoefInvNorm(0:63), pi
   real*4 Tail, aa, phi, StDev
   real*16 ch
   real rate
   
   real*4, dimension(:,:), allocatable :: Corr, Diag, Rnd, Corr2, xxmyym
   real*4, dimension(:,:,:), allocatable :: Realiz
   real*16, dimension(:,:), allocatable :: L, CoVar
   real*4, dimension(:), allocatable :: Mean, StErr, T, Rand, StErr2, xxm, xxm2, yym, yym2
   real*8, dimension(:), allocatable :: Mean2
   character*10, dimension(:), allocatable :: ParName

   !call CPU_TIME(start)
   call System_Clock(start, rate, mode)
   
   ProgName = 'StochHyProp'
   InpFileExt = '.stc'
   Version = '2.01'
   write (*,*) trim(ProgName) // ' v. ', trim(Version)

   Fail = .false.
   Extern = .false.
   Merged = .false.
   
!Calculate coefficients for inverse normal distribution equation (erf-1)
   pi = 4*atan(1.)
   CoefInvNorm = 0
   CoefInvNorm(0) = 1
   do i = 1,63
      do j = 0,i-1
         CoefInvNorm(i) = CoefInvNorm(i) + CoefInvNorm(j)*CoefInvNorm(i-1-j)/(j+1)/(2*j+1) 
      enddo 
   enddo    
   do i = 0,63
      CoefInvNorm(i) = (CoefInvNorm(i) / (2.*i+1)) * pi**i * sqrt(pi) / (2.**(2.*i+1))
   enddo     
   
   
!   Path and filename of executable through argument command line
    Call GetArg (1,BaseName,NumChar)
    do i=1,len(BaseName)
        asc = iachar(BaseName(i:i))
        if ((asc.ge.65).and.(asc.le.90)) then
          BaseName(i:i) = char(asc+32)  ! transform to lower case
        endif 
    enddo

    if(NumChar.gt.3) then
        if (BaseName(NumChar-3:NumChar).eq.InpFileExt) then
           InFile = trim(BaseName)
        else 
           InFile = trim(BaseName)//InpFileExt
        endif
    else
        InFile = trim(ProgName) //InpFileExt
    end if
    

    
! Read input file    
   OPEN (unit=3, FILE = InFile , STATUS='old')  ! Main information file

   write (*,*) 'Reading input file ' // InFile
   do i=1,3   
     read (3,*) dumchar
   enddo
   read (3,*) NSamples
   read (3,*) n
   if (n.eq.1) then
       CreateMerged = .true.
   else
       CreateMerged = .false.
   endif 
  
   read (3,*) NRealizations
   read (3,*) Tail
   read (3,*) n
   if (n.eq.1) then
       ConsStDev = .true.
   else
       ConsStDev = .false.
   endif 
   read (3,*) n
   if (n.eq.1) then
       Key = .true.
   else
       Key = .false.
   endif 
   if (Tail.gt..49) then
       write (*,*) 'Tail should not exceed 0.49.'
       goto 999
   endif 

   read (3,*) OutFile
   if (trim(OutFile).eq.'*') then
    i = index(InFile,".")
	OutFile = InFile(1:i) // 'out'
   endif
   i = index(OutFile,".")
   CsvFile = OutFile(1:i) // 'csv'
   MrgFile = OutFile(1:i) // 'rep'

   OPEN (unit=55, FILE = MrgFile , STATUS='replace')  ! MrgOutput file
      write (55, *) '*  ' // trim(ProgName) // ' v. ', trim(Version)      

   read (3,*) n
   if (n.eq.1) then
       if (NSamples.gt.1) then
          write (*,*) 'Use NSamples = 1 with RETC or Hydrus output file.'
          goto 999
       endif   

       Extern = .true.
   endif
   
   read (3,*) ExtFile
   i = index(ExtFile,".")
   if (i.eq.0) then
       ExtFile = trim(ExtFile)//'.out'
   endif
   !write(*,*) Hydrus, ExtFile
   !read*
   
   
   
 if (Extern) then
     
   OPEN (unit=6, FILE = ExtFile , STATUS='old')  ! Main information file
   read (6,'(A80)') Tmp1
   n=len(trim(Tmp1))  !n=73 for Hydrus, 72 for RETC
   
   if (n.eq.73) then
      write (*,*) 'Reading Hydrus output file ' // ExtFile
   else   
      write (*,*) 'Reading RETC output file ' // ExtFile
   endif
   
   do while (dumchar.ne.'Correlatio')  
     read (6,*) dumchar
   enddo
   read (6,*) dumchar
   if (n.ne.73) then
       read (6,*) dumchar
   endif    
   read (6,'(A80)') Tmp1
   do i=1,80
      Tmp2(81-i:81-i)=Tmp1(i:i)
   enddo
   read (Tmp2,*) NParm
   
  
 else    
   
   read (3,*) dumchar
   read (3,*) dumchar
   read (3,*) dumchar
   read (3,*) dumchar
   read (3,'(A160)') Tmp3
   do i=1,160
      Tmp4(161-i:161-i)=Tmp3(i:i)
   enddo
   read (Tmp4,*) NParm

 endif  


 
   Allocate (Corr(NParm,NParm),L(NParm,NParm))
   Allocate (CoVar(NParm,NParm), Diag(NParm,NParm), Corr2(NParm,NParm), xxmyym(NParm,NParm))
   Allocate (Mean(NParm), Mean2(NParm), StErr(NParm), StErr2(NParm), ParName(NParm))
   Allocate (xxm(NParm), xxm2(NParm), yym(NParm), yym2(NParm))
   Allocate (Rnd(NParm,NRealizations), T(NParm), Realiz(NParm,NRealizations,NSamples), Rand(NParm))
   Corr=0

! ************************************************************************************************
 
   write (*,'(A11, I8, A18, I2, A15)') 'Processing ', NRealizations, ' realizations for ', NParm, ' parameters.'
   if (NSamples.gt.1) then
      write (*,*) NSamples, 'samples'
   endif    
   
   write (*,*) 

444 Continue

 do s = 1, NSamples  
     mean2 = 0
     xxm = 0
     xxm2 = 0
     xxmyym = 0
     Corr2 = 1

   if ((NSamples.gt.1).or.(Merged)) then
      do i = 6, 24
         call CursPos(i,0)
         write (*,*) '                                                                                                                                                                                                           '
      enddo    
      call CursPos(6,1)
      if (.not.Merged) then
            write (*,'(a,I3.3,a,I3.3)') 'Processing Sample ', s, '/', NSamples
      else      
            write (*,*) 'Processing Merged Sample Set                                        '               
      endif
      
      i = index(CsvFile,".")
      CsvOpen = ''
      OutOpen = ''
      st = ''
      if (Merged) then
          st = 'mrg'
      else
          write (st,'(I3.3)') s
      endif
      write (CsvOpen, *) trim(CsvFile(1:i-1)), '_', trim(st), '.csv' 
      write (OutOpen, *) trim(OutFile(1:i-1)), '_', trim(st), '.out' 
   else 
      CsvOpen = CsvFile
      OutOpen = OutFile
   endif
   
      OPEN (unit=4, FILE = OutOpen , STATUS='replace')  ! Output file
      OPEN (unit=5, FILE = CsvOpen , STATUS='replace')  ! CsvOutput file
      
      
      write (4,*) '*  ' // trim(ProgName) // ' v. ', trim(Version)
      write (5,*) trim(ProgName) // ' v. ', trim(Version)


   
 if (Extern) then  
    
   do i = 1,NParm
      read (6,*) dumchar, (Corr(i,j), j = 1,i)
   enddo
   do while (trim(dumchar).ne.'Variable')  
     read (6,*) dumchar
   enddo

   do i = 1,NParm
      read (6,*) ParName(i), Mean(i), StErr(i)
   enddo
   
   close(6)
   

 elseif (.not.Merged) then
        
    read (3,*) dumchar
    
    do i = 1,NParm
      read (3,*) ParName(i), Mean(i), StErr(i), (Corr(i,j), j = 1,i)
      if (StErr(i).eq.0) then
          StErr(i) = Mean(i)/1.e8
      endif    
      if (StErr(i).eq.0) then
          StErr(i) = 1.e-8
      endif    
      if (Corr(i,i).ne.1.) then
         write (*,*) 'Autocorrelation of parameter ' // trim(ParName(i)) // ' should be equal to 1.'
         goto 999
      endif    
    enddo
   
endif  
 
!   write (*,*) ((trim(ParName(i)), ', '), i=1,NParm)
   
   Diag=0
   CoVar=0
   L=0
   T=0
   
   call random_seed()
   call random_number(Rnd)
   Rnd = Tail + Rnd*(1-2*Tail)


   

   do i = 1,NParm
      do j = 1,i-1
         Corr(j,i) = Corr(i,j)
      enddo
      Diag(i,i) = StErr(i)
   enddo

   
  
! Choleski factorisation (https://moonbooks.org/Codes/Cholesky-decomposition-in-fortran-90/)

   call CursPos(7,1)
   write (*,*) '. Performing Cholesky factorization'


   CoVar = (matmul (matmul (Diag, Corr), Diag))
   L = 0.

   L(1,1) = sqrt( CoVar(1,1) )
   do i = 2 , NParm
	 L(i,1) = CoVar(1,i) / L(1,1)
   enddo

   do j = 2 , NParm
	  do i = 1 , NParm
		if (i.lt.j) then
		   L(i,j) = 0.0
		else
		   if (i.eq.j) then
			ch = 0.
			do k = 1 , i -1
			   ch = ch + L(i,k)**2
            enddo
            if (CoVar(i,i).lt.ch) then
              Fail = .true.
            endif
			L(i,j) = sqrt(CoVar(i,i) - ch)
		   else
			ch = 0.
			do k = 1 , i - 1
			   ch = ch + L(j,k)*L(i,k)
			enddo
			L(i,j) = (CoVar(j,i) - ch) / L(j,j)
		   endif
		endif
	  enddo
   enddo   
   
   
     write (4,'(A1,20A16)') '*', ParName
     write (5,*) ((trim(ParName(i)),','), i = 1,NParm)

   if ((ConsStDev).and.(Tail.gt.0)) then
      call InvNorm(aa, 0., 1., Tail, CoefInvNorm)
      phi = exp(-0.5*(aa**2)) / sqrt(2*pi)
      StDev = 1. / sqrt(1 + 2*aa*phi/(1.-2.*Tail))
   else
      StDev = 1.
   endif
   
   call CursPos(7,1)
   write (*,*) 'X Performing Cholesky factorization'
   call CursPos(8,1)
   write (*,*) '. Generating inverse normal randoms for realization '

   do i=1,NRealizations
     if (mod(i,1000).eq.0) then
        call CursPos(8,53)
        write (*,*) i, '/', NRealizations
     endif    
     do j=1,NParm  
       call InvNorm(Rand(j), 0., StDev, Rnd(j,i), CoefInvNorm)
     enddo 
     T = matmul(L, Rand)
      do j = 1,NParm
         Realiz(j,i,s) = T(j) + Mean(j)
      enddo    
      
      
      write(4,'(20F16.8)') ((Realiz(j,i,s)), j = 1,NParm)
      write (5,'(*(F16.8,","))') ((Realiz(j,i,s)), j = 1,NParm)
      
      
   enddo
   call CursPos(8,1)
   write (*,*) 'X Generating inverse normal randoms                                                                 '

   write (4,*) '* End of file'

   write(*,*)

   if (Fail) then
     write (*,*) ''
     write (*,*) 'Sorry, Cholesky decomposition FAILED.'
     write (*,*) 'Check consistency of input parameters.'
   else
   call CursPos(9,1)
   write (*,*) '. Generating realizations ...'

     mean2 = 0
     xxm = 0
     xxm2 = 0
     xxmyym = 0
     
     do j=1,NParm  
        do i=1,NRealizations
           mean2(j) = mean2(j)+ Realiz(j,i,s)
 
        enddo
        mean2(j) = mean2(j) / NRealizations
     enddo
     
     do j=1,NParm  
        do i=1,NRealizations
           xxm(j) = xxm(j) + Realiz(j,i,s)-mean2(j)
           xxm2(j) = xxm2(j) + (Realiz(j,i,s)-mean2(j))**2
        enddo
        StErr2(j) = sqrt(xxm2(j) / NRealizations)
     enddo
     
     do i=1,NParm  
       do j=i+1,NParm
          do k=1,NRealizations
            xxmyym(i,j) = xxmyym(i,j) + (Realiz(i,k,s)-mean2(i))*(Realiz(j,k,s)-mean2(j))  
          enddo  
       enddo
     enddo  

     do i=1,NParm  
       do j=i+1,NParm
         Corr2(i,j) = xxmyym(i,j) / sqrt(xxm2(i)*xxm2(j))   
       enddo
     enddo  
   call CursPos(9,1)
   write (*,*) 'X Generating realizations     '

   endif
    
   Form = '(A12, 2F12.5, F12.4, 20F8.4)'
     if (.not.merged) then
      
      write (55,*) 'Input values sample', s
      do i=1,NParm  
       
       write (55, Form) ParName(i), Mean(i), StErr(i), (Corr(j,i), j = 1,i)
      enddo
     endif
     
        
       write (55,*) 
       
       write (55,*) 'Properties of the set of', NRealizations, ' generated realizations'
       do i=1,NParm  
         
         write (55, Form) ParName(i), Mean2(i), StErr2(i), (Corr2(j,i), j = 1,i)
       enddo
       
       write (*,*) 
       write (*,*) 
       write (55,*) 
       write (55,*) 
       
     close(4)
     close(5)

     
 enddo ! s - NSamples    - Calculate composed sample paramters

 if ((NSamples.gt.1).and.(CreateMerged)) then
     
     mean2 = 0
     xxm = 0
     xxm2 = 0
     xxmyym = 0
     Corr2 = 1

     do j=1,NParm  
       do s=1,NSamples
         do i=1,NRealizations
            mean2(j) = mean2(j)+ Realiz(j,i,s)
         enddo
        enddo  
        mean2(j) = mean2(j) / NRealizations / NSamples
     enddo
     
     do j=1,NParm  
       do s=1,NSamples  
        do i=1,NRealizations
           xxm(j) = xxm(j) + Realiz(j,i,s)-mean2(j)
           xxm2(j) = xxm2(j) + (Realiz(j,i,s)-mean2(j))**2
        enddo
       enddo 
       StErr2(j) = sqrt(xxm2(j) / NRealizations / NSamples)
     enddo
     
     do i=1,NParm  
       do j=i+1,NParm
         do s=1,NSamples  
          do k=1,NRealizations
            xxmyym(i,j) = xxmyym(i,j) + (Realiz(i,k,s)-mean2(i))*(Realiz(j,k,s)-mean2(j))  
          enddo  
         enddo 
       enddo
     enddo  

     do i=1,NParm  
       do j=i+1,NParm
         Corr2(i,j) = xxmyym(i,j) / sqrt(xxm2(i)*xxm2(j))   
       enddo
     enddo  

 
     Form = '(A12, 2F12.5, F12.4, 20F8.4)'
     
     do i=1,3
       write (*,*) 
       write (55,*)
     enddo  
     
     write (55,*) 'MERGED SAMPLE PROPERTIES: '
       do i=1,NParm  
         
         write (55, Form) ParName(i), Mean2(i), StErr2(i), (Corr2(j,i), j = 1,i)
       enddo
      
     if (CreateMerged) then
         Merged = .true.  
     else    
         Merged = .false.  
     endif
 else
     Merged = .false.  
     
 endif
 
     If (Merged) then 
         
        NSamples = 1
        Extern = .false.

        Mean = Mean2
        StErr = StErr2
        Corr = Corr2
        do i = 1,NParm
          do j = 1,i-1
            Corr(i,j) = Corr(j,i)
          enddo
        enddo

       goto 444  
     endif
         
         
     close(55)
      !Call CPU_TIME(finish)
     call System_Clock(finish, rate, mode)
     call CursPos(6,1)
     write (*,*) 'X All samples processed                                                                                                                                        '
     do i=8,9
         call CursPos(i,1)
         write (*,*) '                                                                                                                                           '
     enddo
     call CursPos(8,1)

          write (*,*) 
          write (*,*) 'Processing completed.'
          write (*,*) 'Output files created.'

          write (*,*) 
          write (*,'(F7.2, A30)') (finish-start)/rate, 'seconds of processing time.'

   
999 continue
    
   if (Key) then
     write (*,*) 'Any key to exit'
     read*
   endif
   
    end program

    
    
! Subroutine InvNorm returns value x for which cumulative density function equals cdf,
! given mean and stdev
! Quirijn and Matias de Jong van Lier, 2017.
    subroutine InvNorm(x,mean,stdev,cdf,CoefInvNorm)
    implicit none
    real mean, stdev, x, cdf, c
    real*8 CoefInvNorm(0:63)
    integer i, j
    
    c = 2 * cdf - 1
    
    x=0
    do i=0,63
       x = x + CoefInvNorm(i)* (c**(2*i+1)) 
    enddo    
    x = x * stdev * sqrt(2.) + mean
    
    return
    end
    
    subroutine CursPos(irow,icol)
    use dflib
    use dfwin
    !implicit none
    integer fhandle
    logical lstat
    Type(T_COORD) wpos
    fhandle = GetStdHandle(STD_OUTPUT_HANDLE)
    if (icol.ge.0) then
       wpos.x = icol
    endif
    if (irow.ge.0) then
       wpos.y = irow
    endif
    lstat = SetConsoleCursorPosition(fhandle, wpos)

    return
    end
