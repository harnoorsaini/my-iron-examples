
!!!!!!!!!!!!!!!!!_______                     HEADER COMMENTS.             ___________!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! THE MODULE ENCAPSULATED THE FOLLOWING INFORMATION: ..... !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! _________________________________________________________!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE parsing



  USE OpenCMISS
  USE OpenCMISS_Iron

  INTEGER                                      :: i,FileStat,n,Stat,POS1,POS2                                             !! These variables are used as indices in the parsing algorithm...

  CHARACTER(LEN=100)                           :: Rdline                                                                        !! This variable stores the lines READ  Statement.
  CHARACTER(LEN=100),ALLOCATABLE               :: EquationSetKeyWords(:),ProblemKeyWords(:),MaterialFieldKeyWords(:)            !! Data Structures holding keywords that the parsing algorithm...
  CHARACTER(LEN=100),ALLOCATABLE               :: FibreFieldKeyWords(:),SolverKeywords(:),CoordinateSystemKeywords(:)
  CHARACTER(LEN=100),ALLOCATABLE               :: MeshKeyWords(:),DecompositionKeyWords(:),DependentFieldKeyWords(:)
  CHARACTER(LEN=100),ALLOCATABLE               :: BasisKeyWords(:),OutputKeyWords(:),BoundaryConditionKeywords(:)
  CHARACTER(LEN=100),ALLOCATABLE               :: RegionKeywords(:),ControlLoopKeyWords(:),FieldKeyWords(:),FunctionKeyWords(:)
  CHARACTER(LEN=100),ALLOCATABLE               :: GeometricFieldKeyWords(:),GeneratedMeshKeyWords(:),FieldsKeywords(:)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! DEFINING DERIVED TYPES TO STORE USER INPUTS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  CHARACTER(LEN=100),ALLOCATABLE               :: BoundaryConditionId(:,:),BoundaryDirichelet(:,:), &
                                                    & BoundaryTractionNeumann(:,:), BoundaryPressureNeumann(:,:)
  CONTAINS


  SUBROUTINE   GenericParsing(KeywordArray, BlockToBeParsed,NumOfArguments,Id,Rdline, Argument1,Argument2,Argument3, &
    & Argument4,Argument5,Argument6,Argument7,Argument8,Argument9,Argument10,Argument11,Argument12,Argument13, &
      & Argument14, Argument15)

    !! These variables are used as indices in the parsing algorithm.
    INTEGER                                                    :: counter,FileStat,n,Stat,k,POS1,POS2,NumOfArguments,Id
    !! Arguments where user inputs  are stored.
    CHARACTER(LEN=*)    , OPTIONAL   , DIMENSION(:)            :: Argument1,Argument2,Argument3,Argument4,Argument5
    CHARACTER(LEN=*)    , OPTIONAL   , DIMENSION(:)            :: Argument6,Argument7,Argument8,Argument9,Argument10
    CHARACTER(LEN=*)    , OPTIONAL   , DIMENSION(:)            :: Argument11,Argument12,Argument13,Argument14,Argument15
    ! Variable that temporarily store the lines of input file
    CHARACTER(LEN=*)                                           :: Rdline
    ! THis is the data structure that contains the keywords
    CHARACTER(LEN=*)    , DIMENSION(:)                         :: KeywordArray
    ! Block to be parsed in the input file.
    CHARACTER(LEN=*)                                           :: BlockToBeParsed
    ! Data structure where the user inputs are temporarily stored
    CHARACTER(LEN=100)  , ALLOCATABLE                          :: StoreArguments(:,:)
    LOGICAL                                                    :: Flag
    ! Allocating the size of StoreArgument data structure that temprarily stores the input parameters

    ALLOCATE(StoreArguments(15,15))

    !!!! STORING DEFAULT VALUE

       IF  (present(Argument1)) THEN
          StoreArguments(1,:) = Argument1
       END IF

       IF  (present(Argument2) ) THEN
         StoreArguments(2,:) = Argument2
       END IF

       IF  (present(Argument3) ) THEN
         StoreArguments(3,:) = Argument3
       END IF

       IF  (present(Argument4) ) THEN
         StoreArguments(4,:) = Argument4
       END IF
       IF  (present(Argument5)) THEN
         StoreArguments(5,:) = Argument5
       END IF

       IF  (present(Argument6)) THEN
         StoreArguments(6,:) =  Argument6
       END IF

       IF  (present(Argument7)) THEN
         StoreArguments(7,:) = Argument7
       END IF

       IF  (present(Argument8) ) THEN
         StoreArguments(8,:) = Argument8
       END IF

       IF  (present(Argument9)) THEN
         StoreArguments(9,:) =  Argument9
       END IF

       IF  (present(Argument10) ) THEN
         StoreArguments(10,:) = Argument10
       END IF
       IF  (present(Argument11)) THEN
         StoreArguments(11,:) = Argument11
       END IF

       IF  (present(Argument12)) THEN
         StoreArguments(12,:) = Argument12
       END IF

       IF  (present(Argument13)) THEN
         StoreArguments(13,:)  = Argument13
       END IF

       IF  (present(Argument14)) THEN
         StoreArguments(14,:) = Argument14
       END IF

       IF  (present(Argument15)) THEN
         StoreArguments(15,:) =  Argument15
       END IF

    Rdline = to_upper(Rdline)
                                        !! converting the string to upper case
    IF (TRIM(Rdline)== "START_" // TRIM(BlockToBeParsed)) THEN        !! Start of parsing the  block in the input file.

       Id = Id + 1

       Readline:DO

          counter=0

          Flag = .TRUE.

          READ(12,'(A)',IOSTAT=FileStat), Rdline
          Rdline = to_upper(Rdline)

          Rdline =  strip_space(Rdline)

          CALL remove_CHARACTER(Rdline,"!")

          IF (TRIM(Rdline)=="END_"//TRIM(BlockToBeParsed)) THEN

            EXIT

          END IF

          DO WHILE (Flag .EQV. .TRUE. .AND. counter .LE. NumOfArguments)

            counter = counter + 1

            IF (counter ==  NumOfArguments) FLAG = .FALSE.  ! to stop memory leakage
            IF (Rdline.EQ.KeywordArray(counter)) THEN

              Flag = .FALSE.

              READ(12,'(A)',IOSTAT=FileStat), Rdline

              Rdline = to_upper(Rdline)

              CALL skip_blank_lines(Rdline)

              CALL remove_CHARACTER(Rdline,"!")

              Rdline =  strip_space(Rdline)

              CALL split_inline_argument_new(Rdline,StoreArguments(counter,:))



            END IF

          END DO

       !! storing new values
      END DO Readline

        IF  (present(Argument1)) THEN
           Argument1 = StoreArguments(1,:)
       END IF

       IF  (present(Argument2) ) THEN
          Argument2 = StoreArguments(2,:)
       END IF

       IF  (present(Argument3) ) THEN
          Argument3 = StoreArguments(3,:)
       END IF

       IF  (present(Argument4) ) THEN
          Argument4 = StoreArguments(4,:)
       END IF

       IF  (present(Argument5)) THEN
         Argument5 = StoreArguments(5,:)
       END IF

       IF  (present(Argument6)) THEN
         Argument6 = StoreArguments(6,:)
       END IF

       IF  (present(Argument7)) THEN
         Argument7 = StoreArguments(7,:)
       END IF

       IF  (present(Argument8) ) THEN
         Argument8 = StoreArguments(8,:)
       END IF

       IF  (present(Argument9)) THEN
         Argument9 = StoreArguments(9,:)
       END IF

       IF  (present(Argument10) ) THEN
         Argument10 = StoreArguments(10,:)
       END IF

       IF  (present(Argument11)) THEN
         Argument11 = StoreArguments(11,:)
       END IF

       IF  (present(Argument12)) THEN
         Argument12 = StoreArguments(12,:)
       END IF

       IF  (present(Argument13)) THEN
         Argument13 = StoreArguments(13,:)
       END IF

       IF  (present(Argument14)) THEN
         Argument14 = StoreArguments(14,:)
       END IF

       IF  (present(Argument15)) THEN
         Argument15 = StoreArguments(15,:)
       END IF



    END IF


  END SUBROUTINE GenericParsing


  !!!!!!!!!!!!!!!!!!!!!!! FOLLOWING SUBROUTINE PARSE "BOUDARY CONDITIONS"  IN INPUT FILE. !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!! FOR COMMENTS PLEASE REFER TO THE "Field_parsing_subroutine."  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE   BC_parsing_subroutine(String1,boundary_conditions_arg1,boundary_conditions_arg2,boundary_conditions_arg3, &
                 & boundary_conditions_arg4, Rdline,num_of_dirichelet)


    INTEGER                                            :: i,FileStat,n,Stat,k,POS1,POS2,num_of_dirichelet
    CHARACTER(LEN=100) , dimension(:,:),ALLOCATABLE    :: boundary_conditions_arg1,boundary_conditions_arg2, &
                                                            & boundary_conditions_arg4,boundary_conditions_arg3
    CHARACTER(LEN=100)                                 :: Rdline,g
    CHARACTER(LEN=100),   dimension(:)                 :: String1


    Rdline = to_upper(Rdline)
    IF (trim(Rdline)=="START_BOUNDARY_CONDITIONS") then
      num_of_dirichelet= 0
      num_of_traction_neumann=0
      num_of_pressure_neumann=0

      Readline: DO

        READ(12,'(A)'), Rdline
        Rdline = to_upper(Rdline)
        Rdline =  strip_space(Rdline)
        CALL remove_CHARACTER(Rdline,"!")
        IF (trim(Rdline)=="END_BOUNDARY_CONDITIONS") EXIT

          ! parse ID
          IF (Rdline.EQ.String1(1)) then
             READ(12,'(A)',IOSTAT=FileStat), Rdline
             Rdline = to_upper(Rdline)
             CALL skip_blank_lines(Rdline)
             CALL remove_CHARACTER(Rdline,"!")
             Rdline =  strip_space(Rdline)
             boundary_conditions_arg1(:,1) = Rdline
          END IF
          ! parse Dirichlet BC
          IF (Rdline.EQ.String1(2)) then
            ! get the complete line
            READ(12,'(A)',IOSTAT=FileStat), Rdline
            ! skip it if the line is empty
            CALL skip_blank_lines(Rdline)
            Rdline = to_upper(Rdline)
            Rdline =  strip_space(Rdline)
            CALL remove_CHARACTER(Rdline,"!")
            ! get: [ TYPE , LOCATION , COMPONENTS , VALUE, ADD|SET]
            CALL split_inline_argument(Rdline,boundary_conditions_arg2,num_of_dirichelet)
          ENDIF

          IF (Rdline.EQ.String1(3)) then
            READ(12,'(A)',IOSTAT=FileStat), Rdline

            Rdline = to_upper(Rdline)
            CALL skip_blank_lines(Rdline)
            CALL remove_CHARACTER(Rdline,"!")
            Rdline =  strip_space(Rdline)

            CALL split_inline_argument(Rdline,boundary_conditions_arg3,num_of_traction_neumann)
          END IF

          IF (Rdline.EQ.String1(4)) then
            READ(12,'(A)',IOSTAT=FileStat), Rdline
            Rdline = to_upper(Rdline)
            CALL skip_blank_lines(Rdline)
            CALL remove_CHARACTER(Rdline,"!")
            Rdline =  strip_space(Rdline)
            CALL split_inline_argument(Rdline,boundary_conditions_arg4,num_of_pressure_neumann)
          END IF


      END DO Readline

    END IF   !! IF (trim(Rdline)=="END_BOUNDARY_CONDITIONS")

  END SUBROUTINE BC_parsing_subroutine

  !!!!!!!!!!!!!!!!!!!!! THE FOLLOWING SUBROUTINE COUNTS THE NUMBER OF TIMES A KEYWORD EXIST IN THE INPUT FILE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE  searching(Rdline,search_string,counter)

    !!!! initializing values


    CHARACTER(LEN=100)                         :: Rdline
    CHARACTER(LEN=*)                           :: search_string
    INTEGER                                    :: counter

    Rdline = to_upper(Rdline)
    Rdline =  strip_space(Rdline)
    IF (trim(Rdline)==search_string) counter = counter+1

  END SUBROUTINE  searching




  !!!!!!!!!!!REMOVE PROCEEDING , RECEEDING SPACES AND SPACES INBETWEEN CHARACTERS INSIDE A STRING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  FUNCTION strip_space( input_string ) RESULT ( output_string )

    ! -- Arguments
    CHARACTER( * ), INTENT( IN )  :: input_string
    INTEGER                       :: n

    ! -- Function result
    CHARACTER( LEN( input_string ) ) :: output_string

    ! -- Local parameters
    INTEGER,        PARAMETER :: IACHAR_SPACE = 32, IACHAR_TAB   = 9


    ! -- Local variables
    INTEGER :: i
    INTEGER :: iachar_CHARACTER

    ! -- Initialise output string
    output_string = ' '

    ! -- Initialise output string "useful" LENgth counter
    n = 0

    ! -- Loop over string elements
    DO i = 1, LEN( input_string )

      ! -- Convert the current CHARACTER to its position
      ! -- in the ASCII collating sequence
      iachar_CHARACTER = IACHAR( input_string( i:i ) )

      ! -- IF the CHARACTER is NOT a space ' ' or a tab '->|'
      ! -- copy it to the output string.
      IF ( iachar_CHARACTER /= IACHAR_SPACE .AND. &
           iachar_CHARACTER /= IACHAR_TAB) THEN
        n = n + 1
        output_string( n:n ) = input_string( i:i )
      END IF

    END DO


  END FUNCTION strip_space



  !!!!!!!!!!!!!!!!!!!!! THE FOLLOWING FUNCTION... !!!!!!!!!!!!!!!!!!!!!!
  !!!!!         REMOVE SPECIFIED CHARACTER FROM A STRING!!!!!!!!!!!!!!!!



  SUBROUTINE remove_CHARACTER(string,CHARACTER)
    implicit none
    INTEGER :: ind
    CHARACTER(LEN=100)  :: string
    CHARACTER           :: CHARACTER
    ind = SCAN(string, CHARACTER)
    IF (ind .NE.0) then
      string = string(:ind-1)
    END IF
  END SUBROUTINE remove_CHARACTER




  !!!!!CONVERT STRING TO REAL NUMBER !!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!



  REAL(CMISSRP) function str2real(str_2)
    implicit none
    ! Arguments
    CHARACTER(LEN=*)      :: str_2
    str_2 = to_upper(str_2)
    READ(str_2,*,IOSTAT=Stat)  str2real
  END function str2real


  !!!!!CONVERT STRING TO INTEGER NUMBER !!!!!!!!!!!!!! !!!!!!!!!!!!!!!!

  INTEGER function str2int(str)
    implicit none
    ! Arguments
    CHARACTER(LEN=*)      :: str
    str = to_upper(str)
    READ(str,*,IOSTAT=Stat)  str2int
  END function str2int

  !!!!!Skip blank lines in the input file !!!!!!!!!!!!!! !!!!!!!!!!!!!!!!

  SUBROUTINE skip_blank_lines(Rdline)


    INTEGER :: FileStat
    CHARACTER(LEN=100)              :: Rdline

    DO
      IF (Rdline .NE."") EXIT
      READ(12,*) Rdline
    END DO

  END SUBROUTINE skip_blank_lines

    !!!!!!!!!!!!!!!!!!!!!!!!!! The following SUBROUTINE splits a comma separated string !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE split_inline_argument_new(Rdline,array1,i)


    INTEGER :: k,POS2,POS1
    INTEGER,optional::i
    CHARACTER(LEN=100)  :: Rdline
    CHARACTER(LEN=*) , dimension(:) :: array1
    k=0
    DO

        POS1 = 1
        DO
            POS2 = INDEX(Rdline(POS1:), ",")
            IF (POS2 == 0) THEN
                k = k + 1
                array1(k) = Rdline(POS1:)
                EXIT
            END IF
            k = k + 1
            array1(k) = Rdline(POS1:POS1+POS2-2)
            POS1 = POS2+POS1
        END DO

        IF (Rdline == "") EXIT
        READ(12,'(A)') Rdline
        Rdline = to_upper(Rdline)
        Rdline =  strip_space(Rdline)
        CALL remove_CHARACTER(Rdline,"!")
        IF(present(i)) i=i+1

   END DO

  END SUBROUTINE split_inline_argument_new


  !!!!!!!!!!!!!!!!!!!!!!!!!! The following SUBROUTINE splits a comma separated string !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE split_inline_argument(Rdline,array1,i)


    INTEGER :: k,POS2,POS1
    INTEGER,optional::i
    CHARACTER(LEN=100)  :: Rdline
    CHARACTER(LEN=100) , dimension(:,:),ALLOCATABLE :: array1
    k=0
    DO

        POS1 = 1
        DO
            POS2 = INDEX(Rdline(POS1:), ",")
            IF (POS2 == 0) THEN
                k = k + 1
                array1(k,1) = Rdline(POS1:)
                EXIT
            END IF
            k = k + 1
            array1(k,1) = Rdline(POS1:POS1+POS2-2)
            POS1 = POS2+POS1
        END DO

        IF (Rdline == "") EXIT
        READ(12,'(A)') Rdline
        Rdline = to_upper(Rdline)
        Rdline =  strip_space(Rdline)
        CALL remove_CHARACTER(Rdline,"!")
        IF(present(i)) i=i+1
   END DO

  END SUBROUTINE split_inline_argument


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! THe following function converts lower case string to upper case !!!!!!!!!!!!!!!!!!!!!!!!

  function to_upper(strIn) result(strOut)
  ! Adapted from http://www.star.le.ac.uk/~cgp/fortran.html (25 May 2012)
  ! Original author: Clive Page

     implicit none

     character(len=*), intent(in) :: strIn
     character(len=len(strIn)) :: strOut
     integer :: i,j

     do i = 1, len(strIn)
          j = iachar(strIn(i:i))
          if (j>= iachar("a") .and. j<=iachar("z") ) then
               strOut(i:i) = achar(iachar(strIn(i:i))-32)
          else
               strOut(i:i) = strIn(i:i)
          end if
     end do

  end function to_upper

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! THE FOLLOWING FUNCTION  !!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!..LINK OBJECTS OF DIFFERENT !!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!....DERIVED TYPES WITH SAME IDs!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE MATCH_IDs(argument1,argument2,index1,string1,string2)

    CHARACTER(LEN=*)            :: argument1,argument2(:)
    INTEGER                     :: index1  , counter
    LOGICAL                     :: FLAG
    CHARACTER(LEN=*) , optional :: string1,string2

    FLAG = .TRUE.
    DO counter = 1,SIZE(argument2)

      IF (argument1 == argument2(counter)) then
        index1 = counter
        FLAG = .FALSE.
      END IF

   END DO
   !! if there are no similar ID found, then throw an error
   IF (FLAG) THEN
     CALL HANDLE_ERROR(string1 //"  and  "//string2//" does not contain objects with similar IDs")
   END IF
  END SUBROUTINE MATCH_IDs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! THE FOLLOWING FUNCTION  !!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!..used for error handling !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 SUBROUTINE HANDLE_ERROR(ERROR_STRING)

    CHARACTER(LEN=*), INTENT(IN) :: ERROR_STRING
    WRITE(*,'(">>ERROR: ",A)') ERROR_STRING(1:LEN_TRIM(ERROR_STRING))
    STOP

 END SUBROUTINE HANDLE_ERROR

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE parsing
