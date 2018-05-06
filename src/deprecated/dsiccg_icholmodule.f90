MODULE dsiccg_icholmodule
! Malak Olamaie
! 28/01/2013


!USE  params
USE  constants

CONTAINS
!DECK DCHKW
      SUBROUTINE DCHKW (NAME, LOCIW, LENIW, LOCW, LENW, IERR, ITER, ERR)
!***BEGIN PROLOGUE  DCHKW
!***SUBSIDIARY
!***PURPOSE  SLAP WORK/IWORK Array Bounds Checker.
!            This routine checks the work array lengths and interfaces
!            to the SLATEC error handler if a problem is found.
!***LIBRARY   SLATEC (SLAP)
!***CATEGORY  R2
!***TYPE      DOUBLE PRECISION (SCHKW-S, DCHKW-D)
!***KEYWORDS  ERROR CHECKING, SLAP, WORKSPACE CHECKING
!***AUTHOR  Seager, Mark K., (LLNL)
!             Lawrence Livermore National Laboratory
!             PO BOX 808, L-60
!             Livermore, CA 94550 (510) 423-3141
!             seager@llnl.gov
!***DESCRIPTION
!
! *Usage:
!     CHARACTER*(*) NAME
!     INTEGER LOCIW, LENIW, LOCW, LENW, IERR, ITER
!     DOUBLE PRECISION ERR
!
!     CALL DCHKW( NAME, LOCIW, LENIW, LOCW, LENW, IERR, ITER, ERR )
!
! *Arguments:
! NAME   :IN       Character*(*).
!         Name of the calling routine.  This is used in the output
!         message, if an error is detected.
! LOCIW  :IN       Integer.
!         Location of the first free element in the integer workspace
!         array.
! LENIW  :IN       Integer.
!         Length of the integer workspace array.
! LOCW   :IN       Integer.
!         Location of the first free element in the double precision
!         workspace array.
! LENRW  :IN       Integer.
!         Length of the double precision workspace array.
! IERR   :OUT      Integer.
!         Return error flag.
!               IERR = 0 => All went well.
!               IERR = 1 => Insufficient storage allocated for
!                           WORK or IWORK.
! ITER   :OUT      Integer.
!         Set to zero on return.
! ERR    :OUT      Double Precision.
!         Set to the smallest positive magnitude if all went well.
!         Set to a very large number if an error is detected.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  D1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   880225  DATE WRITTEN
!   881213  Previous REVISION DATE
!   890915  Made changes requested at July 1989 CML Meeting.  (MKS)
!   890922  Numerous changes to prologue to make closer to SLATEC
!           standard.  (FNF)
!   890929  Numerous changes to reduce SP/DP differences.  (FNF)
!   900805  Changed XERRWV calls to calls to XERMSG.  (RWC)
!   910411  Prologue converted to Version 4.0 format.  (BAB)
!   910502  Corrected XERMSG calls to satisfy Section 6.2.2 of ANSI
!           X3.9-1978.  (FNF)
!   910506  Made subsidiary.  (FNF)
!   920511  Added complete declaration section.  (WRB)
!   921015  Added code to initialize ITER and ERR when IERR=0.  (FNF)
!***END PROLOGUE  DCHKW
!     .. Scalar Arguments ..
      DOUBLE PRECISION ERR
      INTEGER IERR, ITER, LENIW, LENW, LOCIW, LOCW
      CHARACTER NAME*(*)
!     .. Local Scalars ..
      CHARACTER XERN1*8, XERN2*8, XERNAM*8
!     .. External Functions ..
 !     DOUBLE PRECISION D1MACH
!      EXTERNAL D1MACH
!     .. External Subroutines ..
!      EXTERNAL XERMSG
!***FIRST EXECUTABLE STATEMENT  DCHKW
!
!         Check the Integer workspace situation.
!
      IERR = 0
      ITER = 0
      ERR = D1MACH(1)
      IF( LOCIW.GT.LENIW ) THEN
         IERR = 1
         ERR = D1MACH(2)
         XERNAM = NAME
         WRITE (XERN1, '(I8)') LOCIW
         WRITE (XERN2, '(I8)') LENIW
         CALL XERMSG ('SLATEC', 'DCHKW',&
           'In ' // XERNAM // ', INTEGER work array too short.  ' //&
           'IWORK needs ' // XERN1 // '; have allocated ' // XERN2,&
           1, 1)
      ENDIF
!
!         Check the Double Precision workspace situation.
      IF( LOCW.GT.LENW ) THEN
         IERR = 1
         ERR = D1MACH(2)
         XERNAM = NAME
         WRITE (XERN1, '(I8)') LOCW
         WRITE (XERN2, '(I8)') LENW
         CALL XERMSG ('SLATEC', 'DCHKW',&
           'In ' // XERNAM // ', DOUBLE PRECISION work array too ' //&
           'short.  RWORK needs ' // XERN1 // '; have allocated ' //&
           XERN2, 1, 1)
      ENDIF
      RETURN
      END SUBROUTINE DCHKW
!===================================================================
!DECK D1MACH
        FUNCTION D1MACH(I)
!
! D1MACH returns double precision real machine constants.
!
!  Discussion:
!
!    Assuming that the internal representation of a double precision real
!    number is in base B, with T the number of base-B digits in the mantissa,
!    and EMIN the smallest possible exponent and EMAX the largest possible 
!    exponent, then
!
!      D1MACH(1) = B^(EMIN-1), the smallest positive magnitude.
!      D1MACH(2) = B^EMAX*(1-B^(-T)), the largest magnitude.
!      D1MACH(3) = B^(-T), the smallest relative spacing.
!      D1MACH(4) = B^(1-T), the largest relative spacing.
!      D1MACH(5) = log10(B).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 April 2007
!
!  Author:
!
!    Original FORTRAN77 version by Phyllis Fox, Andrew Hall, Norman Schryer.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Phyllis Fox, Andrew Hall, Norman Schryer,
!    Algorithm 528:
!    Framework for a Portable Library,
!    ACM Transactions on Mathematical Software,
!    Volume 4, Number 2, June 1978, page 176-188.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, chooses the parameter to be returned.
!    1 <= I <= 5.
!
!    Output, real ( kind = 8 ) D1MACH, the value of the chosen parameter.
!
! Now uses fortran90 functions, Malak Olamaie January 2013
         IMPLICIT NONE
       
         INTEGER I
         double precision   D1MACH
         double precision   X,XX

         X=1.d0
        
         IF ( I < 1 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'D1MACH - Fatal error!'
            write ( *, '(a)' ) '  The input argument I is out of bounds.'
            write ( *, '(a)' ) '  Legal values satisfy 1 <= I <= 5.'
            write ( *, '(a,i12)' ) '  I = ', I
            D1MACH = 0.0D+00
            STOP
         ELSEIF( I == 1 ) THEN
             D1MACH = TINY(X)
             
         ELSEIF( I == 2 ) THEN
             D1MACH = HUGE(X)
             
         ELSEIF( I == 3 ) THEN
             D1MACH = EPSILON(X)/RADIX(X)
         ELSEIF( I == 4 ) THEN
             D1MACH = EPSILON(X) 
         ELSEIF( I == 5 ) THEN
             XX=RADIX(X)
             D1MACH = DLOG10(XX)
         ELSEIF( 5 < I ) THEN
             write ( *, '(a)' ) ' '
             write ( *, '(a)' ) 'D1MACH - Fatal error!'
             write ( *, '(a)' ) '  The input argument I is out of bounds.'
             write ( *, '(a)' ) '  Legal values satisfy 1 <= I <= 5.'
             write ( *, '(a,i12)' ) '  I = ', i
             D1MACH = 0.0D+00
             STOP
         ENDIF

        RETURN
        END FUNCTION D1MACH
!=================================================================================================
!DECK XERMSG
      SUBROUTINE XERMSG (LIBRAR, SUBROU, MESSG, NERR, LEVEL)
! ***BEGIN PROLOGUE  XERMSG
! ***PURPOSE  Process error messages for SLATEC and other libraries.
! ***LIBRARY   SLATEC (XERROR)
! ***CATEGORY  R3C
! ***TYPE      ALL (XERMSG-A)
! ***KEYWORDS  ERROR MESSAGE, XERROR
! ***AUTHOR  Fong, Kirby, (NMFECC at LLNL)
! ***DESCRIPTION
! 
!   XERMSG processes a diagnostic message in a manner determined by the
!   value of LEVEL and the current value of the library error control
!   flag, KONTRL.  See subroutine XSETF for details.
! 
!    LIBRAR   A character constant (or character variable) with the name
!             of the library.  This will be 'SLATEC' for the SLATEC
!             Common Math Library.  The error handling package is
!             general enough to be used by many libraries
!             simultaneously, so it is desirable for the routine that
!             detects and reports an error to identify the library name
!             as well as the routine name.
! 
!    SUBROU   A character constant (or character variable) with the name
!             of the routine that detected the error.  Usually it is the
!             name of the routine that is calling XERMSG.  There are
!             some instances where a user callable library routine calls
!             lower level subsidiary routines where the error is
!             detected.  In such cases it may be more informative to
!             supply the name of the routine the user called rather than
!             the name of the subsidiary routine that detected the
!             error.
! 
!    MESSG    A character constant (or character variable) with the text
!             of the error or warning message.  In the example below,
!             the message is a character constant that contains a
!             generic message.
! 
!                   CALL XERMSG ('SLATEC', 'MMPY',
!                  *'THE ORDER OF THE MATRIX EXCEEDS THE ROW DIMENSION',
!                  *3, 1)
! 
!             It is possible (and is sometimes desirable) to generate a
!             specific message--e.g., one that contains actual numeric
!             values.  Specific numeric values can be converted into
!             character strings using formatted WRITE statements into
!             character variables.  This is called standard Fortran
!             internal file I/O and is exemplified in the first three
!             lines of the following example.  You can also catenate
!             substrings of characters to construct the error message.
!             Here is an example showing the use of both writing to
!             an internal file and catenating character strings.
! 
!                   CHARACTER*5 CHARN, CHARL
!                   WRITE (CHARN,10) N
!                   WRITE (CHARL,10) LDA
!                10 FORMAT(I5)
!                   CALL XERMSG ('SLATEC', 'MMPY', 'THE ORDER'//CHARN//
!                  *   ' OF THE MATRIX EXCEEDS ITS ROW DIMENSION OF'//
!                  *   CHARL, 3, 1)
! 
!             There are two subtleties worth mentioning.  One is that
!             the // for character catenation is used to construct the
!             error message so that no single character constant is
!             continued to the next line.  This avoids confusion as to
!             whether there are trailing blanks at the end of the line.
!             The second is that by catenating the parts of the message
!             as an actual argument rather than encoding the entire
!             message into one large character variable, we avoid
!             having to know how long the message will be in order to
!             declare an adequate length for that large character
!             variable.  XERMSG calls XERPRN to print the message using
!             multiple lines if necessary.  If the message is very long,
!             XERPRN will break it into pieces of 72 characters (as
!             requested by XERMSG) for printing on multiple lines.
!             Also, XERMSG asks XERPRN to prefix each line with ' *  '
!             so that the total line length could be 76 characters.
!             Note also that XERPRN scans the error message backwards
!             to ignore trailing blanks.  Another feature is that
!             the substring '$$' is treated as a new line sentinel
!             by XERPRN.  If you want to construct a multiline
!             message without having to count out multiples of 72
!             characters, just use '$$' as a separator.  '$$'
!             obviously must occur within 72 characters of the
!             start of each line to have its intended effect since
!             XERPRN is asked to wrap around at 72 characters in
!             addition to looking for '$$'.
! 
!    NERR     An integer value that is chosen by the library routine's
!             author.  It must be in the range -99 to 999 (three
!             printable digits).  Each distinct error should have its
!             own error number.  These error numbers should be described
!             in the machine readable documentation for the routine.
!             The error numbers need be unique only within each routine,
!             so it is reasonable for each routine to start enumerating
!             errors from 1 and proceeding to the next integer.
! 
!    LEVEL    An integer value in the range 0 to 2 that indicates the
!             level (severity) of the error.  Their meanings are
! 
!            -1  A warning message.  This is used if it is not clear
!                that there really is an error, but the user's attention
!                may be needed.  An attempt is made to only print this
!                message once.
! 
!             0  A warning message.  This is used if it is not clear
!                that there really is an error, but the user's attention
!                may be needed.
! 
!             1  A recoverable error.  This is used even if the error is
!                so serious that the routine cannot return any useful
!                answer.  If the user has told the error package to
!                return after recoverable errors, then XERMSG will
!                return to the Library routine which can then return to
!                the user's routine.  The user may also permit the error
!                package to terminate the program upon encountering a
!                recoverable error.
! 
!             2  A fatal error.  XERMSG will not return to its caller
!                after it receives a fatal error.  This level should
!                hardly ever be used; it is much better to allow the
!                user a chance to recover.  An example of one of the few
!                cases in which it is permissible to declare a level 2
!                error is a reverse communication Library routine that
!                is likely to be called repeatedly until it integrates
!                across some interval.  If there is a serious error in
!                the input such that another step cannot be taken and
!                the Library routine is called again without the input
!                error having been corrected by the caller, the Library
!                routine will probably be called forever with improper
!                input.  In this case, it is reasonable to declare the
!                error to be fatal.
! 
!    Each of the arguments to XERMSG is input; none will be modified by
!    XERMSG.  A routine may make multiple calls to XERMSG with warning
!    level messages; however, after a call to XERMSG with a recoverable
!    error, the routine should return to the user.  Do not try to call
!    XERMSG with a second recoverable error after the first recoverable
!    error because the error package saves the error number.  The user
!    can retrieve this error number by calling another entry point in
!    the error handling package and then clear the error number when
!    recovering from the error.  Calling XERMSG in succession causes the
!    old error number to be overwritten by the latest error number.
!    This is considered harmless for error numbers associated with
!    warning messages but must not be done for error numbers of serious
!    errors.  After a call to XERMSG with a recoverable error, the user
!    must be given a chance to call NUMXER or XERCLR to retrieve or
!    clear the error number.
! ***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
!                 Error-handling Package, SAND82-0800, Sandia
!                 Laboratories, 1982.
! ***ROUTINES CALLED  FDUMP, J4SAVE, XERCNT, XERHLT, XERPRN, XERSVE
! ***REVISION HISTORY  (YYMMDD)
!   880101  DATE WRITTEN
!   880621  REVISED AS DIRECTED AT SLATEC CML MEETING OF FEBRUARY 1988.
!           THERE ARE TWO BASIC CHANGES.
!           1.  A NEW ROUTINE, XERPRN, IS USED INSTEAD OF XERPRT TO
!               PRINT MESSAGES.  THIS ROUTINE WILL BREAK LONG MESSAGES
!               INTO PIECES FOR PRINTING ON MULTIPLE LINES.  '$$' IS
!               ACCEPTED AS A NEW LINE SENTINEL.  A PREFIX CAN BE
!               ADDED TO EACH LINE TO BE PRINTED.  XERMSG USES EITHER
!               ' ***' OR ' *  ' AND LONG MESSAGES ARE BROKEN EVERY
!               72 CHARACTERS (AT MOST) SO THAT THE MAXIMUM LINE
!               LENGTH OUTPUT CAN NOW BE AS GREAT AS 76.
!           2.  THE TEXT OF ALL MESSAGES IS NOW IN UPPER CASE SINCE THE
!               FORTRAN STANDARD DOCUMENT DOES NOT ADMIT THE EXISTENCE
!               OF LOWER CASE.
!   880708  REVISED AFTER THE SLATEC CML MEETING OF JUNE 29 AND 30.
!           THE PRINCIPAL CHANGES ARE
!           1.  CLARIFY COMMENTS IN THE PROLOGUES
!           2.  RENAME XRPRNT TO XERPRN
!           3.  REWORK HANDLING OF '$$' IN XERPRN TO HANDLE BLANK LINES
!               SIMILAR TO THE WAY FORMAT STATEMENTS HANDLE THE /
!               CHARACTER FOR NEW RECORDS.
!   890706  REVISED WITH THE HELP OF FRED FRITSCH AND REG CLEMENS TO
!           CLEAN UP THE CODING.
!   890721  REVISED TO USE NEW FEATURE IN XERPRN TO COUNT CHARACTERS IN
!           PREFIX.
!   891013  REVISED TO CORRECT COMMENTS.
!   891214  Prologue converted to Version 4.0 format.  (WRB)
!   900510  Changed test on NERR to be -9999999 < NERR < 99999999, but
!           NERR .ne. 0, and on LEVEL to be -2 < LEVEL < 3.  Added
!           LEVEL=-1 logic, changed calls to XERSAV to XERSVE, and
!           XERCTL to XERCNT.  (RWC)
!   920501  Reformatted the REFERENCES section.  (WRB)
! ***END PROLOGUE  XERMSG
      INTEGER NERR,LEVEL
      CHARACTER*(*) LIBRAR, SUBROU, MESSG
      CHARACTER*8 XLIBR, XSUBR
      CHARACTER*72  TEMP
      CHARACTER*20  LFIRST
      INTEGER LKNTRL, MAXMES, KDUMMY, I, KOUNT, LERR, LLEVEL, MKNTRL, &
          LTEMP, NUNIT, KUNIT
! ***FIRST EXECUTABLE STATEMENT  XERMSG
      LKNTRL = J4SAVE (2, 0, .FALSE.)
      MAXMES = J4SAVE (4, 0, .FALSE.)
! 
!       LKNTRL IS A LOCAL COPY OF THE CONTROL FLAG KONTRL.
!       MAXMES IS THE MAXIMUM NUMBER OF TIMES ANY PARTICULAR MESSAGE
!          SHOULD BE PRINTED.
! 
!       WE PRINT A FATAL ERROR MESSAGE AND TERMINATE FOR AN ERROR IN
!          CALLING XERMSG.  THE ERROR NUMBER SHOULD BE POSITIVE,
!          AND THE LEVEL SHOULD BE BETWEEN 0 AND 2.
! 

         IF (NERR.LT.-9999999 .OR. NERR.GT.99999999 .OR. NERR.EQ.0 .OR.&
             LEVEL.LT.-1 .OR. LEVEL.GT.2) THEN
         CALL XERPRN (' ***', -1, 'FATAL ERROR IN...$$ ' //&
           'XERMSG -- INVALID ERROR NUMBER OR LEVEL$$ '//&
           'JOB ABORT DUE TO FATAL ERROR.', 72)
         CALL XERSVE (' ', ' ', ' ', 0, 0, 0, KDUMMY)
         CALL XERHLT (' ***XERMSG -- INVALID INPUT')
         RETURN
      ENDIF
! 
!       RECORD THE MESSAGE.
! 
      I = J4SAVE (1, NERR, .TRUE.)
      CALL XERSVE (LIBRAR, SUBROU, MESSG, 1, NERR, LEVEL, KOUNT)
! 
!       HANDLE PRINT-ONCE WARNING MESSAGES.
! 
      IF (LEVEL.EQ.-1 .AND. KOUNT.GT.1) RETURN
! 
!       ALLOW TEMPORARY USER OVERRIDE OF THE CONTROL FLAG.
! 
      XLIBR  = LIBRAR
      XSUBR  = SUBROU
      LFIRST = MESSG
      LERR   = NERR
      LLEVEL = LEVEL
      CALL XERCNT (XLIBR, XSUBR, LFIRST, LERR, LLEVEL, LKNTRL)
! 
      LKNTRL = MAX(-2, MIN(2,LKNTRL))
      MKNTRL = ABS(LKNTRL)
! 
!       SKIP PRINTING IF THE CONTROL FLAG VALUE AS RESET IN XERCNT IS
!       ZERO AND THE ERROR IS NOT FATAL.
! 
      IF (LEVEL.LT.2 .AND. LKNTRL.EQ.0) GO TO 30
      IF (LEVEL.EQ.0 .AND. KOUNT.GT.MAXMES) GO TO 30
      IF (LEVEL.EQ.1 .AND. KOUNT.GT.MAXMES .AND. MKNTRL.EQ.1) GO TO 30
      IF (LEVEL.EQ.2 .AND. KOUNT.GT.MAX(1,MAXMES)) GO TO 30
! 
!       ANNOUNCE THE NAMES OF THE LIBRARY AND SUBROUTINE BY BUILDING A
!       MESSAGE IN CHARACTER VARIABLE TEMP (NOT EXCEEDING 66 CHARACTERS)
!       AND SENDING IT OUT VIA XERPRN.  PRINT ONLY IF CONTROL FLAG
!       IS NOT ZERO.
! 
      IF (LKNTRL .NE. 0) THEN
         TEMP(1:21) = 'MESSAGE FROM ROUTINE '
         I = MIN(LEN(SUBROU), 16)
         TEMP(22:21+I) = SUBROU(1:I)
         TEMP(22+I:33+I) = ' IN LIBRARY '
         LTEMP = 33 + I
         I = MIN(LEN(LIBRAR), 16)
         TEMP(LTEMP+1:LTEMP+I) = LIBRAR (1:I)
         TEMP(LTEMP+I+1:LTEMP+I+1) = '.'
         LTEMP = LTEMP + I + 1
         CALL XERPRN (' ***', -1, TEMP(1:LTEMP), 72)
      ENDIF
! 
!       IF LKNTRL IS POSITIVE, PRINT AN INTRODUCTORY LINE BEFORE
!       PRINTING THE MESSAGE.  THE INTRODUCTORY LINE TELLS THE CHOICE
!       FROM EACH OF THE FOLLOWING THREE OPTIONS.
!       1.  LEVEL OF THE MESSAGE
!              'INFORMATIVE MESSAGE'
!              'POTENTIALLY RECOVERABLE ERROR'
!              'FATAL ERROR'
!       2.  WHETHER CONTROL FLAG WILL ALLOW PROGRAM TO CONTINUE
!              'PROG CONTINUES'
!              'PROG ABORTED'
!       3.  WHETHER OR NOT A TRACEBACK WAS REQUESTED.  (THE TRACEBACK
!           MAY NOT BE IMPLEMENTED AT SOME SITES, SO THIS ONLY TELLS
!           WHAT WAS REQUESTED, NOT WHAT WAS DELIVERED.)
!              'TRACEBACK REQUESTED'
!              'TRACEBACK NOT REQUESTED'
!       NOTICE THAT THE LINE INCLUDING FOUR PREFIX CHARACTERS WILL NOT
!       EXCEED 74 CHARACTERS.
!       WE SKIP THE NEXT BLOCK IF THE INTRODUCTORY LINE IS NOT NEEDED.
! 
      IF (LKNTRL .GT. 0) THEN
! 
!       THE FIRST PART OF THE MESSAGE TELLS ABOUT THE LEVEL.
! 
         IF (LEVEL .LE. 0) THEN
            TEMP(1:20) = 'INFORMATIVE MESSAGE,'
            LTEMP = 20
         ELSEIF (LEVEL .EQ. 1) THEN
            TEMP(1:30) = 'POTENTIALLY RECOVERABLE ERROR,'
            LTEMP = 30
         ELSE
            TEMP(1:12) = 'FATAL ERROR,'
            LTEMP = 12
         ENDIF
! 
!       THEN WHETHER THE PROGRAM WILL CONTINUE.
! 
         IF ((MKNTRL.EQ.2 .AND. LEVEL.GE.1) .OR.&
            (MKNTRL.EQ.1 .AND. LEVEL.EQ.2)) THEN
            TEMP(LTEMP+1:LTEMP+14) = ' PROG ABORTED,'
            LTEMP = LTEMP + 14
         ELSE
            TEMP(LTEMP+1:LTEMP+16) = ' PROG CONTINUES,'
            LTEMP = LTEMP + 16
         ENDIF
! 
!       FINALLY TELL WHETHER THERE SHOULD BE A TRACEBACK.
! 
         IF (LKNTRL .GT. 0) THEN
            TEMP(LTEMP+1:LTEMP+20) = ' TRACEBACK REQUESTED'
            LTEMP = LTEMP + 20
         ELSE
            TEMP(LTEMP+1:LTEMP+24) = ' TRACEBACK NOT REQUESTED'
            LTEMP = LTEMP + 24
         ENDIF
         CALL XERPRN (' ***', -1, TEMP(1:LTEMP), 72)
      ENDIF
! 
!       NOW SEND OUT THE MESSAGE.
! 
      CALL XERPRN (' *  ', -1, MESSG, 72)
! 
!       IF LKNTRL IS POSITIVE, WRITE THE ERROR NUMBER AND REQUEST A
!          TRACEBACK.
! 
      IF (LKNTRL .GT. 0) THEN
         WRITE (TEMP, '(''ERROR NUMBER = '', I8)') NERR
         DO 10 I=16,22
            IF (TEMP(I:I) .NE. ' ') GO TO 20
   10    CONTINUE
! 
   20    CALL XERPRN (' *  ', -1, TEMP(1:15) // TEMP(I:23), 72)
         CALL FDUMP
      ENDIF
! 
!       IF LKNTRL IS NOT ZERO, PRINT A BLANK LINE AND AN END OF MESSAGE.
! 
      IF (LKNTRL .NE. 0) THEN
         CALL XERPRN (' *  ', -1, ' ', 72)
         CALL XERPRN (' ***', -1, 'END OF MESSAGE', 72)
         CALL XERPRN ('    ',  0, ' ', 72)
      ENDIF
! 
!       IF THE ERROR IS NOT FATAL OR THE ERROR IS RECOVERABLE AND THE
!       CONTROL FLAG IS SET FOR RECOVERY, THEN RETURN.
! 
   30 IF (LEVEL.LE.0 .OR. (LEVEL.EQ.1 .AND. MKNTRL.LE.1)) RETURN
! 
!       THE PROGRAM WILL BE STOPPED DUE TO AN UNRECOVERED ERROR OR A
!       FATAL ERROR.  PRINT THE REASON FOR THE ABORT AND THE ERROR
!       SUMMARY IF THE CONTROL FLAG AND THE MAXIMUM ERROR COUNT PERMIT.
! 
      IF (LKNTRL.GT.0 .AND. KOUNT.LT.MAX(1,MAXMES)) THEN
         IF (LEVEL .EQ. 1) THEN
            CALL XERPRN &
              (' ***', -1, 'JOB ABORT DUE TO UNRECOVERED ERROR.', 72)
         ELSE
            CALL XERPRN(' ***', -1, 'JOB ABORT DUE TO FATAL ERROR.', 72)
         ENDIF
         CALL XERSVE (' ', ' ', ' ', -1, 0, 0, KDUMMY)
         CALL XERHLT (' ')
      ELSE
         CALL XERHLT (MESSG)
      ENDIF
      RETURN
      END SUBROUTINE XERMSG
!=====================================================================================
!DECK J4SAVE
      INTEGER FUNCTION J4SAVE (IWHICH, IVALUE, ISET)
! ***BEGIN PROLOGUE  J4SAVE
! ***SUBSIDIARY
! ***PURPOSE  Save or recall global variables needed by error
!            handling routines.
! ***LIBRARY   SLATEC (XERROR)
! ***TYPE      INTEGER (J4SAVE-I)
! ***KEYWORDS  ERROR MESSAGES, ERROR NUMBER, RECALL, SAVE, XERROR
! ***AUTHOR  Jones, R. E., (SNLA)
! ***DESCRIPTION
! 
!      Abstract
!        J4SAVE saves and recalls several global variables needed
!        by the library error handling routines.
! 
!     Description of Parameters
!      --Input--
!        IWHICH - Index of item desired.
!                = 1 Refers to current error number.
!                = 2 Refers to current error control flag.
!                = 3 Refers to current unit number to which error
!                    messages are to be sent.  (0 means use standard.)
!                = 4 Refers to the maximum number of times any
!                     message is to be printed (as set by XERMAX).
!                = 5 Refers to the total number of units to which
!                     each error message is to be written.
!                = 6 Refers to the 2nd unit for error messages
!                = 7 Refers to the 3rd unit for error messages
!                = 8 Refers to the 4th unit for error messages
!                = 9 Refers to the 5th unit for error messages
!        IVALUE - The value to be set for the IWHICH-th parameter,
!                 if ISET is .TRUE. .
!        ISET   - If ISET=.TRUE., the IWHICH-th parameter will BE
!                 given the value, IVALUE.  If ISET=.FALSE., the
!                 IWHICH-th parameter will be unchanged, and IVALUE
!                 is a dummy parameter.
!      --Output--
!        The (old) value of the IWHICH-th parameter will be returned
!        in the function value, J4SAVE.
! 
! ***SEE ALSO  XERMSG
! ***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
!                 Error-handling Package, SAND82-0800, Sandia
!                 Laboratories, 1982.
! ***ROUTINES CALLED  (NONE)
! ***REVISION HISTORY  (YYMMDD)
!   790801  DATE WRITTEN
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900205  Minor modifications to prologue.  (WRB)
!   900402  Added TYPE section.  (WRB)
!   910411  Added KEYWORDS section.  (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
! ***END PROLOGUE  J4SAVE
      INTEGER IWHICH,IVALUE
      LOGICAL ISET
      INTEGER IPARAM(9)
      SAVE IPARAM
      DATA IPARAM(1),IPARAM(2),IPARAM(3),IPARAM(4)/0,2,0,10/
      DATA IPARAM(5)/1/
      DATA IPARAM(6),IPARAM(7),IPARAM(8),IPARAM(9)/0,0,0,0/
! ***FIRST EXECUTABLE STATEMENT  J4SAVE
      J4SAVE = IPARAM(IWHICH)
      IF (ISET) IPARAM(IWHICH) = IVALUE
      RETURN
      END FUNCTION J4SAVE
!============================================================================
!DECK XERPRN
      SUBROUTINE XERPRN (PREFIX, NPREF, MESSG, NWRAP)
! ***BEGIN PROLOGUE  XERPRN
! ***SUBSIDIARY
! ***PURPOSE  Print error messages processed by XERMSG.
! ***LIBRARY   SLATEC (XERROR)
! ***CATEGORY  R3C
! ***TYPE      ALL (XERPRN-A)
! ***KEYWORDS  ERROR MESSAGES, PRINTING, XERROR
! ***AUTHOR  Fong, Kirby, (NMFECC at LLNL)
! ***DESCRIPTION
! 
! This routine sends one or more lines to each of the (up to five)
! logical units to which error messages are to be sent.  This routine
! is called several times by XERMSG, sometimes with a single line to
! print and sometimes with a (potentially very long) message that may
! wrap around into multiple lines.
! 
! PREFIX  Input argument of type CHARACTER.  This argument contains
!         characters to be put at the beginning of each line before
!         the body of the message.  No more than 16 characters of
!         PREFIX will be used.
! 
! NPREF   Input argument of type INTEGER.  This argument is the number
!         of characters to use from PREFIX.  If it is negative, the
!         intrinsic function LEN is used to determine its length.  If
!         it is zero, PREFIX is not used.  If it exceeds 16 or if
!         LEN(PREFIX) exceeds 16, only the first 16 characters will be
!         used.  If NPREF is positive and the length of PREFIX is less
!         than NPREF, a copy of PREFIX extended with blanks to length
!         NPREF will be used.
! 
! MESSG   Input argument of type CHARACTER.  This is the text of a
!         message to be printed.  If it is a long message, it will be
!         broken into pieces for printing on multiple lines.  Each line
!         will start with the appropriate prefix and be followed by a
!         piece of the message.  NWRAP is the number of characters per
!         piece; that is, after each NWRAP characters, we break and
!         start a new line.  In addition the characters '$$' embedded
!         in MESSG are a sentinel for a new line.  The counting of
!         characters up to NWRAP starts over for each new line.  The
!         value of NWRAP typically used by XERMSG is 72 since many
!         older error messages in the SLATE! Library are laid out to
!         rely on wrap-around every 72 characters.
! 
! NWRAP   Input argument of type INTEGER.  This gives the maximum size
!         piece into which to break MESSG for printing on multiple
!         lines.  An embedded '$$' ends a line, and the count restarts
!         at the following character.  If a line break does not occur
!         on a blank (it would split a word) that word is moved to the
!         next line.  Values of NWRAP less than 16 will be treated as
!         16.  Values of NWRAP greater than 132 will be treated as 132.
!         The actual line length will be NPREF + NWRAP after NPREF has
!         been adjusted to fall between 0 and 16 and NWRAP has been
!         adjusted to fall between 16 and 132.
! 
! ***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
!                 Error-handling Package, SAND82-0800, Sandia
!                 Laboratories, 1982.
! ***ROUTINES CALLED  I1MACH, XGETUA
! ***REVISION HISTORY  (YYMMDD)
!   880621  DATE WRITTEN
!   880708  REVISED AFTER THE SLATEC CML SUBCOMMITTEE MEETING OF
!           JUNE 29 AND 30 TO CHANGE THE NAME TO XERPRN AND TO REWORK
!           THE HANDLING OF THE NEW LINE SENTINEL TO BEHAVE LIKE THE
!           SLASH CHARACTER IN FORMAT STATEMENTS.
!   890706  REVISED WITH THE HELP OF FRED FRITSCH AND REG CLEMENS TO
!           STREAMLINE THE CODING AND FIX A BUG THAT CAUSED EXTRA BLANK
!           LINES TO BE PRINTED.
!   890721  REVISED TO ADD A NEW FEATURE.  A NEGATIVE VALUE OF NPREF
!           CAUSES LEN(PREFIX) TO BE USED AS THE LENGTH.
!   891013  REVISED TO CORRECT ERROR IN CALCULATING PREFIX LENGTH.
!   891214  Prologue converted to Version 4.0 format.  (WRB)
!   900510  Added code to break messages between words.  (RWC)
!   920501  Reformatted the REFERENCES section.  (WRB)
! ***END PROLOGUE  XERPRN
      CHARACTER*(*) PREFIX, MESSG
      INTEGER NPREF, NWRAP , N
      CHARACTER*148 CBUFF
      INTEGER IU(5), NUNIT
      CHARACTER*2 NEWLIN
      PARAMETER (NEWLIN = '$$')
      INTEGER I, LPREF, LWRAP, LENMSG, NEXTC, LPIECE, IDELTA
! ***FIRST EXECUTABLE STATEMENT  XERPRN
      CALL XGETUA(IU,NUNIT)
! 
!       A ZERO VALUE FOR A LOGICAL UNIT NUMBER MEANS TO USE THE STANDARD
!       ERROR MESSAGE UNIT INSTEAD.  I1MACH(4) RETRIEVES THE STANDARD
!       ERROR MESSAGE UNIT.
! 
      N = I1MACH(4)
      DO 10 I=1,NUNIT
         IF (IU(I) .EQ. 0) IU(I) = N
   10 CONTINUE
! 
!       LPREF IS THE LENGTH OF THE PREFIX.  THE PREFIX IS PLACED AT THE
!       BEGINNING OF CBUFF, THE CHARACTER BUFFER, AND KEPT THERE DURING
!       THE REST OF THIS ROUTINE.
! 
      IF ( NPREF .LT. 0 ) THEN
         LPREF = LEN(PREFIX)
      ELSE
         LPREF = NPREF
      ENDIF
      LPREF = MIN(16, LPREF)
      IF (LPREF .NE. 0) CBUFF(1:LPREF) = PREFIX
! 
!       LWRAP IS THE MAXIMUM NUMBER OF CHARACTERS WE WANT TO TAKE AT ONE
!       TIME FROM MESSG TO PRINT ON ONE LINE.
! 
      LWRAP = MAX(16, MIN(132, NWRAP))
! 
!       SET LENMSG TO THE LENGTH OF MESSG, IGNORE ANY TRAILING BLANKS.
! 
      LENMSG = LEN(MESSG)
      N = LENMSG
      DO 20 I=1,N
         IF (MESSG(LENMSG:LENMSG) .NE. ' ') GO TO 30
         LENMSG = LENMSG - 1
   20 CONTINUE
   30 CONTINUE
! 
!       IF THE MESSAGE IS ALL BLANKS, THEN PRINT ONE BLANK LINE.
! 
      IF (LENMSG .EQ. 0) THEN
         CBUFF(LPREF+1:LPREF+1) = ' '
         DO 40 I=1,NUNIT
            WRITE(IU(I), '(A)') CBUFF(1:LPREF+1)
   40    CONTINUE
         RETURN
      ENDIF
! 
!       SET NEXTC TO THE POSITION IN MESSG WHERE THE NEXT SUBSTRING
!       STARTS.  FROM THIS POSITION WE SCAN FOR THE NEW LINE SENTINEL.
!       WHEN NEXTC EXCEEDS LENMSG, THERE IS NO MORE TO PRINT.
!       WE LOOP BACK TO LABEL 50 UNTIL ALL PIECES HAVE BEEN PRINTED.
! 
!       WE LOOK FOR THE NEXT OCCURRENCE OF THE NEW LINE SENTINEL.  THE
!       INDEX INTRINSIC FUNCTION RETURNS ZERO IF THERE IS NO OCCURRENCE
!       OR IF THE LENGTH OF THE FIRST ARGUMENT IS LESS THAN THE LENGTH
!       OF THE SECOND ARGUMENT.
! 
!       THERE ARE SEVERAL CASES WHICH SHOULD BE CHECKED FOR IN THE
!       FOLLOWING ORDER.  WE ARE ATTEMPTING TO SET LPIECE TO THE NUMBER
!       OF CHARACTERS THAT SHOULD BE TAKEN FROM MESSG STARTING AT
!       POSITION NEXTC.
! 
!       LPIECE .EQ. 0   THE NEW LINE SENTINEL DOES NOT OCCUR IN THE
!                       REMAINDER OF THE CHARACTER STRING.  LPIECE
!                       SHOULD BE SET TO LWRAP OR LENMSG+1-NEXTC,
!                       WHICHEVER IS LESS.
! 
!       LPIECE .EQ. 1   THE NEW LINE SENTINEL STARTS AT MESSG(NEXTC:
!                       NEXTC).  LPIECE IS EFFECTIVELY ZERO, AND WE
!                       PRINT NOTHING TO AVOID PRODUCING UNNECESSARY
!                       BLANK LINES.  THIS TAKES CARE OF THE SITUATION
!                       WHERE THE LIBRARY ROUTINE HAS A MESSAGE OF
!                       EXACTLY 72 CHARACTERS FOLLOWED BY A NEW LINE
!                       SENTINEL FOLLOWED BY MORE CHARACTERS.  NEXTC
!                       SHOULD BE INCREMENTED BY 2.
! 
!       LPIECE .GT. LWRAP+1  REDUCE LPIECE TO LWRAP.
! 
!       ELSE            THIS LAST CASE MEANS 2 .LE. LPIECE .LE. LWRAP+1
!                       RESET LPIECE = LPIECE-1.  NOTE THAT THIS
!                       PROPERLY HANDLES THE END CASE WHERE LPIECE .EQ.
!                       LWRAP+1.  THAT IS, THE SENTINEL FALLS EXACTLY
!                       AT THE END OF A LINE.
! 
      NEXTC = 1
   50 LPIECE = INDEX(MESSG(NEXTC:LENMSG), NEWLIN)
      IF (LPIECE .EQ. 0) THEN
! 
!       THERE WAS NO NEW LINE SENTINEL FOUND.
! 
         IDELTA = 0
         LPIECE = MIN(LWRAP, LENMSG+1-NEXTC)
         IF (LPIECE .LT. LENMSG+1-NEXTC) THEN
            DO 52 I=LPIECE+1,2,-1
               IF (MESSG(NEXTC+I-1:NEXTC+I-1) .EQ. ' ') THEN
                  LPIECE = I-1
                  IDELTA = 1
                  GOTO 54
               ENDIF
   52       CONTINUE
         ENDIF
   54    CBUFF(LPREF+1:LPREF+LPIECE) = MESSG(NEXTC:NEXTC+LPIECE-1)
         NEXTC = NEXTC + LPIECE + IDELTA
      ELSEIF (LPIECE .EQ. 1) THEN
! 
!       WE HAVE A NEW LINE SENTINEL AT MESSG(NEXTC:NEXTC+1).
!       DON'T PRINT A BLANK LINE.
! 
         NEXTC = NEXTC + 2
         GO TO 50
      ELSEIF (LPIECE .GT. LWRAP+1) THEN
! 
!       LPIECE SHOULD BE SET DOWN TO LWRAP.
! 
         IDELTA = 0
         LPIECE = LWRAP
         DO 56 I=LPIECE+1,2,-1
            IF (MESSG(NEXTC+I-1:NEXTC+I-1) .EQ. ' ') THEN
               LPIECE = I-1
               IDELTA = 1
               GOTO 58
            ENDIF
   56    CONTINUE
   58    CBUFF(LPREF+1:LPREF+LPIECE) = MESSG(NEXTC:NEXTC+LPIECE-1)
         NEXTC = NEXTC + LPIECE + IDELTA
      ELSE
! 
!       IF WE ARRIVE HERE, IT MEANS 2 .LE. LPIECE .LE. LWRAP+1.
!       WE SHOULD DECREMENT LPIECE BY ONE.
! 
         LPIECE = LPIECE - 1
         CBUFF(LPREF+1:LPREF+LPIECE) = MESSG(NEXTC:NEXTC+LPIECE-1)
         NEXTC  = NEXTC + LPIECE + 2
      ENDIF
! 
!       PRINT
! 
      DO 60 I=1,NUNIT
         WRITE(IU(I), '(A)') CBUFF(1:LPREF+LPIECE)
   60 CONTINUE
! 
      IF (NEXTC .LE. LENMSG) GO TO 50
      RETURN
      END SUBROUTINE XERPRN
!=============================================================================
!DECK XERSVE
      SUBROUTINE XERSVE (LIBRAR, SUBROU, MESSG, KFLAG, NERR, LEVEL, &
        ICOUNT)
! ***BEGIN PROLOGUE  XERSVE
! ***SUBSIDIARY
! ***PURPOSE  Record that an error has occurred.
! ***LIBRARY   SLATEC (XERROR)
! ***CATEGORY  R3
! ***TYPE      ALL (XERSVE-A)
! ***KEYWORDS  ERROR, XERROR
! ***AUTHOR  Jones, R. E., (SNLA)
! ***DESCRIPTION
! 
! *Usage:
! 
!        INTEGER  KFLAG, NERR, LEVEL, ICOUNT
!        CHARACTER * (len) LIBRAR, SUBROU, MESSG
!  
!        CALL XERSVE (LIBRAR, SUBROU, MESSG, KFLAG, NERR, LEVEL, ICOUNT)
! 
! *Arguments:
! 
!        LIBRAR :IN    is the library that the message is from.
!        SUBROU :IN    is the subroutine that the message is from.
!        MESSG  :IN    is the message to be saved.
!        KFLAG  :IN    indicates the action to be performed.
!                      when KFLAG > 0, the message in MESSG is saved.
!                      when KFLAG=0 the tables will be dumped and
!                      cleared.
!                      when KFLAG < 0, the tables will be dumped and
!                      not cleared.
!        NERR   :IN    is the error number.
!        LEVEL  :IN    is the error severity.
!        ICOUNT :OUT   the number of times this message has been seen,
!                      or zero if the table has overflowed and does not
!                      contain this message specifically.  When KFLAG=0,
!                      ICOUNT will not be altered.
! 
! *Description:
! 
!   Record that this error occurred and possibly dump and clear the
!   tables.
! 
! ***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
!                 Error-handling Package, SAND82-0800, Sandia
!                 Laboratories, 1982.
! ***ROUTINES CALLED  I1MACH, XGETUA
! ***REVISION HISTORY  (YYMMDD)
!   800319  DATE WRITTEN
!   861211  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900413  Routine modified to remove reference to KFLAG.  (WRB)
!   900510  Changed to add LIBRARY NAME and SUBROUTINE to calling
!           sequence, use IF-THEN-ELSE, make number of saved entries
!           easily changeable, changed routine name from XERSAV to
!           XERSVE.  (RWC)
!   910626  Added LIBTAB and SUBTAB to SAVE statement.  (BKS)
!   920501  Reformatted the REFERENCES section.  (WRB)
! ***END PROLOGUE  XERSVE
      INTEGER KFLAG, NERR, LEVEL, ICOUNT
      INTEGER, PARAMETER :: LENTAB=10
      INTEGER LUN(5)
      CHARACTER*(*) LIBRAR, SUBROU, MESSG
      CHARACTER*8  LIBTAB(LENTAB), SUBTAB(LENTAB), LIB, SUB
      CHARACTER*20 MESTAB(LENTAB), MES
      INTEGER NERTAB(LENTAB), LEVTAB(LENTAB), KOUNT(LENTAB)
      INTEGER KOUNTX, NMSG
      SAVE LIBTAB, SUBTAB, MESTAB, NERTAB, LEVTAB, KOUNT, KOUNTX, NMSG
      DATA KOUNTX/0/, NMSG/0/
      INTEGER I,NUNIT,KUNIT
! ***FIRST EXECUTABLE STATEMENT  XERSVE
! 
      IF (KFLAG.LE.0) THEN
! 
!        Dump the table.
! 
         IF (NMSG.EQ.0) RETURN
! 
!        Print to each unit.
! 
         CALL XGETUA (LUN, NUNIT)
         DO 20 KUNIT = 1,NUNIT
            IUNIT = LUN(KUNIT)
            IF (IUNIT.EQ.0) IUNIT = I1MACH(4)
! 
!           Print the table header.
! 
            WRITE (IUNIT,9000)
! 
!           Print body of table.
! 
            DO 10 I = 1,NMSG
               WRITE (IUNIT,9010) LIBTAB(I), SUBTAB(I), MESTAB(I), &
                 NERTAB(I),LEVTAB(I),KOUNT(I)
   10       CONTINUE
! 
!           Print number of other errors.
! 
            IF (KOUNTX.NE.0) WRITE (IUNIT,9020) KOUNTX
            WRITE (IUNIT,9030)
   20    CONTINUE
! 
!        Clear the error tables.
! 
         IF (KFLAG.EQ.0) THEN
            NMSG = 0
            KOUNTX = 0
         ENDIF
      ELSE
! 
!        PROCESS A MESSAGE...
!        SEARCH FOR THIS MESSG, OR ELSE AN EMPTY SLOT FOR THIS MESSG,
!        OR ELSE DETERMINE THAT THE ERROR TABLE IS FULL.
! 
         LIB = LIBRAR
         SUB = SUBROU
         MES = MESSG
         DO 30 I = 1,NMSG
            IF (LIB.EQ.LIBTAB(I) .AND. SUB.EQ.SUBTAB(I) .AND.&
              MES.EQ.MESTAB(I) .AND. NERR.EQ.NERTAB(I) .AND.&
              LEVEL.EQ.LEVTAB(I)) THEN
                  KOUNT(I) = KOUNT(I) + 1
                  ICOUNT = KOUNT(I)
                  RETURN
            ENDIF
   30    CONTINUE
! 
         IF (NMSG.LT.LENTAB) THEN
! 
!           Empty slot found for new message.
! 
            NMSG = NMSG + 1
            LIBTAB(I) = LIB
            SUBTAB(I) = SUB
            MESTAB(I) = MES
            NERTAB(I) = NERR
            LEVTAB(I) = LEVEL
            KOUNT (I) = 1
            ICOUNT    = 1
         ELSE
! 
!           Table is full.
! 
            KOUNTX = KOUNTX+1
            ICOUNT = 0
         ENDIF
      ENDIF
      RETURN
! 
!     Formats.
! 
 9000 FORMAT ('0          ERROR MESSAGE SUMMARY' / &
        ' LIBRARY    SUBROUTINE MESSAGE START             NERR', &
        '     LEVEL     COUNT')
 9010 FORMAT (1X,A,3X,A,3X,A,3I10)
 9020 FORMAT ('0OTHER ERRORS NOT INDIVIDUALLY TABULATED = ', I10)
 9030 FORMAT (1X)
      END SUBROUTINE XERSVE
!=======================================================================
!DECK XERHLT
      SUBROUTINE XERHLT (MESSG)
! ***BEGIN PROLOGUE  XERHLT
! ***SUBSIDIARY
! ***PURPOSE  Abort program execution and print error message.
! ***LIBRARY   SLATEC (XERROR)
! ***CATEGORY  R3C
! ***TYPE      ALL (XERHLT-A)
! ***KEYWORDS  ABORT PROGRAM EXECUTION, ERROR, XERROR
! ***AUTHOR  Jones, R. E., (SNLA)
! ***DESCRIPTION
! 
!     Abstract
!        ***Note*** machine dependent routine
!        XERHLT aborts the execution of the program.
!        The error message causing the abort is given in the calling
!        sequence, in case one needs it for printing on a dayfile,
!        for example.
! 
!     Description of Parameters
!        MESSG is as in XERMSG.
! 
! ***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
!                 Error-handling Package, SAND82-0800, Sandia
!                 Laboratories, 1982.
! ***ROUTINES CALLED  (NONE)
! ***REVISION HISTORY  (YYMMDD)
!   790801  DATE WRITTEN
!   861211  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900206  Routine changed from user-callable to subsidiary.  (WRB)
!   900510  Changed calling sequence to delete length of character
!           and changed routine name from XERABT to XERHLT.  (RWC)
!   920501  Reformatted the REFERENCES section.  (WRB)
! ***END PROLOGUE  XERHLT
      CHARACTER*(*) MESSG
! ***FIRST EXECUTABLE STATEMENT  XERHLT
      STOP
      END SUBROUTINE XERHLT

!=======================================================================
!DECK XGETUA
      SUBROUTINE XGETUA (IUNITA, N)
! ***BEGIN PROLOGUE  XGETUA
! ***PURPOSE  Return unit number(s) to which error messages are being
!            sent.
! ***LIBRARY   SLATEC (XERROR)
! ***CATEGORY  R3C
! ***TYPE      ALL (XGETUA-A)
! ***KEYWORDS  ERROR, XERROR
! ***AUTHOR  Jones, R. E., (SNLA)
! ***DESCRIPTION
! 
!     Abstract
!        XGETUA may be called to determine the unit number or numbers
!        to which error messages are being sent.
!        These unit numbers may have been set by a call to XSETUN,
!        or a call to XSETUA, or may be a default value.
! 
!     Description of Parameters
!      --Output--
!        IUNIT - an array of one to five unit numbers, depending
!                on the value of N.  A value of zero refers to the
!                default unit, as defined by the I1MACH machine
!                constant routine.  Only IUNIT(1),...,IUNIT(N) are
!                defined by XGETUA.  The values of IUNIT(N+1),...,
!                IUNIT(5) are not defined (for N .LT. 5) or altered
!                in any way by XGETUA.
!        N     - the number of units to which copies of the
!                error messages are being sent.  N will be in the
!                range from 1 to 5.
! 
! ***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
!                 Error-handling Package, SAND82-0800, Sandia
!                 Laboratories, 1982.
! ***ROUTINES CALLED  J4SAVE
! ***REVISION HISTORY  (YYMMDD)
!   790801  DATE WRITTEN
!   861211  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
! ***END PROLOGUE  XGETUA
      INTEGER IUNITA(5),N
      INTEGER I,INDEX
! ***FIRST EXECUTABLE STATEMENT  XGETUA
      N = J4SAVE(5,0,.FALSE.)
      DO 30 I=1,N
         INDEX = I+4
         IF (I.EQ.1) INDEX = 3
         IUNITA(I) = J4SAVE(INDEX,0,.FALSE.)
   30 CONTINUE
      RETURN
      END SUBROUTINE XGETUA
!========================================================================
!DECK XERCNT
      SUBROUTINE XERCNT (LIBRAR, SUBROU, MESSG, NERR, LEVEL, KONTRL)
! ***BEGIN PROLOGUE  XERCNT
! ***SUBSIDIARY
! ***PURPOSE  Allow user control over handling of errors.
! ***LIBRARY   SLATE! (XERROR)
! ***CATEGORY  R3C
! ***TYPE      ALL (XERCNT-A)
! ***KEYWORDS  ERROR, XERROR
! ***AUTHOR  Jones, R. E., (SNLA)
! ***DESCRIPTION
! 
!     Abstract
!        Allows user control over handling of individual errors.
!        Just after each message is recorded, but before it is
!        processed any further (i.e., before it is printed or
!        a decision to abort is made), a call is made to XERCNT.
!        If the user has provided his own version of XERCNT, he
!        can then override the value of KONTROL used in processing
!        this message by redefining its value.
!        KONTRL may be set to any value from -2 to 2.
!        The meanings for KONTRL are the same as in XSETF, except
!        that the value of KONTRL changes only for this message.
!        If KONTRL is set to a value outside the range from -2 to 2,
!        it will be moved back into that range.
! 
!     Description of Parameters
! 
!      --Input--
!        LIBRAR - the library that the routine is in.
!        SUBROU - the subroutine that XERMSG is being called from
!        MESSG  - the first 20 characters of the error message.
!        NERR   - same as in the call to XERMSG.
!        LEVEL  - same as in the call to XERMSG.
!        KONTRL - the current value of the control flag as set
!                 by a call to XSETF.
! 
!      --Output--
!        KONTRL - the new value of KONTRL.  If KONTRL is not
!                 defined, it will remain at its original value.
!                 This changed value of control affects only
!                 the current occurrence of the current message.
! 
! ***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
!                 Error-handling Package, SAND82-0800, Sandia
!                 Laboratories, 1982.
! ***ROUTINES CALLED  (NONE)
! ***REVISION HISTORY  (YYMMDD)
!   790801  DATE WRITTEN
!   861211  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900206  Routine changed from user-callable to subsidiary.  (WRB)
!   900510  Changed calling sequence to include LIBRARY and SUBROUTINE
!           names, changed routine name from XERCTL to XERCNT.  (RWC)
!   920501  Reformatted the REFERENCES section.  (WRB)
! ***END PROLOGUE  XERCNT
      INTEGER  NERR, LEVEL, KONTRL
      CHARACTER*(*) LIBRAR, SUBROU, MESSG
! ***FIRST EXECUTABLE STATEMENT  XERCNT
      RETURN
      END SUBROUTINE XERCNT



!========================================================================
!DECK FDUMP
      SUBROUTINE FDUMP
! ***BEGIN PROLOGUE  FDUMP
! ***PURPOSE  Symbolic dump (should be locally written).
! ***LIBRARY   SLATEC (XERROR)
! ***CATEGORY  R3
! ***TYPE      ALL (FDUMP-A)
! ***KEYWORDS  ERROR, XERMSG
! ***AUTHOR  Jones, R. E., (SNLA)
! ***DESCRIPTION
! 
!        ***Note*** Machine Dependent Routine
!        FDUMP is intended to be replaced by a locally written
!        version which produces a symbolic dump.  Failing this,
!        it should be replaced by a version which prints the
!        subprogram nesting list.  Note that this dump must be
!        printed on each of up to five files, as indicated by the
!        XGETUA routine.  See XSETUA and XGETUA for details.
! 
!     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
! 
! ***REFERENCES  (NONE)
! ***ROUTINES CALLED  (NONE)
! ***REVISION HISTORY  (YYMMDD)
!   790801  DATE WRITTEN
!   861211  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
! ***END PROLOGUE  FDUMP
! ***FIRST EXECUTABLE STATEMENT  FDUMP
      RETURN
      END SUBROUTINE FDUMP



!========================================================================
!DECK I1MACH
      INTEGER FUNCTION I1MACH (I)
      IMPLICIT NONE
      INTEGER :: I
      double precision :: X
      DOUBLE PRECISION :: XX
!***BEGIN PROLOGUE  I1MACH
!***PURPOSE  Return integer machine dependent constants.
!***LIBRARY   SLATEC
!***CATEGORY  R1
!***TYPE      INTEGER (I1MACH-I)
!***KEYWORDS  MACHINE CONSTANTS
!***AUTHOR  Fox, P. A., (Bell Labs)
!           Hall, A. D., (Bell Labs)
!           Schryer, N. L., (Bell Labs)
!***DESCRIPTION
!
!   I1MACH can be used to obtain machine-dependent parameters for the
!   local machine environment.  It is a function subprogram with one
!   (input) argument and can be referenced as follows:
!
!        K = I1MACH(I)
!
!   where I=1,...,16.  The (output) value of K above is determined by
!   the (input) value of I.  The results for various values of I are
!   discussed below.
!
!   I/O unit numbers:
!     I1MACH( 1) = the standard input unit.
!     I1MACH( 2) = the standard output unit.
!     I1MACH( 3) = the standard punch unit.
!     I1MACH( 4) = the standard error message unit.
!
!   Words:
!     I1MACH( 5) = the number of bits per integer storage unit.
!     I1MACH( 6) = the number of characters per integer storage unit.
!
!   Integers:
!     assume integers are represented in the S-digit, base-A form
!
!                sign ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) )
!
!                where 0 .LE. X(I) .LT. A for I=0,...,S-1.
!     I1MACH( 7) = A, the base.
!     I1MACH( 8) = S, the number of base-A digits.
!     I1MACH( 9) = A**S - 1, the largest magnitude.
!
!   Floating-Point Numbers:
!     Assume floating-point numbers are represented in the T-digit,
!     base-B form
!                sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
!
!                where 0 .LE. X(I) .LT. B for I=1,...,T,
!                0 .LT. X(1), and EMIN .LE. E .LE. EMAX.
!     I1MACH(10) = B, the base.
!
!   Single-Precision:
!     I1MACH(11) = T, the number of base-B digits.
!     I1MACH(12) = EMIN, the smallest exponent E.
!     I1MACH(13) = EMAX, the largest exponent E.
!
!   Double-Precision:
!     I1MACH(14) = T, the number of base-B digits.
!     I1MACH(15) = EMIN, the smallest exponent E.
!     I1MACH(16) = EMAX, the largest exponent E.
!
!   To alter this function for a particular environment, the desired
!   set of DATA statements should be activated by removing the C from
!   column 1.  Also, the values of I1MACH(1) - I1MACH(4) should be
!   checked for consistency with the local operating system.
!
!***REFERENCES  P. A. Fox, A. D. Hall and N. L. Schryer, Framework for
!                 a portable library, ACM Transactions on Mathematical
!                 Software 4, 2 (June 1978), pp. 177-188.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   750101  DATE WRITTEN
!   960411  Modified for Fortran 90 (BE after suggestions by EHG).   
!   980727  Modified value of I1MACH(6) (BE after suggestion by EHG).   
!***END PROLOGUE  I1MACH
!
      X  = 1.0      
      XX = 1.0D0

      SELECT CASE (I)
        CASE (1)
          I1MACH = 5 ! Input unit
        CASE (2)
          I1MACH = 6 ! Output unit
        CASE (3)
          I1MACH = 0 ! Punch unit is no longer used
        CASE (4)
          I1MACH = 0 ! Error message unit
        CASE (5)
          I1MACH = BIT_SIZE(I)
        CASE (6)
          I1MACH = 4            ! Characters per integer is hopefully no
                                ! longer used. 
                                ! If it is used it has to be set manually.
                                ! The value 4 is correct on IEEE-machines.
        CASE (7)
          I1MACH = RADIX(1)
        CASE (8)
          I1MACH = BIT_SIZE(I) - 1
        CASE (9)
          I1MACH = HUGE(1)
        CASE (10)
          I1MACH = RADIX(X)
        CASE (11)
          I1MACH = DIGITS(X)
        CASE (12)
          I1MACH = MINEXPONENT(X)
        CASE (13)
          I1MACH = MAXEXPONENT(X)
        CASE (14)
          I1MACH = DIGITS(XX)
        CASE (15)
          I1MACH = MINEXPONENT(XX)
        CASE (16)
          I1MACH = MAXEXPONENT(XX) 
        CASE DEFAULT
          WRITE (*, FMT = 9000)
 9000     FORMAT ('1ERROR    1 IN I1MACH - I OUT OF BOUNDS')
          STOP
        END SELECT
      RETURN
      END FUNCTION I1MACH
!=======================================================================
!DECK DSICS
      SUBROUTINE DSICS (N, NELT, IA, JA, A, ISYM, NEL, IEL, JEL, EL, D, R, IWARN)
! ***BEGIN PROLOGUE  DSICS
! ***PURPOSE  Incompl. Cholesky Decomposition Preconditioner SLAP Set Up.
!            Routine to generate the Incomplete Cholesky decomposition,
!            L*D*L-trans, of a symmetric positive definite matrix, A,
!            which is stored in SLAP Column format.  The unit lower
!            triangular matrix L is stored by rows, and the inverse of
!            the diagonal matrix D is stored.
! ***LIBRARY   SLATEC (SLAP)
! ***CATEGORY  D2E
! ***TYPE      DOUBLE PRECISION (SSICS-S, DSICS-D)
! ***KEYWORDS  INCOMPLETE CHOLESKY FACTORIZATION,
!             ITERATIVE PRECONDITION, LINEAR SYSTEM, SLAP SPARSE
! ***AUTHOR  Greenbaum, Anne, (Courant Institute)
!           Seager, Mark K., (LLNL)
!             Lawrence Livermore National Laboratory
!             PO BOX 808, L-60
!             Livermore, CA 94550 (510) 423-3141
!             seager@llnl.gov
! ***DESCRIPTION
! 
! *Usage:
!     INTEGER N, NELT, IA(NELT), JA(NELT), ISYM
!     INTEGER NEL, IEL(NEL), JEL(NEL), IWARN
!     DOUBLE PRECISION A(NELT), EL(NEL), D(N), R(N)
! 
!     CALL DSICS( N, NELT, IA, JA, A, ISYM, NEL, IEL, JEL, EL, D, R,
!    $    IWARN )
! 
! *Arguments:
! N      :IN       Integer.
!         Order of the Matrix.
! NELT   :IN       Integer.
!         Number of elements in arrays IA, JA, and A.
! IA     :INOUT    Integer IA(NELT).
! JA     :INOUT    Integer JA(NELT).
! A      :INOUT    Double Precision A(NELT).
!         These arrays should hold the matrix A in the SLAP Column
!         format.  See "Description", below.
! ISYM   :IN       Integer.
!         Flag to indicate symmetric storage format.
!         If ISYM=0, all non-zero entries of the matrix are stored.
!         If ISYM=1, the matrix is symmetric, and only the lower
!         triangle of the matrix is stored.
! NEL    :OUT      Integer.
!         Number of non-zeros in the lower triangle of A.   Also
!         corresponds to the length of the IEL, JEL, EL arrays.
! IEL    :OUT      Integer IEL(NEL).
! JEL    :OUT      Integer JEL(NEL).
! EL     :OUT      Double Precision EL(NEL).
!         IEL, JEL, EL contain the unit lower triangular factor  of the
!         incomplete decomposition   of the A  matrix  stored  in  SLAP
!         Row format.   The Diagonal of   ones   *IS*   stored.     See
!         "Description", below for more details about the SLAP Row fmt.
! D      :OUT      Double Precision D(N)
!         Upon return this array holds D(I) = 1./DIAG(A).
! R      :WORK     Double Precision R(N).
!         Temporary double precision workspace needed for the
!         factorization.
! IWARN  :OUT      Integer.
!         This is a warning variable and is zero if the IC factoriza-
!         tion goes well.  It is set to the row index corresponding to
!         the last zero pivot found.  See "Description", below.
! 
! *Description
!       =================== S L A P Column format ==================
!       This routine  requires that  the matrix A  be stored in  the
!       SLAP Column format.  In this format the non-zeros are stored
!       counting down columns (except for  the diagonal entry, which
!       must appear first in each  "column")  and are stored  in the
!       double precision array A.   In other words,  for each column
!       in the matrix put the diagonal entry in  A.  Then put in the
!       other non-zero  elements going down  the column (except  the
!       diagonal) in order.   The  IA array holds the  row index for
!       each non-zero.  The JA array holds the offsets  into the IA,
!       A arrays  for  the  beginning  of each   column.   That  is,
!       IA(JA(ICOL)),  A(JA(ICOL)) points   to the beginning  of the
!       ICOL-th   column    in    IA and   A.      IA(JA(ICOL+1)-1),
!       A(JA(ICOL+1)-1) points to  the  end of the   ICOL-th column.
!       Note that we always have  JA(N+1) = NELT+1,  where N is  the
!       number of columns in  the matrix and NELT  is the number  of
!       non-zeros in the matrix.
! 
!       Here is an example of the  SLAP Column  storage format for a
!       5x5 Matrix (in the A and IA arrays '|'  denotes the end of a
!       column):
! 
!           5x5 Matrix      SLAP Column format for 5x5 matrix on left.
!                              1  2  3    4  5    6  7    8    9 10 11
!       |11 12  0  0 15|   A: 11 21 51 | 22 12 | 33 53 | 44 | 55 15 35
!       |21 22  0  0  0|  IA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3
!       | 0  0 33  0 35|  JA:  1  4  6    8  9   12
!       | 0  0  0 44  0|
!       |51  0 53  0 55|
! 
!       ==================== S L A P Row format ====================
! 
!       This routine requires  that the matrix A  be  stored  in the
!       SLAP  Row format.   In this format  the non-zeros are stored
!       counting across  rows (except for the diagonal  entry, which
!       must  appear first  in each  "row")  and  are stored  in the
!       double precision  array A.  In other words, for each row  in
!       the matrix  put the diagonal  entry in A.   Then put in  the
!       other  non-zero elements  going across  the row  (except the
!       diagonal) in order.  The JA array holds the column index for
!       each non-zero.  The IA array holds the offsets  into the JA,
!       A  arrays  for  the   beginning  of  each  row.    That  is,
!       JA(IA(IROW)),A(IA(IROW)) are the first elements of the IROW-
!       th row in  JA and A,  and  JA(IA(IROW+1)-1), A(IA(IROW+1)-1)
!       are  the last elements  of the  IROW-th row.   Note  that we
!       always have  IA(N+1) = NELT+1, where N is the number of rows
!       in the matrix  and  NELT is the  number of non-zeros  in the
!       matrix.
! 
!       Here is an example of the SLAP Row storage format for a  5x5
!       Matrix (in the A and JA arrays '|' denotes the end of a row):
! 
!           5x5 Matrix         SLAP Row format for 5x5 matrix on left.
!                              1  2  3    4  5    6  7    8    9 10 11
!       |11 12  0  0 15|   A: 11 12 15 | 22 21 | 33 35 | 44 | 55 51 53
!       |21 22  0  0  0|  JA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3
!       | 0  0 33  0 35|  IA:  1  4  6    8  9   12
!       | 0  0  0 44  0|
!       |51  0 53  0 55|
! 
!       With the SLAP  format some  of  the   "inner  loops" of this
!       routine should vectorize  on  machines with hardware support
!       for vector   gather/scatter  operations.  Your compiler  may
!       require a compiler directive to  convince it that  there are
!       no  implicit  vector  dependencies.  Compiler directives for
!       the Alliant    FX/Fortran and CRI   CFT/CFT77 compilers  are
!       supplied with the standard SLAP distribution.
! 
!       The IC factorization does not always exist for SPD matrices.
!       In the event that a zero pivot is found it is set  to be 1.0
!       and the factorization proceeds.   The integer variable IWARN
!       is set to the last row where the Diagonal was fudged.  This
!       eventuality hardly ever occurs in practice.
! 
! ***SEE ALSO  DCG, DSICCG
! ***REFERENCES  1. Gene Golub and Charles Van Loan, Matrix Computations,
!                  Johns Hopkins University Press, Baltimore, Maryland,
!                  1983.
! ***ROUTINES CALLED  XERMSG
! ***REVISION HISTORY  (YYMMDD)
!   890404  DATE WRITTEN
!   890404  Previous REVISION DATE
!   890915  Made changes requested at July 1989 CML Meeting.  (MKS)
!   890922  Numerous changes to prologue to make closer to SLATEC
!           standard.  (FNF)
!   890929  Numerous changes to reduce SP/DP differences.  (FNF)
!   900805  Changed XERRWV calls to calls to XERMSG.  (RWC)
!   910411  Prologue converted to Version 4.0 format.  (BAB)
!   920511  Added complete declaration section.  (WRB)
!   920929  Corrected format of reference.  (FNF)
!   930701  Updated CATEGORY section.  (FNF, WRB)
! ***END PROLOGUE  DSICS
!     .. Scalar Arguments ..
      INTEGER ISYM, IWARN, N, NEL, NELT
!     .. Array Arguments ..
      DOUBLE PRECISION A(NELT), D(N), EL(NEL), R(N)
      INTEGER IA(NELT), IEL(NEL), JA(NELT), JEL(NEL)
!     .. Local Scalars ..
      DOUBLE PRECISION ELTMP
      INTEGER I, IBGN, IC, ICBGN, ICEND, ICOL, IEND, IR, IRBGN, IREND, IROW, IRR, J, JBGN, JELTMP, JEND
      CHARACTER XERN1*8
!     .. External Subroutines ..
!      EXTERNAL XERMSG
! ***FIRST EXECUTABLE STATEMENT  DSICS
! 
!         Set the lower triangle in IEL, JEL, EL
! 
      IWARN = 0
! 
!         All matrix elements stored in IA, JA, A.  Pick out the lower
!         triangle (making sure that the Diagonal of EL is one) and
!         store by rows.
!
      NEL = 1
      IEL(1) = 1
      JEL(1) = 1
      EL(1) = 1
      D(1) = A(1)
!CVD$R NOCONCUR

      DO 30 IROW = 2, N
!         Put in the Diagonal.

         NEL = NEL + 1
         IEL(IROW) = NEL
         JEL(NEL) = IROW
         EL(NEL) = 1
         D(IROW) = A(JA(IROW))      
! 
!         Look in all the lower triangle columns for a matching row.
!         Since the matrix is symmetric, we can look across the
!         IROW-th row by looking down the IROW-th column (if it is
!         stored ISYM=0)...
 
         IF( ISYM.EQ.0 ) THEN
            ICBGN = JA(IROW)
            ICEND = JA(IROW+1)-1
         ELSE
            ICBGN = 1
            ICEND = IROW-1
         ENDIF
 

         DO 20 IC = ICBGN, ICEND

            IF( ISYM.EQ.0 ) THEN
               ICOL = IA(IC)
               IF( ICOL.GE.IROW ) GOTO 20
            ELSE
               ICOL = IC
            ENDIF
 
            JBGN = JA(ICOL)+1
            JEND = JA(ICOL+1)-1
 
            IF( JBGN.LE.JEND .AND. IA(JEND).GE.IROW ) THEN

!CVD$ NOVECTOR
 
               DO 10 J = JBGN, JEND

                  IF( IA(J).EQ.IROW ) THEN
                     NEL = NEL + 1
                     JEL(NEL) = ICOL
                     EL(NEL)  = A(J)

                     GOTO 20
                  ENDIF

 10            CONTINUE
            ENDIF

 20      CONTINUE
 30   CONTINUE


      IEL(N+1) = NEL+1
 

       
! 
!         Sort ROWS of lower triangle into descending order (count out
!         along rows out from Diagonal).
! 
      DO 60 IROW = 2, N
         IBGN = IEL(IROW)+1

         IEND = IEL(IROW+1)-1

         IF( IBGN.LT.IEND ) THEN

            DO 50 I = IBGN, IEND-1

!CVD$ NOVECTOR


               DO 40 J = I+1, IEND
                   IF( JEL(I) .GT. JEL(J)) THEN
 
                     JELTMP = JEL(J)
                     JEL(J) = JEL(I)
                     JEL(I) = JELTMP
                     ELTMP = EL(J)
                     EL(J) = EL(I)
                     EL(I) = ELTMP
                  ENDIF

 40            CONTINUE
 50         CONTINUE
         ENDIF
 60   CONTINUE
! 
!         Perform the Incomplete Cholesky decomposition by looping
!         over the rows.
!         Scale the first column.  Use the structure of A to pick out
!         the rows with something in column 1.
! 

      IRBGN = JA(1)+1
      IREND = JA(2)-1

      DO 65 IRR = IRBGN, IREND
         IR = IA(IRR)
!         Find the index into EL for EL(1,IR).
!         Hint: it's the second entry.
         I = IEL(IR)+1
         EL(I) = EL(I)/D(1)
  65   CONTINUE
!



      DO 110 IROW = 2, N
! 
!         Update the IROW-th diagonal.
! 
         DO 66 I = 1, IROW-1
            R(I) = 0
 66      CONTINUE

         IBGN = IEL(IROW)+1
         IEND = IEL(IROW+1)-1
          IF( IBGN.LE.IEND ) THEN

!CLLL. OPTION ASSERT (NOHAZARD)
!CDIR$ IVDEP
!CVD$ NODEPCHK

            DO 70 I = IBGN, IEND
                 R(JEL(I)) = EL(I)*D(JEL(I))
                 D(IROW) = D(IROW) - EL(I)*R(JEL(I))
              ! write(*,*)'EL(I)=',EL(I), 'D(JEL(I))',D(JEL(I)),'D(IROW)=',D(IROW),'R(JEL(I))=',R(JEL(I))

  70         CONTINUE


! 
!         Check to see if we have a problem with the diagonal.
!
            IF( D(IROW).LE.0.0D0 ) THEN
               IF( IWARN.EQ.0 ) IWARN = IROW
               D(IROW) = 1
            ENDIF
         ENDIF
! 
!         Update each EL(IROW+1:N,IROW), if there are any.
!         Use the structure of A to determine the Non-zero elements
!         of the IROW-th column of EL.
! 


         IRBGN = JA(IROW)
         IREND = JA(IROW+1)-1
         DO 100 IRR = IRBGN, IREND
            IR = IA(IRR)
            IF( IR.LE.IROW ) GOTO 100

!         Find the index into EL for EL(IR,IROW)

            IBGN = IEL(IR)+1
            IEND = IEL(IR+1)-1
            IF( JEL(IBGN).GT.IROW ) GOTO 100


            DO 90 I = IBGN, IEND
               IF( JEL(I).EQ.IROW ) THEN
                  ICEND = IEND
 91               IF( JEL(ICEND).GE.IROW ) THEN
                     ICEND = ICEND - 1
                     GOTO 91
                  ENDIF

!Sum up the EL(IR,1:IROW-1)*R(1:IROW-1) contributions.
!CLLL. OPTION ASSERT (NOHAZARD)
!CDIR$ IVDEP
!CVD$ NODEPCHK

                  DO 80 IC = IBGN, ICEND
                     EL(I) = EL(I) - EL(IC)*R(JEL(IC))
 80               CONTINUE
                  EL(I) = EL(I)/D(IROW)                
                  GOTO 100
               ENDIF
 90         CONTINUE
! 
!         If we get here, we have real problems...
            WRITE (XERN1, '(I8)') IROW
            CALL XERMSG ('SLATEC', 'DSICS',&
              'A and EL data structure mismatch in row '// XERN1, 1, 2)
 100     CONTINUE
 110  CONTINUE

 
! 
!         Replace diagonals by their inverses.
! 
!CVD$ CONCUR
      DO 120 I =1, N
         D(I) = 1.0D0/D(I)
   120  CONTINUE
     
      RETURN
! ------------- LAST LINE OF DSICS FOLLOWS ----------------------------
      END  SUBROUTINE DSICS 
! ----------------------------------------------------------------------------------
!*DECK DS2Y
      SUBROUTINE DS2Y (N, NELT, IA, JA, A, ISYM)
!***BEGIN PROLOGUE  DS2Y
!***PURPOSE  SLAP Triad to SLAP Column Format Converter.
!            Routine to convert from the SLAP Triad to SLAP Column
!            format.
!***LIBRARY   SLATEC (SLAP)
!***CATEGORY  D1B9
!***TYPE      DOUBLE PRECISION (SS2Y-S, DS2Y-D)
!***KEYWORDS  LINEAR SYSTEM, SLAP SPARSE
!***AUTHOR  Seager, Mark K., (LLNL)
!              Lawrence Livermore National Laboratory
!             PO BOX 808, L-60
!             Livermore, CA 94550 (510) 423-3141
!             seager@llnl.gov
!***DESCRIPTION
! 
! *Usage:
!     INTEGER N, NELT, IA(NELT), JA(NELT), ISYM
!     DOUBLE PRECISION A(NELT)
! 
!     CALL DS2Y( N, NELT, IA, JA, A, ISYM )
! 
! *Arguments:
! N      :IN       Integer
!         Order of the Matrix.
! NELT   :IN       Integer.
!         Number of non-zeros stored in A.
! IA     :INOUT    Integer IA(NELT).
! JA     :INOUT    Integer JA(NELT).
! A      :INOUT    Double Precision A(NELT).
!         These arrays should hold the matrix A in either the SLAP
!         Triad format or the SLAP Column format.  See "Description",
!         below.  If the SLAP Triad format is used, this format is
!         translated to the SLAP Column format by this routine.
! ISYM   :IN       Integer.
!         Flag to indicate symmetric storage format.
!         If ISYM=0, all non-zero entries of the matrix are stored.
!         If ISYM=1, the matrix is symmetric, and only the lower
!         triangle of the matrix is stored.
! 
! *Description:
!       The Sparse Linear Algebra Package (SLAP) utilizes two matrix
!       data structures: 1) the  SLAP Triad  format or  2)  the SLAP
!       Column format.  The user can hand this routine either of the
!       of these data structures.  If the SLAP Triad format is give
!       as input then this routine transforms it into SLAP Column
!       format.  The way this routine tells which format is given as
!       input is to look at JA(N+1).  If JA(N+1) = NELT+1 then we
!       have the SLAP Column format.  If that equality does not hold
!       then it is assumed that the IA, JA, A arrays contain the
!       SLAP Triad format.
!
!       =================== S L A P Triad format ===================
!       This routine requires that the  matrix A be   stored in  the
!       SLAP  Triad format.  In  this format only the non-zeros  are
!       stored.  They may appear in  *ANY* order.  The user supplies
!       three arrays of  length NELT, where  NELT is  the number  of
!       non-zeros in the matrix: (IA(NELT), JA(NELT), A(NELT)).  For
!       each non-zero the user puts the row and column index of that
!       matrix element  in the IA and  JA arrays.  The  value of the
!       non-zero   matrix  element is  placed  in  the corresponding
!       location of the A array.   This is  an  extremely  easy data
!       structure to generate.  On  the  other hand it   is  not too
!       efficient on vector computers for  the iterative solution of
!       linear systems.  Hence,   SLAP changes   this  input    data
!       structure to the SLAP Column format  for  the iteration (but
!       does not change it back).
! 
!       Here is an example of the  SLAP Triad   storage format for a
!       5x5 Matrix.  Recall that the entries may appear in any order.
! 
!           5x5 Matrix      SLAP Triad format for 5x5 matrix on left.
!                              1  2  3  4  5  6  7  8  9 10 11
!       |11 12  0  0 15|   A: 51 12 11 33 15 53 55 22 35 44 21
!       |21 22  0  0  0|  IA:  5  1  1  3  1  5  5  2  3  4  2
!       | 0  0 33  0 35|  JA:  1  2  1  3  5  3  5  2  5  4  1
!       | 0  0  0 44  0|
!       |51  0 53  0 55|
! 
!       =================== S L A P Column format ==================
! 
!       This routine  requires that  the matrix A  be stored in  the
!       SLAP Column format.  In this format the non-zeros are stored
!       counting down columns (except for  the diagonal entry, which
!       must appear first in each  "column")  and are stored  in the
!       real array A.  In other words, for each column in the matrix
!       put the diagonal entry in A.  Then put in the other non-zero
!       elements going down   the  column (except  the diagonal)  in
!       order.  The IA array holds the row  index for each non-zero.
!       The JA array holds the offsets into the IA, A arrays for the
!       beginning of   each    column.    That  is,    IA(JA(ICOL)),
!       A(JA(ICOL)) points to the beginning of the ICOL-th column in
!       IA and  A.  IA(JA(ICOL+1)-1),  A(JA(ICOL+1)-1) points to the
!       end  of   the ICOL-th  column.  Note   that  we  always have
!       JA(N+1) = NELT+1, where  N  is the number of columns in  the
!       matrix and  NELT   is the number of non-zeros in the matrix.
! 
!       Here is an example of the  SLAP Column  storage format for a
!       5x5 Matrix (in the A and IA arrays '|'  denotes the end of a
!       column):
! 
!           5x5 Matrix      SLAP Column format for 5x5 matrix on left.
!                              1  2  3    4  5    6  7    8    9 10 11
!       |11 12  0  0 15|   A: 11 21 51 | 22 12 | 33 53 | 44 | 55 15 35
!       |21 22  0  0  0|  IA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3
!       | 0  0 33  0 35|  JA:  1  4  6    8  9   12
!       | 0  0  0 44  0|
!       |51  0 53  0 55|
! 
!***REFERENCES  (NONE)
!***ROUTINES CALLED  QS2I1D
!***REVISION HISTORY  (YYMMDD)
!   871119  DATE WRITTEN
!   881213  Previous REVISION DATE
!   890915  Made changes requested at July 1989 CML Meeting.  (MKS)
!   890922  Numerous changes to prologue to make closer to SLATEC
!           standard.  (FNF)
!   890929  Numerous changes to reduce SP/DP differences.  (FNF)
!   910411  Prologue converted to Version 4.0 format.  (BAB)
!   910502  Corrected C***FIRST EXECUTABLE STATEMENT line.  (FNF)
!   920511  Added complete declaration section.  (WRB)
!   930701  Updated CATEGORY section.  (FNF, WRB)
!***END PROLOGUE  DS2Y
!     .. Scalar Arguments ..
      INTEGER ISYM, N, NELT
!     .. Array Arguments ..
      DOUBLE PRECISION A(NELT)
      INTEGER IA(NELT), JA(NELT)
!     .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER I, IBGN, ICOL, IEND, ITEMP, J
!     .. External Subroutines ..
!      EXTERNAL QS2I1D
!***FIRST EXECUTABLE STATEMENT  DS2Y
! 
!         Check to see if the (IA,JA,A) arrays are in SLAP Column
!         format.  If it's not then transform from SLAP Triad.
! 
      IF( JA(N+1).EQ.NELT+1 ) RETURN
! 
!         Sort into ascending order by COLUMN (on the ja array).
!         This will line up the columns.
! 
      CALL QS2I1D( JA, IA, A, NELT, 1 )
! 
!         Loop over each column to see where the column indices change
!         in the column index array ja.  This marks the beginning of the
!         next column.
! 
! VD$R NOVECTOR
      JA(1) = 1
      DO 20 ICOL = 1, N-1
         DO 10 J = JA(ICOL)+1, NELT
            IF( JA(J).NE.ICOL ) THEN
               JA(ICOL+1) = J
               GOTO 20
            ENDIF
 10      CONTINUE
 20   CONTINUE
      JA(N+1) = NELT+1
! 
!         Mark the n+2 element so that future calls to a SLAP routine
!         utilizing the YSMP-Column storage format will be able to tell.
! 
          !JA(N+2) = 0

                 DO I=N+2,NELT
                      JA(I)=0
                 ENDDO


!         Now loop through the IA array making sure that the diagonal
!         matrix element appears first in the column.  Then sort the
!         rest of the column in ascending order.
! 
      DO 70 ICOL = 1, N
         IBGN = JA(ICOL)
         IEND = JA(ICOL+1)-1
         DO 30 I = IBGN, IEND
            IF( IA(I).EQ.ICOL ) THEN
! 
!              Swap the diagonal element with the first element in the
!              column.
! 
               ITEMP = IA(I)
               IA(I) = IA(IBGN)
               IA(IBGN) = ITEMP
               TEMP = A(I)
               A(I) = A(IBGN)
               A(IBGN) = TEMP
               GOTO 40
            ENDIF
 30      CONTINUE
 40      IBGN = IBGN + 1
         IF( IBGN.LT.IEND ) THEN
            DO 60 I = IBGN, IEND
               DO 50 J = I+1, IEND
                  IF( IA(I).GT.IA(J) ) THEN
                     ITEMP = IA(I)
                     IA(I) = IA(J)
                     IA(J) = ITEMP
                     TEMP = A(I)
                     A(I) = A(J)
                     A(J) = TEMP
                  ENDIF
 50            CONTINUE
 60         CONTINUE
         ENDIF
 70   CONTINUE
      RETURN
! ------------- LAST LINE OF DS2Y FOLLOWS ----------------------------
      END SUBROUTINE DS2Y
! -------------------------------------------------------------------------------------------------
!*DECK QS2I1D
      SUBROUTINE QS2I1D (IA, JA, A, N, KFLAG)
!***BEGIN PROLOGUE  QS2I1D
!***SUBSIDIARY
!***PURPOSE  Sort an integer array, moving an integer and DP array.
!            This routine sorts the integer array IA and makes the same
!            interchanges in the integer array JA and the double pre-
!            cision array A.  The array IA may be sorted in increasing
!            order or decreasing order.  A slightly modified QUICKSORT
!            algorithm is used.
!***LIBRARY   SLATEC (SLAP)
!***CATEGORY  N6A2A
!***TYPE      DOUBLE PRECISION (QS2I1R-S, QS2I1D-D)
!***KEYWORDS  SINGLETON QUICKSORT, SLAP, SORT, SORTING
!***AUTHOR  Jones, R. E., (SNLA)
!           Kahaner, D. K., (NBS)
!           Seager, M. K., (LLNL) seager@llnl.gov
!           Wisniewski, J. A., (SNLA)
!***DESCRIPTION
!     Written by Rondall E Jones
!     Modified by John A. Wisniewski to use the Singleton QUICKSORT
!     algorithm. date 18 November 1976.
! 
!     Further modified by David K. Kahaner
!     National Bureau of Standards
!     August, 1981
! 
!     Even further modification made to bring the code up to the
!     Fortran 77 level and make it more readable and to carry
!     along one integer array and one double precision array during
!     the sort by
!     Mark K. Seager
!     Lawrence Livermore National Laboratory
!     November, 1987
!     This routine was adapted from the ISORT routine.
! 
!     ABSTRACT
!         This routine sorts an integer array IA and makes the same
!         interchanges in the integer array JA and the double precision
!         array A.
!         The array IA may be sorted in increasing order or decreasing
!         order.  A slightly modified quicksort algorithm is used.
! 
!     DESCRIPTION OF PARAMETERS
!        IA - Integer array of values to be sorted.
!        JA - Integer array to be carried along.
!         A - Double Precision array to be carried along.
!         N - Number of values in integer array IA to be sorted.
!     KFLAG - Control parameter
!           = 1 means sort IA in INCREASING order.
!           =-1 means sort IA in DECREASING order.
! 
!***SEE ALSO  DS2Y
!***REFERENCES  R. C. Singleton, Algorithm 347, An Efficient Algorithm
!                 for Sorting With Minimal Storage, Communications ACM
!                 12:3 (1969), pp.185-7.
!***ROUTINES CALLED  XERMSG
!***REVISION HISTORY  (YYMMDD)
!   761118  DATE WRITTEN
!   890125  Previous REVISION DATE
!   890915  Made changes requested at July 1989 CML Meeting.  (MKS)
!   890922  Numerous changes to prologue to make closer to SLATEC
!           standard.  (FNF)
!   890929  Numerous changes to reduce SP/DP differences.  (FNF)
!   900805  Changed XERROR calls to calls to XERMSG.  (RWC)
!   910411  Prologue converted to Version 4.0 format.  (BAB)
!   910506  Made subsidiary to DS2Y and corrected reference.  (FNF)
!   920511  Added complete declaration section.  (WRB)
!   920929  Corrected format of reference.  (FNF)
!   921012  Corrected all f.p. constants to double precision.  (FNF)
!***END PROLOGUE  QS2I1D
! VD$R NOVECTOR
! VD$R NOCONCUR
!     .. Scalar Arguments ..
      INTEGER KFLAG, N
!     .. Array Arguments ..
      DOUBLE PRECISION A(N)
      INTEGER IA(N), JA(N)
!     .. Local Scalars ..
      DOUBLE PRECISION R, TA, TTA
      INTEGER I, IIT, IJ, IT, J, JJT, JT, K, KK, L, M, NN
!     .. Local Arrays ..
      INTEGER IL(21), IU(21)
!     .. External Subroutines ..
!      EXTERNAL XERMSG
!     .. Intrinsic Functions ..
      INTRINSIC ABS, INT
!***FIRST EXECUTABLE STATEMENT  QS2I1D
      NN = N
      IF (NN.LT.1) THEN
         CALL XERMSG ('SLATEC', 'QS2I1D',&
           'The number of values to be sorted was not positive.', 1, 1)
         RETURN
      ENDIF
      IF( N.EQ.1 ) RETURN
      KK = ABS(KFLAG)
      IF ( KK.NE.1 ) THEN
         CALL XERMSG ('SLATEC', 'QS2I1D',&
           'The sort control parameter, K, was not 1 or -1.', 2, 1)
         RETURN
      ENDIF
! 
!     Alter array IA to get decreasing order if needed.
! 
      IF( KFLAG.LT.1 ) THEN
         DO 20 I=1,NN
            IA(I) = -IA(I)
 20      CONTINUE
      ENDIF
! 
!     Sort IA and carry JA and A along.
!     And now...Just a little black magic...
      M = 1
      I = 1
      J = NN
      R = .375D0
 210  IF( R.LE.0.5898437D0 ) THEN
         R = R + 3.90625D-2
      ELSE
         R = R-.21875D0
      ENDIF
 225  K = I
! 
!     Select a central element of the array and save it in location
!     it, jt, at.
! 
      IJ = I + INT ((J-I)*R)
      IT = IA(IJ)
      JT = JA(IJ)
      TA = A(IJ)
! 
!     If first element of array is greater than it, interchange with it.
! 
      IF( IA(I).GT.IT ) THEN
         IA(IJ) = IA(I)
         IA(I)  = IT
         IT     = IA(IJ)
         JA(IJ) = JA(I)
         JA(I)  = JT
         JT     = JA(IJ)
         A(IJ)  = A(I)
         A(I)   = TA
         TA     = A(IJ)
      ENDIF
      L=J
! 
!     If last element of array is less than it, swap with it.
! 
      IF( IA(J).LT.IT ) THEN
         IA(IJ) = IA(J)
         IA(J)  = IT
         IT     = IA(IJ)
         JA(IJ) = JA(J)
         JA(J)  = JT
         JT     = JA(IJ)
         A(IJ)  = A(J)
         A(J)   = TA
         TA     = A(IJ)
! 
!     If first element of array is greater than it, swap with it.
! 
         IF ( IA(I).GT.IT ) THEN
            IA(IJ) = IA(I)
            IA(I)  = IT
            IT     = IA(IJ)
            JA(IJ) = JA(I)
            JA(I)  = JT
            JT     = JA(IJ)
            A(IJ)  = A(I)
            A(I)   = TA
            TA     = A(IJ)
         ENDIF
      ENDIF
! 
!     Find an element in the second half of the array which is
!     smaller than it.
! 
  240 L=L-1
      IF( IA(L).GT.IT ) GO TO 240
! 
!     Find an element in the first half of the array which is
!     greater than it.
! 
  245 K=K+1
      IF( IA(K).LT.IT ) GO TO 245
! 
!     Interchange these elements.
! 
      IF( K.LE.L ) THEN
         IIT   = IA(L)
         IA(L) = IA(K)
         IA(K) = IIT
         JJT   = JA(L)
         JA(L) = JA(K)
         JA(K) = JJT
         TTA   = A(L)
         A(L)  = A(K)
         A(K)  = TTA
         GOTO 240
      ENDIF
! 
!     Save upper and lower subscripts of the array yet to be sorted.
! 
      IF( L-I.GT.J-K ) THEN
         IL(M) = I
         IU(M) = L
         I = K
         M = M+1
      ELSE
         IL(M) = K
         IU(M) = J
         J = L
         M = M+1
      ENDIF
      GO TO 260
! 
!     Begin again on another portion of the unsorted array.
! 
  255 M = M-1
      IF( M.EQ.0 ) GO TO 300
      I = IL(M)
      J = IU(M)
  260 IF( J-I.GE.1 ) GO TO 225
      IF( I.EQ.J ) GO TO 255
      IF( I.EQ.1 ) GO TO 210
      I = I-1
  265 I = I+1
      IF( I.EQ.J ) GO TO 255
      IT = IA(I+1)
      JT = JA(I+1)
      TA =  A(I+1)
      IF( IA(I).LE.IT ) GO TO 265
      K=I
  270 IA(K+1) = IA(K)
      JA(K+1) = JA(K)
      A(K+1)  =  A(K)
      K = K-1
      IF( IT.LT.IA(K) ) GO TO 270
      IA(K+1) = IT
      JA(K+1) = JT
      A(K+1)  = TA
      GO TO 265
! 
!     Clean up, if necessary.
! 
  300 IF( KFLAG.LT.1 ) THEN
         DO 310 I=1,NN
            IA(I) = -IA(I)
 310     CONTINUE
      ENDIF
      RETURN
! ------------- LAST LINE OF QS2I1D FOLLOWS ----------------------------
      END SUBROUTINE QS2I1D 

! -------------------------------------------------------------------------------------------------
!DECK DCG
      SUBROUTINE DCG (N, B, X, NELT, IA, JA, A, ISYM, MATVEC, MSOLVE,&
        ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, R, Z, P, DZ, RWORK,&
        IWORK)
!***BEGIN PROLOGUE  DCG
!***PURPOSE  Preconditioned Conjugate Gradient Sparse Ax=b Solver.
!            Routine to solve a symmetric positive definite linear
!            system  Ax = b  using the Preconditioned Conjugate
!            Gradient method.
!***LIBRARY   SLATEC (SLAP)
!***CATEGORY  D2B4
!***TYPE      DOUBLE PRECISION (SCG-S, DCG-D)
!***KEYWORDS  ITERATIVE PRECONDITION, SLAP, SPARSE,
!             SYMMETRIC LINEAR SYSTEM
!***AUTHOR  Greenbaum, Anne, (Courant Institute)
!           Seager, Mark K., (LLNL)
!             Lawrence Livermore National Laboratory
!             PO BOX 808, L-60
!             Livermore, CA 94550 (510) 423-3141
!             seager@llnl.gov
!***DESCRIPTION
! 
! *Usage:
!     INTEGER  N, NELT, IA(NELT), JA(NELT), ISYM, ITOL, ITMAX
!     INTEGER  ITER, IERR, IUNIT, IWORK(USER DEFINED)
!     DOUBLE PRECISION B(N), X(N), A(NELT), TOL, ERR, R(N), Z(N)
!     DOUBLE PRECISION P(N), DZ(N), RWORK(USER DEFINED)
!     EXTERNAL MATVEC, MSOLVE
! 
!     CALL DCG(N, B, X, NELT, IA, JA, A, ISYM, MATVEC, MSOLVE,
!    $     ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, R, Z, P, DZ,
!    $     RWORK, IWORK )
! 
! *Arguments:
! N      :IN       Integer.
!         Order of the Matrix.
! B      :IN       Double Precision B(N).
!         Right-hand side vector.
! X      :INOUT    Double Precision X(N).
!         On input X is your initial guess for solution vector.
!         On output X is the final approximate solution.
! NELT   :IN       Integer.
!         Number of Non-Zeros stored in A.
! IA     :IN       Integer IA(NELT).
! JA     :IN       Integer JA(NELT).
! A      :IN       Double Precision A(NELT).
!         These arrays contain the matrix data structure for A.
!         It could take any form.  See "Description", below,
!         for more details.
! ISYM   :IN       Integer.
!         Flag to indicate symmetric storage format.
!         If ISYM=0, all non-zero entries of the matrix are stored.
!         If ISYM=1, the matrix is symmetric, and only the upper
!         or lower triangle of the matrix is stored.
! MATVEC :EXT      External.
!         Name of a routine which performs the matrix vector multiply
!         Y = A*X given A and X.  The name of the MATVEC routine must
!         be declared external in the calling program.  The calling
!         sequence to MATVEC is:
! 
! 
!         Where N is the number of unknowns, Y is the product A*X
!         upon return X is an input vector, NELT is the number of
!         non-zeros in the SLAP IA, JA, A storage for the matrix A.
!         ISYM is a flag which, if non-zero, denotest that A is
!         symmetric and only the lower or upper triangle is stored.
! MSOLVE :EXT      External.
!         Name of a routine which solves a linear system MZ = R for
!         Z given R with the preconditioning matrix M (M is supplied via
!         RWORK and IWORK arrays).  The name of the MSOLVE routine must
!         be declared external in the calling program.  The calling
!         sequence to MSOLVE is:
! 
!             CALL MSOLVE(N, R, Z, NELT, IA, JA, A, ISYM, RWORK, IWORK)
! 
!         Where N is the number of unknowns, R is the right-hand side
!         vector and Z is the solution upon return.  NELT, IA, JA, A and
!         ISYM are defined as above.  RWORK is a double precision array
!         that can be used to pass necessary preconditioning information
!         and/or workspace to MSOLVE.  IWORK is an integer work array
!         for the same purpose as RWORK.
! ITOL   :IN       Integer.
!         Flag to indicate type of convergence criterion.
!         If ITOL=1, iteration stops when the 2-norm of the residual
!         divided by the 2-norm of the right-hand side is less than TOL.
!         If ITOL=2, iteration stops when the 2-norm of M-inv times the
!         residual divided by the 2-norm of M-inv times the right hand
!         side is less than TOL, where M-inv is the inverse of the
!         diagonal of A.
!         ITOL=11 is often useful for checking and comparing different
!         routines.  For this case, the user must supply the "exact"
!         solution or a very accurate approximation (one with an error
!         much less than TOL) through a common block,
!             COMMON /DSLBLK/ SOLN( )
!         If ITOL=11, iteration stops when the 2-norm of the difference
!         between the iterative approximation and the user-supplied
!         solution divided by the 2-norm of the user-supplied solution
!         is less than TOL.  Note that this requires the user to set up
!         the "COMMON /DSLBLK/ SOLN(LENGTH)" in the calling routine.
!         The routine with this declaration should be loaded before the
!         stop test so that the correct length is used by the loader.
!         This procedure is not standard Fortran and may not work
!         correctly on your system (although it has worked on every
!         system the authors have tried).  If ITOL is not 11 then this
!         common block is indeed standard Fortran.
! TOL    :INOUT    Double Precision.
!         Convergence criterion, as described above.  (Reset if IERR=4.)
! ITMAX  :IN       Integer.
!         Maximum number of iterations.
! ITER   :OUT      Integer.
!         Number of iterations required to reach convergence, or
!         ITMAX+1 if convergence criterion could not be achieved in
!         ITMAX iterations.
! ERR    :OUT      Double Precision.
!         Error estimate of error in final approximate solution, as
!         defined by ITOL.
! IERR   :OUT      Integer.
!         Return error flag.
!           IERR = 0 => All went well.
!           IERR = 1 => Insufficient space allocated for WORK or IWORK.
!           IERR = 2 => Method failed to converge in ITMAX steps.
!           IERR = 3 => Error in user input.
!                       Check input values of N, ITOL.
!           IERR = 4 => User error tolerance set too tight.
!                       Reset to 500*D1MACH(3).  Iteration proceeded.
!           IERR = 5 => Preconditioning matrix, M, is not positive
!                       definite.  (r,z) < 0.
!           IERR = 6 => Matrix A is not positive definite.  (p,Ap) < 0.
! IUNIT  :IN       Integer.
!         Unit number on which to write the error at each iteration,
!         if this is desired for monitoring convergence.  If unit
!         number is 0, no writing will occur.
! R      :WORK     Double Precision R(N).
! Z      :WORK     Double Precision Z(N).
! P      :WORK     Double Precision P(N).
! DZ     :WORK     Double Precision DZ(N).
!         Double Precision arrays used for workspace.
! RWORK  :WORK     Double Precision RWORK(USER DEFINED).
!         Double Precision array that can be used by  MSOLVE.
! IWORK  :WORK     Integer IWORK(USER DEFINED).
!         Integer array that can be used by  MSOLVE.
! 
! *Description
!       This routine does  not care  what matrix data   structure is
!       used for  A and M.  It simply   calls  the MATVEC and MSOLVE
!       routines, with  the arguments as  described above.  The user
!       could write any type of structure and the appropriate MATVEC
!       and MSOLVE routines.  It is assumed  that A is stored in the
!       IA, JA, A  arrays in some fashion and  that M (or INV(M)) is
!       stored  in  IWORK  and  RWORK   in  some fashion.   The SLAP
!       routines DSDCG and DSICCG are examples of this procedure.
! 
!       Two  examples  of  matrix  data structures  are the: 1) SLAP
!       Triad  format and 2) SLAP Column format.
! 
! *Cautions:
!     This routine will attempt to write to the Fortran logical output
!     unit IUNIT, if IUNIT .ne. 0.  Thus, the user must make sure that
!     this logical unit is attached to a file or terminal before calling
!     this routine with a non-zero value for IUNIT.  This routine does
!     not check for the validity of a non-zero IUNIT unit number.
! 
!***SEE ALSO  DSDCG, DSICCG
!***REFERENCES  1. Louis Hageman and David Young, Applied Iterative
!                  Methods, Academic Press, New York, 1981.
!               2. Concus, Golub and O'Leary, A Generalized Conjugate
!                  Gradient Method for the Numerical Solution of
!                  Elliptic Partial Differential Equations, in Sparse
!                  Matrix Computations, Bunch and Rose, Eds., Academic
!                  Press, New York, 1979.
!               3. Mark K. Seager, A SLAP for the Masses, in
!                  G. F. Carey, Ed., Parallel Supercomputing: Methods,
!                  Algorithms and Applications, Wiley, 1989, pp.135-155.
!***ROUTINES CALLED  D1MACH, DAXPY, DCOPY, DDOT, ISDCG
!***REVISION HISTORY  (YYMMDD)
!   890404  DATE WRITTEN
!   890404  Previous REVISION DATE
!   890915  Made changes requested at July 1989 CML Meeting.  (MKS)
!   890921  Removed TeX from comments.  (FNF)
!   890922  Numerous changes to prologue to make closer to SLATEC
!           standard.  (FNF)
!   890929  Numerous changes to reduce SP/DP differences.  (FNF)
!   891004  Added new reference.
!   910411  Prologue converted to Version 4.0 format.  (BAB)
!   910502  Removed MATVEC and MSOLVE from ROUTINES CALLED list.  (FNF)
!   920407  COMMON BLOCK renamed DSLBLK.  (WRB)
!   920511  Added complete declaration section.  (WRB)
!   920929  Corrected format of references.  (FNF)
!   921019  Changed 500.0 to 500 to reduce SP/DP differences.  (FNF)
!***END PROLOGUE  DCG
!     .. Scalar Arguments ..
      DOUBLE PRECISION ERR, TOL
      INTEGER IERR, ISYM, ITER, ITMAX, ITOL, IUNIT, N, NELT
!     .. Array Arguments ..
      DOUBLE PRECISION A(NELT), B(N), DZ(N), P(N), R(N), RWORK(*), X(N),Z(N)
      INTEGER IA(NELT), IWORK(*), JA(NELT)
!     .. Subroutine Arguments ..
!      EXTERNAL MATVEC, MSOLVE

  INTERFACE 
          SUBROUTINE MSOLVE(N, R, Z, NELT, IA, JA, A, ISYM, RWORK, IWORK)
      INTEGER N, ISYM, NELT
       DOUBLE PRECISION A(NELT) , RWORK(*) ,Z(N), R(N)  
      INTEGER IA(NELT), IWORK(*), JA(NELT)

          END SUBROUTINE MSOLVE
  END INTERFACE

  INTERFACE 
          SUBROUTINE MATVEC(N, X, R, NELT, IA, JA, A, ISYM)
      INTEGER N, ISYM, NELT
      DOUBLE PRECISION A(NELT), R(N), X(N)
      INTEGER IA(NELT),  JA(NELT)
     
          END SUBROUTINE MATVEC
  END INTERFACE


!     .. Local Scalars ..
      DOUBLE PRECISION AK, AKDEN, BK, BKDEN, BKNUM, BNRM, SOLNRM, TOLMIN
      INTEGER I, K
!     .. External Functions ..
!      DOUBLE PRECISION D1MACH, DDOT
!      INTEGER ISDCG
!      EXTERNAL D1MACH, DDOT, ISDCG
!     .. External Subroutines ..
!      EXTERNAL DAXPY, DCOPY
!***FIRST EXECUTABLE STATEMENT  DCG
! 
!         Check some of the input data.
! 
      ITER = 0
      IERR = 0
      IF( N.LT.1 ) THEN
         IERR = 3
         RETURN
      ENDIF
      TOLMIN = 500*D1MACH(3)
      IF( TOL.LT.TOLMIN ) THEN
         TOL = TOLMIN
         IERR = 4
      ENDIF
! 
!         Calculate initial residual and pseudo-residual, and check
!         stopping criterion.
      CALL MATVEC(N, X, R, NELT, IA, JA, A, ISYM)
      DO 10 I = 1, N
         R(I) = B(I) - R(I)
 10   CONTINUE
      CALL MSOLVE(N, R, Z, NELT, IA, JA, A, ISYM, RWORK, IWORK)
! 
      IF( ISDCG(N, B, X, NELT, IA, JA, A, ISYM, MSOLVE, ITOL, TOL,&
          ITMAX, ITER, ERR, IERR, IUNIT, R, Z, P, DZ,&
          RWORK, IWORK, AK, BK, BNRM, SOLNRM) .NE. 0 ) GO TO 200
      IF( IERR.NE.0 ) RETURN
! 
!         ***** Iteration loop *****
! 
      DO 100 K=1,ITMAX
         ITER = K
! 
!         Calculate coefficient bk and direction vector p.
         BKNUM = DDOT(N, Z, 1, R, 1)
         IF( BKNUM.LE.0.0D0 ) THEN
            IERR = 5
            RETURN
         ENDIF
         IF(ITER .EQ. 1) THEN
            CALL DCOPY(N, Z, 1, P, 1)
         ELSE
            BK = BKNUM/BKDEN
            DO 20 I = 1, N
               P(I) = Z(I) + BK*P(I)
 20         CONTINUE
         ENDIF
         BKDEN = BKNUM
! 
!         Calculate coefficient ak, new iterate x, new residual r,
!         and new pseudo-residual z.
         CALL MATVEC(N, P, Z, NELT, IA, JA, A, ISYM)
         AKDEN = DDOT(N, P, 1, Z, 1)
         IF( AKDEN.LE.0.0D0 ) THEN
            IERR = 6
            RETURN
         ENDIF
         AK = BKNUM/AKDEN
         CALL DAXPY(N, AK, P, 1, X, 1)
         CALL DAXPY(N, -AK, Z, 1, R, 1)
         CALL MSOLVE(N, R, Z, NELT, IA, JA, A, ISYM, RWORK, IWORK)
! 
!         check stopping criterion.
         IF( ISDCG(N, B, X, NELT, IA, JA, A, ISYM, MSOLVE, ITOL, TOL,&
             ITMAX, ITER, ERR, IERR, IUNIT, R, Z, P, DZ, RWORK,&
             IWORK, AK, BK, BNRM, SOLNRM) .NE. 0 ) GO TO 200
! 
 100  CONTINUE
! 
!         *****   end of loop  *****
! 
!         stopping criterion not satisfied.
      ITER = ITMAX + 1
      IERR = 2
! 
 200  RETURN
! ------------- LAST LINE OF DCG FOLLOWS -----------------------------
      END SUBROUTINE DCG 
!-----------------------------------------------------------------------
!DECK DSLLTI----->MSOLVE
      SUBROUTINE  DSLLTI(N, B, X, NELT, IA, JA, A, ISYM, RWORK, IWORK)
!  ***BEGIN PROLOGUE  DSLLTI
! ***PURPOSE  SLAP MSOLVE for LDL' (IC) Factorization.
!             This routine acts as an interface between the SLAP generic
!             MSOLVE calling convention and the routine that actually
!                            -1
!            computes (LDL')  B = X.
!***LIBRARY   SLATEC (SLAP)
!***CATEGORY  D2E
!***TYPE      DOUBLE PRECISION (SSLLTI-S, DSLLTI-D)
!***KEYWORDS  ITERATIVE PRECONDITION, LINEAR SYSTEM SOLVE, SLAP, SPARSE
!***AUTHOR  Greenbaum, Anne, (Courant Institute)
!            Seager, Mark K., (LLNL)
!             Lawrence Livermore National Laboratory
!             PO BOX 808, L-60
!              Livermore, CA 94550 (510) 423-3141
!              seager@llnl.gov
!***DESCRIPTION
!       It is assumed that RWORK and IWORK have initialized with
!       the information required for DLLTI2:
!          IWORK(1) = NEL
!          IWORK(2) = Starting location of IEL in IWORK.
!          IWORK(3) = Starting location of JEL in IWORK.
!          IWORK(4) = Starting location of EL in RWORK.
!          IWORK(5) = Starting location of DINV in RWORK.
!       See the DESCRIPTION of DLLTI2 for details.
!***REFERENCES  (NONE)
!***ROUTINES CALLED  DLLTI2
!***REVISION HISTORY  (YYMMDD)
!   871119  DATE WRITTEN
!   881213  Previous REVISION DATE
!   890915  Made changes requested at July 1989 CML Meeting.  (MKS)
!   890922  Numerous changes to prologue to make closer to SLATEC
!           standard.  (FNF)
!   890929  Numerous changes to reduce SP/DP differences.  (FNF)
!   910411  Prologue converted to Version 4.0 format.  (BAB)
!   910502  Corrected conversion error.  (FNF)
!   920511  Added complete declaration section.  (WRB)
!   921113  Corrected C***CATEGORY line.  (FNF)
!   930701  Updated CATEGORY section.  (FNF, WRB)
!***END PROLOGUE  DSLLTI
!     .. Scalar Arguments ..
      INTEGER ISYM, N, NELT
!     .. Array Arguments ..
      DOUBLE PRECISION A(NELT), B(N), RWORK(*), X(N)
      INTEGER IA(NELT), IWORK(*), JA(NELT)
!     .. Local Scalars ..
      INTEGER LOCDIN, LOCEL, LOCIEL, LOCJEL, NEL
!     .. External Subroutines ..
!     EXTERNAL DLLTI2
! ***FIRST EXECUTABLE STATEMENT  DSLLTI
      NEL = IWORK(1)
      LOCIEL = IWORK(3)
      LOCJEL = IWORK(2)
      LOCEL  = IWORK(4)
      LOCDIN = IWORK(5)
      CALL DLLTI2(N, B, X, NEL, IWORK(LOCIEL),IWORK(LOCJEL),RWORK(LOCEL), RWORK(LOCDIN))
! 
      RETURN
! ------------- LAST LINE OF DSLLTI FOLLOWS ----------------------------
      END SUBROUTINE DSLLTI
!---------------------------------------------------------------------
!DECK DCOPY
      SUBROUTINE DCOPY (N, DX, INCX, DY, INCY)
!***BEGIN PROLOGUE  DCOPY
!***PURPOSE  Copy a vector.
!***LIBRARY   SLATEC (BLAS)
!***CATEGORY  D1A5
!***TYPE      DOUBLE PRECISION (SCOPY-S, DCOPY-D, CCOPY-C, ICOPY-I)
!***KEYWORDS  BLAS, COPY, LINEAR ALGEBRA, VECTOR
!***AUTHOR  Lawson, C. L., (JPL)
!           Hanson, R. J., (SNLA)
!           Kincaid, D. R., (U. of Texas)
!           Krogh, F. T., (JPL)
!***DESCRIPTION
! 
!                B L A S  Subprogram
!    Description of Parameters
! 
!     --Input--
!        N  number of elements in input vector(s)
!       DX  double precision vector with N elements
!     INCX  storage spacing between elements of DX
!       DY  double precision vector with N elements
!     INCY  storage spacing between elements of DY
! 
!     --Output--
!       DY  copy of vector DX (unchanged if N .LE. 0)
! 
!     Copy double precision DX to double precision DY.
!     For I = 0 to N-1, copy DX(LX+I*INCX) to DY(LY+I*INCY),
!     where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is
!     defined in a similar way using INCY.
! 
!***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
!                 Krogh, Basic linear algebra subprograms for Fortran
!                 usage, Algorithm No. 539, Transactions on Mathematical
!                 Software 5, 3 (September 1979), pp. 308-323.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   791001  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920310  Corrected definition of LX in DESCRIPTION.  (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DCOPY
      INTEGER N, INCX, INCY
      DOUBLE PRECISION DX(*), DY(*)
      INTEGER IX, IY, I, M, MP1, NS
!***FIRST EXECUTABLE STATEMENT  DCOPY
      IF (N .LE. 0) RETURN
      IF (INCX .EQ. INCY) IF (INCX-1) 5,20,60
! 
!     Code for unequal or nonpositive increments.
! 
    5 IX = 1
      IY = 1
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
      IF (INCY .LT. 0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DY(IY) = DX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
! 
!     Code for both increments equal to 1.
! 
!     Clean-up loop so remaining vector length is a multiple of 7.
! 
   20 M = MOD(N,7)
      IF (M .EQ. 0) GO TO 40
      DO 30 I = 1,M
        DY(I) = DX(I)
   30 CONTINUE
      IF (N .LT. 7) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,7
        DY(I) = DX(I)
        DY(I+1) = DX(I+1)
        DY(I+2) = DX(I+2)
        DY(I+3) = DX(I+3)
        DY(I+4) = DX(I+4)
        DY(I+5) = DX(I+5)
        DY(I+6) = DX(I+6)
   50 CONTINUE
      RETURN
! 
!     Code for equal, positive, non-unit increments.
! 
   60 NS = N*INCX
      DO 70 I = 1,NS,INCX
        DY(I) = DX(I)
   70 CONTINUE
      RETURN
      END SUBROUTINE DCOPY
!-------------------------------------------------------------------
!DECK DDOT
       FUNCTION DDOT (N, DX, INCX, DY, INCY)
!***BEGIN PROLOGUE  DDOT
!***PURPOSE  Compute the inner product of two vectors.
!***LIBRARY   SLATEC (BLAS)
!***CATEGORY  D1A4
!***TYPE      DOUBLE PRECISION (SDOT-S, DDOT-D, CDOTU-C)
!***KEYWORDS  BLAS, INNER PRODUCT, LINEAR ALGEBRA, VECTOR
!***AUTHOR  Lawson, C. L., (JPL)
!           Hanson, R. J., (SNLA)
!           Kincaid, D. R., (U. of Texas)
!           Krogh, F. T., (JPL)
!***DESCRIPTION
! 
!                B L A S  Subprogram
!    Description of Parameters
! 
!     --Input--
!        N  number of elements in input vector(s)
!       DX  double precision vector with N elements
!     INCX  storage spacing between elements of DX
!       DY  double precision vector with N elements
!     INCY  storage spacing between elements of DY
! 
!     --Output--
!     DDOT  double precision dot product (zero if N .LE. 0)
! 
!     Returns the dot product of double precision DX and DY.
!     DDOT = sum for I = 0 to N-1 of  DX(LX+I*INCX) * DY(LY+I*INCY),
!     where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is
!     defined in a similar way using INCY.
! 
!***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
!                 Krogh, Basic linear algebra subprograms for Fortran
!                 usage, Algorithm No. 539, Transactions on Mathematical
!                 Software 5, 3 (September 1979), pp. 308-323.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   791001  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920310  Corrected definition of LX in DESCRIPTION.  (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DDOT
      INTEGER N, INCX, INCY
       double precision DDOT
      DOUBLE PRECISION DX(*), DY(*)
      INTEGER IX, IY, I, M, MP1, NS

!***FIRST EXECUTABLE STATEMENT  DDOT
      DDOT = 0.0D0
      IF (N .LE. 0) RETURN
      IF (INCX .EQ. INCY) IF (INCX-1) 5,20,60
! 
!     Code for unequal or nonpositive increments.
! 
    5 IX = 1
      IY = 1
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
      IF (INCY .LT. 0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DDOT = DDOT + DX(IX)*DY(IY)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
! 
!     Code for both increments equal to 1.
! 
!     Clean-up loop so remaining vector length is a multiple of 5.
! 
   20 M = MOD(N,5)
      IF (M .EQ. 0) GO TO 40
      DO 30 I = 1,M
         DDOT = DDOT + DX(I)*DY(I)
   30 CONTINUE
      IF (N .LT. 5) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
      DDOT = DDOT + DX(I)*DY(I) + DX(I+1)*DY(I+1) + DX(I+2)*DY(I+2) +&
                   DX(I+3)*DY(I+3) + DX(I+4)*DY(I+4)
   50 CONTINUE
      RETURN
! 
!     Code for equal, positive, non-unit increments.
! 
   60 NS = N*INCX
      DO 70 I = 1,NS,INCX
        DDOT = DDOT + DX(I)*DY(I)
   70 CONTINUE
      RETURN
      END FUNCTION DDOT
!----------------------------------------------------------------------
!DECK DLLTI2
      SUBROUTINE DLLTI2 (N, B, X, NEL, IEL, JEL, EL, DINV)
!***BEGIN PROLOGUE  DLLTI2
!***PURPOSE  SLAP Backsolve routine for LDL' Factorization.
!            Routine to solve a system of the form  L*D*L' X = B,
!            where L is a unit lower triangular matrix and D is a
!            diagonal matrix and ' means transpose.
!***LIBRARY   SLATEC (SLAP)
!***CATEGORY  D2E
!***TYPE      DOUBLE PRECISION (SLLTI2-S, DLLTI2-D)
!***KEYWORDS  INCOMPLETE FACTORIZATION, ITERATIVE PRECONDITION, SLAP,
!             SPARSE, SYMMETRIC LINEAR SYSTEM SOLVE
!***AUTHOR  Greenbaum, Anne, (Courant Institute)
!           Seager, Mark K., (LLNL)
!             Lawrence Livermore National Laboratory
!             PO BOX 808, L-60
!             Livermore, CA 94550 (510) 423-3141
!             seager@llnl.gov
! ***DESCRIPTION
! 
! *Usage:
!     INTEGER N, NEL, IEL(NEL), JEL(NEL)
!     DOUBLE PRECISION B(N), X(N), EL(NEL), DINV(N)
!
!     CALL DLLTI2( N, B, X, NEL, IEL, JEL, EL, DINV )
! 
! *Arguments:
! N      :IN       Integer
!         Order of the Matrix.
! B      :IN       Double Precision B(N).
!         Right hand side vector.
! X      :OUT      Double Precision X(N).
!         Solution to L*D*L' x = b.
! NEL    :IN       Integer.
!         Number of non-zeros in the EL array.
! IEL    :IN       Integer IEL(NEL).
! JEL    :IN       Integer JEL(NEL).
! EL     :IN       Double Precision     EL(NEL).
!         IEL, JEL, EL contain the unit lower triangular factor   of
!         the incomplete decomposition   of the A  matrix  stored in
!         SLAP Row format.   The diagonal of ones *IS* stored.  This
!         structure can be set  up  by  the DS2LT routine.  See  the
!         "Description", below for more details about the  SLAP  Row
!         format.
! DINV   :IN       Double Precision DINV(N).
!         Inverse of the diagonal matrix D.
! 
! *Description:
!       This routine is supplied with  the SLAP package as a routine
!       to perform the MSOLVE operation in the SCG iteration routine
!       for  the driver  routine DSICCG.   It must be called via the
!       SLAP  MSOLVE calling sequence  convention  interface routine
!       DSLLI.
!         **** THIS ROUTINE ITSELF DOES NOT CONFORM TO THE ****
!               **** SLAP MSOLVE CALLING CONVENTION ****
! 
!       IEL, JEL, EL should contain the unit lower triangular factor
!       of  the incomplete decomposition of  the A matrix  stored in
!       SLAP Row format.   This IC factorization  can be computed by
!       the  DSICS routine.  The  diagonal  (which is all one's) is
!       stored.
! 
!       With  the SLAP  Row format  the "inner loop" of this routine
!       should vectorize   on machines with   hardware  support  for
!       vector gather/scatter operations.  Your compiler may require
!       a  compiler directive  to  convince   it that there  are  no
!       implicit vector  dependencies.  Compiler directives  for the
!       Alliant FX/Fortran and CRI CFT/CFT77 compilers  are supplied
!       with the standard SLAP distribution.
! 
!***SEE ALSO  DSICCG, DSICS
!***REFERENCES  (NONE)
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   871119  DATE WRITTEN
!   881213  Previous REVISION DATE
!   890915  Made changes requested at July 1989 CML Meeting.  (MKS)
!   890922  Numerous changes to prologue to make closer to SLATEC
!           standard.  (FNF)
!   890929  Numerous changes to reduce SP/DP differences.  (FNF)
!   910411  Prologue converted to Version 4.0 format.  (BAB)
!   920511  Added complete declaration section.  (WRB)
!   921113  Corrected C***CATEGORY line.  (FNF)
!   930701  Updated CATEGORY section.  (FNF, WRB)
! ***END PROLOGUE  DLLTI2
!     .. Scalar Arguments ..
      INTEGER N, NEL
!     .. Array Arguments ..
      DOUBLE PRECISION B(N), DINV(N), EL(NEL), X(N)
      INTEGER IEL(NEL), JEL(NEL)
!     .. Local Scalars ..
      INTEGER I, IBGN, IEND, IROW
! ***FIRST EXECUTABLE STATEMENT  DLLTI2
! 
!         Solve  L*y = b,  storing result in x.
! 
      DO 10 I=1,N
         X(I) = B(I)
 10   CONTINUE
      DO 30 IROW = 1, N
         IBGN = IEL(IROW) + 1
         IEND = IEL(IROW+1) - 1
         IF( IBGN.LE.IEND ) THEN
! LLL. OPTION ASSERT (NOHAZARD)
! DIR$ IVDEP
! VD$ NOCONCUR
! VD$ NODEPCHK
            DO 20 I = IBGN, IEND
               X(IROW) = X(IROW) - EL(I)*X(JEL(I))
 20         CONTINUE
         ENDIF
 30   CONTINUE
! 
!         Solve  D*Z = Y,  storing result in X.
! 
      DO 40 I=1,N
         X(I) = X(I)*DINV(I)
 40   CONTINUE
! 
!         Solve  L-trans*X = Z.
! 
      DO 60 IROW = N, 2, -1
         IBGN = IEL(IROW) + 1
         IEND = IEL(IROW+1) - 1
         IF( IBGN.LE.IEND ) THEN
! LLL. OPTION ASSERT (NOHAZARD)
! DIR$ IVDEP
! VD$ NOCONCUR
! VD$ NODEPCHK
            DO 50 I = IBGN, IEND
               X(JEL(I)) = X(JEL(I)) - EL(I)*X(IROW)
 50         CONTINUE
         ENDIF
 60   CONTINUE
! 
      RETURN
! ------------- LAST LINE OF DLLTI2 FOLLOWS ----------------------------
      END SUBROUTINE DLLTI2
!---------------------------------------------------------------------
!DECK DNRM2
       FUNCTION DNRM2 (N, DX, INCX)
!***BEGIN PROLOGUE  DNRM2
!***PURPOSE  Compute the Euclidean length (L2 norm) of a vector.
!***LIBRARY   SLATEC (BLAS)
!***CATEGORY  D1A3B
!***TYPE      DOUBLE PRECISION (SNRM2-S, DNRM2-D, SCNRM2-C)
!***KEYWORDS  BLAS, EUCLIDEAN LENGTH, EUCLIDEAN NORM, L2,
!             LINEAR ALGEBRA, UNITARY, VECTOR
!***AUTHOR  Lawson, C. L., (JPL)
!           Hanson, R. J., (SNLA)
!           Kincaid, D. R., (U. of Texas)
!           Krogh, F. T., (JPL)
!***DESCRIPTION
! 
!                B L A S  Subprogram
!    Description of parameters
! 
!     --Input--
!        N  number of elements in input vector(s)
!       DX  double precision vector with N elements
!     INCX  storage spacing between elements of DX
! 
!     --Output--
!    DNRM2  double precision result (zero if N .LE. 0)
! 
!     Euclidean norm of the N-vector stored in DX with storage
!     increment INCX.
!     If N .LE. 0, return with result = 0.
!     If N .GE. 1, then INCX must be .GE. 1
! 
!     Four phase method using two built-in constants that are
!     hopefully applicable to all machines.
!         CUTLO = maximum of  SQRT(U/EPS)  over all known machines.
!         CUTHI = minimum of  SQRT(V)      over all known machines.
!     where
!         EPS = smallest no. such that EPS + 1. .GT. 1.
!         U   = smallest positive no.   (underflow limit)
!         V   = largest  no.            (overflow  limit)
! 
!     Brief outline of algorithm.
! 
!     Phase 1 scans zero components.
!     move to phase 2 when a component is nonzero and .LE. CUTLO
!     move to phase 3 when a component is .GT. CUTLO
!     move to phase 4 when a component is .GE. CUTHI/M
!     where M = N for X() real and M = 2*N for complex.
! 
!     Values for CUTLO and CUTHI.
!     From the environmental parameters listed in the IMSL converter
!     document the limiting values are as follows:
!     CUTLO, S.P.   U/EPS = 2**(-102) for  Honeywell.  Close seconds are
!                   Univac and DEC at 2**(-103)
!                   Thus CUTLO = 2**(-51) = 4.44089E-16
!     CUTHI, S.P.   V = 2**127 for Univac, Honeywell, and DEC.
!                   Thus CUTHI = 2**(63.5) = 1.30438E19
!     CUTLO, D.P.   U/EPS = 2**(-67) for Honeywell and DEC.
!                   Thus CUTLO = 2**(-33.5) = 8.23181D-11
!     CUTHI, D.P.   same as S.P.  CUTHI = 1.30438D19
!     DATA CUTLO, CUTHI /8.232D-11,  1.304D19/
!     DATA CUTLO, CUTHI /4.441E-16,  1.304E19/
! 
!***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
!                 Krogh, Basic linear algebra subprograms for Fortran
!                 usage, Algorithm No. 539, Transactions on Mathematical
!                 Software 5, 3 (September 1979), pp. 308-323.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   791001  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DNRM2
      INTEGER N, INCX
      INTEGER NEXT
      DOUBLE PRECISION DX(*), CUTLO, CUTHI, HITEST, SUM, XMAX, ZERO,ONE
      double precision     DNRM2
      SAVE CUTLO, CUTHI, ZERO, ONE
      DATA ZERO, ONE /0.0D0, 1.0D0/
! 
      DATA CUTLO, CUTHI /8.232D-11,  1.304D19/
      INTEGER NN, I, J
!***FIRST EXECUTABLE STATEMENT  DNRM2
      IF (N .GT. 0) GO TO 10
         DNRM2  = ZERO
         GO TO 300
! 
   10 ASSIGN 30 TO NEXT
      SUM = ZERO
      NN = N * INCX
! 
!                                                 BEGIN MAIN LOOP
! 
      I = 1
   20    GO TO NEXT,(30, 50, 70, 110)
   30 IF (ABS(DX(I)) .GT. CUTLO) GO TO 85
      ASSIGN 50 TO NEXT
      XMAX = ZERO
! 
!                        PHASE 1.  SUM IS ZERO
! 
   50 IF (DX(I) .EQ. ZERO) GO TO 200
      IF (ABS(DX(I)) .GT. CUTLO) GO TO 85
! 
!                                PREPARE FOR PHASE 2.
! 
      ASSIGN 70 TO NEXT
      GO TO 105
! 
!                                PREPARE FOR PHASE 4.
! 
  100 I = J
      ASSIGN 110 TO NEXT
      SUM = (SUM / DX(I)) / DX(I)
  105 XMAX = ABS(DX(I))
      GO TO 115
! 
!                   PHASE 2.  SUM IS SMALL.
!                             SCALE TO AVOID DESTRUCTIVE UNDERFLOW.
! 
   70 IF (ABS(DX(I)) .GT. CUTLO) GO TO 75
! 
!                     COMMON CODE FOR PHASES 2 AND 4.
!                     IN PHASE 4 SUM IS LARGE.  SCALE TO AVOID OVERFLOW.
! 
  110 IF (ABS(DX(I)) .LE. XMAX) GO TO 115
         SUM = ONE + SUM * (XMAX / DX(I))**2
         XMAX = ABS(DX(I))
         GO TO 200
! 
  115 SUM = SUM + (DX(I)/XMAX)**2
      GO TO 200
! 
!                  PREPARE FOR PHASE 3.
! 
   75 SUM = (SUM * XMAX) * XMAX
! 
!     FOR REAL OR D.P. SET HITEST = CUTHI/N
!     FOR COMPLEX      SET HITEST = CUTHI/(2*N)
! 
   85 HITEST = CUTHI / N
! 
!                   PHASE 3.  SUM IS MID-RANGE.  NO SCALING.
! 
      DO 95 J = I,NN,INCX
      IF (ABS(DX(J)) .GE. HITEST) GO TO 100
   95    SUM = SUM + DX(J)**2
      DNRM2 = SQRT(SUM)
      GO TO 300
! 
  200 CONTINUE
      I = I + INCX
      IF (I .LE. NN) GO TO 20
! 
!              END OF MAIN LOOP.
! 
!              COMPUTE SQUARE ROOT AND ADJUST FOR SCALING.
! 
      DNRM2 = XMAX * SQRT(SUM)
  300 CONTINUE
      RETURN
      END FUNCTION DNRM2
!--------------------------------------------------------------------
!DECK DSMV------>MATVEC
      SUBROUTINE DSMV (N, X, Y, NELT, IA, JA, A, ISYM)
!***BEGIN PROLOGUE  DSMV
!***PURPOSE  SLAP Column Format Sparse Matrix Vector Product.
!            Routine to calculate the sparse matrix vector product:
!            Y = A*X.
!***LIBRARY   SLATEC (SLAP)
!***CATEGORY  D1B4
!***TYPE      DOUBLE PRECISION (SSMV-S, DSMV-D)
!***KEYWORDS  MATRIX VECTOR MULTIPLY, SLAP, SPARSE
!***AUTHOR  Greenbaum, Anne, (Courant Institute)
!           Seager, Mark K., (LLNL)
!             Lawrence Livermore National Laboratory
!             PO BOX 808, L-60
!             Livermore, CA 94550 (510) 423-3141
!             seager@llnl.gov
!***DESCRIPTION
! 
! *Usage:
!     INTEGER  N, NELT, IA(NELT), JA(NELT), ISYM
!     DOUBLE PRECISION X(N), Y(N), A(NELT)
! 
!     CALL DSMV(N, X, Y, NELT, IA, JA, A, ISYM )
! 
! *Arguments:
! N      :IN       Integer.
!         Order of the Matrix.
! X      :IN       Double Precision X(N).
!         The vector that should be multiplied by the matrix.
! Y      :OUT      Double Precision Y(N).
!         The product of the matrix and the vector.
! NELT   :IN       Integer.
!         Number of Non-Zeros stored in A.
! IA     :IN       Integer IA(NELT).
! JA     :IN       Integer JA(NELT).
! A      :IN       Double Precision A(NELT).
!         These arrays should hold the matrix A in the SLAP Column
!         format.  See "Description", below.
! ISYM   :IN       Integer.
!         Flag to indicate symmetric storage format.
!         If ISYM=0, all non-zero entries of the matrix are stored.
!         If ISYM=1, the matrix is symmetric, and only the upper
!         or lower triangle of the matrix is stored.
! 
! *Description
! 
!       With  the SLAP  format  the "inner  loops" of  this  routine
!       should vectorize   on machines with   hardware  support  for
!       vector gather/scatter operations.  Your compiler may require
!       a  compiler directive  to  convince   it that there  are  no
!       implicit vector  dependencies.  Compiler directives  for the
!       Alliant FX/Fortran and CRI CFT/CFT77 compilers  are supplied
!       with the standard SLAP distribution.
! 
! *Cautions:
!     This   routine   assumes  that  the matrix A is stored in SLAP
!     Column format.  It does not check  for  this (for  speed)  and
!     evil, ugly, ornery and nasty things  will happen if the matrix
!     data  structure  is,  in fact, not SLAP Column.  Beware of the
!     wrong data structure!!!
! 
!***SEE ALSO  DSMTV
!***REFERENCES  (NONE)
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   871119  DATE WRITTEN
!   881213  Previous REVISION DATE
!   890915  Made changes requested at July 1989 CML Meeting.  (MKS)
!   890922  Numerous changes to prologue to make closer to SLATEC
!           standard.  (FNF)
!   890929  Numerous changes to reduce SP/DP differences.  (FNF)
!   910411  Prologue converted to Version 4.0 format.  (BAB)
!   920511  Added complete declaration section.  (WRB)
!   930701  Updated CATEGORY section.  (FNF, WRB)
!***END PROLOGUE  DSMV
!     .. Scalar Arguments ..
      INTEGER ISYM, N, NELT
!     .. Array Arguments ..
      DOUBLE PRECISION A(NELT), X(N), Y(N)
      INTEGER IA(NELT), JA(NELT)
!     .. Local Scalars ..
      INTEGER I, IBGN, ICOL, IEND, IROW, J, JBGN, JEND
!***FIRST EXECUTABLE STATEMENT  DSMV
! 
!         Zero out the result vector.
! 
      DO 10 I = 1, N
         Y(I) = 0
 10   CONTINUE
! 
!         Multiply by A.
! 
! VD$R NOCONCUR
      DO 30 ICOL = 1, N
         IBGN = JA(ICOL)
         IEND = JA(ICOL+1)-1
! LLL. OPTION ASSERT (NOHAZARD)
! DIR$ IVDEP
! VD$ NODEPCHK
         DO 20 I = IBGN, IEND
            Y(IA(I)) = Y(IA(I)) + A(I)*X(ICOL)
 20      CONTINUE
 30   CONTINUE
! 
      IF( ISYM.EQ.1 ) THEN
! 
!         The matrix is non-symmetric.  Need to get the other half in...
!         This loops assumes that the diagonal is the first entry in
!         each column.
! 
         DO 50 IROW = 1, N
            JBGN = JA(IROW)+1
            JEND = JA(IROW+1)-1
            IF( JBGN.GT.JEND ) GOTO 50
            DO 40 J = JBGN, JEND
               Y(IROW) = Y(IROW) + A(J)*X(IA(J))
 40         CONTINUE
 50      CONTINUE
      ENDIF
      RETURN
! ------------- LAST LINE OF DSMV FOLLOWS ----------------------------
      END SUBROUTINE DSMV
!-----------------------------------------------------------------------
!DECK ISDCG
       FUNCTION ISDCG (N, B, X, NELT, IA, JA, A, ISYM, MSOLVE,&
        ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, R, Z, P, DZ, RWORK,&
        IWORK, AK, BK, BNRM, SOLNRM)
!***BEGIN PROLOGUE  ISDCG
!***SUBSIDIARY
!***PURPOSE  Preconditioned Conjugate Gradient Stop Test.
!            This routine calculates the stop test for the Conjugate
!            Gradient iteration scheme.  It returns a non-zero if the
!            error estimate (the type of which is determined by ITOL)
!            is less than the user specified tolerance TOL.
!***LIBRARY   SLATEC (SLAP)
!***CATEGORY  D2B4
!***TYPE      DOUBLE PRECISION (ISSCG-S, ISDCG-D)
!***KEYWORDS  LINEAR SYSTEM, SLAP, SPARSE, STOP TEST
!***AUTHOR  Greenbaum, Anne, (Courant Institute)
!           Seager, Mark K., (LLNL)
!             Lawrence Livermore National Laboratory
!             PO BOX 808, L-60
!             Livermore, CA 94550 (510) 423-3141
!             seager@llnl.gov
!***DESCRIPTION
! 
! *Usage:
!     INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, ITOL, ITMAX, ITER
!     INTEGER IERR, IUNIT, IWORK(USER DEFINED)
!     DOUBLE PRECISION B(N), X(N), A(N), TOL, ERR, R(N), Z(N)
!     DOUBLE PRECISION P(N), DZ(N), RWORK(USER DEFINED), AK, BK
!     DOUBLE PRECISION BNRM, SOLNRM
!     EXTERNAL MSOLVE
! 
!     IF( ISDCG(N, B, X, NELT, IA, JA, A, ISYM, MSOLVE, ITOL, TOL,
!    $     ITMAX, ITER, ERR, IERR, IUNIT, R, Z, P, DZ, RWORK, IWORK,
!    $     AK, BK, BNRM, SOLNRM) .NE. 0 ) THEN ITERATION DONE
! 
! *Arguments:
! N      :IN       Integer.
!         Order of the Matrix.
! B      :IN       Double Precision B(N).
!         Right-hand side vector.
! X      :IN       Double Precision X(N).
!         The current approximate solution vector.
! NELT   :IN       Integer.
!         Number of Non-Zeros stored in A.
! IA     :IN       Integer IA(NELT).
! JA     :IN       Integer JA(NELT).
! A      :IN       Double Precision A(NELT).
!         These arrays should hold the matrix A in either the SLAP
!         Triad format or the SLAP Column format.  See "Description"
!         in the DCG, DSDCG or DSICCG routines.
! ISYM   :IN       Integer.
!         Flag to indicate symmetric storage format.
!         If ISYM=0, all non-zero entries of the matrix are stored.
!         If ISYM=1, the matrix is symmetric, and only the upper
!         or lower triangle of the matrix is stored.
! MSOLVE :EXT      External.
!         Name of a routine which solves a linear system MZ = R for
!         Z given R with the preconditioning matrix M (M is supplied via
!         RWORK and IWORK arrays).  The name of the MSOLVE routine must
!         be declared external in the calling program.  The calling
!         sequence to MSOLVE is:
!             CALL MSOLVE(N, R, Z, NELT, IA, JA, A, ISYM, RWORK, IWORK)
!         Where N is the number of unknowns, R is the right-hand side
!         vector and Z is the solution upon return.  NELT, IA, JA, A and
!         ISYM are defined as above.  RWORK is a double precision array
!         that can be used to pass necessary preconditioning information
!         and/or workspace to MSOLVE.  IWORK is an integer work array
!         for the same purpose as RWORK.
! ITOL   :IN       Integer.
!         Flag to indicate type of convergence criterion.
!         If ITOL=1, iteration stops when the 2-norm of the residual
!         divided by the 2-norm of the right-hand side is less than TOL.
!         If ITOL=2, iteration stops when the 2-norm of M-inv times the
!         residual divided by the 2-norm of M-inv times the right hand
!         side is less than TOL, where M-inv is the inverse of the
!         diagonal of A.
!         ITOL=11 is often useful for checking and comparing different
!         routines.  For this case, the user must supply the "exact"
!         solution or a very accurate approximation (one with an error
!         much less than TOL) through a common block,
!             COMMON /DSLBLK/ SOLN( )
!         If ITOL=11, iteration stops when the 2-norm of the difference
!         between the iterative approximation and the user-supplied
!         solution divided by the 2-norm of the user-supplied solution
!         is less than TOL.  Note that this requires the user to set up
!         the "COMMON /DSLBLK/ SOLN(LENGTH)" in the calling routine.
!         The routine with this declaration should be loaded before the
!         stop test so that the correct length is used by the loader.
!         This procedure is not standard Fortran and may not work
!         correctly on your system (although it has worked on every
!         system the authors have tried).  If ITOL is not 11 then this
!         common block is indeed standard Fortran.
! TOL    :IN       Double Precision.
!         Convergence criterion, as described above.
! ITMAX  :IN       Integer.
!         Maximum number of iterations.
! ITER   :IN       Integer.
!         Current iteration count.  (Must be zero on first call.)
! ERR    :OUT      Double Precision.
!         Error estimate of error in the X(N) approximate solution, as
!         defined by ITOL.
! IERR   :OUT      Integer.
!         Error flag.  IERR is set to 3 if ITOL is not one of the
!         acceptable values, see above.
! IUNIT  :IN       Integer.
!         Unit number on which to write the error at each iteration,
!         if this is desired for monitoring convergence.  If unit
!         number is 0, no writing will occur.
! R      :IN       Double Precision R(N).
!         The residual R = B-AX.
! Z      :WORK     Double Precision Z(N).
!         Workspace used to hold the pseudo-residual M Z = R.
! P      :IN       Double Precision P(N).
!         The conjugate direction vector.
! DZ     :WORK     Double Precision DZ(N).
!         Workspace used to hold temporary vector(s).
! RWORK  :WORK     Double Precision RWORK(USER DEFINED).
!         Double Precision array that can be used by MSOLVE.
! IWORK  :WORK     Integer IWORK(USER DEFINED).
!         Integer array that can be used by MSOLVE.
! AK     :IN       Double Precision.
! BK     :IN       Double Precision.
!         Current conjugate gradient parameters alpha and beta.
! BNRM   :INOUT    Double Precision.
!         Norm of the right hand side.  Type of norm depends on ITOL.
!         Calculated only on the first call.
! SOLNRM :INOUT    Double Precision.
!         2-Norm of the true solution, SOLN.  Only computed and used
!         if ITOL = 11.
! 
! *Function Return Values:
!       0 : Error estimate (determined by ITOL) is *NOT* less than the
!           specified tolerance, TOL.  The iteration must continue.
!       1 : Error estimate (determined by ITOL) is less than the
!           specified tolerance, TOL.  The iteration can be considered
!           complete.
! 
! *Cautions:
!     This routine will attempt to write to the Fortran logical output
!     unit IUNIT, if IUNIT .ne. 0.  Thus, the user must make sure that
!     this logical unit is attached to a file or terminal before calling
!     this routine with a non-zero value for IUNIT.  This routine does
!     not check for the validity of a non-zero IUNIT unit number.
! 
!***SEE ALSO  DCG, DSDCG, DSICCG
!***ROUTINES CALLED  D1MACH, DNRM2
!***COMMON BLOCKS    DSLBLK
!***REVISION HISTORY  (YYMMDD)
!   890404  DATE WRITTEN
!   890404  Previous REVISION DATE
!   890915  Made changes requested at July 1989 CML Meeting.  (MKS)
!   890921  Removed TeX from comments.  (FNF)
!   890922  Numerous changes to prologue to make closer to SLATEC
!           standard.  (FNF)
!   890929  Numerous changes to reduce SP/DP differences.  (FNF)
!   891003  Removed C***REFER TO line, per MKS.
!   910411  Prologue converted to Version 4.0 format.  (BAB)
!   910502  Removed MSOLVE from ROUTINES CALLED list.  (FNF)
!   910506  Made subsidiary to DCG.  (FNF)
!   920407  COMMON BLOCK renamed DSLBLK.  (WRB)
!   920511  Added complete declaration section.  (WRB)
!   920930  Corrected to not print AK,BK when ITER=0.  (FNF)
!   921026  Changed 1.0E10 to D1MACH(2) and corrected D to E in
!           output format.  (FNF)
! ***END PROLOGUE  ISDCG
!     .. Scalar Arguments ..
      INTEGER  ISDCG 
      DOUBLE PRECISION AK, BK, BNRM, ERR, SOLNRM, TOL
      INTEGER IERR, ISYM, ITER, ITMAX, ITOL, IUNIT, N, NELT
!     .. Array Arguments ..
      DOUBLE PRECISION A(NELT), B(N), DZ(N), P(N), R(N), RWORK(*), X(N), Z(N)
      INTEGER IA(NELT), IWORK(*), JA(NELT)
!     .. Subroutine Arguments ..
!      EXTERNAL MSOLVE  

      INTERFACE 
          SUBROUTINE MSOLVE(N, R, Z, NELT, IA, JA, A, ISYM, RWORK, IWORK)
      INTEGER N, ISYM, NELT
       DOUBLE PRECISION A(NELT) , RWORK(*) ,Z(N), R(N)  
      INTEGER IA(NELT), IWORK(*), JA(NELT)

          END SUBROUTINE MSOLVE
     END INTERFACE


!     .. Arrays in Common ..
      DOUBLE PRECISION SOLN(1)
!     .. Local Scalars ..
      INTEGER I
!     .. External Functions ..
!      DOUBLE PRECISION D1MACH, DNRM2
!      EXTERNAL D1MACH, DNRM2
!     .. Common blocks ..
!      COMMON /DSLBLK/ SOLN
! ***FIRST EXECUTABLE STATEMENT  ISDCG
      ISDCG = 0
! 
      IF( ITOL.EQ.1 ) THEN
!         err = ||Residual||/||RightHandSide|| (2-Norms).
         IF(ITER .EQ. 0) BNRM = DNRM2(N, B, 1)
         ERR = DNRM2(N, R, 1)/BNRM
      ELSE IF( ITOL.EQ.2 ) THEN
!                  -1              -1
!         err = ||M  Residual||/||M  RightHandSide|| (2-Norms).
         IF(ITER .EQ. 0) THEN
            CALL MSOLVE(N, R, Z, NELT, IA, JA, A, ISYM, RWORK, IWORK)
            BNRM = DNRM2(N, DZ, 1)
         ENDIF
         ERR = DNRM2(N, Z, 1)/BNRM
      ELSE IF( ITOL.EQ.11 ) THEN
!         err = ||x-TrueSolution||/||TrueSolution|| (2-Norms).
         IF(ITER .EQ. 0) SOLNRM = DNRM2(N, SOLN, 1)
         DO 10 I = 1, N
            DZ(I) = X(I) - SOLN(I)
 10      CONTINUE
         ERR = DNRM2(N, DZ, 1)/SOLNRM
      ELSE
! 
!         If we get here ITOL is not one of the acceptable values.
         ERR = D1MACH(2)
         IERR = 3
      ENDIF
! 
      IF(IUNIT .NE. 0) THEN
         IF( ITER.EQ.0 ) THEN
            WRITE(IUNIT,1000) N, ITOL
            WRITE(IUNIT,1010) ITER, ERR
         ELSE
            WRITE(IUNIT,1010) ITER, ERR, AK, BK
         ENDIF
      ENDIF
      IF(ERR .LE. TOL) ISDCG = 1
      RETURN
 1000 FORMAT(' Preconditioned Conjugate Gradient for ',&
          'N, ITOL = ',I5, I5,&
          /' ITER','   Error Estimate','            Alpha',&
          '             Beta')
 1010 FORMAT(1X,I4,1X,D16.7,1X,D16.7,1X,D16.7)
! ------------- LAST LINE OF ISDCG FOLLOWS ------------------------------
      END FUNCTION ISDCG
!-----------------------------------------------------------------------
!DECK DAXPY
      SUBROUTINE DAXPY (N, DA, DX, INCX, DY, INCY)
!***BEGIN PROLOGUE  DAXPY
!***PURPOSE  Compute a constant times a vector plus a vector.
!***LIBRARY   SLATEC (BLAS)
!***CATEGORY  D1A7
!***TYPE      DOUBLE PRECISION (SAXPY-S, DAXPY-D, CAXPY-C)
!***KEYWORDS  BLAS, LINEAR ALGEBRA, TRIAD, VECTOR
!***AUTHOR  Lawson, C. L., (JPL)
!           Hanson, R. J., (SNLA)
!           Kincaid, D. R., (U. of Texas)
!           Krogh, F. T., (JPL)
!***DESCRIPTION
! 
!                B L A S  Subprogram
!    Description of Parameters
! 
!     --Input--
!        N  number of elements in input vector(s)
!       DA  double precision scalar multiplier
!       DX  double precision vector with N elements
!     INCX  storage spacing between elements of DX
!       DY  double precision vector with N elements
!     INCY  storage spacing between elements of DY
! 
!     --Output--
!       DY  double precision result (unchanged if N .LE. 0)
! 
!     Overwrite double precision DY with double precision DA*DX + DY.
!     For I = 0 to N-1, replace  DY(LY+I*INCY) with DA*DX(LX+I*INCX) +
!       DY(LY+I*INCY),
!     where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is
!     defined in a similar way using INCY.
! 
! ***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
!                 Krogh, Basic linear algebra subprograms for Fortran
!                 usage, Algorithm No. 539, Transactions on Mathematical
!                 Software 5, 3 (September 1979), pp. 308-323.
! ***ROUTINES CALLED  (NONE)
! ***REVISION HISTORY  (YYMMDD)
!   791001  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920310  Corrected definition of LX in DESCRIPTION.  (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
! ***END PROLOGUE  DAXPY
      INTEGER N, INCX, INCY
      DOUBLE PRECISION DX(*), DY(*), DA
      INTEGER IX, IY, I, M, MP1, NS
! ***FIRST EXECUTABLE STATEMENT  DAXPY
      IF (N.LE.0 .OR. DA.EQ.0.0D0) RETURN
      IF (INCX .EQ. INCY) IF (INCX-1) 5,20,60
! 
!     Code for unequal or nonpositive increments.
! 
    5 IX = 1
      IY = 1
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
      IF (INCY .LT. 0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DY(IY) = DY(IY) + DA*DX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
! 
!     Code for both increments equal to 1.
! 
!     Clean-up loop so remaining vector length is a multiple of 4.
! 
   20 M = MOD(N,4)
      IF (M .EQ. 0) GO TO 40
      DO 30 I = 1,M
        DY(I) = DY(I) + DA*DX(I)
   30 CONTINUE
      IF (N .LT. 4) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,4
        DY(I) = DY(I) + DA*DX(I)
        DY(I+1) = DY(I+1) + DA*DX(I+1)
        DY(I+2) = DY(I+2) + DA*DX(I+2)
        DY(I+3) = DY(I+3) + DA*DX(I+3)
   50 CONTINUE
      RETURN
! 
!     Code for equal, positive, non-unit increments.
! 
   60 NS = N*INCX
      DO 70 I = 1,NS,INCX
        DY(I) = DA*DX(I) + DY(I)
   70 CONTINUE
      RETURN
      END SUBROUTINE DAXPY






END MODULE dsiccg_icholmodule
