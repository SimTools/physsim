CC********************************************************************CC
C*                                                                    *C
C*=================                                                   *C
C* Program PRMAIN                                                     *C
C*=================                                                   *C
C*(Purpose)                                                           *C
C*       Main program of TOPAZ Production job.                        *C
C*(Relation)                                                          *C
C*       Calls the following subroutines;                             *C
C*           PRJOBI : initializes this job.                           *C
C*           PRRUN  : analyzes a run in run loop.                     *C
C*           PRJOBT : terminates this job.                            *C
C*(Author)                                                            *C
C*   14-Oct-1985  A. Miyamoto  Original Version.                      *C
C*   25-Apr-1987  A. Miyamoto  For Full system                        *C
C*                                                                    *C
CC********************************************************************CC
C
C ======================
      PROGRAM PRMAIN
C ======================
C
      EXTERNAL PRBLK1, UCONST, BS51DT
C
C=========< Entry Point >===============================================
C
C -------------------------------------------------------------------- C
C     Initialize this job.                                             C
C -------------------------------------------------------------------- C
C
      CALL PRINIT
C
C -------------------------------------------------------------------- C
C     Describe Analysis routine.                                       C
C -------------------------------------------------------------------- C
C
      CALL FLANLI
C
C -------------------------------------------------------------------- C
C     Do job initialization task.                                      C
C -------------------------------------------------------------------- C
C
      CALL PRJOBI
C
C -------------------------------------------------------------------- C
C     Process Event.
C -------------------------------------------------------------------- C
C
      CALL PRCEVT
C
C -------------------------------------------------------------------- C
C     Terminate Production job.
C -------------------------------------------------------------------- C
C
      CALL PRJOBT
C
      PRINT *,' Production Job normal end.'
      STOP
C
      END
