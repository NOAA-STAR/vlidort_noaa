MODULE PCRTM_Interp_Utility

  IMPLICIT NONE

  CONTAINS

    SUBROUTINE Interp_Linear_Dble(nref, xref, yref, n, x, y)

    INTEGER,               INTENT(IN) :: nref
    REAL*8, DIMENSION(nref), INTENT(IN) :: xref,yref
    INTEGER,               INTENT(IN) :: n

    REAL*8, DIMENSION(n),    INTENT(IN) :: x
    REAL*8, DIMENSION(n),    INTENT(out):: y

    INTEGER :: i, r, rlo
    REAL*8    :: slope, dx

    rlo = 1
    DO i = 1, n
      IF (x(i) <= xref(1)) THEN
        y(i) = yref(1)
      ELSEIF (x(i) >= xref(nref)) THEN
        y(i) = yref(nref)
      ELSE

        DO r =1, nref-1
          IF ((xref(r) <= x(i)) .AND. (x(i) < xref(r+1))) THEN
            rlo = r
            EXIT
          END IF
        END DO

        slope = (yref(r+1) - yref(r)) / (xref(r+1) - xref(r))
        dx =  x(i) - xref(r)
        y(i) = yref(r) + slope*dx
      END IF
    END DO

    IF (.false.) write(*,*) rlo
  END SUBROUTINE interp_linear_Dble

  SUBROUTINE Interp_Linear(nref, xref, yref, n, x, y)

    INTEGER,               INTENT(IN) :: nref
    REAL, DIMENSION(nref), INTENT(IN) :: xref,yref
    INTEGER,               INTENT(IN) :: n

    REAL, DIMENSION(n),    INTENT(IN) :: x
    REAL, DIMENSION(n),    INTENT(out):: y

    INTEGER :: i, r, rlo
    REAL    :: slope, dx

    rlo = 1
    DO i = 1, n
      IF (x(i) <= xref(1)) THEN
        y(i) = yref(1)
      ELSEIF (x(i) >= xref(nref)) THEN
        y(i) = yref(nref)
      ELSE

        DO r =1, nref-1
          IF ((xref(r) <= x(i)) .AND. (x(i) < xref(r+1))) THEN
            rlo = r
            EXIT
          END IF
        END DO

        slope = (yref(r+1) - yref(r)) / (xref(r+1) - xref(r))
        dx =  x(i) - xref(r)
        y(i) = yref(r) + slope*dx
      END IF
    ENDDO

    IF (.false.) write(*,*) rlo

  END SUBROUTINE Interp_Linear


  SUBROUTINE quadterp (xtab, ytab, xint, yint)
    IMPLICIT NONE
    REAL, DIMENSION(0:2)   :: xtab, ytab
    REAL                   :: xint,yint
    REAL                   :: x0,x1,x2
    REAL                   :: y0,y1,y2

    yint = ( (xint-xtab(1)) * ( xint-xtab(2)) )/ ( (xtab(0)-xtab(1)) * (xtab(0)-xtab(2)) ) * ytab(0) + &
         ( (xint-xtab(0)) * ( xint-xtab(2)) )/ ( (xtab(1)-xtab(0)) * (xtab(1)-xtab(2)) ) * ytab(1) + &
         ( (xint-xtab(0)) * ( xint-xtab(1)) )/ ( (xtab(2)-xtab(0)) * (xtab(2)-xtab(1)) ) * ytab(2)


  END SUBROUTINE quadterp


  SUBROUTINE quadterp2 (xtab, ytab, xint, yint)
    IMPLICIT NONE
    REAL*8, DIMENSION(0:2)   :: xtab, ytab
    REAL*8                   :: xint,yint
    REAL*8                   :: x0,x1,x2
    REAL*8                   :: y0,y1,y2

    yint = ( (xint-xtab(1)) * ( xint-xtab(2)) )/ ( (xtab(0)-xtab(1)) * (xtab(0)-xtab(2)) ) * ytab(0) + &
         ( (xint-xtab(0)) * ( xint-xtab(2)) )/ ( (xtab(1)-xtab(0)) * (xtab(1)-xtab(2)) ) * ytab(1) + &
         ( (xint-xtab(0)) * ( xint-xtab(1)) )/ ( (xtab(2)-xtab(0)) * (xtab(2)-xtab(1)) ) * ytab(2)


  END SUBROUTINE quadterp2



END MODULE PCRTM_Interp_Utility
