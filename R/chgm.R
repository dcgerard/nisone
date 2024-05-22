## CHGM computes the confluent hypergeometric function M(a,b,x).
#
#  Licensing:
#
#    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However,
#    they give permission to incorporate this routine into a user program
#    provided that the copyright is acknowledged.
#
#  Modified:
#
#    27 July 2012
#
#  Author:
#
#    Shanjie Zhang, Jianming Jin
#
#  Reference:
#
#    Shanjie Zhang, Jianming Jin,
#    Computation of Special Functions,
#    Wiley, 1996,
#    ISBN: 0-471-11963-6,
#    LC: QA351.C45.
chgm <- function(a, b, x) {
  a0 <- a
  a1 <- a
  x0 <- x
  hg <- 0.0

  if (b == 0.0 || b == -abs(as.integer(b))) {
    hg <- 1.0e300
  } else if (a == 0.0 || x == 0.0) {
    hg <- 1.0
  } else if (a == -1.0) {
    hg <- 1.0 - x / b
  } else if (a == b) {
    hg <- exp(x)
  } else if (a - b == 1.0) {
    hg <- (1.0 + x / b) * exp(x)
  } else if (a == 1.0 && b == 2.0) {
    hg <- (exp(x) - 1.0) / x
  } else if (a == as.integer(a) && a < 0.0) {
    m <- as.integer(-a)
    r <- 1.0
    hg <- 1.0
    for (k in 1:m) {
      r <- r * (a + k - 1.0) / k / (b + k - 1.0) * x
      hg <- hg + r
    }
  }

  if (hg != 0.0) {
    return(hg)
  }

  if (x < 0.0) {
    a <- b - a
    a0 <- a
    x <- abs(x)
  }

  nl <- if (a < 2.0) 0 else 1

  if (a >= 2.0) {
    la <- as.integer(a)
    a <- a - la - 1.0
  }

  for (n in 0:nl) {
    if (a0 >= 2.0) {
      a <- a + 1.0
    }

    if (x <= 30.0 + abs(b) || a < 0.0) {
      hg <- 1.0
      rg <- 1.0
      for (j in 1:500) {
        rg <- rg * (a + j - 1.0) / (j * (b + j - 1.0)) * x
        hg <- hg + rg
        if (abs(rg / hg) < 1.0e-15) {
          break
        }
      }
    } else {
      ta <- gamma(a)
      tb <- gamma(b)
      xg <- b - a
      tba <- gamma(xg)
      sum1 <- 1.0
      sum2 <- 1.0
      r1 <- 1.0
      r2 <- 1.0
      for (i in 1:8) {
        r1 <- -r1 * (a + i - 1.0) * (a - b + i) / (x * i)
        r2 <- -r2 * (b - a + i - 1.0) * (a - i) / (x * i)
        sum1 <- sum1 + r1
        sum2 <- sum2 + r2
      }
      hg1 <- tb / tba * x^(-a) * cos(pi * a) * sum1
      hg2 <- tb / ta * exp(x) * x^(a - b) * sum2
      hg <- hg1 + hg2
    }

    if (n == 0) {
      y0 <- hg
    } else if (n == 1) {
      y1 <- hg
    }
  }

  if (a0 >= 2.0) {
    for (i in 1:(la - 1)) {
      hg <- ((2.0 * a - b + x) * y1 + (b - a) * y0) / a
      y0 <- y1
      y1 <- hg
      a <- a + 1.0
    }
  }

  if (x0 < 0.0) {
    hg <- hg * exp(x0)
  }

  a <- a1
  x <- x0

  return(hg)
}

## Original fortran code:
#     subroutine chgm ( a, b, x, hg ) bind(C, name = "chgm_")
#
#     !*****************************************************************************80
#     !
#     !! CHGM computes the confluent hypergeometric function M(a,b,x).
#     !
#     !  Licensing:
#     !
#     !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However,
#     !    they give permission to incorporate this routine into a user program
#     !    provided that the copyright is acknowledged.
#     !
#     !  Modified:
#     !
#     !    27 July 2012
#     !
#     !  Author:
#     !
#     !    Shanjie Zhang, Jianming Jin
#     !
#     !  Reference:
#     !
#     !    Shanjie Zhang, Jianming Jin,
#     !    Computation of Special Functions,
#     !    Wiley, 1996,
#     !    ISBN: 0-471-11963-6,
#     !    LC: QA351.C45.
#     !
#     !  Parameters:
#     !
#     !    Input, real ( kind = rk ) A, B, parameters.
#     !
#     !    Input, real ( kind = rk ) X, the argument.
#     !
#     !    Output, real ( kind = rk ) HG, the value of M(a,b,x).
#     !
#       integer, parameter :: rk = kind ( 1.0D+00 )
#
#
#       real(c_double) :: a
#       real ( kind = rk ) a0
#       real ( kind = rk ) a1
#       real(c_double) :: b
#       real(c_double) :: hg
#       real ( kind = rk ) hg1
#       real ( kind = rk ) hg2
#       integer i
#       integer j
#       integer k
#       integer la
#       integer m
#       integer n
#       integer nl
#       real ( kind = rk ) pi
#       real ( kind = rk ) r
#       real ( kind = rk ) r1
#       real ( kind = rk ) r2
#       real ( kind = rk ) rg
#       real ( kind = rk ) sum1
#       real ( kind = rk ) sum2
#       real ( kind = rk ) ta
#       real ( kind = rk ) tb
#       real ( kind = rk ) tba
#       real(c_double) :: x
#       real ( kind = rk ) x0
#       real ( kind = rk ) xg
#       real ( kind = rk ) y0
#       real ( kind = rk ) y1
#
#       pi = 3.141592653589793D+00
#       a0 = a
#       a1 = a
#       x0 = x
#       hg = 0.0D+00
#
#       if ( b == 0.0D+00 .or. b == - abs ( int ( b ) ) ) then
#         hg = 1.0D+300
#       else if ( a == 0.0D+00 .or. x == 0.0D+00 ) then
#         hg = 1.0D+00
#       else if ( a == -1.0D+00 ) then
#         hg = 1.0D+00 - x / b
#       else if ( a == b ) then
#         hg = exp ( x )
#       else if ( a - b == 1.0D+00 ) then
#         hg = ( 1.0D+00 + x / b ) * exp ( x )
#       else if ( a == 1.0D+00 .and. b == 2.0D+00 ) then
#         hg = ( exp ( x ) - 1.0D+00 ) / x
#       else if ( a == int ( a ) .and. a < 0.0D+00 ) then
#         m = int ( - a )
#         r = 1.0D+00
#         hg = 1.0D+00
#         do k = 1, m
#           r = r * ( a + k - 1.0D+00 ) / k / ( b + k - 1.0D+00 ) * x
#           hg = hg + r
#         end do
#       end if
#
#       if ( hg /= 0.0D+00 ) then
#         return
#       end if
#
#       if ( x < 0.0D+00 ) then
#         a = b - a
#         a0 = a
#         x = abs ( x )
#       end if
#
#       if ( a < 2.0D+00 ) then
#         nl = 0
#       end if
#
#       if ( 2.0D+00 <= a ) then
#         nl = 1
#         la = int ( a )
#         a = a - la - 1.0D+00
#       end if
#
#       do n = 0, nl
#
#         if ( 2.0D+00 <= a0 ) then
#           a = a + 1.0D+00
#         end if
#
#         if ( x <= 30.0D+00 + abs ( b ) .or. a < 0.0D+00 ) then
#
#           hg = 1.0D+00
#           rg = 1.0D+00
#           do j = 1, 500
#             rg = rg * ( a + j - 1.0D+00 ) &
#               / ( j * ( b + j - 1.0D+00 ) ) * x
#             hg = hg + rg
#             if ( abs ( rg / hg ) < 1.0D-15 ) then
#               exit
#             end if
#           end do
#
#         else
#
#           call gamma ( a, ta )
#           call gamma ( b, tb )
#           xg = b - a
#           call gamma ( xg, tba )
#           sum1 = 1.0D+00
#           sum2 = 1.0D+00
#           r1 = 1.0D+00
#           r2 = 1.0D+00
#           do i = 1, 8
#             r1 = - r1 * ( a + i - 1.0D+00 ) * ( a - b + i ) / ( x * i )
#             r2 = - r2 * ( b - a + i - 1.0D+00 ) * ( a - i ) / ( x * i )
#             sum1 = sum1 + r1
#             sum2 = sum2 + r2
#           end do
#           hg1 = tb / tba * x ** ( - a ) * cos ( pi * a ) * sum1
#           hg2 = tb / ta * exp ( x ) * x ** ( a - b ) * sum2
#           hg = hg1 + hg2
#
#         end if
#
#         if ( n == 0 ) then
#           y0 = hg
#         else if ( n == 1 ) then
#           y1 = hg
#         end if
#
#       end do
#
#       if ( 2.0D+00 <= a0 ) then
#         do i = 1, la - 1
#           hg = ( ( 2.0D+00 * a - b + x ) * y1 + ( b - a ) * y0 ) / a
#           y0 = y1
#           y1 = hg
#           a = a + 1.0D+00
#         end do
#       end if
#
#       if ( x0 < 0.0D+00 ) then
#         hg = hg * exp ( x0 )
#       end if
#
#       a = a1
#       x = x0
#
#       return
#     end subroutine chgm
#
# end module fortfuncs
