! --------------------
! this function implements the CD for lasso
! algorithm for the logistic regression problem
! para.in includes: epsilon, maxIter, lambda, 


      subroutine cd_lasso_bin(X, y, p, n, beta, nonPen,  paraIn)
       
      integer :: i, l, p, n, breakInd
      integer, dimension(p) :: nonPen
      double precision, dimension(n,p) :: X
      double precision, dimension(n) :: y, r, Xl, W, eta, pi
      double precision, dimension(p) :: blast, beta
      double precision, dimension(4) :: paraIn
      double precision :: lambda, v, z, thresh
     
      epsilon = paraIn(1)
      maxIter = INT(paraIn(2))
      lambda = paraIn(3)
      eta = MATMUL(X, beta)
      pi = 1/(1+EXP(-eta)) 
      r = y - pi 
      breakInd = 0
      !print *, SUM(ABS(r))
      !print *, SUM(pi)/n
      !print *,pi
      !print *,r
      do  i = 1, maxIter
      
           blast = beta          
          
           DO   l = 1, p
                  Xl = X(:,l)
                  W = pi*(1-pi)
                  v=SUM(Xl * W *  Xl)/n
                  z=SUM(Xl * (y-pi))/n
                  z=z + v * beta(l)
                  thresh = MAX(dble(0), ABS(z)-lambda)
                  
                  IF(nonPen(l) == 0) THEN
                      beta(l)=SIGN(dble(1),z)*thresh/v
                  ELSE 
                      beta(l)=z/v
                  END IF
                  eta = eta + (beta(l)-blast(l))*Xl
                  pi = 1/(1+EXP(-eta))
                  r = y - pi
                  !if(MOD(l,10)==0) THEN
                  !print *, SUM(ABS(r))
                  !END IF
                  if(SUM(ABS(r))<epsilon*p) THEN
                  ! print *, SUM(ABS(r))
                  breakInd = 1
                  EXIT
                  END IF
                                     
          END DO 
          !print *, SUM(ABS(r))     
          IF (breakInd ==1) THEN 
          		EXIT                   
          END IF
          IF (maxIter /= 1) THEN
          !print *, beta
             IF(SUM(ABS(blast-beta))<epsilon*p) THEN
                  EXIT                   
              END IF
          END IF 
              
      end do

      end
