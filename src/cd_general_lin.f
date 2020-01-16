! --------------------
! this function implements the CD for lasso
! algorithm for the following optimization problem
! paraIn includes: epsilon, maxIter, lambda, gamma
! ptype = 1: lasso, 2: mcp, 3: scad


      subroutine cd_general_lin(X, y, p, n, beta, pf, ptype, paraIn)
       
      integer :: i, l, p, n, ptype
      double precision, dimension(n,p) :: X
      double precision, dimension(n) :: y, r, Xl
      double precision, dimension(p) :: blast, beta, pf
      double precision, dimension(4) :: paraIn
      double precision :: lam, lambda, v, z, thresh, gamma
      double precision :: tmp, epsilon
     
      epsilon = paraIn(1)
      maxIter = INT(paraIn(2))
      lam = paraIn(3)
      gamma = paraIn(4)
       
      r = y - MATMUL(X, beta) 
      do  i = 1, maxIter
      
           blast = beta          
          
           DO   l = 1, p
                  lambda = lam*pf(l)
                  Xl = X(:,l)
                  v=SUM(Xl * Xl)/n
                  z=SUM(Xl *r)/n
                  z=z + v * beta(l)
                  thresh = MAX(dble(0), ABS(z)-lambda)
                  
                  IF(ptype == 1) THEN
                      beta(l)=SIGN(dble(1),z)*thresh/v
                  ELSE IF(ptype == 2) THEN
                  IF(ABS(z)<lambda*gamma) THEN
                  beta(l)=SIGN(dble(1),z)*thresh/v/(1-1/gamma)
                  ELSE
                  beta(l)=z/v
                  END IF
                  ELSE IF(ptype == 3) THEN
                  IF(ABS(z) < 2 * lambda) THEN
                      beta(l)=SIGN(dble(1),z)*thresh/v
                  ELSE IF(ABS(z)<lambda*gamma) THEN
                  tmp = ABS(z)-lambda*gamma/(gamma-1)
                  thresh = MAX(dble(0), tmp)/(1-1/(gamma-1))
                      beta(l)=SIGN(dble(1),z)*thresh/v
                  ELSE
                      beta(l)=z/v
                  END IF
                  END IF
                  r = r - (beta(l)-blast(l))*Xl
                   
          END DO      
          
          IF (maxIter /= 1) THEN
             IF(SUM(ABS(blast-beta))<epsilon*p) THEN
                  !print *, "fortran loop =         ", i, beta(1:5)
                  EXIT                   
              END IF
          END IF 
          
              
      end do

     
      end
