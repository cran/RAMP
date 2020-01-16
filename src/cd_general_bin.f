! --------------------
! this function implements the CD for lasso
! algorithm for the following optimization problem
! paraIn includes: epsilon, maxIter, lambda, gamma
! ptype = 1: lasso, 2: mcp, 3: scad


      subroutine cd_general_bin(X, y, p, n, beta, pf, ptype, paraIn)
       
      integer :: i, l, p, n, ptype  
      double precision, dimension(n,p) :: X
      double precision, dimension(n) :: y, r, Xl, W, eta, pi
      double precision, dimension(p) :: blast, beta, pf
      double precision, dimension(4) :: paraIn
      double precision :: lambda, v, z, thresh, gamma
      double precision :: tmp, lam, epsilon
     
      epsilon = paraIn(1)
      maxIter = INT(paraIn(2))
      lam = paraIn(3)
      gamma = paraIn(4)
      eta = MATMUL(X, beta)
      pi = 1/(1+EXP(-eta)) 
      r = y - pi 
      breakInd = 0 

      do  i = 1, maxIter
      
           blast = beta          
          
           DO   l = 1, p
                  lambda = lam*pf(l)
                  Xl = X(:,l)
                  W = pi*(1-pi)
                  v=SUM(Xl * W *  Xl)/n
                  z=SUM(Xl * (y-pi))/n
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
                  
                  
                  eta = eta + (beta(l)-blast(l))*Xl
                  pi = 1/(1+EXP(-eta))
                  r = y - pi 
                  IF(SUM(ABS(r))<epsilon*p) THEN
                  ! print *, SUM(ABS(r))
                  breakInd = 1
                  EXIT
                  END IF
                   
          END DO      
          IF (breakInd ==1) THEN 
          		EXIT                   
          END IF
          IF (maxIter /= 1) THEN
             IF(SUM(ABS(blast-beta))<epsilon*p) THEN
                  !print *, "fortran loop =         ", i, beta(1:5)
                  EXIT                   
              END IF
          END IF 
          
              
      end do

       
      end
