! --------------------
! this function implements the CD for lasso
! algorithm for the following optimization problem
! para.in includes: epsilon, maxIter, lambda, 


      subroutine cd_lasso_lin(X, y, p, n, beta, nonPen,  paraIn)
       
      integer :: i, l, p, n     
      integer, dimension(p) :: nonPen
      double precision, dimension(n,p) :: X
      double precision, dimension(n) :: y, r, Xl
      double precision, dimension(p) :: blast, beta
      double precision, dimension(4) :: paraIn
      double precision :: lambda, v, z, thresh
     
      epsilon = paraIn(1)
      maxIter = INT(paraIn(2))
      lambda = paraIn(3)
       
      r = y - MATMUL(X, beta) 
      do  i = 1, maxIter
      
           blast = beta          
          
           DO   l = 1, p
                  Xl = X(:,l)
                  v=SUM(Xl * Xl)/n
                  z=SUM(Xl *r)/n
                  z=z + v * beta(l)
                  thresh = MAX(dble(0), ABS(z)-lambda)
                  
                  IF(nonPen(l) == 0) THEN
                      beta(l)=SIGN(dble(1),z)*thresh/v
                  ELSE 
                      beta(l)=z/v
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
