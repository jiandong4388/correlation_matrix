library(expm)


       GFT_inverse_mapping <- function(gamma_in, tol_value){
         C <- matrix(0,nrow=0,ncol=0)
         iter_number <- -1
         
         # Check if input is of proper format : gamma is of suitable length
         # and tolerance value belongs to a proper interval
         n <- 0.5*(1+sqrt(1+8*length(gamma_in)))
         if (!(is.vector(gamma_in) && n%%1==0)){
           (stop("Dimension of 'gammar is incorrect") )
            }else if (!(tol_value>=1e-14 && tol_value<=1e-4)){
              (stop("Incorrect tolerance value")) }
         else{
           
           # Place elements from gamma into off-diagonal parts
           # and put zeros on the main diagonal of nxn symmetric matrix A
           A <- matrix(0,nrow=n,ncol=n)
           A[upper.tri(A, diag=FALSE)] <- gamma_in
           A <- A + t(A)
           
           # Read properties ofthe input matrix
           diag_vec <- diag(A)
           
           # Iterative algorithm to get the proper diagonal vector
           dist <- sqrt(n)
           while(dist > sqrt(n)*tol_value ){
             diag_delta <- log(diag(expm(A)))
             diag_vec <- diag_vec - diag_delta
             diag(A) <- diag_vec
             dist <- norm(diag_delta, type="2")
             iter_number<- iter_number + 1
         }
         # Get a unique reciprocal correlation matrix
         C <- expm(A)
         diag(C)<- rep(1,n)
           
       } 
       return(list(C=C, iter_number=iter_number))

   }