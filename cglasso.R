######   Solve Covariance Graphical Lasso Model using Coordinate Descent algorithms
CovGlassoCD= function(Q,Rho,Sig,outertol,innertol,outermax,innermax)
{
# CovGlassoCD solves the covariance graphical lasso problem (Bien and Tibshirani 2011, Biometrika):
# Sigma_cd  = argmin { logdet(Sigma)+trace(Q/Sig)+sum(sum(Rho.*abs(Sig))) }
# using coordinate descent algorithm (Wang 2012)
# # Input:
#   Q:  p*p sample covariance matrix  Q = Y*Y'/n; 
#   Rho: p*p shrinkage parameter matrix
#   Sig: inital value
#   outertol: threshold for the change of objective functions for terminating outer iterations
#   innertol: threshold for the maximum change of beta for terminating inner loop
#   outermax: maximum number of itermations of the outer loop
#   innermax: maximum number of iterations of the inner loop for the lasso problem
#   
#   Output:
#   Sig: clustering point of Sig
#   C:   inv(Sig)
#   obj: sequence of the objective function at each iteration: 
#                 obj = logdet(Sigma)+trace(Q/Sig)+sum(sum(Rho.*abs(Sig)));
#   niter: number of iterations


# Written by Hao Wang @ U of South Carolina

 p = nrow(Q);
   C = solve(Sig); # Initial value

   
ind.noi.all = mat.or.vec(p-1,p);
for( i in 1:p)
  {

    ind.noi = seq(1,p);
    ind.noi = ind.noi[-i];
    ind.noi.all[,i] = ind.noi; 
   
  }

  ob = -determinant(C,logarithm = TRUE)$modulus+sum(C*Q)+sum(Rho*abs(Sig));

  obj <- ob;
   
#    print(paste('iter = ',0,'ob=',ob));

for(iter1 in 1:outermax)
  {


    ob.old = ob;
    Sig.old = abs(Sig);


    for(i in 1:p)
      {

        ind.noi <- ind.noi.all[,i];

        C11 <- C[ind.noi,ind.noi]; C12 = C[ind.noi,i];
        S11 <- Q[ind.noi,ind.noi]; S12 = Q[ind.noi,i];
        Sig11 <- Sig[ind.noi,ind.noi];

        invSig11 <- C11-C12%*%t(C12)/C[i,i];
        invSig11S12 <- invSig11%*%S12;

        S11invSig11 <- S11%*%invSig11;
        W1 <- invSig11%*%S11invSig11;

####### Update gam
        beta <- Sig[ind.noi,i];
        a <- 1/2; b <- Rho[i,i]/2; c <- (t(beta)%*%W1%*%beta-2*t(beta)%*%invSig11S12+Q[i,i])/2;

        {if (b==0)
          {gam <- c/a;}
        else
          {gam <- (-a+sqrt(a*a+4*b*c))/(2*b);}
         }


        
####### Update beta          
        V <- W1/as.numeric(gam)+Rho[i,i]*invSig11;
        u = invSig11S12/as.numeric(gam);


        Vbeta = V%*%beta;
     for(iter2 in 1:innermax)
       {
        beta.old <- beta;
        
        for(j in 1:(p-1))
        {
          betaj.old <- beta[j];
          x = u[j]-Vbeta[j]+V[j,j]*beta[j];
          
#          {
#          if(x == 0)
#            {signx = runif(1);}
#          else
#            {signx = sign(x)}
#           }
          beta[j] <- max(0,abs(x)-Rho[i,ind.noi[j]])*sign(x)/V[j,j];

        beta.change = beta[j]-betaj.old;
       if(beta.change!=0)
         {
           Vbeta = Vbeta+beta.change*V[,j];
         }
          
  
        }

         if(max(abs(beta.old-beta))<innertol)
           {break;}
       
      }

        Sig[ind.noi,i] <- beta;
        Sig[i,ind.noi] <- beta;
        Sig[i,i] <- gam+t(beta)%*%invSig11%*%beta;

        
        invSig11beta <-  invSig11%*%beta;
        C[ind.noi,ind.noi] <- invSig11+invSig11beta%*%t(invSig11beta)/as.numeric(gam);
       
        C12 <- -invSig11beta/as.numeric(gam);

        C[ind.noi,i]=C12;
        C[i,ind.noi]=C12;

        C[i,i] = 1/gam
        
      }
    
      ob <- -determinant(C,logarithm = TRUE)$modulus+sum(C*Q)+sum(Rho*abs(Sig));
    obj <- c(obj,ob);
#    print(paste('iter = ', iter1,'ob=',ob));
      if(abs(ob-ob.old)<outertol)
        {break;}

    
  }

  return(list(Sig = Sig, C = C, obj = obj, niter = iter1))

}






######   Solve Covariance Graphical Lasso Model using ECM algorithms
CovGlassoECM = function(Q,Rho,Sig,tol,itermax)
{
  
# CovGlassoECM solves the covariance graphical lasso problem (Bien and Tibshirani 2011, Biometrika):
# Sigma_cd  = argmin { logdet(Sigma)+trace(Q/Sig)+sum(sum(Rho.*abs(Sig))) }
# using ECM algorithm (Wang 2012)
# # Input:
#   Q:  p*p sample covariance matrix  Q = Y*Y'/n; 
#   Rho: p*p shrinkage parameter matrix
#   Sig: inital value
#   tol: threshold for the change of objective functions for terminating outer iterations
#   itermax: maximum number of iterations
#   
#   Output:
#   Sig: clustering point of Sig
#   C:   inv(Sig)
#   obj: sequence of the objective function at each iteration: 
#                 obj = logdet(Sigma)+trace(Q/Sig)+sum(sum(Rho.*abs(Sig)));
#   niter: number of iterations

   p = nrow(Q);
   S = Q;
   C = solve(Sig);
   n = 1;
   Lambda = Rho*n;
   
ind.noi.all = mat.or.vec(p-1,p);
for( i in 1:p)
  {

    ind.noi = seq(1,p);
    ind.noi = ind.noi[-i];
    ind.noi.all[,i] = ind.noi; 
   
  }

  ob = -determinant(C,logarithm = TRUE)$modulus+sum(C*Q)+sum(Rho*abs(Sig));

  obj <- ob;
   
#    print(paste('iter = ',0,'ob=',ob));

for(iter in 1:itermax)
  {


    ob.old = ob;
    Sig.old = abs(Sig);


    for(i in 1:p)
      {

        ind.noi <- ind.noi.all[,i];

        C11 <- C[ind.noi,ind.noi]; C12 = C[ind.noi,i];
        S11 <- S[ind.noi,ind.noi];S12 = S[ind.noi,i];
        Sig11 <- Sig[ind.noi,ind.noi];

        invSig11 <- C11-C12%*%t(C12)/C[i,i];
        invSig11S12 <- invSig11%*%S12;

        S11invSig11 <- S11%*%invSig11;
        W1 <- invSig11%*%S11invSig11;

####### Update gam
        beta <- Sig[ind.noi,i];
        a <- n/2; b <- Lambda[i,i]/2; c <- (t(beta)%*%W1%*%beta-2*t(beta)%*%invSig11S12+S[i,i])/2;

        {if (b==0)
          {gam <- c/a;}
        else
          {gam <- (-a+sqrt(a*a+4*b*c))/(2*b);}
         }


        
####### Update beta
        A1 <- (S11invSig11/as.numeric(gam)+Lambda[i,i]*diag(p-1));
        A  <- solve(A1,Sig11); 
        B  <- diag(Sig.old[ind.noi,i]/Lambda[ind.noi,i]);        
        M <- B%*%solve(A+B,A);

        beta <- M%*%invSig11S12/as.numeric(gam);
        Sig[ind.noi,i] <- beta;
        Sig[i,ind.noi] <- beta;
        Sig[i,i] <- gam+t(beta)%*%invSig11%*%beta;


        invSig11beta <-  invSig11%*%beta;
        C[ind.noi,ind.noi] <- invSig11+invSig11beta%*%t(invSig11beta)/as.numeric(gam);
       
        C12 <- -invSig11beta/as.numeric(gam);

        C[ind.noi,i]=C12;
        C[i,ind.noi]=C12;

        C[i,i] = 1/gam
        
      }
    
      ob <- -determinant(C,logarithm = TRUE)$modulus+sum(C*Q)+sum(Rho*abs(Sig));
    obj <- c(obj,ob);
     # dif <- mean(abs(Sig-Sig.old))

     # print(paste('iter = ', iter,'ob=',ob,'ave. dif=',dif));
#    print(paste('iter = ', iter,'ob=',ob));
      if(abs(ob-ob.old)<tol)
        {break;}
    
  }


  return(list(Sig = Sig, C = C, obj = obj, niter = iter))

}



