# removed debugger in spcov package so that simulations don't hang
GGDescent <- function (Sigma, Sigma0, S, lambda, del, nsteps, step.size, nesterov = FALSE, 
    backtracking = FALSE, tol = 0.001, trace = 0) 
{
    if (backtracking) {
        beta <- backtracking
        if (beta <= 0 | beta >= 1) 
            stop("Backtracking parameter beta must be in (0,1).")
    }
    tt <- step.size
    converged <- FALSE
    exit <- FALSE
    obj.starting <- spcov:::ComputeObjective(Sigma, S, lambda)
    Sigma.starting <- Sigma
    Omega <- Sigma
    Sigma.last <- Sigma
    ttts <- NULL
    ttt <- tt
    for (i in seq(nsteps)) {
        inv.Sigma0 <- solve(Sigma0)
        log.det.Sigma0 <- spcov:::LogDet(Sigma0)
        grad.g <- spcov:::ComputeGradientOfg(Omega, S, Sigma0, inv.Sigma0 = inv.Sigma0)
        grad.g <- (grad.g + t(grad.g))/2
        g.omega <- spcov:::g(Omega, S, Sigma0, inv.Sigma0 = inv.Sigma0, 
            log.det.Sigma0 = log.det.Sigma0)
        while (backtracking) {
            soft.thresh <- ProxADMM(Omega - ttt * grad.g, del, 
                1, P = lambda * ttt, rho = 0.1)$X
            gen.grad.g <- (Omega - soft.thresh)/ttt
            left <- spcov:::g(soft.thresh, S, Sigma0, inv.Sigma0 = inv.Sigma0, 
                log.det.Sigma0 = log.det.Sigma0)
            right <- g.omega - ttt * sum(grad.g * gen.grad.g) + 
                ttt * sum(gen.grad.g^2)/2
            if (is.na(left) || is.na(right)) {
                print("left or right is NA.")
#                browser()
            }
            if (left <= right) {
                Sigma <- soft.thresh
                ttts <- c(ttts, ttt)
                if (mean(abs(Sigma - Sigma.last)) < tol) {
                  converged <- TRUE
                  break
                }
                if (nesterov) 
                  Omega <- Sigma + (i - 1)/(i + 2) * (Sigma - 
                    Sigma.last)
                else Omega <- Sigma
                Sigma.last <- Sigma
                if (trace > 0) 
                  cat("--true objective:", spcov:::ComputeObjective(Sigma, 
                    S, lambda), fill = T)
                if (trace > 0) 
                  cat(i, ttt, " ")
                break
            }
            ttt <- beta * ttt
            if (ttt < 1e-15) {
                cat("Step size too small: no step taken", fill = T)
                exit <- TRUE
                break
            }
        }
        if (!backtracking) {
            Sigma <- ProxADMM(Sigma - ttt * grad.g, del, 1, P = lambda * 
                ttt, rho = 0.1)$X
            if (mean(abs(Sigma - Sigma.last)) < tol) 
                converged <- TRUE
            if (nesterov) 
                Omega <- Sigma + (i - 1)/(i + 2) * (Sigma - Sigma.last)
            else Omega <- Sigma
            Sigma.last <- Sigma
        }
        if (converged) {
            if (trace > 0) {
                cat("--GG converged in", i, "steps!")
                if (backtracking) 
                  cat(" (last step size:", ttt, ")", fill = T)
                else cat(fill = T)
            }
            break
        }
        if (exit) {
            break
        }
    }
    obj.end <- spcov:::ComputeObjective(Sigma, S, lambda)
    if (obj.starting < obj.end) {
        if (nesterov) {
            cat("Objective rose with Nesterov.  Using generalized gradient instead.", 
                fill = T)
            return(GGDescent(Sigma = Sigma.starting, Sigma0 = Sigma0, 
                S = S, lambda = lambda, del = del, nsteps = nsteps, 
                step.size = step.size, nesterov = FALSE, backtracking = backtracking, 
                tol = tol, trace = trace))
        }
#        browser()
        cat("--Returning initial Sigma since GGDescent/Nesterov did not decrease objective", 
            fill = T)
        Sigma <- Sigma.starting
    }
    list(Sigma = Sigma, niter = i, step.sizes = ttts)
}
assignInNamespace('GGDescent', GGDescent, 'spcov')
