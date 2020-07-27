##VRE Outbreak in Hospital - Deterministic Compartment Model
#Replication of Wolkewitz et al. 2008

library(EpiModel)
library(ggplot2)

hai_fun <- function(t, t0, parms) {
  with(as.list(c(t0, parms)), {
    ##main ODEs
    
      #Population sizes
      n.p <- c.p + u.p
      n.s <- c.s + u.s
      n.e <- c.e + u.e
      
      #Contaminated Patients
      dc.p <- (lam.u*n.p+(lam.c-lam.u)*c.p)*phi + (beta.sp*c.s/n.s+beta.ep*c.e/n.e)*(n.p-c.p) - lam.c*c.p
      
      #Uncontaminated Patient
      du.p <- (lam.u*n.p+(lam.c-lam.u)*c.p)*(1-phi) - lam.u*(n.p-c.p) - (beta.sp*c.s/n.s+beta.ep*c.e/n.e)*(n.p-c.p)
      
      #Contaminated HCW
      dc.s <- (beta.ps*c.p/n.p+beta.es*c.e/n.e)*(n.s-c.s) - mu*c.s
      
      #Uncontaminated HCW
      du.s <- mu*c.s - (beta.ps*c.p/n.p+beta.es*c.e/n.e)*(n.s-c.s)
      
      #Contaminated Environment
      dc.e <- (beta.se*c.s/n.s+beta.pe*c.p/n.p)*(n.e-c.e) - kappa*c.e - alpha*c.e
      
      #Uncontamintated Environment
      du.e <- kappa*c.e + alpha*c.e - (beta.se*c.s/n.s+beta.pe*c.p/n.p)*(n.e-c.e)
    
    ##output
    
      list(c(dc.p, du.p, dc.s, du.s, dc.e, du.e), n.p=n.p, Prevalence=(c.p/n.p))
    
  })
}

param <- param.dcm(phi=0.1, lam.u=0.1/24, lam.c=0.05/24, mu=24/24, kappa=1/24, alpha=c(0,1/24,5/24),
                   beta.sp=0.3/24, beta.se=2/24, beta.ps=2/24, beta.pe=2/24, beta.es=2/24, beta.ep=0.3/24)

init <- init.dcm(c.p=1, u.p=19, c.s=0, u.s=5, c.e=0, u.e=100)

control <- control.dcm(nsteps=1200, new.mod=hai_fun)

mod_hai <- dcm(param, init, control)

mod_hai

mod_output <- as.data.frame(mod_hai)
mod_output$Day <- (mod_output$time+24)/24 

#Plot of Prevalence of VRE Colonization
ggplot(data=mod_output, aes(x=Day, y=Prevalence)) + geom_point(aes(color=as.factor(run))) + ggtitle("Prevalence of VRE Colonization") + 
  scale_color_manual(name="Decontamination Rate", labels=c("None", "Low Daylight", "High Daylight"), values=c("#000000", "#EBB172", "#EB9234"))
