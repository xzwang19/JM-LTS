
func_sum <- function(fit, dat, meth, random){
  
  loglik <- extract_log_lik(fit)
  
  ### calculate WAIC ###
  waic <- waic(loglik)
  waic <- waic$waic
  
  ### model fit ###
  fit <- summary(fit)$summary
  fit <- as.data.frame(fit)
  fit$param <- rownames(fit)
  setDT(fit)
  fit <- fit[,.(param, mean, `2.5%`, `97.5%`, n_eff, Rhat)]
  fit$mean <- round(fit$mean, 4)
  fit$`2.5%` <- round(fit$`2.5%`,4)
  fit$`97.5%` <- round(fit$`97.5%`, 4)
  
  fitsum <- rbind(fit[grep('mu', fit$param),],
                    fit[grep('fshift', fit$param),],
                    fit[grep('shape', fit$param),],
                    fit[grep('alpha', fit$param),],
                    fit[grep('sigma_rshift', fit$param),],
                    fit[grep('sigma_e', fit$param),],
                    fit[grep('sigma_b', fit$param),][1:2],
                    fit[grep('CORR_b', fit$param),],
                    fit[grep('^rshift', fit$param),])
  
  fitsum$meth <- meth
  fitsum$waic <- waic
  fitsum$random <- random
  return(fitsum)
  
}