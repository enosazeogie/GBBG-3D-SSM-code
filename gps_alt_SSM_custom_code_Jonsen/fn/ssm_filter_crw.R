##' @title fit a CRW state-space model to data - adapted from `aniMotum`
##'
##' @description generates initial values for model parameters and unobserved states;
##' structures data and initial values for C++ \code{TMB} template;
##' fits state-space model; minimizes the joint log-likelihood via the selected
##' optimizer (\code{nlminb} or \code{optim}); structures and passes output
##' object to \code{fit_ssm}
##'
##' @details called by \code{fit_ssm}, not intended for general use. \code{sfilter} can only fit to an
##' individual track, use \code{fit_ssm} to fit to multiple tracks (see ?fit_ssm).
##'
##' @param x Argos data passed through prefilter()
##' @param time.step the regular time interval, in hours, to predict to.
##' Alternatively, a vector of prediction times, possibly not regular, must be
##' specified as a data.frame with id and POSIXt dates.
##' @param omega (fixed) standard deviation representing average drone altitude error
##' @param df the fixed degrees of freedom for the t-distribution for altitude errors
##' @param verbose is deprecated, use ssm_control(verbose = 1) instead, see \code{ssm_control} for details
##'
##' @importFrom TMB MakeADFun sdreport newtonOption FreeADFun
##' @importFrom stats approx cov sd predict nlminb optim na.omit
##' @importFrom utils flush.console
##' @importFrom tibble as_tibble
##' @importFrom sf st_crs st_coordinates st_geometry<- st_as_sf st_set_crs


ssm_filter_crw <-
  function(x,
           time.step = 6,
           omega = 6.1297, ## average drone error (galt - lalt) - GPS alt uncorrected for geoid
           df = 4,
           verbose = TRUE) {
    
    if(all(c("lon","lat") %in% names(x))){
      x <- sf::st_as_sf(x, coords = c("lon","lat"), crs = 4326) |>
        sf::st_transform(crs = "+proj=merc +units=km +datum=WGS84")
      
      loc <- as.data.frame(st_coordinates(x))
      names(loc) <- c("x","y")
      # loc <- loc |> 
      #   mutate(x = x - mean(x), 
      #         y = y - mean(y))
      st_geometry(x) <- NULL
      d <- cbind(x, loc)
    } else {
      d <- x
    }
    
    d$isd <- TRUE
    
    ## generate prediction intervals
    if (!inherits(time.step, "data.frame") & all(!is.na(time.step))) {
      ## prediction times - assume on time.step-multiple of the hour
      tsp <- time.step * 3600
      tms <- (as.numeric(d$date) - as.numeric(d$date[1])) / tsp
      index <- floor(tms)
      
      if (time.step > 1) {
        ## time.step >= 1 h
        ## trunc so predictions start on the hour immediately prior to 1st obs
        ts <-
          data.frame(id = d$id[1],
                     date = seq(
                       trunc(d$date[1], "hour"),
                       by = tsp,
                       length.out = max(index) + 2
                     ))
      } else {
        ## time.step < 1 h
        ## trunc so predictions start on the time.step immediately prior to 1st obs
        ts1 <- trunc(d$date[1] - tsp, "mins") + tsp
        
        if(ts1 <= d$date[1]) {
          ts <-
            data.frame(id = d$id[1],
                       date = seq(ts1,
                                  by = tsp,
                                  length.out = max(index) + 2))
        } else {
          ts <-
            data.frame(id = d$id[1],
                       date = seq(ts1 - tsp,
                                  by = tsp,
                                  length.out = max(index) + 2))
        }
      }
      
    } else if (inherits(time.step, "data.frame") & all(!is.na(time.step))) {
      ts <- subset(time.step, id %in% unique(d$id))
      
    } else if (inherits(time.step, "data.frame") & any(is.na(time.step))) {
      stop("NA's are not allowed in user-supplied prediction times data.frame")
    }
    
    if (all(!is.na(time.step))) {
      ## add 0.5 s to observation time(s) that exactly match prediction time(s)
      if (sum(d$date %in% ts$date) > 0) {
        o.times <- which(d$date %in% ts$date)
        d[o.times, "date"] <- d[o.times, "date"] + 0.5
      }
      
      ## merge data and prediction times
      ## add is.data flag (distinguish obs from reg states)
      d.all <- full_join(d, ts, by = c("id", "date"))
      d.all <- d.all[order(d.all$date), ]
      d.all$isd <- with(d.all, ifelse(is.na(isd), FALSE, isd))
      d.all$id <- with(d.all, ifelse(is.na(id), na.omit(unique(id))[1], id))
    } else {
      d.all <- d
    }
    
    ## calc delta times in hours for observations & interpolation points (states)
    dt <- as.numeric(difftime(d.all$date, c(as.POSIXct(NA),d.all$date[-nrow(d.all)]), 
                              units = "hours"))
    dt[1] <- 0.000001 # - 0 causes numerical issues in CRW model
    
    ## use approx & MA filter to obtain state initial values
    x.init1 <-
      approx(x = select(d, date, x),
             xout = d.all$date,
             rule = 2)$y
    x.init <-
      as.numeric(stats::filter(x.init1, rep(1, 5) / 5))
    x.na <- which(is.na(x.init))
    x.init[x.na] <- x.init1[x.na]
    
    y.init1 <-
      approx(x = select(d, date, y),
             xout = d.all$date,
             rule = 2)$y
    y.init <-
      as.numeric(stats::filter(y.init1, rep(1, 5) / 5))
    y.na <- which(is.na(y.init))
    y.init[y.na] <- y.init1[y.na]
    
    z.init1 <-
      approx(x = select(d, date, z),
             xout = d.all$date,
             rule = 2)$y
    z.init <-
      as.numeric(stats::filter(z.init1, rep(1, 5) / 5))
    z.na <- which(is.na(z.init))
    z.init[z.na] <- z.init1[z.na]
    
    xs <- cbind(x.init, y.init, z.init)
    
    state0 <- c(xs[1,1], xs[1,2], 2, 0, 0, 0)
    attributes(state0) <- NULL
    
    # ## Estimate stochastic innovations
    # es <- (xs[-1,] - xs[-nrow(xs),])/dt
    # 
    # ## Estimate components of variance
    # v <- rbind(rep(0,3), es)
    
    v <- cbind(c(0,diff(x.init)), c(0,diff(y.init)), c(0,diff(z.init))) / dt
    
    parameters <- list(
      X = t(xs),
      V = t(v),
      log_sigma = c(1,1,1), 
      log_tau = c(1,1,1),
      l_rho_p = 0,
      l_rho_o = 0)
    
    ## TMB - data list
    ## convert from string to integer for C++
    ## calculate fitted & predicted indices + delta t for proper speed estimates
    
    data <- list(
      Y = rbind(d.all$x, d.all$y, d.all$z), 
      tdiff = dt,
      state0 = state0,
      isd = as.integer(d.all$isd),
      omega = omega, 
      df_z = df)
    

    ## TMB - create objective function
    inner.control <- list(smartsearch = TRUE)
    rnd <- c("X", "V")
    
    try(dyn.unload(dynlib("src/crwz")), silent = TRUE)
    suppressMessages(compile("src/crwz.cpp"))
    dyn.load(dynlib("src/crwz"))
    
    obj <-
      MakeADFun(
        data,
        parameters,
        random = rnd,
        hessian = TRUE,
        method = "L-BFGS-B",
        DLL = "crwz",
        silent = TRUE
      )
    
    #    obj$env$tracemgc <- verbose
    
    ## add par values to trace if control$verbose = TRUE
    myfn <- function(x) {
      cat("\r", "pars:  ", round(x, 5), "     ")
      flush.console()
      obj$fn(x)
    }
    
    
    
    ## Set parameter bounds - most are -Inf, Inf
    L = c(log_sigma=c(-Inf,-Inf,-Inf),
                    l_rho_p=-7,
          log_tau=c(-Inf,-Inf,-Inf),
              l_rho_o=-7)
    #          log_df_z = log(1)) 
    
    U = c(log_sigma=c(Inf,Inf,Inf),
                    l_rho_p=7,
          l_tau=c(Inf,Inf,Inf),
              l_rho_o=7)
    #          log_df_z = 3.912)
    
    #    names(L)[c(1:3,5:7)] <- c("log_sigma", "log_sigma", "log_sigma", "log_tau", "log_tau", "log_tau")
    #    names(U)[c(1:3,5:7)] <- c("log_sigma", "log_sigma", "log_sigma", "log_tau", "log_tau", "log_tau")
    
    ## Minimize objective function
    oldw <- getOption("warn")
    options(warn = -1)  ## turn warnings off but check if optimizer crashed & return warning at end
    # opt <- try(do.call(optim,
    #                             args = list(
    #                               par = obj$par,
    #                               fn = ifelse(verbose, myfn, obj$fn),
    #                               gr = obj$gr,
    #                               method = "L-BFGS-B",
    #                               lower = L,
    #                               upper = U
    #                             )
    #                           ), silent = TRUE)
    
    opt <- try(nlminb(obj$par,
                      ifelse(verbose, myfn, obj$fn),
                      #obj$fn,
                      obj$gr,
                      control = list(eval.max = 3000, iter.max = 2000, rel.tol = 0.001, x.tol = 0.015),
                      lower = L,
                      upper = U
    ))
    cat("\n")
    
    
    FreeADFun(obj)
    
    ## if error then exit with limited output to aid debugging
    ## check if pdHess is FALSE at end and return warning
    rep <- try(sdreport(obj)) #, skip.delta.method = !control$se
    
    options(warn = oldw) ## turn warnings back on
    
    ## Parameters, states and the fitted values
    fxd <- summary(rep, "report")
    if("sigma" %in% row.names(fxd)) {
      row.names(fxd)[which(row.names(fxd) == "sigma")] <- c("sigma_x","sigma_y","sigma_z")
    } 
    if("tau" %in% row.names(fxd)) {
      row.names(fxd)[which(row.names(fxd) == "tau")] <- c("tau_x","tau_y","tau_z")
    }
    
    tmp <- summary(rep, "random")
    loc <- tmp[rownames(tmp) == "X",]
    vel <- tmp[rownames(tmp) == "V",]
    
    loc <- cbind(loc[seq(1, dim(loc)[1], by = 3),],
                 loc[seq(2, dim(loc)[1], by = 3),],
                 loc[seq(3, dim(loc)[1], by = 3),]) 
    loc <- as.data.frame(loc, row.names = 1:nrow(loc))
    names(loc) <- c("x", "x.se", "y", "y.se", "z", "z.se")
    
    vel <- cbind(vel[seq(1, dim(vel)[1], by = 3),],
                 vel[seq(2, dim(vel)[1], by = 3),],
                 vel[seq(3, dim(vel)[1], by = 3),])
    vel <- as.data.frame(vel, row.names = 1:nrow(vel))
    names(vel) <- c("u", "u.se", "v", "v.se", "w", "w.se")
    
    rdm <- cbind(loc, vel)
    rdm$id <- unique(d.all$id)
    rdm$date <- d.all$date
    rdm$isd <- d.all$isd
    
    rdm <- rdm[, c("id", "date", "x", "y", "z",
                   "x.se", "y.se", "z.se", "u", "v", "w", "u.se", "v.se", "w.se", "isd")]
    
    fv <- subset(rdm, isd) |>
      select(-isd)
    fv <- as_tibble(fv)
    
    pv <- subset(rdm, !isd) |>
      select(-isd)
    pv <- as_tibble(pv)
    
    npar <- length(opt[["par"]])
    n <- nrow(fv)
    AICc <- 2 * npar + 2 * opt[["objective"]] + 2 * npar * (npar + 1)/(n - npar - 1)
    
    out <- list(
      call = call,
      predicted = pv,
      fitted = fv,
      par = fxd,
      data = x,
      isd = d.all$isd,
      inits = parameters,
      ts = time.step,
      opt = opt,
      tmb = obj,
      rep = rep,
      AICc = AICc,
      optimiser = "nlminb"
    )
    
    out
  }
