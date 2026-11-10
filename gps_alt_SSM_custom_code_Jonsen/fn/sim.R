sim <- function(n=1000, 
                sigma = c(2, 2, 1.5), 
                tau = c(10, 10, 150), 
                df = 2,
                p = 0.2) {

d <- data.frame(u=rnorm(n, 0.001, sigma[1]), 
                v=rnorm(n, 0, sigma[2]), 
                z=rnorm(n, 0.01, sigma[3]))

mu <- y1 <- matrix(0, nrow = n, ncol = 3)

## simulate process
for(i in 2:n) {
  mu[i,1] <- mu[i-1, 1] + d$u[i]
  mu[i,2] <- mu[i-1, 2] + d$v[i]
  mu[i,3] <- mu[i-1, 3] + d$z[i]
  if(mu[i,3] < 0) mu[i,3] <- 0
}

## simulate measurement error
y1[,1] <- rnorm(n, mu[,1], tau[1])
y1[,2] <- rnorm(n, mu[,2], tau[2])
y1[,3] <- rt(n, df) * tau[3] + mu[,3]


## generate time
dts <- seq(ISOdatetime(2024,07,22,00,00,00, tz = "UTC"), by = "1 mins", l = n)

d <- tibble(date = dts, 
                x = y1[,1], 
                y = y1[,2], 
                z = y1[,3], 
                xt = mu[,1], 
                yt = mu[,2], 
                zt = mu[,3])

keep <- rbinom(n, 1, p)
keep <- ifelse(keep == 1, TRUE, FALSE)

return(d[keep, ])

}