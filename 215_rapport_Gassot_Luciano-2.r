set.seed(1)

# --- EXERCICE 1 ---

rm(list = ls())

# - Question 2 -


# METHODE 1

# Simulation de la densite g avec deux lois uniformes 
rgen_g1 <- function(n) {
  x <- runif(n, - pi / 2, pi / 2)
  y <- runif(n, -1, 1)
  return(cbind(x,y))
}

# Constante m
m1 <- (5+2/sqrt(2))*2*pi*exp(pi)

rho1<- function(x,y){
  return ((abs(sin((2/pi)*x*x-pi/4))+4*cos(x)*cos(x)+y^4)*exp(-2*(x+abs(y)))/((5+sqrt(2)/2)*exp(pi)))
}

# Algo de Rejet
rgen_f1 <- function(n){     
  ans <- c() #Echantillon de sortie 
  rho.ans <- c()
  nl <- 0
  nt <- 0
  ratio <- 1
  k <- n
  while (k>0) {
    nl <- nl + 1
    nt <- nt + floor(k/ratio)
    sim_g <- rgen_g1(floor(k/ratio))
    x <- sim_g[,1]
    y <- sim_g[,2]
    rho <- rho1(x,y)
    w <- which(runif(floor(k/ratio)) <= rho)
    rho.ans <- append(rho.ans,rho)
    ans  <- rbind(ans, sim_g[w, ])
    ratio <- mean(rho)
    k <- n - nrow(ans)
  }
  return(list("Simulation1"=ans, "Rho1"=rho.ans, "nl1"=nl, "nt1"=nt))
}

# METHODE 2

# Simulation selon la densite g1 avec la methode de la fonction inverse
fct_inv_x <- function(x){
  return (-(log(1-x*(1-exp(-2*pi)))+pi)/2)
}

# Simulation selon la densite g2 avec la methode de la fonction inverse
fct_inv_y <- function(x){
  return ((log(2*x*(exp(2)-1)+1)-2)/2*(x<1/2)-(log(1-(2*x-1)*(1-exp(-2))))/2*(x>=1/2))
}


# Simulation de n realisations suivant la densite g
rgen_g <- function(n) {  
  u <- runif(n)
  x <- fct_inv_x(u)
  y <- fct_inv_y(u)
  return(cbind(x,y))
}

# Constante m 
m <- ((5+sqrt(2)/2)*(1-exp(-2))*(exp(pi)-exp(-pi)))/2    



# Fonction rho
rho<- function(x,y){    
  return ((abs(sin((2/pi)*x^2-pi/4))+4*cos(x)^2+y^4)/(5+sqrt(2)/2))
}


# Simulation de n realisations suivant la densite f par methode du rejet
rgen_f <- function(n){     
  ans <- c() #Echantillon de sortie 
  rho.ans <- c()
  nl <- 0
  nt <- 0
  ratio <- 1
  k <- n
  while (k>0) {
    nl <- nl + 1
    nt <- nt + floor(k/ratio)
    sim_g <- rgen_g(floor(k/ratio))
    x <- sim_g[,1]
    y <- sim_g[,2]
    rho <- rho(x,y)
    w <- which(runif(floor(k/ratio)) <= rho)
    rho.ans <- append(rho.ans,rho)
    ans  <- rbind(ans, sim_g[w, ])
    ratio <- mean(rho)
    k <- n - nrow(ans)
  }
  return(list("Simulation"=ans[1:n,], "Rho"=rho.ans, "nl"=nl, "nt"=nt))
}


n <- 10000


# Autoevaluation

n_test <- 100

nl <- 0
nt <- 0

for (i in 1:n_test){
  test <- rgen_f(n)
  nl <- nl + test$nl
  nt <- nt + test$nt
}

nl <- nl/n_test
nt <- nt/n_test


# - Question 3 -

# Echantillon z de densite f 
f <- rgen_f(n)
z <- f$Simulation 

# - Question 4.b) -

# Estimateur b_n et intervalle de confiance associ?
bn_estim <- function(y, level) {
  delta <- 1/mean(y)
  s2 <- var(y)*(delta^4)
  eps <- (pnorm(0.5 * (1 + level))*sqrt(s2/length(y)))
  return(data.frame(
    value = delta,
    var = s2,
    ic_inf = delta - eps, 
    ic_sup = delta + eps
  ))
}

# Variable dans l'estimateur b_n
y_b <- m*f$Rho

b_n <- bn_estim(y_b, 0.95)


# - Question 4.c) -

# Methode Bootstrap

# On g?n?re ? partir de k ?chantillons de taille 10000 k estimateurs afin 
# d'estimer l'esperance de l'estimateur

k <- 10000 #nombre d'estimateurs
v_boot <- c()
for (i in 1:k) {
  # x_boot : vecteur aleatoire de taille 1000 des valeurs de rho
  x_boot <- sample(x = m*f$Rho, size = 1000, replace = TRUE) 
  # b_n_boot : estimateur de x_boot
  b_n_boot <- bn_estim(x_boot, 0.95)
  # v_boot : vecteur de k estimateurs
  v_boot <- append(v_boot,b_n_boot$value) #vecteur de k estimateurs 
}

# Calcul du biais estim?
biais_estim  <- (mean(v_boot) - b_n$value)



# - Question 5.b) -

# Estimateur a_n et intervalle de confiance associ?
an_estim <- function(y, level) {
  delta <- mean(y)
  s2 <- var(y)
  eps <- (pnorm(0.5 * (1 + level)) * sqrt(s2 / length(y)))
  return(data.frame(
    value = delta, 
    var = s2,
    ic_inf = delta - eps, 
    ic_sup = delta + eps
  ))
}

# Variable dans l'estimateur a_n
y_a <- function(x,y){
  return((1/(m*rho(x,y))))
}

a_n <- an_estim(y_a(z[,1],z[,2]), 0.95)


# - Question 6 -

# Rapport des couts 
st_1 <-  a_n$var/b_n$var


# - Question 7 -

psi_x <- function(x){
  return(exp(-2*x)*(((exp(2)-1)/exp(2))*(abs(sin((2/pi)*x*x-(pi/4)))+4*cos(x)*cos(x))+3*(exp(2)-7)/(4*exp(2)))*( abs(x) <= (pi/2)))
}

# Estimateur f_(X,n) de la densite marginale de X 
f_x_n <- function(a, x) {
  return(a*psi_x(x))
}

# Intervalle d'?valuation
s <- seq((-pi/2), (pi/2), length.out = 10000) 




# Comparaison graphique de la distribution marginale de l'echantillon z a l'estimateur f_x_n

par(mfrow = c(1, 1))

hist(
  x = z[, 1], freq = F, main = "Histogramme de X", xlab = "", ylab = "Frequences",
  col = "grey70" , breaks = 60, ylim = c(0,1))

lines( s, f_x_n(a_n$value, s), col = "orange", lwd = 3)

legend("topright", "Estimateur de la densite marginale X", col = "orange",
       lty = c(1, 4), lwd = c(2, 3), box.lty = 0, bg = "gray95", inset = .05
)

# - Question 10 -

# Estimateur w_n et intervalle de confiance associ?
wn_estim <- function(t, y, level) {
  delta <- mean(y(t))
  s2 <- var(y(t))
  eps <- (pnorm(0.5 * (1 + level)) * sqrt(s2 / length(y(t))))
  return(data.frame(
    value = delta, 
    var = s2,
    ic_inf = delta - eps, 
    ic_sup = delta + eps
  ))
}

# Fonction psi  
psi <- function(x,y){
  return(((abs(sin((2/pi)*x*x-pi/4))+4*cos(x)*cos(x)+y^4)*exp(-2*(x+abs(y))))*(abs(x) <= (pi/2))*(abs(y) <= 1))
}

# Variable dans l'estimateur w_n 
y_w <- function(t){
  return(((psi(t, z[, 2])*f_x_n(a_n$value,z[, 1]))/(psi(z[, 1], z[, 2]))))
}

w_n <- wn_estim(-1, y_w, 0.95)

# - Question 11 - 

# Rapport des couts 
st_2 <-  w_n$var/var(y_a(z[,1],z[,2])*psi_x(-1)^2)


# --- EXERCICE 2 ---

rm(list = ls())

# - Question 1 -
rmvnorm <- function(n, mu, sigma) { # Simulation de X
  Z <- matrix(rnorm(3 * n), nrow = 3, ncol = n)
  X <- t(chol(sigma))
  return(mu + X %*% Z)
}

# On utilise les donnees de l'exercice pour simuler 10000 simulations de X
n <- 10000
mu <- c(0.1, 0, 0.1)
sigma <- matrix(c(0.047, 0, 0.0117, 0, 0.047, 0, 0.0117, 0, 0.047), 3, 3)

x <- rmvnorm(n, mu, sigma)


#  - Question 2.b) -

h <- function(x, n) { # Fonction h utilisant pmin pour prendre le min de 2 vecteurs
  return(pmin(rep(3, n), colMeans(exp(-x))))
}

h.x <- h(x, n)
delta <- mean(h.x)
delta

var.h <- var(h.x)
erreur <- var.h / n
erreur

# - Question 3. b) -

# On introduit la variable antith?tique A= 2*mu - X
a <- -x + 2 * mu

par(mfrow = c(1, 2))

q_a <- quantile(a)
q_x <- quantile(x)

qqplot(q_x, q_a, xlab = "Quantiles de X",ylab = "Quantiles de A(X)", main = "QQplot de A(X) et X", col = "blue")
abline(a = 0, b = 1)


hist(a, freq = F,main = "Repartition de l'enchantillon A(X)",xlab = "A(X)", col = "grey70")
lines(density(x), col = "red")
legend("topright","Densit? de X",col = "red",lwd = 1,box.lty = 1)

# Calcul de l'estimateur
h.a <- h(a, n)
delta.ant <- mean(c(h.x, h.a))
delta.ant

# Calcul de l'erreur quadratique moyenne
rho <- cor(h.x, h.a)

erreur.ant <- var.h * (1 + rho) / (2 * n)
erreur.ant

# Calcul du  facteur de r?duction de variance th?orique

library(microbenchmark)
test <- microbenchmark::microbenchmark(h(x, n), rmvnorm(n, mu, sigma), times = 1000)
print(test)

C_h <- mean(test$time[which(test$expr == "h(x, n)")])
C_X <- mean(test$time[which(test$expr == "rmvnorm(n, mu, sigma)")])


R1 <- 2 * (C_h + C_X) / ((C_X + 2 * C_h) * (1 + rho))
R1


# - Question 4. a) -


h0 <- function(x) {
  return(exp(-colMeans(x)))
}

h0.x <- h0(x)

rho.control <- cor(h0.x, h.x)
rho.control

t1 <- c(-1 / 3, -1 / 3, -1 / 3)
m <- exp(t(mu) %*% t1 + 1 / 2 * t(t1) %*% sigma %*% t1) [1, 1]

essai <- seq(from = 0, to = 1000, by = 1)
err <- c()
del <- c()
b <- c()
for (l in essai) {
  h.x.l <- h.x[1:l]
  h.x.nl <- h.x[l:n]
  
  h0.x.l <- h0.x[1:l]
  h0.x.nl <- h0.x[l:n]
  bl <- sum((h0.x.l - m) * (h.x.l - mean(h.x.l))) / sum((h0.x.l - m) * (h0.x.l - m))
  
  d <- mean(h.x.nl - bl * (h0.x.nl - m))
  
  e <- var(h.x.nl) / (n - l) + (bl * bl * var(h0.x) - 2 * bl * cov(h.x.nl, h0.x.nl)) / (n - l)
  
  err <- append(err, e)
  b <- append(b, bl)
  del <- append(del, d)
}

par(mfrow = c(1, 2))

&


w <- which(err == min(err))
essai[w] #Valeur de l
err[w] #Erreur minimale
b[w] #Valeur de b
del[w] #Valeur de l'estimateur


# --- EXERCICE 3 ---

rm(list = ls())

n <- 10000

m <- 2
p <- 0.2
theta <- 2

# Creation de formules adaptees a l'enonce
rgeom_modif <- function(n, p = 0.2) {
  return(rgeom(n, p) + 1)
}

dgeom_modif <- function(n, p = 0.2) {
  return(dgeom(n - 1, p))
}

pgeom_modif <- function(n, p = 0.2) {
  return(pgeom(n - 1, p))
}

estim_MC <- function(n) {
  ans <- c()
  y <- rgeom_modif(n, p)
  for (i in 1:n) {
    ans <- append(ans, sum(log(1 + rgamma(y[i], m, theta))))
  }
  return(list("delta" = mean(ans), "Var" = var(ans) / n, "Erreur" = var(ans) / n^2))
}


MC <- estim_MC(n)
MC

# Construction de l'estimateur stratifie
estim_Strat <- function(n, p = 0.2) {
  L <- 15
  # Allocation proportionnelle
  nk <- c()
  pk <- dgeom_modif(1:(L - 1), p)
  nk <- floor(n * pk)
  pk <- append(pk, 1 - pgeom_modif(L, p))
  nk <- append(nk, n - sum(nk))
  estim <- c()
  erreur <- c()
  for (i in 1:(L - 1)) { # Simulation des Si(k)
    sk <- replicate(nk[i], sum(log(1 + rgamma(i, theta, m))))
    estim <- append(estim, pk[i] * mean(sk))
    erreur <- append(erreur, pk[i] * var(sk))
  }
  U <- runif(nk[L])
  P_Y15 <- pgeom_modif(L, p)
  # Simulation de n15 Y suivant la loi (Y >= 15)
  Y_cond <- qgeom(P_Y15 + (1 - P_Y15) * U, p)
  
  sL <- c()
  # Simulation de n15 Si(15)
  for (i in Y_cond) {
    sk <- sum(log(1 + rgamma(i, theta, m)))
    sL <- append(sL, sk)
  }
  estim <- append(estim, pk[L] * mean(sL))
  erreur <- append(erreur, pk[L] * var(sL))
  sum(estim)
  sum(erreur) / n
  return(list("delta" = sum(estim), "Var" = sum(erreur) / n, "Erreur" = sum(erreur) / n^2))
}

Strat <- estim_Strat(n)
Strat

# Calcul du facteur de r?duction de variance
library(microbenchmark)
test <- microbenchmark(estim_MC(n), estim_Strat(n))
print(test)

C_MC <- mean(test$time[which(test$expr == "estim_MC(n)")])
C_Strat <- mean(test$time[which(test$expr == "estim_Strat(n)")])

R <- (MC$Var * C_MC) / (C_Strat * Strat$Var)
R
