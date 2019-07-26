# Libreria ----------------------------------------------------------------

require(gamlss)

# Lectura de datos --------------------------------------------------------

datos <- read.csv2("SCD.csv", header=T)
datos1 <- datos[-c(6,10), ]
datos1$peso <- c(1,1,1,1,1,1,1,0.5,1,1,1,1,1,1,1,1)

# Análisis de posibles distribuciones -------------------------------------
mgs <- fitDist(GS,type= "realplus", data=datos1, k=log(16))
mgs$fits

# Histogramas de distribuciones
mod1<-histDist(datos1$GS, family=EXP, main="Distribución EXP", xlab="Grado de secado",ylab="Densidad")
mod2<-histDist(datos1$GS, family=WEI3, main="Distribución WEI3", xlab="Grado de secado",ylab="Densidad")
mod3<-histDist(datos1$GS, family=GA, main="Distribución GA", xlab="Grado de secado",ylab="Densidad")
mod4<-histDist(datos1$GS, family=GG, main="Distribución GG", xlab="Grado de secado",ylab="Densidad")


hist(datos1$GS, freq=FALSE, col='gray', las=1,
     ylab="Densidad", xlab="Grado de secado",
     main="Distribución EXP")
curve(dEXP(x, mu=mod1$mu), lwd=3, add=TRUE)

hist(datos1$GS, freq=FALSE, col='gray', las=1,
     ylab="Densidad", xlab="Grado de secado",
     main="Distribución WEI3")
curve(dWEI3(x, mu=mod2$mu, sigma=mod2$sigma),
      lwd=3, add=TRUE)

hist(datos1$GS, freq=FALSE, col='gray', las=1,
     ylab="Densidad", xlab="Grado de secado",
     main="Distribución GA")
curve(dGA(x, mu=mod3$mu, sigma=mod3$sigma),
      lwd=3, add=TRUE)

hist(datos1$GS, freq=FALSE, col='gray', las=1,
     ylab="Densidad", xlab="Grado de secado",
     main="Distribución GG")
curve(dGG(x, mu=mod4$mu, sigma=mod4$sigma, nu=mod4$nu),
      lwd=3, add=TRUE)

# Modelo de referencia ----------------------------------------------------

mod0 <- lm(GS ~ TE+C+N+I(TE^2)+TE*C+TE*N+I(C^2)+C*N, data=datos1)
# Resultados del modelo
summary(mod0)
# Medidadas de ajuste
AIC(mod0)
cor(datos1$GS, fitted(mod0))
plot(mod0)

# Ajustado mod0 con gamlss para obtener pseudo R2
mod0 <- gamlss(GS~TE+C+N+I(TE^2)+TE*C+TE*N+I(C^2)+C*N, data=datos1)
Rsq(mod0)

# Modelos alternativo 1 ---------------------------------------------------

mod1 <- gamlss(sqrt(GS) ~ TE+C+N+I(TE^2)+TE*N+C*N, data=datos1)
summary(mod1)
AIC(mod1)
cor(datos1$GS, fitted(mod1))
plot(mod1)
Rsq(mod1)

# Modelo alternativo 2 ----------------------------------------------------

mod2 <- gamlss(GS~TE+C+N, family=EXP, data=datos1, k=log(16), 
               weights=datos1$peso)

mod2.2 <- stepGAICAll.A(mod2, 
                        scope=list(lower=~1, 
                                   upper=~TE+N+C+I(TE^2)+I(C^2)+I(N^2)+TE*N*C), 
                        data=datos1)
summary(mod2.2)
GAIC(mod2.2)
cor(datos1$GS, fitted(mod2.2))
wp(mod2.2)
Rsq(object = mod2.2, type = "both")


# Modelo alternativo 3 ----------------------------------------------------

mod3 <- gamlss(GS~TE+C+N, family=GA, data=datos1, k=log(16), weights = datos1$peso)

mod3.3 <- stepGAICAll.A(mod3, scope=list(lower=~1, upper=~TE+N+C+I(TE^2)+I(C^2)+I(N^2)
                                         +TE*N*C), data=datos1)
summary(mod3.3)
GAIC(mod3.3)
cor(datos1$GS, fitted(mod3.3))
wp(mod3.3)
Rsq(object = mod3.3, type = "both")

# Modelo alternativo 4 ----------------------------------------------------

mod4 <- gamlss(GS~TE+C+N, family=GG, data=datos1, k=log(16), weights = datos1$peso)

mod4.4 <- stepGAICAll.A(mod4, scope=list(lower=~1, upper=~TE+N+C+I(TE^2)+I(C^2)+I(N^2)
                                         +TE*N*C), data=datos1)
summary(mod4.4)
GAIC(mod4.4)
cor(datos1$GS, fitted(mod4.4))
wp(mod4.4)
Rsq(object = mod4.4, type = "both")

# Modelo alternativo 5 ----------------------------------------------------

mod5 <- gamlss(GS~1, family=WEI3, data=datos1, k=log(16), weights = datos1$peso)
mod5.5 <- stepGAICAll.A(mod5, scope=list(lower=~1, upper=~TE+N+C+I(TE^2)+I(C^2)+I(N^2)
                                         +TE*N*C), data=datos1)
summary(mod5.5)
GAIC(mod5.5)
cor(datos1$GS, fitted(mod5.5))
wp(mod5.5)
Rsq(object = mod4.4, type = "both")

# Superficie de respuesta -------------------------------------------------

funcion <-function(x){
  TE <-x[1]
  N <-x[2]
  C <-x[3]
  -exp(predict(object=mod5.5,
               newdata=data.frame(TE=TE,N=N,C=C)))
}

optim(par=c(50,4,8),fn=funcion,method="L-BFGS-B",
      lower=c(45,3,6),upper=c(61,7,10))

#fijamos N=3
N <-3
TEF <- seq(from=45, to=61, length.out=20)
C <- seq(from=6, to=10, length.out=20)

Rend <- function(TEF, C) {
  res <- coef(mod5.5) * c(1, N, TEF^2, TEF, C, TEF*N, N*C)
  exp(sum(res))
}

Rend <- Vectorize(Rend)
GS1 <-outer(TEF, C, Rend)

persp(x=C, y=TEF, z=GS1, xlab="C", ylab="T",
      theta=135, phi=40, ticktype = "detailed", col='gray',
      zlab="Grado de secado")


