# Libreria ----------------------------------------------------------------

require(gamlss)

# Lectura de datos --------------------------------------------------------

datos <- read.csv2("SCD.csv", header=T)
datos1 <- datos[-c(6,10),]

# Modelo de referencia ----------------------------------------------------

modgs <- lm(GS~TE+C+N+I(TE^2)+TE*C+TE*N+I(C^2)+C*N, data=datos1)
# Resultados del modelo
summary(modgs)
# Medidadas de ajuste
AIC(modgs)
cor(datos1$GS, fitted(modgs))
plot(modgs)


# Análisis de posibles distribuciones -------------------------------------

# Ajuste de distribuciones
mgs <- fitDist(GS,type= "realplus", data=datos1, k=log(16))
mgs$fits

histDist(datos1$GS, family=EXP, main="Distribución EXP", xlab="",ylab="Densidad")
histDist(datos1$GS, family=WEI3, main="Distribución WEI3", xlab="",ylab="Densidad")
histDist(datos1$GS, family=GA, main="Distribución GA", xlab="",ylab="Densidad")
histDist(datos1$GS, family=GG, main="Distribución GG", xlab="",ylab="Densidad")


# Modelos alternativos  ---------------------------------------------------

# Modelos alternativo 1 ---------------------------------------------------
mod1 <- lm(sqrt(GS)~TE+C+N+I(TE^2)+TE*N+C*N, data=datos1)
summary(mod1)
AIC(mod1)
cor(datos1$GS, fitted(mod1))
plot(mod1)

# Modelo alternativo 2 ----------------------------------------------------

mod2 <- gamlss(GS~TE+C+N, family=EXP, data=datos1, k=log(16))

mod2.2 <- stepGAICAll.A(mod2, scope=list(lower=~1, upper=~TE+N+C+I(TE^2)+I(C^2)+I(N^2)
                                         +TE*N*C), data=datos1)
summary(mod2.2)
GAIC(mod2.2)
cor(datos1$GS, fitted(mod2.2))
wp(mod2.2)


# Modelo alternativo 3 ----------------------------------------------------

mod3 <- gamlss(GS~TE+C+N, family=GA, data=datos1, k=log(16))

mod3.3 <- stepGAICAll.A(mod3, scope=list(lower=~1, upper=~TE+N+C+I(TE^2)+I(C^2)+I(N^2)
                                         +TE*N*C), data=datos1)
summary(mod3.3)
GAIC(mod3.3)
cor(datos1$GS, fitted(mod3.3))
wp(mod3.3)


# Modelo alternativo 4 ----------------------------------------------------

mod4 <- gamlss(GS~TE+C+N, family=GG, data=datos1, k=log(16))

mod4.4 <- stepGAICAll.A(mod4, scope=list(lower=~1, upper=~TE+N+C+I(TE^2)+I(C^2)+I(N^2)
                                         +TE*N*C), data=datos1)
summary(mod4.4)
GAIC(mod4.4)
cor(datos1$GS, fitted(mod4.4))
wp(mod4.4)


# Modelo alternativo 5 ----------------------------------------------------

mod5 <- gamlss(GS~1, family=WEI3, data=datos1, k=log(16))
mod5.5 <- stepGAICAll.A(mod5, scope=list(lower=~1, upper=~TE+N+C+I(TE^2)+I(C^2)+I(N^2)
                                         +TE*N*C), data=datos1)
summary(mod5.5)
GAIC(mod5.5)
cor(datos1$GS, fitted(mod5.5))
wp(mod5.5)


# Superficie de respuesta -------------------------------------------------

funcion<-function(x){
  TE<-x[1]
  N<-x[2]
  C<-x[3]
  -exp(predict(object=mod5.5,
               newdata=data.frame(TE=TE,N=N,C=C)))
}

optim(par=c(50,4,8),fn=funcion,method="L-BFGS-B",
      lower=c(45,3,6),upper=c(61,7,10))

# fijamos TE=45
TEF <- 45
N <- seq(from=3, to=7, length.out=20)
C <- seq(from=6, to=10, length.out=20)

Rend <- function(N, C) {
  res <- coef(mod5.5) * c(1, N, TEF^2, TEF, C, C^2,TEF*N,N*C,TEF*C)
  exp(sum(res))
}


Rend <- Vectorize(Rend)
GS1 <-outer(N, C, Rend)

persp(x=N, y=C, z=GS1,
      theta=125, phi=30, ticktype = "detailed", col='light blue',
      zlab="Grado de secado")
