#### Script para visualizar, describir y modelar curvas de progresos de enfermedades  ####
#### by JRE https://github.com/jrodriguez88/Model_Diseases #####

##### Cargar los Datos : t (dias) -- incidencia (%) 

rm(list = ls())
#data <- read.table("clipboard", header = T)  


#Crea datos ejemplo Paper Nutter https://www.researchgate.net/publication/286088063_Visualizing_Describing_and_Modeling_Disease_Progress_Curves_using_Epimodel 

t <- c(21,33,36,40,43)
i <- c(29,59,80,91,96)
data <- cbind(t,i)
y <- data[,2]/100                  # Convierte incidencia(%) a fraccion-y-: incidencia/100
y0 <- 0.05                         # Opcional, definir incidencia inicial (y0) para intercomparar modelos

#Cada Modelo tiene su propia ecuacion y requiere transformarse a su forma lineal

mono <- log(1/(1-y))               # Transformacion linear para el modelo Monomolecular
expo <- log(y)                     # Transformacion linear para el modelo Exponencial
logi <- log(y/(1-y))               # Transformacion linear para el modelo Logistico
gomp <- -log(-log(y))              # Transformacion linear para el modelo Gompertz
weib <- log(log(1/(1-y)))          # Transformacion linear para el modelo Weibull


# Consolida datos observados con las transformaciones lineales para cada modelo

data <- cbind(data,y, mono, expo, logi, gomp, weib)  # Crea matrix con datos y transformaciones
View(data)
data
dia<- max(data[,1])              # Ultimo dia de monitoreo, punto de referencia
ciclo <- "100"                   # El ciclo de analisis depende de cada enfermedad y el tiempo de interaccion con el cultivo


## Genera Regresiones Lineales para cada transformacion. La regresion genera una pendiente (r) y un intercepto (yo)
rmono<-lm(mono~data[,1])           # Regresion lineal para el modelo Monomolecular
rexpo<-lm(expo~data[,1])           # Regresion lineal para el modelo Exponencial
rlogi<-lm(logi~data[,1])           # Regresion lineal para el modelo Logistico
rgomp<-lm(gomp~data[,1])           # Regresion lineal para el modelo Gompertz
rweib<-lm(weib~data[,1])           # Regresion lineal para el modelo Weibull

# Calculo de las tasa de progreso de la epidemis (r) por cada modelo
rm<-rmono[[1]][[2]]
re<-rexpo[[1]][[2]] 
rl<-rlogi[[1]][[2]] 
rg<-rgomp[[1]][[2]] 
rw<-rweib[[1]][[2]] 

#Calculo del Inoculo Inicial (y0) por cada modelo
ym<-rmono[[1]][[1]]
ye<-rexpo[[1]][[1]] 
yl<-rlogi[[1]][[1]] 
yg<-rgomp[[1]][[1]] 
yw<-rweib[[1]][[1]]

#Calculo de estadisticos R? y Residuales
r1<- summary(lm(mono~t))$r.squared
r2<- summary(lm(expo~t))$r.squared
r3<- summary(lm(logi~t))$r.squared
r4<- summary(lm(gomp~t))$r.squared
r5<- summary(lm(weib~t))$r.squared

s1<- summary(lm(mono~t))$sigma
s2<- summary(lm(expo~t))$sigma
s3<- summary(lm(logi~t))$sigma
s4<- summary(lm(gomp~t))$sigma
s5<- summary(lm(weib~t))$sigma

 
#Genera Tabla de Estadisticas
Modelo <- c("Monomolecular","Exponencial", "Logistico", "Gompertz","Weibull")
Pendiente <- as.vector(rbind(rm, re, rl, rg, rw))
Intercepto <- as.vector(rbind(ym,ye,yl,yg,yw))
R2 <- as.vector(rbind(r1,r2,r3,r4,r5))
SD_Residual <- as.vector(rbind(s1,s2,s3,s4,s5))
estadisticas_regresion <- cbind(Modelo, Intercepto, Pendiente, R2, SD_Residual)
estadisticas_regresion
View(estadisticas_regresion)


##Curva de Tasa
#plot(t, mono)
#(x[2]-x[1])/(t[2]-t[1])
#(t[2]+t[1])/2




##### Cargar Funciones graficas para los modelos

## Modelo Exponencial 
plotexp <- function(y0,r,maxt){
    curve(
        y0*exp(r*x),
        from=0,
        to=maxt,
        xlab='Tiempo (dias)',
        ylab='Incidencia',
        main='Modelo Exponencial',
        col='mediumblue',
        ylim = c(0, 1))
    abline(v=dia, col="red")
    points(data[,1],data[,3], col="red")
}


## Modelo Monomolecular  
plotmono <- function(y0,r,maxt){
    curve(
        1-(1-y0)*exp(-r*x),
        from=0,
        to=maxt,
        xlab='Tiempo (dias)',
        ylab='Incidencia',
        main='Modelo Monomolecular',
        col='mediumblue',
        ylim = c(0, 1))
    abline(v=dia, col="red")
    points(data[,1],data[,3], col="red")
}


## Modelo Logistico 
plotlog <- function(y0,r,maxt){
    curve(
        1/(1+(1-y0)/y0*exp(-r*x)),
        from=0,
        to=maxt,
        xlab='Tiempo (dias)',
        ylab='Incidencia',
        main='Modelo Logistic',
        col='mediumblue',
        ylim = c(0, 1))
    abline(v=dia, col="red")
    points(data[,1],data[,3], col="red")
}


## Modelo Gompertz 
plotgomp <- function(y0,r,maxt){
    curve(
        exp(log(y0)*exp(-r*x)),
        from=0, to=maxt, 
        xlab='Tiempo (dias)',
        ylab='Incidencia',
        main='Modelo Gompertz',
        col='mediumblue',
        ylim = c(0, 1))
    abline(v=dia, col="red")
    points(data[,1],data[,3], col="red")
}


## Modelo Weibull
##   a - Unidades de Tiempo, indica el dia de inicio de la epidemis, a=1
##   b - Parametro de escala. Es inversamente proporcional a "r", b=1/rw
##   c - Parametro que controla la forma de la curva. Es variable, si c=1 -> Weibull=Monomolecular
plotweib <- function(a,b,c,maxt){
    curve(
        1-exp(-((x-a)/b)^c),
        from=0,
        to=maxt,
        xlab='Tiempo (dias)',
        ylab='Incidencia',
        main='Modelo Weibull',
        col='mediumblue',
        ylim = c(0, 1))
    abline(v=dia, col="red")
    points(data[,1],data[,3], col="red")
}



########### GENERACION DE GRAFICOS

### Graficos de Progreso de la Enfermedad

#Grafico con Incidencia Inicial de 5% , y0=0.05 // Permite comparacion de Modelos

par(oma=c(0,0,2,0), mfrow=c(2,3))
plot(t,y, xlab = "Tiempo(dias)", ylab = "Incidencia", main = "Incidencia Observada") 
lines(t,y, col="red")
plotexp(y0, re, ciclo)
plotmono(y0, rm, ciclo)
plotlog(y0, rl, ciclo)
plotgomp(y0,rg, ciclo)
plotweib(1, 1/rw, 1, ciclo)
mtext(paste0("Modelos de Progreso de la Epidemia, y0=", y0*100, "%"), outer = TRUE, cex = 1.5)


#Grafico con Incidencia Inicial derivada de Regresion Lineal= Intercepto

#par(oma=c(0,0,2,0), mfrow=c(2,3))
#plot(t,y, xlab = "Tiempo(dias)", ylab = "Incidencia", main = "Incidencia Observada") 
#lines(t,y, col="red")
#plotexp(ye, re, ciclo)
#plotmono(ym, rm, ciclo)
#plotlog(yl, rl, ciclo)
#plotgomp(yg,rg, ciclo)
#plotweib(1, 1/rw, 1, ciclo)
#mtext("Modelos de Progreso de la Epidemia, y0=Intercepto de Regresion", outer = TRUE, cex = 1.5)

