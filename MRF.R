# ----- BIBLIOTEKI -----
library('XML')
library("RCurl")
library('rvest')
library('data.table')
library('xts')
library('dygraphs')
library('htmltools')
library('car')
library('FinTS')
library('fGarch')
library('tseries')
library('moments')
library('rugarch')

# ----- FUNKCJE -----
# METRYKI FUNKCJE 
ME <- function(h,s){
  me=mean(h-s, na.rm = TRUE)
  return(me)
}
MAE <- function(h,s){
  mae=mean(abs(h-s), na.rm = TRUE)
  return(mae)
}
RMSE <- function(h,s){
  rmse=sqrt(mean((h-s)^2, na.rm = TRUE))
  return(rmse)
}
AMAPE <- function(h,s){
  amape=mean((abs(h-s))/(abs(h+s)), na.rm = TRUE)
  return(amape)
}
TIC <- function(h,s){
  tic=sqrt(mean((h-s)^2))/(sqrt(mean(h^2))+sqrt(mean(s^2)))
  return(tic)
}
MMEU <- function(h,s){
  u <- vector(mode = "numeric")
  s_u <- vector(mode = "numeric")
  o <- vector(mode = "numeric")
  s_o <- vector(mode = "numeric")
  for ( i in 1:length(h)){
    if (h[i] < s[i]) {
      u <- append(u, h[i])
      s_u <- append(s_u, s[i])
    }
    else if (h[i] > s[i]) {
      o <- append(o, h[i])
      s_o <- append(s_o, s[i])
    }
  }
  mmeu=(sum(abs(o-s_o))+sum(sqrt(abs(u-s_u))))/length(h)
  return(mmeu)
}
MMEO <- function(h,s){
  u <- vector(mode = "numeric")
  s_u <- vector(mode = "numeric")
  o <- vector(mode = "numeric")
  s_o <- vector(mode = "numeric")
  for ( i in 1:length(h)){
    if (h[i] < s[i]) {
      u <- append(u, h[i])
      s_u <- append(s_u, s[i])
    }
    else if (h[i] > s[i]) {
      o <- append(o, h[i])
      s_o <- append(s_o, s[i])
    }
  }
  
  mmeo=(sum(sqrt(abs(o-s_o)))+sum(abs(u-s_u)))/length(h)
  return(mmeo)
} 


# ----- OBRÓBKA DANYCH -----
data = read.csv("data.csv")
# sprawdzenie formatu
str(data)
# zmiana formatu daty
data$date = as.Date(data$date, "%Y-%m-%d")
# selekcja okna czasowego (2010-2020) 
data = data[4233:6980,]
# reset indeksów
row.names(data) <- NULL 

# ----- ANALIZA DANYCH -----
# obiekt pomocniczy
share_xts = xts(data[, c("rate")], order.by = data$date)
# wykres WIG
dygraph(share_xts, main = "WIG") %>%
  dyRangeSelector(height = 50)
# INTERPRETACJA:
#   grupowanie wariancji

# wykres ACF
acf(data$rate, lag.max = 36, na.action = na.pass,
    ylim = c(-0.3, 0.3), 
    col = "darkblue", lwd = 7,
    main = "Wykres ACF stóp zwrotu")
pacf(data$rate, lag.max = 36, na.action = na.pass,
    ylim = c(-0.3, 0.3), 
    col = "darkblue", lwd = 7,
    main = "Wykres PACF stóp zwrotu")
# INTERPRETACJA:
#   bardzo istotne pierwsze opóźnienie, potem stablizacja i niektóre opóźnienia są istotne

# test Durbina Watsona na autokorelację reszt (wykrywanie procesu MA)
durbinWatsonTest(lm(data$rate ~ 1),
                 max.lag = 5) # 5 pierwszych opóźnień
# INTERPRETACJA:
#   H0: corr = 0 => istotne pierwsze opóźnienie -> potencjał na uwzględnienie procesu MA(1)

# Wykres ACF^2 dla stóp zwrotu (efekty ARCH)
acf(data$rate^2, lag.max = 100, na.action = na.pass,
    ylim = c(-0.3, 0.3),
    col = "darkblue", lwd = 7,
    main = "Wykres ACF kwadratów stóp zwrotu")
# INTERPRETACJA:
#   P-O-T-Ę-Ż-N-E efekty ARCH

# test ARCH
ArchTest(data$rate, lags = 11)
# ITNERPRETACJA:
#   test arch dla 11 opóźnień silnie odrzuca H0 o braku ARCH - formalne potwierdzenie ACF^2,
#   nie ma sensu przeprowadzania testu durbina watsona dla poszczególnych opóźnień, wykres i archtest
#   pokazują, że dane wykazują silne efekty arch - autoregresyjna warunkowa heteroskedastycznosć
#   (test korzysta z mnożnikóW Lagrangea')

# Histogram + gęstość dopasowanego rozkładu normalnego
hist(data$rate, prob = T, breaks = 100, main = "Histogram log-stóp zwrotu WIG20", xlab="returns", col="skyblue1")
curve(dnorm(x, mean = mean(data$rate, na.rm = T),
            sd  = sd(data$rate, na.rm = T)),
      col = "darkblue", lwd = 2, add = TRUE,
)
# INTERPRETACJA:
#   minimalnie grubsze ogony większa kurtoza w porównaniu z rozkładem normalnym o empirycznych parametrach, lewy ogon wydaje się dłuższy
#   co oznacza lewostronną asymetrię rozkładu

# statystyki opisowe
empstats <- basicStats(data$rate)
# INTERPRETACJA:
#   Wysoka kurtoza (ok 12), ujemna skośność (ok -1.1) - asymetria lewostronna (dłuższy lewy ogon rozkładu)
#   sugeruje to możliwość lepszych wyników dla modeli uwzględniających asymetrię

# formalny test Jarque-Bera na normalność rozkładu
jarque.bera.test(data$rate)
# INTERPRETACJA:
#   silnie odrzucona H0 o normalności rozkładu -> dotychczasowe wnioski są zatem poprawne

# DYGRESJA:
#   Modelujemy wariancje i sprawdzamy jakość jej prognoz. Równanie wariancji zawiera w sobie parametr dla
#   reszty z równania średniej -> uwzględnienie procesu MA(q) może poprawić modelowanie


# ---- MODELOWANIE I PREDYKCJA JEDNODNIOWA ----
# przygotowanie danych do modelowania
data$obs = 1:length(data$rate)
start <- data$obs[data$date == as.Date("2020-01-02")] # pierwszy dzien out of the sample
stop <- data$obs[data$date == as.Date("2020-12-30")] # ostatni dzien out of the sample
dta = data[data$obs <= (start-1), ]

# specyfikacja
# *** MA(1)-GARCH(1,1) e~Norm***
spec.garch <- ugarchspec(variance.model = list(model = "sGARCH",
                                               garchOrder = c(1, 1)),
                         mean.model = list(armaOrder = c(0, 1),
                                           include.mean = F),
                         distribution.model = "norm")
# estymacja
garch11ma1 <- ugarchfit(spec = spec.garch, data = dta$rate)
garch11ma1

# INTERPRETACJA WYDRUKU:
# test Ljunga-Boxa dla standaryzowanych reszt (H0: brak korelacji - przyjęta zawsze) =>
#     autokorelacja reszt została wyjaśniona

# test Ljunga-Boxa dla kwadratów standaryzowanych reszt (H0: brak korelacji - przyjęta zawsze) =>
#     autokorelacja kwadratów reszt została wyjaśniona

# test ważony ARCH LM (H0: brak efektów ARCH - przyjęta zawsze) =>
#     efekt ARCH został wyjaśniony

# tutaj są wydruki do pojedynczych estymacji, tak żeby popatrzeć w parametry
plot(garch11ma1@fit$residuals, type="l")

# *** GARCH(1,1) e~Norm ***
spec.garch <- ugarchspec(variance.model = list(model = "sGARCH",
                                               garchOrder = c(1, 1)),
                         mean.model = list(armaOrder = c(0, 0),
                                           include.mean = F),
                         distribution.model = "norm")
# estymacja
garch11 <- ugarchfit(spec = spec.garch, data = dta$rate)
garch11

# Testy Ljunga-Boxa na biały szum
#    zmienna
Box.test(data$rate, type = "Ljung-Box", lag = 25)
#    reszty z modelu MA(1)-GARCH(1,1)
Box.test(garch11ma1@fit$residuals, type = "Ljung-Box", lag = 25)
#    reszty z modelu GARCH(1,1)
Box.test(garch11@fit$residuals, type = "Ljung-Box", lag = 25)
# INTERPRETACJA:
#    analizowana zmienna nie jest białym szumem, tak samo reszty z modelu GARCH(1,1) e~N
#    ale reszty z modelu MA(1)-GARCH(1,1), e~N są (H0: biały szum)
#    Uwzględnienie komponentu MA(1) w modelu GARCH jest dobrym pomysłem
{
sigma_forecast.garch11 <- ugarchforecast(garch11, n.ahead = 1)
# technikalia - zapis wartosci predykcji wariancji jako numeric
# GARCH
ME(sigma_forecast.garch11@forecast$sigmaFor[1,1]^2, data[data$obs == start,]$rate^2)
plot(dta$date,residuals_std, type = "l")
plot(garch11, which = 11)
plot(garch11, which = 10)


# *** GARCH(1,1) e~snorm ***
spec.garch <- ugarchspec(variance.model = list(model="sGARCH",
                                               garchOrder = c(1,1)),
                         mean.model = list(armaOrder=c(0,0),
                                           include.mean = F),
                         distribution.model = "snorm")
# estymacja
garch11.snorm <- ugarchfit(spec = spec.garch, data = dta$rate)
garch11.snorm
# predykcja
sigma_forecast.garch11.snorm <- ugarchforecast(garch11.snorm, n.ahead = 1)
ME(sigma_forecast.garch11.snorm@forecast$sigmaFor[1,1]^2, data[data$obs == start,]$rate^2)


# *** GARCH(1,1) e~t-student ***
spec.garch <- ugarchspec(variance.model = list(model="sGARCH",
                                               garchOrder = c(1,1)),
                         mean.model = list(armaOrder=c(1,1),
                                           include.mean = F),
                         distribution.model = "std")
# estymacja
garch11.std <- ugarchfit(spec = spec.garch, data = dta$rate)
garch11.std
# predykcja
sigma_forecast.garch11.std <- ugarchforecast(garch11.std, n.ahead = 1)
ME(sigma_forecast.garch11.std@forecast$sigmaFor[1,1]^2, data[data$obs == start,]$rate^2)
# mniejszy błąd niż GARCH(1,1) norm, wyraźnie niższe też kryteria informacyjne, reszta bez zmian
plot(dta$date,garch11.std@fit$residuals, type = "l")
plot(garch11.std, which = 11)
plot(garch11.std, which = 10)


# *** GARCH(1,1) e~t-student skośny ***
spec.garch <- ugarchspec(variance.model = list(model="sGARCH",
                                               garchOrder = c(1,1)),
                         mean.model = list(armaOrder=c(0,1),
                                           include.mean = F),
                         distribution.model = "sstd")
# estymacja
garch11.sstd <- ugarchfit(spec = spec.garch, data = dta$rate)
garch11.sstd
# predykcja
sigma_forecast.garch11.sstd <- ugarchforecast(garch11.sstd, n.ahead = 1)
ME(sigma_forecast.garch11.sstd@forecast$sigmaFor[1,1]^2, data[data$obs == start,]$rate^2)
# mniejszy błąd niż GARCH(1,1) norm, wyraźnie niższe też kryteria informacyjne, reszta bez zmian
plot(dta$date,garch11.sstd@fit$residuals, type = "l")
plot(garch11.sstd, which = 11)
plot(garch11.sstd, which = 10)


# *** GARCH(1,1) e~GED ***
spec.garch <- ugarchspec(variance.model = list(model="sGARCH",
                                               garchOrder = c(1,1)),
                         mean.model = list(armaOrder=c(0,0),
                                           include.mean = F),
                         distribution.model = "ged")
# estymacja
garch11.ged <- ugarchfit(spec = spec.garch, data = dta$rate)
garch11.ged
# predykcja
sigma_forecast.garch11.ged <- ugarchforecast(garch11.ged, n.ahead = 1)
ME(sigma_forecast.garch11.ged@forecast$sigmaFor[1,1]^2, data[data$obs == start,]$rate^2)
# mniejszy błąd niż GARCH(1,1) norm, większy niż std, kryteria wyższe niż std, reszta bez zmian

# *** GARCH(1,1) e~sGED ***
spec.garch <- ugarchspec(variance.model = list(model="sGARCH",
                                               garchOrder = c(1,1)),
                         mean.model = list(armaOrder=c(0,0),
                                           include.mean = F),
                         distribution.model = "sged")
# estymacja
garch11.sged <- ugarchfit(spec = spec.garch, data = dta$rate)
garch11.sged
# predykcja
sigma_forecast.garch11.sged <- ugarchforecast(garch11.sged, n.ahead = 1)
ME(sigma_forecast.garch11.sged@forecast$sigmaFor[1,1]^2, data[data$obs == start,]$rate^2)

# *** E-GARCH(1,1) e~norm ***
spec.garch <- ugarchspec(variance.model = list(model = "eGARCH",
                                               garchOrder = c(1, 1)),
                         mean.model = list(armaOrder = c(0, 0),
                                           include.mean = F),
                         distribution.model = "norm")
# estymacja
egarch11.norm <- ugarchfit(spec = spec.garch, data = dta$rate)
egarch11.norm
# predykcja
sigma_forecast.egarch11.norm <- ugarchforecast(egarch11.norm, n.ahead = 1)
ME(sigma_forecast.egarch11.norm@forecast$sigmaFor[1,1]^2, data[data$obs == start,]$rate^2)
# alpha < 0 -> nie wiem czy to problem,
# gamma > 0 -> brak asymetrii, wyższe kryteria nawet niż GARCH(1,1) norm, reszty niewyjaśnione, kwadraty reszt
# wyjaśnione, choć mniej pewnie -> ten model na razie jest najgorszy

# *** E-GARCH(1,1) e~std ***
spec.garch <- ugarchspec(variance.model = list(model = "eGARCH",
                                               garchOrder = c(1, 1)),
                         mean.model = list(armaOrder = c(0, 0),
                                           include.mean = F),
                         distribution.model = "std")
# estymacja
egarch11.std <- ugarchfit(spec = spec.garch, data = dta$rate)
egarch11.std
# predykcja
sigma_forecast.egarch11.std <- ugarchforecast(egarch11.std, n.ahead = 1)
ME(sigma_forecast.egarch11.std@forecast$sigmaFor[1,1]^2, data[data$obs == start,]$rate^2)
# alpha < 0 -> nie wiem czy to problem,
# gamma > 0 -> brak asymetrii, najniższe narazie kryteria, reszty niewyjaśnione, kwadraty reszt
# wyjaśnione, choć mniej pewnie -> predykcje lepsze niż egarch(1,1) norm, ale gorsze niż garch, kiepski model

# *** E-GARCH(1,1) e~std ***
spec.garch <- ugarchspec(variance.model = list(model = "eGARCH",
                                               garchOrder = c(1, 1)),
                         mean.model = list(armaOrder = c(0, 0),
                                           include.mean = F),
                         distribution.model = "sstd")
# estymacja
egarch11.sstd <- ugarchfit(spec = spec.garch, data = dta$rate)
egarch11.sstd
# predykcja
sigma_forecast.egarch11.sstd <- ugarchforecast(egarch11.sstd, n.ahead = 1)
ME(sigma_forecast.egarch11.sstd@forecast$sigmaFor[1,1]^2, data[data$obs == start,]$rate^2)


# *** E-GARCH(1,1) e~ged ***
spec.garch <- ugarchspec(variance.model = list(model = "eGARCH",
                                               garchOrder = c(1, 1)),
                         mean.model = list(armaOrder = c(0, 0),
                                           include.mean = F),
                         distribution.model = "ged")
# estymacja
egarch11.ged <- ugarchfit(spec = spec.garch, data = dta$rate)
egarch11.ged
# predykcja
sigma_forecast.egarch11.ged <- ugarchforecast(egarch11.ged, n.ahead = 1)
ME(sigma_forecast.egarch11.ged@forecast$sigmaFor[1,1]^2, data[data$obs == start,]$rate^2)
# alpha < 0 -> nie wiem czy to problem,
# gamma > 0 -> brak asymetrii, kryteria wyższe niż egarch.std, niższe niż garch11.std, reszty niewyjaśnione, kwadraty reszt
# wyjaśnione, choć mniej pewnie -> predykcje lepsze niż egarch(1,1) norm, ale gorsze niż garch, kiepski model


# *** T-GARCH(1,1) e~norm ***
spec.garch <- ugarchspec(variance.model = list(model = "fGARCH",
                                               submodel = "TGARCH",
                                               garchOrder = c(1, 1)),
                         mean.model = list(armaOrder = c(0, 0),
                                           include.mean = F),
                         distribution.model = "norm")
# estymacja
tgarch11.norm <- ugarchfit(spec = spec.garch, data = dta$rate)
tgarch11.norm
# predykcja
sigma_forecast.tgarch11.norm <- ugarchforecast(tgarch11.norm, n.ahead = 1)
ME(sigma_forecast.tgarch11.norm@forecast$sigmaFor[1,1]^2, data[data$obs == start,]$rate^2)
# ogólnie parametry przynajmniej są spoko, ale jakość modelu gorsza niż zwykły garch 

# *** T-GARCH(1,1) e~std ***
spec.garch <- ugarchspec(variance.model = list(model = "fGARCH",
                                               submodel = "TGARCH",
                                               garchOrder = c(1, 1)),
                         mean.model = list(armaOrder = c(0, 0),
                                           include.mean = F),
                         distribution.model = "std")
# estymacja
tgarch11.std <- ugarchfit(spec = spec.garch, data = dta$rate)
tgarch11.std
# predykcja
sigma_forecast.tgarch11.std <- ugarchforecast(tgarch11.std, n.ahead = 1)
ME(sigma_forecast.tgarch11.std@forecast$sigmaFor[1,1]^2, data[data$obs == start,]$rate^2)
# ogólnie parametry przynajmniej są spoko, parametr Tresholdu jest istotny czyli ta asymetria jednak jest XD,
# ale jakość modelu gorsza niż zwykły garch , kryteria też są jedne z niższych

# *** T-GARCH(1,1) e~sstd ***
spec.garch <- ugarchspec(variance.model = list(model = "fGARCH",
                                               submodel = "TGARCH",
                                               garchOrder = c(1, 1)),
                         mean.model = list(armaOrder = c(0, 0),
                                           include.mean = F),
                         distribution.model = "sstd")
# estymacja
tgarch11.sstd <- ugarchfit(spec = spec.garch, data = dta$rate)
tgarch11.sstd
# predykcja
sigma_forecast.tgarch11.sstd <- ugarchforecast(tgarch11.sstd, n.ahead = 1)
ME(sigma_forecast.tgarch11.sstd@forecast$sigmaFor[1,1]^2, data[data$obs == start,]$rate^2)

# *** T-GARCH(1,1) e~ged ***
spec.garch <- ugarchspec(variance.model = list(model = "fGARCH",
                                               submodel = "TGARCH",
                                               garchOrder = c(1, 1)),
                         mean.model = list(armaOrder = c(0, 0),
                                           include.mean = F),
                         distribution.model = "ged")
# estymacja
tgarch11.ged <- ugarchfit(spec = spec.garch, data = dta$rate)
tgarch11.ged
# predykcja
sigma_forecast.tgarch11.ged <- ugarchforecast(tgarch11.ged, n.ahead = 1)
ME(sigma_forecast.tgarch11.ged@forecast$sigmaFor[1,1]^2, data[data$obs == start,]$rate^2)
# ogólnie parametry przynajmniej są spoko, parametr Tresholdu jest istotny czyli ta asymetria jednak jest XD,
# ale jakość modelu gorsza niż zwykły garch , kryteria gorsze niż dla std

# *** C-GARCH(1,1) e~norm ***
spec.garch <- ugarchspec(variance.model = list(model = "csGARCH",
                                               garchOrder = c(1, 1)),
                         mean.model = list(armaOrder = c(0, 0),
                                           include.mean = F),
                         distribution.model = "norm")
# estymacja
cgarch11.norm <- ugarchfit(spec = spec.garch, data = dta$rate)
cgarch11.norm
# predykcja
sigma_forecast.cgarch11.norm <- ugarchforecast(cgarch11.norm, n.ahead = 1)
ME(sigma_forecast.cgarch11.norm@forecast$sigmaFor[1,1]^2, data[data$obs == start,]$rate^2)
# kryteria kiepsko, predykcje kiepsko, na parametrach się nie znam, model lipa

# *** C-GARCH(1,1) e~std ***
spec.garch <- ugarchspec(variance.model = list(model = "csGARCH",
                                               garchOrder = c(1, 1)),
                         mean.model = list(armaOrder = c(0, 0),
                                           include.mean = F),
                         distribution.model = "std")
# estymacja
cgarch11.std <- ugarchfit(spec = spec.garch, data = dta$rate)
cgarch11.std
# predykcja
sigma_forecast.cgarch11.std <- ugarchforecast(cgarch11.std, n.ahead = 1)
ME(sigma_forecast.cgarch11.std@forecast$sigmaFor[1,1]^2, data[data$obs == start,]$rate^2)
# lepszy niż wersja z rozkładem normalnym

# *** C-GARCH(1,1) e~ged ***
spec.garch <- ugarchspec(variance.model = list(model = "csGARCH",
                                               garchOrder = c(1, 1)),
                         mean.model = list(armaOrder = c(0, 0),
                                           include.mean = F),
                         distribution.model = "ged")
# estymacja
cgarch11.ged <- ugarchfit(spec = spec.garch, data = dta$rate)
cgarch11.ged
# predykcja
sigma_forecast.cgarch11.ged <- ugarchforecast(cgarch11.ged, n.ahead = 1)
ME(sigma_forecast.cgarch11.ged@forecast$sigmaFor[1,1]^2, data[data$obs == start,]$rate^2)
# o lepsza predykcja niż e garch, prawie tak dobra jak garch ged
}



# ---- PĘTLA PREDYKCYJNA ----
# INSTRUKCJA:
#   dodając nowy model należy:
#     1. utworzyć wektor do przechowywania predykcji
#     2. utworzyć specyfikację
#     3. oszacować model w pętli
#     4. zrobić predykcje w pętli
#     5. zapisać wyniki predykcji do wektora wpętli
#     6. zapisać wektor predykcji (podniesiony do kwadratu) jako kolumnę data2

data2 <- data[start:stop, ] # df do zapisywania predykcji wariancji

# pętla była już wykonana, wyniki można wczytać
data2 <- read.csv("GarchOOSresults.csv")
data2 <- data2[,2:dim(data2)[2]]

# wektory do przechowywania rezultatów (każdy model musi mieć swój wektor)
sigma_preds.garch11ma1.norm <- rep(NA, times = stop-start+1) # MA(1)-GARCH(1,1) normalny
sigma_preds.garch11ma1.std <- rep(NA, times = stop-start+1) # MA(1)-GARCH(1,1) t-student
sigma_preds.garch11ma1.ged <- rep(NA, times = stop-start+1) # MA(1)-GARCH(1,1) ged
sigma_preds.garch11ma1.sstd <- rep(NA, times = stop-start+1) # MA(1)-GARCH(1,1) t-student skośny

sigma_preds.egarch11ma1.norm <- rep(NA, times = stop-start+1) # MA(1)-EGARCH(1,1) normalny
sigma_preds.egarch11ma1.std <- rep(NA, times = stop-start+1) # MA(1)-EGARCH(1,1) t-student
sigma_preds.egarch11ma1.ged <- rep(NA, times = stop-start+1) # MA(1)-EGARCH(1,1) ged
sigma_preds.egarch11ma1.sstd <- rep(NA, times = stop-start+1) # MA(1)-EGARCH(1,1) t-student skośny

sigma_preds.tgarch11ma1.norm <- rep(NA, times = stop-start+1) # MA(1)-TGARCH(1,1) normalny
sigma_preds.tgarch11ma1.std <- rep(NA, times = stop-start+1) # MA(1)-TGARCH(1,1) t-student
sigma_preds.tgarch11ma1.ged <- rep(NA, times = stop-start+1) # MA(1)-TGARCH(1,1) ged
sigma_preds.tgarch11ma1.sstd <- rep(NA, times = stop-start+1) # MA(1)-TGARCH(1,1) t-student skośny

sigma_preds.cgarch11ma1.norm <- rep(NA, times = stop-start+1) # MA(1)-CGARCH(1,1) normalny
sigma_preds.cgarch11ma1.std <- rep(NA, times = stop-start+1) # MA(1)-CGARCH(1,1) t-student
sigma_preds.cgarch11ma1.ged <- rep(NA, times = stop-start+1) # MA(1)-CGARCH(1,1) ged
sigma_preds.cgarch11ma1.sstd <- rep(NA, times = stop-start+1) # MA(1)-CGARCH(1,1) t-student skośny

# specyfikacje modeli
spec.garch11ma1.norm <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                                   mean.model = list(armaOrder = c(0, 1), include.mean = F),
                                   distribution.model = "norm")
spec.garch11ma1.std <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                                  mean.model = list(armaOrder = c(0, 1), include.mean = F),
                                  distribution.model = "std")
spec.garch11ma1.ged <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                                  mean.model = list(armaOrder = c(0, 1), include.mean = F),
                                  distribution.model = "ged")
spec.garch11ma1.sstd <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                                   mean.model = list(armaOrder = c(0, 1), include.mean = F),
                                   distribution.model = "sstd")

spec.egarch11ma1.norm <- ugarchspec(variance.model = list(model = "eGARCH", garchOrder = c(1, 1)),
                                    mean.model = list(armaOrder = c(0, 1), include.mean = F),
                                    distribution.model = "norm")
spec.egarch11ma1.std <- ugarchspec(variance.model = list(model = "eGARCH", garchOrder = c(1, 1)),
                                   mean.model = list(armaOrder = c(0, 1), include.mean = F),
                                   distribution.model = "std")
spec.egarch11ma1.ged <- ugarchspec(variance.model = list(model = "eGARCH", garchOrder = c(1, 1)),
                                   mean.model = list(armaOrder = c(0, 1), include.mean = F),
                                   distribution.model = "ged")
spec.egarch11ma1.sstd <- ugarchspec(variance.model = list(model = "eGARCH", garchOrder = c(1, 1)),
                                    mean.model = list(armaOrder = c(0, 1), include.mean = F),
                                    distribution.model = "sstd")

spec.tgarch11ma1.norm <- ugarchspec(variance.model = list(model = "fGARCH", submodel = "TGARCH", garchOrder = c(1, 1)),
                                    mean.model = list(armaOrder = c(0, 1), include.mean = F),
                                    distribution.model = "norm")
spec.tgarch11ma1.std <- ugarchspec(variance.model = list(model = "fGARCH", submodel = "TGARCH", garchOrder = c(1, 1)),
                                   mean.model = list(armaOrder = c(0, 1), include.mean = F),
                                   distribution.model = "std")
spec.tgarch11ma1.ged <- ugarchspec(variance.model = list(model = "fGARCH", submodel = "TGARCH", garchOrder = c(1, 1)),
                                   mean.model = list(armaOrder = c(0, 1), include.mean = F),
                                   distribution.model = "ged")
spec.tgarch11ma1.sstd <- ugarchspec(variance.model = list(model = "fGARCH", submodel = "TGARCH", garchOrder = c(1, 1)),
                                    mean.model = list(armaOrder = c(0, 1), include.mean = F),
                                    distribution.model = "sstd")

spec.cgarch11ma1.norm <- ugarchspec(variance.model = list(model = "csGARCH", garchOrder = c(1, 1)),
                                    mean.model = list(armaOrder = c(0, 1), include.mean = F),
                                    distribution.model = "norm")
spec.cgarch11ma1.std <- ugarchspec(variance.model = list(model = "csGARCH", garchOrder = c(1, 1)),
                                   mean.model = list(armaOrder = c(0, 1), include.mean = F),
                                   distribution.model = "std")
spec.cgarch11ma1.ged <- ugarchspec(variance.model = list(model = "csGARCH", garchOrder = c(1, 1)),
                                   mean.model = list(armaOrder = c(0, 1), include.mean = F),
                                   distribution.model = "ged")
spec.cgarch11ma1.sstd <- ugarchspec(variance.model = list(model = "csGARCH", garchOrder = c(1, 1)),
                                    mean.model = list(armaOrder = c(0, 1), include.mean = F),
                                    distribution.model = "sstd")


time1 = Sys.time()
for(k in start:stop){
  # dane dla bieżącej iteracji (chcę przewidzieć bieżące k)
  tmp.data <- data[data$obs <= (k-1), ] 
  
  # estymacje modeli
  tmp.garch11ma1.norm <- ugarchfit(spec = spec.garch11ma1.norm, data = na.omit(tmp.data$rate), solver = "hybrid")
  tmp.garch11ma1.std <- ugarchfit(spec = spec.garch11ma1.std, data = na.omit(tmp.data$rate), solver = "hybrid")
  tmp.garch11ma1.ged <- ugarchfit(spec = spec.garch11ma1.ged, data = na.omit(tmp.data$rate), solver = "hybrid")
  tmp.garch11ma1.sstd <- ugarchfit(spec = spec.garch11ma1.sstd, data = na.omit(tmp.data$rate), solver = "hybrid")
  
  tmp.egarch11ma1.norm <- ugarchfit(spec = spec.egarch11ma1.norm, data = na.omit(tmp.data$rate), solver = "hybrid")
  tmp.egarch11ma1.std <- ugarchfit(spec = spec.egarch11ma1.std, data = na.omit(tmp.data$rate), solver = "hybrid")
  tmp.egarch11ma1.ged <- ugarchfit(spec = spec.egarch11ma1.ged, data = na.omit(tmp.data$rate), solver = "hybrid")
  tmp.egarch11ma1.sstd <- ugarchfit(spec = spec.egarch11ma1.sstd, data = na.omit(tmp.data$rate), solver = "hybrid")
  
  tmp.tgarch11ma1.norm <- ugarchfit(spec = spec.tgarch11ma1.norm, data = na.omit(tmp.data$rate), solver = "hybrid")
  tmp.tgarch11ma1.std <- ugarchfit(spec = spec.tgarch11ma1.std, data = na.omit(tmp.data$rate), solver = "hybrid")
  tmp.tgarch11ma1.ged <- ugarchfit(spec = spec.tgarch11ma1.ged, data = na.omit(tmp.data$rate), solver = "hybrid")
  tmp.tgarch11ma1.sstd <- ugarchfit(spec = spec.tgarch11ma1.sstd, data = na.omit(tmp.data$rate), solver = "hybrid")
  
  tmp.cgarch11ma1.norm <- ugarchfit(spec = spec.cgarch11ma1.norm, data = na.omit(tmp.data$rate), solver = "hybrid")
  tmp.cgarch11ma1.std <- ugarchfit(spec = spec.cgarch11ma1.std, data = na.omit(tmp.data$rate), solver = "hybrid")
  tmp.cgarch11ma1.ged <- ugarchfit(spec = spec.cgarch11ma1.ged, data = na.omit(tmp.data$rate), solver = "hybrid")
  tmp.cgarch11ma1.sstd <- ugarchfit(spec = spec.cgarch11ma1.sstd, data = na.omit(tmp.data$rate), solver = "hybrid")
  
  # predykcje modeli (sigma - odchylenie standardowe)
  sigma_forecast.garch11ma1.norm <- ugarchforecast(tmp.garch11ma1.norm, n.ahead = 1)
  sigma_forecast.garch11ma1.std <- ugarchforecast(tmp.garch11ma1.std, n.ahead = 1)
  sigma_forecast.garch11ma1.ged <- ugarchforecast(tmp.garch11ma1.ged, n.ahead = 1)
  sigma_forecast.garch11ma1.sstd <- ugarchforecast(tmp.garch11ma1.sstd, n.ahead = 1)
  
  sigma_forecast.egarch11ma1.norm <- ugarchforecast(tmp.egarch11ma1.norm, n.ahead = 1)
  sigma_forecast.egarch11ma1.std <- ugarchforecast(tmp.egarch11ma1.std, n.ahead = 1)
  sigma_forecast.egarch11ma1.ged <- ugarchforecast(tmp.egarch11ma1.ged, n.ahead = 1)
  sigma_forecast.egarch11ma1.sstd <- ugarchforecast(tmp.egarch11ma1.sstd, n.ahead = 1)
  
  sigma_forecast.tgarch11ma1.norm <- ugarchforecast(tmp.tgarch11ma1.norm, n.ahead = 1)
  sigma_forecast.tgarch11ma1.std <- ugarchforecast(tmp.tgarch11ma1.std, n.ahead = 1)
  sigma_forecast.tgarch11ma1.ged <- ugarchforecast(tmp.tgarch11ma1.ged, n.ahead = 1)
  sigma_forecast.tgarch11ma1.sstd <- ugarchforecast(tmp.tgarch11ma1.sstd, n.ahead = 1)
  
  sigma_forecast.cgarch11ma1.norm <- ugarchforecast(tmp.cgarch11ma1.norm, n.ahead = 1)
  sigma_forecast.cgarch11ma1.std <- ugarchforecast(tmp.cgarch11ma1.std, n.ahead = 1)
  sigma_forecast.cgarch11ma1.ged <- ugarchforecast(tmp.cgarch11ma1.ged, n.ahead = 1)
  sigma_forecast.cgarch11ma1.sstd <- ugarchforecast(tmp.cgarch11ma1.sstd, n.ahead = 1)
  
  # zapis wyników do wektora
  sigma_preds.garch11ma1.norm[k-start+1] <- sigma_forecast.garch11ma1.norm@forecast$sigmaFor[1,1]
  sigma_preds.garch11ma1.std[k-start+1] <- sigma_forecast.garch11ma1.std@forecast$sigmaFor[1,1]
  sigma_preds.garch11ma1.ged[k-start+1] <- sigma_forecast.garch11ma1.ged@forecast$sigmaFor[1,1]
  sigma_preds.garch11ma1.sstd[k-start+1] <- sigma_forecast.garch11ma1.sstd@forecast$sigmaFor[1,1]
  
  sigma_preds.egarch11ma1.norm[k-start+1] <- sigma_forecast.egarch11ma1.norm@forecast$sigmaFor[1,1]
  sigma_preds.egarch11ma1.std[k-start+1] <- sigma_forecast.egarch11ma1.std@forecast$sigmaFor[1,1]
  sigma_preds.egarch11ma1.ged[k-start+1] <- sigma_forecast.egarch11ma1.ged@forecast$sigmaFor[1,1]
  sigma_preds.egarch11ma1.sstd[k-start+1] <- sigma_forecast.egarch11ma1.sstd@forecast$sigmaFor[1,1]
  
  sigma_preds.tgarch11ma1.norm[k-start+1] <- sigma_forecast.tgarch11ma1.norm@forecast$sigmaFor[1,1]
  sigma_preds.tgarch11ma1.std[k-start+1] <- sigma_forecast.tgarch11ma1.std@forecast$sigmaFor[1,1]
  sigma_preds.tgarch11ma1.ged[k-start+1] <- sigma_forecast.tgarch11ma1.ged@forecast$sigmaFor[1,1]
  sigma_preds.tgarch11ma1.sstd[k-start+1] <- sigma_forecast.tgarch11ma1.sstd@forecast$sigmaFor[1,1]
  
  sigma_preds.cgarch11ma1.norm[k-start+1] <- sigma_forecast.cgarch11ma1.norm@forecast$sigmaFor[1,1]
  sigma_preds.cgarch11ma1.std[k-start+1] <- sigma_forecast.cgarch11ma1.std@forecast$sigmaFor[1,1]
  sigma_preds.cgarch11ma1.ged[k-start+1] <- sigma_forecast.cgarch11ma1.ged@forecast$sigmaFor[1,1]
  sigma_preds.cgarch11ma1.sstd[k-start+1] <- sigma_forecast.cgarch11ma1.sstd@forecast$sigmaFor[1,1]
}


stop
time2 = Sys.time()
# czas obliczen
print(time2 - time1)
# dodanie kolumn ze zmiennością zrealizowaną (kwadrat stóp zwrotu) i przewidywaną (sigma^2)
data2$r_sq <- data2$rate^2 # zrealizowana
data2$sigma_preds_garch11ma1.norm <- sigma_preds.garch11ma1.norm^2
data2$sigma_preds_garch11ma1.std <- sigma_preds.garch11ma1.std^2
data2$sigma_preds_garch11ma1.ged <- sigma_preds.garch11ma1.ged^2
data2$sigma_preds_garch11ma1.sstd <- sigma_preds.garch11ma1.sstd^2

data2$sigma_preds_egarch11ma1.norm <- sigma_preds.egarch11ma1.norm^2
data2$sigma_preds_egarch11ma1.std <- sigma_preds.egarch11ma1.std^2
data2$sigma_preds_egarch11ma1.ged <- sigma_preds.egarch11ma1.ged^2
data2$sigma_preds_egarch11ma1.sstd <- sigma_preds.egarch11ma1.sstd^2

data2$sigma_preds_tgarch11ma1.norm <- sigma_preds.tgarch11ma1.norm^2
data2$sigma_preds_tgarch11ma1.std <- sigma_preds.tgarch11ma1.std^2
data2$sigma_preds_tgarch11ma1.ged <- sigma_preds.tgarch11ma1.ged^2
data2$sigma_preds_tgarch11ma1.sstd <- sigma_preds.tgarch11ma1.sstd^2

data2$sigma_preds_cgarch11ma1.norm <- sigma_preds.cgarch11ma1.norm^2
data2$sigma_preds_cgarch11ma1.std <- sigma_preds.cgarch11ma1.std^2
data2$sigma_preds_cgarch11ma1.ged <- sigma_preds.cgarch11ma1.ged^2
data2$sigma_preds_cgarch11ma1.sstd <- sigma_preds.cgarch11ma1.sstd^2

# zapis wyników do .csv
#write.csv(data2, "GarchOOSresults.csv")

mtr <- as.data.frame(matrix(nrow = 16, ncol = 7))
for(i in 1:dim(mtr)[1]){
  mtr[i,1] <- ME(na.omit(data2[,i+5]), na.omit(data2$r_sq))
  mtr[i,2] <- MAE(na.omit(data2[,i+5]),na.omit(data2$r_sq))
  mtr[i,3] <- AMAPE(na.omit(data2[,i+5]), na.omit(data2$r_sq))
  mtr[i,4] <- RMSE(na.omit(data2[,i+5]), na.omit(data2$r_sq))
  mtr[i,5] <- MMEO(na.omit(data2[,i+5]), na.omit(data2$r_sq))
  mtr[i,6] <- MMEU(na.omit(data2[,i+5]), na.omit(data2$r_sq))
  mtr[i,7] <- TIC(na.omit(data2[,i+5]), na.omit(data2$r_sq))
}
colnames(mtr) <- c("ME", "MAE", "AMAPE", "RMSE", "MMEO", "MMEU", "TIC")
rownames(mtr) <- c("GARCH-N", "GARCH-T", "GARCH-G", "GARCH-ST",
                   "EGARCH-N", "EGARCH-T", "EGARCH-G", "EGARCH-ST",
                   "TGARCH-N", "TGARCH-T", "TGARCH-G", "TGARCH-ST",
                   "CGARCH-N", "CGARCH-T", "CGARCH-G", "CGARCH-ST")

mtr

# *** Random walk model *** 
data2$sigma_preds_RW <- data$rate[(start-1):(stop-1)]^2
# Zapis wyników RW
mtr["RW",1] <- ME(na.omit(data2$sigma_preds_RW), na.omit(data2$r_sq))
mtr["RW",2] <- MAE(na.omit(data2$sigma_preds_RW),na.omit(data2$r_sq))
mtr["RW",3] <- AMAPE(na.omit(data2$sigma_preds_RW), na.omit(data2$r_sq))
mtr["RW",4] <- RMSE(na.omit(data2$sigma_preds_RW), na.omit(data2$r_sq))
mtr["RW",5] <- MMEO(na.omit(data2$sigma_preds_RW), na.omit(data2$r_sq))
mtr["RW",6] <- MMEU(na.omit(data2$sigma_preds_RW), na.omit(data2$r_sq))
mtr["RW",7] <- TIC(na.omit(data2$sigma_preds_RW), na.omit(data2$r_sq))

# *** Historical Average model ***
sigma_preds.HA <- rep(NA, times = stop-start+1)
for(i in start:stop){
  sigma_preds.HA[i-start+1] <- sum(data$rate[1:(i-1)]^2)/length(data$rate[1:(i-2)])
}
data2$sigma_preds_HA <- sigma_preds.HA
# Zapis wyników RW
mtr["HA",1] <- ME(na.omit(data2$sigma_preds_HA), na.omit(data2$r_sq))
mtr["HA",2] <- MAE(na.omit(data2$sigma_preds_HA),na.omit(data2$r_sq))
mtr["HA",3] <- AMAPE(na.omit(data2$sigma_preds_HA), na.omit(data2$r_sq))
mtr["HA",4] <- RMSE(na.omit(data2$sigma_preds_HA), na.omit(data2$r_sq))
mtr["HA",5] <- MMEO(na.omit(data2$sigma_preds_HA), na.omit(data2$r_sq))
mtr["HA",6] <- MMEU(na.omit(data2$sigma_preds_HA), na.omit(data2$r_sq))
mtr["HA",7] <- TIC(na.omit(data2$sigma_preds_HA), na.omit(data2$r_sq))
mtr

# *** Moving 1-year-Average model ***
# jako długość roku przyjęto liczbę obserewacji z okresu predykcyjnego (mogą wystąpić różnice, ale pomijalne)
sigma_preds.MA1 <- rep(NA, times = stop-start+1)
for(i in start:stop){
  sigma_preds.MA1[i-start+1] <- mean(data$rate[(i-dim(data2)[1]-1):(i-1)]^2)
}
data2$sigma_preds_MA1 <- sigma_preds.MA1
# Zapis wyników MA1
mtr["MA-1year",1] <- ME(na.omit(data2$sigma_preds_MA1), na.omit(data2$r_sq))
mtr["MA-1year",2] <- MAE(na.omit(data2$sigma_preds_MA1),na.omit(data2$r_sq))
mtr["MA-1year",3] <- AMAPE(na.omit(data2$sigma_preds_MA1), na.omit(data2$r_sq))
mtr["MA-1year",4] <- RMSE(na.omit(data2$sigma_preds_MA1), na.omit(data2$r_sq))
mtr["MA-1year",5] <- MMEO(na.omit(data2$sigma_preds_MA1), na.omit(data2$r_sq))
mtr["MA-1year",6] <- MMEU(na.omit(data2$sigma_preds_MA1), na.omit(data2$r_sq))
mtr["MA-1year",7] <- TIC(na.omit(data2$sigma_preds_MA1), na.omit(data2$r_sq))
mtr

# *** Moving 3-year-Average model ***
# jako długość roku przyjęto liczbę obserewacji z okresu predykcyjnego (mogą wystąpić różnice, ale pomijalne)
sigma_preds.MA3 <- rep(NA, times = stop-start+1)
for(i in start:stop){
  sigma_preds.MA3[i-start+1] <- mean(data$rate[(i-3*dim(data2)[1]-1):(i-1)]^2)
}
data2$sigma_preds_MA3 <- sigma_preds.MA3
# Zapis wyników MA1
mtr["MA-3year",1] <- ME(na.omit(data2$sigma_preds_MA3), na.omit(data2$r_sq))
mtr["MA-3year",2] <- MAE(na.omit(data2$sigma_preds_MA3),na.omit(data2$r_sq))
mtr["MA-3year",3] <- AMAPE(na.omit(data2$sigma_preds_MA3), na.omit(data2$r_sq))
mtr["MA-3year",4] <- RMSE(na.omit(data2$sigma_preds_MA3), na.omit(data2$r_sq))
mtr["MA-3year",5] <- MMEO(na.omit(data2$sigma_preds_MA3), na.omit(data2$r_sq))
mtr["MA-3year",6] <- MMEU(na.omit(data2$sigma_preds_MA3), na.omit(data2$r_sq))
mtr["MA-3year",7] <- TIC(na.omit(data2$sigma_preds_MA3), na.omit(data2$r_sq))
mtr

######################################################
#                DANE TYGODNIOWE                     #
######################################################

# ----- OBRÓBKA DANYCH -----
data_w <- read.csv("wig_w.csv")
data_w <- data_w[517:1091,c("Data", "Zamkniecie")]
colnames(data_w) <- c("date", "price")
data_w$date <- as.Date(data_w$date, "%Y-%m-%d") 
data_w$rate <- diff.xts(log(data_w$price))
data_w <- na.omit(data_w)
#write.csv(data_w, "dane_tyg.csv", row.names = F)
# ----- ANALIZA DANYCH -----
# obiekt pomocniczy
share_xts_w = xts(data_w[, c("rate")], order.by = data_w$date)
# wykres WIG
dygraph(share_xts_w, main = "WIG tygodniowy") %>%
  dyRangeSelector(height = 50)
# INTERPRETACJA:
#   grupowanie wariancji

# wykres ACF
acf(data_w$rate, lag.max = 36, na.action = na.pass,
    ylim = c(-0.3, 0.3), 
    col = "darkblue", lwd = 7,
    main = "Wykres ACF stóp zwrotu (tygodniowe)")
pacf(data_w$rate, lag.max = 36, na.action = na.pass,
     ylim = c(-0.3, 0.3), 
     col = "darkblue", lwd = 7,
     main = "Wykres PACF stóp zwrotu (tygodniowe)")
# INTERPRETACJA:
#   raczej nieistotne opóźnienia, może będzie MA(1)

# test Durbina Watsona na autokorelację reszt (wykrywanie procesu MA)
durbinWatsonTest(lm(data_w$rate ~ 1),
                 max.lag = 5) # 5 pierwszych opóźnień
# INTERPRETACJA:
#   H0: corr = 0 => istotne drugie opóźnienie -> potencjał na uwzględnienie procesu MA(2)

# Wykres ACF^2 dla stóp zwrotu (efekty ARCH)
acf(data_w$rate^2, lag.max = 100, na.action = na.pass,
    ylim = c(-0.3, 0.3),
    col = "darkblue", lwd = 7,
    main = "Wykres ACF kwadratów stóp zwrotu (tygodniowe)")
# INTERPRETACJA:
#   no tak nie bardzo są te efekty ARCH 

# test ARCH
ArchTest(data_w$rate, lags = 3)
# ITNERPRETACJA:
#   test arch dla 3 opóźnień silnie odrzuca H0 o braku ARCH - formalne potwierdzenie ACF^2,
#   Warto sprawdzić czy opóźnienie 1 jednak nie jest istotne

durbinWatsonTest(lm(data_w$rate^2 ~ 1),
                 max.lag = 5) # 5 pierwszych opóźnień
# INTERPRETACJA:
#     Tylko opóźnienie 2 jest istotne, formalne potwierdzenie słabego efektu ARCH

# Histogram + gęstość dopasowanego rozkładu normalnego
hist(data_w$rate, prob = T, breaks = 100, main = "Histogram log-stóp zwrotu WIG (tygodniowe)", xlab="returns", col="skyblue1")
curve(dnorm(x, mean = mean(data_w$rate, na.rm = T),
            sd  = sd(data_w$rate, na.rm = T)),
      col = "darkblue", lwd = 2, add = TRUE,
)
# INTERPRETACJA:
#   minimalnie grubsze ogony większa kurtoza w porównaniu z rozkładem normalnym o empirycznych parametrach, lewy ogon wydaje się dłuższy
#   co oznacza lewostronną asymetrię rozkładu

# statystyki opisowe
empstats_w <- basicStats(data_w$rate)
# INTERPRETACJA:
#   Wysoka kurtoza (ok 19), ujemna skośność (ok -2.25) - asymetria lewostronna (dłuższy lewy ogon rozkładu)
#   sugeruje to możliwość lepszych wyników dla modeli uwzględniających asymetrię

# formalny test Jarque-Bera na normalność rozkładu
jarque.bera.test(data_w$rate)
# INTERPRETACJA:
#   silnie odrzucona H0 o normalności rozkładu -> dotychczasowe wnioski są zatem poprawne

# ---- MODELOWANIE I PREDYKCJA JEDNODNIOWA ----
# przygotowanie danych do modelowania
data_w$obs = 1:length(data_w$rate)
start <- data_w$obs[data_w$date == as.Date("2020-01-05")] # pierwszy dzien out of the sample
stop <- data_w$obs[data_w$date == as.Date("2020-12-27")] # ostatni dzien out of the sample
dta = data_w[data_w$obs <= (start-1), ]

# specyfikacja
# *** MA(2)-GARCH(1,1) e~Norm***
spec.garch <- ugarchspec(variance.model = list(model = "sGARCH",
                                               garchOrder = c(1, 1)),
                         mean.model = list(armaOrder = c(0, 2),
                                           include.mean = F),
                         distribution.model = "norm")
# estymacja
garch11ma2 <- ugarchfit(spec = spec.garch, data = dta$rate)
garch11ma2

# INTERPRETACJA WYDRUKU:
# test Ljunga-Boxa dla standaryzowanych reszt (H0: brak korelacji - przyjęta zawsze) =>
#     autokorelacja reszt została wyjaśniona

# test Ljunga-Boxa dla kwadratów standaryzowanych reszt (H0: brak korelacji - przyjęta zawsze) =>
#     autokorelacja kwadratów reszt została wyjaśniona

# test ważony ARCH LM (H0: brak efektów ARCH - przyjęta zawsze) =>
#     efekt ARCH został wyjaśniony
# MA nieistotne
# tutaj są wydruki do pojedynczych estymacji, tak żeby popatrzeć w parametry
plot(garch11ma2@fit$residuals, type="l")

# *** GARCH(1,1) e~Norm ***
spec.garch <- ugarchspec(variance.model = list(model = "sGARCH",
                                               garchOrder = c(1, 1)),
                         mean.model = list(armaOrder = c(0, 0),
                                           include.mean = F),
                         distribution.model = "norm")
# estymacja
garch11 <- ugarchfit(spec = spec.garch, data = dta$rate)
garch11
# Testy Ljunga-Boxa na biały szum
#    zmienna
Box.test(data_w$rate, type = "Ljung-Box", lag = 25)
#    reszty z modelu MA(2)-GARCH(1,1)
Box.test(garch11ma2@fit$residuals, type = "Ljung-Box", lag = 25)
#    reszty z modelu GARCH(1,1)
Box.test(garch11@fit$residuals, type = "Ljung-Box", lag = 25)
# INTERPRETACJA:
#    analizowana zmienna nie jest białym szumem, ale reszty z GARCH(1,1) już tak
#    nie ma co kombinować przy równaniu średniej

# ---- PĘTLA PREDYKCYJNA ----
# INSTRUKCJA:
#   dodając nowy model należy:
#     1. utworzyć wektor do przechowywania predykcji
#     2. utworzyć specyfikację
#     3. oszacować model w pętli
#     4. zrobić predykcje w pętli
#     5. zapisać wyniki predykcji do wektora w pętli
#     6. zapisać wektor predykcji (podniesiony do kwadratu) jako kolumnę data3

data3 <- data_w[start:stop, ] # df do zapisywania predykcji wariancji


# wektory do przechowywania rezultatów (każdy model musi mieć swój wektor)
sigma_preds.garch11.norm <- rep(NA, times = stop-start+1) # GARCH(1,1) normalny
sigma_preds.garch11.std <- rep(NA, times = stop-start+1) # GARCH(1,1) t-student
sigma_preds.garch11.ged <- rep(NA, times = stop-start+1) # GARCH(1,1) ged
sigma_preds.garch11.sstd <- rep(NA, times = stop-start+1) # GARCH(1,1) t-student skośny

sigma_preds.egarch11.norm <- rep(NA, times = stop-start+1) # EGARCH(1,1) normalny
sigma_preds.egarch11.std <- rep(NA, times = stop-start+1) # EGARCH(1,1) t-student
sigma_preds.egarch11.ged <- rep(NA, times = stop-start+1) # EGARCH(1,1) ged
sigma_preds.egarch11.sstd <- rep(NA, times = stop-start+1) # EGARCH(1,1) t-student skośny

sigma_preds.tgarch11.norm <- rep(NA, times = stop-start+1) # TGARCH(1,1) normalny
sigma_preds.tgarch11.std <- rep(NA, times = stop-start+1) # TGARCH(1,1) t-student
sigma_preds.tgarch11.ged <- rep(NA, times = stop-start+1) # TGARCH(1,1) ged
sigma_preds.tgarch11.sstd <- rep(NA, times = stop-start+1) # TGARCH(1,1) t-student skośny

sigma_preds.cgarch11.norm <- rep(NA, times = stop-start+1) # CGARCH(1,1) normalny
sigma_preds.cgarch11.std <- rep(NA, times = stop-start+1) # CGARCH(1,1) t-student
sigma_preds.cgarch11.ged <- rep(NA, times = stop-start+1) # CGARCH(1,1) ged
sigma_preds.cgarch11.sstd <- rep(NA, times = stop-start+1) # CGARCH(1,1) t-student skośny

# specyfikacje modeli
spec.garch11.norm <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                                   mean.model = list(armaOrder = c(0, 0), include.mean = F),
                                   distribution.model = "norm")
spec.garch11.std <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                                  mean.model = list(armaOrder = c(0, 0), include.mean = F),
                                  distribution.model = "std")
spec.garch11.ged <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                                  mean.model = list(armaOrder = c(0, 0), include.mean = F),
                                  distribution.model = "ged")
spec.garch11.sstd <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                                   mean.model = list(armaOrder = c(0, 0), include.mean = F),
                                   distribution.model = "sstd")

spec.egarch11.norm <- ugarchspec(variance.model = list(model = "eGARCH", garchOrder = c(1, 1)),
                                    mean.model = list(armaOrder = c(0, 0), include.mean = F),
                                    distribution.model = "norm")
spec.egarch11.std <- ugarchspec(variance.model = list(model = "eGARCH", garchOrder = c(1, 1)),
                                   mean.model = list(armaOrder = c(0, 0), include.mean = F),
                                   distribution.model = "std")
spec.egarch11.ged <- ugarchspec(variance.model = list(model = "eGARCH", garchOrder = c(1, 1)),
                                   mean.model = list(armaOrder = c(0, 0), include.mean = F),
                                   distribution.model = "ged")
spec.egarch11.sstd <- ugarchspec(variance.model = list(model = "eGARCH", garchOrder = c(1, 1)),
                                    mean.model = list(armaOrder = c(0, 0), include.mean = F),
                                    distribution.model = "sstd")

spec.tgarch11.norm <- ugarchspec(variance.model = list(model = "fGARCH", submodel = "TGARCH", garchOrder = c(1, 1)),
                                    mean.model = list(armaOrder = c(0, 0), include.mean = F),
                                    distribution.model = "norm")
spec.tgarch11.std <- ugarchspec(variance.model = list(model = "fGARCH", submodel = "TGARCH", garchOrder = c(1, 1)),
                                   mean.model = list(armaOrder = c(0, 0), include.mean = F),
                                   distribution.model = "std")
spec.tgarch11.ged <- ugarchspec(variance.model = list(model = "fGARCH", submodel = "TGARCH", garchOrder = c(1, 1)),
                                   mean.model = list(armaOrder = c(0, 0), include.mean = F),
                                   distribution.model = "ged")
spec.tgarch11.sstd <- ugarchspec(variance.model = list(model = "fGARCH", submodel = "TGARCH", garchOrder = c(1, 1)),
                                    mean.model = list(armaOrder = c(0, 0), include.mean = F),
                                    distribution.model = "sstd")

spec.cgarch11.norm <- ugarchspec(variance.model = list(model = "csGARCH", garchOrder = c(1, 1)),
                                    mean.model = list(armaOrder = c(0, 0), include.mean = F),
                                    distribution.model = "norm")
spec.cgarch11.std <- ugarchspec(variance.model = list(model = "csGARCH", garchOrder = c(1, 1)),
                                   mean.model = list(armaOrder = c(0, 0), include.mean = F),
                                   distribution.model = "std")
spec.cgarch11.ged <- ugarchspec(variance.model = list(model = "csGARCH", garchOrder = c(1, 1)),
                                   mean.model = list(armaOrder = c(0, 0), include.mean = F),
                                   distribution.model = "ged")
spec.cgarch11.sstd <- ugarchspec(variance.model = list(model = "csGARCH", garchOrder = c(1, 1)),
                                    mean.model = list(armaOrder = c(0, 0), include.mean = F),
                                    distribution.model = "sstd")


time1 = Sys.time()
for(k in start:stop){
  # dane dla bieżącej iteracji (chcę przewidzieć bieżące k)
  tmp.data <- data_w[data_w$obs <= (k-1), ] 
  
  # estymacje modeli
  tmp.garch11.norm <- ugarchfit(spec = spec.garch11.norm, data = na.omit(tmp.data$rate), solver = "hybrid")
  tmp.garch11.std <- ugarchfit(spec = spec.garch11.std, data = na.omit(tmp.data$rate), solver = "hybrid")
  tmp.garch11.ged <- ugarchfit(spec = spec.garch11.ged, data = na.omit(tmp.data$rate), solver = "hybrid")
  tmp.garch11.sstd <- ugarchfit(spec = spec.garch11.sstd, data = na.omit(tmp.data$rate), solver = "hybrid")
  
  tmp.egarch11.norm <- ugarchfit(spec = spec.egarch11.norm, data = na.omit(tmp.data$rate), solver = "hybrid")
  tmp.egarch11.std <- ugarchfit(spec = spec.egarch11.std, data = na.omit(tmp.data$rate), solver = "hybrid")
  tmp.egarch11.ged <- ugarchfit(spec = spec.egarch11.ged, data = na.omit(tmp.data$rate), solver = "hybrid")
  tmp.egarch11.sstd <- ugarchfit(spec = spec.egarch11.sstd, data = na.omit(tmp.data$rate), solver = "hybrid")
  
  tmp.tgarch11.norm <- ugarchfit(spec = spec.tgarch11.norm, data = na.omit(tmp.data$rate), solver = "hybrid")
  tmp.tgarch11.std <- ugarchfit(spec = spec.tgarch11.std, data = na.omit(tmp.data$rate), solver = "hybrid")
  tmp.tgarch11.ged <- ugarchfit(spec = spec.tgarch11.ged, data = na.omit(tmp.data$rate), solver = "hybrid")
  tmp.tgarch11.sstd <- ugarchfit(spec = spec.tgarch11.sstd, data = na.omit(tmp.data$rate), solver = "hybrid")
  
  tmp.cgarch11.norm <- ugarchfit(spec = spec.cgarch11.norm, data = na.omit(tmp.data$rate), solver = "hybrid")
  tmp.cgarch11.std <- ugarchfit(spec = spec.cgarch11.std, data = na.omit(tmp.data$rate), solver = "hybrid")
  tmp.cgarch11.ged <- ugarchfit(spec = spec.cgarch11.ged, data = na.omit(tmp.data$rate), solver = "hybrid")
  tmp.cgarch11.sstd <- ugarchfit(spec = spec.cgarch11.sstd, data = na.omit(tmp.data$rate), solver = "hybrid")
  
  # predykcje modeli (sigma - odchylenie standardowe)
  sigma_forecast.garch11.norm <- ugarchforecast(tmp.garch11.norm, n.ahead = 1)
  sigma_forecast.garch11.std <- ugarchforecast(tmp.garch11.std, n.ahead = 1)
  sigma_forecast.garch11.ged <- ugarchforecast(tmp.garch11.ged, n.ahead = 1)
  sigma_forecast.garch11.sstd <- ugarchforecast(tmp.garch11.sstd, n.ahead = 1)
  
  sigma_forecast.egarch11.norm <- ugarchforecast(tmp.egarch11.norm, n.ahead = 1)
  sigma_forecast.egarch11.std <- ugarchforecast(tmp.egarch11.std, n.ahead = 1)
  sigma_forecast.egarch11.ged <- ugarchforecast(tmp.egarch11.ged, n.ahead = 1)
  sigma_forecast.egarch11.sstd <- ugarchforecast(tmp.egarch11.sstd, n.ahead = 1)
  
  sigma_forecast.tgarch11.norm <- ugarchforecast(tmp.tgarch11.norm, n.ahead = 1)
  sigma_forecast.tgarch11.std <- ugarchforecast(tmp.tgarch11.std, n.ahead = 1)
  sigma_forecast.tgarch11.ged <- ugarchforecast(tmp.tgarch11.ged, n.ahead = 1)
  sigma_forecast.tgarch11.sstd <- ugarchforecast(tmp.tgarch11.sstd, n.ahead = 1)
  
  sigma_forecast.cgarch11.norm <- ugarchforecast(tmp.cgarch11.norm, n.ahead = 1)
  sigma_forecast.cgarch11.std <- ugarchforecast(tmp.cgarch11.std, n.ahead = 1)
  sigma_forecast.cgarch11.ged <- ugarchforecast(tmp.cgarch11.ged, n.ahead = 1)
  sigma_forecast.cgarch11.sstd <- ugarchforecast(tmp.cgarch11.sstd, n.ahead = 1)
  
  # zapis wyników do wektora
  sigma_preds.garch11.norm[k-start+1] <- sigma_forecast.garch11.norm@forecast$sigmaFor[1,1]
  sigma_preds.garch11.std[k-start+1] <- sigma_forecast.garch11.std@forecast$sigmaFor[1,1]
  sigma_preds.garch11.ged[k-start+1] <- sigma_forecast.garch11.ged@forecast$sigmaFor[1,1]
  sigma_preds.garch11.sstd[k-start+1] <- sigma_forecast.garch11.sstd@forecast$sigmaFor[1,1]
  
  sigma_preds.egarch11.norm[k-start+1] <- sigma_forecast.egarch11.norm@forecast$sigmaFor[1,1]
  sigma_preds.egarch11.std[k-start+1] <- sigma_forecast.egarch11.std@forecast$sigmaFor[1,1]
  sigma_preds.egarch11.ged[k-start+1] <- sigma_forecast.egarch11.ged@forecast$sigmaFor[1,1]
  sigma_preds.egarch11.sstd[k-start+1] <- sigma_forecast.egarch11.sstd@forecast$sigmaFor[1,1]
  
  sigma_preds.tgarch11.norm[k-start+1] <- sigma_forecast.tgarch11.norm@forecast$sigmaFor[1,1]
  sigma_preds.tgarch11.std[k-start+1] <- sigma_forecast.tgarch11.std@forecast$sigmaFor[1,1]
  sigma_preds.tgarch11.ged[k-start+1] <- sigma_forecast.tgarch11.ged@forecast$sigmaFor[1,1]
  sigma_preds.tgarch11.sstd[k-start+1] <- sigma_forecast.tgarch11.sstd@forecast$sigmaFor[1,1]
  
  sigma_preds.cgarch11.norm[k-start+1] <- sigma_forecast.cgarch11.norm@forecast$sigmaFor[1,1]
  sigma_preds.cgarch11.std[k-start+1] <- sigma_forecast.cgarch11.std@forecast$sigmaFor[1,1]
  sigma_preds.cgarch11.ged[k-start+1] <- sigma_forecast.cgarch11.ged@forecast$sigmaFor[1,1]
  sigma_preds.cgarch11.sstd[k-start+1] <- sigma_forecast.cgarch11.sstd@forecast$sigmaFor[1,1]
}


stop
time2 = Sys.time()
# czas obliczen
print(time2 - time1)
# dodanie kolumn ze zmiennością zrealizowaną (kwadrat stóp zwrotu) i przewidywaną (sigma^2)
data3$r_sq <- data3$rate^2 # zrealizowana
data3$sigma_preds_garch11.norm <- sigma_preds.garch11.norm^2
data3$sigma_preds_garch11.std <- sigma_preds.garch11.std^2
data3$sigma_preds_garch11.ged <- sigma_preds.garch11.ged^2
data3$sigma_preds_garch11.sstd <- sigma_preds.garch11.sstd^2

data3$sigma_preds_egarch11.norm <- sigma_preds.egarch11.norm^2
data3$sigma_preds_egarch11.std <- sigma_preds.egarch11.std^2
data3$sigma_preds_egarch11.ged <- sigma_preds.egarch11.ged^2
data3$sigma_preds_egarch11.sstd <- sigma_preds.egarch11.sstd^2

data3$sigma_preds_tgarch11.norm <- sigma_preds.tgarch11.norm^2
data3$sigma_preds_tgarch11.std <- sigma_preds.tgarch11.std^2
data3$sigma_preds_tgarch11.ged <- sigma_preds.tgarch11.ged^2
data3$sigma_preds_tgarch11.sstd <- sigma_preds.tgarch11.sstd^2

data3$sigma_preds_cgarch11.norm <- sigma_preds.cgarch11.norm^2
data3$sigma_preds_cgarch11.std <- sigma_preds.cgarch11.std^2
data3$sigma_preds_cgarch11.ged <- sigma_preds.cgarch11.ged^2
data3$sigma_preds_cgarch11.sstd <- sigma_preds.cgarch11.sstd^2

pred_tyh <- read.csv("predykcje_tygodniowe.csv")
mtrw <- as.data.frame(matrix(nrow = 16, ncol = 7))
for(i in 1:dim(mtrw)[1]){
  mtrw[i,1] <- ME(na.omit(data3[,i+5]), na.omit(data3$r_sq))
  mtrw[i,2] <- MAE(na.omit(data3[,i+5]),na.omit(data3$r_sq))
  mtrw[i,3] <- AMAPE(na.omit(data3[,i+5]), na.omit(data3$r_sq))
  mtrw[i,4] <- RMSE(na.omit(data3[,i+5]), na.omit(data3$r_sq))
  mtrw[i,5] <- MMEO(na.omit(data3[,i+5]), na.omit(data3$r_sq))
  mtrw[i,6] <- MMEU(na.omit(data3[,i+5]), na.omit(data3$r_sq))
  mtrw[i,7] <- TIC(na.omit(data3[,i+5]), na.omit(data3$r_sq))
}
colnames(mtrw) <- c("ME", "MAE", "AMAPE", "RMSE", "MMEO", "MMEU", "TIC")
rownames(mtrw) <- c("GARCH-N", "GARCH-T", "GARCH-G", "GARCH-ST",
                   "EGARCH-N", "EGARCH-T", "EGARCH-G", "EGARCH-ST",
                   "TGARCH-N", "TGARCH-T", "TGARCH-G", "TGARCH-ST",
                   "CGARCH-N", "CGARCH-T", "CGARCH-G", "CGARCH-ST")

mtrw

# *** Random walk model *** 
data3$sigma_preds_RW <- data_w$rate[(start-1):(stop-1)]^2
# Zapis wyników RW
mtrw["RW",1] <- ME(na.omit(data3$sigma_preds_RW), na.omit(data3$r_sq))
mtrw["RW",2] <- MAE(na.omit(data3$sigma_preds_RW),na.omit(data3$r_sq))
mtrw["RW",3] <- AMAPE(na.omit(data3$sigma_preds_RW), na.omit(data3$r_sq))
mtrw["RW",4] <- RMSE(na.omit(data3$sigma_preds_RW), na.omit(data3$r_sq))
mtrw["RW",5] <- MMEO(na.omit(data3$sigma_preds_RW), na.omit(data3$r_sq))
mtrw["RW",6] <- MMEU(na.omit(data3$sigma_preds_RW), na.omit(data3$r_sq))
mtrw["RW",7] <- TIC(na.omit(data3$sigma_preds_RW), na.omit(data3$r_sq))

# *** Historical Average model ***
sigma_preds.HA <- rep(NA, times = stop-start+1)
for(i in start:stop){
  sigma_preds.HA[i-start+1] <- sum(data_w$rate[1:(i-1)]^2)/length(data_w$rate[1:(i-2)])
}
data3$sigma_preds_HA <- sigma_preds.HA
# Zapis wyników RW
mtrw["HA",1] <- ME(na.omit(data3$sigma_preds_HA), na.omit(data3$r_sq))
mtrw["HA",2] <- MAE(na.omit(data3$sigma_preds_HA),na.omit(data3$r_sq))
mtrw["HA",3] <- AMAPE(na.omit(data3$sigma_preds_HA), na.omit(data3$r_sq))
mtrw["HA",4] <- RMSE(na.omit(data3$sigma_preds_HA), na.omit(data3$r_sq))
mtrw["HA",5] <- MMEO(na.omit(data3$sigma_preds_HA), na.omit(data3$r_sq))
mtrw["HA",6] <- MMEU(na.omit(data3$sigma_preds_HA), na.omit(data3$r_sq))
mtrw["HA",7] <- TIC(na.omit(data3$sigma_preds_HA), na.omit(data3$r_sq))
mtrw

# *** Moving 1-year-Average model ***
# jako długość roku przyjęto liczbę obserewacji z okresu predykcyjnego (mogą wystąpić różnice, ale pomijalne)
sigma_preds.MA1 <- rep(NA, times = stop-start+1)
for(i in start:stop){
  sigma_preds.MA1[i-start+1] <- mean(data_w$rate[(i-dim(data3)[1]-1):(i-1)]^2)
}
data3$sigma_preds_MA1 <- sigma_preds.MA1
# Zapis wyników MA1
mtrw["MA-1year",1] <- ME(na.omit(data3$sigma_preds_MA1), na.omit(data3$r_sq))
mtrw["MA-1year",2] <- MAE(na.omit(data3$sigma_preds_MA1),na.omit(data3$r_sq))
mtrw["MA-1year",3] <- AMAPE(na.omit(data3$sigma_preds_MA1), na.omit(data3$r_sq))
mtrw["MA-1year",4] <- RMSE(na.omit(data3$sigma_preds_MA1), na.omit(data3$r_sq))
mtrw["MA-1year",5] <- MMEO(na.omit(data3$sigma_preds_MA1), na.omit(data3$r_sq))
mtrw["MA-1year",6] <- MMEU(na.omit(data3$sigma_preds_MA1), na.omit(data3$r_sq))
mtrw["MA-1year",7] <- TIC(na.omit(data3$sigma_preds_MA1), na.omit(data3$r_sq))
mtrw

# *** Moving 3-year-Average model ***
# jako długość roku przyjęto liczbę obserewacji z okresu predykcyjnego (mogą wystąpić różnice, ale pomijalne)
sigma_preds.MA3 <- rep(NA, times = stop-start+1)
for(i in start:stop){
  sigma_preds.MA3[i-start+1] <- mean(data_w$rate[(i-3*dim(data3)[1]-1):(i-1)]^2)
}
data3$sigma_preds_MA3 <- sigma_preds.MA3
# Zapis wyników MA1
mtrw["MA-3year",1] <- ME(na.omit(data3$sigma_preds_MA3), na.omit(data3$r_sq))
mtrw["MA-3year",2] <- MAE(na.omit(data3$sigma_preds_MA3),na.omit(data3$r_sq))
mtrw["MA-3year",3] <- AMAPE(na.omit(data3$sigma_preds_MA3), na.omit(data3$r_sq))
mtrw["MA-3year",4] <- RMSE(na.omit(data3$sigma_preds_MA3), na.omit(data3$r_sq))
mtrw["MA-3year",5] <- MMEO(na.omit(data3$sigma_preds_MA3), na.omit(data3$r_sq))
mtrw["MA-3year",6] <- MMEU(na.omit(data3$sigma_preds_MA3), na.omit(data3$r_sq))
mtrw["MA-3year",7] <- TIC(na.omit(data3$sigma_preds_MA3), na.omit(data3$r_sq))
mtrw

#write.csv(data2, "predykcje_dzienne.csv", row.names = F)
#write.csv(data3, "predykcje_tygodniowe.csv",row.names = F)

