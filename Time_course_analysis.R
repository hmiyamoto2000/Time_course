#Time_course_analysis

#Make th folder for Time_course
dir.create("./Result_Time_course", showWarnings = TRUE, recursive = FALSE, mode = "0777")

dat
write.csv(dat, "./Result_Time_course/File_name.csv")

dd_s_t
sink('./Result_Time_course/File_name.txt', append = TRUE)
print (dd_s_t)
sink()


#全体のデータが正規分布かどうかを評価
library(MVN)
DD <- read.csv("C_whole_Med.csv")
input=DD[,2:5]
mvn(input,univariateTest="SW")

mvn_r=mvn(input,univariateTest="SW")
sink('./Result_Time_course/SW_mvn_r.txt', append = TRUE)
print (mvn_r)
sink()


#個体別で評価

library("vars")
C_1=read.csv("C_whole_Med.csv") #中央値　CとT
C_1_nonID=C_1[,2:9]
d_C_1=ts(C_1_nonID,start=c(1),freq=1)
plot(d_C_1)
VARselect(d_C_1, lag.max=4, type="const")$selection

VARselect_r=VARselect(d_C_1, lag.max=4, type="const")$selection
sink('./Result_Time_course/VARselect_r.txt', append = TRUE)
print (VARselect_r)
sink()

#あるいは
lag <- VARselect(d_C_1[,c(5:6)], lag.max=2)$selection[1]
sink('./Result_Time_course/lag_r.txt', append = TRUE)
print (lag)
sink()

#2:3は、この場合は、XB_p__Bacteroidetes, XB_g__Methanobrevibacter 

#あるいはts型に慣れてない場合データフレームにする。canada <- as.data.frame(as.matrix(Canada))を参考として、
AA <- as.data.frame(as.matrix(d_C_1))

split.screen(c(3,1))  # 描画デバイス分割
screen(1)
plot(AA$XB_g__Desulfovibrio, type='l', main='Desulfovibrio')

screen(2)
plot(AA$XB_p__Bacteroidetes, type='l', main='Bacteroidetes')

screen(3)
plot(AA$XB_g__Methanobrevibacter, type='l', main='Methanobrevibacter')

split.screen(c(3,1)) 
screen(1)
plot(AA$XB_g__Pseudobutyrivibrio, type='l', main='Pseudobutyrivibrio')

screen(2)
plot(AA$XB_g__Shuttleworthia, type='l', main='Shuttleworthia')


var.AA <- VAR(AA[,c(5:6)], p=2, type='const') #あるいはvar.3p<-VAR(d_C_1,p=2,type="const")
serial_var.AA_r=serial.test(var.AA, lags.pt = 2, type = "PT.asymptotic")
#serial.test(var.AA, lags.pt = 2, type = "PT.asymptotic") #デフォルトである "PT.asymptotic" 

sink('./Result_Time_course/serial_var.AA_r.txt', append = TRUE)
print (serial_var.AA_r)
sink()

summary(var.AA)
summary_var.AA=summary(var.AA)

sink('./Result_Time_course/summary_var.AA_r.txt', append = TRUE)
print (summary_var.AA)
sink()

plot(var.AA)
par(mfrow = c(2, 2))
plot(var.AA$varresult$XB_p__Bacteroidetes)
par(mfrow = c(2, 2))
plot(var.AA$varresult$XB_g__Methanobrevibacter)

#OLStest
Stability_var=stability(var.AA)
par(col="skyblue",col.axis="mediumpurple",mfrow=c(2,2))
plot(Stability_var$stability$XB_p__Bacteroidetes)
par(col="skyblue",col.axis="mediumpurple",mfrow=c(2,2))
plot(Stability_var$stability$XB_g__Methanobrevibacter)

#色なし
par(col="black",col.axis="black",mfrow=c(2,2))
plot(Stability_var$stability$XB_p__Bacteroidetes)
par(col="black",col.axis="black",mfrow=c(2,2))
plot(Stability_var$stability$XB_g__Methanobrevibacter)

sctest(var.AA) #できない

serial.test(var.AA,lags.pt=6)

serial_test_var.AA=serial.test(var.AA,lags.pt=6)

sink('./Result_Time_course/serial_test_var.AA_r.txt', append = TRUE)
print (serial_test_var.AA)
sink()


#グレンジャー因果の評価
causality(var.AA) #causality(var.AA,cause="XB_p__Bacteroidetes")

Granger=causality(var.AA)

sink('./Result_Time_course/Granger_r.txt', append = TRUE)
print (Granger)
sink()


#予測
var.3p.prd <- predict(var.AA, n.ahead =10, ci = 0.95) #10倍先を予測

sink('./Result_Time_course/var.3p.prd_r.txt', append = TRUE)
print (var.3p.prd)
sink()

od<- options(digits=3)   #表 示する数値の桁数を3とする
var.3p.prd
plot(var.3p.prd)


#ARIMAモデル https://www.rbt.his.u-fukui.ac.jp/~naniwa/pub/R4.html
#https://best-biostatistics.com/toukei-er/entry/basics-of-time-series-analysis/
library("forecast")
library(TSstudio)
library(MLmetrics)

#擬似数字羅列の場合
#set.seed(100)
#d_C_1 <- armaSim(n = 200, model = list(ar = c(0.8, -0.6, 0.7), d = 1,ma = c(-0.5, 0.6)), n.start = 300)

#こちらで読み込みした場合でも同じ名前で
C_1=read.csv("C_whole_Med.csv")
C_1=read.csv("CT_whole_Med.csv")

C_1_nonID=C_1[,2:5]
d_C_1=ts(C_1_nonID,start=c(1),freq=1)
ts.plot(d_C_1)
ts_plot(d_C_1, title = "Week",Ytitle = "Relative abundance")
autoplot(d_C_1)
autoplot(d_C_1,xlab='week',ylab='Relative abundance',main='Control group')
#xlim=c(1,11)とするとその範囲
#scale_x_continuous(breaks=seq(1,11,1))

ggtsdisplay(d_C_1[,4],main='Bacteroidetes')
#ACF(相関係数)、偏相関係数(PACF)を分けたグラフを表示

#ADF検定の例　#https://qiita.com/YM_DSKR/items/2528548913378bfbf9bc
library(tseries)
adf.test(d_C_1[,3])
adf_test=adf.test(d_C_1[,4])
sink('./Result_Time_course/adf_test_T_Bac_r.txt', append = TRUE)
print (adf_test)
sink()
#Augmented Dickey-Fuller Test
#data:  d_C_1[, 2]
#Dickey-Fuller = -3, Lag order = 2, p-value = 0.3
#alternative hypothesis: stationary
#jarque.bera.test(d_C_1[,2:9])

#Jargu=jarque.bera.test(C_1_nonID[,2])
#sink('./Result_Time_course/jarque_bera_test_r.txt', append = TRUE)
#print (Jargu)
#sink()

#pvalue0.3となり帰無仮説は棄却されない→単位根あり→差分とるべき(非定常過程)
#差分をとって再検証
adf.test(diff(d_C_1[,2]))
#Augmented Dickey-Fuller Test
#data:  diff(d_C_1[, 2])
#Dickey-Fuller = -2, Lag order = 2, p-value = 0.5
#alternative hypothesis: stationary
#p0.5となり帰無仮説は棄却されない→単位根あり→差分をとった方が良い
#もしp<0.05なら、帰無仮説は棄却される→単位根なし→差分をとらなくて良い(定常過程)

#差分を取る場合のplot
plot(diff(d_C_1[,2]),main='Bacteroidetes')
autoplot(diff(d_C_1),xlab='week',ylab='Relative abundance',main='Control group')

#KPSS検定
library(urca)
summary(ur.kpss(d_C_1[,2]))
####################### 
# KPSS Unit Root Test # 
####################### 
#Test is of type: mu with 2 lags. 
#Value of test-statistic is: 0.469 
#Critical value for a significance level of: 
#               10pct  5pct 2.5pct  1pct
#critical values 0.347 0.463  0.574 0.739
#有意水準5%とした時の棄却点が0.463であり、棄却点0.469とあまり変わらない
#単位根がないという帰無仮説が棄却されないので、単位根が必ずしもあるとは言えない(差分を取る必要がない)

#もしValue of test-statistic isが、15.4007であれば、単位根がないという帰無仮説が棄却されたので、単位根があるとみなす
#ndiffs関数で何回差分を取るべきかを判定することができる
library(forecast)
ndiffs(d_C_1[,2])
summary(ur.kpss(diff(d_C_1[,2])))
####################### 
# KPSS Unit Root Test # 
####################### 
#Test is of type: mu with 2 lags. 
#Value of test-statistic is: 0.217 
#Critical value for a significance level of: 
#  10pct  5pct 2.5pct  1pct
#critical values 0.347 0.463  0.574 0.739
#一回の差分で0.463を下回る0.217となった。

#二階差分の場合
summary(ur.kpss(diff(diff(d_C_1[,2]))))




ts.plot(C_1_nonID)
acf(d_C_1) #自己相関係数を図示
ts.plot(diff(d_C_1))#前後のデータの差分を求め、その自己相関係数を図示
ts_lags(d_C_1,lags = 1:12)

##Phillips-Perron検定
#時系列データがランダムにふらふらしているだけなのか、ある一定の傾向があるのかを検定する方法がPhillips-Perron検定
PP.test(d_C_1[,2])
Phillips_Perron=PP.test(d_C_1[,2])
sink('./Result_Time_course/Phillips_Perron_r.txt', append = TRUE)
print (Phillips_Perron)
sink()

#p>0.05の時: 結果は以下の通り、ランダムであることは否定できない（帰無仮説はランダムウォークである）
#単位根がない　→ ランダムウォークではない
#adf.test(d_C_1[,2]) #計算できない


##Durbin-Watson統計量の評価
library(car)
t=C_1[,6] #七列目にtの番号順を入れる
lm1 <- lm(d_C_1[,2] ~ t)
plot(t, d_C_1[,2], type = "l")
residualPlot(lm1, variable = "t")
durbinWatsonTest(lm1)

Durbin_Watson=durbinWatsonTest(lm1)
sink('./Result_Time_course/Durbin_Watson_r.txt', append = TRUE)
print (Durbin_Watson)
sink()


#Durbin-Watson統計量は、正の自己相関があるときに0に近い値、
#負の自己相関があるときに4に近い値、自己相関がないときに2となります。
#この場合は0.3727528と、0に近い値となりました。
#自己相関がないことを帰無仮説とした検定ではp値が0となり、帰無仮説は棄却されました。

#Bacteroidetesの評価:
#lag Autocorrelation D-W Statistic p-value
#1      -0.4407488      2.877169   0.242
#Alternative hypothesis: rho != 0
#Durbin-Watson統計量は2.877169と、2に近い値となりました。
#また、自己相関がないことを帰無仮説とした検定では、帰無仮説は棄却されませんでした。

#https://ito4303.sakura.ne.jp/posts/2023-07-09-durbin-watson/
#https://htsuda.net/stats/regression.html



acf(diff(d_C_1), na.action = na.omit)
auto.arima(d_C_1[,2])#AIC （赤池情報量基準）で最適な ARIMA モデルを推定する。
ARIMA=auto.arima(d_C_1[,2])
sink('./Result_Time_course/ARIMA_Bac_r.txt', append = TRUE)
print (ARIMA)
sink()


fit <- armaFit(~ arima(0,0,0), data = d_C_1) #うまくいかない
#推定結果が ARIMA(3,1,3) となったので、パラメータの推定とフィッティングを
#行う。

opar = par(mfrow = c(2,2), cex = 0.7)#うまくいかない
summary(fit)#うまくいかない
par(opar)#うまくいかない
ts.plot(fitted(fit), y, col = c(3, 1))#うまくいかない
dev.off()#うまくいかない

#時系列データがランダムにふらふらしているだけなのか、ある一定の傾向があるのかを検定する方法がPhillips-Perron検定
PP.test(d_C_1[,2])
Phillips_Perron=PP.test(d_C_1[,2])
#p>0.05の時: 結果は以下の通り、ランダムであることは否定できない（帰無仮説はランダムウォークである）
#単位根がない　→ ランダムウォークではない
#adf.test(d_C_1[,2]) #計算できない

#上記の前提の上で実施

ar1 <- ar(d_C_1[,2]) #自己相関 AutoRegressive AR モデル
aa1 <- auto.arima(d_C_1[,2], trace=TRUE, stepwise=FALSE, max.p=10, max.q=10)
sink('./Result_Time_course/ARIMA_r.txt', append = TRUE)
print (aa1)
sink()

fc1 <- forecast(aa1,h=10)
plot(fc1)

lines(d_C_1[,2])
grid()

#Reference
#https://blog.statsbeginner.net/entry/2018/05/04/223605
#https://ito4303.sakura.ne.jp/posts/2023-07-09-durbin-watson/
#https://htsuda.net/stats/regression.html
#ADF検定の例　#https://qiita.com/YM_DSKR/items/2528548913378bfbf9bc
#ARIMAモデル https://www.rbt.his.u-fukui.ac.jp/~naniwa/pub/R4.html
#https://best-biostatistics.com/toukei-er/entry/basics-of-time-series-analysis/
