# load packages] ----------------------------------------------------------

library(tidyr)
library(dplyr)
library(nlme)
library(magrittr)
library(ggplot2)
library(fitdistrplus)


# load data ---------------------------------------------------------------

wilson<-read.table("C:/Users/Farjala/Documents/Manuscritos/FrontiersMicrobiology/Wilson.txt",
                   header=TRUE, stringsAsFactors = FALSE)
str(wilson)



# olhada geral na distribuicao dos dados ----------------------------------

qplot(Production, data=metabolism, geom="histogram")
qplot(Respiration, data=metabolism, geom="histogram")
qplot(BGE, data=metabolism, geom="histogram")

plot(density(Production))
plot(density(Respiration))
plot(density(BGE))

# vai ficando claro que os dados n?o s?o nada normais. A BGE, que ? uma fra??o, foge particularmente da distribui??o normal (esperado!) 
## olhando pelo ajuste da distribui??o gaussiana isto fica claro!
metabolism  <- metabolism %>%
  mutate(z.Production = (Production-mean(Production))/sd(Production))

qqnorm(metabolism$z.Production)
with(metabolism, qqnorm(z.Production))
abline(0,1)

z.Respiration<-(Respiration-mean(Respiration))/sd(Respiration)
qqnorm(z.Respiration)
abline(0,1)

z.BGE<-(BGE-mean(BGE))/sd(BGE)
qqnorm(z.BGE)
abline(0,1)


# data reorganization -----------------------------------------------------

meta <- wilson %>%
  select(Lake:BGE) %>%
  gather(key = "response", value = "value", Respiration, Production, BGE) %>%
  mutate(N = ifelse(Treatment %in% c("N","NP"), yes = "N", no = "Control"),
         P = ifelse(Treatment %in% c("P","NP"), yes = "P", no = "Control"),
         Treatment = ifelse(Treatment == "B", yes = "Control", no = Treatment),
         Treatment = factor(Treatment, levels = c("Control","N","P","NP")),
         Lake = factor(Lake, levels = c("Cab", "Carap", "Comp", "Aboi", "Atol")),
         clarity = ifelse(Lake %in% c("Cab","Carap"), yes = "Clear", no = "Dark"))
head(meta)

## checando se todos os lagos tem nomes
any(is.na(meta$Lake))

# analysis (linear mixed models) ------------------------------------------

#GLMM com temperature como vari?vel cont?nua e o Bloco inserido dentro de Lagoa como vari?veis aleat?rias

PB <- filter (meta, response %in% c("Production"))
lmeProduction<-lme(value ~ (N*P)*Temperature, random = ~1|Lake/Block, data=PB)
summary(lmeProduction)

# opção de figura com regressões
PB <- filter (meta, response %in% c("Production"))
head(PB)
coef(lm(value ~ Temperature, data = PB))
ggplot(PB, aes(x=Temperature, y = value, colour=Lake)) +
  geom_point(position = position_jitter(w = 0.5), size=2.5) +
  stat_summary(fun.y = median, geom="line")+
  geom_abline(intercept = 0.1275176256, slope = -0.0008845793, size = 1.0)+
  scale_color_manual(values = c("Cab"="chartreuse3","Carap"="blue", "Comp"="orange", "Aboi"="red", "Atol"="grey30"))+
  scale_y_continuous(limits = c(0.01, 0.26))+
  xlab("Temperature") + ylab("Bacterial Production") +
  theme_bw() + theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank()) +
  theme(axis.title = element_text(size = rel(1.7))) + theme(axis.text = element_text(size = rel(1.5))) #mudei o tamanho das escalas e legendas


# opção de figura com boxplots
p1 <- ggplot(PB, aes(x = factor(Temperature), y = value)) + geom_boxplot(aes(fill = Treatment)) + scale_y_continuous(limits = c(0, 0.30)) #retirei os outliers da figura e mudei a escala para otimizar o espa?o 
p2 <- p1 + xlab("Temperature") + ylab("Production") # alterei as legendas (n?o foi especialmente necass?rio aqui)
p3 <- p2 + theme_bw() + theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank()) # retirei as linhas de grade e peinel de fundo cinza
p4 <- p3 + theme(axis.title = element_text(size = rel(1.7))) + theme(axis.text = element_text(size = rel(1.5))) #mudei o tamanho das escalas e legendas
p4

RB <- filter (meta, response %in% c("Respiration"))
lmeRespiration<-lme(value~(N*P)*Temperature, random=~1|Lake/Block, data=RB)
summary(lmeRespiration)

RB <- filter (meta, response %in% c("Respiration"))
head(RB)
coef(lm(value ~ Temperature, data = RB))
ggplot(RB, aes(x=Temperature, y = value, colour=Lake))+
  geom_point(position = position_jitter(w = 0.5), size=2.5) +
  geom_abline(intercept = -0.05480916, slope = 0.01438665, size = 1.0)+
  stat_summary(fun.y = median, geom="line")+
  scale_color_manual(values = c("Cab"="chartreuse3","Carap"="blue", "Comp"="orange", "Aboi"="red", "Atol"="grey30"))+
  scale_y_continuous(limits = c(0.01, 0.85))+
  xlab("Temperature") + ylab("Bacterial Respiration") +
  theme_bw() + theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank()) +
  theme(axis.title = element_text(size = rel(1.7))) + theme(axis.text = element_text(size = rel(1.5))) #mudei o tamanho das escalas e legendas

p1 <- ggplot(RB, aes(x = factor(Temperature), y = value)) + geom_boxplot(aes(fill = Treatment)) + scale_y_continuous(limits = c(0, 0.60)) #retirei os outliers da figura e mudei a escala para otimizar o espa?o 
p2 <- p1 + xlab("Temperature") + ylab("Respiration") # alterei as legendas (n?o foi especialmente necass?rio aqui)
p3 <- p2 + theme_bw() + theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank()) # retirei as linhas de grade e peinel de fundo cinza
p4 <- p3 + theme(axis.title = element_text(size = rel(1.7))) + theme(axis.text = element_text(size = rel(1.5))) #mudei o tamanho das escalas e legendas
p4


BGE <- filter (meta, response %in% c("BGE"))
lmeBGE<-lme(value~(N*P)*Temperature, random=~1|Lake/Block, data=BGE)
summary(lmeBGE)

BGE <- filter (meta, response %in% c("BGE"))
head(BGE)
coef(lm(value ~ Temperature, data = BGE))
ggplot(BGE, aes(x=Temperature, y = value, colour=Lake))+
  geom_point(position = position_jitter(w = 0.5), size=2.5) +
  geom_abline(intercept = 0.488336043, slope = -0.01156826, size = 1.0)+
  stat_summary(fun.y = median, geom="line")+
  scale_color_manual(values = c("Cab"="chartreuse3","Carap"="blue", "Comp"="orange", "Aboi"="red", "Atol"="grey30"))+
  scale_y_continuous(limits = c(0, 0.60))+
  xlab("Temperature") + ylab("BGE") +
  theme_bw() + theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank()) +
  theme(axis.title = element_text(size = rel(1.7))) + theme(axis.text = element_text(size = rel(1.5))) #mudei o tamanho das escalas e legendas

p1 <- ggplot(BGE, aes(x = factor(Temperature), y = value)) + geom_boxplot(aes(fill = Treatment)) + scale_y_continuous(limits = c(0, 0.60)) #retirei os outliers da figura e mudei a escala para otimizar o espa?o 
p2 <- p1 + xlab("Temperature") + ylab("BGE") # alterei as legendas (n?o foi especialmente necass?rio aqui)
p3 <- p2 + theme_bw() + theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank()) # retirei as linhas de grade e peinel de fundo cinza
p4 <- p3 + theme(axis.title = element_text(size = rel(1.5))) + theme(axis.text = element_text(size = rel(1.0))) #mudei o tamanho das escalas e legendas
p4

# finalmente, juntados todos os plots em uma Figura apenas para publica??o

ggplot(meta, aes(x = factor(Temperature), y = value)) + geom_boxplot(aes(fill = Treatment)) + #scale_y_continuous(limits = c(0, 0.6)) + #retirei os outliers da figura e mudei a escala para otimizar o espa?o 
  xlab("Temperature") + ylab("response") + # alterei as legendas (n?o foi especialmente necass?rio aqui)
  theme_bw() + theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank()) + # retirei as linhas de grade e peinel de fundo cinza
  theme(axis.title = element_text(size = rel(1.5))) + theme(axis.text = element_text(size = rel(1.0))) + #mudei o tamanho das escalas e legendas
  facet_wrap( ~ response, scales = "free")


# olhando agora dentro de cada lagoa. H? padr?es entre as lagoas em rela??o ? adi??o de nutrientes? ---------
#Cabiunas
PBCab <- filter (meta, response %in% c("Production"), Lake %in% c("Cab"))
lmePBCab <-lme(value~(N*P), random=~1|Temperature/Block, data=PBCab)
summary(lmePBCab)
p1 <- ggplot(PBCab, aes(x = Treatment, y = value)) + geom_boxplot() + scale_y_continuous(limits = c(0, 0.30)) #retirei os outliers da figura e mudei a escala para otimizar o espa?o 
p2 <- p1 + xlab("Treatment") + ylab("Production") # alterei as legendas (n?o foi especialmente necass?rio aqui)
p3 <- p2 + theme_bw() + theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank()) # retirei as linhas de grade e peinel de fundo cinza
p4 <- p3 + theme(axis.title = element_text(size = rel(1.5))) + theme(axis.text = element_text(size = rel(1.0))) #mudei o tamanho das escalas e legendas
p4

RBCab <- filter (meta, response %in% c("Respiration"), Lake %in% c("Cab"))
lmeRBCab <-lme(value~(N*P), random=~1|Temperature/Block, data=RBCab)
summary(lmeRBCab)
p1 <- ggplot(RBCab, aes(x = Treatment, y = value)) + geom_boxplot() + scale_y_continuous(limits = c(0, 0.60)) #retirei os outliers da figura e mudei a escala para otimizar o espa?o 
p2 <- p1 + xlab("Treatment") + ylab("Respiration") # alterei as legendas (n?o foi especialmente necass?rio aqui)
p3 <- p2 + theme_bw() + theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank()) # retirei as linhas de grade e peinel de fundo cinza
p4 <- p3 + theme(axis.title = element_text(size = rel(1.5))) + theme(axis.text = element_text(size = rel(1.0))) #mudei o tamanho das escalas e legendas
p4

BGECab <- filter (meta, response %in% c("BGE"), Lake %in% c("Cab"))
lmeBGECab <-lme(value~(N*P), random=~1|Temperature/Block, data=BGECab)
summary(lmeBGECab)
p1 <- ggplot(BGECab, aes(x = Treatment, y = value)) + geom_boxplot() + scale_y_continuous(limits = c(0, 0.60)) #retirei os outliers da figura e mudei a escala para otimizar o espa?o 
p2 <- p1 + xlab("Treatment") + ylab("BGE") # alterei as legendas (n?o foi especialmente necass?rio aqui)
p3 <- p2 + theme_bw() + theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank()) # retirei as linhas de grade e peinel de fundo cinza
p4 <- p3 + theme(axis.title = element_text(size = rel(1.5))) + theme(axis.text = element_text(size = rel(1.0))) #mudei o tamanho das escalas e legendas
p4

#Carapebus
PBCarap <- filter (meta, response %in% c("Production"), Lake %in% c("Carap"))
lmePBCarap <-lme(value~(N*P), random=~1|Temperature/Block, data=PBCarap)
summary(lmePBCarap)
p1 <- ggplot(PBCarap, aes(x = Treatment, y = value)) + geom_boxplot() + scale_y_continuous(limits = c(0, 0.30)) #retirei os outliers da figura e mudei a escala para otimizar o espa?o 
p2 <- p1 + xlab("Treatment") + ylab("Production") # alterei as legendas (n?o foi especialmente necass?rio aqui)
p3 <- p2 + theme_bw() + theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank()) # retirei as linhas de grade e peinel de fundo cinza
p4 <- p3 + theme(axis.title = element_text(size = rel(1.5))) + theme(axis.text = element_text(size = rel(1.0))) #mudei o tamanho das escalas e legendas
p4

RBCarap <- filter (meta, response %in% c("Respiration"), Lake %in% c("Carap"))
lmeRBCarap <-lme(value~(N*P), random=~1|Temperature/Block, data=RBCarap)
summary(lmeRBCarap)
p1 <- ggplot(RBCarap, aes(x = Treatment, y = value)) + geom_boxplot() + scale_y_continuous(limits = c(0, 0.60)) #retirei os outliers da figura e mudei a escala para otimizar o espa?o 
p2 <- p1 + xlab("Treatment") + ylab("Respiration") # alterei as legendas (n?o foi especialmente necass?rio aqui)
p3 <- p2 + theme_bw() + theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank()) # retirei as linhas de grade e peinel de fundo cinza
p4 <- p3 + theme(axis.title = element_text(size = rel(1.5))) + theme(axis.text = element_text(size = rel(1.0))) #mudei o tamanho das escalas e legendas
p4

BGECarap <- filter (meta, response %in% c("BGE"), Lake %in% c("Carap"))
lmeBGECarap <-lme(value~(N*P), random=~1|Temperature/Block, data=BGECarap)
summary(lmeBGECarap)
p1 <- ggplot(BGECarap, aes(x = Treatment, y = value)) + geom_boxplot() + scale_y_continuous(limits = c(0, 0.60)) #retirei os outliers da figura e mudei a escala para otimizar o espa?o 
p2 <- p1 + xlab("Treatment") + ylab("BGE") # alterei as legendas (n?o foi especialmente necass?rio aqui)
p3 <- p2 + theme_bw() + theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank()) # retirei as linhas de grade e peinel de fundo cinza
p4 <- p3 + theme(axis.title = element_text(size = rel(1.5))) + theme(axis.text = element_text(size = rel(1.0))) #mudei o tamanho das escalas e legendas
p4

#Comprida
PBComp <- filter (meta, response %in% c("Production"), Lake %in% c("Comp"))
lmePBComp <-lme(value~(N*P), random=~1|Temperature/Block, data=PBComp)
summary(lmePBComp)
p1 <- ggplot(PBComp, aes(x = Treatment, y = value)) + geom_boxplot() + scale_y_continuous(limits = c(0, 0.20)) #retirei os outliers da figura e mudei a escala para otimizar o espa?o 
p2 <- p1 + xlab("Treatment") + ylab("Production") # alterei as legendas (n?o foi especialmente necass?rio aqui)
p3 <- p2 + theme_bw() + theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank()) # retirei as linhas de grade e peinel de fundo cinza
p4 <- p3 + theme(axis.title = element_text(size = rel(1.5))) + theme(axis.text = element_text(size = rel(1.0))) #mudei o tamanho das escalas e legendas
p4

RBComp <- filter (meta, response %in% c("Respiration"), Lake %in% c("Comp"))
lmeRBComp <-lme(value~(N*P), random=~1|Temperature/Block, data=RBComp)
summary(lmeRBComp)
p1 <- ggplot(RBComp, aes(x = Treatment, y = value)) + geom_boxplot() + scale_y_continuous(limits = c(0, 0.60)) #retirei os outliers da figura e mudei a escala para otimizar o espa?o 
p2 <- p1 + xlab("Treatment") + ylab("Respiration") # alterei as legendas (n?o foi especialmente necass?rio aqui)
p3 <- p2 + theme_bw() + theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank()) # retirei as linhas de grade e peinel de fundo cinza
p4 <- p3 + theme(axis.title = element_text(size = rel(1.5))) + theme(axis.text = element_text(size = rel(1.0))) #mudei o tamanho das escalas e legendas
p4

BGEComp <- filter (meta, response %in% c("BGE"), Lake %in% c("Comp"))
lmeBGEComp <-lme(value~(N*P), random=~1|Temperature/Block, data=BGEComp)
summary(lmeBGEComp)
p1 <- ggplot(BGEComp, aes(x = Treatment, y = value)) + geom_boxplot() + scale_y_continuous(limits = c(0, 0.60)) #retirei os outliers da figura e mudei a escala para otimizar o espa?o 
p2 <- p1 + xlab("Treatment") + ylab("BGE") # alterei as legendas (n?o foi especialmente necass?rio aqui)
p3 <- p2 + theme_bw() + theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank()) # retirei as linhas de grade e peinel de fundo cinza
p4 <- p3 + theme(axis.title = element_text(size = rel(1.5))) + theme(axis.text = element_text(size = rel(1.0))) #mudei o tamanho das escalas e legendas
p4

#Amarra-Boi
PBAboi <- filter (meta, response %in% c("Production"), Lake %in% c("Aboi"))
lmePBAboi <-lme(value~(N*P), random=~1|Temperature/Block, data=PBAboi)
summary(lmePBAboi)
p1 <- ggplot(PBAboi, aes(x = Treatment, y = value)) + geom_boxplot() + scale_y_continuous(limits = c(0, 0.20)) #retirei os outliers da figura e mudei a escala para otimizar o espa?o 
p2 <- p1 + xlab("Treatment") + ylab("Production") # alterei as legendas (n?o foi especialmente necass?rio aqui)
p3 <- p2 + theme_bw() + theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank()) # retirei as linhas de grade e peinel de fundo cinza
p4 <- p3 + theme(axis.title = element_text(size = rel(1.5))) + theme(axis.text = element_text(size = rel(1.0))) #mudei o tamanho das escalas e legendas
p4

RBAboi <- filter (meta, response %in% c("Respiration"), Lake %in% c("Aboi"))
lmeRBAboi <-lme(value~(N*P), random=~1|Temperature/Block, data=RBAboi)
summary(lmeRBAboi)
p1 <- ggplot(RBAboi, aes(x = Treatment, y = value)) + geom_boxplot() + scale_y_continuous(limits = c(0, 0.60)) #retirei os outliers da figura e mudei a escala para otimizar o espa?o 
p2 <- p1 + xlab("Treatment") + ylab("Respiration") # alterei as legendas (n?o foi especialmente necass?rio aqui)
p3 <- p2 + theme_bw() + theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank()) # retirei as linhas de grade e peinel de fundo cinza
p4 <- p3 + theme(axis.title = element_text(size = rel(1.5))) + theme(axis.text = element_text(size = rel(1.0))) #mudei o tamanho das escalas e legendas
p4

BGEAboi <- filter (meta, response %in% c("BGE"), Lake %in% c("Aboi"))
lmeBGEAboi <-lme(value~(N*P), random=~1|Temperature/Block, data=BGEAboi)
summary(lmeBGEAboi)
p1 <- ggplot(BGEAboi, aes(x = Treatment, y = value)) + geom_boxplot() + scale_y_continuous(limits = c(0, 0.60)) #retirei os outliers da figura e mudei a escala para otimizar o espa?o 
p2 <- p1 + xlab("Treatment") + ylab("BGE") # alterei as legendas (n?o foi especialmente necass?rio aqui)
p3 <- p2 + theme_bw() + theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank()) # retirei as linhas de grade e peinel de fundo cinza
p4 <- p3 + theme(axis.title = element_text(size = rel(1.5))) + theme(axis.text = element_text(size = rel(1.0))) #mudei o tamanho das escalas e legendas
p4

#Atoleiro
PBAtol <- filter (meta, response %in% c("Production"), Lake %in% c("Atol"))
lmePBAtol <-lme(value~(N*P), random=~1|Temperature/Block, data=PBAtol)
summary(lmePBAtol)
p1 <- ggplot(PBAtol, aes(x = Treatment, y = value)) + geom_boxplot() + scale_y_continuous(limits = c(0, 0.20)) #retirei os outliers da figura e mudei a escala para otimizar o espa?o 
p2 <- p1 + xlab("Treatment") + ylab("Production") # alterei as legendas (n?o foi especialmente necass?rio aqui)
p3 <- p2 + theme_bw() + theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank()) # retirei as linhas de grade e peinel de fundo cinza
p4 <- p3 + theme(axis.title = element_text(size = rel(1.5))) + theme(axis.text = element_text(size = rel(1.0))) #mudei o tamanho das escalas e legendas
p4

RBAtol <- filter (meta, response %in% c("Respiration"), Lake %in% c("Atol"))
lmeRBAtol <-lme(value~(N*P), random=~1|Temperature/Block, data=RBAtol)
summary(lmeRBAtol)
p1 <- ggplot(RBAtol, aes(x = Treatment, y = value)) + geom_boxplot() + scale_y_continuous(limits = c(0, 0.60)) #retirei os outliers da figura e mudei a escala para otimizar o espa?o 
p2 <- p1 + xlab("Treatment") + ylab("Respiration") # alterei as legendas (n?o foi especialmente necass?rio aqui)
p3 <- p2 + theme_bw() + theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank()) # retirei as linhas de grade e peinel de fundo cinza
p4 <- p3 + theme(axis.title = element_text(size = rel(1.5))) + theme(axis.text = element_text(size = rel(1.0))) #mudei o tamanho das escalas e legendas
p4

BGEAtol <- filter (meta, response %in% c("BGE"), Lake %in% c("Atol"))
lmeBGEAtol <-lme(value~(N*P), random=~1|Temperature/Block, data=BGEAtol)
summary(lmeBGEAtol)
p1 <- ggplot(BGEAtol, aes(x = Treatment, y = value)) + geom_boxplot() + scale_y_continuous(limits = c(0, 0.60)) #retirei os outliers da figura e mudei a escala para otimizar o espa?o 
p2 <- p1 + xlab("Treatment") + ylab("BGE") # alterei as legendas (n?o foi especialmente necass?rio aqui)
p3 <- p2 + theme_bw() + theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank()) # retirei as linhas de grade e peinel de fundo cinza
p4 <- p3 + theme(axis.title = element_text(size = rel(1.5))) + theme(axis.text = element_text(size = rel(1.0))) #mudei o tamanho das escalas e legendas
p4

## Cabiunas e Carapebus apresentaram composrtamento muito semelhante e diferente de Comprida, Amarra-Boi e Atoleiro, que, por sua vez, apresentaram comportamento parecido.
### Assim, reanalizarei juntando as lagoas "claras" (Cab e Carap) e "escuras" (Comp, Aboi e Atol) 

#escuras
PBdark <- filter (meta, response %in% c("Production"), Lake %in% c("Comp", "Aboi", "Atol"))
lmePBdark <-lme(value~(N*P), random=~1|Temperature*Lake/Block, data=PBdark)
summary(lmePBdark)
p1 <- ggplot(PBdark, aes(x = Treatment, y = value)) + geom_boxplot(aes(fill = Lake)) + scale_y_continuous(limits = c(0, 0.15)) + scale_x_discrete("Treatment", labels = c("B" = "Control")) #retirei os outliers da figura e mudei a escala para otimizar o espa?o 
p2 <- p1 + xlab("Treatment") + ylab("Production") + xlim("B","N","P","NP") # alterei as legendas (n?o foi especialmente necass?rio aqui)
p3 <- p2 + theme_bw() + theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank()) # retirei as linhas de grade e peinel de fundo cinza
p4 <- p3 + theme(axis.title = element_text(size = rel(1.5))) + theme(axis.text = element_text(size = rel(1.0))) #mudei o tamanho das escalas e legendas
p4
anova(lmePBdark)

RBdark <- filter (meta, response %in% c("Respiration"), Lake %in% c("Comp", "Aboi", "Atol"))
lmeRBdark <-lme(value~(N*P), random=~1|Temperature*Lake/Block, data=RBdark)
summary(lmeRBdark)
p1 <- ggplot(RBdark, aes(x = Treatment, y = value)) + geom_boxplot(aes(fill = Lake)) + scale_y_continuous(limits = c(0, 0.60)) #retirei os outliers da figura e mudei a escala para otimizar o espa?o 
p2 <- p1 + xlab("Treatment") + ylab("Respiration") #+ xlim("B","N","P","NP") # alterei as legendas (n?o foi especialmente necass?rio aqui)
p3 <- p2 + theme_bw() + theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank()) # retirei as linhas de grade e peinel de fundo cinza
p4 <- p3 + theme(axis.title = element_text(size = rel(1.5))) + theme(axis.text = element_text(size = rel(1.0))) #mudei o tamanho das escalas e legendas
p4
anova(lmeRBdark)

BGEdark <- filter (meta, response %in% c("BGE"), Lake %in% c("Comp", "Aboi", "Atol"))
lmeBGEdark <-lme(value~(N*P), random=~1|Temperature*Lake/Block, data=BGEdark)
summary(lmeBGEdark)
p1 <- ggplot(BGEdark, aes(x = Treatment, y = value)) + geom_boxplot(aes(fill = Lake)) + scale_y_continuous(limits = c(0, 0.50)) #retirei os outliers da figura e mudei a escala para otimizar o espa?o 
p2 <- p1 + xlab("Treatment") + ylab("BGE") + xlim("B","N","P","NP") # alterei as legendas (n?o foi especialmente necass?rio aqui)
p3 <- p2 + theme_bw() + theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank()) # retirei as linhas de grade e peinel de fundo cinza
p4 <- p3 + theme(axis.title = element_text(size = rel(1.5))) + theme(axis.text = element_text(size = rel(1.0))) #mudei o tamanho das escalas e legendas
p4
anova(lmeBGEdark)

#claras
PBclear <- filter (meta, response %in% c("Production"), Lake %in% c("Cab", "Carap"))
lmePBclear <-lme(value~(N*P), random=~1|Temperature*Lake/Block, data=PBclear)
summary(lmePBclear)
p1 <- ggplot(PBclear, aes(x = Treatment, y = value)) + geom_boxplot(aes(fill = Lake)) + scale_y_continuous(limits = c(0, 0.3)) #retirei os outliers da figura e mudei a escala para otimizar o espa?o 
p2 <- p1 + xlab("Treatment") + ylab("Production") + xlim("B","N","P","NP") # alterei as legendas (n?o foi especialmente necass?rio aqui)
p3 <- p2 + theme_bw() + theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank()) # retirei as linhas de grade e peinel de fundo cinza
p4 <- p3 + theme(axis.title = element_text(size = rel(1.5))) + theme(axis.text = element_text(size = rel(1.0))) #mudei o tamanho das escalas e legendas
p4
anova(lmePBclear)


RBclear <- filter (meta, response %in% c("Respiration"), Lake %in% c("Cab", "Carap"))
lmeRBclear <-lme(value~(N*P), random=~1|Temperature/Block, data=RBclear)
summary(lmeRBclear)
p1 <- ggplot(RBclear, aes(x = Treatment, y = value)) + geom_boxplot(aes(fill = Lake)) + scale_y_continuous(limits = c(0, 1)) #retirei os outliers da figura e mudei a escala para otimizar o espa?o 
p2 <- p1 + xlab("Treatment") + ylab("Production") #+ xlim("B","N","P","NP") # alterei as legendas (n?o foi especialmente necass?rio aqui)
p3 <- p2 + theme_bw() + theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank()) # retirei as linhas de grade e peinel de fundo cinza
p4 <- p3 + theme(axis.title = element_text(size = rel(1.5))) + theme(axis.text = element_text(size = rel(1.0))) #mudei o tamanho das escalas e legendas
p4
anova(lmeRBclear)

BGEclear <- filter (meta, response %in% c("BGE"), Lake %in% c("Cab", "Carap"))
lmeBGEclear <-lme(value~(N*P), random=~1|Temperature/Block, data=BGEclear)
summary(lmeBGEclear)
p1 <- ggplot(BGEclear, aes(x = Treatment, y = value)) + geom_boxplot(aes(fill = Lake)) + scale_y_continuous(limits = c(0, 0.60)) #retirei os outliers da figura e mudei a escala para otimizar o espa?o 
p2 <- p1 + xlab("Treatment") + ylab("BGE") + xlim("Control","N","P","NP") # alterei as legendas (n?o foi especialmente necass?rio aqui)
p3 <- p2 + theme_bw() + theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank()) # retirei as linhas de grade e peinel de fundo cinza
p4 <- p3 + theme(axis.title = element_text(size = rel(1.5))) + theme(axis.text = element_text(size = rel(1.0))) #mudei o tamanho das escalas e legendas
p4
anova(lmeBGEclear)

# finalmente, juntados todos os plots (claras e escuras) em uma Figura apenas para publica??o

ggplot(meta, aes(x = Treatment, y = value)) + 
  geom_boxplot(aes(fill = Lake)) + #scale_y_continuous(limits = c(0, 0.6)) + #retirei os outliers da figura e mudei a escala para otimizar o espa?o 
  xlab("Treatment") + ylab("Response") + # alterei as legendas (n?o foi especialmente necass?rio aqui)
  theme_bw() + theme(panel.grid.major.y = element_blank(), 
                     panel.grid.minor.y = element_blank(), 
                     panel.grid.major.x = element_blank()) + # retirei as linhas de grade e peinel de fundo cinza
  theme(axis.title = element_text(size = rel(1.5))) + 
  theme(axis.text = element_text(size = rel(1.0))) + #mudei o tamanho das escalas e legendas
  facet_grid(response ~ clarity, scales = "free_y")
### vendo como fica 
