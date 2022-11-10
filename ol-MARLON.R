if (!require('nortest')) install.packages('nortest')
if(!require('lawstar')) install.packages('lawstar')
if(!require('MASS')) install.packages('MASS')
if(!require('dplyr')) install.packages('dplyr')
if(!require('reshape')) install.packages('reshape')
if(!require('lawstar')) install.packages('lawstar')
if(!require('nlme')) install.packages('nlme')
if(!require('gamlss')) install.packages('gamlss')
if(!require('esquisse')) install.packages('esquisse')
if(!require('hnp')) install.packages('hnp')
if(!require('lattice')) install.packages('lattice')
if(!require('ggplot2')) install.packages('ggplot2')

# Apontamento do diretório de trabalho
setwd('/home/marlon-rogerio/apps/tese-guiomar')

# Carregamento do OLLST no projeto
source.with.encoding("OLLST-gamlss.R", encoding = 'UTF-8')

# Importacao do dados
dados = read.table('r1.txt', h=T)
dados <- data.frame(dados)

attach(dados)

dados <- dados[order(dados$EP, decreasing=F), ] 
names(dados)
# No caso, ordena de acordo por w em ordem alfabética inversa, x em ordem alfabética e y em ordem decrescente. Se quiser que tudo seja na mesma ordem, basta um TRUE ou FALSE, que valerá para todos.

# Mostragem dos dados
View(dados)

# Verificar a normalidade da variável resposta no conjunto de dados
lillie.test(dados$RESP)


# Geração do boxplot
ggplot(dados) +
  aes(x = "", y = RESP, group = ESP) +
  geom_boxplot(fill = "#112446") +
  theme_bw() +
  facet_wrap(vars(EP))
# Eh possivel inferir uma dependência dos dados pelo boxplot


histDist(dados$RESP, family = "OLLST", xlab = 'Resposta', ylab = 'Frequência')
# Presença de bimodalidade

# Tratamento das vaiáveis de trabalho
especie = as.factor(dados$ESP)
length(especie)
epoca = as.factor(dados$EP)
length(epoca)
resposta = dados$RESP
resposta
length(resposta)

modelo <-gamlss::gamlss(resposta~re(fixed=~especie+epoca:especie, random=~1|epoca), data=dados, family = "OLLST")

summary(modelo)
coef(getSmo(modelo))
ranef(getSmo(modelo))
VarCorr(getSmo(modelo))
summary(getSmo(modelo))
intervals(getSmo(modelo))
fitted(getSmo(modelo))
fixef(getSmo(modelo))

r = residuals(modelo)
my.hnp <- hnp(r,halfnormal = F, print.on=TRUE, plot=FALSE)
plot(my.hnp, main="(a2)", xlab="Half-ollst scores",
     ylab="Resíduos de quantis", legpos="topleft")






############# ----------------------------
library(nlme)
require(esquisse)
require(hnp)
require(gamlss)
library(dplyr)
library(tidyr)
library(ggplot2)

# Para o conjunto de dados r1.txt, considerar as seguintes variáveis:
#   ESP:    relativo a espécie do feijão (== TRATAMENTO)
#   EP:     relativo a época (== TEMPO)
#   RESP:   relativo a resposta (== RESPOSTA)

# Apontamento para o diretório de trabalho e carregamento dos dados
setwd('/home/marlon-rogerio/apps/tese-guiomar')

# Carregamento da distribuição de tratamento dos dados
source.with.encoding("OLLST-gamlss.R", encoding = 'UTF-8')
source.with.encoding("OLLSN-gamlss.R", encoding = 'UTF-8')

# Carregamento dos dados a serem consulmidos
dados = read.table("r5.txt", h = T)
shapiro.test(dados$RESP)
# W = 0.86129, p-value = 1.754e-08

View(dados)


dados %>%
  mutate(id = row_number()) %>%
  gather(hora, y, -id) %>%
  mutate(hora = factor(hora, levels = c("0", "1", "2", "3", "4", "6"))) %>%
  ggplot(aes(x = hora, y = y, group = id)) + geom_line() + geom_point()

# Acessa as variáveis presentes na estrutura dos dados sem referenciar diretamente os dados
attach(dados)

y <- dados$RESP
y
tratamento <- as.factor(dados$ESP)
tratamento
subject <- as.factor(rep(c(1:9),rep(9,9)))
subject
hora<- as.factor(dados$EP)
dados <- data.frame(y,tratamento,subject,hora)

str(dados)

# Gráfico de perfis 
ggplot(dados, aes(x=hora, y=y, group=1)) +
  geom_line() + 
  geom_point() + 
  labs(x="Hora", y="Y")

View(dados)
attach(dados)
# Ajuste do modelo sem interação
modelo_sem_interacao <- gamlss(y ~ hora+hora:e_tratamento, family = "LO", data = dados, n.cyc=1000)
summary(modelo_sem_interacao)
summary(getSmo(modelo_sem_interacao))

# Ajuste do modelo com interação
modelo_com_interacao <- gamlss(y ~ hora*tratamento, family = "OLLSN", data = dados)
summary(modelo_com_interacao)

str(modelo_com_interacao)

# Verificação pela verossimilhança se a interação foi siginificativa ou não 
tcal <- modelo_sem_interacao$P.deviance - modelo_com_interacao$P.deviance
tcal
valorp <- 1-pchisq(tcal, 27)
valorp
p_valor <- LR.test(modelo_sem_interacao, modelo_com_interacao)
p_valor 

# Ajuste do modelo hieráquico com efeito aleatório do intercepto
tratamento
e_tratamento = relevel(tratamento,"GBS22") # Comparação em relação ao tratamento 1
# e_tratamento = relevel(subject,1) # Comparação em relação ao tratamento 1
e_tratamento
hora
#t3 <- gamlss(C ~ re(fixed = ~ trat+epoca:trat, random=~1|trat), data = dados, family = "OLLSN")
t3 <- gamlss(y ~ re(fixed = ~ hora+hora:e_tratamento, random = ~ 1|subject), data = dados, family = "OLLST", n.cyc=1000)
t3 <- gamlss(y ~ re(fixed = ~ hora+hora:e_tratamento, random = ~ 1|subject), data = dados, family = "GU", n.cyc=1000)
summary(t3)
summary(getSmo(t3))



