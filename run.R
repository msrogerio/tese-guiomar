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

dados = read.table("v5.txt", h = T)
shapiro.test(dados$RESP)
# W = 0.86129, p-value = 1.754e-08

View(dados)

y <- dados$RESP
y
tratamento <- as.factor(dados$ESP)
tratamento
subject <- as.factor(rep(c(1:8),rep(27,8)))
subject
hora<- as.factor(dados$EP)
dados <- data.frame(y,tratamento,subject,hora)

str(dados)

ggplot(dados, aes(x=hora, y=y, group=1)) +
  geom_line() + 
    geom_point() + 
    labs(x="Hora", y="Y")

tratamento
e_tratamento = relevel(tratamento,"AGET") # Comparação em relação ao tratamento 1
e_tratamento
hora
t3 <- gamlss(y ~ re(fixed = ~ hora+hora:e_tratamento, random = ~ 1|subject), data = dados, family = "RG", n.cyc=1000)
summary(t3)
return <- summary(getSmo(t3))

# Criar um vetor de níveis
## ! Importante
## Os primeiros valores mostados no tTable fazem referência oa intercepto
## Para esse tipo de estudo, não há interesse em tais valores. 
## Considerando que as linhas do intercpto são em quantidades iguais ao tamanho
## do vetor de níveis da unidade temporal do experimento, nesse caso, para 
## retirar a saída do intercepto do tTable, basta ler tal vetor e deduzir 
## a quantidade de linhas no tTable
tempos <- levels(hora)

tTable <- return$tTable ## Velores de interesse e comparação de cada tratamento dentro de cada tempo 
novo_frame <- 



