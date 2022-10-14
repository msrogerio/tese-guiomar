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

modelo <-gamlss(resposta~re(fixed=~especie+epoca:especie, random=~1|Subject), data=dados, family = "OLLST")

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


