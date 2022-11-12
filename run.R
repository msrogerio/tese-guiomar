library(nlme)
require(esquisse)
require(hnp)
require(gamlss)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

### Biblioteca de comunicação API
if(!require('httr')) {
  install.packages('httr')
  library('httr')
}


# Para o conjunto de dados r1.txt, considerar as seguintes variáveis:
#   ESP:    relativo a espécie do feijão (== TRATAMENTO)
#   EP:     relativo a época (== TEMPO)
#   RESP:   relativo a resposta (== RESPOSTA)

# Apontamento para o diretório de trabalho e carregamento dos dados
setwd('/home/marlon-rogerio/apps/tese-guiomar')
dados = read.table("v5.txt", h = T)
shapiro.test(dados$RESP)
# W = 0.86129, p-value = 1.754e-08

View(dados)


# FUNÇÃO PARA CRIAR O CONTEXTO SUBJECT
# Ler a calcula a quantidade de tratamentos, virifica se há ou não 
# dissonância e cria uma variável do tipo subject, ou seja, um contexto 
# de agrupamento dos tratamentos.
cria.subject <- function(dados) {
  base.subject <- data.frame(dados%>%group_by(ESP)%>%count())
  comparador <- mean(base.subject$n)
  subject <- comparador
  for (i in base.subject$n) {
    if (comparador != i) {
      subject <- NULL
      break
    }
  }
  temp <- length(dados$ESP)/comparador
  subject <- as.factor(rep(c(1:temp),rep(comparador,temp)))
  return(subject)
}

subject <- cria.subject(dados)
subject

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
tratamento.evidencia <- "AGET"
e_tratamento = relevel(tratamento, tratamento_evidencia) # Comparação em relação ao tratamento 1
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
niveis_de_tempo <- levels(hora)

# Captura a saída dos 
tTable <- return$tTable ## Velores de interesse e comparação de cada tratamento dentro de cada tempo
tTable <- slice(data.frame(tTable), -(1:length(niveis_de_tempo)))

## A primeira coluna dos frames retornados são chamadas observações. Não há cabeçalho que idetifique-as 
# nesse caso é necessário capturar os dados capturar os nomes de cada linha (obsevações) para que 
# possamos quebrar a string de dados a função row.names faz essa captura. 
string_bruta_comparacoes <- row.names(tTable[,0])
# string_bruta_comparacoes <- string_bruta_comparacoes[-c(1:length(niveis_de_tempo))]
# string_bruta_comparacoes


# É necessário correr toda a lista de comparações e segmentar a string 
# em colunas para que fique vizivelmente melhor
cont <- 1
tempos <- NULL
tratamentos.comparado <- NULL
temp <- NULL

while(cont <= length(string_bruta_comparacoes)){
  # print(string_bruta_comparacoes[cont])
  temp <- NULL
  temp <- str_split(string_bruta_comparacoes[cont], ':e_tratamento', simplify = TRUE)
  ## Verifica se e a variável está vazia. Se esse for o caso,
  # o primero valor é inserido
  if (is.null(tratamentos.comparado)) {
    tratamentos.comparado <- c(temp[,2])
  } else {
    # Se já houver algum valor na variável um novo valor será acrescentado
    tratamentos.comparado <- append(tratamentos.comparado, temp[,2])
  }
  if (is.null(tempos)) {
    tempos <- str_split(temp, 'hora', simplify = TRUE)[1,2]
    string_bruta_comparacoes
  } else {
    tempos <- append(tempos, c(str_split(temp, 'hora', simplify = TRUE)[1,2]))
  }
  cont <- cont + 1
}

tratamento.evidencia <- rep(tratamento.evidencia, length(tempos))
value <- tTable[,1]
Std.Error<- tTable[,2]
DF <- tTable[,3]
t.value <- tTable[,4]
p.value<- tTable[,5]

significancia <- NULL

###
## Verificando a significância dos dados pelo p.valor
for (j in p.value) {
  if (is.null(significancia)){
    if (j < 0.5) {
      significancia <- c("*")
    } else {
      significancia <- c("SN")
    }
  }else {
    if (j < 0.5) {
      significancia <- append(significancia, "*")
    } else {
      significancia <- append(significancia, "SN")
    }
  }
}

###
## Monta uma nova tTable
novo_frame <- data.frame(
  tempos, tratamento.evidencia, tratamentos.comparado,
  value, Std.Error, DF, t.value, p.value, significancia
)
