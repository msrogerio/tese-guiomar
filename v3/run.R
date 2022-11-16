library(nlme)
require(hnp)
require(gamlss)
library(dplyr)
library(tidyr)
library(ggplot2)
library(httr)
library(sqldf)
library(stringr)



#   Para o conjunto de dados r1.txt, considerar as seguintes variáveis:
#   ESP:    relativo a espécie do feijão (== TRATAMENTO)
#   EP:     relativo a época (== TEMPO)
#   RESP:   relativo a resposta (== RESPOSTA)
#   Apontamento para o diretório de trabalho e carregamento dos dados
setwd('/home/marlon-rogerio/apps/tese-guiomar/v3')
dados = read.table("v3.txt", h = T)

shapiro.test(dados$RESP) # Teste de normailidade
# View(dados)


#   FUNÇÃO PARA CRIAR O CONTEXTO grupo
#   Ler e calcula a quantidade de tratamentos, verifica se há ou não 
#   dissonância e cria uma variável do tipo grupo, ou seja, um contexto 
#   de agrupamento dos tratamentos.
cria.grupo <- function(dados) {
  base.grupo <- data.frame(dados%>%group_by(ESP)%>%count())
  comparador <- mean(base.grupo$n)
  grupo <- comparador
  for (i in base.grupo$n) {
    if (comparador != i) {
      grupo <- NULL
      break
    }
  }
  temp <- length(dados$ESP)/comparador
  grupo <- as.factor(rep(c(1:temp),rep(comparador,temp)))
  return(grupo)
}


#   Função para quebra as string e recuperar o valor do tempo 
#   e tratamento comparado
#
#   : return novas.colunas
#       str(novas.colunas)
#       data.frame':	2 variables:
#       $ tratamentos.comparado
#       $ tempos
tratamentos.tempos <- function(string_bruta_comparacoes){
  # É necessário correr toda a lista de comparações e segmentar a string 
  # em colunas para que fique vizivelmente melhor
  
  # A primeira coluna dos frames retornados são chamadas observações. Não há cabeçalho que idetifique-as 
  # nesse caso é necessário capturar os dados capturar os nomes de cada linha (obsevações) para que 
  # possamos quebrar a string de dados a função row.names faz essa captura. 
  cont <- 1
  tempos <- NULL
  tratamentos.comparado <- NULL
  temp <- NULL
  while(cont <= length(string_bruta_comparacoes)){
    # print(string_bruta_comparacoes[cont])
    temp <- NULL
    temp <- str_split(string_bruta_comparacoes[cont], ':tratamento.removido', simplify = TRUE)
    # print(temp)
    # Verifica se e a variável está vazia. Se esse for o caso,
    # o primero valor é inserido
    if (is.null(tratamentos.comparado)) {
      tratamentos.comparado <- c(temp[,2])
    } else {
      # Se já houver algum valor na variável um novo valor será acrescentado
      tratamentos.comparado <- append(tratamentos.comparado, temp[,2])
    }
    if (is.null(tempos)) {
      tempos <- str_split(temp, 'epoca', simplify = TRUE)[1,2]
      string_bruta_comparacoes
    } else {
      tempos <- append(tempos, c(str_split(temp, 'epoca', simplify = TRUE)[1,2]))
    }
    cont <- cont + 1
  }
  novas.colunas <- data.frame(tratamentos.comparado, tempos)
  return(novas.colunas)
}


#   Função para verificar a significancia dos dados e cria uma nova coluna 
#   para ser posta no novo data frame
verifica.significancia <- function(p.valor){
  significancia <- NULL
  for (j in p.valor) {
    if (is.null(significancia)){
      if (j < 0.5) {
        significancia <- c("*")
      } else {
        significancia <- c("NS")
      }
    }else {
      if (j < 0.5) {
        significancia <- append(significancia, "*")
      } else {
        significancia <- append(significancia, "NS")
      }
    }
  }
  return (significancia)
}


#   Função para tramento da tabela resultado do ajuste do modelo 
ajuste.tabela <- function(t.table, tratamento.evidencia, dados) {
  
  # Criar um vetor de níveis
  # ! Importante
  # Os primeiros valores mostados no tTable fazem referência oa intercepto
  # Para esse tipo de estudo, não há interesse em tais valores. 
  # Considerando que as linhas do intercpto são em quantidades iguais ao tamanho
  # do vetor de níveis da unidade temporal do experimento, nesse caso, para 
  # retirar a saída do intercepto do tTable, basta ler tal vetor e deduzir 
  # a quantidade de linhas no tTable
  
  niveis.tempo <- levels(dados$epoca)
  
  string_bruta_comparacoes <- row.names(t.table[,0])
  novas.colunas <- tratamentos.tempos(string_bruta_comparacoes)
  tempos <- novas.colunas$tempos
  tratamentos.comparado <- novas.colunas$tratamentos.comparado
  tratamento.evidencia <- rep(tratamento.evidencia, length(tempos))
  value <- t.table[,1]
  Std.Error<- t.table[,2]
  DF <- t.table[,3]
  t.value <- t.table[,4]
  p.value<- t.table[,5]

  significancia <- verifica.significancia(p.value)
  
  novo_frame <- data.frame(
    tempos, tratamento.evidencia, tratamentos.comparado,
    value, Std.Error, DF, t.value, p.value, significancia
  )
}


#   Bloco de variáveis
y <- dados$RESP
grupo <- cria.grupo(dados)
tratamento <- as.factor(dados$ESP)
epoca <- as.factor(dados$EP)
dados <- data.frame(y,tratamento,grupo,epoca)
lista.tratamentos <- levels(tratamento)


ggplot(dados, aes(x=epoca, y=y, group=1)) +
  geom_line() + 
  geom_point() + 
  labs(x="ÉPOCA", y="Y")

ggsave('imagem.png')

tabela.a <- NULL
tabela.b <- NULL
tabela.c <- NULL
cont <- 1

#  Persistência do ajuste para cada tratamento 
for (i in lista.tratamentos) {
  
  tratamento.removido <- relevel(dados$tratamento, i)
  ajuste <- gamlss(y ~ re(fixed = ~ epoca+epoca:tratamento.removido, random = ~ 1|grupo), data = dados, family = "RG", n.cyc=100)
  return <- summary(getSmo(ajuste))
  tabela.t <- return$tTable
  niveis.tempo <- levels(epoca)
  tabela.t <- slice(data.frame(tabela.t), -(1:length(niveis.tempo)))
  
  tabela.a <- ajuste.tabela(tabela.t, i, dados)
  nome.arquivo <- sprintf("%s/result/%s.csv", getwd(), cont)
  write.table(tabela.a, file = nome.arquivo, sep=",", na="", quote=TRUE, row.names=FALSE)
  cont <- cont +1
}


