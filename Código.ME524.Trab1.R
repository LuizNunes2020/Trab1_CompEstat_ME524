# Pacotes
library(tidyverse)
library(ggplot2)
library(coda)
library(GGally)
library(reshape2)


# Parâmetros iniciais
set.seed(123456)
lambda <- 30
n_iter <- 10000  # Número de iterações do Gibbs
a <- 2  # Parâmetro a da Beta
b <- 8  # Parâmetro b da Beta

# Criando a matriz de capturas (exemplo com 12 indivíduos e 14 sessões de captura)
capturas <- matrix(c(
  0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1,
  1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1,
  0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1,
  1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1,
  0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
  0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0,
  0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 1,
  0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0
), nrow = 12, byrow = TRUE)

# Calculando o número de capturas em cada sessão (n_i)
n <- colSums(capturas)
n

# Calculando o número de recapturas (m_j) para as sessões subsequentes
m <- numeric(length(n) - 1)
for (j in 2:length(n)) {
  if (j == 2) {
    m[j - 1] <- sum(capturas[, 1] > 0 & capturas[, j] == 1)
  } else {
    m[j - 1] <- sum(rowSums(capturas[, 1:(j - 1)]) > 0 & capturas[, j] == 1)
  }
}
m

# r é o número total de capturas menos recapturas
r <- sum(n) - sum(m)
r

# Função de log-verossimilhança ajustada para garantir que N seja um valor escalar
log_vero_N <- function(N, p, n, r, lambda) {
  if (!is.numeric(N) || length(N) != 1 || N < r || is.na(N)) return(-Inf)
  
  log_prior <- (N * log(lambda) - lambda)  # Componente do log da distribuição Poisson
  log_factorial_terms <- lgamma(N + 1) - lgamma(N - r + 1)  # lgamma é usada para N!
  log_likelihood <- sum((N - n) * log(1 - p))  # Componente do log da verossimilhança
  
  return(log_prior + log_factorial_terms + log_likelihood)
}
# Função para amostrar os p_i condicionalmente a N
sample_p <- function(N, n) {
  p <- numeric(length(n))
  for (i in 1:length(n)) {
    p[i] <- rbeta(1, n[i] + 1, N - n[i] + 1)
  }
  return(p)
}

# Implementação do Metropolis-Hastings para amostrar N
mh_N <- function(p, n, r, lambda, N_init, n_iter = 1) {
  N_samples <- numeric(n_iter)
  N_samples[1] <- N_init
  
  for (t in 2:n_iter) {
    N_current <- max(N_samples[t - 1], r)  # Garante que N_current seja pelo menos r
    N_proposal <- rpois(1, lambda = N_current)  # Proposta a partir de Poisson centrada em N^(t-1)
    N_proposal <- max(N_proposal, r)  # Garantir que N >= r
    
    # Calcula o log de alpha para a aceitação
    log_alpha <- log_vero_N(N_proposal, p, n, r, lambda) - log_vero_N(N_current, p, n, r, lambda)
    alpha <- exp(log_alpha)
    
    if (runif(1) < alpha) {
      N_samples[t] <- N_proposal
    } else {
      N_samples[t] <- N_current
    }
  }
  
  return(N_samples[n_iter])  # Retorna apenas o último valor amostrado
}

# Implementação do Amostrador de Gibbs
gibbs_sampler <- function(n_iter = 1000, lambda = 30, n, r, N_init = 12) {
  N_samples <- numeric(n_iter)
  p_samples <- matrix(0, nrow = n_iter, ncol = length(n))
  
  # Inicializações
  N_samples[1] <- N_init
  p_samples[1, ] <- runif(length(n), 0, 1)
  
  # Iterações do Gibbs
  for (i in 2:n_iter) {
    # Amostragem de N usando o Metropolis-Hastings
    N_samples[i] <- mh_N(p_samples[i - 1, ], n, r, lambda, N_samples[i - 1])
    
    # Amostragem de p_i para cada sessão de captura
    for (j in 1:length(n)) {
      p_samples[i, j] <- rbeta(1, n[j] + 1, N_samples[i] - n[j] + 1)
    }
  }
  
  list(N_samples = N_samples, p_samples = p_samples)
}

# Testando o algoritmo do Amostrador de Gibbs
resultados_gibbs_caso1 <- gibbs_sampler(n_iter = n_iter, lambda = lambda, n = n, r = r, N_init = 12)

# Separar as cadeias para o cálculo de R-hat
n_chains <- 2
n_iter_per_chain <- length(resultados_gibbs_caso1$N_samples) / n_chains

# Dividir as amostras em duas cadeias
N_samples_chain1 <- resultados_gibbs_caso1$N_samples[1:n_iter_per_chain]
N_samples_chain2 <- resultados_gibbs_caso1$N_samples[(n_iter_per_chain + 1):(2 * n_iter_per_chain)]

# Criar objetos MCMC para as duas cadeias
mcmc_chain1 <- mcmc(N_samples_chain1)
mcmc_chain2 <- mcmc(N_samples_chain2)

# Combinar em um objeto mcmc.list para calcular R-hat
mcmc_list <- mcmc.list(mcmc_chain1, mcmc_chain2)

# Calcular R-hat para N
rhat_N <- gelman.diag(mcmc_list)$psrf[1]
print(paste("R-hat for N:", rhat_N))

# Trace plot para N 
ggplot(data.frame(iter = 1:length(resultados_gibbs_caso1$N_samples), N = resultados_gibbs_caso1$N_samples), aes(x = iter, y = N)) +
  geom_line(color = "blue") +
  theme_minimal() +
  labs(title = "Trace Plot para N (Amostrador de Gibbs)",
       x = "Iterações",
       y = "N") +
  theme(plot.title = element_text(hjust = 0.5))

# Histograma da distribuição posterior de N com destaque em N = 12
ggplot(data.frame(N_samples = resultados_gibbs_caso1$N_samples), aes(x = N_samples)) +
  geom_histogram(binwidth = 1, fill = "blue", color = "black", alpha = 0.7) +
  geom_vline(xintercept = 12, color = "red", linetype = "dashed", size = 1) +
  theme_minimal() +
  labs(title = "Distribuição Posterior de N",
       x = "N",
       y = "Frequência") +
  theme(plot.title = element_text(hjust = 0.5))


# Calcular R-hat para cada p_i e gerar trace plots
rhat_p <- numeric(ncol(resultados_gibbs_caso1$p_samples))
for (i in 1:ncol(resultados_gibbs_caso1$p_samples)) {
  p_samples_chain1 <- resultados_gibbs_caso1$p_samples[1:n_iter_per_chain, i]
  p_samples_chain2 <- resultados_gibbs_caso1$p_samples[(n_iter_per_chain + 1):(2 * n_iter_per_chain), i]
  
  mcmc_list_p <- mcmc.list(mcmc(p_samples_chain1), mcmc(p_samples_chain2))
  rhat_p[i] <- gelman.diag(mcmc_list_p)$psrf[1]
  
  # Trace plot para p_i
  traceplot(mcmc_list_p, main = paste("Trace Plot for p", i))
  
  # Histograma para p_i
  ggplot(data.frame(p_samples = resultados_gibbs_caso1$p_samples[, i]), aes(x = p_samples)) +
    geom_histogram(binwidth = 0.05, fill = "blue", color = "black", alpha = 0.7) +
    theme_minimal() +
    labs(title = paste("Posterior Distribution of p", i), x = paste("p", i), y = "Frequency")
}

# Mostrar os valores de R-hat para cada p_i
print(rhat_p)

# Criar um data frame combinando todos os p_i para uma densidade conjunta
df_p <- data.frame(iter = rep(1:n_iter, ncol(resultados_gibbs_caso1$p_samples)),
                   p = as.vector(resultados_gibbs_caso1$p_samples),
                   variable = rep(paste0("p", 1:ncol(resultados_gibbs_caso1$p_samples)), each = n_iter))

# Gráfico de densidade combinado para todos os p_i
ggplot(df_p, aes(x = p, color = variable)) +
  geom_density() +
  theme_minimal() +
  labs(title = "Distribuições de Densidade dos p_i",
       x = "Valor de p",
       y = "Densidade") +
  theme(legend.position = "right")


# Supondo que os dados das amostras de p_i do Caso 1 estejam no objeto `resultados_gibbs_caso1$p_samples`
p_samples_caso1 <- as.data.frame(resultados_gibbs_caso1$p_samples)

# Ajustar os nomes das colunas para representar os p_i
colnames(p_samples_caso1) <- paste0("p", 1:ncol(p_samples_caso1))

# Calcular a matriz de correlação
cor_matrix_caso1 <- cor(p_samples_caso1)

# Transformar a matriz de correlação em um formato longo para o ggplot
cor_data_caso1 <- melt(cor_matrix_caso1)

# Criar o heatmap usando ggplot2
ggplot(cor_data_caso1, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  labs(title = "Heatmap de Correlação dos p_i (Caso 1)",
       x = "p_i",
       y = "p_i",
       fill = "Correlação") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Estimativa do tamanho populacional (a média da amostra posterior de N)
N_estimated <- mean(resultados_gibbs_caso1$N_samples)

# Área efetivamente amostrada (exemplo: 524 km² conforme o estudo de Silveira)
A <- 524

# Densidade de onças por 100 km²
densidade <- (N_estimated / A) * 100

print(paste("Densidade estimada de onças-pintadas:", round(densidade, 2), "indivíduos por 100 km²"))

















# Função para amostrar os p_i condicionalmente a N, considerando a priori Beta(a, b)
sample_p <- function(N, n, a, b) {
  p <- numeric(length(n))
  for (i in 1:length(n)) {
    p[i] <- rbeta(1, a + n[i], b + N - n[i])
  }
  return(p)
}
# Implementação do Amostrador de Gibbs para o Caso 2
gibbs_sampler <- function(n_iter = 10000, lambda = 30, n, r, a, b, N_init = 12) {
  N_samples <- numeric(n_iter)
  p_samples <- matrix(0, nrow = n_iter, ncol = length(n))
  
  # Inicializações
  N_samples[1] <- N_init
  p_samples[1, ] <- runif(length(n), 0, 1)
  
  # Iterações do Gibbs
  for (i in 2:n_iter) {
    # Amostragem de N usando o Metropolis-Hastings
    N_samples[i] <- mh_N(p_samples[i - 1, ], n, r, lambda, N_samples[i - 1])
    
    # Amostragem de p_i para cada sessão de captura
    p_samples[i, ] <- sample_p(N_samples[i], n, a, b)
  }
  
  list(N_samples = N_samples, p_samples = p_samples)
}

# Testando o algoritmo do Amostrador de Gibbs para o Caso 2
resultados_gibbs_caso2 <- gibbs_sampler(n_iter = n_iter, lambda = lambda, n = n, r = r, a = a, b = b, N_init = 12)

# Separar as cadeias para o cálculo de R-hat
n_chains <- 2
n_iter_per_chain <- length(resultados_gibbs_caso2$N_samples) / n_chains

# Dividir as amostras em duas cadeias
N_samples_chain1 <- resultados_gibbs_caso2$N_samples[1:n_iter_per_chain]
N_samples_chain2 <- resultados_gibbs_caso2$N_samples[(n_iter_per_chain + 1):(2 * n_iter_per_chain)]

# Criar objetos MCMC para as duas cadeias
mcmc_chain1 <- mcmc(N_samples_chain1)
mcmc_chain2 <- mcmc(N_samples_chain2)

# Combinar em um objeto mcmc.list para calcular R-hat
mcmc_list <- mcmc.list(mcmc_chain1, mcmc_chain2)

# Calcular R-hat para N
rhat_N <- gelman.diag(mcmc_list)$psrf[1]
print(paste("R-hat for N:", rhat_N))

# Trace plot para N
ggplot(data.frame(iter = 1:length(resultados_gibbs_caso2$N_samples), N = resultados_gibbs_caso2$N_samples), aes(x = iter, y = N)) +
  geom_line(color = "blue") +
  theme_minimal() +
  labs(title = "Trace Plot para N (Amostrador de Gibbs)",
       x = "Iterações",
       y = "N") +
  theme(plot.title = element_text(hjust = 0.5))

# Histograma para N
ggplot(data.frame(N_samples = resultados_gibbs_caso2$N_samples), aes(x = N_samples)) +
  geom_histogram(binwidth = 1, fill = "blue", color = "black", alpha = 0.7) +
  theme_minimal() +
  labs(title = "Posterior Distribution of N (Caso 2)", x = "N", y = "Frequency")

# Calcular R-hat para cada p_i e gerar trace plots
rhat_p <- numeric(ncol(resultados_gibbs_caso2$p_samples))
for (i in 1:ncol(resultados_gibbs_caso2$p_samples)) {
  p_samples_chain1 <- resultados_gibbs_caso2$p_samples[1:n_iter_per_chain, i]
  p_samples_chain2 <- resultados_gibbs_caso2$p_samples[(n_iter_per_chain + 1):(2 * n_iter_per_chain), i]
  
  mcmc_list_p <- mcmc.list(mcmc(p_samples_chain1), mcmc(p_samples_chain2))
  rhat_p[i] <- gelman.diag(mcmc_list_p)$psrf[1]
  
  # Trace plot para p_i
  traceplot(mcmc_list_p, main = paste("Trace Plot for p", i, "(Caso 2)"))
  
  # Histograma para p_i
  ggplot(data.frame(p_samples = resultados_gibbs_caso2$p_samples[, i]), aes(x = p_samples)) +
    geom_histogram(binwidth = 0.05, fill = "blue", color = "black", alpha = 0.7) +
    theme_minimal() +
    labs(title = paste("Posterior Distribution of p", i, "(Caso 2)"), x = paste("p", i), y = "Frequency")
}

# Mostrar os valores de R-hat para cada p_i
print(rhat_p)



# Supondo que os dados das amostras de p_i do Caso 2 estejam no objeto `resultados_gibbs_caso2$p_samples`
p_samples_caso2 <- as.data.frame(resultados_gibbs_caso2$p_samples)

# Ajustar os nomes das colunas para representar os p_i
colnames(p_samples_caso2) <- paste0("p", 1:ncol(p_samples_caso2))

# Transformar os dados em formato longo para facilitar a criação do gráfico
p_samples_long <- p_samples_caso2 %>%
  pivot_longer(cols = everything(), names_to = "variable", values_to = "value")

# Criar o gráfico de densidade conjunto dos p_i
ggplot(p_samples_long, aes(x = value, color = variable)) +
  geom_density() +
  theme_minimal() +
  labs(title = "Distribuições de Densidade dos p_i (Caso 2)",
       x = "Valor de p",
       y = "Densidade") +
  theme(legend.title = element_blank())


# Converter as amostras de p_i para um data frame
p_samples_df <- as.data.frame(resultados_gibbs_caso2$p_samples)
colnames(p_samples_df) <- paste0("p", 1:14)  # Nomear as colunas como p1, p2, ..., p14

# Calcular a matriz de correlação entre os p_i
correlation_matrix <- cor(p_samples_df)

# Plotar a matriz de correlação como um heatmap
ggplot(melt(correlation_matrix), aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, 
                       limit = c(-1, 1), space = "Lab", 
                       name = "Correlação") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Heatmap de Correlação dos p_i (Caso 2)",
       x = "p_i",
       y = "p_i")




# Estimativa do tamanho populacional (a média da amostra posterior de N)
N_estimated_caso2 <- mean(resultados_gibbs_caso2$N_samples)

# Área efetivamente amostrada (exemplo: 524 km² conforme o estudo de Silveira)
A <- 524
# Densidade de onças por 100 km²
densidade_2 <- (N_estimated_caso2 / A) * 100

print(paste("Densidade estimada de onças-pintadas:", round(densidade_2, 2), "indivíduos por 100 km²"))











# Supondo que as amostras de N para os dois casos estejam nos seguintes objetos:
N_samples_caso1 <- resultados_gibbs_caso1$N_samples
N_samples_caso2 <- resultados_gibbs_caso2$N_samples

# Definindo o nível do HDI (por exemplo, 95%)
nivel_hdi <- 0.95

# Calculando os HDIs para N nos dois casos
hdi_caso1 <- HPDinterval(mcmc(N_samples_caso1), prob = nivel_hdi)
hdi_caso2 <- HPDinterval(mcmc(N_samples_caso2), prob = nivel_hdi)

# Imprimindo os resultados
cat("HDI para N no Caso 1 (", nivel_hdi * 100, "%):\n")
print(hdi_caso1)
cat("\nHDI para N no Caso 2 (", nivel_hdi * 100, "%):\n")
print(hdi_caso2)

# Criando um data frame com os resultados para visualização
hdi_data <- data.frame(
  Caso = c("Caso 1", "Caso 2"),
  Lower = c(hdi_caso1[1], hdi_caso2[1]),
  Upper = c(hdi_caso1[2], hdi_caso2[2])
)

# Gráfico comparativo dos HDIs
ggplot(hdi_data, aes(x = Caso, y = Lower, ymax = Upper, ymin = Lower)) +
  geom_linerange(size = 2, color = "blue") +
  geom_point(aes(y = (Lower + Upper) / 2), size = 3, color = "red") +
  theme_minimal() +
  labs(title = paste("Comparação dos HDIs de N para", nivel_hdi * 100, "% de Credibilidade"),
       x = "Caso",
       y = "N") +
  coord_flip()

# Criando uma tabela com os resultados dos HDIs
hdi_table <- data.frame(
  Caso = c("Caso 1", "Caso 2"),
  HDI_Lower = c(hdi_caso1[1], hdi_caso2[1]),
  HDI_Upper = c(hdi_caso1[2], hdi_caso2[2]),
  Intervalo = paste0("[", round(hdi_caso1[1], 2), ", ", round(hdi_caso1[2], 2), "]", collapse = " "),
  Intervalo_Caso2 = paste0("[", round(hdi_caso2[1], 2), ", ", round(hdi_caso2[2], 2), "]", collapse = " ")
)
hdi_table











boxplot(N_samples_caso1, N_samples_caso2, names = c("Caso 1", "Caso 2"),
        main = "Boxplot das Distribuições de N para os Casos",
        ylab = "N", col = c("blue", "green"))

