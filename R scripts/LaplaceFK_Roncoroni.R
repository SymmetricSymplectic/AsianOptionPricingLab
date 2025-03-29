#### Reproducir Roncoroni-Fusai ####

# Función para calcular la transformada de Laplace en el modelo de raíz cuadrada
ltsr <- function(spot, expiry, sg, rf, gamma) {
  lambda <- sqrt(rf^2 + 2*gamma*sg*sg)
  numerator <- 2*gamma*(exp(expiry*lambda) - 1)
  denominator <- lambda + rf + (lambda - rf)*exp(expiry*lambda)
  lt <- exp(-spot*numerator/denominator)/(gamma^2)
  return(lt)
}

# Función para realizar la inversión numérica de la transformada de Laplace
# usando el algoritmo de Abate-Whitt con la aceleración de Euler
AWSR <- function(spot, strike, expiry, sg, rf, aa = 18.4, terms = 20, extraterms = 10) {
  X <- strike*expiry
  
  # Calcular la transformada en gamma = aa/(2*strike)
  lt <- ltsr(spot, expiry, sg, rf, aa/(2*X))
  sum <- lt*exp(aa/2)/(2*X)
  
  # Aplicar el algoritmo de Euler
  k <- 1:(terms + extraterms)
  arg_real <- aa/(2*X)
  arg_imag <- pi*k/X
  
  term <- numeric(length(k))
  for(i in 1:length(k)) {
    arg <- complex(real = arg_real, imaginary = arg_imag[i])
    term[i] <- ((-1)^k[i]) * ltsr(spot, expiry, sg, rf, arg) * exp(aa/2)/X
  }
  
  csum <- sum + cumsum(term)
  sumr <- Re(csum[terms:(terms+extraterms)])
  
  # Calcular coeficientes binomiales para la extrapolación
  j <- 0:extraterms
  bincoeff <- factorial(extraterms)/(factorial(j)*factorial(extraterms-j))
  
  # Resultado extrapolado
  euler <- sum(bincoeff*sumr) * (2^(-extraterms))
  
  # Aplicar la fórmula final
  optprice <- exp(-rf*expiry) * (euler + spot*(exp(rf*expiry)-1)/rf - X)/expiry
  
  return(optprice)
}

# Función para generar la tabla de precios de opciones asiáticas
generate_asian_option_table <- function() {
  # Parámetros
  spot <- 1  # El spot parece ser 1 en los ejemplos
  expiry <- 1 # Asumimos vencimiento de 1 año
  rf <- 0.05 # Tasa libre de riesgo (r en el texto)
  
  # Valores de volatilidad (sigma) y strikes a evaluar
  sigmas <- c(0.1, 0.3, 0.5)
  strikes <- c(0.9, 0.95, 1.0, 1.05, 1.1)
  
  # Crear tabla para almacenar resultados
  results <- matrix(0, nrow = length(strikes), ncol = length(sigmas) + 1)
  colnames(results) <- c("K", paste("σ =", sigmas))
  results[,1] <- strikes
  
  # Calcular precios de opciones para cada combinación
  for (i in 1:length(strikes)) {
    for (j in 1:length(sigmas)) {
      results[i, j+1] <- AWSR(spot, strikes[i], expiry, sigmas[j], rf)
    }
  }
  
  # Convertir a data frame y formatear
  results_df <- as.data.frame(results)
  
  # Imprimir tabla de resultados
  print(results_df, digits = 5)
  
  return(results_df)
}

# Ejecutar la función para generar la tabla
asian_option_prices <- generate_asian_option_table()


# Crear una función para comparar los resultados con los valores de referencia
compare_asian_option_prices <- function() {
  # Parámetros
  spot <- 1
  expiry <- 1
  rf <- 0.05
  
  # Valores de volatilidad (sigma) y strikes a evaluar
  sigmas <- c(0.1, 0.3, 0.5)
  strikes <- c(0.9, 0.95, 1.0, 1.05, 1.1)
  
  # Crear matriz para los resultados calculados
  calculated_results <- matrix(0, nrow = length(strikes), ncol = length(sigmas) + 1)
  colnames(calculated_results) <- c("K", paste("σ =", sigmas))
  calculated_results[,1] <- strikes
  
  # Valores de referencia de la tabla LaTeX proporcionada
  reference_values <- matrix(c(
    0.9, 0.13745, 0.15384, 0.18691,
    0.95, 0.09294, 0.12001, 0.15821,
    1.0, 0.05258, 0.09075, 0.13253,
    1.05, 0.02268, 0.06640, 0.10987,
    1.1, 0.00687, 0.04696, 0.09016
  ), nrow = 5, byrow = TRUE)
  
  reference_results <- as.data.frame(reference_values)
  colnames(reference_results) <- c("K", paste("σ =", sigmas))
  
  # Calcular precios con nuestra implementación
  for (i in 1:length(strikes)) {
    for (j in 1:length(sigmas)) {
      calculated_results[i, j+1] <- AWSR(spot, strikes[i], expiry, sigmas[j], rf)
    }
  }
  
  calculated_df <- as.data.frame(calculated_results)
  
  # Calcular diferencias absolutas y porcentuales
  diff_abs <- calculated_results
  diff_pct <- calculated_results
  
  for (i in 1:length(strikes)) {
    for (j in 1:length(sigmas)) {
      diff_abs[i, j+1] <- calculated_results[i, j+1] - reference_values[i, j+1]
      diff_pct[i, j+1] <- (calculated_results[i, j+1] / reference_values[i, j+1] - 1) * 100
    }
  }
  
  diff_abs_df <- as.data.frame(diff_abs)
  diff_pct_df <- as.data.frame(diff_pct)
  
  # Nombrar adecuadamente las tablas de diferencias
  colnames(diff_abs_df) <- c("K", paste("σ =", sigmas))
  colnames(diff_pct_df) <- c("K", paste("σ =", sigmas))
  
  # Crear tabla final con todos los resultados
  cat("Precios calculados con el método de Laplace:\n")
  print(calculated_df, digits = 5)
  
  cat("\nValores de referencia de los autores:\n")
  print(reference_results, digits = 5)
  
  cat("\nDiferencias absolutas (Calculado - Referencia):\n")
  print(diff_abs_df, digits = 5)
  
  cat("\nDiferencias porcentuales (%):\n")
  print(diff_pct_df, digits = 3)
  
  # Crear una tabla formateada para presentación
  comparison_table <- data.frame(
    Strike = strikes,
    Calc_Sigma_01 = calculated_results[, 2],
    Ref_Sigma_01 = reference_values[, 2],
    Diff_Sigma_01 = diff_pct[, 2],
    Calc_Sigma_03 = calculated_results[, 3],
    Ref_Sigma_03 = reference_values[, 3],
    Diff_Sigma_03 = diff_pct[, 3],
    Calc_Sigma_05 = calculated_results[, 4],
    Ref_Sigma_05 = reference_values[, 4],
    Diff_Sigma_05 = diff_pct[, 4]
  )
  
  # Crear gráfico de comparación
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    library(ggplot2)
    library(reshape2)
    
    # Preparar datos para graficar
    plot_data <- data.frame(
      Strike = rep(strikes, 2 * length(sigmas)),
      Sigma = rep(rep(sigmas, each = length(strikes)), 2),
      Value = c(as.vector(as.matrix(calculated_df[, -1])), 
                as.vector(as.matrix(reference_results[, -1]))),
      Type = rep(c("Calculado", "Referencia"), each = length(strikes) * length(sigmas))
    )
    
    # Crear gráfico
    p <- ggplot(plot_data, aes(x = Strike, y = Value, color = Type, shape = as.factor(Sigma))) +
      geom_point(size = 3) +
      geom_line(aes(linetype = Type)) +
      scale_color_manual(values = c("blue", "red")) +
      labs(title = "Comparación de precios de opciones asiáticas",
           subtitle = "Método de Laplace vs. Valores de referencia",
           x = "Strike (K)",
           y = "Precio de la opción",
           color = "Fuente",
           shape = "Volatilidad (σ)",
           linetype = "Fuente") +
      theme_minimal() +
      theme(legend.position = "bottom")
    
    print(p)
  } else {
    cat("\nPara generar gráficos, instale el paquete 'ggplot2':\n install.packages('ggplot2')\n")
  }
  
  return(list(
    calculated = calculated_df,
    reference = reference_results,
    diff_abs = diff_abs_df,
    diff_pct = diff_pct_df,
    comparison = comparison_table
  ))
}

# Ejecutar la comparación
results <- compare_asian_option_prices()

#### Comparacion con Euler-Maruyama ####

# Método de Euler-Maruyama para simular el proceso de raíz cuadrada
# dX(t) = rX(t)dt + σ√X(t)dW(t)

# Función para simular el proceso de raíz cuadrada usando Euler-Maruyama
simulate_square_root_process <- function(X0, r, sigma, T, n_steps, n_paths) {
  dt <- T/n_steps
  sqrt_dt <- sqrt(dt)
  
  # Matriz para almacenar todas las trayectorias
  X <- matrix(0, nrow = n_steps + 1, ncol = n_paths)
  X[1, ] <- X0  # Valor inicial
  
  # Simular cada trayectoria
  for (i in 1:n_paths) {
    for (j in 1:n_steps) {
      # Generar incremento browniano
      dW <- rnorm(1, mean = 0, sd = sqrt_dt)
      
      # Aplicar Euler-Maruyama
      # Para evitar valores negativos, usamos max(X[j, i], 0)
      X[j+1, i] <- X[j, i] + r * X[j, i] * dt + sigma * sqrt(max(X[j, i], 0)) * dW
      
      # Asegurar que el proceso no tome valores negativos
      X[j+1, i] <- max(X[j+1, i], 0)
    }
  }
  
  return(X)
}

# Función para calcular el precio de una opción asiática usando Monte Carlo
price_asian_option_mc <- function(spot, strike, T, r, sigma, n_steps = 1000, n_paths = 10000, control_variate = FALSE) {
  # Simular el proceso
  X <- simulate_square_root_process(spot, r, sigma, T, n_steps, n_paths)
  
  # Calcular la media aritmética para cada trayectoria
  # Para una integración más precisa, usamos la regla del trapecio
  avg_values <- numeric(n_paths)
  dt <- T/n_steps
  
  for (i in 1:n_paths) {
    # Aproximación de la integral usando la regla del trapecio
    # Corregimos los índices y aseguramos que todas las operaciones sean numéricas
    avg_values[i] <- (X[1, i]/2 + sum(X[2:n_steps, i]) + X[n_steps+1, i]/2) * dt / T
  }
  
  # El resto de la función sigue igual...
  # Calcular payoff de la opción asiática
  payoffs <- pmax(avg_values - strike, 0)
  
  # Aplicar descuento
  option_price <- exp(-r * T) * mean(payoffs)
  
  # Calcular error estándar Monte Carlo
  se <- sd(payoffs) / sqrt(n_paths)
  
  # Si se solicita, implementar control variate para reducir varianza
  if (control_variate) {
    # Usamos el valor esperado analítico del subyacente como variable de control
    # E[X(t)] = X0 * exp(rt)
    expected_avg <- spot * (exp(r * T) - 1) / (r * T)
    
    # Calcular media aritmética teórica
    MC_avg <- mean(avg_values)
    
    # Ajustar el precio usando la técnica de control variate
    beta <- -1  # Coeficiente óptimo simplificado
    option_price_cv <- option_price + beta * (MC_avg - expected_avg)
    
    return(list(price = option_price_cv, se = se))
  }
  
  return(list(price = option_price, se = se))
}

# Función para generar tabla comparativa
compare_pricing_methods <- function(laplace_fn, monte_carlo_fn) {
  # Parámetros
  spot <- 1
  T <- 1
  r <- 0.05
  
  # Valores de volatilidad y strikes
  sigmas <- c(0.1, 0.3, 0.5)
  strikes <- c(0.9, 0.95, 1.0, 1.05, 1.1)
  
  # Crear tabla para almacenar resultados
  results <- data.frame(
    K = strikes,
    sigma_0.1_laplace = 0,
    sigma_0.1_mc = 0,
    sigma_0.1_diff = 0,
    sigma_0.3_laplace = 0,
    sigma_0.3_mc = 0,
    sigma_0.3_diff = 0,
    sigma_0.5_laplace = 0,
    sigma_0.5_mc = 0,
    sigma_0.5_diff = 0
  )
  
  # Calcular precios con ambos métodos
  for (i in 1:length(strikes)) {
    for (j in 1:length(sigmas)) {
      # Calcular precio con transformada de Laplace
      laplace_price <- laplace_fn(spot, strikes[i], T, sigmas[j], r)
      
      # Calcular precio con Monte Carlo (con control variate)
      mc_result <- monte_carlo_fn(spot, strikes[i], T, r, sigmas[j], 
                                  control_variate = TRUE,
                                  n_steps = 100, 
                                  n_paths = 5000)
      mc_price <- mc_result$price
      
      # Calcular diferencia relativa
      diff_pct <- 100 * (mc_price - laplace_price) / laplace_price
      
      # Guardar resultados
      col_idx <- 3 * (j - 1) + 2
      results[i, col_idx] <- laplace_price
      results[i, col_idx + 1] <- mc_price
      results[i, col_idx + 2] <- diff_pct
    }
  }
  
  return(results)
}

# Ejemplo de uso
# Asumiendo que ya tenemos la función AWSR de la solución anterior
compare_results <- compare_pricing_methods(AWSR, price_asian_option_mc)

# Formatear resultados para mejor visualización
format_comparison_table <- function(results) {
  formatted <- results
  
  # Columnas para cada sigma
  for (j in 1:3) {
    base_col <- 3 * (j - 1) + 2
    
    # Formatear valores numéricos
    formatted[, base_col] <- sprintf("%.5f", results[, base_col])      # Laplace
    formatted[, base_col + 1] <- sprintf("%.5f", results[, base_col + 1])  # Monte Carlo
    formatted[, base_col + 2] <- sprintf("%.2f%%", results[, base_col + 2]) # Diferencia
  }
  
  # Renombrar columnas para mejor visualización
  colnames(formatted) <- c(
    "K", 
    "σ=0.1 (Laplace)", "σ=0.1 (MC)", "σ=0.1 (Diff)",
    "σ=0.3 (Laplace)", "σ=0.3 (MC)", "σ=0.3 (Diff)",
    "σ=0.5 (Laplace)", "σ=0.5 (MC)", "σ=0.5 (Diff)"
  )
  
  return(formatted)
}

# Generar tabla formateada
formatted_comparison <- format_comparison_table(compare_results)
print(formatted_comparison)

# Función para analizar la convergencia de Monte Carlo
analyze_mc_convergence <- function(spot, strike, T, r, sigma, max_paths = 1000, steps = 10) {
  path_counts <- seq(100, max_paths, length.out = steps)
  results <- data.frame(
    paths = path_counts,
    price = numeric(steps),
    se = numeric(steps)
  )
  
  # Precio de referencia (Laplace) g
  reference_price <- AWSR(spot, strike, T, sigma, r)
  
  for (i in 1:steps) {
    mc_result <- price_asian_option_mc(spot, strike, T, r, sigma, 
                                       n_steps = 500, 
                                       n_paths = path_counts[i],
                                       control_variate = TRUE)
    
    results$price[i] <- mc_result$price
    results$se[i] <- mc_result$se
  }
  
  # Añadir error relativo
  results$rel_error <- abs(results$price - reference_price) / reference_price * 100
  
  return(results)
}

# Ejemplo: analizar convergencia para K=1.0 y σ=0.3
convergence_data <- analyze_mc_convergence(1, 1.0, 1, 0.05, 0.3)

#### Comparación de los tres métodos ####
# Función para comparar los tres métodos: Valores de referencia, Laplace y Monte Carlo
compare_three_methods <- function() {
  # Parámetros
  spot <- 1
  expiry <- 1
  rf <- 0.05
  
  # Valores de volatilidad (sigma) y strikes a evaluar
  sigmas <- c(0.1, 0.3, 0.5)
  strikes <- c(0.9, 0.95, 1.0, 1.05, 1.1)
  
  # Valores de referencia de la tabla LaTeX proporcionada
  reference_values <- matrix(c(
    0.9, 0.13745, 0.15384, 0.18691,
    0.95, 0.09294, 0.12001, 0.15821,
    1.0, 0.05258, 0.09075, 0.13253,
    1.05, 0.02268, 0.06640, 0.10987,
    1.1, 0.00687, 0.04696, 0.09016
  ), nrow = 5, byrow = TRUE)
  
  reference_df <- as.data.frame(reference_values)
  colnames(reference_df) <- c("K", paste("σ =", sigmas))
  
  # Crear dataframe para resultados
  results <- data.frame(
    Strike = strikes,
    matrix(0, nrow = length(strikes), ncol = 9)  # 3 métodos x 3 volatilidades
  )
  
  colnames(results) <- c(
    "Strike",
    "Ref_σ=0.1", "Laplace_σ=0.1", "MC_σ=0.1",
    "Ref_σ=0.3", "Laplace_σ=0.3", "MC_σ=0.3",
    "Ref_σ=0.5", "Laplace_σ=0.5", "MC_σ=0.5"
  )
  
  # Calcular precios con todos los métodos
  for (i in 1:length(strikes)) {
    # Valores de referencia
    for (j in 1:length(sigmas)) {
      # Índice de columna para cada volatilidad en el dataframe de resultados
      col_idx_ref <- 3*j - 1   # 2, 5, 8
      col_idx_lap <- 3*j       # 3, 6, 9
      col_idx_mc <- 3*j + 1    # 4, 7, 10
      
      # Valor de referencia
      results[i, col_idx_ref] <- reference_values[i, j+1]
      
      # Calcular precio con transformada de Laplace
      results[i, col_idx_lap] <- AWSR(spot, strikes[i], expiry, sigmas[j], rf)
      
      # Calcular precio con Monte Carlo
      mc_result <- price_asian_option_mc(
        spot, strikes[i], expiry, rf, sigmas[j], 
        n_steps = 100, 
        n_paths = 5000, 
        control_variate = TRUE
      )
      results[i, col_idx_mc] <- mc_result$price
    }
  }
  
  # Calcular diferencias porcentuales
  diff_results <- data.frame(
    Strike = strikes,
    matrix(0, nrow = length(strikes), ncol = 6)  # 2 comparaciones x 3 volatilidades
  )
  
  colnames(diff_results) <- c(
    "Strike",
    "Laplace_vs_Ref_σ=0.1", "MC_vs_Ref_σ=0.1",
    "Laplace_vs_Ref_σ=0.3", "MC_vs_Ref_σ=0.3",
    "Laplace_vs_Ref_σ=0.5", "MC_vs_Ref_σ=0.5"
  )
  
  for (i in 1:length(strikes)) {
    for (j in 1:length(sigmas)) {
      col_idx_ref <- 3*j - 1   # 2, 5, 8 en results
      col_idx_lap <- 3*j       # 3, 6, 9 en results
      col_idx_mc <- 3*j + 1    # 4, 7, 10 en results
      
      diff_col_lap <- 2*j      # 2, 4, 6 en diff_results
      diff_col_mc <- 2*j + 1   # 3, 5, 7 en diff_results
      
      # Diferencia porcentual: (método - referencia) / referencia * 100
      ref_value <- results[i, col_idx_ref]
      lap_value <- results[i, col_idx_lap]
      mc_value <- results[i, col_idx_mc]
      
      diff_results[i, diff_col_lap] <- (lap_value - ref_value) / ref_value * 100
      diff_results[i, diff_col_mc] <- (mc_value - ref_value) / ref_value * 100
    }
  }
  
  # Formatear resultados para mejor visualización
  formatted_results <- results
  formatted_diff <- diff_results
  
  # Formatear valores de precios con 5 decimales
  for (j in 2:ncol(formatted_results)) {
    formatted_results[, j] <- sprintf("%.5f", results[, j])
  }
  
  # Formatear valores de diferencias con 2 decimales y símbolo de porcentaje
  for (j in 2:ncol(formatted_diff)) {
    formatted_diff[, j] <- sprintf("%.2f%%", diff_results[, j])
  }
  
  # Crear gráficos comparativos si está disponible ggplot2
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    library(ggplot2)
    library(reshape2)
    
    # Preparar datos para graficar precios
    plot_data <- data.frame(
      Strike = rep(strikes, 9),
      Method = rep(c(rep("Referencia", length(strikes)), 
                     rep("Laplace", length(strikes)), 
                     rep("Monte Carlo", length(strikes))), 3),
      Sigma = c(rep("0.1", 3*length(strikes)), 
                rep("0.3", 3*length(strikes)), 
                rep("0.5", 3*length(strikes))),
      Price = c(
        results[, 2], results[, 3], results[, 4],  # σ=0.1
        results[, 5], results[, 6], results[, 7],  # σ=0.3
        results[, 8], results[, 9], results[, 10]  # σ=0.5
      )
    )
    
    # Gráfico de precios
    p1 <- ggplot(plot_data, aes(x = Strike, y = Price, color = Method, shape = Method)) +
      geom_point(size = 3) +
      geom_line() +
      facet_wrap(~Sigma, scales = "free_y", ncol = 3, 
                 labeller = labeller(Sigma = function(x) paste("σ =", x))) +
      scale_color_manual(values = c("red", "blue", "green4")) +
      labs(title = "Comparación de precios de opciones asiáticas",
           subtitle = "Todos los métodos por nivel de volatilidad",
           x = "Strike (K)",
           y = "Precio de la opción") +
      theme_minimal()
    
    # Preparar datos para graficar diferencias porcentuales
    diff_plot_data <- data.frame(
      Strike = rep(strikes, 6),
      Method = rep(c(rep("Laplace vs Ref", length(strikes)), 
                     rep("Monte Carlo vs Ref", length(strikes))), 3),
      Sigma = c(rep("0.1", 2*length(strikes)), 
                rep("0.3", 2*length(strikes)), 
                rep("0.5", 2*length(strikes))),
      Diff = c(
        diff_results[, 2], diff_results[, 3],  # σ=0.1
        diff_results[, 4], diff_results[, 5],  # σ=0.3
        diff_results[, 6], diff_results[, 7]   # σ=0.5
      )
    )
    
    # Gráfico de diferencias
    p2 <- ggplot(diff_plot_data, aes(x = Strike, y = Diff, color = Method, shape = Method)) +
      geom_point(size = 3) +
      geom_line() +
      geom_hline(yintercept = 0, linetype = "dashed", color = "darkgrey") +
      facet_wrap(~Sigma, scales = "free_y", ncol = 3,
                 labeller = labeller(Sigma = function(x) paste("σ =", x))) +
      scale_color_manual(values = c("blue", "green4")) +
      labs(title = "Diferencias porcentuales respecto a valores de referencia",
           subtitle = "Por nivel de volatilidad",
           x = "Strike (K)",
           y = "Diferencia (%)") +
      theme_minimal()
    
    # Imprimir gráficos
    print(p1)
    print(p2)
  } else {
    cat("\nPara generar gráficos, instale el paquete 'ggplot2':\n install.packages('ggplot2')\n")
  }
  
  # Mostrar tablas
  cat("\nComparación de precios de opciones asiáticas:\n")
  print(formatted_results, row.names = FALSE)
  
  cat("\nDiferencias porcentuales respecto a los valores de referencia:\n")
  print(formatted_diff, row.names = FALSE)
  
  # Calcular estadísticas de error
  error_stats <- data.frame(
    Method = c("Laplace", "Monte Carlo"),
    RMSE_01 = c(
      sqrt(mean((results[, 3] - results[, 2])^2)),
      sqrt(mean((results[, 4] - results[, 2])^2))
    ),
    RMSE_03 = c(
      sqrt(mean((results[, 6] - results[, 5])^2)),
      sqrt(mean((results[, 7] - results[, 5])^2))
    ),
    RMSE_05 = c(
      sqrt(mean((results[, 9] - results[, 8])^2)),
      sqrt(mean((results[, 10] - results[, 8])^2))
    ),
    MAE_01 = c(
      mean(abs(results[, 3] - results[, 2])),
      mean(abs(results[, 4] - results[, 2]))
    ),
    MAE_03 = c(
      mean(abs(results[, 6] - results[, 5])),
      mean(abs(results[, 7] - results[, 5]))
    ),
    MAE_05 = c(
      mean(abs(results[, 9] - results[, 8])),
      mean(abs(results[, 10] - results[, 8]))
    ),
    MAPE_01 = c(
      mean(abs(diff_results[, 2])),
      mean(abs(diff_results[, 3]))
    ),
    MAPE_03 = c(
      mean(abs(diff_results[, 4])),
      mean(abs(diff_results[, 5]))
    ),
    MAPE_05 = c(
      mean(abs(diff_results[, 6])),
      mean(abs(diff_results[, 7]))
    )
  )
  
  # Formatear estadísticas de error
  formatted_stats <- error_stats
  for (j in 2:ncol(formatted_stats)) {
    if (j <= 7) {
      # RMSE y MAE con 5 decimales
      formatted_stats[, j] <- sprintf("%.5f", error_stats[, j])
    } else {
      # MAPE con 2 decimales y símbolo de porcentaje
      formatted_stats[, j] <- sprintf("%.2f%%", error_stats[, j])
    }
  }
  
  cat("\nEstadísticas de error por método y volatilidad:\n")
  print(formatted_stats, row.names = FALSE)
  
  # Retornar todos los resultados
  return(list(
    prices = results,
    differences = diff_results,
    error_stats = error_stats,
    formatted_prices = formatted_results,
    formatted_differences = formatted_diff,
    formatted_stats = formatted_stats
  ))
}

# Ejemplo de uso
comparison_results <- compare_three_methods()



#### Comparación del desempeño computacional de ambos métodos ####
# Función para medir el rendimiento de ambos métodos
compare_performance <- function() {
  # Parámetros comunes para ambos métodos
  spot <- 1
  strike <- 1.0
  T <- 1
  r <- 0.05
  sigma <- 0.3
  
  # Parámetros para escalabilidad (probaremos diferentes configuraciones)
  mc_paths_list <- c(100, 500, 1000, 5000)
  mc_steps_list <- c(100, 500, 1000)
  
  laplace_terms_list <- c(10, 20, 50)
  laplace_extraterms_list <- c(5, 10, 15)
  
  # Tabla para resultados
  results <- data.frame(
    Method = character(),
    Config = character(),
    Time_ms = numeric(),
    Memory_MB = numeric(),
    Price = numeric()
  )
  
  # Benchmark para Monte Carlo
  for (paths in mc_paths_list) {
    for (steps in mc_steps_list) {
      # Limpiar memoria antes de cada prueba
      gc()
      
      # Inicio de medición
      config <- paste0("paths=", paths, ", steps=", steps)
      start_time <- Sys.time()
      mem_start <- gc(reset = TRUE)
      
      # Ejecutar Monte Carlo
      result <- price_asian_option_mc(
        spot = spot, 
        strike = strike, 
        T = T, 
        r = r, 
        sigma = sigma, 
        n_steps = steps, 
        n_paths = paths,
        control_variate = TRUE
      )
      
      # Fin de medición
      end_time <- Sys.time()
      mem_end <- gc(reset = FALSE)
      
      # Calcular métricas
      time_ms <- as.numeric(difftime(end_time, start_time, units = "secs")) * 1000
      mem_used <- sum(mem_end[,6] - mem_start[,6]) # Medir MB
      
      # Guardar resultados
      results <- rbind(results, data.frame(
        Method = "Monte Carlo",
        Config = config,
        Time_ms = time_ms,
        Memory_MB = mem_used,
        Price = result$price
      ))
    }
  }
  
  # Benchmark para Laplace
  for (terms in laplace_terms_list) {
    for (extraterms in laplace_extraterms_list) {
      # Limpiar memoria antes de cada prueba
      gc()
      
      # Inicio de medición
      config <- paste0("terms=", terms, ", extraterms=", extraterms)
      start_time <- Sys.time()
      mem_start <- gc(reset = TRUE)
      
      # Ejecutar Laplace
      result <- AWSR(
        spot = spot, 
        strike = strike, 
        expiry = T, 
        sg = sigma, 
        rf = r,
        aa = 18.4,
        terms = terms,
        extraterms = extraterms
      )
      
      # Fin de medición
      end_time <- Sys.time()
      mem_end <- gc(reset = FALSE)
      
      # Calcular métricas
      time_ms <- as.numeric(difftime(end_time, start_time, units = "secs")) * 1000
      mem_used <- sum(mem_end[,6] - mem_start[,6]) # / 1024^2  # Convertir a MB
      
      # Guardar resultados
      results <- rbind(results, data.frame(
        Method = "Laplace",
        Config = config,
        Time_ms = time_ms,
        Memory_MB = mem_used,
        Price = result
      ))
    }
  }
  
  return(results)
}

# Función para visualizar los resultados de rendimiento
plot_performance <- function(perf_results) {
  # Convertir Method a factor para ordenar en gráficos
  perf_results$Method <- factor(perf_results$Method, 
                                levels = c("Laplace", "Monte Carlo"))
  
  # Crear resumen por método
  method_summary <- aggregate(
    cbind(Time_ms, Memory_MB) ~ Method, 
    data = perf_results, 
    FUN = function(x) c(mean = mean(x), 
                        median = median(x), 
                        min = min(x), 
                        max = max(x))
  )
  
  # Formato para mostrar resumen
  method_summary_formatted <- data.frame(
    Method = method_summary$Method,
    Mean_Time_ms = method_summary$Time_ms[,"mean"],
    Median_Time_ms = method_summary$Time_ms[,"median"],
    Min_Time_ms = method_summary$Time_ms[,"min"],
    Max_Time_ms = method_summary$Time_ms[,"max"],
    Mean_Memory_MB = method_summary$Memory_MB[,"mean"],
    Median_Memory_MB = method_summary$Memory_MB[,"median"]
  )
  
  # Si está disponible el paquete ggplot2, crear visualizaciones
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    library(ggplot2)
    
    # Gráfico de tiempo de ejecución
    p1 <- ggplot(perf_results, aes(x = Method, y = Time_ms, fill = Method)) +
      geom_boxplot() +
      scale_y_log10() +
      labs(title = "Comparación de Tiempo de Ejecución",
           x = "Método",
           y = "Tiempo (ms, escala log)") +
      theme_minimal()
    
    # Gráfico de uso de memoria
    p2 <- ggplot(perf_results, aes(x = Method, y = Memory_MB, fill = Method)) +
      geom_boxplot() +
      labs(title = "Comparación de Uso de Memoria",
           x = "Método",
           y = "Memoria (MB)") +
      theme_minimal()
    
    # Imprimir gráficos
    print(p1)
    print(p2)
  }
  
  # Devolver el resumen
  return(method_summary_formatted)
}

# Función para medir el tiempo de convergencia
compare_convergence_speed <- function() {
  # Parámetros fijos
  spot <- 1
  strike <- 1.0
  T <- 1
  r <- 0.05
  sigma <- 0.3
  
  # Valor de referencia (considerar Laplace con términos altos como "exacto")
  reference_price <- AWSR(spot, strike, T, sigma, r, aa = 18.4, terms = 50, extraterms = 20)
  
  # Error objetivo
  target_error <- 0.001  # 0.1%
  
  # Para Monte Carlo - aumentar gradualmente número de simulaciones
  mc_paths_seq <- seq(100, 1000, by = 100)
  mc_results <- data.frame(
    paths = integer(),
    time_ms = numeric(),
    price = numeric(),
    rel_error = numeric()
  )
  
  for (paths in mc_paths_seq) {
    gc()
    start_time <- Sys.time()
    
    result <- price_asian_option_mc(
      spot = spot, 
      strike = strike, 
      T = T, 
      r = r, 
      sigma = sigma, 
      n_steps = 500,  # fijo para este test
      n_paths = paths,
      control_variate = TRUE
    )
    
    end_time <- Sys.time()
    time_ms <- as.numeric(difftime(end_time, start_time, units = "secs")) * 1000
    
    rel_error <- abs(result$price - reference_price) / reference_price
    
    mc_results <- rbind(mc_results, data.frame(
      paths = paths,
      time_ms = time_ms,
      price = result$price,
      rel_error = rel_error
    ))
    
    # Terminar si alcanzamos el error objetivo
    if (rel_error <= target_error) {
      break
    }
  }
  
  # Para Laplace - aumentar gradualmente número de términos
  lap_terms_seq <- seq(5, 50, by = 5)
  lap_results <- data.frame(
    terms = integer(),
    time_ms = numeric(),
    price = numeric(),
    rel_error = numeric()
  )
  
  for (terms in lap_terms_seq) {
    gc()
    start_time <- Sys.time()
    
    result <- AWSR(
      spot = spot, 
      strike = strike, 
      expiry = T, 
      sg = sigma, 
      rf = r,
      aa = 18.4,
      terms = terms,
      extraterms = round(terms/2)  # proporcional a terms
    )
    
    end_time <- Sys.time()
    time_ms <- as.numeric(difftime(end_time, start_time, units = "secs")) * 1000
    
    rel_error <- abs(result - reference_price) / reference_price
    
    lap_results <- rbind(lap_results, data.frame(
      terms = terms,
      time_ms = time_ms,
      price = result,
      rel_error = rel_error
    ))
    
    # Terminar si alcanzamos el error objetivo
    if (rel_error <= target_error) {
      break
    }
  }
  
  # Resumen
  mc_conv_time <- mc_results[mc_results$rel_error <= target_error, "time_ms"][1]
  if (is.na(mc_conv_time)) mc_conv_time <- max(mc_results$time_ms)
  
  lap_conv_time <- lap_results[lap_results$rel_error <= target_error, "time_ms"][1]
  if (is.na(lap_conv_time)) lap_conv_time <- max(lap_results$time_ms)
  
  summary <- data.frame(
    Method = c("Monte Carlo", "Laplace"),
    Convergence_Time_ms = c(mc_conv_time, lap_conv_time),
    Final_Error = c(min(mc_results$rel_error), min(lap_results$rel_error))
  )
  
  return(list(
    summary = summary,
    mc_results = mc_results,
    lap_results = lap_results
  ))
}

# Ejecutar pruebas de rendimiento
perf_results <- compare_performance()
print(perf_results)

# Visualizar resultados
summary <- plot_performance(perf_results)
print(summary)

# Análisis de convergencia
conv_results <- compare_convergence_speed()
print(conv_results$summary)



