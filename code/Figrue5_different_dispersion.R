library(ggplot2)

## -----------------------------
## 1. parameters
## -----------------------------
mu1 <- -0.5
sd1 <- 1.2  # smaller dispersion

mu2 <- 0.5
sd2 <- 2.8   # larger dispersion

x <- seq(-7, 8, length.out = 1000)

## -----------------------------
## 2. data for pdf and Sf
## -----------------------------
df_pdf <- data.frame(
  x = rep(x, 2),
  y = c(dnorm(x, mean = mu1, sd = sd1),
        dnorm(x, mean = mu2, sd = sd2)),
  group = factor(rep(c("1", "2"), each = length(x)))
)

df_sf <- data.frame(
  x = rep(x, 2),
  y = c(pnorm(x, mean = mu1, sd = sd1, lower.tail = FALSE),
        pnorm(x, mean = mu2, sd = sd2, lower.tail = FALSE)),
  group = factor(rep(c("1", "2"), each = length(x)))
)

## -----------------------------
## 3. pdf plot
## -----------------------------
p_pdf <- ggplot(df_pdf, aes(x = x, y = y, linetype = group)) +
  geom_line(linewidth = 1.2) +
  geom_hline(yintercept = 0, linewidth = 0.7) + 
  scale_linetype_manual(values = c("solid", "twodash")) +
  coord_cartesian(
    xlim = c(-6.5, 7.2),
    ylim = c(0, 0.35)
  ) +
  #labs(title = "Underlying Density Functions") +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "none"
  )

## -----------------------------
## 4. Sf plot
## -----------------------------
p_sf <- ggplot(df_sf, aes(x = x, y = y, linetype = group)) +
  geom_line(linewidth = 1.2) +
  geom_hline(yintercept = 0, linewidth = 0.7) + 
  scale_linetype_manual(values = c("solid", "twodash")) +
  coord_cartesian(
    xlim = c(-5.8, 6.2),
    ylim = c(0, 1.05)
  ) +
  #labs(title = "Underlying Survival Functions") +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "none"
  )

## -----------------------------
## 5. print separately
## -----------------------------
p_pdf
p_sf

