library(ggplot2)
library(latex2exp)

t <- seq(0, 1, length.out = 2001)

# baselines
S1 <- 1 - t
S2 <- (1 - t)^2
S3 <- (1 - t)^(1/2)
S4  <- plogis(-2   * qlogis(t))
S5  <- plogis(-0.5 * qlogis(t))

# helper: long format via stack()
to_long <- function(t, mat){
  st <- stack(as.data.frame(mat))  # values, ind
  data.frame(
    t  = rep(t, times = ncol(mat)),
    id = factor(st$ind, levels = colnames(mat)),
    S  = st$values
  )
}

# ---- Figure 1: S1, S2, S3 ----
df1 <- to_long(t, cbind(S1 = S1, S2 = S2, S3 = S3))

lab_vec <- c(
  S1 = expression(S[1](t)),
  S2 = expression(S[2](t)),
  S3 = expression(S[3](t))
)

p1 <- ggplot(df1, aes(x = t, y = S, color = id, linetype = id, group = id)) +
  geom_line(linewidth = 1.1) +
  coord_cartesian(xlim = c(0, 1.05), ylim = c(0, 1.05), expand = FALSE) +
  scale_x_continuous(
    limits = c(0, 1.05),
    breaks = seq(0, 1, 0.2),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(0, 1.05),
    breaks = c(0.2, 0.4, 0.6, 0.8, 1),
    expand = c(0, 0)
  ) +
  scale_color_manual(
    name   = NULL,
    values = c(S1 = "#000000", S2 = "#CC3333", S3 = "#0072B2"), # #FF6A98 #F8766D #EA8331
    breaks = c("S1", "S2", "S3"),
    labels = lab_vec
  ) +
  scale_linetype_manual(
    name   = NULL,
    values = c(S1 = "solid", S2 = "22", S3 = "42"),
    breaks = c("S1", "S2", "S3"),
    labels = lab_vec
  ) +
  labs(
    x = TeX(r"($t$)"),
    y = TeX(r"(Baseline Sf  $S(t)$)")
  ) +
  theme_classic(base_size = 12,base_family = "serif") +
  theme(
    legend.position = c(0.8, 0.85),
    legend.title = element_blank(),
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 17),
    legend.text = element_text(size = 17),
    legend.key.width  = grid::unit(1.3, "cm"),
    legend.key.height = grid::unit(1.3, "cm")
  )

p1

# ---- Figure 2: S1, S4, S5 ----
df1 <- to_long(t, cbind(S1 = S1, S4 = S4, S5 = S5))

lab_vec <- c(
  S1 = expression(S[1](t)),
  S4 = expression(S[4](t)),
  S5 = expression(S[5](t))
)

p1 <- ggplot(df1, aes(x = t, y = S, color = id, linetype = id, group = id)) +
  geom_line(linewidth = 1.1) +
  coord_cartesian(xlim = c(0, 1.05), ylim = c(0, 1.05), expand = FALSE) +
  scale_x_continuous(
    limits = c(0, 1.05),
    breaks = seq(0, 1, 0.2),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(0, 1.05),
    breaks = c(0.2, 0.4, 0.6, 0.8, 1),
    expand = c(0, 0)
  ) +
  scale_color_manual(
    name   = NULL,
    values = c(S1 = "#000000", S4 = "#7CAE00", S5 = "#B983FF"), # #E58700 #00BA38 #7CAE00
    breaks = c("S1", "S4", "S5"), 
    labels = lab_vec
  ) +
  scale_linetype_manual(
    name   = NULL,
    values = c(S1 = "solid", S4 = "44", S5 = "2262"),
    breaks = c("S1", "S4", "S5"),
    labels = lab_vec
  ) +
  labs(
    x = TeX(r"($t$)"),
    y = TeX(r"(Baseline Sf  $S(t)$)")
  ) +
  theme_classic(base_size = 12,base_family = "serif") +
  theme(
    legend.position = c(0.8, 0.85),
    legend.title = element_blank(),
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 17),
    legend.text = element_text(size = 17),
    legend.key.width  = grid::unit(1.3, "cm"),
    legend.key.height = grid::unit(1.2, "cm")
  )

p1
