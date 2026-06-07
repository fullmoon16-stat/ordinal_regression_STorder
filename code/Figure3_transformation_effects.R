library(ggplot2)
library(latex2exp)
library(gridExtra)

# grid
tt  <- seq(0, 1, length.out = 1000)

# baseline
S1 <- function(t) 1 - t

# H families
H1 <- function(t, k) 1 - t^k
H2 <- function(t, k) (1 - t)^k
H3 <- function(t, k) plogis(k - qlogis(t))
H4 <- function(t, k, s) plogis(k - s * qlogis(t))

# -----------------------------
# parameter choices
# -----------------------------
k12_vals <- c(0.3, 0.5, 2, 4)
k3_vals  <- c(-2, -1, 1, 2)
s4_vals  <- c(2)
k4_fix   <- c(-2,-1,0,1,2)

# -----------------------------
# helper plotting function
# -----------------------------
make_plot <- function(df, y_label) {
  
  # baseline / transformed 분리
  df_base <- subset(df, type == "Baseline")
  df_tr   <- subset(df, type == "Transformed")
  df_tr$group <- factor(df_tr$group, levels = unique(df_tr$group))
  
  ggplot() +
    geom_line(
      data = df_base,
      aes(x = t, y = y),
      color = "black",
      linewidth = 1.2,
      linetype = "solid",
      show.legend = FALSE
    ) +
    # transformed curves: 색 + linetype 모두 group에 매핑
    geom_line(
      data = df_tr,
      aes(x = t, y = y, color = group, linetype = group),
      linewidth = 1.2
    ) + 
    annotate(
      "text",
      x = 1.05, y = 0,
      label = "t",
      size = 7,family = "serif"
    ) +
    annotate(
      "text",
      x = 0.03, y = 1.07,
      label = y_label,
      parse = TRUE,
      size = 7, family = "serif"
    ) +
    coord_cartesian(xlim = c(0, 1.02), ylim = c(0, 1.02), expand = FALSE, clip = "off" ) +
    scale_x_continuous(
      breaks = seq(0, 1, by = 0.2),
      labels = function(x) sprintf("%.1f", x),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      breaks = c(0.2,0.4,0.6,0.8,1),
      labels = function(y) sprintf("%.1f", y),
      expand = c(0, 0)
    ) +
    theme_classic(base_size = 12,base_family = "serif") +
    theme(
      legend.position = c(0.88, 0.86),
      legend.title = element_blank(),
      axis.text = element_text(size = 14),
      axis.title = element_blank(),
      legend.text = element_text(size = 16),
      
      legend.key.width  = grid::unit(1.1, "cm"),  
      legend.key.height = grid::unit(1, "cm"),
      plot.margin = margin(30, 25, 10, 10)
    )
}

# -----------------------------
# Figure 3(a): S1(rho_11)
# -----------------------------
df1 <- data.frame(
  t = tt,
  y = S1(tt),
  group = "1-t",
  type = "Baseline"
)

for (k in k12_vals) {
  df1 <- rbind(
    df1,
    data.frame(
      t = tt,
      y = H1(tt, k),
      group = paste0("k = ", k),
      type = "Transformed"
    )
  )
}

p1 <- make_plot(df1, y_label = "S[1](rho[1*','*1](t))")


# -----------------------------
# Figure 3(b): S1(rho_12)
# -----------------------------
df2 <- data.frame(
  t = tt,
  y = S1(tt),
  group = "1-t",
  type = "Baseline"
)

for (k in k12_vals) {
  df2 <- rbind(
    df2,
    data.frame(
      t = tt,
      y = H2(tt, k),
      group = paste0("k = ", k),
      type = "Transformed"
    )
  )
}

p2 <- make_plot(df2, y_label = "S[1](rho[1*','*2](t))")

# -----------------------------
# Figure 3(c): S1(rho_13)
# -----------------------------
df3 <- data.frame(
  t = tt,
  y = S1(tt),
  group = "1-t",
  type = "Baseline"
)

for (k in k3_vals) {
  df3 <- rbind(
    df3,
    data.frame(
      t = tt,
      y = H3(tt, k),
      group = paste0("k = ", k),
      type = "Transformed"
    )
  )
}

p3 <- make_plot(df3, y_label = "S[1](rho[1*','*3](t))")

# -----------------------------
# Figure 3(d): S1(rho_14)
# -----------------------------
df4 <- data.frame(
  t = tt,
  y = S1(tt),
  group = "1-t",
  type = "Baseline"
)

for (s in s4_vals) {
  for (k in k4_fix) {
    df4 <- rbind(
      df4,
      data.frame(
        t = tt,
        y = H4(tt, k = k, s = s),
        group = paste0("k = ", k,", s = ", s),
        type = "Transformed"
      )
    )
  }
  
}

p4 <- make_plot(df4, y_label = "S[1](rho[1*','*4](t))")


##################################


# baseline
S1 <- function(t) t

# H families
H1 <- function(t, k) t^k
H2 <- function(t, k) 1-(1 - t)^k
H3 <- function(t, k) plogis(-k + qlogis(t))
H4 <- function(t, k, s) plogis(-k + s * qlogis(t))

# -----------------------------
# parameter choices
# -----------------------------
k12_vals <- c(0.3, 0.5, 2, 4)
k3_vals  <- c(-2, -1, 1, 2)
s4_vals  <- c(2)
k4_fix   <- c(-2,-1,0,1,2)

# -----------------------------
# helper plotting function
# -----------------------------
make_plot <- function(df, y_label) {
  
  # baseline / transformed 분리
  df_base <- subset(df, type == "Baseline")
  df_tr   <- subset(df, type == "Transformed")
  df_tr$group <- factor(df_tr$group, levels = unique(df_tr$group))
  
  ggplot() +
    geom_line(
      data = df_base,
      aes(x = t, y = y),
      color = "black",
      linewidth = 1.2,
      linetype = "solid",
      show.legend = FALSE
    ) +
    # transformed curves: 색 + linetype 모두 group에 매핑
    geom_line(
      data = df_tr,
      aes(x = t, y = y, color = group, linetype = group),
      linewidth = 1.2
    ) + 
    annotate(
      "text",
      x = 1.05, y = 0,
      label = "t",
      size = 7,family = "serif"
    ) +
    annotate(
      "text",
      x = 0.01, y = 1.07,
      label = y_label,
      parse = TRUE,
      size = 7,family = "serif"
    ) +
    coord_cartesian(xlim = c(0, 1.02), ylim = c(0, 1.02), expand = FALSE, clip = "off" ) +
    scale_x_continuous(
      breaks = seq(0, 1, by = 0.2),
      labels = function(x) sprintf("%.1f", x),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      breaks = c(0.2,0.4,0.6,0.8,1),
      labels = function(y) sprintf("%.1f", y),
      expand = c(0, 0)
    ) +
    theme_classic(base_size = 12) + 
    theme(
      legend.position = "none",
      axis.text = element_text(size = 14),
      axis.title = element_blank(),
      plot.margin = margin(30, 25, 10, 10)
    ) 
}

# -----------------------------
# Figure 3(a): S1(rho_11)
# -----------------------------
df1 <- data.frame(
  t = tt,
  y = S1(tt),
  group = "1-t",
  type = "Baseline"
)

for (k in k12_vals) {
  df1 <- rbind(
    df1,
    data.frame(
      t = tt,
      y = H1(tt, k),
      group = paste0("k = ", k),
      type = "Transformed"
    )
  )
}

p11 <- make_plot(df1, y_label = "rho[1*','*1](t)")

# -----------------------------
# Figure 3(b): S1(rho_12)
# -----------------------------
df2 <- data.frame(
  t = tt,
  y = S1(tt),
  group = "1-t",
  type = "Baseline"
)

for (k in k12_vals) {
  df2 <- rbind(
    df2,
    data.frame(
      t = tt,
      y = H2(tt, k),
      group = paste0("k = ", k),
      type = "Transformed"
    )
  )
}

p22 <- make_plot(df2, y_label = "rho[1*','*2](t)")

# -----------------------------
# Figure 3(c): S1(rho_13)
# -----------------------------
df3 <- data.frame(
  t = tt,
  y = S1(tt),
  group = "1-t",
  type = "Baseline"
)

for (k in k3_vals) {
  df3 <- rbind(
    df3,
    data.frame(
      t = tt,
      y = H3(tt, k),
      group = paste0("k = ", k),
      type = "Transformed"
    )
  )
}

p33 <- make_plot(df3, y_label = "rho[1*','*3](t)")

# -----------------------------
# Figure 3(d): S1(rho_14)
# -----------------------------
df4 <- data.frame(
  t = tt,
  y = S1(tt),
  group = "1-t",
  type = "Baseline"
)

for (s in s4_vals) {
  for (k in k4_fix) {
    df4 <- rbind(
      df4,
      data.frame(
        t = tt,
        y = H4(tt, k = k, s = s),
        group = paste0("k = ", k,", s = ", s),
        type = "Transformed"
      )
    )
  }
  
}

p44 <- make_plot(df4, y_label = "rho[1*','*4](t)")






########### 한꺼번에 그리기 #########
grid.arrange(p11, p1, ncol = 2)
grid.arrange(p22, p2, ncol = 2)
grid.arrange(p33, p3, ncol = 2)
grid.arrange(p44, p4, ncol = 2)

