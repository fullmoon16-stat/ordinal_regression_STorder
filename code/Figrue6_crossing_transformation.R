library(ggplot2)

## -----------------------------
## 1. fixed parameters
## -----------------------------
phi    <- 0.8
psi    <- 3.5
m      <- 6
j_star <- 3

## -----------------------------
## 2. rho functions
## -----------------------------
rho1 <- function(j) {
  j + phi * sin(pi / psi * j)
}

rho2 <- function(j) {
  j - phi * sin(pi / (m - psi) * (j - psi))
}

## -----------------------------
## 3. data for discontinuous curve
##    first piece: 0 <= j <= 3
##    second piece: 3 < j <= 6
## -----------------------------
seg1 <- data.frame(
  j   = seq(0, j_star, length.out = 400)
)
seg1$rho <- rho1(seg1$j)

seg2 <- data.frame(
  j   = seq(j_star + 1e-4, m, length.out = 400)
)
seg2$rho <- rho2(seg2$j)

## integer points
pts <- data.frame(
  j   = 0:m,
  rho = c(rho1(0:j_star), rho2((j_star + 1):m))
)

## -----------------------------
## 4. plot
## -----------------------------
p <- ggplot() +
  ## identity line: rho(j) = j
  geom_segment(
    aes(x = 0, y = 0, xend = m, yend = m),
    linewidth = 1.1,
    color = "black",
    linetype = "solid"
  ) +
  
  ## piecewise rho(j), drawn as two separate segments
  geom_line(
    data = seg1,
    aes(x = j, y = rho),
    linewidth = 1.1,
    color = "black",
    linetype = "twodash"
  ) +
  geom_line(
    data = seg2,
    aes(x = j, y = rho),
    linewidth = 1.1,
    color = "black",
    linetype = "twodash"
  ) +
  
  ## integer points
  geom_point(
    data = pts,
    aes(x = j, y = rho),
    size = 4.5,
    shape = 16,
    color = "black"
  ) +
  coord_cartesian(
    expand = FALSE,
    clip = "off"
  ) +
  scale_x_continuous(
    breaks = 0:m,
    limits = c(0, m+ 0.25),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = 1:m,
    limits = c(0, m + 0.25),
    expand = c(0, 0)
  ) +
  
  labs(
    x = expression(j),
    y = expression(rho(j))
  ) +
  
  theme_classic(base_size = 12) +
  theme(
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 15),
    legend.position = "none"
  )

p

