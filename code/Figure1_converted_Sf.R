library(ggplot2)
library(latex2exp)

df <- data.frame(
  t  = c(0, 1, 2, 3, 4),
  Sf = c(1, 0.75, 0.5, 0.25, 0)
)

hseg <- data.frame(
  x    = c(0, 0, 0),
  xend = c(1, 2, 3),
  y    = c(0.75, 0.50, 0.25),
  yend = c(0.75, 0.50, 0.25)
)

vseg <- data.frame(
  x    = c(1, 2, 3),
  xend = c(1, 2, 3),
  y    = c(0, 0, 0),
  yend = c(0.75, 0.50, 0.25)
)

ggplot(df, aes(x = t, y = Sf)) +
  geom_segment(
    data = hseg,
    aes(x = x, xend = xend, y = y, yend = yend),
    inherit.aes = FALSE,
    linetype = "dotted",
    linewidth = 0.7
  ) +
  geom_segment(
    data = vseg,
    aes(x = x, xend = xend, y = y, yend = yend),
    inherit.aes = FALSE,
    linetype = "dotted",
    linewidth = 0.7
  ) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 4, stroke = 0) +
  coord_cartesian(
    xlim = c(0, 4.1),
    ylim = c(0, 1.05),
    expand = FALSE,
    clip = "off"
  ) +
  scale_x_continuous(
    breaks = 0:4,
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = c(0.25, 0.50, 0.75, 1.00),
    expand = c(0, 0)
  ) +
  labs(
    x = TeX(r'($t$)'),
    y = TeX(r'(Converted Sf  $\bar{G}_1^*(t)$)')
  ) +
  theme_classic(base_size = 12,base_family = "serif") +
  theme(
    axis.title.x = element_text(size = 17),
    axis.title.y = element_text(size = 17),
    axis.text.x  = element_text(size = 15),
    axis.text.y  = element_text(size = 15)
  )
