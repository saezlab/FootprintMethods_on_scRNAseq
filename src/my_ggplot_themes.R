my_theme = function() {
  background_grid(major = "xy", minor = "none", size.major = 0.4) +
    theme(title = element_text(size=16),
          axis.text = element_text(size=14),
          legend.text = element_text(size=14),
          legend.title = element_blank()
    )
}