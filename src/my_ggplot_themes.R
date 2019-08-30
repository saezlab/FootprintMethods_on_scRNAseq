my_theme = function() {
  background_grid(major = "xy", minor = "none", size.major = 0.4) +
    theme(title = element_text(size=14),
          axis.text = element_text(size=14),
          legend.text = element_text(size=14),
          legend.title = element_blank(),
          strip.background = element_rect(colour = "white", fill="white"),
          strip.text.x  = element_text(size=14)
    )
}

my_poster_theme = function() {
  background_grid(major = "xy", minor = "none", size.major = 0.4) +
    theme(title = element_text(size=24),
          axis.text = element_text(size=24),
          legend.text = element_text(size=24),
          legend.title = element_blank()
    )
}
