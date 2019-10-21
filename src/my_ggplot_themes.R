my_theme = function(grid = NULL) {
  
  if (is.null(grid)) {
    p = background_grid(major = "xy", minor = "none", size.major = 0.4) 
  } else if (grid == "y") {
    p = background_grid(major = "y", minor = "none", size.major = 0.4) 
  } else if (grid == "x") {
    p = background_grid(major = "x", minor = "none", size.major = 0.4) 
  } else if (grid == "no") {
    p = theme()
  }
  p = p + 
    theme(title = element_text(size=14),
          axis.text = element_text(size=14),
          legend.text = element_text(size=14),
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
