#' list of colorblind-friendly colors, developed by James and client
cb_friendly_colors <- c(
  `bright_purple`   = "#2759F6",
  `dark_purple`     = "#402999",
  `purple`          = "#9176C8",
  `blue`            = "#4FAEEB",
  `blue_grey`       = "#92A6BC",
  `forest_green`    = "#3C877B",
  `pink`            = "#E93380",
  `olive_green`     = "olivedrab3",
  `yellow`          = "yellow",
  `light_orange`    = "#FFD37D",
  `dark_orange`     = "#D5392C",
  `army_green`      = "#C3C380",
  `black`           = "black",
  `dark_grey`       = "darkgrey",
  `light_blue`      = "lightblue",
  `brown`           = "#661100"   
)

#' fetch color from list by name
cb_friendly_cols <- function(...) {
  cols <- c(...)
  
  if (is.null(cols))
    return (cb_friendly_colors)
  
  cb_friendly_colors[cols]
}

#' define main colorblind-friendly palette as well as sub-palettes
cb_friendly_palettes <- list(
  `main`  = cb_friendly_cols("bright_purple", "dark_purple", "purple", "blue", 
                               "blue_grey", "forest_green", "pink", "olive_green", 
                               "yellow", "light_orange", "dark_orange",
                               "army_green", "black", "dark_grey", "light_blue", 
                               "brown"),
  `cool`  = cb_friendly_cols("bright_purple", "dark_purple", "purple", "blue"),
  `hot`   = cb_friendly_cols("yellow", "light_orange", "dark_orange"),
  `grey`  = cb_friendly_cols("black", "dark_grey", "blue_grey")
)

#' access cb friendly palette by name, reversing if necessary
#' 
#' @param palette name of the palette to be returned
#' @param reverse boolean, reverse order of colors in palette
cb_friendly_pal <-function(palette = 'main', reverse = F, ...){
  pal <- cb_friendly_palettes[[palette]]
  if (reverse) pal <- rev(pal)
  colorRampPalette(pal, ...)
}

#' use cb friendly colors as color aesthetic with ggplot
#' 
#' @param palette name of the palette to be returned
#' @param discrete boolean, whether to make palette discretely divided into colors or continuous
#' @param reverse boolean, reverse order of colors in palette
#' @export
scale_color_cb_friendly <- function(palette = "main", discrete = TRUE, reverse = FALSE, ...) {
  pal <- cb_friendly_pal(palette = palette, reverse = reverse)
  
  if (discrete) {
    discrete_scale("colour", paste0("cb_friendly_", palette), palette = pal, ...)
  } else {
    scale_color_gradientn(colours = pal(256), ...)
  }
}

#' use cb friendly colors as fill aesthetic with ggplot
#' 
#' @param palette name of the palette to be returned
#' @param discrete boolean, whether to make palette discretely divided into colors or continuous
#' @param reverse boolean, reverse order of colors in palette
#' @export
scale_fill_cb_friendly <- function(palette = "main", discrete = TRUE, reverse = FALSE, ...) {
  pal <- cb_friendly_pal(palette = palette, reverse = reverse)
  
  if (discrete) {
    discrete_scale("fill", paste0("cb_friendly_", palette), palette = pal, ...)
  } else {
    scale_fill_gradientn(colours = pal(256), ...)
  }
}