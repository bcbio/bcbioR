#' list of colorblind-friendly colors, developed by James and client
cb_friendly_colors <- c(
  `blue`            = "#2759F6",
  `light_orange`    = "#FFD37D",
  `olive_green`     = "olivedrab3",
  `purple`          = "#9176C8",
  `pink`            = "#E93380",
  `sky_blue`        = "#4FAEEB",
  `blue_grey`       = "#92A6BC",
  `forest_green`    = "#3C877B",
  `yellow`          = "yellow",
  `dark_purple`     = "#402999",
  `dark_orange`     = "#D5392C",
  `army_green`      = "#C3C380",
  `black`           = "black",
  `dark_grey`       = "darkgrey",
  `light_blue`      = "lightblue",
  `brown`           = "#661100",
  `white`           = "white"
)

#' list cb friendly colors
#' @export
list_cb_friendly_cols <- function(){
  return(cb_friendly_colors)
}

#' fetch color from list by name
#'
#' @param ... list of color names
#'
#' @export
cb_friendly_cols <- function(...) {
  cols <- c(...)

  if (is.null(cols))
    return (cb_friendly_colors)

  cb_friendly_colors[cols]
}

#' define main colorblind-friendly palette as well as sub-palettes
cb_friendly_palettes <- list(
  `main`  = cb_friendly_cols("blue", "purple", "sky_blue",
                               "blue_grey", "forest_green", "pink", "olive_green",
                               "yellow", "dark_purple", "dark_orange",
                               "army_green", "black", "dark_grey", "light_blue",
                               "brown","light_orange"),
  `cool`  = cb_friendly_cols("blue", "dark_purple", "purple", "sky_blue"),
  `hot`   = cb_friendly_cols("yellow", "light_orange", "dark_orange"),
  `grey`  = cb_friendly_cols("black", "dark_grey", "blue_grey"),
  `heatmap` = cb_friendly_cols("blue", "white", "brown")
)

#' access cb friendly palette by name, reversing if necessary
#'
#' @param palette name of the palette to be returned
#' @param reverse boolean, reverse order of colors in palette
#' @param ... pass to ggplot
#' @export
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
#' @param ... pass to ggplot
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
#' @param ... pass to ggplot
#' @export
scale_fill_cb_friendly <- function(palette = "main", discrete = TRUE, reverse = FALSE, ...) {
  pal <- cb_friendly_pal(palette = palette, reverse = reverse)

  if (discrete) {
    discrete_scale("fill", paste0("cb_friendly_", palette), palette = pal, ...)
  } else {
    scale_fill_gradientn(colours = pal(256), ...)
  }
}
