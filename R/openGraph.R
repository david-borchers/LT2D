#' Function to open an operating system specific graphics device.
#'
#' @param width, height the (nominal) width and height of the canvas of the plotting window in inches
#' @param ... other arguments passed into graphics device
#' @author  John K. Kruschke
#' @note function downloaded from \link{http://doingbayesiandataanalysis.blogspot.com.au/2013/01/uniform-r-code-for-opening-saving.html}
#' @export
openGraph = function( width=7 , height=7 , ... ) {
  if ( .Platform$OS.type != "windows" ) { # Mac OS, Linux
    X11( width=width , height=height , type="cairo" , ... ) 
  } else { # Windows OS
    windows( width=width , height=height , ... )
  }
}
