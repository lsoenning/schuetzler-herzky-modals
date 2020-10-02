
my_settings <- list()
my_settings$axis.line$col <- NA
my_settings$strip.border$col <- NA
my_settings$strip.background$col <- NA
my_settings$layout.heights$strip = 1.4
my_settings$layout.heights$axis.xlab.padding <- 0
my_settings$layout.widths$ylab.axis.padding <- 0
my_settings$clip$panel = "off"

my_settings$layout.heights$axis.xlab.padding <- 0
my_settings$layout.heights$key.axis.padding <- 0
my_settings$layout.heights$xlab <- 0
my_settings$layout.heights$bottom.padding <- 0.2
my_settings$layout.heights$top.padding <- 1
my_settings$layout.heights$key.padding <- 0
my_settings$layout.heights$axis.top <- 0
my_settings$layout.heights$main <- 0
my_settings$layout.heights$main.key.padding <- 0

my_settings$layout.widths$ylab.axis.padding <- 0 
my_settings$layout.widths$ylab <- 0
my_settings$layout.widths$axis.panel <- 10
my_settings$layout.widths$left.padding <- 1.2
my_settings$layout.widths$right.padding <- 0
my_settings$layout.widths$key.right <- 0
my_settings$layout.widths$axis.right <- 0
my_settings$layout.widths$ylab.right <- 0

my_settings$axis.text$col <- "grey30"

my_settings$axis.components$left$tck <- .6
my_settings$axis.components$left$pad1 <- 0.5
my_settings$axis.components$left$pad2 <- 2
my_settings$axis.components$bottom$tck <- .6
my_settings$axis.components$bottom$pad1 <- 0.5
my_settings$axis.components$bottom$pad2 <- 1.4

my_settings$box.rectangle$col = 1
my_settings$box.umbrella$col = 1
my_settings$box.dot$col = 1
my_settings$plot.symbol$col = 1


# function for text with edges

panel.text_halo <- 
	function(x, y=NULL, labels, col='black', bg='white', 
					 theta= seq(0, 2*pi, length.out=50), r=0.1, ... ) {
		
		xy <- xy.coords(x,y)
		xo <- r*strwidth('A')
		yo <- r*strheight('A')
		
		# draw background text with small shift in x and y in background colour
		for (i in theta) {
			panel.text(x= xy$x + cos(i)*xo, 
								 y= xy$y + sin(i)*yo, labels, col=bg, ... )
		}
		# draw actual text in exact xy position in foreground colour
		panel.text(xy$x, xy$y, labels, col=col, ... )
	}


# functions for axes

axis_L <-
	function(side, ..., line.col)
	{
		if (side %in% c("left", "bottom")) {
			col <- trellis.par.get("axis.text")$col
			axis.default(side, ..., line.col = col)
			if (side == "bottom")
				grid::grid.lines(y = 0)
			if (side == "left")
				grid::grid.lines(x = 0)
		}
	}


axis_left = function(side, ..., line.col)
{
	if (side %in% "left") {
		col <- trellis.par.get("axis.text")$col
		axis.default(side, ..., line.col = col)
		if (side == "left")
			grid::grid.lines(x = 0)
	}
}

axis_bottom = function(side, ..., line.col)
{
	if (side %in% "bottom") {
		col <- trellis.par.get("axis.text")$col
		axis.default(side, ..., line.col = col)
		if (side == "bottom")
			grid::grid.lines(y = 0)
	}
}

axis_top = function(side, ..., line.col)
{
	if (side %in% "top") {
		col <- trellis.par.get("axis.text")$col
		axis.default(side, ..., line.col = col)
		if (side == "top")
			grid::grid.lines(y = 1)
	}
}

