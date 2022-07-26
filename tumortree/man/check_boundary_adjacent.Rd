% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_spatial_info.R
\name{check_boundary_adjacent}
\alias{check_boundary_adjacent}
\title{Check if cell is adjacent boundary}
\usage{
check_boundary_adjacent(x_loc, y_loc, alive_cells)
}
\arguments{
\item{x_loc}{numeric}

\item{y_loc}{numeric}

\item{alive_cells}{data.frame}
}
\value{
integer
}
\description{
Check if cell is next to an extreme x or y location in the tumor to determine edge.
}
