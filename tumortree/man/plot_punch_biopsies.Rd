% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bulk_sample_sim.R
\name{plot_punch_biopsies}
\alias{plot_punch_biopsies}
\title{Plot punch biopsies}
\usage{
plot_punch_biopsies(
  alive_cells,
  punched_cell_indices_list,
  color_by_edge = FALSE,
  color_by_punch_id = FALSE,
  size = 2
)
}
\arguments{
\item{alive_cells}{data.frame all alive cells in simulated tumor}

\item{color_by_edge}{logical if TRUE, color by edge and center classification}

\item{color_by_punch_id}{logical if TRUE, color by punch id}

\item{size}{numeric size scale for scatter plot points}
}
\value{
ggplot
}
\description{
Visualize punch biopsies on simulated tumor
}
