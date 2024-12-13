% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/update_multiple_methods_with_shared_args.R
\name{update_multiple_methods_with_shared_args}
\alias{update_multiple_methods_with_shared_args}
\title{Update multiple methods with shared arguments}
\usage{
update_multiple_methods_with_shared_args(
  multiple_methods,
  shared_args_list,
  type = c("pt", "ci")
)
}
\arguments{
\item{multiple_methods}{List of lists, each containing method as a character and the method arguments as a list (although names can differ, e.g. \code{"pt_method"} and \code{"pt_method_args"})}

\item{shared_args_list}{List of shared arguments to be assigned to each element in multiple_methods (if not locally defined)}
}
\value{
Updated multiple_methods list
}
\description{
Update multiple methods with shared arguments
}
\examples{
multiple_methods <- list(m1 = list(pt_method = "ML", pt_method_args = NULL),
                         m2 = list(pt_method = "Dirichlet", pt_method_args = NULL))
print(multiple_methods)
shared_args_list <- list(Dir_alpha = 10)
multiple_methods <- update_multiple_methods_with_shared_args(multiple_methods, shared_args_list)
print(multiple_methods)

multiple_methods <- list(m1 = list(pt_method = "ML", pt_method_args = NULL),
                         m2 = list(pt_method = "Dirichlet", pt_method_args = list(Dir_alpha = 5)))
print(multiple_methods)
shared_args_list <- list(Dir_alpha = 10)
multiple_methods <- update_multiple_methods_with_shared_args(multiple_methods, shared_args_list)
print(multiple_methods)
}
