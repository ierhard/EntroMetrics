#' Convert entropy units (based on )
#'
#' @param entropy Entropy value (numeric)
#' @param current_unit Current unit of entropy (`"log2"` = bits, `"ln"` = nats)
#' @param new_unit Desired unit of entropy (`"log2"` for bits, `"ln"` for nats, and `"normalize"` for normalized entropy/evenness)
#' @param alphabet_size Size of the alphabet of the distribution whose entropy is in question (only needed if `new_unit == "normalize"`)
#'
#' @return Entropy in desired units
#'
#' @examples
#' convert_entropy_units(entropy = 1.1, current_unit = "log2", new_unit = "ln")
convert_entropy_units <- function(entropy,
                                  current_unit = c("log2", "ln"),
                                  new_unit = c("log2", "ln", "normalize"),
                                  alphabet_size) {

  current_unit <- match.arg(current_unit)
  new_unit <- match.arg(new_unit)

  current_unit_base <- switch(
    current_unit,
    "log2" = 2,
    "ln" = exp(1)
  )

  new_unit_base <- switch(
    new_unit,
    "log2" = 2,
    "ln" = exp(1),
    "normalize" = "normalize"
  )

  if (new_unit_base != "normalize") {
    entropy_new <- entropy * log(x = current_unit_base, base = new_unit_base)
  } else {
    entropy_new <- entropy / log(x = alphabet_size, base = current_unit_base)
  }

  return(entropy_new)
}

