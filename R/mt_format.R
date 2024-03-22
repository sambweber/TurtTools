

# -------------------------------------------------------------------------------------------------------------------------
# mt_format
# -------------------------------------------------------------------------------------------------------------------------

#' Puts data in right format for a MoniTool model run
#
#' @param data is a dataframe which must have a column called 'date' which was the date on which
#' the count was done, along with at least one of columns entitled'nests' or 'activities' (can be both).
#' May also optionally have columns called 'beach' to differentiate sites, and 'season' to name nesting
#' seasons. If 'season' is missing it will be guessed from the dates provided.

#' @param max.days is a number giving the maximum number of days over which the
#' activities counted in each survey may have occurred. So `max.days = 1` means all activities occurred on the
#' the night before survey was done (for morning counts). For `max.days > 1` counted
#' activities are assumed to have occurred over all days since the previous count, up to
#' a maximum of `date - max.days`. This requires you to make some informed assumption about
#' how many days activities remain visible for. For more complex monitoring designs it is
#' also possible to include a column in `data` called `datestart` which contains the earliest dates
#' when counts in each survey could have occurred. For example, if tracks were obliterated two days
#' before the survey on some occassions, then `datestart` could be set to
#' `date - 2 days` for these observations and to `date - max.days` for the remainder.
#'
#' @reference.date is a character vector giving a nominal start date for the nesting season.
#' Should be in format %d%B or %B%d (e.g. '1 November', 'Dec 31', '12 Jul')
#' *Note* This is not the date when counts begin, but a date that defines where one nesting season ends and another begins.
#' In year round nesting it will typically be the date of a seasonal low point. In other cases it will occur in a
#' the middle of a period of zero counts (or assumed zeros) separating seasonal peaks. The value
#' is used to convert dates in `data` to the number of days since the start of the nesting season (`day`) which is
#' need for analysis.

#' @param min.obs is a single numeric giving the minimum number of observations that are needed
#' for a given beach in a given season to fit a model. Beaches with `< min.obs` counts will
#' be removed from the data
#'
#' @param sites.together A logical indicating whether counts on beaches in each season should be modelled together.
#' This is intended to allow users to estimate some joint parameters across sites (e.g. peak date) but the
#' initialization process is currently too slow with multiple beaches so it is recommended to leave sites.together = F.
#' @export

as_phenology = function(data, reference.date, max.days = 1, min.obs = 10, sites.together = F){

  if(!has_name(data,'date')) stop('Data must contain a column called date')
  data = tidyr::drop_na(data,date)
  if(!is(data$date,'Date')){
    data$date = suppressWarnings(as.Date(lubridate::parse_date_time(data$date,orders=c("Ymd","dmY"))))
    if(any_na(data$date)) stop("column 'date' is not of class 'Date' and format cannot be guessed")
  }

  data %<>% mutate(reference_date = set_ref_date(date,reference.date),
                   day = as.numeric(date - reference_date))

  if(has_name(data,'beach')) data$beach <- factor(data$beach)

  if(!has_name(data,'season')) {

    data %<>% group_by(reference_date) %>%
      mutate(season = max(year(date))) %>% rowwise()  %>%
      mutate(season = factor(paste(unique(c(year(reference_date),season)),collapse='-')))

  }

  # Order observations by date and check if any duplicates
  data %<>%
    group_by(across(any_of(c('season','beach')))) %>%
    dplyr::arrange(date)

  dups = subset(dplyr::count(data,date),n>1,select=-n)

  if(nrow(dups)){
    stop(paste(c('Duplicate survey dates found',
                 capture.output(print(data.frame(dups), row.names = FALSE))),
               collapse = "\n"))
  }

  # Calculate start of monitoring windows if not explicitly given
  if(!has_name(data,'datestart')) {

    data %<>%
      mutate(window = map_dbl(diff(c(-1,day)),~min(.x,max.days)), datestart = date-window)

  } else data %<>% mutate(window = date-datestart)

  # Remove beaches that have too few counts to be effectively modelled
  data %<>%
    mutate(too_few = n()<=min.obs) %>% ungroup()

  if(any(data$too_few)){
    warning(paste(c('The following seasons/beaches have insufficient data and have been removed:',
                    capture.output(print(data.frame(distinct(subset(data,too_few),season,beach)), row.names = FALSE))),
                  collapse = "\n"))
  }

  nest.vars = c('season','reference_date')
  if(!sites.together) nest.vars = c(nest.vars,'beach')

  # Prep final output
  data %<>% subset(!too_few) %>%
    droplevels() %>%
    dplyr::select(any_of(c('season','beach','day','datestart','date','window')),everything(),-too_few) %>%
    nest(data = -any_of(nest.vars)) %>%
    mutate(data = map2(data,reference_date, ~set_attr(.x,'ref_date',.y))) %>%
    dplyr::select(-reference_date)

  class(data) <- c('phenology_df',class(data))

  return(data)
}

# -----------------------------------------------------------------------------------
# set_ref_date
# -----------------------------------------------------------------------------------

#' Used internally by `MT_prep` to adjust dates using a reference date.

#' @param x A vector of class `Date` to be adjusted
#' @param ref_date A character vector giving the reference date in format %d%B or %B%d (e.g. '1 November', 'Dec 31')

set_ref_date =

  function(x,ref_date){

    stopifnot('First variable should be a Date vector' = is(x,'Date'))

    ref_date = as.Date(lubridate::parse_date_time(paste(ref_date,lubridate::year(x)),c('dBY','BdY')))

    dplyr::if_else(x >= ref_date,ref_date,ref_date - lubridate::years(1))

  }
