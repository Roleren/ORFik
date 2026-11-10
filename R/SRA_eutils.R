#' Read an NCBI E-utilities URL as XML, within NCBI's usage rules
#'
#' All ORFik calls to NCBI E-utilities (esearch, efetch, esummary, elink) go through this function.
#' It follows the E-utilities usage guidelines (\url{https://www.ncbi.nlm.nih.gov/books/NBK25497/}):
#' \itemize{
#'  \item adds \code{tool=ORFik}, the contact e-mail and the NCBI API key when they are set;
#'  \item spaces requests: at most 3 per second without an API key, 10 per second with one (per R process);
#'  \item retries on HTTP 429 (too many requests), 5xx and connection failures with exponential back-off,
#'  honouring a \code{Retry-After} header;
#'  \item sends long requests (for example efetch with hundreds of ids) as HTTP POST instead of an over-long URL.
#' }
#' Set the key and e-mail once per session or in ~/.Renviron:
#' \code{Sys.setenv(ENTREZ_KEY = "<key>", ENTREZ_EMAIL = "you@example.org")} (\code{NCBI_API_KEY} is also read).
#' A free API key is created under "API Key Management" in the settings of an NCBI account.
#' @param url character, a full E-utilities URL
#' @param max_tries integer, attempts before giving up
#' @return an xml2 document
#' @noRd
eutils_read_xml <- function(url, max_tries = 6) {
  key <- Sys.getenv("ENTREZ_KEY", Sys.getenv("NCBI_API_KEY"))
  email <- getOption("ORFik.entrez_email", Sys.getenv("ENTREZ_EMAIL"))
  base <- sub("\\?.*", "", url)
  query <- if (grepl("\\?", url)) sub("^[^?]*\\?", "", url) else ""
  params <- strsplit(query, "&", fixed = TRUE)[[1]]
  params <- params[nzchar(params)]
  fields <- setNames(as.list(vapply(gsub("+", " ", sub("^[^=]*=?", "", params), fixed = TRUE), utils::URLdecode, "")), sub("=.*", "", params))
  if (is.null(fields$tool)) fields$tool <- "ORFik"
  if (nzchar(email) && is.null(fields$email)) fields$email <- email
  if (nzchar(key) && is.null(fields$api_key)) fields$api_key <- key
  post <- nchar(url) > 1500
  min_interval <- if (nzchar(key)) 0.11 else 0.34
  for (i in seq_len(max_tries)) {
    wait <- min_interval - (as.numeric(Sys.time()) - eutils_state$last)
    if (wait > 0) Sys.sleep(wait)
    eutils_state$last <- as.numeric(Sys.time())
    res <- tryCatch(if (post) httr::POST(base, body = fields, encode = "form", httr::timeout(120))
                    else httr::GET(base, query = fields, httr::timeout(120)),
                    error = function(e) e)
    status <- if (inherits(res, "error")) NA_integer_ else httr::status_code(res)
    if (!is.na(status) && status == 200) {
      txt <- httr::content(res, as = "raw")
      doc <- tryCatch(xml2::read_xml(txt), error = function(e) NULL)
      # NCBI sometimes answers 200 with an <ERROR>API rate limit exceeded</ERROR> body
      if (!is.null(doc) && !grepl("rate limit exceeded", as.character(doc), ignore.case = TRUE)) return(doc)
    } else if (!is.na(status) && !(status == 429 || status >= 500)) {
      stop("NCBI E-utilities returned HTTP ", status, " for ", base)
    }
    if (i == max_tries) break
    retry_after <- if (!inherits(res, "error")) suppressWarnings(as.numeric(httr::headers(res)[["retry-after"]])) else NA
    delay <- if (length(retry_after) == 1 && !is.na(retry_after)) retry_after else min(2^i, 60) + stats::runif(1)
    message("NCBI E-utilities busy (", if (is.na(status)) conditionMessage(res) else paste("HTTP", status),
            "), retry ", i, "/", max_tries - 1, " in ", round(delay, 1), " s",
            if (!nzchar(key)) " (an NCBI API key raises the limit: see ?eutils_read_xml)")
    Sys.sleep(delay)
  }
  stop("NCBI E-utilities did not answer after ", max_tries, " attempts: ", base)
}

eutils_state <- new.env(parent = emptyenv())
eutils_state$last <- 0
