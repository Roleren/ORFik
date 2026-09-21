# Private cache ownership and storage-only retries. Never infer ownership from
# a directory name alone, and never recursively remove the output directory.
.ofst_storage_abort <- function(...) {
  stop(structure(list(message = paste0("OFST merge: ", ...), call = NULL),
                 class = c("ofst_scratch_error", "error", "condition")))
}

.ofst_memory_error <- function(e) {
  grepl("cannot allocate|bad_alloc|vector memory|long vectors not supported|negative length vectors",
        conditionMessage(e), ignore.case = TRUE)
}

.ofst_read_line <- function(path) {
  tryCatch({
    x <- suppressWarnings(readLines(path, n = 1L, warn = FALSE))
    if (length(x) == 1L && nzchar(x)) x else NULL
  }, error = function(e) NULL)
}

.ofst_owner_text <- function(x, pattern = ".+") {
  is.character(x) && is.null(dim(x)) && length(x) == 1L && !is.na(x) && grepl(pattern, x)
}

.ofst_process_start <- function(pid) {
  stat <- .ofst_read_line(file.path("/proc", as.character(pid), "stat"))
  if (is.null(stat)) return(NULL)
  # The process name in parentheses can itself contain spaces/parentheses.
  fields <- strsplit(sub("^.*\\) ", "", stat), " +")[[1L]]
  if (length(fields) >= 20L && .ofst_owner_text(fields[20L], "^[0-9]+$")) fields[20L] else NULL
}

.ofst_process_owner <- function() {
  machine <- .ofst_read_line("/etc/machine-id")
  list(pid = Sys.getpid(), start = .ofst_process_start(Sys.getpid()),
       host = unname(Sys.info()["nodename"]), user = unname(Sys.info()["effective_user"]),
       machine = if (is.null(machine)) NULL else digest::digest(machine, algo = "sha256", serialize = FALSE),
       boot = .ofst_read_line("/proc/sys/kernel/random/boot_id"),
       pid_namespace = suppressWarnings(Sys.readlink("/proc/self/ns/pid")))
}

.ofst_owner_alive <- function(owner) {
  current <- .ofst_process_owner()
  if (!is.list(owner) || !identical(owner$host, current$host) ||
      !identical(owner$user, current$user)) return(NA)
  if (!is.null(owner$machine) && !is.null(current$machine) &&
      !identical(owner$machine, current$machine)) return(NA)
  uuid <- "^[[:xdigit:]]{8}-[[:xdigit:]]{4}-[[:xdigit:]]{4}-[[:xdigit:]]{4}-[[:xdigit:]]{12}$"
  if (!.ofst_owner_text(owner$boot, uuid) || !.ofst_owner_text(current$boot, uuid)) return(NA)
  if (!identical(owner$boot, current$boot)) {
    # A reboot proves death only when a stable machine identity also matches.
    if (!is.null(owner$machine) && identical(owner$machine, current$machine)) return(FALSE)
    return(NA)
  }
  if (!.ofst_owner_text(owner$pid_namespace, "^pid:\\[[0-9]+\\]$") ||
      !identical(owner$pid_namespace, current$pid_namespace) || is.null(current$start)) return(NA)
  if (!is.numeric(owner$pid) || is.object(owner$pid) || !is.null(dim(owner$pid)) ||
      length(owner$pid) != 1L || is.na(owner$pid) || !is.finite(owner$pid) ||
      owner$pid < 1 || owner$pid != trunc(owner$pid) ||
      !.ofst_owner_text(owner$start, "^[0-9]+$")) return(NA)
  if (!dir.exists(file.path("/proc", as.character(owner$pid)))) return(FALSE)
  observed <- .ofst_process_start(owner$pid)
  if (is.null(observed)) return(NA)
  identical(observed, owner$start) # Handles PID reuse, not just PID existence.
}

.ofst_cache_owner <- function(path) {
  if (!dir.exists(path) || !identical(Sys.readlink(path), "")) return(NULL)
  marker <- file.path(path, ".ofst-cache-owner.rds")
  if (!file.exists(marker) || !identical(Sys.readlink(marker), "")) return(NULL)
  owner <- tryCatch(readRDS(marker), error = function(e) NULL)
  canonical <- normalizePath(path, mustWork = TRUE)
  if (!is.list(owner) || !identical(owner$format, "ORFik-private-ofst-cache-v1") ||
      !identical(owner$id, basename(canonical)) || !identical(owner$parent, dirname(canonical)) ||
      !grepl("^orfik-ofst-(filter|anchor|output|tables)-", owner$id)) return(NULL)
  owner
}

.ofst_write_cache_owner <- function(path, owner) {
  temporary <- tempfile(".owner-", tmpdir = path)
  on.exit(unlink(temporary), add = TRUE)
  tryCatch({
    saveRDS(owner, temporary)
    if (!file.rename(temporary, file.path(path, ".ofst-cache-owner.rds")))
      stop("could not publish the ownership record")
  }, error = function(e) .ofst_storage_abort("cannot record cache ownership in '", path,
                                             "'. Check disk space/permissions. ", conditionMessage(e)))
  invisible(NULL)
}

.ofst_cleanup_cache <- function(path, current_run = FALSE, expected_id = NULL,
                                 discard_recovery = FALSE) {
  suspendInterrupts(tryCatch({
    if (!dir.exists(path)) return(TRUE)
    owner <- .ofst_cache_owner(path)
    if (is.null(owner) || (!is.null(expected_id) && !identical(owner$id, expected_id))) return(FALSE)
    if (current_run) {
      if (!identical(owner$process, .ofst_process_owner())) return(FALSE)
    } else if (!identical(.ofst_owner_alive(owner$process), FALSE)) return(FALSE)
    if (!identical(owner$state, "temporary") && !(current_run && discard_recovery)) {
      message("OFST cache: preserving output recovery files in '", path,
              "' (publication was interrupted); inspect these before removing the cache.")
      return(FALSE)
    }
    # The durable output-side record can find primary scratch even when R's
    # next session has a different tempdir. Only mutually linked owned caches
    # qualify; arbitrary paths in a malformed marker are never deleted.
    for (child in owner$children) {
      if (!is.list(child) || !is.character(child$path) || length(child$path) != 1L) return(FALSE)
      if (!dir.exists(child$path)) next
      child_owner <- .ofst_cache_owner(child$path)
      if (is.null(child_owner) || !identical(child_owner$id, child$id) ||
          !grepl("^orfik-ofst-(filter|tables)-", child_owner$id) || length(child_owner$children) ||
          !identical(child_owner$anchor_id, owner$id) ||
          !identical(child_owner$process, owner$process)) return(FALSE)
      if (!.ofst_cleanup_cache(child$path, current_run, child$id)) return(FALSE)
    }
    removed <- unlink(path, recursive = TRUE)
    if (removed != 0L || dir.exists(path)) {
      message("OFST cache cleanup could not remove '", path, "'; remove this owned cache later when storage is accessible.")
      return(FALSE)
    }
    if (!current_run) message("OFST cache: removed abandoned cache '", path, "'.")
    TRUE
  }, error = function(e) {
    message("OFST cache cleanup failed for '", path, "': ", conditionMessage(e))
    FALSE
  }))
}

.ofst_cleanup_stale_caches <- function(parent) {
  if (!dir.exists(parent)) return(invisible(NULL))
  paths <- list.files(parent, pattern = "^orfik-ofst-(filter|anchor|output|tables)-", full.names = TRUE)
  for (path in paths) .ofst_cleanup_cache(path)
  invisible(NULL)
}

.ofst_new_cache <- function(parent, kind = "filter") {
  if (!dir.exists(parent) || file.access(parent, 2L) != 0L)
    .ofst_storage_abort("scratch directory '", parent, "' does not exist or is not writable. Set filter_tmpdir to a disk with free space.")
  .ofst_cleanup_stale_caches(parent)
  path <- tempfile(paste0("orfik-ofst-", kind, "-"), tmpdir = normalizePath(parent, mustWork = TRUE))
  if (!dir.create(path, showWarnings = FALSE)) .ofst_storage_abort("could not create scratch directory '", path, "'.")
  registered <- FALSE
  # This exact directory was just created by this function, before any payload
  # was written. Clean even a failed ownership-record write (e.g. full disk).
  on.exit(if (!registered) suspendInterrupts(unlink(path, recursive = TRUE)), add = TRUE)
  owner <- list(format = "ORFik-private-ofst-cache-v1", id = basename(path),
                parent = dirname(path), process = .ofst_process_owner(),
                created = as.character(Sys.time()), state = "temporary", children = list())
  .ofst_write_cache_owner(path, owner)
  registered <- TRUE
  path
}

.ofst_link_cache <- function(anchor, child) {
  parent_owner <- .ofst_cache_owner(anchor)
  child_owner <- .ofst_cache_owner(child)
  if (is.null(parent_owner) || is.null(child_owner)) .ofst_storage_abort("cache ownership record is missing; cannot safely register scratch cleanup.")
  child_owner$anchor_id <- parent_owner$id
  .ofst_write_cache_owner(child, child_owner)
  parent_owner$children <- c(parent_owner$children, list(list(path = child, id = child_owner$id)))
  .ofst_write_cache_owner(anchor, parent_owner)
}

.ofst_with_scratch <- function(primary, fallback, fun) {
  anchor <- NULL
  anchor_error <- NULL
  same_directory <- !is.null(fallback) && identical(normalizePath(primary, mustWork = FALSE),
                                                   normalizePath(fallback, mustWork = FALSE))
  if (!is.null(fallback) && !same_directory)
    anchor <- tryCatch(.ofst_new_cache(fallback, "anchor"),
                       ofst_scratch_error = function(e) { anchor_error <<- e; NULL })
  on.exit(if (!is.null(anchor)) .ofst_cleanup_cache(anchor, current_run = TRUE), add = TRUE)
  attempt <- function() {
    scratch <- .ofst_scratch_dir(primary)
    on.exit(.ofst_cleanup_cache(scratch, current_run = TRUE), add = TRUE)
    if (!is.null(anchor)) .ofst_link_cache(anchor, scratch)
    fun(scratch)
  }
  result <- tryCatch(attempt(), ofst_scratch_error = function(e) e)
  if (!inherits(result, "ofst_scratch_error")) return(result)
  if (is.null(fallback) || same_directory)
    .ofst_storage_abort(conditionMessage(result), " No distinct filter_fallback_dir is configured; choose a disk with free space.")
  if (is.null(anchor))
    .ofst_storage_abort("primary scratch failed: ", conditionMessage(result),
                        " Fallback directory '", fallback, "' is unavailable: ", conditionMessage(anchor_error))
  message("OFST rescue: primary scratch storage failed: ", conditionMessage(result))
  message("OFST rescue: retrying once in output-side cache '", anchor,
          "'. Input files will be read again; completed outputs are not deleted.")
  tryCatch(fun(anchor), ofst_scratch_error = function(e)
    .ofst_storage_abort("both scratch locations failed. Primary: ", conditionMessage(result),
                        " Fallback: ", conditionMessage(e),
                        " Free disk space/check quota and inodes before retrying."))
}
