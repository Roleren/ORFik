% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shift_footprints_helpers.R
\name{shifts_save}
\alias{shifts_save}
\title{Save shifts for Ribo-seq}
\usage{
shifts_save(shifts, folder)
}
\arguments{
\item{shifts}{a list of data.table/data.frame objects.
Must be named with the full path to ofst/bam files that defines the shifts.}

\item{folder}{directory to save file,
Usually: file.path(libFolder(df), "pshifted"), where df is the ORFik
experiment / or your path of default file types.
It will be named file.path(folder, "shifting_table.rds").
For ORFik to work optimally,
the folder should be the /pshifted/ folder relative to default files.}
}
\value{
invisible(NULL), file saved to disc as "shifting_table.rds".
}
\description{
Should be stored in pshifted folder relative to default files
}
\examples{
df <- ORFik.template.experiment.zf()
shifts <- shifts_load(df)
original_shifts <- file.path(libFolder(df), "pshifted", "shifting_table.rds")
# Move to temp
new_shifts_path <- file.path(tempdir(), "shifting_table.rds")
new_shifts <- c(shifts, shifts)
names(new_shifts)[2] <- file.path(tempdir(), "RiboSeqTemp.ofst")
saveRDS(new_shifts, new_shifts_path)
new_shifts[[1]][1,2] <- -10
# Now update the new shifts, here we input only first
shifts_save(new_shifts[1], tempdir())
readRDS(new_shifts_path) # You still get 2 outputs

}
\seealso{
Other pshifting: 
\code{\link{changePointAnalysis}()},
\code{\link{detectRibosomeShifts}()},
\code{\link{shiftFootprints}()},
\code{\link{shiftFootprintsByExperiment}()},
\code{\link{shiftPlots}()},
\code{\link{shifts_load}()}
}
\concept{pshifting}
