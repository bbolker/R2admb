## more problems:
##    parameter definitions (int number) --> %_% dummy?
##    are semicolons preserved?
##    C++-style output:  a << 1 << 2  --> %<% ?
##    for loops
##   this may just not be worth it ..

##'Write model to text file
##'
##'@usage write_model(model,con)
##'@param model An R function defining the model (see Details)
##'@param file name to write to (without extension)
##'@return Nothing: creates file in current working directory
##'@keywords misc
##'@note Based on \code{R2WinBUGS::write.model}
##'@details For convenience, you can write out your TPL file
##'in the form of an R function, with a couple of exceptions aimed
##'at reconciling C++ and R syntax:
##'\itemize{
##' \item{section headings such as PARAMETER_SECTION, LOCAL_CALCS, etc.
##'       should be preceded by two comment characters (##), at the
##'       beginning of a line (R will assume that any line beginning
##'       with two comment characters followed by a capital letter
##'       is a section heading)}
##' \item{comments should be otherwise be specified with comment
##'     characters rather than C++-style \code{//} comments (C-style
##'     \code{/* */} comments are not allowed)}
##'}
write.model <- function (model, file = "model", digits = 5) 
{
    if (tolower(substr(file,nchar(file)-4,nchar(file)))==".tpl") {
        con <- file
    } else con <- paste0(file,".tpl")
    model.text <- as.character(body(model))
    model.text <- gsub("#+","//",gsub("^##([A-Z]+)","\\1",model.text))
        model.text <- c("model", replaceScientificNotationR(body(model), 
            digits = digits))
    }
    writeLines(model.text, con = file)
}
