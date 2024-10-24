#' A set of fetal movements data collected before and during the Covid-19 pandemic
#'
#' The data were recorded via 4D ultrasound scans from 40 fetuses (20 before Covid and 20 during Covid) at 32 weeks gestation,
#' and consist of the number of movements each fetus carries out in relation to the recordable scan length.
#'
#' @docType data
#'
#' @usage data(fetal_covid_data)
#'
#' @format An object of class \code{"data.frame"}
#' \describe{
#'  \item{UpperFaceMovements}{Inner Brow Raiser, Outer Brow Raiser, Brow Lower, Cheek Raiser, Nose Wrinkle.}
#'  \item{Headmovements}{Turn Right, Turn Left, Up, Down.}
#'  \item{MouthMovements}{Upper Lip Raiser, Nasolabial Furrow, Lip Puller,
#'      Lower Lip Depressor, Lip Pucker, Tongue Show, Lip Stretch, Lip Presser, Lip Suck, Lips Parting, Jaw Drop, Mouth Stretch.}
#'  \item{TouchMovements}{Upper Face, Side Face, Lower Face, Mouth Area.}
#'  \item{EyeBlink}{All scans were coded for eye blink.}
#'  \item{status_bi}{"during the pandemic" is coded by 1, "before the pandemic" is coded by 0.}
#'  \item{status}{specifies whether it is during or before the pandemic.}
#' }
#' @references Reissland, N., Ustun, B. and Einbeck, J. (2024).
#'             The effects of lockdown during the COVID-19 pandemic on fetal movement profiles.
#'             BMC Pregnancy and Childbirth, 24(1), 1-7.
#' @keywords datasets
#' @examples
#'
#' data(fetal_covid_data)
#' head(fetal_covid_data)
"fetal_covid_data"
