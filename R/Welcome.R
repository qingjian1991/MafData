#' Welcome
#'
#' This is my First R Pakcges, The Input is your names.
#'
#'@author Qingjian Chen (Sun Yat-sen University)
#'@details This is my first R Package
#'@param name user names
#'@return no return information.
#'@export

Welcome = function(name =""){
  message(sprintf("Happy %s", name))
}


