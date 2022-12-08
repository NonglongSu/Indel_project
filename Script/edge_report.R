suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(tidyverse))

# setwd("~/Dropbox (ASU)/Indel_project/Script")



markup = function(x, y){
  
  x1 = str_replace_all(x, fixed(" "), "")
  y1 = str_replace_all(y, fixed(" "), "")
  
  x2 = str_split(x1, boundary("character"))[[1]]
  y2 = str_split(y1, boundary("character"))[[1]]
  
  pos = which(x2 != y2 & y2 != '-' & y2 != '*')   
  
  x2[pos] = str_c("<span style=\"color:red;\">",x2[pos],"</span>")
  y2[pos] = str_c("<span style=\"color:red;\">",y2[pos],"</span>")
  
  x3 = str_c(x2, collapse = "")
  y3 = str_c(y2, collapse = "")
  
  return(c(x3, y3))
}

##################################################################################
# file.1 = "../test_human_mouse_rat/Data_6/Figure/mafft_edge/align.ori.l.txt"
# file.2 = "../test_human_mouse_rat/Data_6/Figure/mafft_edge/align.better.l.txt"
# file.3 = "../test_human_mouse_rat/Data_6/Figure/mafft_edge//align.ori.r.txt"
# file.4 = "../test_human_mouse_rat/Data_6/Figure/mafft_edge/align.better.r.txt"
main = function(file.1, file.2, file.3, file.4, ouFile){
  
  data.1 = read_delim(file.1, "\t", col_names = TRUE)
  data.2 = read_delim(file.2, "\t", col_names = TRUE)
  data.3 = read_delim(file.3, "\t", col_names = TRUE)
  data.4 = read_delim(file.4, "\t", col_names = TRUE)
  
  chain.1.1 = c()
  chain.1.2 = c()
  chain.2.1 = c()
  chain.2.2 = c()
  
  chain.3.1 = c()
  chain.3.2 = c()
  chain.4.1 = c()
  chain.4.2 = c()
  
  for (i in 1:nrow(data.1)) {
    chain.1.1 = c(chain.1.1, markup(data.1$amino_acid.1[i], data.1$amino_acid.2[i]))
    chain.1.2 = c(chain.1.2, markup(data.1$wid.1[i], data.1$wid.2[i]))
    
    chain.2.1 = c(chain.2.1, markup(data.2$amino_acid.1[i], data.2$amino_acid.2[i]))
    chain.2.2 = c(chain.2.2, markup(data.2$wid.1[i], data.2$wid.2[i]))
  }
  
  for (i in 1:nrow(data.3)) {
    chain.3.1 = c(chain.3.1, markup(data.3$amino_acid.1[i], data.3$amino_acid.2[i]))
    chain.3.2 = c(chain.3.2, markup(data.3$wid.1[i], data.3$wid.2[i]))
    
    chain.4.1 = c(chain.4.1, markup(data.4$amino_acid.1[i], data.4$amino_acid.2[i]))
    chain.4.2 = c(chain.4.2, markup(data.4$wid.1[i], data.4$wid.2[i]))
  }
  
  
# Embed HTML via concatenation
sink(file = ouFile, type = c("output", "message"), append = FALSE)
cat("<html xmlns=\"http://www.w3.org/1999/xhtml\">
  <head>
    <style>
      td { font-family: monospace; }
      table,td { border: 1px solid black; }
      table {border-collapse: collapse;}
      .td3 {border-left: 5px double;}
      
    </style>
  </head>
  <body>
    <table align=\"center\">
      <tr>
        <th></th>
        <th>Original Alignment</th>
        <th></th>
        <th>Better Alignment</th>
      </tr>","\n")

for(i in seq(1, nrow(data.1) * 2, 2) ){
cat("      <tr>
       <td>"                  , chain.1.1[i], "<br/>", chain.1.1[i+1], "</td>", "\n",
    "      <td>"              , chain.1.2[i], "<br/>", chain.1.2[i+1], "</td>", "\n",
    "      <td class=\"td3\">", chain.2.1[i], "<br/>", chain.2.1[i+1], "</td>", "\n",
    "      <td>"              , chain.2.2[i], "<br/>", chain.2.2[i+1], "</td>", "\n",
    "      </tr>","\n")
}

cat("   </table>
    <div style=\"clear:both; margin: 50px\"></div>   
    <table align=\"center\">
      <tr>
        <th></th>
        <th>Original Alignment</th>
        <th></th>
        <th>Better Alignment</th>
      </tr>","\n")

for(i in seq(1, nrow(data.3) * 2, 2) ){
cat("      <tr>
       <td>"                  , chain.3.1[i], "<br/>", chain.3.1[i+1], "</td>", "\n",
      "    <td>"              , chain.3.2[i], "<br/>", chain.3.2[i+1], "</td>", "\n",
      "    <td class=\"td3\">", chain.4.1[i], "<br/>", chain.4.1[i+1], "</td>", "\n",
      "    <td>"              , chain.4.2[i], "<br/>", chain.4.2[i+1], "</td>", "\n",
      "    </tr>","\n")
}


cat("    </table>
  </body>
</html>")
sink(NULL)

closeAllConnections()

}

args = commandArgs(trailingOnly = TRUE)
main(args[1], args[2], args[3], args[4], args[5])



