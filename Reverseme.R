reverseComplement = function(input, rev = T){
  input = chartr("ATGCatgc","TACGtacg", input)
  if (rev==T){
    reversed = strsplit(as.character(input), split ="")
    input = reversed[[1]][nchar(input):1]
  }
  return(paste(input,collapse = ""))
}
reverseComplement("gatggcggacttcgggttatcta",T)
TAGAAGATAAGTCGTAGGGATGC
hello
