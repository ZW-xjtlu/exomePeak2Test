#assay_dds <- matrix(rpois(100, lambda = 10),nrow = 10,ncol = 10)
#IP_vec = assay_dds[,1]
#input_vec = assay_dds[,2]
mayers_hypertest <- function(IP_vec, input_vec, correct_IP = 0, correct_input = 0){

  IP_vec =  round(IP_vec + correct_IP)
  input_vec = round(input_vec + correct_input)
  IP_vec = pmax(IP_vec,0)
  input_vec = pmax(input_vec,0)
  total_white_balls <- sum(IP_vec)
  total_black_balls <- sum(input_vec)
  balls_draw <- IP_vec + input_vec
  white_balls_drawn <- IP_vec
  pvalue <- phyper(q = white_balls_drawn,
         m = total_white_balls,
         n = total_black_balls,
         k = balls_draw, lower.tail = F)
  p.adjust(pvalue, method = "BH")
}
