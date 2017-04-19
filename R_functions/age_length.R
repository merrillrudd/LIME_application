age_length <- function(highs, lows, L_a, CVlen){

	lbprobs <- function(mnl, sdl) return(pnorm(highs, mnl, sdl) - pnorm(lows, mnl, sdl))
	vlprobs <- Vectorize(lbprobs, vectorize.args = c("mnl", "sdl"))
	plba <- t(vlprobs(L_a, L_a * CVlen))
	plba <- plba/rowSums(plba)
	return(plba)
}