AMCMC.update <- function(draw,cur.mn,cur.var,cur.it){
	if(cur.it >0){
		mn <- ((cur.it-1)*cur.mn+draw)/cur.it
		if(cur.it==1){
			v <- matrix(0,nrow=length(draw),ncol=length(draw))
		} else {
			v <- (cur.it-2)*cur.var+(cur.it-1)*(cur.mn%*%t(cur.mn))+draw%*%t(draw)
			v <- (v-cur.it*(mn%*%t(mn)))/(cur.it-1)
		}
	} else {
		mn <- matrix(0,nrow=nrow(cur.mn),ncol=1)
		v <- matrix(0,nrow=nrow(draw),ncol=nrow(draw))
	}
	return(list(mn=mn,var=v))
}

