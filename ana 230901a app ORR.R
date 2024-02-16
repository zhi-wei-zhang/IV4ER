# simulation study for estimating causal exposure-response relationship (23d)

# unadjusted analysis with no adjustment for confounding
est.ua = function(x, y, family=gaussian(), xg=NULL, ng=99) {
	mod = glm(y~x, family=family)
	if (is.null(xg)) xg = (1:ng)/(ng+1) else ng = length(xg)
	new.data = data.frame(x=xg)
	predict(mod, newdata=new.data, type="response")
}

# G-computation adjustment for measured confounding
est.gc = function(w, x, y, family=gaussian(), xg=NULL, ng=99) {
	mod = glm(y~w+x, family=family)
	if (is.null(xg)) xg = (1:ng)/(ng+1) else ng = length(xg)
	est = rep(NA, ng)
	for (g in 1:ng) {
		new.data = data.frame(w=w, x=rep(xg[g],length(w)))
		est[g] = mean(predict(mod, newdata=new.data, type="response"))
	}
	est
}

# IV analysis adjusting for all (measured and unmeasured) confounding
est.iv = function(z, x, y, family=gaussian(), xg=NULL, ng=99) {
	n = length(z)
	v = rep(NA, n)
	K = max(z)
	for (k in 1:K) {
		i.k = which(z==k)
		n.k = length(i.k)
		if (n.k>0) v[i.k] = rank(x[i.k])/(n.k+1)
	}
	v1 = qnorm(v)
	mod = glm(y~v1+x, family=family)
	if (is.null(xg)) xg = (1:ng)/(ng+1) else ng = length(xg)
	vg = (1:ng)/(ng+1)
	est = rep(NA, ng)
	for (g in 1:ng) {
		new.data = data.frame(v1=qnorm(vg), x=rep(xg[g],ng))
		est[g] = mean(predict(mod, newdata=new.data, type="response"))
	}
	est
}

dat = read.csv("Frey22 reformatted.csv", header=TRUE)
dim(dat); names(dat); summary(dat)
z = as.numeric(dat$Dose=="High")+1
x = log10(dat$Exposure)
y = 1-as.numeric(dat$Response=="NR") # as.numeric(dat$Response=="CR") # 
n = length(y)

x.lo = 7.05; x.hi = 8.65
x.rge = x.hi-x.lo
ng = 99
xg = x.lo+x.rge*(1:ng)/(ng+1)
nr = 1000
nm = 2

pt.est = matrix(NA, ng, nm)
pt.est[,1] = est.ua(x, y, family = binomial(), xg=xg)
pt.est[,2] = est.iv(z, x, y, family = binomial(), xg=xg)
pt.est

bt.est = array(NA, dim=c(nr, nm, ng))
bt.se = matrix(NA, ng, nm)
bt.good = bt.se
for (r in 1:nr) {
	if (r%%100==0) print(paste(r, date()))
	ind = sample(1:n, n, replace=TRUE)
	zb = z[ind]; xb = x[ind]; yb = y[ind]
	bt.est[r,1,] = est.ua(xb, yb, family = binomial(), xg=xg)
	bt.est[r,2,] = est.iv(zb, xb, yb, family = binomial(), xg=xg)
}

for (m in 1:nm) {
	for (g in 1:ng) {
		good = !is.na(bt.est[,m,g])
		bt.good[g,m] = sum(good)
		bt.se[g,m] = sd(bt.est[good,m,g])
	}
}

est.smry = cbind(xg, pt.est, bt.se, bt.good)
est.smry

file.out = "ana rst 230901a app ORR.csv"
write.csv(est.smry, file=file.out)

# sim.smry = read.csv(file.out)
za = qnorm(0.975)

plot(xg, pt.est[,1], type="l", lty=1, lwd=2, col=1, xlab="x", ylab="Overall response rate", ylim=c(0,0.8))
lines(xg, pt.est[,1]-za*bt.se[,1], type="l", lty=1, lwd=1, col=1)
lines(xg, pt.est[,1]+za*bt.se[,1], type="l", lty=1, lwd=1, col=1)
lines(xg, pt.est[,2], type="l", lty=3, lwd=2, col=3)
lines(xg, pt.est[,2]-za*bt.se[,2], type="l", lty=3, lwd=1, col=3)
lines(xg, pt.est[,2]+za*bt.se[,2], type="l", lty=3, lwd=1, col=3)
legend("top", lty=c(1,3), lwd=2, col=c(1,3), legend=c("Unadjusted", "IV"))

# end of program