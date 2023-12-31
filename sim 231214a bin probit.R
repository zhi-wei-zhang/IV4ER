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

ed.50 = 1
dose.rge = 1L # 2L
dd = ed.50*(2^((-dose.rge):dose.rge))
nd = length(dd)
pz = rep(1/nd, nd)
beta = c(0, 1, 1, 0.5) # coefficients of (1,X,W,U) in E(Y|X,W,U)
gamma = c(0.5, 1) # coefficients of (W,U) in multiplier to d: eta = exp((W,U)gamma)
ng = 99; x.lo = -2; x.hi = 2
xg = x.lo+(x.hi-x.lo)*(1:ng)/(ng+1)
yg.tru = pnorm((beta[1]+beta[2]*xg)/sqrt(1+beta[3]^2+beta[4]^2))
n = 30 # 120 # 60 # 
nr = 1000
nm = 3
one = rep(1, n)

all.est = array(NA, dim=c(nr, nm, ng))
est.bias = matrix(NA, ng, nm)
est.sd = est.bias
est.mse = est.bias
est.good = est.bias
for (r in 1:nr) {
	if (r%%100==0) print(paste(r, date()))
	w = rnorm(n); u = rnorm(n)
	eta = as.numeric(cbind(w,u)%*%gamma)
	z = sample(1:nd, n, replace=TRUE)
	d = dd[z]
	x = z-2+0.5*eta # as.vector(eta*d/(ed.50+(eta*d)))
	mu.y = pnorm(cbind(one,x,w,u)%*%beta)
	y = rbinom(n, 1, mu.y) # +rnorm(n)
	all.est[r,1,] = est.ua(x, y, family = binomial(link="probit"), xg=xg)
	all.est[r,2,] = est.gc(w, x, y, family = binomial(link="probit"), xg=xg)
	all.est[r,3,] = est.iv(z, x, y, family = binomial(link="probit"), xg=xg)
}

for (m in 1:nm) {
	for (g in 1:ng) {
		good = !is.na(all.est[,m,g])
		est.good[g,m] = sum(good)
		est.bias[g,m] = mean(all.est[good,m,g])-yg.tru[g]
		est.sd[g,m] = sd(all.est[good,m,g])
		est.mse[g,m] = (est.bias[g,m]^2)+(est.sd[g,m]^2)
	}
}

sim.smry = cbind(yg.tru, est.bias, est.sd, est.mse, est.good)
sim.smry = rbind(colMeans(sim.smry), sim.smry)
sim.smry

file.out = paste("sim rst 231214a bin probit n=", n, ".csv", sep="")
write.csv(sim.smry, file=file.out)

# sim.smry = read.csv(file.out)

# bias plot
# est.bias = sim.smry[-1,2:4]
plot(xg, est.bias[,1], type="l", xlab="x", ylab="Bias", main=paste("n =", n))
for (i in 2:3) lines(xg, est.bias[,i], type="l", lty=i, col=i)
legend("top", lty=1:3, col=1:3, legend=c("UA", "GC", "IV"))

# sd plot
plot(xg, est.sd[,1], type="l", xlab="x", ylab="SD", ylim=c(min(est.sd), max(est.sd)))
for (i in 2:3) lines(xg, est.sd[,i], type="l", lty=i, col=i)
legend("top", lty=1:3, col=1:3, legend=c("UA", "GC", "IV"))

# rmse plot
est.rmse = sqrt(est.mse)
plot(xg, est.rmse[,1], type="l", xlab="x", ylab="RMSE", ylim=c(min(est.rmse), max(est.rmse)))
for (i in 2:3) lines(xg, est.rmse[,i], type="l", lty=i, col=i)
legend("top", lty=1:3, col=1:3, legend=c("UA", "GC", "IV"))


# end of program