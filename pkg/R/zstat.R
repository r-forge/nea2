zstat <-
function(BSTAR)
{
 bstar.mean = apply(BSTAR,c(1,2),mean)
 dev = c(BSTAR) - c(bstar.mean)
    dev = array(dev, dim(BSTAR))
 stat.p = apply(dev,c(1,3),max)
 stat.m = apply(dev,c(1,3),function(x) -max(-x))
  stat = ifelse(stat.p>stat.m, stat.p, -stat.m)
  mstat = apply(stat,1,mean)
  sdstat= apply(stat,1,sd)
  Z = (stat - mstat)/sdstat
  return(Z)
}

