# procedure to obtain an optimal one break partitions for a segment that
# starts at date start and ends at date last. It returns the optimal break
# and the associated SSR.
# Procedure used with the sequential method whenthe T*(T+1)/2 vector of SSR
# is not computed.
# start: beginning of the segment considered
# b1: first possible break date
# b2: last possible break date
# last end of segment considered
partione = function(b1,b2,last,vssr,vssrev){

dvec = matrix(0L,nrow = last, ncol = 1)

j = b1;

while (j <= b2) {
dvec[j,1]=vssr[j,1]+vssrev[last-j,1];
j=j+1;
}

ssrmin=min(dvec[seq(b1,b2),1]);
dx=(b1-1)+which.min(dvec[seq(b1,b2),1]);
out = list('ssrmin' = ssrmin, 'dx' = dx)
return(out)
}
