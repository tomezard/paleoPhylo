as.paleoPhylo <- function (nm, pn, st, en, xx = NA, label = nm, grp=NA) 
  {
  if(length(nm)!=length(unique(nm))) stop (paste("ID codes are not unique, specifically",nm[duplicated(nm)]))
  if(!is.na(pn[st==max(st)])) 
    {
    warning(paste("The oldest species' ancestor has been changed to 'NA' from", pn[st==max(st)], "."))
    pn[st==max(st)] <- NA
    }

  getX <- length(xx)==1 & is.na(xx[1])
  dat<- data.frame(nm,pn,st,en,xx,label,grp)
  dat <- dat[rev(order(dat$st, dat$en)),]
  pP <- list(nm = as.character(dat$nm), pn = as.character(dat$pn), 
    st = dat$st, en = dat$en,
    xx = dat$xx, label = as.character(dat$label), grp=dat$grp)
  class(pP) <- "paleoPhylo"
 
  #if I need to find the x locations
  if(getX) 
    {
    xxloc <- getXloc(pP)[, c(1, 6)]
    wxloc <- merge(dat, xxloc, by=c("nm"))	#with xlocs in
    dat <- wxloc
    dat <- with(wxloc, data.frame(nm, pn, st, en, xx=xx.y, label, grp))
    dat <- dat[rev(order(dat$st, dat$en)),]
    pP <- list(nm = as.character(dat$nm), pn = as.character(dat$pn), 
      st = dat$st, en = dat$en,
      xx = dat$xx, label = as.character(dat$label), grp=dat$grp)
    class(pP) <- "paleoPhylo"
   }
  
  return(pP)
  }
    