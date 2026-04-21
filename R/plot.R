# ---------------------------------------------------------------------------------------------------------------------------------------
# plot.MTpred: plot method for MTpred
# --------------------------------------------------------------------------------------------------------------------------------------- 
  
plot.phenology_df <- function(object,colour = c(nests='orange',activities='blue'),shape=21, x.scale = c('day','date'),
                              nrow = NULL,ncol = NULL){
 
raw = dplyr::select(object,any_of(c('season','beach','data'))) %>% 
      unnest(data) %>%
      gather(key = 'y.var',value='N.obs',activities,nests)

x.scale = match.arg(x.scale)  
  
pl = 
  ggplot(raw,aes(x = .data[[x.scale]])) + 
  scale_colour_manual(values = colour, name = NULL) + 
  labs(y = 'Count', x = x.scale)

  #pl = pl + 
  #geom_ribbon(aes(ymin = .lower,ymax = .upper,group = y.var),alpha=.3) +
  #geom_line(aes(colour=y.var)) + 

pl = pl + geom_point(data=raw,aes(y = N.obs,colour=y.var),shape=shape) 
  
vars = c('beach','season')
vars = vars[map_lgl(vars,~length(unique(object[[.x]]))>1)]
if(!length(vars)) return(pl)
xfree = ifelse('season' %in% vars & x.scale == 'date', 'free_x', 'fixed')
vars = paste('~',paste(vars,collapse='+'))

# Options for pagination if too many plots to view of one page
nplots = nrow(object)
if(!missing(ncol) & !missing(nrow)) pages = nplots(ncol*nrow) else pages = 1

if(pages<=1) print(pl + facet_wrap(as.formula(vars), scales = xfree)) else {
  
  for(i in 1:pages) print(pl+ ggforce::facet_wrap_paginate(as.formula(vars),scales = 'free_y',ncol=ncol,nrow=nrow,page=i)) 

  }
  
}
