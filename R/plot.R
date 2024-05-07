# ---------------------------------------------------------------------------------------------------------------------------------------
# plot.MTpred: plot method for MTpred
# --------------------------------------------------------------------------------------------------------------------------------------- 
  
plot.phenology_df <- function(object,colour = c(nests='orange',activities='blue'),shape=21, x.scale = c('day','date')){
 
raw = dplyr::select(object,any_of(c('season','beach','data'))) %>% 
      unnest(data) %>%
      gather(key = 'y.var',value='N.obs',activities,nests)
  
pl = 
  ggplot(raw,aes(x = .data[[x.scale]])) + 
  scale_colour_manual(values = colour, name = NULL) + 
  labs(y = 'Count', x = x.scale)

  #pl = pl + 
  #geom_ribbon(aes(ymin = .lower,ymax = .upper,group = y.var),alpha=.3) +
  #geom_line(aes(colour=y.var)) + 

vars = c('beach','season')
vars = vars[map_lgl(vars,~length(unique(object[[.x]]))>1)]
xfree = ifelse('season' %in% vars & x.scale == 'date', 'free_x', 'fixed')
vars = paste('~',paste(vars,collapse='+'))

pl + 
geom_point(data=raw,aes(y = N.obs,colour=y.var),shape=shape) +
facet_wrap(vars, scales = xfree) 
  
}
