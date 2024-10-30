setwd("C:/Users/Lu Mingyuan/Documents/2. NUS Mathematics")
setwd("./9. MA3220 Ordinary Differential Equations")
options(scipen=999);rm(list=ls(all.names=T))
opar=par(mar=c(bottom=0,left=2,top=0,right=2),mfrow=c(1,1),no.readonly=T) # Save init par
plot.new();set.seed(1)
## Method 1
library(ggplot2)
matrix_DE=function(t,state,params){
  x=state[1];y=state[2];
  dx=-2*x+y;dy=-x;
  return(list(c(dx,dy)));
}
# library(phaseR)
# phase_plot=flowField(deriv=matrix_DE,xlim=c(-5,5),ylim=c(-5,5), # args
#                      points=35,add=FALSE,system="two.dim", # kwargs
#                      xlab="u",ylab="v"); # kwargs
## Method 2
library(deSolve);library(comprehenr)
system=function(t,state,params){ # Modifiable based directly on matrix equation
  x=state[1];y=state[2];
  dx=-2*x+y;dy=-x;return(list(c(dx,dy)));
}
solution_path=function(init_conds,time_span){ # Path traced
  path=lapply(init_conds,function(init){ # init_conds arg for init param
    return(ode(y=init,times=time_span,func=system,parms=NULL)%>%
      as.data.frame()%>%rename(temp_x='1',temp_y='2')%>%
        mutate(init_x=init[1],init_y=init[2]))
  });
  return(bind_rows(path));
}
INIT_CONDITIONS=c(to_list(for (init_y in -4:4) c(-5,init_y)),
                  to_list(for (init_y in -4:4) c(5,init_y))) # from comprehenr pkg
TIME_SPAN=seq(from=0,to=10,by=0.1)
path_traced=solution_path(init_conds=INIT_CONDITIONS,time_span=TIME_SPAN)

FLOW_FIELD_WIDTH=1000;RESOLUTION_FACTOR=0.05;
resolution=FLOW_FIELD_WIDTH*RESOLUTION_FACTOR
colCount=FLOW_FIELD_WIDTH%/%resolution;rowCount=colCount;

phase_plane_grid=expand.grid(x=seq(-5,5,length.out=colCount),
                             y=seq(-5,5,length.out=rowCount))%>%
  mutate(x_prime=-2*x+y,y_prime=-x,
         angle=atan2(y_prime,x_prime),magnitude=sqrt(x_prime^2+y_prime^2))%>%
  mutate(xend=x+cos(angle)*0.2,yend=y+sin(angle)*0.5)#%>%
  #mutate(across(c(x,y,xend,yend),~.x*resolution))
plot2=ggplot(data=phase_plane_grid)+
  geom_segment(aes(x=x,y=y,xend=xend,yend=yend),
               arrow=arrow(length=unit(0.1,units="inches")),linewidth=0.3)+
  geom_path(data=path_traced,aes(x=temp_x,y=temp_y,group=interaction(init_x,init_y)),
                                 color="blue",linewidth=0.5)+
  geom_abline(intercept=0,slope=1,linetype="dashed",color="red")+
  scale_color_viridis_c(option="plasma",end=0.9)+
  theme_minimal(base_size=14)+
  theme(panel.border=element_rect(color="black",fill=NA),
        panel.grid=element_blank(),
        panel.background=element_rect(fill="white"))+
  labs(title="Phase Plane Diagram for u'(t)=-2u(t)+v(t), v'(t)=-u(t)",x="u(t)",y="v(t)")
plot(plot2)