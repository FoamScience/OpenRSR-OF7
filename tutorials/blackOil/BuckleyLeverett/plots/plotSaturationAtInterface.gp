# Needs case index to be passed from CMD
if (!exists("case")) case = 1
# Needs solver name to passed from CMD
if (!exists("solver")) solver = "Impes"
# Needs time at which the frontal saturation is at domain1 boundary
if (!exists("time")) time = 86400*0.55

field = "water.alpha"
set term pngcairo enhanced
set out sprintf("%s-%s-profile-at-%.2f-days.png", solver, field, time/86400.0)

set title sprintf("%s profile from injection point", field)

set xlabel "Distance [m]"
set ylabel field
set xrange [0:12]
set yrange [0:1]

set key outside vertical font key_font

tf=sprintf("theory/case%d/%s/%s-%.2f-days.csv", case, field, field, time/86400.0) 
df=sprintf("postProcessing/singleGraph/%d/line_%s_p.xy", time, field)
pti(name)=sprintf("%s (%.2f Days)", name, time/86400.0)
plot tf t pti("Theory") w l dt 4, df t pti(solver) w l
