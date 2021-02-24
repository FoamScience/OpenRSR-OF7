# Needs case index to be passed from CMD
if (!exists("case")) case = 1
# Needs solver name to be passed from CMD
if (!exists("solver")) solver = "Impes"
# Needs solver name to be passed from CMD
if (!exists("start")) start = 0.0

field = "p"
key_font = "Iosevka Term, 9"
set term pngcairo enhanced
set out sprintf("%s-%s-profile-%dTo%d.png", solver, field, (start+0.1)*10, (start+0.5)*10)

set title sprintf("%s profile from injection point", field)

set xlabel "Distance [m]"
set ylabel field
set xrange [0:12]
set yrange [0:8e7]

set key outside vertical font key_font
array times[5]
do for [i=1:5] { times[i] = (start+0.1*i)*86400 }

tf(n)=sprintf("theory/case%d/%s/%s-%.2f-days.csv", case, field, field, times[n]/86400.0) 
df(n)=sprintf("postProcessing/singleGraph/%d/line_water.alpha_%s.xy", times[n], field)
pti(name,n)=sprintf("%s (%.2f Days)", name, times[n]/86400.0)
plot for [i=1:5] tf(i) notitle w l dt 4, \
     for [i=1:5] df(i) u 1:3 notitle w l, \
     1/0 t "Theory" dt 4 lc -1, 1/0 t solver lc -1, \
     for [i=1:5] 1/0 ls (i+5) t sprintf("%.2f Days", times[i]/86400.0)


