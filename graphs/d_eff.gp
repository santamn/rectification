set term lua tikz
set output "d_eff.tex"

set key top left offset 3, -2
set logscale x                 

set xlabel "$|f|$"    
set xrange [1:150]         
set ylabel "$D_{\\mathrm{eff}}$"
set yrange [0:10]

plot "../data/mono/d_eff_150_x10.dat" using 1:2 with points pt 7  title "Forward", \
      "../data/mono/d_eff_150_x10.dat" using 1:3 with points pt 12 title "Reverse"

unset terminal
