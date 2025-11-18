set term lua tikz
set output "mu.tex"

set key top left offset 3, -2         
set logscale x                 

set xlabel "$|f|$"             
set ylabel "$\\mu(f)$"         

plot "../data/mu.dat" using 1:2 with points pt 7  title "Forward", \
      "../data/mu.dat" using 1:3 with points pt 12 title "Reverse"

unset terminal
