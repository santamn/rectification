set term lua tikz
set output "alpha.tex"

unset key

set xlabel "$f$"             
set ylabel "$\\alpha$"

plot "../data/alpha.dat" using 1:2 with points pt 7 

unset terminal
