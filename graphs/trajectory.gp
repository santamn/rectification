# gnuplot script: animated trajectory + moving rod + channel boundaries
# Data format (columns): step  x  y  theta[radian]
# Example line: 123  0.0123  1.0456  -0.37

file = "data/di/trajectory_f1_seed0.dat"
out  = "data/di/trajectory.gif"


# rod length
L = 0.01

# animation controls
speed_factor = 10       # 10x faster playback (time advances 10x per frame)
frame_stride = 200 * speed_factor      # larger -> fewer frames -> faster to render
frame_delay  = 3        # 1/100 sec per frame

# plot window controls
x_window = 1.0          # show x-range with width = 1 (centered on data)

# channel boundary: omega(x)
omega(x) = sin(2.0*pi*x) + 0.25*sin(4.0*pi*x) + 1.12

set term gif animate optimize delay frame_delay size 900,700
set output out

set size ratio -1
set key off
set samples 2000
set tics out

# determine x-range from data (with margin)

stats file using 2 name "X" nooutput
if (!exists("X_records") || X_records < 1) {
    print sprintf("ERROR: failed to read data file: %s", file)
    exit
}
xmin_data = X_min
xmax_data = X_max

# Force the displayed x-range width to x_window, centered on the data.
xmid = 0.5*(xmin_data + xmax_data)
xmin = xmid - 0.5*x_window
xmax = xmid + 0.5*x_window
xspan = xmax - xmin
set xrange [xmin : xmax]

# determine y-range so that ±omega(x) fit in the current xrange
# (simple scan, avoids needing analytic max)
ymax = 0
nscan = 2000
xlo = xmin
xhi = xmax
xw  = xhi - xlo
if (xw == 0) { xw = 1 }

do for [k=0:nscan] {
    xx = xlo + xw*k/nscan
    yy = omega(xx)
    if (yy > ymax) { ymax = yy }
}
set yrange [-1.05*ymax : 1.05*ymax]

# total records (use the same column as xrange stats; avoid xrange filtering issues)
N = X_records

# styles
set style line 1 lw 2 lc rgb "#000000"   # boundaries
set style line 2 lw 1 lc rgb "#1f77b4"   # trajectory
set style line 3 lw 2 lc rgb "#d62728"   # rod

# animation loop
# - boundaries: y=±omega(x)
# - trajectory: data up to i-th row
# - rod: vector centered at (x,y) with angle theta

do for [i=1:N:frame_stride] {
    set title sprintf("trajectory (row=%d/%d)", i, N)

    plot \
        omega(x)   w l ls 1, \
        -omega(x)  w l ls 1, \
        file u 2:3 every ::1::i w l ls 2, \
        file u ($2-0.5*L*cos($4)):($3-0.5*L*sin($4)):(L*cos($4)):(L*sin($4)) every ::i::i w vectors nohead ls 3
}

set output
