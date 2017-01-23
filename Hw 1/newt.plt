set terminal svg enhanced size 1000 1000 fname "Times" fsize 36
set output "newt.svg"
set title "Sequence Pn of Newton Method"
set xlabel "i"
set ylabel "P_i"
plot "./graph_newt.dat" using 1:2 title ""