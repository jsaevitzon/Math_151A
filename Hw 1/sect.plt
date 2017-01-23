set terminal svg enhanced size 1000 1000 fname "Times" fsize 36
set output "sect.svg"
set title "Sequence Pn of Secant Method"
set xlabel "i"
set ylabel "P_i"
plot "./graph_secant.dat" using 1:2 title "" with points