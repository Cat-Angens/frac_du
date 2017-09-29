set terminal png size 1280,720 linewidth 2
set output "result.png"


set multiplot layout 2,2 title "Решение при постоянной равномерной скорости"
set tmargin 2
set yrange [-0.1:1.1]

set title "Начальное распределение"
plot "tvd\\field_0001.txt" w lp pt 6 title "tvd-схема", "notvd\\field_0001.txt" w lp pt 6 title "схема 2-го порядка"

unset title
plot "tvd\\field_0034.txt" w lp pt 6 title "tvd-схема", "notvd\\field_0034.txt" w lp pt 6 title "схема 2-го порядка"

plot "tvd\\field_0067.txt" w lp pt 6 title "tvd-схема", "notvd\\field_0067.txt" w lp pt 6 title "схема 2-го порядка"

plot "tvd\\field_0100.txt" w lp pt 6 title "tvd-схема", "notvd\\field_0100.txt" w lp pt 6 title "схема 2-го порядка"