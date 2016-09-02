binwidth = 5
set boxwidth binwidth
bin(x, width) = width*floor(x/width) + binwidth/2.0
plot 'Amplitude_488.txt' u (bin($1,binwidth)):(1.0) smooth freq with boxes

