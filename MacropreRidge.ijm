run("8-bit");
run("Convolve...", "text1=[-1 -1 -1 -1 -1\n-1 5 5 5 -1\n-1 5 5 5 -1\n-1 5 5 5 -1\n-1 -1 -1 -1 -1\n] normalize stack");
run("Subtract Background...", "rolling=10 stack");
run("Ridge Detection");
