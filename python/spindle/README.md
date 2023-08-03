# Spindle analysis:


# Plots

    get_parameters.py run* > parameters.txt
    plot_fiber_length.py run* > fiber_length.txt
    plot_spindle_length.py run* > spindle_length.txt


    plot_spindles.py




# Sort runs (optional)

    get_parameters run* > x
    sort -k 2 -n x > z
    mkdir tmp
    mv r* tmp
    cut -c 1-7 z | xargs -n 1 -I {} collect.py run%%%% tmp/{}
    rmdir tmp
    rm x y z

# Underlying reports

	report time spindle:length > spindle_length.txt
	report3 fiber:mark > fiber_length.txt
