
gammas_test=../data/dl3/predictions_rta_gammas_test_sample.csv

dir = build/

resolution = $(dir)/angular_resolution.pdf
resolution_vs_energy = $(dir)/angular_resolution_vs_energy.pdf

all: $(dir) $(resolution) $(resolution_vs_energy)

$(dir):
	mkdir -p  $(dir)

$(resolution): plot_angular_resolution.py $(dir)
	python plot_angular_resolution.py $(gammas_test) $(resolution)

$(resolution_vs_energy): plot_angular_resolution_vs_energy.py $(dir)
	python plot_angular_resolution_vs_energy.py $(gammas_test) $(resolution_vs_energy)
