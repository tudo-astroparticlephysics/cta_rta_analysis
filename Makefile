SHELL=/bin/bash
build_dir = build
data_dir = data

plot_overview = $(build_dir)/separator_performance.pdf
plot_overview = $(build_dir)/separator_performance.pdf
plot_overview_diffuse = $(build_dir)/separator_performance_diffuse.pdf
plot_overview_regressor = $(build_dir)/regressor_performance.pdf

plot_effective_area_pointlike = $(build_dir)/effective_area_pointlike.pdf
plot_effective_area = $(build_dir)/effective_area.pdf

plot_angular_resolution = $(build_dir)/angular_resolution.pdf
plot_angular_resolution_pointlike = $(build_dir)/angular_resolution_pointlike.pdf

plot_energy_resolution = $(build_dir)/energy_resolution.pdf
plot_energy_resolution_pointlike = $(build_dir)/energy_resolution_pointlike.pdf

plot_theta_square = $(build_dir)/theta_square.pdf

plot_auc_per_type = $(build_dir)/auc_per_type.pdf
plot_auc = $(build_dir)/auc.pdf

predicted_gammas = $(build_dir)/gamma_dl2.h5
predicted_protons = $(build_dir)/proton_dl2.h5

proton_test = $(build_dir)/proton_test.h5
proton_test_small = $(build_dir)/proton_test_small.h5
gamma_test = $(build_dir)/gamma_test.h5

electron_test = $(build_dir)/electron_test.h5

proton_train = $(build_dir)/proton_train.h5
gamma_train = $(build_dir)/gamma_train.h5

gamma_pointlike = $(build_dir)/gamma_pointlike.h5

predictions_separator = $(build_dir)/aict_predictions_separation.h5
model_separator =  $(build_dir)/separator.pkl

predictions_regressor = $(build_dir)/aict_predictions_regression.h5
model_regressor = $(build_dir)/regressor.pkl

sensitivity = $(build_dir)/sensitivity.pdf 
sensitivity_extrapolate = $(build_dir)/sensitivity_extrapolate.pdf 
sensitivity_exact = $(build_dir)/sensitivity_exact.pdf 
sensitivity_simple = $(build_dir)/sensitivity_simple.pdf 

config = configs/iact_config.yaml

PLOTS := $(plot_angular_resolution) $(plot_angular_resolution_pointlike) $(plot_overview) $(plot_overview_regressor)
PLOTS += $(sensitivity) $(sensitivity_extrapolate) $(sensitivity_exact) $(sensitivity_simple)
PLOTS += $(plot_effective_area_pointlike) $(plot_effective_area) 
PLOTS += $(plot_energy_resolution_pointlike) $(plot_energy_resolution) $(plot_theta_square) $(plot_auc) $(plot_auc_per_type)

all:  $(build_dir)/APPLICATION_DONE_BACKGROUND $(build_dir)/APPLICATION_DONE_SIGNAL $(PLOTS)

sensitivity: $(plot_sensitivity_all)

sensitivity_all: $(plot_sensitivity_all)  $(plot_sensitivity_sst) $(plot_sensitivity_lst) $(plot_sensitivity_mst)

clean:
	rm -rf $(build_dir)

$(build_dir):
	mkdir -p $(build_dir)

$(proton_test) $(proton_test_small) $(proton_train): $(data_dir)/protons.h5 | $(build_dir)
	# aict_apply_cuts $(config) $(data_dir)/dl2/protons.hdf5 $(build_dir)/protons_cutted.hdf5 -k telescope_events
	aict_split_data $(data_dir)/protons.h5 $(build_dir)/proton -n train -f 0.005 -n test_small -f 0.005 -n test  -f 0.99  -t cta
$(gamma_test) $(gamma_train): $(data_dir)/gammas_diffuse.h5 | $(build_dir)
	# aict_apply_cuts $(config) $(data_dir)/dl2/gammas.hdf5 $(build_dir)/gammas_cutted.hdf5 -k telescope_events
	aict_split_data $(data_dir)/gammas_diffuse.h5 $(build_dir)/gamma -n train -f 0.2 -n test -f 0.8  -t cta
$(electron_test): | $(build_dir)
	# aict_apply_cuts $(config) $(data_dir)/dl2/gammas.hdf5 $(build_dir)/gammas_cutted.hdf5 -k telescope_events
	cp $(data_dir)/electrons.h5 $(electron_test)
$(gamma_pointlike): | $(build_dir)
	# aict_apply_cuts $(config) $(data_dir)/dl2/gammas.hdf5 $(build_dir)/gammas_cutted.hdf5 -k telescope_events
	cp $(data_dir)/gammas.h5 $(gamma_pointlike)


$(model_separator) $(predictions_separator): $(proton_train) $(gamma_train) $(config)
	aict_train_separation_model $(config) $(gamma_train) $(proton_train) $(predictions_separator) $(model_separator)

$(model_regressor) $(predictions_regressor): $(gamma_train) $(config)
	aict_train_energy_regressor $(config) $(gamma_train) $(predictions_regressor) $(model_regressor)
#

$(build_dir)/APPLICATION_DONE_BACKGROUND: $(model_separator) $(config) $(proton_test) $(model_regressor) $(electron_test)
	aict_apply_separation_model $(config) $(proton_test) $(model_separator) --yes --chunksize 1000000
	aict_apply_separation_model $(config) $(proton_test_small) $(model_separator) --yes --chunksize 100000
	aict_apply_separation_model $(config) $(electron_test) $(model_separator) --yes --chunksize 100000
	aict_apply_energy_regressor $(config) $(proton_test) $(model_regressor) --yes --chunksize 1000000
	aict_apply_energy_regressor $(config) $(proton_test_small) $(model_regressor) --yes --chunksize 100000
	aict_apply_energy_regressor $(config) $(electron_test) $(model_regressor) --yes --chunksize 100000
	touch $(build_dir)/APPLICATION_DONE_BACKGROUND

$(build_dir)/APPLICATION_DONE_SIGNAL: $(model_separator) $(config) $(gamma_test) $(model_regressor) $(gamma_pointlike)
	aict_apply_separation_model $(config) $(gamma_test) $(model_separator) --yes --chunksize 400000
	aict_apply_separation_model $(config) $(gamma_pointlike) $(model_separator) --yes --chunksize 400000
	aict_apply_energy_regressor $(config) $(gamma_test) $(model_regressor) --yes --chunksize 400000
	aict_apply_energy_regressor $(config) $(gamma_pointlike) $(model_regressor) --yes --chunksize 400000
	touch $(build_dir)/APPLICATION_DONE_SIGNAL

# plot a bunch of performance values for the classifier
$(plot_overview): $(model_separator) $(predictions_separator) matplotlibrc $(config)
	aict_plot_separator_performance $(config) $(predictions_separator) $(model_separator) -o $(plot_overview)
$(plot_overview_regressor): $(model_regressor) $(predictions_regressor) matplotlibrc $(config)
	aict_plot_regressor_performance $(config) $(predictions_regressor) $(model_regressor) -o $(plot_overview_regressor)

$(plot_effective_area_pointlike): $(gamma_pointlike) $(build_dir)/APPLICATION_DONE_SIGNAL 
	cta_plot_effective_area $(gamma_pointlike) -o $(plot_effective_area_pointlike) -t 0 -t 0.7
$(plot_effective_area): $(gamma_test)	$(build_dir)/APPLICATION_DONE_SIGNAL 
	cta_plot_effective_area $(gamma_test) -o $(plot_effective_area) -t 0 -t 0.7

$(plot_angular_resolution): $(gamma_test) $(build_dir)/APPLICATION_DONE_SIGNAL
	cta_plot_angular_resolution $(gamma_test) -o $(plot_angular_resolution) -m 4 --reference --plot_e_reco
$(plot_angular_resolution_pointlike): $(gamma_pointlike) $(build_dir)/APPLICATION_DONE_SIGNAL
	cta_plot_angular_resolution $(gamma_pointlike) -o $(plot_angular_resolution_pointlike) -m 4 --reference --plot_e_reco

$(plot_energy_resolution_pointlike): $(gamma_pointlike) $(build_dir)/APPLICATION_DONE_SIGNAL
	cta_plot_energy_resolution $(gamma_pointlike) -o $(plot_energy_resolution_pointlike) --reference
$(plot_energy_resolution): $(gamma_test) $(build_dir)/APPLICATION_DONE_SIGNAL
	cta_plot_energy_resolution $(gamma_test) -o $(plot_energy_resolution) --reference

$(plot_theta_square): $(gamma_pointlike) $(electron_test) $(proton_test) $(build_dir)/APPLICATION_DONE_SIGNAL $(build_dir)/APPLICATION_DONE_BACKGROUND
	cta_plot_theta_square $(gamma_pointlike) $(proton_test) $(electron_test) -j 1 -o $(plot_theta_square)

$(plot_auc_per_type): $(gamma_test) $(proton_test_small) $(build_dir)/APPLICATION_DONE_SIGNAL $(build_dir)/APPLICATION_DONE_BACKGROUND
	cta_plot_auc_per_type $(gamma_test) $(proton_test_small) -o $(plot_auc_per_type)
$(plot_auc): $(gamma_test) $(proton_test_small) $(build_dir)/APPLICATION_DONE_SIGNAL $(build_dir)/APPLICATION_DONE_BACKGROUND
	cta_plot_auc $(gamma_test) $(proton_test_small) -o $(plot_auc)

$(sensitivity): $(gamma_pointlike) $(proton_test) $(electron_test) $(build_dir)/APPLICATION_DONE_SIGNAL $(build_dir)/APPLICATION_DONE_BACKGROUND | $(build_dir)
	cta_plot_sensitivity $(gamma_pointlike) $(proton_test) $(electron_test) -o $(sensitivity) --reference --requirement -m 4 --method histogram
$(sensitivity_extrapolate): $(gamma_pointlike) $(proton_test) $(electron_test) $(build_dir)/APPLICATION_DONE_SIGNAL $(build_dir)/APPLICATION_DONE_BACKGROUND | $(build_dir)
	cta_plot_sensitivity $(gamma_pointlike) $(proton_test) $(electron_test) -o $(sensitivity_extrapolate) --reference --requirement -m 4 
$(sensitivity_simple): $(gamma_pointlike) $(proton_test) $(electron_test) $(build_dir)/APPLICATION_DONE_SIGNAL $(build_dir)/APPLICATION_DONE_BACKGROUND | $(build_dir)
	cta_plot_sensitivity $(gamma_pointlike) $(proton_test) $(electron_test) -o $(sensitivity_simple) --reference --requirement -m 4 --method simple
$(sensitivity_exact): $(gamma_pointlike) $(proton_test) $(electron_test) $(build_dir)/APPLICATION_DONE_SIGNAL $(build_dir)/APPLICATION_DONE_BACKGROUND | $(build_dir)
	cta_plot_sensitivity $(gamma_pointlike) $(proton_test) $(electron_test) -o $(sensitivity_exact) --reference --requirement -m 4 --method exact

# $(sensitivity_sst_fits): $(gamma_test_sst) $(proton_test_sst)
# 	python effective_area/calculate_sensitivity.py $(gamma_test_sst) $(proton_test_sst) $(sensitivity_sst_fits) -n 20
# $(plot_sensitivity_sst): $(sensitivity_sst_fits)
# 	python effective_area/plot_sensitivity.py   $(sensitivity_sst_fits)  -o $(plot_sensitivity_sst)


# $(sensitivity_mst_fits): $(gamma_test_mst) $(proton_test_mst)
# 	python effective_area/calculate_sensitivity.py $(gamma_test_mst) $(proton_test_mst) $(sensitivity_mst_fits) -n 20
# $(plot_sensitivity_mst): $(sensitivity_mst_fits)
# 	python effective_area/plot_sensitivity.py   $(sensitivity_mst_fits)  -o $(plot_sensitivity_mst)

# $(sensitivity_lst_fits): $(gamma_test_lst) $(proton_test_lst)
# 	python effective_area/calculate_sensitivity.py $(gamma_test_lst) $(proton_test_lst) $(sensitivity_lst_fits) -n 20
# $(plot_sensitivity_lst): $(sensitivity_lst_fits)
# 	python effective_area/plot_sensitivity.py   $(sensitivity_lst_fits)  -o $(plot_sensitivity_lst)


# $(plot_sensitivity_combined): $(sensitivity_all_fits) $(sensitivity_lst_fits) $(sensitivity_mst_fits) $(sensitivity_sst_fits)
# 	python effective_area/plot_sensitivity.py   $(sensitivity_sst_fits) $(sensitivity_mst_fits) $(sensitivity_lst_fits) $(sensitivity_all_fits) -l 'SST' -l 'MST' -l 'LST' -l 'all' -c 'lightslategray' -c 'silver' -c 'gray' -c 'crimson' -o $(plot_sensitivity_combined)

# $(plot_sensitivity_sst): $(gamma_test_sst) $(proton_test_sst)
# 	python effective_area/plot_sensitivity.py  -g $(gamma_test_sst) -p $(proton_test_sst) -n 20 -o $(plot_sensitivity_sst)

# $(plot_sensitivity_mst): $(gamma_test_mst) $(proton_test_mst)
# 	python effective_area/plot_sensitivity.py  -g $(gamma_test_mst) -p $(proton_test_mst) -n 20 -o $(plot_sensitivity_mst)

# $(plot_sensitivity_lst): $(gamma_test_lst) $(proton_test_lst)
# 	python effective_area/plot_sensitivity.py  -g $(gamma_test_lst) -p $(proton_test_lst) -n 20 -o $(plot_sensitivity_lst)

# $(plot_sensitivity_combined): $(gamma_test_lst) $(proton_test_lst) $(gamma_test_mst) $(proton_test_mst) $(gamma_test_sst) $(proton_test_sst)
# 	python effective_area/plot_sensitivity.py  -g $(gamma_test_sst) -g $(gamma_test_mst) -g $(gamma_test_lst) -p $(proton_test_sst)  -p $(proton_test_mst) -p $(proton_test_lst) -l 'SST' -l 'MST' -l 'LST' -n 20 -o $(plot_sensitivity_combined)
