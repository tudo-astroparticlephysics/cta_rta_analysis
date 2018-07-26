SHELL=/bin/bash
build_dir = build
data_dir = data

plot_overview = $(build_dir)/separator_performance.pdf
plot_overview_diffuse = $(build_dir)/separator_performance_diffuse.pdf
plot_overview_regressor = $(build_dir)/regressor_performance.pdf

sensitivity_all_fits = $(build_dir)/sensitivity.fits
plot_sensitivity_all = $(build_dir)/sensitivity.pdf

sensitivity_sst_fits = $(build_dir)/sensitivity_sst.fits
plot_sensitivity_sst = $(build_dir)/sensitivity_sst.pdf

sensitivity_mst_fits = $(build_dir)/sensitivity_mst.fits
plot_sensitivity_mst = $(build_dir)/sensitivity_mst.pdf

sensitivity_lst_fits = $(build_dir)/sensitivity_lst.fits
plot_sensitivity_lst = $(build_dir)/sensitivity_lst.pdf

plot_sensitivity_combined = $(build_dir)/sensitivity_combined.pdf

gamma_output = $(data_dir)/dl2/gammas.hdf5
gamma_output_diffuse = $(data_dir)/dl2/gammas_diffuse.hdf5
proton_output = $(data_dir)/dl2/protons.hdf5

proton_test = $(build_dir)/protons_test.hdf5
proton_diffuse_test = $(build_dir)/protons_diffuse_test.hdf5
gamma_test = $(build_dir)/gammas_test.hdf5
gamma_diffuse_test = $(build_dir)/gammas_diffuse_test.hdf5

electron_test = $(build_dir)/electrons_test.hdf5

proton_train = $(build_dir)/protons_train.hdf5
gamma_train = $(build_dir)/gammas_train.hdf5
gamma_diffuse_train = $(build_dir)/gammas_diffuse_train.hdf5

predictions_separator = $(build_dir)/aict_predictions_separation.hdf5
predictions_separator_diffuse = $(build_dir)/aict_predictions_separation_diffuse.hdf5
model_separator =  $(build_dir)/separator.pkl
model_separator_diffuse =  $(build_dir)/separator_diffuse.pkl

predictions_regressor = $(build_dir)/aict_predictions_regression.hdf5
model_regressor = $(build_dir)/regressor.pkl
model_regressor_diffuse =  $(build_dir)/regressor_diffuse.pkl

config = configs/iact_config.yaml

gamma_test_sst = $(build_dir)/gammas_test_sst.hdf5
proton_test_sst = $(build_dir)/protons_test_sst.hdf5

gamma_test_mst = $(build_dir)/gammas_test_mst.hdf5
proton_test_mst = $(build_dir)/protons_test_mst.hdf5

gamma_test_lst = $(build_dir)/gammas_test_lst.hdf5
proton_test_lst = $(build_dir)/protons_test_lst.hdf5


all: $(plot_overview) $(plot_overview_regressor) $(plot_overview_diffuse) $(build_dir)/APPLICATION_DONE

sensitivity: $(plot_sensitivity_all)

sensitivity_all: $(plot_sensitivity_all)  $(plot_sensitivity_sst) $(plot_sensitivity_lst) $(plot_sensitivity_mst)

clean:
	rm -rf $(build_dir)

$(build_dir):
	mkdir -p $(build_dir)


$(gamma_test): $(config) | $(build_dir)


$(proton_test) $(proton_train): $(data_dir)/dl2/protons.hdf5 | $(build_dir)
	# aict_apply_cuts $(config) $(data_dir)/dl2/protons.hdf5 $(build_dir)/protons_cutted.hdf5 -k telescope_events
	aict_split_data $(data_dir)/dl2/protons.hdf5 $(build_dir)/protons -n test -f 0.98 -n train -f 0.02  -t cta

$(proton_diffuse_test): $(proton_test)| $(build_dir)
	cp $(proton_test) $(proton_diffuse_test)


$(gamma_test) $(gamma_train): $(data_dir)/dl2/gammas.hdf5 | $(build_dir)
	# aict_apply_cuts $(config) $(data_dir)/dl2/gammas.hdf5 $(build_dir)/gammas_cutted.hdf5 -k telescope_events
	aict_split_data $(data_dir)/dl2/gammas.hdf5 $(build_dir)/gammas -n test -f 0.97 -n train -f 0.02  -t cta
#

$(electron_test): | $(build_dir)
	# aict_apply_cuts $(config) $(data_dir)/dl2/gammas.hdf5 $(build_dir)/gammas_cutted.hdf5 -k telescope_events
	cp $(data_dir)/dl2/electrons.hdf5 $(electron_test)
#
$(gamma_diffuse_test) $(gamma_diffuse_train): $(data_dir)/dl2/gammas_diffuse.hdf5 | $(build_dir)
	aict_split_data $(gamma_output_diffuse) $(build_dir)/gammas_diffuse -n test -f 0.85 -n train -f 0.15  -t cta
#
# $(proton_diffuse_test): $(proton_test)
# 	cp $(proton_test) $(proton_diffuse_test)

$(model_separator) $(predictions_separator): $(proton_train) $(gamma_train) $(config)
	aict_train_separation_model $(config) $(gamma_train) $(proton_train) $(predictions_separator) $(model_separator)

$(model_regressor) $(predictions_regressor): $(gamma_train) $(config)
	aict_train_energy_regressor $(config) $(gamma_train) $(predictions_regressor) $(model_regressor)
#
$(model_separator_diffuse) $(predictions_separator_diffuse): $(proton_train) $(gamma_diffuse_train) $(config)
	aict_train_separation_model $(config) $(gamma_diffuse_train) $(proton_train) $(predictions_separator_diffuse) $(model_separator_diffuse)

# $(model_regressor_diffuse) $(predictions): $(gamma_diffuse_train) $(config)
# 	aict_train_energy_regressor $(config) $(gamma_diffuse_train) $(predictions_regressor) $(model_regressor_diffuse)


$(build_dir)/APPLICATION_DONE: $(model_separator) $(config) $(gamma_test) $(proton_test) $(model_regressor) $(electron_test)
	aict_apply_separation_model $(config) $(gamma_test) $(model_separator) --yes --chunksize 400000
	aict_apply_separation_model $(config) $(proton_test) $(model_separator) --yes --chunksize 400000
	aict_apply_separation_model $(config) $(electron_test) $(model_separator) --yes --chunksize 400000
	aict_apply_energy_regressor $(config) $(gamma_test) $(model_regressor) --yes --chunksize 400000
	aict_apply_energy_regressor $(config) $(proton_test) $(model_regressor) --yes --chunksize 400000
	aict_apply_energy_regressor $(config) $(electron_test) $(model_regressor) --yes --chunksize 400000
	touch $(build_dir)/APPLICATION_DONE

#
$(build_dir)/APPLICATION_DIFFUSE_DONE: $(model_separator_diffuse) $(config) $(proton_diffuse_test) $(gamma_diffuse_test)
	aict_apply_separation_model $(config) $(gamma_diffuse_test) $(model_separator_diffuse) --yes --chunksize 100000
	aict_apply_separation_model $(config) $(proton_diffuse_test) $(model_separator_diffuse) --yes --chunksize 400000
# 	aict_apply_energy_regressor $(config) $(gamma_diffuse_test) $(model_regressor_diffuse) --yes --chunksize 100000
	touch $(build_dir)/APPLICATION_DIFFUSE_DONE


# plot a bunch of performance values for the classifier
$(plot_overview): $(model_separator) $(predictions_separator) matplotlibrc $(config)
	aict_plot_separator_performance $(config) $(predictions_separator) $(model_separator) -o $(plot_overview)


$(plot_overview_diffuse): $(model_separator_diffuse) $(predictions_separator_dffuse) matplotlibrc $(config)
	aict_plot_separator_performance $(config) $(predictions_separator) $(model_separator_diffuse) -o $(plot_overview_diffuse)


$(plot_overview_regressor): $(model_regressor) $(predictions_regressor) matplotlibrc $(config)
	aict_plot_regressor_performance $(config) $(predictions_regressor) $(model_regressor) -o $(plot_overview_regressor)


# $(plot_effective_area): $(gamma_test) effective_area/plot_effective_area.py
# 	python effective_area/plot_effective_area.py -g $(gamma_test) -o $(plot_effective_area)


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


$(sensitivity_all_fits): $(gamma_test) $(proton_test)
	python effective_area/calculate_sensitivity.py $(gamma_test) $(proton_test) $(sensitivity_all_fits) -n 20
$(plot_sensitivity_all): $(sensitivity_all_fits)
	python effective_area/plot_sensitivity.py   $(sensitivity_all_fits)  -o $(plot_sensitivity_all)

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
