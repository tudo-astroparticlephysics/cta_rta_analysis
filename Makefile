SHELL=/bin/bash
build_dir = build
data_dir = data

plot_overview = $(build_dir)/separator_performance.pdf
plot_overview_regressor = $(build_dir)/regressor_performance.pdf
plot_roc = $(build_dir)/roc.pdf
plot_roc_per_telescope = $(build_dir)/roc_per_telescope.pdf
plot_auc_vs_energy = $(build_dir)/auc_vs_energy.pdf
plot_hists = $(build_dir)/hists.pdf

plot_angular_resolution = $(build_dir)/angular_resolution_1d.pdf
plot_ang_res_energy = $(build_dir)/angular_resolution_energy.png

plot_effective_area = $(build_dir)/effective_area_all.pdf

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
proton_output = $(data_dir)/dl2/protons.hdf5

proton_test = $(build_dir)/protons_test.hdf5
gamma_test = $(build_dir)/gammas_test.hdf5

proton_train = $(build_dir)/protons_train.hdf5
gamma_train = $(build_dir)/gammas_train.hdf5

predictions_separator = $(build_dir)/klaas_predictions_separation.hdf5
model_separator =  $(build_dir)/separator.pkl
config_separator = configs/separator.yaml

predictions_regressor = $(build_dir)/klaas_predictions_regression.hdf5
model_regressor = $(build_dir)/regressor.pkl
config_regressor = configs/regressor.yaml

gamma_test_sst = $(build_dir)/gammas_test_sst.hdf5
proton_test_sst = $(build_dir)/protons_test_sst.hdf5

gamma_test_mst = $(build_dir)/gammas_test_mst.hdf5
proton_test_mst = $(build_dir)/protons_test_mst.hdf5

gamma_test_lst = $(build_dir)/gammas_test_lst.hdf5
proton_test_lst = $(build_dir)/protons_test_lst.hdf5


all: $(plot_overview) $(plot_overview_regressor) $(plot_roc) $(plot_hists) $(plot_roc_per_telescope) $(plot_auc_vs_energy) $(plot_angular_resolution) $(plot_ang_res_energy) $(plot_effective_area) $(plot_sensitivity_all) $(plot_sensitivity_sst) $(plot_sensitivity_lst) $(plot_sensitivity_mst)

clean:
	rm -rf $(build_dir)

$(build_dir):
	mkdir -p $(build_dir)


$(proton_test) $(proton_train): $(proton_output) | $(build_dir)
	klaas_split_data $(proton_output) $(build_dir)/protons -n test -f 0.92 -n train -f 0.08  -t cta

$(gamma_test) $(gamma_train): $(gamma_output) | $(build_dir)
	klaas_split_data $(gamma_output) $(build_dir)/gammas -n test -f 0.92 -n train -f 0.08  -t cta

$(model_separator) $(predictions): $(proton_train) $(gamma_train) $(config_separator)
	klaas_train_separation_model $(config_separator) $(gamma_train) $(proton_train) $(predictions_separator) $(model_separator)

$(model_regressor) $(predictions): $(gamma_train) $(config_regressor)
	klaas_train_energy_regressor $(config_regressor) $(gamma_train) $(predictions_regressor) $(model_regressor)


$(build_dir)/APPLICATION_DONE: $(proton_train) $(gamma_train) $(model_separator) $(config_separator) $(gamma_test) $(proton_test) $(model_regressor)
	klaas_apply_separation_model $(config_separator) $(gamma_test) $(model_separator) --yes --chunksize 100000
	klaas_apply_separation_model $(config_separator) $(proton_test) $(model_separator) --yes --chunksize 100000
	klaas_apply_energy_regressor $(config_regressor) $(gamma_test) $(model_regressor) --yes --chunksize 100000
	touch $(build_dir)/APPLICATION_DONE



# plot a bunch of performance values for the classifier
$(plot_overview): $(model_separator) $(predictions_separator) matplotlibrc configs/separator.yaml
	klaas_plot_separator_performance $(config_separator) $(predictions_separator) $(model_separator) -o $(plot_overview)

$(plot_overview_regressor): $(model_regressor) $(predictions_regressor) matplotlibrc configs/regressor.yaml
	klaas_plot_regressor_performance $(config_regressor) $(predictions_regressor) $(model_regressor) -o $(plot_overview_regressor)


# plot a roc curve
$(plot_roc): $(proton_test) $(gamma_test) matplotlibrc ml/plot_multi_tel_auc.py $(build_dir)/APPLICATION_DONE
	python ml/plot_multi_tel_auc.py $(gamma_test) $(proton_test) -o $(plot_roc)
# plot a roc curve per telescope type
$(plot_roc_per_telescope): $(proton_test) $(gamma_test) matplotlibrc ml/plot_auc_per_type.py $(build_dir)/APPLICATION_DONE
	python ml/plot_auc_per_type.py $(gamma_test) $(proton_test) -o $(plot_roc_per_telescope)
# plot a auc vs energy
$(plot_auc_vs_energy): $(proton_test) $(gamma_test) matplotlibrc ml/plot_auc_vs_energy.py $(build_dir)/APPLICATION_DONE
	python ml/plot_auc_vs_energy.py $(gamma_test) $(proton_test) -o $(plot_auc_vs_energy)

# plot a prediction hists
$(plot_hists):$(gamma_test) matplotlibrc ml/plot_prediction_hists.py $(build_dir)/APPLICATION_DONE
	python ml/plot_prediction_hists.py $(gamma_test) $(proton_test) -o $(plot_hists)


$(gamma_test_sst) : matplotlibrc processing/reconstruct_direction.py $(gamma_test) $(build_dir)/APPLICATION_DONE
	python processing/reconstruct_direction.py $(gamma_test) $(gamma_test_sst) ./processing/instrument_description.pkl -t SST -y
$(proton_test_sst) : matplotlibrc processing/reconstruct_direction.py $(proton_test) $(build_dir)/APPLICATION_DONE
	python processing/reconstruct_direction.py $(proton_test) $(proton_test_sst) ./processing/instrument_description.pkl -t SST -y

$(gamma_test_mst) : matplotlibrc processing/reconstruct_direction.py $(gamma_test) $(build_dir)/APPLICATION_DONE
	python processing/reconstruct_direction.py $(gamma_test) $(gamma_test_mst) ./processing/instrument_description.pkl -t MST -y
$(proton_test_mst) : matplotlibrc processing/reconstruct_direction.py $(proton_test) $(build_dir)/APPLICATION_DONE
	python processing/reconstruct_direction.py $(proton_test) $(proton_test_mst) ./processing/instrument_description.pkl -t MST -y

$(gamma_test_lst) : matplotlibrc processing/reconstruct_direction.py $(gamma_test) $(build_dir)/APPLICATION_DONE
	python processing/reconstruct_direction.py $(gamma_test) $(gamma_test_lst) ./processing/instrument_description.pkl -t LST -y
$(proton_test_lst) : matplotlibrc processing/reconstruct_direction.py $(proton_test) $(build_dir)/APPLICATION_DONE
	python processing/reconstruct_direction.py $(proton_test) $(proton_test_lst) ./processing/instrument_description.pkl -t LST -y


$(plot_ang_res_energy): $(gamma_test) angular_resolution/plot_angular_resolution_vs_energy.py
	python angular_resolution/plot_angular_resolution_vs_energy.py $(gamma_test) -o $(plot_ang_res_energy)

$(plot_angular_resolution): $(gamma_test) angular_resolution/plot_angular_resolution_1d.py
	python angular_resolution/plot_angular_resolution_1d.py $(gamma_test) -o $(plot_angular_resolution)

$(plot_effective_area): $(gamma_test) effective_area/plot_effective_area.py
	python effective_area/plot_effective_area.py -g $(gamma_test) -o $(plot_effective_area)


$(sensitivity_sst_fits): $(gamma_test_sst) $(proton_test_sst)
	python effective_area/calculate_sensitivity.py $(gamma_test_sst) $(proton_test_sst) $(sensitivity_sst_fits) -n 20
$(plot_sensitivity_sst): $(sensitivity_sst_fits)
	python effective_area/plot_sensitivity.py   $(sensitivity_sst_fits)  -o $(plot_sensitivity_sst)

$(sensitivity_mst_fits): $(gamma_test_mst) $(proton_test_mst)
	python effective_area/calculate_sensitivity.py $(gamma_test_mst) $(proton_test_mst) $(sensitivity_mst_fits) -n 20
$(plot_sensitivity_mst): $(sensitivity_mst_fits)
	python effective_area/plot_sensitivity.py   $(sensitivity_mst_fits)  -o $(plot_sensitivity_mst)

$(sensitivity_lst_fits): $(gamma_test_lst) $(proton_test_lst)
	python effective_area/calculate_sensitivity.py $(gamma_test_lst) $(proton_test_lst) $(sensitivity_lst_fits) -n 20
$(plot_sensitivity_lst): $(sensitivity_lst_fits)
	python effective_area/plot_sensitivity.py   $(sensitivity_lst_fits)  -o $(plot_sensitivity_lst)


$(sensitivity_all_fits): $(gamma_test) $(proton_test)
	python effective_area/calculate_sensitivity.py $(gamma_test) $(proton_test) $(sensitivity_all_fits) -n 20
$(plot_sensitivity_all): $(sensitivity_all_fits)
	python effective_area/plot_sensitivity.py   $(sensitivity_all_fits)  -o $(plot_sensitivity_all)

# $(plot_sensitivity_sst): $(gamma_test_sst) $(proton_test_sst)
# 	python effective_area/plot_sensitivity.py  -g $(gamma_test_sst) -p $(proton_test_sst) -n 20 -o $(plot_sensitivity_sst)
#
# $(plot_sensitivity_mst): $(gamma_test_mst) $(proton_test_mst)
# 	python effective_area/plot_sensitivity.py  -g $(gamma_test_mst) -p $(proton_test_mst) -n 20 -o $(plot_sensitivity_mst)
#
# $(plot_sensitivity_lst): $(gamma_test_lst) $(proton_test_lst)
# 	python effective_area/plot_sensitivity.py  -g $(gamma_test_lst) -p $(proton_test_lst) -n 20 -o $(plot_sensitivity_lst)
#
# $(plot_sensitivity_combined): $(gamma_test_lst) $(proton_test_lst) $(gamma_test_mst) $(proton_test_mst) $(gamma_test_sst) $(proton_test_sst)
# 	python effective_area/plot_sensitivity.py  -g $(gamma_test_sst) -g $(gamma_test_mst) -g $(gamma_test_lst) -p $(proton_test_sst)  -p $(proton_test_mst) -p $(proton_test_lst) -l 'SST' -l 'MST' -l 'LST' -n 20 -o $(plot_sensitivity_combined)
