SHELL=/bin/bash
build_dir = build
data_dir = data

plot_overview = $(build_dir)/separator_performance.pdf
plot_overview_regressor = $(build_dir)/regressor_performance.pdf
plot_roc = $(build_dir)/roc.pdf
plot_hists = $(build_dir)/hists.pdf
plot_angular_resolution = $(build_dir)/angular_resolution_1d.pdf

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

gamma_dl3 = $(build_dir)/gamma_dl3.hdf5

# all: $(proton_test) $(proton_train) $(gamma_test) $(gamma_train) ../build/ANGULAR ../build/COLLECTION ../build/ML_PERF ../build/SENSITIVITY

all: $(plot_overview) $(plot_overview_regressor) $(plot_roc) $(plot_hists) $(plot_angular_resolution)

clean:
	rm -rf $(build_dir)

$(build_dir):
	mkdir -p $(build_dir)


$(proton_test) $(proton_train): $(proton_output) | $(build_dir)
	klaas_split_data $(proton_output) $(build_dir)/protons -n test -f 0.7 -n train -f 0.3  -t cta

$(gamma_test) $(gamma_train): $(gamma_output) | $(build_dir)
	klaas_split_data $(gamma_output) $(build_dir)/gammas -n test -f 0.7 -n train -f 0.3  -t cta

$(model_separator) $(predictions): $(proton_train) $(gamma_train) $(config_separator)
	klaas_train_separation_model $(config_separator) $(gamma_train) $(proton_train) $(predictions_separator) $(model_separator)

$(model_regressor) $(predictions): $(gamma_train) $(config_regressor)
	klaas_train_energy_regressor $(config_regressor) $(gamma_train) $(predictions_regressor) $(model_regressor)


$(build_dir)/APPLICATION_DONE: $(proton_train) $(gamma_train) $(model_separator) $(config_separator)
	klaas_apply_separation_model $(config_separator) $(gamma_test) $(model_separator) --yes
	klaas_apply_separation_model $(config_separator) $(proton_test) $(model_separator) --yes
	klaas_apply_energy_regressor $(config_regressor) $(gamma_test) $(model_regressor) --yes
	touch $(build_dir)/APPLICATION_DONE



# plot a bunch of performance values for the classifier
$(plot_overview): $(model_separator) $(predictions_separator) matplotlibrc configs/separator.yaml
	klaas_plot_separator_performance $(config_separator) $(predictions_separator) $(model_separator) -o $(plot_overview)

$(plot_overview_regressor): $(model_regressor) $(predictions_regressor) matplotlibrc configs/regressor.yaml
	klaas_plot_regressor_performance $(config_regressor) $(predictions_regressor) $(model_regressor) -o $(plot_overview_regressor)


 # plot a roc curve
$(plot_roc): $(model_separator) matplotlibrc ml/plot_multi_tel_auc.py $(build_dir)/APPLICATION_DONE
	python ml/plot_multi_tel_auc.py $(gamma_test) $(proton_test) -o $(plot_roc)
	# plot a prediction hists
$(plot_hists): $(model_separator) matplotlibrc ml/plot_prediction_hists.py $(build_dir)/APPLICATION_DONE
	python ml/plot_prediction_hists.py $(gamma_test) $(proton_test) -o $(plot_hists)


	# reconstruct direction
$(gamma_dl3): $(model_regressor) matplotlibrc processing/reconstruct_direction.py $(build_dir)/APPLICATION_DONE
	python processing/reconstruct_direction.py $(gamma_test) ./processing/instrument_description.pkl $(gamma_dl3)
$(plot_angular_resolution): $(gamma_dl3)
	python angular_resolution/plot_angular_resolution_1d.py $(gamma_dl3) -o $(plot_angular_resolution)
