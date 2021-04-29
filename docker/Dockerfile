FROM mambaorg/micromamba:0.12.0

ARG BASE_DATA_URL=https://portal.nersc.gov/cfs/m2650/metatlas/test_data/ci01
ARG REFS_DIR=/global/project/projectdirs/metatlas/projects/spectral_libraries
ARG H5_DIR=/project/projectdirs/metatlas/raw_data/akuftin/20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583

ENV METATLAS_LOCAL=True

EXPOSE 8888

#RUN micromamba install -y -n base -c conda-forge \
#	dill=0.3.3 \
#        gspread=3.7.0 \
#	hdf5=1.10.6 \
#	ipywidgets=7.6.3 \
#	jupyterlab=3.0.14 \
#	matplotlib=3.4.1 \
#	oauth2client=4.1.3 \
#	pandas=1.2.4 \
#	papermill=2.3.3 \
#	pip=21.1 \
#	pymysql=1.0.2 \
#	pytables=3.6.1 \
#	pyyaml=5.4.1 \
#	rdkit=2021.03.1 \
#	scikit-learn=0.24.1 \
#	sqlalchemy=1.4.11 \
#	tabulate=0.8.9 \
#	xlsxwriter=1.4.0 && \
#	rm /opt/conda/pkgs/cache/*

COPY metatlas_env.yaml /metatlas_env.yaml

RUN micromamba install -y -n base -f /metatlas_env.yaml && \
    micromamba clean --all --yes

RUN mkdir -p /io /src /work $REFS_DIR $H5_DIR

ADD $BASE_DATA_URL/msms_refs_v3.tab $REFS_DIR/

ADD $BASE_DATA_URL/20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5 $H5_DIR/
ADD $BASE_DATA_URL/20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_57_Cone-S2_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run40.h5 $H5_DIR/
ADD $BASE_DATA_URL/20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_65_Cone-S3_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run16.h5 $H5_DIR/
ADD $BASE_DATA_URL/20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_73_Cone-S4_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run31.h5 $H5_DIR/

ADD $BASE_DATA_URL/meta_atlas.sqlite3 /work/root_workspace.db

WORKDIR /work

CMD ["/opt/conda/bin/jupyter", "nbclassic", "--ip=0.0.0.0", "--allow-root", "--ServerApp.token=''", "--ServerApp.root_dir=/"]
