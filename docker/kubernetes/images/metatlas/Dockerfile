FROM mambaorg/micromamba:0.12.0

RUN micromamba install -y -n base -c conda-forge \
	dill \
        gspread \
	hdf5 \
	ipywidgets \
	jupyterlab=3.0.14 \
	matplotlib \
	oauth2client \
	pandas \
	papermill \
	pymysql \
	pytables \
	pyyaml \
	rdkit \
	scikit-learn \
	sqlalchemy \
	tabulate \
	xlsxwriter && \
	rm /opt/conda/pkgs/cache/*

RUN pip install dataset  # not availble in conda

ENV METATLAS_LOCAL=True

EXPOSE 8888

RUN mkdir -p /io /src /work

WORKDIR /work

COPY meta_atlas.sqlite3 /work/root_workspace.db

RUN mkdir -p /global/project/projectdirs/metatlas/projects/spectral_libraries

COPY msms_refs_v3.tab /global/project/projectdirs/metatlas/projects/spectral_libraries/

RUN mkdir -p /project/projectdirs/metatlas/raw_data/akuftin/20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583/

COPY *.h5 /project/projectdirs/metatlas/raw_data/akuftin/20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583/

CMD ["/opt/conda/bin/jupyter", "nbclassic", "--ip=0.0.0.0", "--allow-root", "--ServerApp.token=''", "--ServerApp.root_dir=/"]
