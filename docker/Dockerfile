FROM python:3.8-slim-buster

ARG BASE_DATA_URL=https://portal.nersc.gov/cfs/m2650/metatlas/test_data/ci01
ARG REFS_DIR=/global/project/projectdirs/metatlas/projects/spectral_libraries
ARG H5_DIR=/project/projectdirs/metatlas/raw_data/akuftin/20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583

ENV METATLAS_LOCAL=True

EXPOSE 8888

RUN apt-get update && apt-get install -y jq && \
    rm -rf /var/lib/apt/lists/*

COPY requirements.txt /requirements.txt

RUN pip install --quiet -r requirements.txt

RUN mkdir -p /io /src /work $REFS_DIR $H5_DIR

ADD $BASE_DATA_URL/msms_refs_v3.tab $REFS_DIR/

ADD $BASE_DATA_URL/20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_49_Cone-S1_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run34.h5 $H5_DIR/
ADD $BASE_DATA_URL/20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_57_Cone-S2_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run40.h5 $H5_DIR/
ADD $BASE_DATA_URL/20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_65_Cone-S3_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run16.h5 $H5_DIR/
ADD $BASE_DATA_URL/20201106_JGI-AK_PS-KM_505892_OakGall_final_QE-HF_HILICZ_USHXG01583_POS_MSMS_73_Cone-S4_1_Rg70to1050-CE102040-QlobataAkingi-S1_Run31.h5 $H5_DIR/

ADD $BASE_DATA_URL/meta_atlas.sqlite3 /work/root_workspace.db

WORKDIR /work

RUN apt-get update && apt-get install -y libxrender1 && \
    rm -rf /var/lib/apt/lists/*

CMD ["/usr/local/bin/jupyter", "nbclassic", "--ip=0.0.0.0", "--allow-root", "--ServerApp.token=''", "--ServerApp.root_dir=/"]
