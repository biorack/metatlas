FROM continuumio/miniconda2

#
# pactolus requirements
#
RUN conda install -y -c anaconda pandas==0.24.2
RUN conda install -y -c anaconda configparser==3.7.4
RUN conda install -y -c bioconda pymzml==0.7.8
RUN conda install -y -c anaconda pytables==3.4.2
RUN conda install -y -c conda-forge matplotlib-venn==0.11.5
RUN conda install -y -c anaconda pyqt==5.9.2
RUN apt-get update && \
    apt install -y libgl1-mesa-glx libxrender-dev

RUN conda install -y -c anaconda h5py==2.7.0
RUN conda install -y -c rdkit rdkit==2018.09.3.0

ENV HDF5_USE_FILE_LOCKING FALSE


#
# spectral hits requirements
#
RUN conda install -y scipy==1.2.1
RUN apt-get update && \
	apt install -y libgl1-mesa-glx

#
# mzml convert requirements
#
RUN conda install -y -c bioconda pymzml=0.7.8
RUN conda install -y pytables==3.4.2

# 
# Copy all scripts and set env vars
# 
COPY metatlas /usr/local/metatlas
COPY pactolus /root/pactolus
ENV PATH /usr/local/metatlas:/root/pactolus:$PATH
ENV PYTHONPATH /usr/local:/usr/local/metatlas:/usr/local/metatlas/metatlas/helpers:/root

COPY spectral_hits_jaws.py /usr/local/bin/spectral_hits_jaws.py
COPY mzml_loader_jaws.py /usr/local/bin/mzml_loader_jaws.py
