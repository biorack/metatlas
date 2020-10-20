FROM ubuntu:16.04

RUN apt-get update && \
    apt-get -y install unzip bzip2
RUN mkdir -p /wineprefix64/drive_c/pwiz/skyline
ADD pwiz-bin-windows-*.tar.bz2 /wineprefix64/drive_c/pwiz
ADD SkylineTester.zip /
RUN unzip SkylineTester.zip && mv /SkylineTester\ Files/* /wineprefix64/drive_c/pwiz/skyline && rm -fr /wineprefix64/drive_c/pwiz/skyline/TestZipFiles


FROM chambm/wine-dotnet:4.7-x64
COPY --from=0 /wineprefix64/drive_c/pwiz /wineprefix64/drive_c/pwiz

ENV CONTAINER_GITHUB=https://github.com/ProteoWizard/container

LABEL description="Convert MS RAW vendor files to open formats or analyze them with Skyline."
LABEL website=https://github.com/ProteoWizard/container
LABEL documentation=https://github.com/ProteoWizard/container
LABEL license=https://github.com/ProteoWizard/container
LABEL tags="Metabolomics,Proteomics,MassSpectrometry"

ENV WINEDEBUG -all,err+all
ENV WINEPATH "C:\pwiz;C:\pwiz\skyline"

# sudo needed to run wine when container is run as a non-default user (e.g. -u 1234)
# wine*_anyuser scripts are convenience scripts that work like wine/wine64 no matter what user calls them
RUN apt-get update && \
    apt-get -y install sudo && \
    apt-get -y clean && \
    echo "ALL     ALL=NOPASSWD:  ALL" >> /etc/sudoers && \
    echo '#!/bin/sh\nsudo -E -u root wine64 "$@"' > /usr/bin/wine64_anyuser && \
    echo '#!/bin/sh\nsudo -E -u root wine "$@"' > /usr/bin/wine_anyuser && \
    chmod ugo+rx /usr/bin/wine*anyuser && \
    rm -rf \
      /var/lib/apt/lists/* \
      /usr/share/doc \
      /usr/share/doc-base \
      /usr/share/man \
      /usr/share/locale \
      /usr/share/zoneinfo

# create UIDs that Galaxy uses in default configs to launch docker containers; the UID must exist for sudo to work
RUN groupadd -r galaxy -g 1450 && \
    useradd -u 1450 -r -g galaxy -d /home/galaxy -c "Galaxy user" galaxy && \
    useradd -u 1000 -r -g galaxy -d /home/galaxy -c "Galaxy docker user" galaxy_docker && \
    useradd -u 2000 -r -g galaxy -d /home/galaxy -c "Galaxy Travis user" galaxy_travis && \
    useradd -u 999 -r -g galaxy -d /home/galaxy -c "usegalaxy.eu user" galaxy_eu
    
# Set up working directory and permissions to let user xclient save data
RUN mkdir /data
WORKDIR /data

CMD ["wine64_anyuser", "msconvert" ]

## If you need a proxy during build, don't put it into the Dockerfile itself:
## docker build --build-arg http_proxy=http://proxy.example.com:3128/  -t repo/image:version .

ADD mywine /usr/bin/
RUN chmod ugo+rx /usr/bin/mywine
