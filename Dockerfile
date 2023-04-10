###################
# STAGE 1: builder
###################

# Build currently doesn't work on > Java 11 (i18n utils are busted) so build on 8 until we fix this
FROM adoptopenjdk/openjdk8:x86_64-debianslim-jre8u345-b01 as builder

WORKDIR /app/source

ENV PATH="$PATH:/opt/conda/bin:/opt/conda/envs/venv/bin"
ENV FC_LANG en-US
ENV LC_CTYPE en_US.UTF-8

# For a set of softwares which can be installed by conda
RUN apt-get update && apt-get install -y coreutils bash git wget make gettext

# For quartet-dseqc-report
## lein:    backend dependencies and building
ADD ./bin/lein /usr/local/bin/lein
RUN chmod 744 /usr/local/bin/lein
RUN lein upgrade

## install dependencies before adding the rest of the source to maximize caching

## backend dependencies
ADD project.clj .
RUN lein deps

## add the rest of the source
ADD . .
## Fetch all submodule
RUN git submodule update --init --recursive

## build the app
RUN lein uberjar

# For VBT
RUN apt-get update \
 && apt-get install -y --force-yes --no-install-recommends\
      groff \
      g++ \
      wget \
      build-essential \
      zlib1g-dev \
      libbz2-dev \
      liblzma-dev \
 && rm -rf /var/lib/apt/lists/*;

## Add htslib-1.6
WORKDIR /home
RUN wget https://github.com/samtools/htslib/releases/download/1.6/htslib-1.6.tar.bz2 --no-check-certificate
RUN tar xvjf htslib-1.6.tar.bz2
WORKDIR /home/htslib-1.6

## Compile htslib-1.6
RUN ./configure
RUN make
RUN make install

## Add vbt folder
RUN git clone https://github.com/sbg/VBT-TrioAnalysis.git /home/varbenchtools
WORKDIR /home/varbenchtools

RUN cp -R /home/htslib-1.6/htslib /home/varbenchtools/
RUN cp /home/htslib-1.6/libhts.a /home/varbenchtools/lib/
RUN cp /home/htslib-1.6/libhts.so /home/varbenchtools/lib/
RUN cp /home/htslib-1.6/libhts.so.2 /home/varbenchtools/lib/

## Compile vbt
RUN make
RUN ldconfig

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py37_22.11.1-1-Linux-x86_64.sh -O miniconda.sh && bash miniconda.sh -b -p /opt/conda
RUN /opt/conda/bin/conda install -c conda-forge -c bioconda -c anaconda mamba blas lapack cxx-compiler conda-pack gfortran_linux-64
## Note: cromwell==83 must not deleted.
## hap.py must be ran in python2.7
RUN /opt/conda/bin/mamba create -n venv -c bioconda -c conda-forge rtg-tools==3.12.1 hap.py==0.3.14 bedtools==2.27.1 picard==2.25.4 fastqc==0.11.8 fastq-screen==0.13.0 qualimap==2.1.1 cromwell==83

# Pack the conda environment
RUN conda-pack -n venv -o /tmp/env.tar && \
  mkdir /venv && cd /venv && tar xf /tmp/env.tar && \
  rm /tmp/env.tar

RUN /venv/bin/conda-unpack

# Sentieon Genomics
RUN wget https://s3.amazonaws.com/sentieon-release/software/sentieon-genomics-201911.01.tar.gz && tar xzvf sentieon-genomics-201911.01.tar.gz && mv sentieon-genomics-201911.01 /opt/sentieon-genomics

# ###################
# # STAGE 2: runner
# ###################

# FROM adoptopenjdk/openjdk8:x86_64-debianslim-jre8u345-b01 as runner
FROM adoptopenjdk/openjdk11:x86_64-debianslim-jre-11.0.18_10 as runner

LABEL org.opencontainers.image.source https://github.com/chinese-quartet/quartet-dseqc-report.git

# hap.py need python2.7, so we must place /venv/bin before /opt/conda/bin
ENV PATH="/venv/bin:/opt/conda/bin:/varbenchtools:/opt/sentieon-genomics/bin:$PATH"
ENV LD_LIBRARY_PATH="/varbenchtools/lib/:$LD_LIBRARY_PATH"
ENV PYTHONDONTWRITEBYTECODE=1
ENV FC_LANG en-US
ENV LC_CTYPE en_US.UTF-8

RUN apt-get update && apt-get install -y coreutils bash git wget make gettext
RUN echo "**** Install dev packages ****" && \
    apt-get update && \
    apt-get install -y curl && \
    \
    echo "*** Install common development dependencies" && \
    apt-get install -y libmariadb-dev libxml2-dev libcurl4-openssl-dev libssl-dev && \
    \
    echo "**** Cleanup ****" && \
    apt-get clean

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py37_22.11.1-1-Linux-x86_64.sh -O miniconda.sh && bash miniconda.sh -b -p /opt/conda
RUN /opt/conda/bin/conda install -y python=3.9 multiqc==1.8
## For app render.
RUN /opt/conda/bin/pip install git+https://github.com/yjcyxky/biominer-app-util.git
ADD ./resources/requirements.txt /data/requirements.txt
ADD ./bin/dseqc.py /opt/conda/bin/dseqc.py
ADD ./bin/quartet-dseqc-report /opt/conda/bin/quartet-dseqc-report
RUN /opt/conda/bin/pip install -r /data/requirements.txt

WORKDIR /data

COPY --from=builder /app/source/target/uberjar/quartet-dseqc-report*.jar /quartet-dseqc-report.jar
COPY --from=builder /venv /venv
COPY --from=builder /app/source/wes-workflow /venv/wes-workflow
COPY --from=builder /app/source/wgs-workflow /venv/wgs-workflow
COPY --from=builder /home/varbenchtools /varbenchtools
COPY --from=builder /opt/sentieon-genomics /opt/sentieon-genomics

## Config file for cromwell instance
COPY --from=builder /app/source/build/cromwell-local.conf /venv/cromwell-local.conf

# Run it
ENTRYPOINT ["/opt/conda/bin/dseqc.py"]