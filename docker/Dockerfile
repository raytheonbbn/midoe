ARG UBUNTU_VERSION=focal
FROM ubuntu:${UBUNTU_VERSION} AS builder

#################################################
# Install system dependencies used for building #
#################################################
RUN apt-get -y update && \
    DEBIAN_FRONTEND=noninteractive \
    apt-get -y --no-install-recommends install \
    build-essential \
    python3 \
    python3-venv

######################################################
# Create python virtual env and install dependencies #
######################################################
ENV VIRTUAL_ENV=/opt/venv
RUN python3 -m venv $VIRTUAL_ENV
ENV PATH="$VIRTUAL_ENV/bin:$PATH"
COPY ./docker/requirements.txt .
RUN pip install --no-cache-dir wheel && pip install --no-cache-dir -r requirements.txt


################################################
# Finished building. Create runnable container #
################################################
FROM ubuntu:${UBUNTU_VERSION}


###############################
# Install system dependencies #
###############################
RUN apt-get -y update && \
    DEBIAN_FRONTEND=noninteractive \
    apt-get -y --no-install-recommends install \
    python3


###########################
# Copy files from builder #
###########################
COPY --from=builder /opt/venv /opt/venv
COPY --from=builder /usr/lib /usr/lib
COPY --from=builder /usr/local /usr/local


###########################################
# Install detection schema python scripts #
###########################################
WORKDIR /opt/felix/detection-schema
COPY ./schema ./schema
COPY ./gff_to_json ./gff_to_json
COPY ./feature_extraction ./feature_extraction
COPY ./fa_to_json ./fa_to_json
COPY ./evidenceintegration ./evidenceintegration
COPY ./*.py ./

# Ensure entrypoint executable
RUN chmod +x /opt/felix/detection-schema/evidenceintegration/DetectionsToCSV.py


################
# Set ENV vars #
################
ENV VIRTUAL_ENV=/opt/venv
ENV PATH="$VIRTUAL_ENV/bin:$PATH"

# Run setup script
RUN python3 setup.py install

# Set workdir
WORKDIR /outdir

# Set entrypoint
ENTRYPOINT ["/opt/felix/detection-schema/evidenceintegration/DetectionsToCSV.py"]
CMD ["--help"]
