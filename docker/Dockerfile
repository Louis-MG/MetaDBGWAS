FROM ubuntu:latest
ENV DEBIAN_FRONTEND=noninteractive

COPY ./install_docker.sh ./
RUN chmod +x ./install_docker.sh && sh ./install_docker.sh

ENTRYPOINT ["/Metadbgwas/metadbgwas.sh"]
ENV PATH="/Metadbgwas/:${PATH}"
ENV OPENSSL_CONF=/dev/null
ENV HDF5_USE_FILE_LOCKING=FALSE
