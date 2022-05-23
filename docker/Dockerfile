FROM ubuntu:latest
ENV DEBIAN_FRONTEND=noninteractive
USER root
COPY ./install_docker.sh ./
RUN chmod +x ./install_docker.sh && sh ./install_docker.sh
ENV PATH="/MetaDBGWAS/:${PATH}"
