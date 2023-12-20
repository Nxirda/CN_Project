# Utilisation de l'image Ubuntu
FROM ubuntu:20.04

# Mise à jour des packages et installation des dépendances
RUN apt-get update
RUN apt-get install -y gcc 
RUN apt-get install libblas-dev -y
RUN apt-get install build-essential -y
RUN apt-get install liblapacke-dev -y 
RUN apt-get install liblapack-dev -y 

VOLUME /root/env
WORKDIR /root/env
ENV HOSTNAME=docker
