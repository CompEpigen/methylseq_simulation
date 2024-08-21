# syntax=docker/dockerfile:1
FROM ubuntu:20.04
ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y python3.9 python3-pip

RUN mkdir requirements/
RUN mkdir src/

COPY requirements.txt requirements/
COPY *.py src/

RUN pip install -r requirements/requirements.txt
WORKDIR src/
