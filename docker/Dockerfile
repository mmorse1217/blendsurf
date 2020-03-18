# defines base image containing dependencies
FROM ubuntu:18.04 as blendsurf-deps
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \ 
    cmake \ 
    extra-cmake-modules \
    freeglut3-dev \
    libblas-dev \
    liblapack-dev &&\
    apt-get purge -y curl && \
    apt-get autoremove -y && \
    apt-get purge -y ca-certificates &&\
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*


# defines CI build: checks that blendsurf core and renderer still compile
FROM blendsurf-deps as blendsurf-build
COPY . blendsurf/
COPY docker/entrypoint.sh /entrypoint.sh

ENTRYPOINT ["/entrypoint.sh"]

FROM blendsurf-deps as blendsurf-dev
CMD ["/bin/bash"]



