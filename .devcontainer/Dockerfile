
# See here for image contents: https://github.com/microsoft/vscode-dev-containers/tree/v0.154.2/containers/cpp/.devcontainer/base.Dockerfile

# [Choice] Debian / Ubuntu version: debian-10, debian-9, ubuntu-20.04, ubuntu-18.04
ARG VARIANT="buster"
FROM mcr.microsoft.com/vscode/devcontainers/cpp:0-${VARIANT}

# [Optional] Uncomment this section to install additional packages.
# RUN apt-get update && export DEBIAN_FRONTEND=noninteractive \
#     && apt-get -y install --no-install-recommends <your-package-list-here>

ARG vtk_bin="/tmp/vtk_bin"
ARG vtk_src="/tmp/vtk_src"

COPY ./thirdparty/vtk ${vtk_src}

RUN mkdir ${vtk_bin}
RUN  cmake -S ${vtk_src} -B ${vtk_bin} \
    -DCMAKE_BUILD_TYPE=Release \
    -DBUILD_SHARED_LIBS=ON \
    -DVTK_BUILD_TESTING=OFF \
    -DVTK_WRAP_PYTHON=OFF \
    -DVTK_WRAP_JAVA=OFF \
    -DVTK_GROUP_ENABLE_Qt=DONT_WANT \
    -DVTK_GROUP_ENABLE_Web=DONT_WANT \
    -DVTK_SMP_IMPLEMENTATION_TYPE=OpenMP \
    -DVTK_GROUP_ENABLE_Rendering=DONT_WANT

RUN cd ${vtk_bin} && make . install -j4
RUN rm -d ${vtk_bin} -r
RUN rm -d ${vtk_src} -r

