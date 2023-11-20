docker run \
--rm \
-e JUPYTER_ENABLE_LAB=yes \
-v $PWD:/home/jovyan \
-p 8888:8888 \
geertvangeest/ngs-introduction-jupyter:2021.5 \
start-notebook.sh

docker run \
--rm \
-d \
-p 8787:8787 \
-e PASSWORD=test \
-v $PWD:/home/rstudio \
rocker/rstudio

docker run \
--rm \
-p 8443:8443 \
-e PASSWORD=test \
-e PUID=1000 \
-e PGID=1000 \
-e DEFAULT_WORKSPACE=/config/workspace \
-v $PWD:/config/workspace \
geertvangeest/ngs-introduction-vscode:2023.11

docker run \
--rm \
-p 3000:3000 \
-e PUID=1000 \
-e PGID=1000 \
-e DEFAULT_WORKSPACE=/config/workdir \
-v $PWD:/config/workdir \
openvscode

docker run \
--rm \
-p 3000:3000 \
--name vscode_test \
-e PASSWORD=test \
-e PUID=1000 \
-e PGID=1000 \
-e DEFAULT_WORKSPACE=/config/workdir \
--mount type=bind,source=$PWD,target=/config/workdir \
openvscode