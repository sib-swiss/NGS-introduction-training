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
