# BUILD
docker container rm qenum-container 
docker rmi qenum-image
docker build --tag qenum-image -f Dockerfile .
