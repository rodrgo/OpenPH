
# Javaplex is needed to create simplicial complexes 
git clone --depth=1 https://github.com/rodrgo/javaplex numerics/dependencies/javaplex
ant -buildfile numerics/dependencies/javaplex/build.xml

# PHAT is needed for benchmarking (this is our baseline)
git clone https://bitbucket.org/phat-code/phat.git numerics/dependencies/phat
cd numerics/dependencies/phat
cmake .
make
git reset --hard e9800d103fcdd19f0417e89781761f7f85d8ec9b 
cd ../..

# Large benchmarking datasets
bash ./numerics/datasets/wget_pointcloud_datasets.sh


