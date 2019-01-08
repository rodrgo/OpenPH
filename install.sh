
# Dependencies live here
mkdir dependencies && cd "$_"

# Javaplex is needed to create simplicial complexes 
git clone https://github.com/rodrgo/javaplex
cd javaplex
ant

# PHAT is needed for benchmarking (this is our baseline)
git clone https://bitbucket.org/phat-code/phat.git
cd phat
cmake .
make
git reset --hard e9800d103fcdd19f0417e89781761f7f85d8ec9b 

# Large benchmarking datasets
sh datasets/wget_pointcloud_datasets.sh


