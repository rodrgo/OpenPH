#!/usr/bin/env bash

# --------
# Javaplex
# --------

git clone --depth=1 https://github.com/rodrgo/javaplex dependencies/javaplex
ant -buildfile dependencies/javaplex/build.xml

# --------
# PHAT 
# --------

git clone https://bitbucket.org/phat-code/phat.git dependencies/phat
cd dependencies/phat
cmake .
make
git reset --hard e9800d103fcdd19f0417e89781761f7f85d8ec9b 
cd ../..

# --------
# DIPHA
# --------

git clone https://github.com/DIPHA/dipha.git dependencies/dipha
cd dependencies/dipha
cmake .
make
cd ../..

# --------
# Datasets
# --------

#fileid="0B7Yz6TPEpiGERGZFbjlXaUt1ZWM"
#filename="datasets/benchmark_data.zip"
#curl -c ./cookie -s -L "https://drive.google.com/uc?export=download&id=${fileid}" > /dev/null
#curl -Lb ./cookie "https://drive.google.com/uc?export=download&confirm=`awk '/download/ {print $NF}' ./cookie`&id=${fileid}" -o ${filename}
#unzip "${filename}" -d benchmark_data
#rm "${filename}"

#DIR="$(cd "$(dirname "$1")"; pwd -P)/$(basename "$1")"
#POINTCLOUDS="${DIR}datasets/pointclouds"
#mkdir "${POINTCLOUDS}"
#
#arr=("celegans_weighted_undirected_reindexed_for_matlab.txt_maxdist_2.6429_SP_distmat.txt_point_cloud.txt" 
#        "dragon_vrip.ply.txt_1000_.txt" 
#        "dragon_vrip.ply.txt_2000_.txt" 
#        "fractal_9_5_2_linear_edge_list.txt_1866.1116_point_cloud.txt" 
#        "fractal_9_5_2_random_edge_list.txt_0.19795_point_cloud.txt" 
#        "fractal_9_5_2_weight_one_edge_list.txt_2_point_cloud.txt" 
#        "H3N2.all.nt.concat.fa_hdm.txt_point_cloud.txt" 
#        "HIV1_2011.all.nt.concat.fa_hdm.txt_point_cloud.txt" 
#        "house104_edge_list.txt_0.72344_point_cloud.txt" 
#        "human_gene2_sampled_reindexed_for_matlab.txt_maxdist_91.9097_SP_distmat.txt_point_cloud.txt" 
#        "klein_bottle_pointcloud_new_400.txt" 
#        "klein_bottle_pointcloud_new_900.txt" 
#        "network379_edge_list.txt_38.3873_point_cloud.txt" 
#        "random_point_cloud_1000_8_.txt" 
#        "random_point_cloud_100_4_.txt" 
#        "random_point_cloud_50_16_.txt" 
#        "senate104_edge_list.txt_0.68902_point_cloud.txt" 
#        "Vicsek__particles_300_distance_1_noise_0.1_v0_0.03_box_25_timestep_1500_of_3000.txt" 
#        "Vicsek__particles_300_distance_1_noise_0.1_v0_0.03_box_25_timestep_3000_of_3000.txt" 
#        "Vicsek__particles_300_distance_1_noise_0.1_v0_0.03_box_5_timestep_150_of_300.txt" 
#        "Vicsek__particles_300_distance_1_noise_0.1_v0_0.03_box_5_timestep_300_of_300.txt" 
#        "Vicsek__particles_300_distance_1_noise_2_v0_0.03_box_7_timestep_300_of_600.txt" 
#        "Vicsek__particles_300_distance_1_noise_2_v0_0.03_box_7_timestep_600_of_600.txt" 
#        )
#for p in "${arr[@]}"
#do
#        echo "$p"
#        wget -q --show-progress "https://github.com/n-otter/PH-roadmap/raw/master/data_sets/roadmap_datasets_point_cloud/$p" -O "${POINTCLOUDS}/$p"
#done
#
