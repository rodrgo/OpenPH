while read p; do
      wget "https://github.com/n-otter/PH-roadmap/raw/master/data_sets/roadmap_datasets_point_cloud/$p" -O "pointclouds/$p"
done <pointcloud_datasets.list
