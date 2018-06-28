# Project 1 - Planar Segmentation
3D Photography 

CUNY Graduate Center

Professor Stamos

# Task
 Given a range image R, segment the image into planar components

# Input
 Range image R in a .ptx format. The first 10 lines are the header lines that include the number of scan lines, pounts per scan line, and the transformations. Each point has x, y, z coordinates, with respect to the scanner, and the intensity of the point. Range image has a grid structure. 

# Region-growing algorithm
 Sequentially label the points on the same plane with the same label by traversing the cloud using its grid structure. Whether two points lie on the same plane is decided by comparing the normals of the points. Normal for a point is computed by first taking the sum of the covariance matrices of the points in a small neighborhood around the point. Then eigenvalue decomposition is done and the eigenvector corresponding to the smallest eigenvalue is chosen as the normal. 

![Alt text](./Seq_label.jpeg?raw=true "Planar segmentation via sequential labeling")

# RANSAC

- Apply RANSAC by selecting 3 points and then defining a plane using the points.Vote each plane by computing the distance between all the points and the plane. Points with a distance larger than the distance_threshold do not give the plane a vote. Planes with the highest votes are chosen as planar regions.  

![Alt text](./3-pt-plane.jpeg?raw=true "Planar segmentation via defining a plane using 3 points")


- Apply RANSAC by selecting a point and its normal to define the plane. Once again, vote the planes by computing the point to plane distance and giving votes for points that are within the distance threshold. Chose planes with the highest votes as the planar region.

![Alt text](./1-pt.jpeg?raw=true "Planar segmentation via defining a plane using 1 point and its normal")


Both RANSAC algoritms were ran for 500 iterations (on average took 18 hours). 
For the examples with complicated structures, more iterations are required to achieve reasonable planar segmentation. 

Make file included.
 
Need to install the eigen library to run the code: http://eigen.tuxfamily.org/index.php?title=Main_Page

After compiling, run ./p1 input_file.ptx output_1.ptx output_2.ptx output_3.ptx
