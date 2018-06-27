# Project 1 - Planar Segmentation
# 3D Photography 
# CUNY Graduate Center
# Prof Stamos

# Task
 Given a range image R, segment the image into planar components

# Input
 Range image R in a .ptx format. The first 10 lines are the header lines that include the number of scan lines, pounts per scan line, and the transformations. Each point has x, y, z coordinates, with respect to the scanner, and the intensity of the point. Range image has a grid structure. 

# Region-growing algorithm
 Sequentially label the points on the same plane with the same label by traversing the cloud using its grid structure. Whether two points lie on the same plane is dicided by comparing the normals of the points. Normal for a point is computed by first taking the sum of the covariance matrices of the points in a small neighborhood around the point. Then eigenvalue decomposition is done and the eigenvector corresponding to the smallest eigenvalue is chosen as the normal. 

--Insert images for Sequntial labeling. 

# RANSAC

- Apply RANSAC by selecting 3 points and then defining a plane using the points.Vote each plane by computing the distance between all the points and the plane. Points with a distance larger than the distance_threshold do not give the plane a vote. Planes with the highest votes are chosen as planar regions.  

--Insert images for 3 point ransac

- Apply RANSAC by selecting a point and its normal to define the plane. Once again, vote the planes by computing the point to plane distance and giving votes for points that are within the distance threshold. Chose planes with the highest votes as the planar region.
