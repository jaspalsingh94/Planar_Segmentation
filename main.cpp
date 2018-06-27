#include <iostream>
#include <fstream>
#include <string>
#include <eigen3/Eigen/Dense>
#include <typeinfo>
#include <unordered_map>
#include <DisjSets.h>
#include <cmath>
#include <cstdlib>
#include <ctime>

using namespace Eigen;
using namespace std;

#define EIGEN_DONT_VECTORIZE
#define EIGEN_DISABLE_UNALIGNED_ARRAY_ASSERT

bool equal_normals(Vector3f &norm_1, Vector3f &norm_2);
bool not_zero_vector(Vector3f &vec);
void set_equivalence(int label1, int label2, DisjSets &eqSet);

int main(int argc, char **argv){

  if(argc != 5){
    printf("Usage: %s input.ptx out_1.ptx out_2.ptx out_3.ptx\n", argv[0]);
    return 0;
  }
  const string input_file(argv[1]);
  const string out_1(argv[2]);
  const string out_2(argv[3]);
  const string out_3(argv[4]);
  float dist_threshold_1 = 0.1;
  int num_of_iter = 50;
  //Read the point cloud file.
  fstream pclfile;
  pclfile.open(input_file);
  if (!pclfile){
    cerr << "Unable to open the Point cloud file" << endl;
  }

  // Get the # of scanlines and points per scanline
  size_t scan_lines;
  size_t pts_per_scanline;
  pclfile >> scan_lines;
  pclfile >> pts_per_scanline;

  string transform_lines[8]; //skip the 8 lines of transformations
  for(int i =0; i < 8; i++)
  { getline(pclfile, transform_lines[i]); }

  //Declare a 3-dimensional dynamic array to store the PCL data. ->For large ex: stack overflows.
  float*** pclarray = new float**[scan_lines];
  for(size_t i = 0; i < scan_lines; i++)
  {
    pclarray[i] = new float*[pts_per_scanline];
    for(size_t j = 0; j < pts_per_scanline; j++)
    {
      pclarray[i][j] = new float[4];
    }
  }

  // Declare a 2-dimensional dynamic array of Vectors to store the normals.
  Vector3f** normal_vectors = new Vector3f*[scan_lines];
  for(size_t i = 0; i < scan_lines; i++)
  { normal_vectors[i] = new Vector3f[pts_per_scanline];}


  Matrix3f cov_mat; //Declare the covariance matrix for normal computation.
  Vector3f eigen_vals; //Vector for storing eigenvalues
  Vector3f point_array[3][3]; //3x3 matrix of vectors to compute the normals.
  Matrix3f eigen_vecs; // Matrix for storing eigenvectors.

  cov_mat << 0, 0, 0, 0, 0, 0, 0, 0, 0;

  float dummy_1;
  for(size_t i = 0; i < scan_lines; i++){
    for(size_t j = 0; j < pts_per_scanline; j++){
      normal_vectors[i][j] << 0, 0, 0; //set the initial normal vector matrix to 0's

      for(int k = 0; k < 4; k++){
        if(k != 3)
        { pclfile >> pclarray[i][j][k]; }
        else {
          pclfile >> dummy_1;
          pclarray[i][j][k] = 0; //set the intensities to 0.
        }
      }
    }
  }
  cout << endl << "Computing normals... " << endl;
  size_t number_of_pts = 0;
  size_t normals_not_computed = 0;
  for(size_t i = 1; i < (scan_lines -1); i++){
  	for(size_t j = 1; j < (pts_per_scanline -1); j++){
  		/*For each point, compute the normal by first computing the covariance matrix of each point in a 3x3 neighborhood.
		  Then add the covariance matrices together for that neighborhood. Then compute the eigenvalues and the eigenvector: eigenvector
		  corresponding to the smalles eigenvalue is the normal we are going to keep.*/
  		cov_mat << 0, 0, 0, 0, 0, 0, 0, 0, 0;

  		if((pclarray[i][j][0] != 0) || (pclarray[i][j][1] != 0) || (pclarray[i][j][2]) != 0){
        number_of_pts++;
  			int m = i-1;
  			int n = j-1;
        //compute the center
        float center[3] = {0, 0, 0};
        for(int k = 0; k < 3; k++){
          for(int l = 0; l < 3; l++){
              center[0] += pclarray[m+k][n+l][0];
              center[1] += pclarray[m+k][n+l][1];
              center[2] += pclarray[m+k][n+l][2];
          }
        }
        center[0] /= 9;
        center[1] /= 9;
        center[2] /= 9;
  			for(int k =0; k < 3; k++){
  				for(int l =0; l< 3; l++){
            // Subtract each point from its center
  					point_array[k][l] << (pclarray[m+k][n+l][0] - center[0]),
  										 (pclarray[m+k][n+l][1] - center[1]),
  										 (pclarray[m+k][n+l][2] - center[2]);
            // Compute the covariance matrix by P * P.Transpose() -> symmetric matrix
  					cov_mat += point_array[k][l] * point_array[k][l].transpose();
  				}
  			} //end of double for loop for computing the covariance.
        // out = 0;

        /* Compute the eigenvalues and the eigenvectors. The eigenvalues come in sorted order since we have
           a symmetric matrix. We take the eigenvector corresponding to the smalles eigenvalue if one exists.
           We throw away the eigenvalues larger than some threshold and keep those normals = 0. */
        SelfAdjointEigenSolver<Matrix3f> eg(cov_mat);
        /* EigenSolver does not always converge, so we check to see if it did. If it doesnt, the normal stays 0.*/
        if (eg.info() == Success){
          eigen_vals = eg.eigenvalues();
          eigen_vecs = eg.eigenvectors();
          // Now we take the non-zero eigenvalue that is atleast less than 100 and keep the corresponding eigenvector.
          if ((eigen_vals(0) > 0) && (eigen_vals(0) < 100))
          { normal_vectors[i][j] << eigen_vecs(0, 0), eigen_vecs(1, 0), eigen_vecs(2, 0); }
          else if ((eigen_vals(1) > 0) && (eigen_vals(1) < 100))
          { normal_vectors[i][j] << eigen_vecs(0, 1), eigen_vecs(1, 1), eigen_vecs(2, 1); }
          else if ((eigen_vals(2) > 0) && (eigen_vals(2) < 100))
          { normal_vectors[i][j] << eigen_vecs(0, 2), eigen_vecs(1, 2), eigen_vecs(2, 2); }
          else {normals_not_computed++;}
          // Done computing the normals.
        }
        else {normals_not_computed++;}
  		}
  	}
  }
  cout << (number_of_pts - normals_not_computed) << " normals computed out of " << number_of_pts << " points." << endl;

  // Sequential labeling
  float label = 1;
  Vector3f this_normal;
  Vector3f up;
  Vector3f diagonal;
  Vector3f left;
  DisjSets EqualLabels(pts_per_scanline*scan_lines*100);
  unordered_map<int, int> equivalence_Fixer;
  //Taking care of the boundry conditions. Second row and second column
  bool first_vec = true;
  for(size_t j = 1; j < (pts_per_scanline-1); j++) //Second row
  {
    this_normal = normal_vectors[1][j];
    if(not_zero_vector(this_normal)){
      if(!first_vec)
      {
        left = normal_vectors[1][j-1];
        if ((not_zero_vector(left)) && (equal_normals(this_normal, left)))
        { pclarray[1][j][3] = pclarray[1][j-1][3];}
        else
        {
          pclarray[1][j][3] = label;
          label++;
        }
      }
      else
      {
        pclarray[1][j][3] = label;
        label++;
        first_vec = false;
      }
    }
  }

  first_vec = true;
  for(size_t i = 0; i < (scan_lines -1); i++)
  {
    this_normal = normal_vectors[i][1];
    if(not_zero_vector(this_normal))
    {
      if(!first_vec)
      {
        up = normal_vectors[i-1][1];
        if ((not_zero_vector(up)) && (equal_normals(this_normal, up)))
        { pclarray[i][1][3] = pclarray[i-1][1][3];}
        else
        {
          pclarray[i][1][3] = label;
          label++;
        }
      }
      else
      {
        pclarray[i][1][3] = label;
        label++;
        first_vec = false;
      }
    }
  }//Done with boundary conditions.

  int left_label;
  int diagonal_label;
  int up_label;
  for(size_t i = 2; i < (scan_lines -1); i++){
    for(size_t j = 2; j< (pts_per_scanline -1); j++){

      this_normal = normal_vectors[i][j];
      bool no_label = true;
      if (not_zero_vector(this_normal))
      {
        up = normal_vectors[i-1][j];
        diagonal = normal_vectors[i-1][j-1];
        left = normal_vectors[i][j-1];
        up_label = pclarray[i-1][j][3];
        diagonal_label = pclarray[i-1][j-1][3];
        left_label = pclarray[i][j-1][3];
        if((not_zero_vector(up)) && (equal_normals(this_normal, up)))
        {
          pclarray[i][j][3] = pclarray[i-1][j][3];
          no_label = false;
        }
        if((not_zero_vector(diagonal)) && (equal_normals(this_normal, diagonal)) && (diagonal_label != up_label))
        {
          if(no_label)
          {
            pclarray[i][j][3] = pclarray[i-1][j-1][3];
            no_label = false;
          }
          else
          { set_equivalence(diagonal_label, up_label, EqualLabels); }
        }
        if((not_zero_vector(left)) && (equal_normals(this_normal, left)) && (diagonal_label != left_label))
        {
          if(no_label)
          {
            pclarray[i][j][3] = pclarray[i][j-1][3];
            no_label = false;
          }
          else
          { set_equivalence(left_label, pclarray[i][j][3], EqualLabels); } //Need to set equivalence between left and this label.
        }
        if(no_label)
        {
          pclarray[i][j][3] = label;
          label++;
        }
      }
    }
  }//Done labeling. Need to clean up the labels next.


  /*Need to do a second pass through the scan to update all the pontss to their label values.*/
  int diff_normals = 1;
  int update_label;
  for(size_t j = 1; j < (pts_per_scanline - 1); j++){
    for(size_t i = 1; i < (scan_lines - 1); i++)
    {
      label = pclarray[i][j][3];
      if(label != 0)
      {
        update_label = EqualLabels.find(label);
        unordered_map<int, int>::iterator itr = equivalence_Fixer.find(update_label);
        //If label does not exist in the hash table, we add it.
        if(itr == equivalence_Fixer.end())
        {
          diff_normals++;
          equivalence_Fixer[update_label] = diff_normals;
          update_label = diff_normals;
        }
        else
        {
          update_label = itr->second;
        }
        pclarray[i][j][3] = update_label;
      }
    }
  }
  cout << "Number of labels for sequential labeling : " << update_label << endl;
  //write the pointcloud data for sequential labeling to the first file.
  ofstream pcl_out_file_1(out_1);
  // pcl_out_file_1.open(out_1);
  if(pcl_out_file_1.is_open())
  {
    pcl_out_file_1 << scan_lines << endl;
    pcl_out_file_1 << pts_per_scanline;
    for(int i = 0; i < 8; i++)
    { pcl_out_file_1 << transform_lines[i] << endl;}

    for(size_t i = 0; i < scan_lines; i++){
      for(size_t j = 0; j < pts_per_scanline; j++){
        pcl_out_file_1 << pclarray[i][j][0] << ' ' << pclarray[i][j][1] << ' ' << pclarray[i][j][2] << ' ' << (pclarray[i][j][3]) << endl;
      }
    }
    pcl_out_file_1.close();
    cout << "Done writing sequentialy labeled points to " << out_1 << endl;
  }
  else cout << "Unable to open the first output file." << endl;

  // Finished with part (a) : Computation of normals and sequential labeling.
  // Starting part (b) : Define a plane by selecting 3 points and apply RANSAC.


  cout << "Starting to compute planes by picking 3 points and labeling points on planes. This may take a while..." << endl;
  Vector3f point_1 , point_2, point_3, plane_vec_1, plane_vec_2, normal, plane_to_pt_vec, this_vector, point_on_plane;
  float distance;
  float percent_labeled = 0.0;
  float planes[num_of_iter][8]; // array of planes where each plane holds the value of a, b, c, d, votes, point_x, point_y, point_z -> ax + by + cz + d = 0

  // for(size_t i = 0; i < scan_lines; i++){
  //   for(size_t j = 0; j < pts_per_scanline; j++)
  //   { pclarray[i][j][3] = 0; } } //Reset the labels to 0's.

  srand(time(0));
  int labeled_pts = 0;
  int plane;
  float votes = 0;
  bool first = true;
  for(int iter = 0; iter < num_of_iter; iter++)
  {
    point_1 << 0, 0, 0;
    point_2 << 0, 0, 0;
    point_3 << 0, 0, 0;

    size_t row, col;

    while(!not_zero_vector(point_1))
    { row = rand() % scan_lines;
      col = rand() % pts_per_scanline;
      point_1 << pclarray[row][col][0], pclarray[row][col][1], pclarray[row][col][2];
    }

    while(!not_zero_vector(point_2))
    { row = rand() % scan_lines;
      col = rand() % pts_per_scanline;
      point_2 << pclarray[row][col][0], pclarray[row][col][1], pclarray[row][col][2];
    }

    while(!not_zero_vector(point_3))
    { row = rand() % scan_lines;
      col = rand() % pts_per_scanline;
      point_3 << pclarray[row][col][0], pclarray[row][col][1], pclarray[row][col][2];
    } // So now we have our 3 non-zero points.

    plane_vec_1 = point_2 - point_1;
    plane_vec_2 = point_3 - point_1;
    normal = plane_vec_1.cross(plane_vec_2);
    normal.normalize();

    // We have the normal now -> eq of plane: n1*x + ny*y + nz* z + d = 0
    // d = - (normal * point_1)
    planes[iter][0] = normal(0);
    planes[iter][1] = normal(1);
    planes[iter][2] = normal(2);
    planes[iter][3] = normal.dot(point_1);
    planes[iter][3] *= -1;
    planes[iter][4] = 0;
    planes[iter][5] = point_1(0);
    planes[iter][6] = point_1(1);
    planes[iter][7] = point_1(2);

    // Now we have the plane defined. Need to go through each point and give planes votes for points that are within a threshold.
    for(size_t i = 0; i < scan_lines; i++){
      for(size_t j = 0; j< pts_per_scanline; j++){
        this_vector << pclarray[i][j][0], pclarray[i][j][1], pclarray[i][j][2];
        if(not_zero_vector(this_vector)){
          if(first){
            pclarray[i][j][3] = 0;
          }
          plane_to_pt_vec = this_vector - point_1;
          distance = abs(normal.dot(plane_to_pt_vec));
          if(distance < dist_threshold_1){
            planes[iter][4]++;
          }
        }
      }
    }
    first = false;
    if(votes < planes[iter][4]){
      votes = planes[iter][4];
      plane = iter;
    }
    // cout << "This plane has " << planes[iter][4] << " votes." << endl;
  } //Done computing planes and votes

  // Now we will take the top 1/3 of the planes that were computed and set the points to the corresponding label
  int num_of_planes = (num_of_iter/3) + 1;
  for(int label_of_plane = 1; label_of_plane < num_of_planes; label_of_plane++){

    point_on_plane << planes[plane][5], planes[plane][6], planes[plane][7];
    normal << planes[plane][0], planes[plane][1], planes[plane][2];
    // Now take the plane with the highest votes and set the closest points to it.
    for(size_t i = 0; i < scan_lines; i++){
      for(size_t j = 0; j < pts_per_scanline; j++){
        this_vector << pclarray[i][j][0], pclarray[i][j][1], pclarray[i][j][2];
        //We only do this for points that dont have a label.
        if((not_zero_vector(this_vector)) && (pclarray[i][j][3] == 0))
        {
          plane_to_pt_vec = this_vector - point_on_plane;
          distance = abs(normal.dot(plane_to_pt_vec));
          if(distance < dist_threshold_1){
            pclarray[i][j][3] = label_of_plane;
            labeled_pts++;
          }
        }
      }
    }
    planes[plane][4] = 0; //Set the votes for this plane to 0.
    votes = 0;
    if((label_of_plane-1) % 10 == 0)
    { cout << labeled_pts << " points labeled on " << label_of_plane << " planes." << endl;
      percent_labeled += ((float)(labeled_pts))/((float)(number_of_pts));
      cout << 100 - (percent_labeled*100) << " percent points left for labeling." << endl;
      labeled_pts = 0; }

    // The labels for the plane is set. Now revote the planes without the labeled points.
    for(int iter = 0; iter < num_of_iter; iter++){
      //Set the number of votes to zero before revoting
      if(planes[iter][4] < 55)
        { continue; }
      else
        { planes[iter][4] = 0;}
      normal << planes[iter][0], planes[iter][1], planes[iter][2];
      point_on_plane << planes[iter][5], planes[iter][6], planes[iter][7];

      for(size_t i = 0 ; i < scan_lines; i++){
        for(size_t j = 0; j < pts_per_scanline; j++){
          this_vector << pclarray[i][j][0], pclarray[i][j][1], pclarray[i][j][2];

          if(not_zero_vector(this_vector) && (pclarray[i][j][3] == 0))
          {
            plane_to_pt_vec = this_vector - point_on_plane;
            distance = abs(normal.dot(plane_to_pt_vec));
            if(distance < dist_threshold_1){
              planes[iter][4]++;
            }
          }
        }
      }
      if(votes < planes[iter][4]){
        votes = planes[iter][4];
        plane = iter;
      }
    } // Done revoting each plane. The planes already used will trivially have 0 votes.

    if(votes < 50)
    { break;}
  } //Done laebling the points

  // Write the points and labels to the second output file.
  ofstream pcl_out_file_2(out_2);
  // pcl_out_file_1.open(out_1);
  if(pcl_out_file_2.is_open())
  {
    pcl_out_file_2 << scan_lines << endl;
    pcl_out_file_2 << pts_per_scanline;
    for(int i = 0; i < 8; i++)
    { pcl_out_file_2 << transform_lines[i] << endl;}

    for(size_t i = 0; i < scan_lines; i++){
      for(size_t j = 0; j < pts_per_scanline; j++){
        pcl_out_file_2 << pclarray[i][j][0] << ' ' << pclarray[i][j][1] << ' ' << pclarray[i][j][2] << ' ' << (pclarray[i][j][3]) << endl;
      }
    }
    pcl_out_file_2.close();
    cout << "Done writing planes generated by 3 points to " << out_2 << endl;
  }
  else cout << "Unable to open the first output file." << endl;

  // Done with part b of computing a plane with 3 points and labeling points on the best planes.

  // Starting the computation of generating planes via a normal and 1 point and labeling points on the best plane.
  // the varible normal is already defined so we will reuse it.
  cout << "Starting computation of planes by picking a point and normal and labeling points on plane. This will take some time..." << endl;
  percent_labeled = 0.0;
  num_of_iter = 500;
  dist_threshold_1 -= 0.05;
  labeled_pts = 0;
  size_t row, col;
  votes = 0;
  first = true;
  // for(size_t i = 0; i < scan_lines; i++){
  //   for(size_t j = 0; j < pts_per_scanline; j++)
  //   { pclarray[i][j][3] = 0; } } //Reset the labels to 0's.

  int planes_2[num_of_iter][3];
  for(int iter = 0; iter < num_of_iter; iter++)
  {
    point_on_plane << 0, 0, 0;
    while(!not_zero_vector(point_on_plane)){
      row = rand() % scan_lines;
      col = rand() % pts_per_scanline;
      point_on_plane << pclarray[row][col][0], pclarray[row][col][1], pclarray[row][col][2];
    }
    normal << normal_vectors[row][col];
    planes_2[iter][0] = row;
    planes_2[iter][1] = col;
    planes_2[iter][2] = 0; //planes_2[iter][2] will hold the votes.
    /* Picked a random non-zero point and have its corresponding normal. The equation of a plane is normal * point - d = 0.
    To compute the distance a point is from the plane, we need to project that vector on to the normal. However, both the normal and the vector must have the same tail.
    So the vector needs to be from the point on plane to the point, which means we do not need to compute the d for the equation of the plane, since
    the distance is only dependent on the normal and the vector plane_to_point. */
    // For each point, will compute the distance from plane and vote the plane.

    for(size_t i = 0; i < scan_lines; i++){
      for(size_t j = 0; j < pts_per_scanline; j++){
        this_vector << pclarray[i][j][0], pclarray[i][j][1], pclarray[i][j][2];
        if(not_zero_vector(this_vector))
        {
          if(first){
            pclarray[i][j][3] = 0;
          }
          normal.normalize();
          distance = abs(normal.dot(this_vector - point_on_plane));
          if(distance < dist_threshold_1)
          { planes_2[iter][2]++;}
        }
      }
    }//End of For
    first = false;
    if(votes < planes_2[iter][2]){
      votes = planes_2[iter][2];
      plane = iter;
    }
  }

  // Done computing planes.
  // Label points on the best plane and then revote the planes.
  num_of_planes = (num_of_iter/3) + 1;
  for(int label_of_plane = 1; label_of_plane < num_of_planes; label_of_plane++){

    row = planes_2[plane][0];
    col = planes_2[plane][1];
    point_on_plane << pclarray[row][col][0], pclarray[row][col][1], pclarray[row][col][2];
    normal << normal_vectors[row][col];
    for(size_t i = 0; i < scan_lines; i++){
      for(size_t j = 0; j < pts_per_scanline; j++){
        this_vector << pclarray[i][j][0], pclarray[i][j][1], pclarray[i][j][2];
        if(not_zero_vector(this_vector) && pclarray[i][j][3] == 0)
        { normal.normalize();
          distance = abs(normal.dot(this_vector - point_on_plane));
          if(distance < dist_threshold_1){
            pclarray[i][j][3] = label_of_plane;
            labeled_pts++;
          }
        }
      }
    }//Done labeling points on best plane for this iteration.
    planes_2[plane][2] = 0;
    votes = 0;

    if((label_of_plane-1) % 10 == 0){
      cout << labeled_pts << " pts labeled on " << label_of_plane << " planes. " << endl;
      percent_labeled += ((float)(labeled_pts))/((float)(number_of_pts));
      cout << 100 - (percent_labeled*100) << " percent points left for labeling." << endl;
      labeled_pts = 0;
    }

    //Need to revote the planes.
    for(int iter = 0; iter < num_of_iter; iter++){
      if(planes_2[iter][2] < 55){
        continue; }
      else{
        planes_2[iter][2] = 0;}
      row = planes_2[iter][0];
      col = planes_2[iter][1];
      normal << normal_vectors[row][col];
      point_on_plane << pclarray[row][col][0], pclarray[row][col][1], pclarray[row][col][2];

      for(size_t i = 0; i < scan_lines; i++){
        for(size_t j = 0; j < pts_per_scanline; j++){
          this_vector << pclarray[i][j][0], pclarray[i][j][1], pclarray[i][j][2];

          if(not_zero_vector(this_vector) && (pclarray[i][j][3] == 0)){
            normal.normalize();
            distance = abs(normal.dot(this_vector - point_on_plane));
            if(distance < dist_threshold_1){
              planes_2[iter][2]++;} }
        }
      }
      // Save the best plane after revoting in plane for the next iteration.
      if(votes < planes_2[iter][2]){
        votes = planes_2[iter][2];
        plane = iter;
      }
      // cout << "This plane has " << planes_2[iter][2] << " votes." << endl;
    }//Done revoting the planes.

    if(votes < 50){
      break;}
  }//Done labeling all the points.

  // Write the points and labels to the second output file.
  ofstream pcl_out_file_3(out_3);
  // pcl_out_file_3.open(out_3);
  if(pcl_out_file_3.is_open())
  {
    pcl_out_file_3 << scan_lines << endl;
    pcl_out_file_3 << pts_per_scanline;
    for(int i = 0; i < 8; i++)
    { pcl_out_file_3 << transform_lines[i] << endl;}

    for(size_t i = 0; i < scan_lines; i++){
      for(size_t j = 0; j < pts_per_scanline; j++){
        pcl_out_file_3 << pclarray[i][j][0] << ' ' << pclarray[i][j][1] << ' ' << pclarray[i][j][2] << ' ' << (pclarray[i][j][3]) << endl;
      }
    }
    pcl_out_file_3.close();
    cout << "Done writing planes generated by 3 points to " << out_3 << endl;
  }
  else cout << "Unable to open the first output file." << endl;

  //Free the memory
  for(size_t i = 0; i < scan_lines; i++)
  {
    delete [] normal_vectors[i];
    for(size_t j = 0; j < pts_per_scanline; j++){
      delete [] pclarray[i][j];
    }
    delete [] pclarray[i];
  }
  delete [] pclarray;
  delete [] normal_vectors;

  return 0;
}

//Check for a normal vector being zero
bool not_zero_vector(Vector3f &vec){
  if ((vec(0) == 0) && (vec(1) == 0) && (vec(2) == 0))
  { return false;}
  else
  { return true;}
}

// Compares two normals to see if there dot profuct is within 0.93
bool equal_normals(Vector3f &norm_1, Vector3f &norm_2){
  if (abs(norm_1.dot(norm_2)) > 0.93)
  { return true;}
  else
    return false;
}

// Sets equivalence for 2 labels.
void set_equivalence(int label1, int label2, DisjSets &eqSet)
{
  int first_root = eqSet.find(label1);
  int second_root = eqSet.find(label2);
  if(first_root != second_root)
  {
    if(label1 < label2)
    { eqSet.unionSets(first_root, second_root);}
    else if (label1 > label2)
    { eqSet.unionSets(second_root, first_root);}
  }
}
