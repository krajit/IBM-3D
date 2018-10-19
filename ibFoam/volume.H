scalarField volume_(connectivity_matrix_.size(),0.0) ;	 
{
    

	forAll(volume_,i)
	{
		vector  p0=pointcloud_[connectivity_matrix_[i][0]];
	  	vector  p1=pointcloud_[connectivity_matrix_[i][1]];
	  	vector  p2=pointcloud_[connectivity_matrix_[i][2]];
	  	vector  p3=pointcloud_[connectivity_matrix_[i][3]];

		volume_[i]= (1.0/6.0)*fabs( ((p2-p0)^(p3-p0))& (p1-p0) );
	}
}

