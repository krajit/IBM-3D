/*********************************************3D FORCE**************************************************/

           vectorField force(pointcloud_.size(),Zero);           ///Initialising Force at each vertex

           forAll(volume_,i)                                     ///For loop runs over each element
				{

                    //Info<<"element"<<i<<endl;

                    //Coordinates of reference tetrahedron:
					vector  s0=pointcloud_[connectivity_matrix_[i][0]];   
	  				vector  s1=pointcloud_[connectivity_matrix_[i][1]];
	  				vector  s2=pointcloud_[connectivity_matrix_[i][2]];
                    vector  s3=pointcloud_[connectivity_matrix_[i][3]];
 

                    //Coordinates of deformed tetrahedron:
					vector  x0=current_pointcloud_[connectivity_matrix_[i][0]];   
	  				vector  x1=current_pointcloud_[connectivity_matrix_[i][1]];
	  				vector  x2=current_pointcloud_[connectivity_matrix_[i][2]];
                    vector  x3=current_pointcloud_[connectivity_matrix_[i][3]];

                    //difference between coordinates of reference nodes:
                    vector  ds1=s1-s0;
                    vector  ds2=s2-s0;
                    vector  ds3=s3-s0;


                    //Deformation gradient:

                    //First row:

				    vector  a0=(x1[0]-x0[0])*(ds2^ds3)+(x2[0]-x0[0])*(ds3^ds1)+(x3[0]-x0[0])*(ds1^ds2);
                    vector  a1=(x1[1]-x0[1])*(ds2^ds3)+(x2[1]-x0[1])*(ds3^ds1)+(x3[1]-x0[1])*(ds1^ds2);
                    vector  a2=(x1[2]-x0[2])*(ds2^ds3)+(x2[2]-x0[2])*(ds3^ds1)+(x3[2]-x0[2])*(ds1^ds2);
                    
                    //Info<<"volume"<<volume_[i]<<endl;

                    a0=(1.0/(6.0*volume_[i]))*a0;
                    a1=(1.0/(6.0*volume_[i]))*a1;
                    a2=(1.0/(6.0*volume_[i]))*a2;

                    //Info<<"def_grad"<<a0<<endl;
                    //Info<<"def_grad"<<a1<<endl;
                    //Info<<"def_grad"<<a2<<endl;

                    //Piola Kirchoff Stress:

                    //For simplicity we set deformationgradient=Piola Kirchoff stress

                    vector P0=a0;
                    vector P1=a1;
                    vector P2=a2;
                    
                    //Derivative of deformation gradient wrt position:


	               vector dela_delx0 = -(1.0/(6.0*volume_[i]))*((s1^s2) + (s2^s3) + (s3^s1));
	               vector dela_delx1 =  (1.0/(6.0*volume_[i]))*(ds2^ds3);
	               vector dela_delx2 =  (1.0/(6.0*volume_[i]))*(ds3^ds1);
	               vector dela_delx3 =  (1.0/(6.0*volume_[i]))*(ds1^ds2);

                    
                   vector f0( P0&dela_delx0, P1&dela_delx0, P2&dela_delx0);
                   vector f1( P0&dela_delx1, P1&dela_delx1, P2&dela_delx1);
                   vector f2( P0&dela_delx2, P1&dela_delx2, P2&dela_delx2);
                   vector f3( P0&dela_delx3, P1&dela_delx3, P2&dela_delx3);

                   f0 = -volume_[i]*f0;
                   f1 = -volume_[i]*f1;
                   f2 = -volume_[i]*f2;
                   f3 = -volume_[i]*f3;

                   //Info<<"f0"<<f0<<endl;
                   //Info<<"f1"<<f1<<endl;
                   //Info<<"f2"<<f2<<endl;
                   //Info<<"f3"<<f3<<endl;

                   force[connectivity_matrix_[i][0]]=force[connectivity_matrix_[i][0]] + f0;
                   force[connectivity_matrix_[i][1]]=force[connectivity_matrix_[i][1]] + f1;
				   force[connectivity_matrix_[i][2]]=force[connectivity_matrix_[i][2]] + f2;
                   force[connectivity_matrix_[i][3]]=force[connectivity_matrix_[i][3]] + f3;

 
				}  


