/****************************FORCE SPREADING OPERATION************************************/
  

    //FOR LOOP RUNNING OVER ALL CLOUD POINTS

		forAll(temporary_pointcloud_,i)
			{

			 
			 

			   vector probePoint=current_pointcloud_[i];               
			   vector probePoint_force=force[i];

			    
			   labelHashSet cell_neighbors=findNeighbourCells(probePoint);

			  	
							forAllConstIter(labelHashSet, cell_neighbors,iter) 
					   	 		{
									  label celli = iter.key();
									  vector cellicentre=mesh_.cellCentres()[celli];
									  scalar W=weight(cellicentre,probePoint);


                                      //Info<<"Weight of neighbour cell"<<celli<<"is"<<W<<endl;
								 	    
                                      f[celli]=f[celli]+W*probePoint_force;
                                      //Info<<"f of neighbour cell"<<celli<<"is"<<f[celli]<<endl;
							
						 		}
 

             }

//Info<<"I am running till this point"<<endl;
//Info<<"grid_force"<<f<<endl;

/***********************************************************************************************/
