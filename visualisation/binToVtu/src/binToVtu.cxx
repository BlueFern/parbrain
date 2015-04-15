#include <stdio.h>

#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>

#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkHexahedron.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>
#include <vtkUnstructuredGrid.h>

#include <vtkDataSet.h>  // really needed?
#include <vtkDataSetAttributes.h>  // really needed?
#include <vtkProperty.h>  // really needed?

#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkXMLPolyDataWriter.h>

#include <vtkLine.h>
#include <vtkTubeFilter.h>
#include <vtkAppendPolyData.h>
#include <vtkLineSource.h>  // really needed?


int main(int argc, char *argv[]) {
	// general parameters
	#define BLOCK_LENGTH 4e-4	// length of one tissue block (m) 4e-4
    int time_fin = 777;  // final time 400
    //TODO: write/read final time in/from configuration file!

//  read in configuration file
    std::ifstream conf_file("/hpc/home/kdo40/Frontiers_in_Physiology/parbrain/ATP_input/np64_nlev13_sbtr03/info.dat", std::ifstream::binary);

    std::string header;
    std::getline(conf_file, header); // skip header line

    int temp_array[7];
    int b;
    int i = 0;
    while (conf_file >> b) {
    	temp_array[i] = b;
    	i++;
    }
    int n_procs  = temp_array[0];
    int n_blocks_per_rank = temp_array[1];
    int n_state_vars = temp_array[2];
    int m_local = temp_array[3];
    int n_local = temp_array[4];
    int m_global = temp_array[5];
    int n_global = temp_array[6];
    int n_blocks = n_procs * n_blocks_per_rank;
    int n_cols = n_local * n_global;  				// number of columns of tissue blocks (j)
    int n_rows = m_local * m_global;  				// number of rows of tissue blocks (i)

	//read tissue block state binary file
//	std::ifstream is("/home/katharina/workspace_GitHub/parbrain/np02_nlev07_sbtr03_3/tissueBlocks.dat", std::ifstream::binary);
	std::ifstream is("/hpc/home/kdo40/Frontiers_in_Physiology/parbrain/ATP_input/np64_nlev13_sbtr03/tissueBlocks.dat", std::ifstream::binary);
	if (!is) {
		std::cerr << "Cannot read tissueBlock.dat file." << std::endl;
		return EXIT_FAILURE;
	}

	//read flow binary file
//	std::ifstream is_flow("/home/katharina/workspace_GitHub/parbrain/np02_nlev07_sbtr03_3/flow.dat", std::ifstream::binary);
	std::ifstream is_flow("/hpc/home/kdo40/Frontiers_in_Physiology/parbrain/ATP_input/np64_nlev13_sbtr03/flow.dat", std::ifstream::binary);
	if (!is_flow) {
		std::cerr << "Cannot read flow.dat file." << std::endl;
		return EXIT_FAILURE;
	}

    // Read x and y coordinates from tissue_block binary file
	double xCoord[n_blocks];
	double yCoord[n_blocks];

	is.read((char *) xCoord, sizeof(xCoord));
	is.read((char *) yCoord, sizeof(yCoord));

	for (int i = 0; i < n_blocks; i++) {
	    	std::cout << "Block: " << i << "\t x-coord: " << xCoord[i] << "  \t y-coord: " << yCoord[i] << std::endl;
	    }

// Create a vtkPoints object and store the points in it
	// Tissue blocks:
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

	double points_xorigin = xCoord[0] - (BLOCK_LENGTH / 2);
	double points_yorigin = yCoord[0] - (BLOCK_LENGTH / 2);

	//lower "layer" of points
	for (int j = 0; j <= n_cols; j++) {
		double points_xcoord = points_xorigin + j * BLOCK_LENGTH;
		for (int i = 0; i <= n_rows; i++) {
			double points_ycoord = points_yorigin + i * BLOCK_LENGTH;
			points->InsertNextPoint(points_xcoord, points_ycoord, 0);
		}
	}
	//upper "layer" of points
	for (int j = 0; j <= n_cols; j++) {
		double points_xcoord = points_xorigin + j * BLOCK_LENGTH;
		for (int i = 0; i <= n_rows; i++) {
			double points_ycoord = points_yorigin + i * BLOCK_LENGTH;
			points->InsertNextPoint(points_xcoord, points_ycoord, BLOCK_LENGTH);
		}
	}

	// H-Tree: (algorithm from adjacency.c)
	vtkSmartPointer<vtkPoints> points2 = vtkSmartPointer<vtkPoints>::New();
    int n_rows_h = n_rows;
    int n_cols_h = n_cols;
	int n_bifr = log2(n_blocks) + 1;      // number of bifurcations
	int n_branches = (1 << n_bifr) - 1;    // number of branches
	int n_nodes = (1 << (n_bifr - 1)) - 1; // number of nodes  ?

    double xpoints_tree[n_branches];    // x-coord for points for branches & pressure
    double ypoints_tree[n_branches];    // y-coord for points for branches & pressure

    int a, k1, k2, row = 0, col = (1 << (n_bifr-1));
        int xbranch = 0;
        int offset = n_nodes + 1; // because we have to add the leaf nodes
        int nx[n_branches], ny[n_branches], nz1[n_nodes], nz2[n_nodes], h1[n_nodes], h2[n_nodes];
        for (int i = 0; i < col; i++) {
    	ny[i] = i;
        }

        for (int L = n_bifr - 1; L > 0; L--) {
            a = (1 << n_bifr) - (1 << (L+1));

            if (xbranch) {
                for (int j = 0; j < n_cols_h; j+=2) {
                    for (int i = 0; i < n_rows_h; i++) {
                        k1 = a + i + j*n_rows_h;
                        k2 = a + i + (j+1)*n_rows_h;
                        nx[k1] = row + n_blocks; nx[k2] = row + n_blocks; ny[col] = row + n_blocks;
    		    h1[row] = k1;
    		    h2[row] = k2;
    		    row++; col++;
                    }
                }
                n_cols_h /= 2;
            }
            else {
                for (int j = 0; j < n_cols_h; j++) {
                    for (int i = 0; i < n_rows_h; i+=2) {
                        k1 = a + i + j*n_rows_h;
                        k2 = k1 + 1;
    		    nx[k1] = row + n_blocks; nx[k2] = row + n_blocks; ny[col] = row + n_blocks;
    		    h1[row] = k1;
    		    h2[row] = k2;
    		    row++; col++;
                    }
                }
                n_rows_h /= 2;
            }
            xbranch = !xbranch;
        } // L loop: from bottom level up to the top of the tree



        for (int i = 0; i < n_nodes; i++) {
            nz1[i] = ny[h1[i]];
            nz2[i] = ny[h2[i]];
        }
        for (int i = 0; i < n_blocks; i++) {
            xpoints_tree[i] = xCoord[i];
            ypoints_tree[i] = yCoord[i];
            points2->InsertNextPoint(xCoord[i], yCoord[i] , BLOCK_LENGTH);
        }

        for (int i = n_blocks; i < n_branches; i++) {
            xpoints_tree[i] = (xpoints_tree[nz1[i-n_blocks]] + xpoints_tree[nz2[i-n_blocks]]) / 2;
            ypoints_tree[i] = (ypoints_tree[nz1[i-n_blocks]] + ypoints_tree[nz2[i-n_blocks]]) / 2;
            points2->InsertNextPoint((xpoints_tree[nz1[i-n_blocks]] + xpoints_tree[nz2[i-n_blocks]]) / 2 , (ypoints_tree[nz1[i-n_blocks]] + ypoints_tree[nz2[i-n_blocks]]) / 2 , BLOCK_LENGTH);
        }
        points2->InsertNextPoint(0, 0, 7*BLOCK_LENGTH);        // root branch


// Create cell array
    // Tissue blocks:
	vtkSmartPointer<vtkCellArray> cellArray = vtkSmartPointer<vtkCellArray>::New();
	vtkSmartPointer<vtkHexahedron> hexahedr = vtkSmartPointer<vtkHexahedron>::New(); 	// Create hexahedrons (3D blocks)
	double npoints_offset = (n_rows + 1) * (n_cols + 1); // number of points - offset for 2nd "layer"
	for (int block_id = 0; block_id < n_blocks; block_id++) {
		int col = block_id / n_rows;     						// which column are we in?
		int i = block_id + col;
		hexahedr->GetPointIds()->SetId(0, i);					// lower SW corner of hexahedron
		int j = i + n_rows + 1;
		hexahedr->GetPointIds()->SetId(1, j); 					// lower SE corner of hexahedron
		int k = j + 1;
		hexahedr->GetPointIds()->SetId(2, k);					// lower NE corner of hexahedron
		int l = k - (n_rows + 1);
		hexahedr->GetPointIds()->SetId(3, l); 					// lower NW corner of hexahedron
		hexahedr->GetPointIds()->SetId(4, i+npoints_offset); 	// upper SW corner of hexahedron
		hexahedr->GetPointIds()->SetId(5, j+npoints_offset); 	// upper SE corner of hexahedron
		hexahedr->GetPointIds()->SetId(6, k+npoints_offset); 	// upper NE corner of hexahedron
		hexahedr->GetPointIds()->SetId(7, l+npoints_offset); 	// upper NW corner of hexahedron

		cellArray->InsertNextCell(hexahedr); // Create a cell array to store hexaedrons in and add them to it
	}

    // H-Tree:
	vtkSmartPointer<vtkLine> lines = vtkSmartPointer<vtkLine>::New();
	vtkSmartPointer<vtkCellArray> cellArray2 = vtkSmartPointer<vtkCellArray>::New();
	for (int line_id = 0; line_id < n_branches-1; line_id++) { // -1, because of root branch
		lines->GetPointIds()->SetId(0, nx[line_id]);
		lines->GetPointIds()->SetId(1, ny[line_id]);
		cellArray2->InsertNextCell(lines);
	}

	lines->GetPointIds()->SetId(0, ny[n_branches-1]); // root branch
	lines->GetPointIds()->SetId(1, ny[n_branches-1]+1);
	cellArray2->InsertNextCell(lines);




// Create unstructured grid
	// Tissue blocks:
	vtkSmartPointer<vtkUnstructuredGrid> uGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
	uGrid->SetPoints(points);
	uGrid->SetCells(VTK_HEXAHEDRON, cellArray);

//	// H-Tree:
//	vtkSmartPointer<vtkUnstructuredGrid> uGrid2 = vtkSmartPointer<vtkUnstructuredGrid>::New();
//	uGrid2->SetPoints(points2);
//	uGrid2->SetCells(VTK_LINE, cellArray2); //     uGrid2->SetLines(cellArray2);




	char *var_names[] = {"radius_terminating_arteriole","R_k","N_Na_k","N_K_k","N_HCO3_k","N_Cl_k","N_Na_s","N_K_s","N_HCO3_s","K_p","w_k","ca_i","ca_sr_i","v_i","w_i","ip3_i","K_i","ca_j","ca_er_j","v_j","ip3_j","Mp","AMp","AM","NOi","NOj","NOn","cGMP","eNOS","nNOS","ca_n","E_b","E_6c","E_5c"};

// Time step loop - read binary data and add as attributes to cells
	for (int i = 0; i < time_fin; i++) {
		double time, time2;
		is.read((char *) &time, sizeof(time));    // read time from both binary files (is not used for anything at the moment...)
		is_flow.read((char *) &time2, sizeof(time2));
		std::cout << "Time: " << time << std::endl;
		// Tissue blocks:
		std::vector<vtkSmartPointer<vtkDoubleArray> > stateVars;  //create array for each state variable
		for (int v = 0; v < n_state_vars; v++) {
			vtkSmartPointer<vtkDoubleArray> array = vtkSmartPointer<vtkDoubleArray>::New();
			array->SetName(var_names[v]);
			stateVars.push_back(array);
		}

		// H-Tree:
		std::vector<vtkSmartPointer<vtkDoubleArray> > flowVar;
		vtkSmartPointer<vtkDoubleArray> array2 = vtkSmartPointer<vtkDoubleArray>::New(); //create one array for flow variables
		array2->SetName("blood_flow");
		flowVar.push_back(array2);

		// Tissue blocks:
		for (int k = 0; k < n_blocks; k++) { 		// read state variables and add as attributes to uGrid
			double tmpArr[n_state_vars];
			is.read((char *) tmpArr, sizeof(tmpArr));  // read state variables from binary file

			for (int v = 0; v < n_state_vars; v++) {
				stateVars[v]->InsertNextValue(tmpArr[v]);
			}
		}

		// H-Tree:
		for (int k2 = 0; k2 < n_branches; k2++) { 		// read flow variables and add as attributes to uGrid2
			double tmpArr2[1];
			is_flow.read((char *) tmpArr2, sizeof(tmpArr2));  // read flow variables from binary file
			flowVar[0]->InsertNextValue(tmpArr2[0]);

		}

		// Tissue blocks:
		for (int v = 0; v < n_state_vars; v++) {
			uGrid->GetCellData()->AddArray(stateVars[v]);
		}

//		// H-Tree:
//	    uGrid2->GetCellData()->AddArray(flowVar[0]);

		// tubeFilter - polyData
		vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
	    polyData->SetPoints(points2);
	    polyData->SetLines(cellArray2);

		// H-Tree(tubeFilter) - try:
	    polyData->GetCellData()->AddArray(flowVar[0]);

//	    // try:changing radius
//	    vtkSmartPointer<vtkDoubleArray> tubeRadius = vtkSmartPointer<vtkDoubleArray>::New();
//	    tubeRadius->SetName("TubeRadius");
//	    tubeRadius->SetNumberOfTuples(polyData->GetNumberOfPoints());
//	    for (int branch_id=0 ;branch_id < polyData->GetNumberOfPoints(); branch_id++){
////	    	tubeRadius->SetTuple1(i,5e-6*(flowVar[0]->GetValue(i))+3e-5);
//	    	tubeRadius->SetTuple1(branch_id,5e-6*(flowVar[0]->GetValue(i))+3e-5);
//	      }
//
//	    polyData->GetPointData()->AddArray(tubeRadius);
//	    polyData->GetPointData()->SetActiveScalars("TubeRadius");

//	    vtkSmartPointer<vtkXMLPolyDataWriter> tmpWriter = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
//	    tmpWriter->SetInputData(polyData);
//	    tmpWriter->SetFileName("tmpPolydata.vtp");
//	    tmpWriter->Write();



#if 1
	    //H-Tree (tubeFilter):
		vtkSmartPointer<vtkAppendPolyData> appendFilter = vtkSmartPointer<vtkAppendPolyData>::New();
//		for (int line_id = 0; line_id < n_branches-1; line_id++) { // -1, because of root branch
		for (int line_id = 0; line_id < n_blocks; line_id++) { // -1, because of root branch
			vtkSmartPointer<vtkLineSource> lineSource = vtkSmartPointer<vtkLineSource>::New();
			vtkSmartPointer<vtkTubeFilter> tubeFilter = vtkSmartPointer<vtkTubeFilter>::New();
			lineSource->SetPoint1(points2->GetPoint(nx[line_id]));
			lineSource->SetPoint2(points2->GetPoint(ny[line_id]));

			lineSource->Update();
			vtkSmartPointer<vtkPolyData> branchPolyData = lineSource->GetOutput();
			vtkSmartPointer<vtkDoubleArray> flowArray = vtkSmartPointer<vtkDoubleArray>::New();
			flowArray->InsertNextValue((flowVar[0]->GetValue(line_id)));
			flowArray->SetName("blood_flow2");

			branchPolyData->GetCellData()->AddArray(flowArray);


	        tubeFilter->SetInputData(branchPolyData);
//			tubeFilter->SetRadius(5e-6*(flowVar[0]->GetValue(line_id))+3e-5);
//			tubeFilter->SetRadius(8.8e-5*(stateVars[0]->GetValue(line_id))-4.3e-5);
			tubeFilter->SetRadius(0.0003*(stateVars[0]->GetValue(line_id))-0.00023);
			tubeFilter->SetNumberOfSides(20);
			tubeFilter->CappingOn();
			appendFilter->AddInputConnection(tubeFilter->GetOutputPort());
		}

		for (int line_id = n_blocks; line_id < n_branches-1; line_id++) { // -1, because of root branch
			vtkSmartPointer<vtkLineSource> lineSource = vtkSmartPointer<vtkLineSource>::New();
			vtkSmartPointer<vtkTubeFilter> tubeFilter = vtkSmartPointer<vtkTubeFilter>::New();
			lineSource->SetPoint1(points2->GetPoint(nx[line_id]));
			lineSource->SetPoint2(points2->GetPoint(ny[line_id]));

			lineSource->Update();
			vtkSmartPointer<vtkPolyData> branchPolyData = lineSource->GetOutput();
			vtkSmartPointer<vtkDoubleArray> flowArray = vtkSmartPointer<vtkDoubleArray>::New();
			flowArray->InsertNextValue((flowVar[0]->GetValue(line_id)));
			flowArray->SetName("blood_flow2");

			branchPolyData->GetCellData()->AddArray(flowArray);


	        tubeFilter->SetInputData(branchPolyData);
			tubeFilter->SetRadius(5e-5);  //TODO: each level include different radius! (5e-5*pow(2,0.5*level))
//			tubeFilter->SetRadius(5e-6*(stateVars[0]->GetValue(line_id))+3e-5);
			tubeFilter->SetNumberOfSides(20);
			tubeFilter->CappingOn();
			appendFilter->AddInputConnection(tubeFilter->GetOutputPort());
		}

		vtkSmartPointer<vtkLineSource> lineSource = vtkSmartPointer<vtkLineSource>::New();
		vtkSmartPointer<vtkTubeFilter> tubeFilter = vtkSmartPointer<vtkTubeFilter>::New();
		lineSource->SetPoint1(points2->GetPoint(ny[n_branches-1]));
		lineSource->SetPoint2(points2->GetPoint(ny[n_branches-1]+1));

		lineSource->Update();
		vtkSmartPointer<vtkPolyData> branchPolyData = lineSource->GetOutput();
		vtkSmartPointer<vtkDoubleArray> flowArray = vtkSmartPointer<vtkDoubleArray>::New();
//		flowArray->InsertNextValue(20);
		flowArray->InsertNextValue((flowVar[0]->GetValue(n_branches-1)));
		flowArray->SetName("blood_flow2");

		branchPolyData->GetCellData()->AddArray(flowArray);


        tubeFilter->SetInputData(branchPolyData);
		tubeFilter->SetRadius(5e-5);
//			tubeFilter->SetRadius(5e-6*(stateVars[0]->GetValue(line_id))+3e-5);
		tubeFilter->SetNumberOfSides(20);
		tubeFilter->SetCapping(3);
		appendFilter->AddInputConnection(tubeFilter->GetOutputPort());




//		//root branch:
//		vtkSmartPointer<vtkLineSource> lineSource = vtkSmartPointer<vtkLineSource>::New();
//		vtkSmartPointer<vtkTubeFilter> tubeFilter = vtkSmartPointer<vtkTubeFilter>::New();
//		lineSource->SetPoint1(points2->GetPoint(ny[n_branches-1]));
//		lineSource->SetPoint2(points2->GetPoint(ny[n_branches-1]+1));
//		lineSource->Update();
//
//		vtkSmartPointer<vtkPolyData> branchPolyData = lineSource->GetOutput();
//		vtkSmartPointer<vtkDoubleArray> flowArray = vtkSmartPointer<vtkDoubleArray>::New();
//		flowArray->InsertNextValue((flowVar[0]->GetValue(n_branches)));
//		branchPolyData->GetCellData()->AddArray(flowArray);
//
//		tubeFilter->SetInputConnection(branchPolyData);
//		tubeFilter->SetRadius(5e-6*(flowVar[0]->GetValue(n_branches-1))+3e-5);
//		tubeFilter->SetNumberOfSides(20);
//		appendFilter->AddInputConnection(tubeFilter->GetOutputPort());
#endif




	    //H-Tree (tubeFilter): - try
		vtkSmartPointer<vtkAppendPolyData> appendFilter2 = vtkSmartPointer<vtkAppendPolyData>::New();
		vtkSmartPointer<vtkTubeFilter> tube = vtkSmartPointer<vtkTubeFilter>::New();
		tube->SetInputData(polyData);
		tube->SetVaryRadiusToVaryRadiusByAbsoluteScalar();
//		tube->SetRadius(5e-5);
		tube->SetNumberOfSides(20);
		appendFilter2->AddInputConnection(tube->GetOutputPort());



// Write to files:
		// Tissue blocks:
		std::ostringstream filename_buffer;
		//filename_buffer << "/home/katharina/workspace_GitHub/parbrain/np02_nlev07_sbtr03_3/paraView_blocks" << i << ".vtu";
		filename_buffer << "/hpc/home/kdo40/Frontiers_in_Physiology/parbrain/ATP_input/np64_nlev13_sbtr03/paraView_blocks" << i << ".vtu";
		vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New(); // save uGrid to file
		writer->SetFileName(filename_buffer.str().c_str());
		writer->SetInputData(uGrid);
		writer->Write();

//		// H-Tree:
//		std::ostringstream filename_buffer2;
//		filename_buffer2 << "/home/katharina/workspace_GitHub/parbrain/np02_nlev07_sbtr03_3/paraView_Htree" << i << ".vtu";
//		vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer2 = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New(); // save uGrid2 to file
//		writer2->SetFileName(filename_buffer2.str().c_str());
//		writer2->SetInputData(uGrid2);
//		writer2->Write();

//		 H-Tree (tubeFilter):
		std::ostringstream filename_buffer2;
//		filename_buffer2 << "/home/katharina/workspace_GitHub/parbrain/np02_nlev07_sbtr03_3/paraView_Htree_tube" << i << ".vtp";
		filename_buffer2 << "/hpc/home/kdo40/Frontiers_in_Physiology/parbrain/ATP_input/np64_nlev13_sbtr03/paraView_Htree_tube" << i << ".vtp";
		vtkSmartPointer<vtkXMLPolyDataWriter> writer2 = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
		writer2->SetFileName(filename_buffer2.str().c_str());
		writer2->SetInputConnection(appendFilter->GetOutputPort());
		writer2->Write();

		// H-Tree (tubeFilter): - try
//		std::ostringstream filename_buffer3;
////		filename_buffer3 << "/home/katharina/workspace_GitHub/parbrain/np02_nlev07_sbtr03_3/paraView_Htree_tube_try" << i << ".vtp";
//		filename_buffer3 << "/home/katharina/power7/Frontiers_in_Physiology/parbrain/np04_nlev07_sbtr03_K/paraView_Htree_tube_try" << i << ".vtp";
//		vtkSmartPointer<vtkXMLPolyDataWriter> writer3 = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
//		writer3->SetFileName(filename_buffer3.str().c_str());
//		writer3->SetInputConnection(appendFilter2->GetOutputPort());
//		writer3->Write();

	} // time loop



	return EXIT_SUCCESS;
}
