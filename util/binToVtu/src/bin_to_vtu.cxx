/*
Converts the data from the simulate script into vtu files for Paraview
*/

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

#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkXMLPolyDataWriter.h>

#include <vtkLine.h>
#include <vtkTubeFilter.h>
#include <vtkAppendPolyData.h>
#include <vtkLineSource.h>

#define POW_OF_2(x) (1 << (x)) // macro for 2^x using bitwise shift

int main(int argc, char *argv[]) {

	// General parameters:
	#define BLOCK_LENGTH 4e-4	// Length of one tissue block [m].
	char Prefix[] = "";

//*****************************************************************************************
// run with no arguments for debugging
#if 0
	char *dirName="np01_nlev03_sbtr01";
	int tf = 10;
#else
	if (argc != 3)
	{
		printf("Uh oh, spaghettio. You have not entered the correct number of arguments.\n");
		std::cerr << "Usage: " << argv[0] << " <Data directory> <Final time>\n";
		exit(EXIT_FAILURE);
	}

	char *dirName=argv[1];       	// First argument: Folder name.
	int tf = atoi(argv[2]); 		// Second argument: Final time (atoi: str -> int).
#endif
//******************************************************************************************

	// Read configuration file:
	char iSuffix[] = "/info.dat";
	char iOutfile[128];

	sprintf(iOutfile, "%s%s%s", Prefix, dirName, iSuffix);
	std::ifstream conf_file(iOutfile, std::ifstream::binary);

	if (!conf_file)
	{
		std::cerr << "Cannot read " << iOutfile << " file." << std::endl;
		return EXIT_FAILURE;
	}

	std::string header;
	std::getline(conf_file, header); // Skip header line.

	int conf_array[7];  // Create temporary array that stores parameters from configuration file (info.dat).
	int b;
	int i = 0;

	while (conf_file >> b)
	{
	    conf_array[i] = b;
	    i++;
	}

	int n_procs  = conf_array[0];
	int n_blocks_per_rank = conf_array[1];
	int n_state_vars = conf_array[2];
	int m_local = conf_array[3];
	int n_local = conf_array[4];
	int m_global = conf_array[5];
	int n_global = conf_array[6];
	int n_blocks = n_procs * n_blocks_per_rank;
	int n_cols = n_local * n_global;  				// Number of columns of tissue blocks (j).
	int n_rows = m_local * m_global;  				// Number of rows of tissue blocks (i).

	// Read tissue block state binary file:
	char tbSuffix[] = "/tissueBlocks.dat";
	char tbOutfile[128];

	sprintf(tbOutfile, "%s%s%s", Prefix, dirName, tbSuffix);
	std::ifstream is_tissue(tbOutfile, std::ifstream::binary);

	if (!is_tissue)
	{
		std::cerr << "Cannot read " << tbOutfile << " file." << std::endl;
		return EXIT_FAILURE;
	}

	// Read flow binary file:
	char fSuffix[] = "/flow.dat";
	char fOutfile[128];

	sprintf(fOutfile, "%s%s%s", Prefix, dirName, fSuffix);
	std::ifstream is_flow(fOutfile, std::ifstream::binary);

	if (!is_flow)
	{
		std::cerr << "Cannot read " << fOutfile << " file." << std::endl;
		return EXIT_FAILURE;
	}

	// Read pressure binary file:
	char pSuffix[] = "/pressure.dat";
	char pOutfile[128];

	sprintf(pOutfile, "%s%s%s", Prefix, dirName, pSuffix);
	std::ifstream is_pressure(pOutfile, std::ifstream::binary);

	if (!is_pressure)
	{
		std::cerr << "Cannot read " << pOutfile << " file." << std::endl;
		return EXIT_FAILURE;
	}

    // Read x and y coordinates from tissue_block binary file:
	double xCoord[n_blocks];
	double yCoord[n_blocks];

	is_tissue.read((char *) xCoord, sizeof(xCoord));
	is_tissue.read((char *) yCoord, sizeof(yCoord));

	for (int i = 0; i < n_blocks; i++)
	{
	    std::cout << "Block: " << i << "\t x-coord: " << xCoord[i] << "  \t y-coord: " << yCoord[i] << std::endl; // Print all coordinates.
	}

// 1. Create a vtkPoints object and store the points in it:

	// 1.1 Tissue blocks:
	vtkSmartPointer<vtkPoints> points_tb = vtkSmartPointer<vtkPoints>::New();

	double points_xorigin = xCoord[0] - (BLOCK_LENGTH / 2); // Origin - x coordinate.
	double points_yorigin = yCoord[0] - (BLOCK_LENGTH / 2); // Origin - y coordinate.

	// Lower "layer" of points: adds points (corners of 3D blocks)
	for (int j = 0; j <= n_cols; j++)
	{
		double points_xcoord = points_xorigin + j * BLOCK_LENGTH;

		for (int i = 0; i <= n_rows; i++)
		{
			double points_ycoord = points_yorigin + i * BLOCK_LENGTH;
			points_tb->InsertNextPoint(points_xcoord, points_ycoord, 0); // x,y,z coordinates
		}
	}

	// Upper "layer" of points:
	for (int j = 0; j <= n_cols; j++)
	{
		double points_xcoord = points_xorigin + j * BLOCK_LENGTH;

		for (int i = 0; i <= n_rows; i++)
		{
			double points_ycoord = points_yorigin + i * BLOCK_LENGTH;
			points_tb->InsertNextPoint(points_xcoord, points_ycoord, BLOCK_LENGTH);
		}
	}

	// 1.2 H-Tree:
	vtkSmartPointer<vtkPoints> points_tree = vtkSmartPointer<vtkPoints>::New();
    int n_rows_h = n_rows;                 		// Row variable - gets updated.
    int n_cols_h = n_cols;                 		// Column variable - gets updated.
	int n_levels = log2(n_blocks) + 1;       	// Number of bifurcations.
	int n_branches = POW_OF_2(n_levels) - 1;    	// Number of branches. 							*** 1 << n_bifr = 2^n_bifr ***
	int n_nodes = POW_OF_2(n_levels - 1) - 1; 	// Number of nodes. (?) - Leaf nodes excluded!

    double xpoints_tree[n_branches];    // x coordinates for points for tree - branches & pressure.
    double ypoints_tree[n_branches];    // y coordinates for points for tree - branches & pressure.

    // Adjacency matrix (algorithm from adjacency.c) - TODO: Don't repeat logic!
    int a, k1, k2;
    int row = 0;
    int col = POW_OF_2(n_levels-1);
    int xbranch = 0;
	int offset = n_nodes + 1; // +1 because we need to add the leaf nodes.

	int nx[n_branches], ny[n_branches], nz1[n_nodes], nz2[n_nodes], h1[n_nodes], h2[n_nodes];

	for (int i = 0; i < col; i++)
	{
		ny[i] = i;
	}

	// L loop: from bottom level up to the top of the tree.
	for (int L = n_levels - 1; L > 0; L--)
	{
		a = POW_OF_2(n_levels) - POW_OF_2(L+1);

		if (xbranch)
		{
			for (int j = 0; j < n_cols_h; j+=2)
			{
				for (int i = 0; i < n_rows_h; i++)
				{
					k1 = a + i + j*n_rows_h;
					k2 = a + i + (j+1)*n_rows_h;
					nx[k1] = row + n_blocks;
					nx[k2] = row + n_blocks;
					ny[col] = row + n_blocks;
					h1[row] = k1;
					h2[row] = k2;
					row++;
					col++;
				}
			}

			n_cols_h /= 2; // Every 2nd level no of columns gets halved.
		}
		else
		{
			for (int j = 0; j < n_cols_h; j++)
			{
				for (int i = 0; i < n_rows_h; i+=2)
				{
					k1 = a + i + j*n_rows_h;
					k2 = k1 + 1;
					nx[k1] = row + n_blocks;
					nx[k2] = row + n_blocks;
					ny[col] = row + n_blocks;
					h1[row] = k1;
					h2[row] = k2;
					row++;
					col++;
				}
			}

			n_rows_h /= 2; // Every 2nd level no of rows gets halved.
		}
		xbranch = !xbranch; // switch xbranch 0 <--> 1 each level
	}

	// Insert points for leaf nodes:
	for (int i = 0; i < n_blocks; i++)
	{
		xpoints_tree[i] = xCoord[i];
		ypoints_tree[i] = yCoord[i];
		points_tree->InsertNextPoint(xCoord[i], yCoord[i], BLOCK_LENGTH);
	}

	// Insert points for all other nodes (not leaf nodes):
	for (int i = 0; i < n_nodes; i++)
	{
		nz1[i] = ny[h1[i]];
		nz2[i] = ny[h2[i]];
	}

	for (int i = n_blocks; i < n_branches; i++)
	{
		xpoints_tree[i] = (xpoints_tree[nz1[i - n_blocks]] + xpoints_tree[nz2[i - n_blocks]]) / 2;
		ypoints_tree[i] = (ypoints_tree[nz1[i - n_blocks]] + ypoints_tree[nz2[i - n_blocks]]) / 2;
		points_tree->InsertNextPoint((xpoints_tree[nz1[i - n_blocks]] + xpoints_tree[nz2[i - n_blocks]]) / 2 , (ypoints_tree[nz1[i - n_blocks]] + ypoints_tree[nz2[i - n_blocks]]) / 2 , BLOCK_LENGTH);
	}

	// Insert point for root branch at (x,y)=(0,0):
	points_tree->InsertNextPoint(0, 0, n_cols*BLOCK_LENGTH);

// 2. Create cell array:

    // 2.1 Tissue blocks:
	vtkSmartPointer<vtkCellArray> cellArray_tb = vtkSmartPointer<vtkCellArray>::New();  // Create a cell array to store hexahedrons in.
	vtkSmartPointer<vtkHexahedron> hexahedron = vtkSmartPointer<vtkHexahedron>::New(); 	// Create hexahedrons (3D blocks).

	double npoints_offset = (n_rows + 1) * (n_cols + 1); // number of points - offset for 2nd "layer"

	for (int block_id = 0; block_id < n_blocks; block_id++)
	{
		int col = block_id / n_rows;     						// which column are we in?
		int i = block_id + col;
		int j = i + n_rows + 1;
		int k = j + 1;
		int l = k - (n_rows + 1);

		hexahedron->GetPointIds()->SetId(0, i);					// lower SW corner of hexahedron
		hexahedron->GetPointIds()->SetId(1, j); 				// lower SE corner of hexahedron
		hexahedron->GetPointIds()->SetId(2, k);					// lower NE corner of hexahedron
		hexahedron->GetPointIds()->SetId(3, l); 				// lower NW corner of hexahedron
		hexahedron->GetPointIds()->SetId(4, i+npoints_offset); 	// upper SW corner of hexahedron
		hexahedron->GetPointIds()->SetId(5, j+npoints_offset); 	// upper SE corner of hexahedron
		hexahedron->GetPointIds()->SetId(6, k+npoints_offset);  // upper NE corner of hexahedron
		hexahedron->GetPointIds()->SetId(7, l+npoints_offset); 	// upper NW corner of hexahedron

		cellArray_tb->InsertNextCell(hexahedron); // Add hexahedrons to cellArray
	}

    // 2.2 H-Tree:
	vtkSmartPointer<vtkLine> lines = vtkSmartPointer<vtkLine>::New();
	vtkSmartPointer<vtkCellArray> cellArray_tree = vtkSmartPointer<vtkCellArray>::New();

	for (int line_id = 0; line_id < n_branches-1; line_id++)  // -1, because of root branch
	{
		lines->GetPointIds()->SetId(0, nx[line_id]); \
		lines->GetPointIds()->SetId(1, ny[line_id]);
		cellArray_tree->InsertNextCell(lines);
	}

	lines->GetPointIds()->SetId(0, ny[n_branches-1]); // root branch
	lines->GetPointIds()->SetId(1, ny[n_branches-1]+1);
	cellArray_tree->InsertNextCell(lines); // Create a cell array to store hexahedrons in and add them to it

// 3. Create unstructured grid

	// 3.1 Tissue blocks:
	vtkSmartPointer<vtkUnstructuredGrid> uGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
	uGrid->SetPoints(points_tb);
	uGrid->SetCells(VTK_HEXAHEDRON, cellArray_tb);

	// 3.3 H-Tree_lines:
	vtkSmartPointer<vtkUnstructuredGrid> uGrid2 = vtkSmartPointer<vtkUnstructuredGrid>::New();
	uGrid2->SetPoints(points_tree);
	uGrid2->SetCells(VTK_LINE, cellArray_tree); 		// uGrid2->SetLines(cellArray2);

// 4. Add binary data as attributes to cells:

	//char *var_names[] = {"radius_coupled","radius_decoupled","R_k","N_Na_k","N_K_k","N_HCO3_k","N_Cl_k","N_Na_s","N_K_s","N_HCO3_s","K_p","w_k","ca_i","ca_sr_i","v_i","w_i","ip3_i","K_i","ca_j","ca_er_j","v_j","ip3_j","Mp","AMp","AM","input_PLC","input_K_df","input_K_flux","NOi","NOj","NOn","cGMP","eNOS","nNOS","ca_n","E_b","E_6c","E_5c"};
//	char const *var_names[] = {"radius_coupled","R_k","N_Na_k","N_K_k","N_HCO3_k","N_Cl_k","N_Na_s","N_K_s","N_HCO3_s","K_p","w_k","ca_i","ca_sr_i","v_i","w_i","ip3_i","K_i","ca_j","ca_er_j","v_j","ip3_j","Mp","AMp","AM","PLC_input","K_input","flux_ft","NOi","NOj","NOn","cGMP","eNOS","nNOS","ca_n","E_b","E_6c","E_5c"};
	//	char const *var_names[] = {"radius_coupled","R_k","N_Na_k","N_K_k","N_HCO3_k","N_Cl_k","N_Na_s","N_K_s","N_HCO3_s","K_p","w_k","ca_i","ca_sr_i","v_i","w_i","ip3_i","K_i","ca_j","ca_er_j","v_j","ip3_j","Mp","AMp","AM","PLC_input","K_input","flux_ft"};
	char const *var_names[] = {"radius_coupled","R_k","N_Na_k","N_K_k","N_HCO3_k","N_Cl_k","N_Na_s","N_K_s","N_HCO3_s","K_p","w_k","ca_i","ca_sr_i","v_i","w_i","ip3_i","K_i","ca_j","ca_er_j","v_j","ip3_j","Mp","AMp","AM","K_e","PLC_input","K_input","flux_ft"};

	// 4.1 Time step loop:
	double time_tb, time_tree;

	for (int i = 0; i <= tf; i++)
	{
		is_tissue.read((char *) &time_tb, sizeof(time_tb)); // Read time from both binary files (is not used for anything at the moment...)
		is_flow.read((char *) &time_tree, sizeof(time_tree));

		std::cout << "Time: " << time_tb << " \t" << time_tree << std::endl; // Print time from both files for comparison.

		// 4.1.1 Tissue blocks:
		std::vector<vtkSmartPointer<vtkDoubleArray> > stateVars; //Create vector of arrays for each state variable.

		// Set state variable names.
		for (int v = 0; v < n_state_vars; v++)
		{
			vtkSmartPointer<vtkDoubleArray> array_tb = vtkSmartPointer<vtkDoubleArray>::New(); 	// Create array.
			array_tb->SetName(var_names[v]);
			stateVars.push_back(array_tb);
		}

		// 4.1.2 H-Tree:
		std::vector<vtkSmartPointer<vtkDoubleArray> > flowVar; 
		vtkSmartPointer<vtkDoubleArray> array_tree = vtkSmartPointer<vtkDoubleArray>::New(); // Create array.
		array_tree->SetName("blood_flow"); // Only one array for flow variables.
		flowVar.push_back(array_tree);

		// 4.1.3 H-Tree - unscaled
		std::vector<vtkSmartPointer<vtkDoubleArray> > flowVar_unscaled; 	
		vtkSmartPointer<vtkDoubleArray> array_tree_unscaled = vtkSmartPointer<vtkDoubleArray>::New(); // Create array.
		array_tree_unscaled->SetName("blood_flow_unscaled"); // Only one array for flow variables.
		flowVar_unscaled.push_back(array_tree_unscaled);


	// 4.2 Read state variables and add as attributes to uGrid.

		// 4.2.1 Tissue blocks:
		for (int k = 0; k < n_blocks; k++)
		{
			double temp_array_tb[n_state_vars];
			is_tissue.read((char *) temp_array_tb, sizeof(temp_array_tb));  // Read state variables from binary file.

			for (int v = 0; v < n_state_vars; v++)
			{
				stateVars[v]->InsertNextValue(temp_array_tb[v]);
			}
		}

		for (int v = 0; v < n_state_vars; v++)
		{
			uGrid->GetCellData()->AddArray(stateVars[v]);
		}

		// 4.2.2 H-Tree:
		for (int level = 0; level < n_levels; level++)
		{
			// Scale flow variables: Each level times by 2^(nlevels-level-1)
			double n_lines = pow(2,(n_levels-level-1));
			for(int i = 0; i < n_lines; i++)
			{
				double temp_array_tree[1];
				is_flow.read((char *) temp_array_tree, sizeof(temp_array_tree));  // Read flow variables from binary file.
				flowVar[0]->InsertNextValue( temp_array_tree[0] * pow(2, ((n_levels - level - 1) )) );
				flowVar_unscaled[0]->InsertNextValue( temp_array_tree[0] );
			}
		}

		vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
	    polyData->SetPoints(points_tree);
	    polyData->SetLines(cellArray_tree);
	    polyData->GetCellData()->AddArray(flowVar[0]);
	    polyData->GetCellData()->AddArray(flowVar_unscaled[0]);

		// 4.2.3 H-Tree_lines:
		uGrid2->GetCellData()->AddArray(flowVar[0]);
		uGrid2->GetCellData()->AddArray(flowVar_unscaled[0]);


	    // H-Tree (tubeFilter):
	    // a) Leaf branches:
		vtkSmartPointer<vtkAppendPolyData> appendFilter = vtkSmartPointer<vtkAppendPolyData>::New();

		for (int line_id = 0; line_id < n_blocks; line_id++)
		{
			vtkSmartPointer<vtkLineSource> lineSource = vtkSmartPointer<vtkLineSource>::New();
			vtkSmartPointer<vtkTubeFilter> tubeFilter = vtkSmartPointer<vtkTubeFilter>::New();
			lineSource->SetPoint1(points_tree->GetPoint(nx[line_id]));
			lineSource->SetPoint2(points_tree->GetPoint(ny[line_id]));

			lineSource->Update();

			vtkSmartPointer<vtkPolyData> branchPolyData = lineSource->GetOutput();
			vtkSmartPointer<vtkDoubleArray> flowArray = vtkSmartPointer<vtkDoubleArray>::New();
			flowArray->InsertNextValue((flowVar[0]->GetValue(line_id)));
			flowArray->SetName("blood_flow");

			branchPolyData->GetCellData()->AddArray(flowArray);

	        tubeFilter->SetInputData(branchPolyData);
			tubeFilter->SetRadius(0.0003*(stateVars[0]->GetValue(line_id))-0.00023);  	// varying radii
			tubeFilter->SetNumberOfSides(20);
			tubeFilter->CappingOn();  												  	// Cap tubes.
			appendFilter->AddInputConnection(tubeFilter->GetOutputPort());
		}

		// b) Branches further up the tree:
		for (int line_id = n_blocks; line_id < n_branches-1; line_id++) // -1, because of root branch
		{
			vtkSmartPointer<vtkLineSource> lineSource = vtkSmartPointer<vtkLineSource>::New();
			vtkSmartPointer<vtkTubeFilter> tubeFilter = vtkSmartPointer<vtkTubeFilter>::New();
			lineSource->SetPoint1(points_tree->GetPoint(nx[line_id]));
			lineSource->SetPoint2(points_tree->GetPoint(ny[line_id]));

			lineSource->Update();

			vtkSmartPointer<vtkPolyData> branchPolyData = lineSource->GetOutput();
			vtkSmartPointer<vtkDoubleArray> flowArray = vtkSmartPointer<vtkDoubleArray>::New();
			flowArray->InsertNextValue((flowVar[0]->GetValue(line_id)));
			flowArray->SetName("blood_flow");

			branchPolyData->GetCellData()->AddArray(flowArray);

	        tubeFilter->SetInputData(branchPolyData);
			tubeFilter->SetRadius(5e-5);  //TODO: each level include different radius! (5e-5*pow(2,0.5*level))
//			tubeFilter->SetRadius(5e-6*(stateVars[0]->GetValue(line_id))+3e-5);
			tubeFilter->SetNumberOfSides(20);
			tubeFilter->CappingOn();
			appendFilter->AddInputConnection(tubeFilter->GetOutputPort());
		}

		// c) Root branch:
		vtkSmartPointer<vtkLineSource> lineSource = vtkSmartPointer<vtkLineSource>::New();
		vtkSmartPointer<vtkTubeFilter> tubeFilter = vtkSmartPointer<vtkTubeFilter>::New();
		lineSource->SetPoint1(points_tree->GetPoint(ny[n_branches-1]));
		lineSource->SetPoint2(points_tree->GetPoint(ny[n_branches-1]+1));

		lineSource->Update();

		vtkSmartPointer<vtkPolyData> branchPolyData = lineSource->GetOutput();
		vtkSmartPointer<vtkDoubleArray> flowArray = vtkSmartPointer<vtkDoubleArray>::New();
		flowArray->InsertNextValue((flowVar[0]->GetValue(n_branches-1)));
		flowArray->SetName("blood_flow");

		branchPolyData->GetCellData()->AddArray(flowArray);

        tubeFilter->SetInputData(branchPolyData);
		tubeFilter->SetRadius(5e-5);
		tubeFilter->SetNumberOfSides(20);
		tubeFilter->SetCapping(3);
		appendFilter->AddInputConnection(tubeFilter->GetOutputPort());

// 5. Write to files:
		// 5.1 Tissue blocks:
		std::ostringstream filename_buffer_tb;
		char tbVtuSuffix[] = "/paraView_blocks";
		char tbVtuOutfile[128];

		sprintf(tbVtuOutfile, "%s%s%s",Prefix,dirName,tbVtuSuffix);

		filename_buffer_tb << tbVtuOutfile << i << ".vtu";

		vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer_tb = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New(); // save uGrid to file
		writer_tb->SetFileName(filename_buffer_tb.str().c_str());
		writer_tb->SetInputData(uGrid);
		writer_tb->Write();

		// 5.2 H-Tree (tubeFilter):
		std::ostringstream filename_buffer_tree;
		char fVtuSuffix[] = "/paraView_Htree_tubes";
		char fVtuOutfile[128];

		sprintf(fVtuOutfile, "%s%s%s",Prefix,dirName,fVtuSuffix);

		filename_buffer_tree << fVtuOutfile << i << ".vtp";

		vtkSmartPointer<vtkXMLPolyDataWriter> writer_tree = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
		writer_tree->SetFileName(filename_buffer_tree.str().c_str());
		writer_tree->SetInputConnection(appendFilter->GetOutputPort());
		writer_tree->Write();

		// 5.3 H-Tree_lines:
		std::ostringstream filename_buffer_tree2;
		char fVtuSuffix2[] = "/paraView_Htree_lines";
		char fVtuOutfile2[128];

		sprintf(fVtuOutfile2, "%s%s%s",Prefix,dirName,fVtuSuffix2);

		filename_buffer_tree2 << fVtuOutfile2 << i << ".vtu";

		vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer2 = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New(); // save uGrid2 to file
		writer2->SetFileName(filename_buffer_tree2.str().c_str());
		writer2->SetInputData(uGrid2);
		writer2->Write();

	}

	return EXIT_SUCCESS;
}