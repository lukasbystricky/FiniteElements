#include "myfunctions.h"

/***************************************************************************************************
 * Steady state Navier-Stokes flow solver
 * 
 * Author: Lukas Bystricky, lb13f@my.fsu.edu
 * 
 * Solves steady state Navier-Stokes equations
 * 
 * Supports adaptive meshing and convergence studies.
 * 
 * Boundary conditions, and exact solution are specified in "myfunctions.h". 
 * Exact solution is used only if errors are to be calculated.
 **************************************************************************************************/
template<int dim>
class NavierStokesSolver
{
	public:
		NavierStokesSolver();
		void run();
			
		Function<dim> *forcing_function;
        Function<dim> *boundary_values;
        Function<dim> *exact_solution;
		
	private:
        void read_inputs();
		void setup_geometry(int cycle);
		void assemble_system(double nu);
		void assemble_stokes_system(double nu);
		void run_newton_loop(int cycle);
		void solve();
		void refine_grid();
		void output_results(int cycle);
        void calculate_error(int cycle);  
        void print_errors();      
        
		struct InputFile
        {
            int          nx, ny, tpa, mina, maxa, refinement_steps;
            double       nu, r_in, r_out, domain_length, domain_height;
            bool         print_error, verbose, slit_flow;
            std::string  output_name, refinement_type, gmesh_file;
        };
        
        ParameterHandler 	   prm;   
        InputFile              input;        
		Triangulation<dim>     mesh;
		FESystem<dim>          fe;
		DoFHandler<dim>        dof_handler;
		
		ConstraintMatrix       constraints; 
		SparsityPattern        sparsity_pattern;
		SparseMatrix<double>   system_matrix;
		Vector<double>         system_rhs;
		Vector<double>         solution;
		Vector<double> 		   previous_newton_step;

		ConvergenceTable	   velocity_convergence_table;
		ConvergenceTable	   pressure_convergence_table;
};

template<int dim>
NavierStokesSolver<dim>::NavierStokesSolver():
	/***********************************************************************************************
	 *  taylor hood element, quadratic basis functions for 2 velocity components, 
	 * linear basis functions for pressure
	 **********************************************************************************************/
    fe(FE_Q<dim>(2), dim, FE_Q<dim>(1), 1),
    dof_handler(mesh)
{}

/***************************************************************************************************
 * Read and store inputs from user parameter file "input.in"
 **************************************************************************************************/
template<int dim>
void NavierStokesSolver<dim>::read_inputs()
{
	/***********************************************************************************************
	 * Set up paramter file patterns
	 **********************************************************************************************/
	prm.enter_subsection("Fluid constants");
	{
		prm.declare_entry("nu", "1", Patterns::Double(1e-7, 1e4), 
					"kinematic viscosity");
	}
	prm.leave_subsection();
	
	prm.enter_subsection("Problem setup");
	{
		prm.declare_entry("mesh file", "auto", Patterns::Anything(),
					"gmesh file, if none availiable enter 'auto'");
					
		prm.declare_entry("length", "1", Patterns::Double(1e-7, 1e7),
					"length of rectangular domain");
		
		prm.declare_entry("height", "1", Patterns::Double(1e-7, 1e7),
					"height of rectangular domain");
								
		prm.declare_entry("nx", "9", Patterns::Double(5, 10000), 
					"number of nodes in x direction");
					
		prm.declare_entry("ny", "9", Patterns::Double(5, 10000), 
					"number of nodes in y direction");

		prm.declare_entry("r inlet", "-0.980655e5", Patterns::Double(-1e7, 1e7),
					"normal stress at inlet (specify if pressure induced)");
	
		prm.declare_entry("r outlet", "-9.8e4", Patterns::Double(-1e7, 1e7),
					"normal stress at outlet (specify if pressure induced)");
	}
	prm.leave_subsection();
	
	prm.enter_subsection("Run options");
	{	
		prm.declare_entry("slit flow", "true", Patterns::Bool(),
						"slit flow model, (true or false)");
						
		prm.declare_entry("refinement type", "adaptive", Patterns::Anything(),
						"refienment type (adaptive or uniform)");
						
		prm.declare_entry("refinement steps", "5", Patterns::Integer(0,20),
						"refienment steps");
		
		prm.declare_entry("minimum refinement level", "0", Patterns::Integer(0, 20),
						"minimum cell refinement level");
		
		prm.declare_entry("maximum refinement level", "3", Patterns::Integer(0, 20),
						"maximum cell refinement level");
	}
	prm.leave_subsection();  
													
	prm.enter_subsection("Output options");
	{
	  prm.declare_entry("output name", "stokes_flow", Patterns::Anything(), 
						"name of the output file (without extension)");
	  
	  prm.declare_entry("print errors", "true", Patterns::Bool(),
						"print errors to console. If exact solution known and in myfunctions.h, "
						"print L2 and H1 errors. For uniform mesh refinement this is essentially "
						"a convergence study. For adaptive mesh refinement, errors are calcualated "
						"after each adpat.");
						
	  prm.declare_entry("verbose", "true", Patterns::Bool(),
						"verbose output to console during solve");
	}
	prm.leave_subsection();	
	
	/***********************************************************************************************
	 * Read and store paramters
	 **********************************************************************************************/
	 prm.read_input("input.in");
	
	 
	 prm.enter_subsection("Fluid constants");
	 input.nu = prm.get_double("nu");
	 prm.leave_subsection();
	 
	 prm.enter_subsection("Problem setup");
	 
	 input.gmesh_file = prm.get("mesh file");
	 input.nx = prm.get_integer("nx");
	 input.ny = prm.get_integer("ny");
	 input.r_in = prm.get_double("r inlet");
	 input.r_out = prm.get_double("r outlet");
	 input.domain_length = prm.get_double("length");
	 input.domain_height = prm.get_double("height");
	 
	 prm.leave_subsection();
	 
	 prm.enter_subsection("Run options");
	 input.slit_flow = prm.get_bool("slit flow");
	 input.refinement_type = prm.get("refinement type");
	 input.refinement_steps = prm.get_double("refinement steps");
	 input.mina = prm.get_double("minimum refinement level");
	 input.maxa = prm.get_double("maximum refinement level");
	 prm.leave_subsection();
	 
	 prm.enter_subsection("Output options");
	 input.output_name = prm.get("output name");
	 input.print_error = prm.get_bool("print errors");
	 input.verbose =  prm.get_bool("verbose");
	 prm.leave_subsection(); 
	 
	 if (input.verbose)
	 {
		prm.print_parameters (std::cout, ParameterHandler::Text);
	 }
}

template<int dim>
void NavierStokesSolver<dim>::setup_geometry (int cycle)
{  		
	if (cycle == 0)
	{
		if (input.gmesh_file.compare("auto") != 0)
		{
			mesh.clear();
			GridIn<dim> grid_in;
			grid_in.attach_triangulation(mesh);
			std::ifstream input_file(input.gmesh_file.c_str());
			grid_in.read_msh(input_file); 
		}
		else
		{
			std::vector<unsigned int> number_elements(2);
			number_elements[0] = input.nx-1;
			number_elements[1] = input.ny-1;
	  
			Point<dim> bottom_left, top_right;
			bottom_left[0] = 0;
			bottom_left[1] = 0;
			top_right[0] = input.domain_length;
			top_right[1] = input.domain_height;
			
			GridGenerator::subdivided_hyper_rectangle(mesh, number_elements,
				bottom_left, top_right, false);
		}
	}
  
    dof_handler.distribute_dofs(fe);  
    
    if (input.verbose)
	{
		std::printf("Number of active cells:%d\n", mesh.n_active_cells());
		std::printf("Number of degrees of freedom:%d\n", dof_handler.n_dofs()); 
	}
		
	FEValuesExtractors::Vector velocities(0);
    FEValuesExtractors::Scalar pressure(dim);
    
    typename DoFHandler<dim>::active_cell_iterator  
			cell = dof_handler.begin_active(), endc = dof_handler.end();
			
	std::vector<bool> boundary_dofs(dof_handler.n_dofs(), false);
	
    constraints.clear(); 
				
    DoFTools::make_hanging_node_constraints (dof_handler, constraints);	    
    
    //set boundary labels:
    // (1) : inlet
    // (2) : outlet
    // (3) : everything else
	for (; cell!=endc; ++cell)
	{
		for (int f=0; f<GeometryInfo<dim>::faces_per_cell; f++)
		{
			if (cell->face(f)->at_boundary())
			{
				if (fabs(cell->face(f)->center()[0]) < 1e-8)
				{
					cell->face(f)->set_boundary_indicator(1);
				}
				else if (fabs(cell->face(f)->center()[0] - input.domain_length) < 1e-8)	
				{
					cell->face(f)->set_boundary_indicator(2);
				}
				else
				{
					cell->face(f)->set_boundary_indicator(3);
				}
			}
		}
	}
		
	if (input.slit_flow)
	{
		//apply velocity boundary conditions to top and bottom only
		//set tangential (vertical) velocity at inlet and outlet
		
		FEValuesExtractors::Scalar velocity_y(1);
		
		VectorTools::interpolate_boundary_values(dof_handler, 3, 
				*boundary_values, constraints, fe.component_mask(velocities));
						
		VectorTools::interpolate_boundary_values(dof_handler, 1, 
				*boundary_values, constraints, fe.component_mask(velocity_y));
		
		VectorTools::interpolate_boundary_values(dof_handler, 2, 
				*boundary_values, constraints, fe.component_mask(velocity_y));
		
	}
	else
	{
		//set Dirichlet boundary conditions on the velocity everywhere
		VectorTools::interpolate_boundary_values(dof_handler, 1, 
				*boundary_values, constraints, fe.component_mask(velocities));
		
		VectorTools::interpolate_boundary_values(dof_handler, 2, 
				*boundary_values, constraints, fe.component_mask(velocities));
		
		VectorTools::interpolate_boundary_values(dof_handler, 3, 
				*boundary_values, constraints, fe.component_mask(velocities));
					
		//constrain first pressure dof to be 0
		DoFTools::extract_boundary_dofs(dof_handler, fe.component_mask(pressure), 
					boundary_dofs);
					
		const unsigned int first_boundary_dof = std::distance(boundary_dofs.begin(),
			std::find (boundary_dofs.begin(), boundary_dofs.end(), true));
			
		constraints.add_line(first_boundary_dof);
	}
    
    constraints.close();
          
	CompressedSparsityPattern c_sparsity(dof_handler.n_dofs(), dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern (dof_handler, c_sparsity);
    
    constraints.condense(c_sparsity);
	sparsity_pattern.copy_from(c_sparsity);
    
    system_matrix.reinit(sparsity_pattern);
    solution.reinit(dof_handler.n_dofs());     
    system_rhs.reinit(dof_handler.n_dofs());
    previous_newton_step.reinit(dof_handler.n_dofs());      
}

template<int dim>
void NavierStokesSolver<dim>::assemble_system(double nu)
{
	if (input.verbose)
	{
		printf("Assembling matrix...");
	}
	
	system_matrix.reinit(sparsity_pattern);    
    system_rhs.reinit(dof_handler.n_dofs());
    
	Timer timer;
	timer.start ();
	
	QGauss<dim>                 quadrature_formula(fe.degree + 1);
	QGauss<dim-1>               face_quadrature_formula(3);
	
	const int                   dofs_per_cell = fe.dofs_per_cell;
    const int                   n_q_points = quadrature_formula.size();	
    const int                   n_q_points_face = face_quadrature_formula.size();
    
    std::vector<Tensor<1,dim> > previous_newton_velocity_values(n_q_points);
	std::vector<Tensor< 2, dim> > previous_newton_velocity_gradients(n_q_points);
    
	std::vector<Vector<double> > rhs_values (n_q_points, Vector<double>(dim+1));	
	std::vector<Tensor<2,dim> > grad_phi_u (dofs_per_cell);
	std::vector<double>         div_phi_u (dofs_per_cell);
	std::vector<double>         phi_p (dofs_per_cell);
	std::vector<Tensor<1,dim> > phi_u (dofs_per_cell);  
	
	Vector<double>              cell_rhs(dofs_per_cell);
	FullMatrix<double>          cell_matrix (dofs_per_cell, dofs_per_cell);
    
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);    
	const FEValuesExtractors::Vector velocities (0);
	const FEValuesExtractors::Scalar pressure (dim);
	
    FEValues<dim> fe_values(fe, quadrature_formula, 
            update_values | update_gradients | update_JxW_values | update_quadrature_points);    
    FEFaceValues<dim> fe_face_values (fe, face_quadrature_formula, update_values 
			| update_quadrature_points | update_gradients | update_JxW_values);
			   
    typename DoFHandler<dim>::active_cell_iterator  
			cell = dof_handler.begin_active(), endc = dof_handler.end();
    
    for (; cell!=endc; ++cell)
    {
        fe_values.reinit(cell);
        cell_matrix = 0;
        cell_rhs = 0;

		/** Calculate velocity values and gradients from previous newton iteration
		* at each quadrature point in cell*********************************/
		fe_values[velocities].get_function_values(previous_newton_step, previous_newton_velocity_values);        
		fe_values[velocities].get_function_gradients(previous_newton_step, previous_newton_velocity_gradients);
		
		forcing_function->vector_value_list(fe_values.get_quadrature_points(), rhs_values);
		
        //calculate cell contribution to system         
        for (int q = 0; q < n_q_points; q++)
        {
			for (int k=0; k<dofs_per_cell; k++)
			{
				grad_phi_u[k] = fe_values[velocities].gradient (k, q);
				div_phi_u[k]  = fe_values[velocities].divergence (k, q);
				phi_p[k]      = fe_values[pressure].value (k, q);
				phi_u[k]      = fe_values[velocities].value (k, q);
			}
			
			for (int i = 0; i < dofs_per_cell; i++)
			{
				for (int j = 0; j < dofs_per_cell; j++)
				{
					if (input.slit_flow)
					{
						cell_matrix(i,j) += 
								((0.5*nu)*
								double_contract(grad_phi_u[i] + transpose(grad_phi_u[i]),
													grad_phi_u[j] + transpose(grad_phi_u[j]))
										+ phi_u[j]*transpose(previous_newton_velocity_gradients[q])*phi_u[i]
										+ previous_newton_velocity_values[q]*transpose(grad_phi_u[j])*phi_u[i]
										- phi_p[j]*div_phi_u[i] 
										- phi_p[i]*div_phi_u[j])
										*fe_values.JxW(q);
					}
					else
					{
					cell_matrix(i,j) += 
								(nu*double_contract(grad_phi_u[i],grad_phi_u[j])
											+ phi_u[j]*transpose(previous_newton_velocity_gradients[q])*phi_u[i]
											+ previous_newton_velocity_values[q]*transpose(grad_phi_u[j])*phi_u[i]
											- phi_p[j]*div_phi_u[i] 
											- phi_p[i]*div_phi_u[j])
										*fe_values.JxW(q);
					}
				}
				
				int equation_i = fe.system_to_component_index(i).first;
				
				cell_rhs[i] += 
							(fe_values.shape_value(i,q)*rhs_values[q](equation_i)
							 + previous_newton_velocity_values[q]*transpose(previous_newton_velocity_gradients[q])*phi_u[i])
							*fe_values.JxW(q);
			}
		}
		
		if (input.slit_flow)
		{
			/**add line integral term to rhs for inlet and outlet normal stress************************/
			for (unsigned int face=0; face < GeometryInfo<dim>::faces_per_cell; face++)
			{						
				if (cell->face(face)->at_boundary())
				{
					fe_face_values.reinit (cell, face);
								
					for (unsigned int q_boundary = 0; q_boundary < n_q_points_face; q_boundary++)
					{
						for (int k=0; k<dofs_per_cell; k++)
						{
							phi_u[k] = fe_face_values[velocities].value (k, q_boundary);
						}
			
						Point<dim> q_point = fe_face_values.quadrature_point(q_boundary);
						double x = q_point[0];
						
						for (int i = 0; i < dofs_per_cell; i++)
						{
							if (fabs(x) < 1e-8)
							{
								cell_rhs[i] -= input.r_in
										*phi_u[i][0]*fe_face_values.JxW(q_boundary);
							}
							else if (fabs(x - input.domain_length) < 1e-8)
							{
								cell_rhs[i] += input.r_out
										*phi_u[i][0]*fe_face_values.JxW(q_boundary);
							}
						}
					}
				}
			}
		}
		
		cell->get_dof_indices(local_dof_indices);
		constraints.distribute_local_to_global(cell_matrix, cell_rhs, local_dof_indices, 
						system_matrix, system_rhs);  
	}	
    
    timer.stop();
    if (input.verbose)
	{
		printf("done (%es)\n", timer());
	}                        
}

template<int dim>
void NavierStokesSolver<dim>::assemble_stokes_system(double nu)
{
	if (input.verbose)
	{
		printf("Assembling Stokes system...");
	}
	Timer timer;
	timer.start ();
	
	QGauss<dim>                 quadrature_formula(fe.degree + 1);
	QGauss<dim-1>               face_quadrature_formula(3);
	
	const int                   dofs_per_cell = fe.dofs_per_cell;
    const int                   n_q_points = quadrature_formula.size();	
    const int                   n_q_points_face = face_quadrature_formula.size();
    
	std::vector<Vector<double> > rhs_values (n_q_points, Vector<double>(dim+1));
	std::vector<Tensor<2,dim> > grad_phi_u (dofs_per_cell);
	std::vector<double>         div_phi_u (dofs_per_cell);
	std::vector<double>         phi_p (dofs_per_cell);
	std::vector<Tensor<1,dim> > phi_u (dofs_per_cell);  
	
	Vector<double>              cell_rhs(dofs_per_cell);
	FullMatrix<double>          cell_matrix (dofs_per_cell, dofs_per_cell);
    
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);    
	const FEValuesExtractors::Vector velocities (0);
	const FEValuesExtractors::Scalar pressure (dim);
	
    FEValues<dim> fe_values(fe, quadrature_formula, 
            update_values | update_gradients | update_JxW_values | update_quadrature_points);    
    FEFaceValues<dim> fe_face_values (fe, face_quadrature_formula, update_values 
			| update_quadrature_points | update_gradients | update_JxW_values);
			   
    typename DoFHandler<dim>::active_cell_iterator  
			cell = dof_handler.begin_active(), endc = dof_handler.end();
    
    for (; cell!=endc; ++cell)
    {
        fe_values.reinit(cell);
        cell_matrix = 0;
        cell_rhs = 0;
		
		forcing_function->vector_value_list(fe_values.get_quadrature_points(), rhs_values);
		
        //calculate cell contribution to system         
        for (int q = 0; q < n_q_points; q++)
        {
			for (int k=0; k<dofs_per_cell; k++)
			{
				grad_phi_u[k] = fe_values[velocities].gradient (k, q);
				div_phi_u[k]  = fe_values[velocities].divergence (k, q);
				phi_p[k]      = fe_values[pressure].value (k, q);
				phi_u[k]      = fe_values[velocities].value (k, q);
			}
			
			for (int i = 0; i < dofs_per_cell; i++)
			{
				for (int j = 0; j < dofs_per_cell; j++)
				{
					if (input.slit_flow)
					{
						cell_matrix(i,j) += 
								(0.5*nu*
								double_contract(grad_phi_u[i] + transpose(grad_phi_u[i]),
										grad_phi_u[j] + transpose(grad_phi_u[j]))
								- phi_p[i]*div_phi_u[j] 
								- phi_p[j]*div_phi_u[i])
								*fe_values.JxW(q);
					}
					else
					{
					cell_matrix(i,j) += 
								(nu*double_contract(grad_phi_u[i],grad_phi_u[j])
											- phi_p[j]*div_phi_u[i] 
											- phi_p[i]*div_phi_u[j])
										*fe_values.JxW(q);
					}
				}
				
				int equation_i = fe.system_to_component_index(i).first;
				
				cell_rhs[i] += 
							fe_values.shape_value(i,q)*rhs_values[q](equation_i)
							*fe_values.JxW(q); 
			}
		}
		
		if (input.slit_flow)
		{
			/**add line integral term to rhs for inlet and outlet normal stress************************/
			for (unsigned int face=0; face < GeometryInfo<dim>::faces_per_cell; face++)
			{						
				if (cell->face(face)->at_boundary())
				{
					fe_face_values.reinit (cell, face);
								
					for (unsigned int q_boundary = 0; q_boundary < n_q_points_face; q_boundary++)
					{
						for (int k=0; k<dofs_per_cell; k++)
						{
							phi_u[k] = fe_face_values[velocities].value (k, q_boundary);
						}
			
						Point<dim> q_point = fe_face_values.quadrature_point(q_boundary);
						double x = q_point[0];
						
						for (int i = 0; i < dofs_per_cell; i++)
						{
							if (fabs(x) < 1e-8)
							{
								cell_rhs[i] -= input.r_in
										*phi_u[i][0]*fe_face_values.JxW(q_boundary);
							}
							else if (fabs(x - input.domain_length) < 1e-8)
							{
								cell_rhs[i] += input.r_out
										*phi_u[i][0]*fe_face_values.JxW(q_boundary);
							}
						}
					}
				}
			}
		}
		
		cell->get_dof_indices(local_dof_indices);
		constraints.distribute_local_to_global(cell_matrix, cell_rhs, local_dof_indices, 
						system_matrix, system_rhs);  
	}	
    
    timer.stop();
    if (input.verbose)
	{
		printf("done (%es)\n", timer());
	}                        
}

template<int dim>
void NavierStokesSolver<dim>::run_newton_loop(int cycle)
{	
	int MAX_ITER = 10;
	int iter = 0;
	double residual = 0;
	
	assemble_stokes_system(input.nu);
	solve();
		
	previous_newton_step =  solution;

	printf("\nStarting Newton iteration for nu=%f\n", input.nu);
	while (iter == 0 || (residual > 1e-8 && iter < MAX_ITER))
	{		
		assemble_system(input.nu); 
		solve();
		
		Vector<double> res_vec = solution;
		res_vec -= previous_newton_step;
		
		residual = res_vec.l2_norm()/(dof_handler.n_dofs());
		previous_newton_step = solution;
		
		iter++;
		
		printf("Residual = %4.10e\n", residual);
	}
	
	if (iter == MAX_ITER)
	{
		printf("WARNING: Newton's method failed to converge for nu=%f\n", input.nu);
	}
		
	if (input.refinement_type.compare("adaptive") == 0)
	{
		output_results(cycle);
		
		if(input.print_error)
		{
			calculate_error(cycle);
		}
	}
}


template<int dim>
void NavierStokesSolver<dim>::solve()
{
	if (input.verbose)
	{
		printf("Solving linear system... ");
	}
	Timer timer;
	timer.start ();
	
	SparseDirectUMFPACK A_direct;
	A_direct.initialize(system_matrix);

	A_direct.vmult(solution, system_rhs);
    constraints.distribute(solution);
    
    timer.stop ();
    if (input.verbose)
	{
		printf("done (%gs)\n",timer());
	}    
}

template<int dim>
void NavierStokesSolver<dim>::refine_grid()
{
	if (input.verbose)
	{
		printf("Refining mesh...");
	}
	Timer timer;
	timer.start ();	
	
	/** Refine based on velocity gradient jump between cells***************************************/
	const FEValuesExtractors::Vector velocities (0);
	const FEValuesExtractors::Scalar pressure (dim);
	
	Vector<float> estimated_error_per_cell (mesh.n_active_cells());
	KellyErrorEstimator<dim>::estimate (dof_handler, QGauss<1>(dim + 1), 
			typename FunctionMap<dim>::type(), solution, estimated_error_per_cell, 
			fe.component_mask(velocities));
			
	GridRefinement::refine_and_coarsen_fixed_fraction (mesh, 
					estimated_error_per_cell, 0.7, 0.1);

	if (mesh.n_levels() > input.maxa)
	{
		for (typename Triangulation<dim>::active_cell_iterator
					cell = mesh.begin_active(input.maxa);
					cell != mesh.end(); ++cell)
		{
			cell->clear_refine_flag ();
		}
	}
	
	mesh.prepare_coarsening_and_refinement();
	mesh.execute_coarsening_and_refinement ();
	
	timer.stop();
	if (input.verbose)
	{
		printf("done (%gs)\n",timer());
	}
}

template<int dim>
void NavierStokesSolver<dim>::output_results(int cycle)
{	
	std::vector<std::string> solution_names;
	solution_names.push_back("u");
	solution_names.push_back("v");
	solution_names.push_back("pressure");
			
	DataOut<2> data_out;
	data_out.attach_dof_handler(dof_handler);
	data_out.add_data_vector(solution, solution_names);
			
	data_out.build_patches ();
	std::ostringstream filename;
	filename << input.output_name << cycle << ".vtk";
	std::ofstream output (filename.str().c_str());
	data_out.write_vtk (output);
}

template<int dim>
void NavierStokesSolver<dim>::calculate_error(int cycle)
{
	const ComponentSelectFunction<dim> velocity_mask(std::make_pair(0,dim),3);
	const ComponentSelectFunction<dim> pressure_mask(dim,3);
	
	/**Calculate velocity errors*******************************************************************/
    Vector<float> difference_per_cell (mesh.n_active_cells());
    VectorTools::integrate_difference (dof_handler, solution, *exact_solution, difference_per_cell,
                                       QGauss<dim>(3), VectorTools::L2_norm, &velocity_mask);
    const double L2_error_velocity = difference_per_cell.l2_norm();
    
    VectorTools::integrate_difference (dof_handler, solution, *exact_solution, difference_per_cell,
										QGauss<dim>(3),	VectorTools::H1_seminorm, &velocity_mask);
	const double H1_error_velocity = difference_per_cell.l2_norm();

	const QTrapez<1> q_trapez;
	const QIterated<dim> q_iterated (q_trapez, 5);
	VectorTools::integrate_difference (dof_handler,	solution, *exact_solution, difference_per_cell,
										q_iterated,	VectorTools::Linfty_norm, &velocity_mask);
	const double Linfty_error_velocity = difference_per_cell.linfty_norm();
	
	/** Calculate pressure errors******************************************************************/
	VectorTools::integrate_difference (dof_handler, solution, *exact_solution, difference_per_cell, 
										QGauss<dim>(3), VectorTools::L2_norm, &pressure_mask);
    const double L2_error_pressure = difference_per_cell.l2_norm();
    
    VectorTools::integrate_difference (dof_handler, solution, *exact_solution, difference_per_cell,
										QGauss<dim>(3), VectorTools::H1_seminorm, &pressure_mask);
	const double H1_error_pressure = difference_per_cell.l2_norm();

	VectorTools::integrate_difference (dof_handler,	solution, *exact_solution, difference_per_cell,
										q_iterated,	VectorTools::Linfty_norm, &pressure_mask);
	const double Linfty_error_pressure = difference_per_cell.linfty_norm();
	

	velocity_convergence_table.add_value("cycle", cycle);
	velocity_convergence_table.add_value("cells", mesh.n_active_cells());
	velocity_convergence_table.add_value("dofs", dof_handler.n_dofs());
	velocity_convergence_table.add_value("L2", L2_error_velocity);
	velocity_convergence_table.add_value("H1", H1_error_velocity);
	velocity_convergence_table.add_value("Linfty", Linfty_error_velocity);
	
	pressure_convergence_table.add_value("cycle", cycle);
	pressure_convergence_table.add_value("cells", mesh.n_active_cells());
	pressure_convergence_table.add_value("dofs", dof_handler.n_dofs());
	pressure_convergence_table.add_value("L2", L2_error_pressure);
	pressure_convergence_table.add_value("H1", H1_error_pressure);
	pressure_convergence_table.add_value("Linfty", Linfty_error_pressure);
}

template<int dim>
void NavierStokesSolver<dim>::run()
{
	read_inputs(); 
	
	if (input.refinement_type.compare("uniform") == 0)
	{
		double h = 1.0/(double) (input.nx-1);
		  
		for (int cycle = 0; cycle <= input.refinement_steps; cycle++)
		{				
			if (cycle!=0)
			{
				mesh.refine_global(1);
				h        /= 2;
			}
			
			if (input.verbose)
			{
				printf("Uniform refinement cycle %i\nh = %f\n", cycle, h);
			}
			
			setup_geometry(cycle);
			
			run_newton_loop(cycle);
			output_results(cycle);
			calculate_error(cycle);
		}		
	}
	else
	{
		setup_geometry(0);
		
		for (int cycle = 0; cycle <= input.refinement_steps; cycle++)
		{
			if (cycle!=0)
			{
				refine_grid();
				setup_geometry(cycle);
			}
			
			if (input.verbose)
			{
				printf("Adaptive refinement cycle %i\n", cycle);
			}
			
			run_newton_loop(cycle);
			solve();
			output_results(cycle);
			calculate_error(cycle);
			
		}
	}
	
	if (input.print_error)
	{
		print_errors();	
	}
}

template<int dim>
void NavierStokesSolver<dim>::print_errors()
{
	velocity_convergence_table.set_precision("L2", 3);
	velocity_convergence_table.set_precision("H1", 3);
	velocity_convergence_table.set_precision("Linfty", 3);
	velocity_convergence_table.set_scientific("L2", true);
	velocity_convergence_table.set_scientific("H1", true);
	velocity_convergence_table.set_scientific("Linfty", true);
	
	velocity_convergence_table.evaluate_convergence_rates("L2", ConvergenceTable::reduction_rate_log2);
	velocity_convergence_table.evaluate_convergence_rates("H1", ConvergenceTable::reduction_rate_log2);

	pressure_convergence_table.set_precision("L2", 3);
	pressure_convergence_table.set_precision("H1", 3);
	pressure_convergence_table.set_precision("Linfty", 3);
	pressure_convergence_table.set_scientific("L2", true);
	pressure_convergence_table.set_scientific("H1", true);
	pressure_convergence_table.set_scientific("Linfty", true);
	
	pressure_convergence_table.evaluate_convergence_rates("L2", ConvergenceTable::reduction_rate_log2);
	pressure_convergence_table.evaluate_convergence_rates("H1", ConvergenceTable::reduction_rate_log2);
	
	printf("\n");
	printf("Velocity error\n");
	velocity_convergence_table.write_text(std::cout);
	
	printf("\n");
	printf("Pressure error\n");
	pressure_convergence_table.write_text(std::cout);
}

int main ()
{
	const int dim = 2;
	
	NavierStokesSolver<dim> navier_stokes_problem;
	
	ExactSolutionForcingFunction<dim>  ff;
	ExactSolutionBoundaryValues<dim>   bv;
	ExactSolution<dim>                 ex;
	
	//MyZeroFunction<dim>                ff;
	//MyZeroFunction<dim>                bv;
	//StepBoundaryValues<dim>            bv;
	//ExactSolution<dim>                 ex;
	//MyZeroFunction<dim>                bv;
	//DrivenCavityBoundaryValues<dim>    bv;
	
	navier_stokes_problem.forcing_function   = &ff;
	navier_stokes_problem.boundary_values    = &bv;
	navier_stokes_problem.exact_solution     = &ex;
    navier_stokes_problem.run();
}
