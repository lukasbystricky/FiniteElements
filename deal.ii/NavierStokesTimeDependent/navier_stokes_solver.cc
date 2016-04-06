#include "myfunctions.h"

/***************************************************************************************************
 * Time dependent Navier-Stokes flow solver
 * 
 * Author: Lukas Bystricky, lb13f@my.fsu.edu
 * 
 * Solves Navier-Stokes equations. Supports adaptive meshing and convergence studies.
 * 
 * Boundary conditions, initial conditions and exact solution are specified in "myfunctions.h". 
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
		Function<dim> *initial_conditions;
		
	private:
        void read_inputs();
		void setup_geometry(int cycle);
		void assemble_system(double nu, double theta, double t);
		void solve();
		void run_time_loop();
		void refine_grid();
		void output_results(int cycle);
        void calculate_error(int cycle, double t);  
        void print_errors();  
        void calculate_lift_and_drag(double t); 
        
        FILE *pFile;
   
		struct InputFile
        {
            int          nx, ny, tpa, mina, maxa;
            double       nu, t0, tf, dt, r_in, r_out, domain_length, domain_height;
            bool         print_error, verbose, slit_flow;
            std::string  output_name, refinement_type, gmesh_file;
        };
        
        ParameterHandler 	   prm;   
        InputFile              input;  
        Triangulation<dim>	   mesh_hole;
        Triangulation<dim>     mesh_rectangle; 
        Triangulation<dim>	   mesh_top1;  
        Triangulation<dim>	   mesh_top2;   
		Triangulation<dim>     mesh;
		FESystem<dim>          fe;
		DoFHandler<dim>        dof_handler;
		
		ConstraintMatrix       constraints; 
		SparsityPattern        sparsity_pattern;
		SparseMatrix<double>   system_matrix;
		SparseMatrix<double>   system_matrix_copy;
		Vector<double>         system_rhs;
		Vector<double>         solution;
		Vector<double>         old_solution;
		Vector<double> 		   previous_newton_step;

		std::vector<double> 		   drags;
		
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
					
		prm.declare_entry("initial time", "0", Patterns::Double(0, 1000));
					
		prm.declare_entry("final time", "1", Patterns::Double(0, 1000));
					
		prm.declare_entry("time step", "0.05", Patterns::Double(1e-8, 1));
		
		prm.declare_entry("r inlet", "-0.980655e5", Patterns::Double(-1e7, 1e7),
					"normal stress at inlet (specify if slit flow)");
	
		prm.declare_entry("r outlet", "-9.8e4", Patterns::Double(-1e7, 1e7),
					"normal stress at outlet (specify if slit flow)");
	}
	prm.leave_subsection();
	
	prm.enter_subsection("Run options");
	{	
		prm.declare_entry("slit flow", "true", Patterns::Bool(),
						"slit flow model, (true or false)");
						
		prm.declare_entry("refinement type", "adaptive", Patterns::Anything(),
						"refienment type (adaptive or uniform)");
		
		prm.declare_entry("minimum refinement level", "0", Patterns::Integer(0, 3),
						"minimum cell refinement level");
		
		prm.declare_entry("maximum refinement level", "3", Patterns::Integer(0, 3),
						"maximum cell refinement level");
						
		prm.declare_entry("time steps per adapt", "5", Patterns::Integer(1, 1000),
						"number of time steps per adaptive mesh refinement");						
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
	 input.dt = prm.get_double("time step");
	 input.t0 = prm.get_double("initial time");
	 input.tf = prm.get_double("final time");
	 input.r_in = prm.get_double("r inlet");
	 input.r_out = prm.get_double("r outlet");
	 input.domain_length = prm.get_double("length");
	 input.domain_height = prm.get_double("height");
	 prm.leave_subsection();
	 
	 prm.enter_subsection("Run options");
	 input.slit_flow = prm.get_bool("slit flow");
	 input.refinement_type = prm.get("refinement type");
	 input.mina = prm.get_double("minimum refinement level");
	 input.maxa = prm.get_double("maximum refinement level");
	 input.tpa = prm.get_double("time steps per adapt");
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
			
			mesh.refine_global(2);
		}
		else
		{
			std::vector<unsigned int> number_elements(2);
			number_elements[0] = input.nx-1;
			number_elements[1] = input.ny-1;
	  
			Point<dim> bottom_left, top_right;
			bottom_left[0] = 0;
			bottom_left[1] = 0;
			top_right[0] = 1;
			top_right[1] = 1;
			
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

    constraints.clear (); 
				
    DoFTools::make_hanging_node_constraints (dof_handler, constraints);	    
    
    //set boundary labels:
    // (1) : inlet
    // (2) : outlet
    // (3) : top and bottom
    // (4) : everything else
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
				else if ((fabs(cell->face(f)->center()[1] - input.domain_height) < 1e-8) || 	
						 (fabs(cell->face(f)->center()[1]) < 1e-8))
				{
					cell->face(f)->set_boundary_indicator(3);
				}						 
				else
				{
					cell->face(f)->set_boundary_indicator(4);
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
				
		VectorTools::interpolate_boundary_values(dof_handler, 4, 
				*boundary_values, constraints, fe.component_mask(velocities));
						
		VectorTools::interpolate_boundary_values(dof_handler, 1, 
				*boundary_values, constraints, fe.component_mask(velocity_y));
				
		VectorTools::interpolate_boundary_values(dof_handler, 1, 
				*boundary_values, constraints, fe.component_mask(velocities));
		
		VectorTools::interpolate_boundary_values(dof_handler, 2, 
				*boundary_values, constraints, fe.component_mask(velocity_y));
	}
	else
	{
		    
		std::vector<bool> boundary_dofs(dof_handler.n_dofs(), false);
		DoFTools::extract_boundary_dofs(dof_handler, fe.component_mask(pressure), 
					boundary_dofs);
					
		//set Dirichlet boundary conditions on the velocity everywhere		
		VectorTools::interpolate_boundary_values(dof_handler, 1, 
				*boundary_values, constraints, fe.component_mask(velocities));
		
		VectorTools::interpolate_boundary_values(dof_handler, 2, 
				*boundary_values, constraints, fe.component_mask(velocities));
		
		VectorTools::interpolate_boundary_values(dof_handler, 3, 
				*boundary_values, constraints, fe.component_mask(velocities));

		VectorTools::interpolate_boundary_values(dof_handler, 4, 
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
    system_matrix_copy.reinit(sparsity_pattern);
    solution.reinit(dof_handler.n_dofs());     
    old_solution.reinit(dof_handler.n_dofs());
    system_rhs.reinit(dof_handler.n_dofs()); 
    previous_newton_step.reinit(dof_handler.n_dofs());  
}

template<int dim>
void NavierStokesSolver<dim>::assemble_system(double nu, double theta, double t)
{
	if (input.verbose)
	{
		printf("Assembling matrix...");
	}
	Timer timer;
	timer.start ();
	
	system_matrix.reinit(sparsity_pattern);    
    system_rhs.reinit(dof_handler.n_dofs());
    
	QGauss<dim>                 quadrature_formula(2);
	QGauss<dim-1>               face_quadrature_formula(3);
	
	const int                   dofs_per_cell = fe.dofs_per_cell;
    const int                   n_q_points = quadrature_formula.size();	
    const int                   n_q_points_face = face_quadrature_formula.size();
	
	std::vector<Tensor<1,dim> > previous_newton_velocity_values(n_q_points);
	std::vector<Tensor< 2, dim> > previous_newton_velocity_gradients(n_q_points);
	std::vector<Tensor<1,dim> > previous_time_velocity_values(n_q_points);
	std::vector<Tensor< 2, dim> > previous_time_velocity_gradients(n_q_points);
	std::vector<double> previous_time_pressure_values(n_q_points);
	
	std::vector<Vector<double> > rhs_values (n_q_points, Vector<double>(dim+1));
	std::vector<Vector<double> > rhs_values_old (n_q_points, Vector<double>(dim+1));
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
		
		fe_values[velocities].get_function_values(old_solution, previous_time_velocity_values);        
		fe_values[velocities].get_function_gradients(old_solution, previous_time_velocity_gradients);
		fe_values[pressure].get_function_values(old_solution, previous_time_pressure_values);
		
		forcing_function->set_time(t - input.dt);
		forcing_function->vector_value_list(fe_values.get_quadrature_points(), rhs_values_old);
		
		forcing_function->set_time(t);
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
								((1.0/input.dt)*phi_u[j]*phi_u[i] + theta*(0.5*nu*
								double_contract(grad_phi_u[i] + transpose(grad_phi_u[i]),
													grad_phi_u[j] + transpose(grad_phi_u[j]))
										+ phi_u[j]*transpose(previous_newton_velocity_gradients[q])*phi_u[i]
										+ previous_newton_velocity_values[q]*transpose(grad_phi_u[j])*phi_u[i]
										- phi_p[j]*div_phi_u[i])
										- phi_p[i]*div_phi_u[j])
										*fe_values.JxW(q);
					}
					else
					{
						cell_matrix(i,j) += 
								(
									phi_u[j]*phi_u[i] 
									+ theta*input.dt*
									(
										nu*double_contract(grad_phi_u[i],grad_phi_u[j])
										+ phi_u[j]*transpose(previous_newton_velocity_gradients[q])*phi_u[i]
										+ previous_newton_velocity_values[q]*transpose(grad_phi_u[j])*phi_u[i]										
									)
									- input.dt*phi_p[j]*div_phi_u[i]
									- phi_p[i]*div_phi_u[j]
								)
								*fe_values.JxW(q);
					}
				}
				
				int equation_i = fe.system_to_component_index(i).first;
				
				cell_rhs[i] += 
							(
								previous_time_velocity_values[q]*phi_u[i]
								+ input.dt*
								(
									theta*fe_values.shape_value(i,q)*rhs_values[q](equation_i)
									+ theta*previous_newton_velocity_values[q]*transpose(previous_newton_velocity_gradients[q])*phi_u[i]
									+ (1-theta)*fe_values.shape_value(i,q)*rhs_values_old[q](equation_i)
									- (1-theta)*
									(
										nu*double_contract(grad_phi_u[i],previous_time_velocity_gradients[q])
										+ previous_time_velocity_values[q]*transpose(previous_time_velocity_gradients[q])*phi_u[i]
										//- previous_time_pressure_values[q]*div_phi_u[i]
									)
								)
							 )
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
		for (int i = 0; i < dofs_per_cell; i++)
		{
			system_rhs[local_dof_indices[i]] += cell_rhs(i);
			
			for (int j = 0; j < dofs_per_cell; j++)
			{
				 system_matrix.add (local_dof_indices[i], local_dof_indices[j], cell_matrix(i,j));
			 }
		 }        
	}	
    
    timer.stop();
    if (input.verbose)
	{
		printf("done (%gs)\n", timer());
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
	A_direct.initialize(system_matrix_copy);

	A_direct.vmult(solution, system_rhs);
	/*SolverControl solver_control (4000, 1e-10);
	SolverGMRES<>  solver (solver_control);
	
	PreconditionJacobi<> preconditioner;
	preconditioner.initialize(system_matrix, 1.0);
  
	solver.solve (system_matrix_copy, solution, system_rhs,
                 PreconditionIdentity());*/
	
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
		
	for (typename Triangulation<dim>::active_cell_iterator
					cell = mesh.begin_active(input.mina);
					cell != mesh.end_active(input.mina); ++cell)
	{
			cell->clear_coarsen_flag ();
	}
	
	/** Move solution to new mesh******************************************************************/
	SolutionTransfer<dim> solution_trans(dof_handler);
	Vector<double> previous_solution;
	previous_solution = solution;
	mesh.prepare_coarsening_and_refinement();
	solution_trans.prepare_for_coarsening_and_refinement(previous_solution);

	mesh.execute_coarsening_and_refinement ();
	
	setup_geometry(1);
	solution_trans.interpolate(previous_solution, solution);
	
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
void NavierStokesSolver<dim>::calculate_error(int cycle, double t)
{
	exact_solution->set_time(t);
	
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
	velocity_convergence_table.add_value("time", t);
	velocity_convergence_table.add_value("cells", mesh.n_active_cells());
	velocity_convergence_table.add_value("dofs", dof_handler.n_dofs());
	velocity_convergence_table.add_value("L2", L2_error_velocity);
	velocity_convergence_table.add_value("H1", H1_error_velocity);
	velocity_convergence_table.add_value("Linfty", Linfty_error_velocity);
	
	pressure_convergence_table.add_value("cycle", cycle);
	pressure_convergence_table.add_value("time", t);
	pressure_convergence_table.add_value("cells", mesh.n_active_cells());
	pressure_convergence_table.add_value("dofs", dof_handler.n_dofs());
	pressure_convergence_table.add_value("L2", L2_error_pressure);
	pressure_convergence_table.add_value("H1", H1_error_pressure);
	pressure_convergence_table.add_value("Linfty", Linfty_error_pressure);
}

template<int dim>
void NavierStokesSolver<dim>::run_time_loop()
{
	FEValuesExtractors::Vector velocities(0);
	FEValuesExtractors::Scalar pressure(dim);
	
	VectorTools::interpolate(dof_handler, *initial_conditions, old_solution);
	solution = old_solution;
	
	if (input.refinement_type.compare("adaptive") == 0)
	{
		output_results(0);
	}
	
	double t = input.t0;
	double dt_tmp = input.dt;
	int    timestep_number = 0;
	int    n_adapts = 0;
	
	Vector<double> tmp;	
	tmp.reinit(solution.size());
	
	while (t < input.tf)
	{		
		t += input.dt;
		timestep_number++;
		
		if (input.verbose)
		{
			printf("Time step %i at t = %f\n", timestep_number, t);	
		}		
		
		boundary_values->set_time(t);
				
		constraints.clear();
		
		DoFTools::make_hanging_node_constraints(dof_handler, constraints);
		
		if (input.slit_flow)
		{
			//apply velocity boundary conditions to top and bottom only
			//set tangential (vertical) velocity at inlet and outlet
			
			FEValuesExtractors::Scalar velocity_y(1);
			
			VectorTools::interpolate_boundary_values(dof_handler, 4, 
					*boundary_values, constraints, fe.component_mask(velocities));
					
			VectorTools::interpolate_boundary_values(dof_handler, 3, 
					*boundary_values, constraints, fe.component_mask(velocities));
							
			VectorTools::interpolate_boundary_values(dof_handler, 1, 
					*boundary_values, constraints, fe.component_mask(velocity_y));
			
			VectorTools::interpolate_boundary_values(dof_handler, 1, 
				*boundary_values, constraints, fe.component_mask(velocities));
			
			VectorTools::interpolate_boundary_values(dof_handler, 2, 
					*boundary_values, constraints, fe.component_mask(velocity_y));
		}
		else
		{
			std::vector<bool> boundary_dofs(dof_handler.n_dofs(), false);
			DoFTools::extract_boundary_dofs(dof_handler, fe.component_mask(pressure), 
						boundary_dofs);
						
			//set Dirichlet boundary conditions on the velocity everywhere
			VectorTools::interpolate_boundary_values(dof_handler, 1, 
					*boundary_values, constraints, fe.component_mask(velocities));
			
			VectorTools::interpolate_boundary_values(dof_handler, 2, 
					*boundary_values, constraints, fe.component_mask(velocities));
			
			VectorTools::interpolate_boundary_values(dof_handler, 3, 
					*boundary_values, constraints, fe.component_mask(velocities));
	
			VectorTools::interpolate_boundary_values(dof_handler, 4, 
					*boundary_values, constraints, fe.component_mask(velocities));
						
			//constrain first pressure dof to be 0
			DoFTools::extract_boundary_dofs(dof_handler, fe.component_mask(pressure), 
						boundary_dofs);
						
			const unsigned int first_boundary_dof = std::distance(boundary_dofs.begin(),
				std::find (boundary_dofs.begin(), boundary_dofs.end(), true));
				
			constraints.add_line(first_boundary_dof);
		}
		
		constraints.close();
		
		/**Begin Newton loop**********************************************************************/
		int MAX_ITER = 10;
		int iter = 0;
		double residual = 0;
		
		previous_newton_step =  old_solution;
		
		printf("\nStarting Newton iteration for nu=%f\n", input.nu);
		while (iter == 0 || (residual > 1e-12 && iter < MAX_ITER))
		{
			assemble_system(input.nu, 0.5, t);

			system_matrix_copy.copy_from(system_matrix);
			constraints.condense(system_matrix_copy, system_rhs);
			
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
				
		calculate_lift_and_drag(t);
		
		if (input.refinement_type.compare("adaptive") == 0)
		{
			//output_results(timestep_number);
			
			/** Adapt if asked *********************************************************************/
			if (timestep_number % input.tpa == 0)
			{
				refine_grid();
				
				if(input.print_error)
				{
					calculate_error(n_adapts, t);
				}
				
				tmp.reinit(solution.size());
				n_adapts++;	
			}
		}
		
		input.dt = dt_tmp;
		old_solution = solution;
	}
}

template<int dim>
void NavierStokesSolver<dim>::run()
{
	read_inputs(); 	
	pFile = fopen ("drag_data_fine_dt_theta_1.txt","w");
	 
	if (input.refinement_type.compare("uniform") == 0)
	{
		double h = 1.0/(double) (input.nx-1);
		//double h = GridTools::maximal_cell_diameter(mesh)
		input.dt = pow(h,1);  
		  
		for (int cycle = 0; cycle <=input.maxa; cycle++)
		{				
			if (cycle!=0)
			{
				mesh.refine_global(1);
				h        /= 2;
				input.dt = pow(h,1);
			}
			
			if (input.verbose)
			{
				printf("Cycle %i\n h = %g, dt = %g\n", cycle, h, input.dt);	
			}
			
			setup_geometry(cycle);
			
			run_time_loop();
			output_results(cycle);
			calculate_error(cycle, input.tf);
		}		
	}
	else
	{
		setup_geometry(0);
		run_time_loop();
	}
	
	if (input.print_error)
	{
		print_errors();	
	}
	
	fclose (pFile);
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

template<int dim>
void NavierStokesSolver<dim>::calculate_lift_and_drag(double t)
{
	QGauss<dim-1>               face_quadrature_formula(3);
	const int                   n_q_points = face_quadrature_formula.size();
	
	const FEValuesExtractors::Vector velocities (0);
	const FEValuesExtractors::Scalar pressure (dim);
	
	std::vector<double>     pressure_values(n_q_points);
	std::vector<Tensor< 2, dim> > velocity_gradients(n_q_points);
	
	
	Tensor<1,dim>    grad_u_tau;	 
	Tensor<1,dim>    normal_vector; 
	Tensor<1,dim>    tangent_vector; 
	
	Tensor<2,dim>    fluid_stress;
	Tensor<2,dim>    fluid_pressure;
	Tensor<1,dim>    forces;
	
	FEFaceValues<dim> fe_face_values (fe, face_quadrature_formula, update_values 
			| update_quadrature_points | update_gradients | update_JxW_values | update_normal_vectors ); 
	
	typename DoFHandler<dim>::active_cell_iterator  
			cell = dof_handler.begin_active(), endc = dof_handler.end();
	
	double drag = 0;
	double lift = 0;
			
	for (; cell!=endc; ++cell)
    {
		for (int face=0; face < GeometryInfo<dim>::faces_per_cell; face++)
		{						
			if (cell->face(face)->at_boundary())
			{
				fe_face_values.reinit (cell, face);
				std::vector<Point<dim> > q_points = fe_face_values.get_quadrature_points();
				
				if (cell->face(face)->boundary_indicator() == 4)
				{	
					fe_face_values[velocities].get_function_gradients(solution, velocity_gradients);
					fe_face_values[pressure].get_function_values(solution, pressure_values);
							
					for (int q = 0; q < n_q_points; q++)
					{						
						normal_vector = -fe_face_values.normal_vector(q);
						
						fluid_pressure[0][0] = pressure_values[q];
						fluid_pressure[1][1] = pressure_values[q];
						 
						fluid_stress = input.nu*velocity_gradients[q] - fluid_pressure;
						
						forces = fluid_stress*normal_vector*fe_face_values.JxW(q);
						
						drag += forces[0];
						lift += forces[1];						
					}
				}
			}
		}
	}
	
	//calculate pressure drop
	Point<dim> p1, p2;
	p1[0] = 0.15;
	p1[1] = 0.2;
	p2[0] = 0.25;
	p2[1] = 0.2;
	Vector<double> solution_values1(dim+1);
	Vector<double> solution_values2(dim+1);
	
	VectorTools::point_value(dof_handler, solution, p1, solution_values1);
	VectorTools::point_value(dof_handler, solution, p2, solution_values2);
	
	double p_diff = solution_values1(dim) - solution_values2(dim);
	fprintf (pFile, "%f, %f,  %f, %f\n", t, 20.0*drag, 20.0*lift, p_diff);
	fflush (pFile);	
}

int main ()
{
	const int dim = 2;
	
	NavierStokesSolver<dim> navier_stokes_problem;
	
	//DrivenCavityBoundaryValues<dim>    bv;
	//InflowBoundaryValues<dim>                bv;
	PeriodicBoundaryValues<dim>                bv;
	//RotatingCircleForcingFunction<dim> ff;
	//MyZeroFunction<dim>                bv;
	MyZeroFunction<dim>                ff;
	MyZeroFunction<dim>                ex;
	MyZeroFunction<dim>                ic;
	
	/*ExactSolutionBoundaryValues<dim>                bv;
	ExactSolutionForcingFunction<dim>                ff;
	ExactSolution<dim>                ex;
	ExactSolutionInitialConditions<dim>                ic;*/
	
	navier_stokes_problem.forcing_function   = &ff;
	navier_stokes_problem.boundary_values    = &bv;
	navier_stokes_problem.exact_solution     = &ex;
	navier_stokes_problem.initial_conditions = &ic;
    navier_stokes_problem.run();
}
