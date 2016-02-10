#include "myfunctions.h"

template<int dim>
class DarcyFlow
{
	public:
		DarcyFlow();
		void run();
			
		Function<dim> *forcing_function;
        Function<dim> *boundary_values;
        Function<dim> *exact_solution;		
		
	private:
        void read_inputs();
		void setup_geometry(int cycle);
		void assemble();
		void solve();
		void refine_grid();
		void move_mesh();
		void calculate_new_top_values();
		double calculate_new_height();
		Point<dim,double> update_vertex(Point <dim,double> old_vertex);
		void output_results(int cycle);
        void calculate_error(int cycle);  
        void print_errors();      
        
		ParameterHandler prm;     
		Triangulation<dim>     mesh;
		FESystem<dim>          fe;
		DoFHandler<dim>        dof_handler;
		
		ConstraintMatrix       constraints;
		SparsityPattern        sparsity_pattern;
		SparseMatrix<double>   system_matrix;
		Vector<double>         system_rhs;
		Vector<double>         solution;		
		
		ConvergenceTable	   velocity_convergence_table;
		ConvergenceTable	   pressure_convergence_table;
		bool pressure_induced;
		bool verbose;
		
		struct PointComparator
		{
			bool operator()(const Point<dim> & left, const Point<dim> & right) const
			{
				return (left[0] > right[0]);
			}
		};
		
		std::map<Point<dim>, double, PointComparator>       new_domain_heights;
};

template<int dim>
DarcyFlow<dim>::DarcyFlow():
	//Raviart-Thomas elements for the velocity, 
	//discontinous elements for the pressure
    fe (FE_RaviartThomas<dim>(1), 1, FE_DGQ<dim>(1), 1)
    dof_handler(mesh)
{}

template<int dim>
void DarcyFlow<dim>::read_inputs()
{
	prm.declare_entry("mu", "0.1", Patterns::Double(1e-7, 10000), 
					"dynamic viscosity of the fluid (in Pa s)");

	prm.declare_entry("permiability", "1e-3", Patterns::Double(1e-7, 1),
					"average permiability (in mm)");

	prm.declare_entry("pressure induced", "true", Patterns::Bool(),
					"pressure induced flow");
					
	prm.declare_entry("r inlet", "-1e5", Patterns::Double(-1e7, 1e7),
					"normal stress at inlet (specify if pressure induced)");
	
	prm.declare_entry("r outlet", "-1e4", Patterns::Double(-1e7, 1e7),
					"normal stress at outlet (specify if pressure induced)");
									
	prm.enter_subsection("Geometry");
	{
	  prm.declare_entry("initial height", "4", Patterns::Double(0.1, 100),
						"initial height of the foam (in mm)");

	  prm.declare_entry("length", "100", Patterns::Double(1, 100000), 
						"length of the foam (in mm)");
						
	  prm.declare_entry("nx", "250", Patterns::Integer(2, 10000),
						"number of mesh points along foam length");

	  prm.declare_entry("ny", "10", Patterns::Integer(2, 10000),
						"number of mesh points along foam height");
	}
	prm.leave_subsection();

	prm.enter_subsection("Output options");
	{
	  prm.declare_entry("output name", "darcy_flow", Patterns::Anything(), 
						"name of the output file (without extension)");
	  
	  prm.declare_entry("verbose", "true", Patterns::Bool(),
						"verbose output to console during solve");
	}
	prm.leave_subsection();
	
	prm.enter_subsection("Run options");
	{
		prm.declare_entry("convergence study", "false", Patterns::Bool(),
						"run convergence study");
		
		prm.declare_entry("number of convergence cycles", "3", Patterns::Integer(1, 10),
						"number of global refinements in convergence study");
		
		prm.declare_entry("adaptive meshing", "true", Patterns::Bool(),
						"adaptive meshing");
		
		prm.declare_entry("number of adaptive cycles", "3", Patterns::Integer(1, 10),
						"number of adaptive refinements to run");
	}
	prm.leave_subsection();  
	
	prm.read_input("input.in");
}

template<int dim>
void DarcyFlow<dim>::setup_geometry (int cycle)
{  			
	if (cycle == 0)
	{
		prm.enter_subsection("Geometry");
		
		std::vector<unsigned int> number_elements(2);
		number_elements[0] = prm.get_integer("nx")-1;
		number_elements[1] = prm.get_integer("ny")-1;
  
		Point<2> bottom_left (0,0);
		Point<2> top_right (prm.get_double("length"), prm.get_double("initial height"));
		
		GridGenerator::subdivided_hyper_rectangle(mesh, number_elements,
			bottom_left, top_right, false);
			
		prm.leave_subsection();	
	}   	
	
    dof_handler.distribute_dofs(fe);  
    
    if (verbose)
    {
		printf("Number of active cells:%d\n", mesh.n_active_cells());
		printf("Number of degrees of freedom:%d\n", dof_handler.n_dofs()); 
	}	
	
	FEValuesExtractors::Vector velocities(0);	
    FEValuesExtractors::Scalar pressure(dim);

    constraints.clear (); 				
    DoFTools::make_hanging_node_constraints (dof_handler, constraints);	       
    std::vector<bool> boundary_dofs(dof_handler.n_dofs(), false);
    
	prm.enter_subsection("Geometry");
	double domain_height = prm.get_double("initial height");
	double domain_length = prm.get_double("length");
	prm.leave_subsection();
		
    typename DoFHandler<dim>::active_cell_iterator  
			cell = dof_handler.begin_active(), endc = dof_handler.end();
    
    //set boundary labels:
    // (1) : bottom/top
    // (2) : inlet/outlet
	for (; cell!=endc; ++cell)
	{
		for (int f=0; f<GeometryInfo<dim>::faces_per_cell; f++)
		{
			if (cell->face(f)->at_boundary())
			{
				if ((fabs(cell->face(f)->center()[1]) < 1e-8) ||
					(fabs(cell->face(f)->center()[1] - domain_height) < 1e-8))
				{
					cell->face(f)->set_boundary_indicator(1);
				}				
				if ((fabs(cell->face(f)->center()[0]) < 1e-8) ||
					(fabs(cell->face(f)->center()[0] - domain_length) < 1e-8))
				{
					cell->face(f)->set_boundary_indicator(2);
				}
			}
		}
	}		
		
	if (pressure_induced)
	{
		//apply velocity boundary conditions to top and bottom only
		//set tangential (vertical) velocity at inlet and outlet
		
		FEValuesExtractors::Scalar velocity_y(1);
		
		VectorTools::interpolate_boundary_values(dof_handler, 1, 
				*boundary_values, constraints, fe.component_mask(velocities));
				
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
						
		DoFTools::extract_boundary_dofs(dof_handler, fe.component_mask(pressure), 
					boundary_dofs);
					
		const unsigned int first_boundary_dof = std::distance(boundary_dofs.begin(),
			std::find (boundary_dofs.begin(), boundary_dofs.end(), true));
			
		constraints.add_line(first_boundary_dof);	
	}

    constraints.close();
            
	//calculate sparsity pattern 
	CompressedSparsityPattern c_sparsity(dof_handler.n_dofs(), dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern (dof_handler, c_sparsity);
    
    constraints.condense(c_sparsity);
	sparsity_pattern.copy_from(c_sparsity);
    
    //initialize vectors and sparse matrices
    system_matrix.reinit(sparsity_pattern);
    solution.reinit(dof_handler.n_dofs());     
    system_rhs.reinit(dof_handler.n_dofs());      
}

template<int dim>
void DarcyFlow<dim>::assemble()
{	
	Timer timer;
	if (verbose)
	{
		printf("Assembling matrix...");		
		timer.start ();
	}
	
	QGauss<dim>                 quadrature_formula(3);
	QGauss<dim-1>               face_quadrature_formula(3);
	
	const int                   dofs_per_cell = fe.dofs_per_cell;
    const int                   n_q_points = quadrature_formula.size();	
    const int                   n_q_points_face = face_quadrature_formula.size();
	double 						mu = prm.get_double("mu");
	double						kappa = prm.get_double("permiability");
	
	prm.enter_subsection("Geometry");
	double domain_length = prm.get_double("length");
	prm.leave_subsection();
	
	std::vector<Tensor<2,dim> >  grad_phi_u (dofs_per_cell);
	std::vector<double>          div_phi_u (dofs_per_cell);
	std::vector<double>          phi_p (dofs_per_cell);
	std::vector<Tensor<1,dim> >  phi_u (dofs_per_cell);
	std::vector<Vector<double> > rhs_values (n_q_points, Vector<double>(dim+1)); 
	
	FullMatrix<double>          cell_matrix (dofs_per_cell, dofs_per_cell);
    Vector<double>              cell_rhs (dofs_per_cell);
    
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
        forcing_function->vector_value_list(fe_values.get_quadrature_points(), rhs_values);	
        
        cell_matrix = 0;
        cell_rhs  = 0;
        
        //calculate cell contribution to system         
        for (int q = 0; q < n_q_points; q++)
        {
			Point<dim> q_point = fe_values.quadrature_point(q);
			
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
						cell_matrix(i,j) += 
									((mu/2.0)*double_contract(grad_phi_u[j] + transpose(grad_phi_u[j]),grad_phi_u[i] + transpose(grad_phi_u[i]))
											+ (mu/(fabs(kappa)))*phi_u[j]*phi_u[i]										
											- phi_p[j]*div_phi_u[i]
											- phi_p[i]*div_phi_u[j])
											*fe_values.JxW(q);	
												
				}
				
				int equation_i = fe.system_to_component_index(i).first;					
								
				cell_rhs[i] += 
							fe_values.shape_value(i,q)*rhs_values[q](equation_i)
							*fe_values.JxW(q);	  
			}
		}
		
		
		//add line integral term to rhs for inlet and outlet pressure
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
							cell_rhs[i] -= prm.get_double("r inlet")
									*phi_u[i][0]*fe_face_values.JxW(q_boundary);
						}
						else if (fabs(x - domain_length) < 1e-8)
						{
							cell_rhs[i] += prm.get_double("r outlet")
									*phi_u[i][0]*fe_face_values.JxW(q_boundary);
						}						
					}
				}
			}
		}
		
		cell->get_dof_indices(local_dof_indices);
		constraints.distribute_local_to_global(cell_matrix, cell_rhs, 
									local_dof_indices, system_matrix, system_rhs);  
	}	
    
    if (verbose)
    {
		timer.stop();
		printf("done (%gs)\n", timer()); 
	}                       
}

template<int dim>
void DarcyFlow<dim>::solve()
{
	Timer timer;
	if (verbose)
	{
		printf("Solving linear system... ");
		
		timer.start ();
	}
	
	//use direct solver
	SparseDirectUMFPACK A_direct;
	A_direct.initialize(system_matrix);

	A_direct.vmult(solution, system_rhs);	
    constraints.distribute(solution);
    
    if (verbose)
	{
		timer.stop ();
		printf("done (%gs)\n",timer());
	}
}

template<int dim>
void DarcyFlow<dim>::refine_grid()
{
	Timer timer;
	if (verbose)
	{
		printf("Refining mesh...");		
		timer.start ();	
	}
	
	//Refine based on either velocity error or pressure error
	const FEValuesExtractors::Vector velocities (0);
	const FEValuesExtractors::Scalar pressure (2);
	
	Vector<float> estimated_error_per_cell (mesh.n_active_cells());
	KellyErrorEstimator<2>::estimate (dof_handler, QGauss<1>(3), 
			FunctionMap<2>::type(), solution, estimated_error_per_cell, 
			fe.component_mask(velocities));
			
	GridRefinement::refine_and_coarsen_fixed_number(mesh,
			estimated_error_per_cell, 0.3, 0.03);
	
	mesh.execute_coarsening_and_refinement ();
	
	if (verbose)
	{
		timer.stop();
		printf("done (%gs)\n",timer());
	}
}

template<int dim>
void DarcyFlow<dim>::move_mesh()
{
	if (verbose)
	{
		printf("Moving mesh...");
	}
	
	calculate_new_top_values();
	
	std::vector<bool> hanging_dofs(dof_handler.n_dofs(), false);
	DoFTools::extract_hanging_node_dofs(dof_handler, hanging_dofs);
					
	std::vector<bool> vertex_touched (mesh.n_vertices(), false);
	typename DoFHandler<dim>::active_cell_iterator  
			cell = dof_handler.begin_active(), endc = dof_handler.end();
			
	for (; cell!=endc; ++cell)
	{
		for (int v=0; v<GeometryInfo<dim>::vertices_per_cell; v++)
		{
			if (vertex_touched[cell->vertex_index(v)] == false)
			{				
				vertex_touched[cell->vertex_index(v)] = true;				
				int dof_number = cell->vertex_dof_index(v, 0);
				
				printf(hanging_dofs[dof_number] ? "true\n" : "false\n");
				if (hanging_dofs[dof_number] == false)
				{
					cell->vertex(v) = update_vertex(cell->vertex(v));
				}
			}
		}
	}
	
	if (verbose)
	{
		std::ofstream out("new_grid.eps");
		GridOut grid_out;
		grid_out.write_eps(mesh, out);
		
		printf("done\n");
	}
}

template<int dim>
Point<dim> DarcyFlow<dim>::update_vertex(Point<dim> old_vertex)
{
		double x_0 = 0;		
		double x_1 = std::numeric_limits<double>::infinity(); 
		double y_0, y_1;
		
		prm.enter_subsection("Geometry");
		double initial_height = prm.get_double("initial height");
		prm.leave_subsection();
	
		for(typename std::map<Point<dim>, double, PointComparator>::iterator 
						ii=new_domain_heights.begin(); ii!=new_domain_heights.end(); ++ii)
		{
			Point<dim> test_point = (*ii).first;
			if (test_point[0] >= x_0 && test_point[0] <= old_vertex[0])
			{
				x_0 = test_point[0];
				y_0 = (*ii).second;
			}
			else if (test_point[0] < x_1 && test_point[0] > old_vertex[0])
			{
				x_1 = test_point[0];
				y_1 = (*ii).second;
			}
		}
		
		double y_bar = y_0 + (old_vertex[0] - x_0)*((y_1-y_0)/(x_1-x_0));
		Point<dim> new_vertex;
		new_vertex[0] = old_vertex[0];	
		new_vertex[1] = y_bar*old_vertex[1]/initial_height;	
		
		return new_vertex;
}

template<int dim>
void DarcyFlow<dim>::calculate_new_top_values()
{
	prm.enter_subsection("Geometry");
	double initial_height = prm.get_double("initial height");
	prm.leave_subsection();
		
	std::vector<bool> vertex_touched (mesh.n_vertices(), false);
	typename DoFHandler<dim>::active_cell_iterator  
			cell = dof_handler.begin_active(), endc = dof_handler.end();
			
	for (; cell!=endc; ++cell)
	{
		for (int v=0; v<GeometryInfo<dim>::vertices_per_cell; v++)
		{
			if (vertex_touched[cell->vertex_index(v)] == false)
			{
				Point<dim> vertex = cell->vertex(v);
				
				if (fabs(vertex[1] - initial_height) < 1e-8)
				{
					vertex_touched[cell->vertex_index(v)] = true;
					new_domain_heights[cell->vertex(v)] = calculate_new_height();
				}
			}
		}
	}
}

template<int dim>
double DarcyFlow<dim>::calculate_new_height()
{
	prm.enter_subsection("Geometry");
	double initial_height = prm.get_double("initial height");
	prm.leave_subsection();

	return initial_height + 0.1*(rand()/((double)RAND_MAX) - 0.5);
}

template<int dim>
void DarcyFlow<dim>::output_results(int cycle)
{	
	std::vector<std::string> solution_names (2, "velocity");
	solution_names.push_back ("pressure");
	
	std::vector<DataComponentInterpretation::DataComponentInterpretation>
		data_component_interpretation(2, 
				DataComponentInterpretation::component_is_part_of_vector);
	data_component_interpretation.push_back(
			DataComponentInterpretation::component_is_scalar);
			
	DataOut<2> data_out;
	data_out.attach_dof_handler(dof_handler);
	data_out.add_data_vector(solution, solution_names, DataOut<2>::type_dof_data,
			data_component_interpretation);
			
	data_out.build_patches ();
	std::ostringstream filename;
	
	prm.enter_subsection("Output options");
	filename << prm.get("output name") << cycle << ".vtk";
	prm.leave_subsection();
	
	std::ofstream output (filename.str().c_str());
	data_out.write_vtk (output);
}

template<int dim>
void DarcyFlow<dim>::calculate_error(int cycle)
{	
	const ComponentSelectFunction<dim> velocity_mask(std::make_pair(0,dim),3);
	const ComponentSelectFunction<dim> pressure_mask(dim,3);
	
	//Calculate velocity errors*************************************************
    Vector<float> difference_per_cell (mesh.n_active_cells());
    VectorTools::integrate_difference (dof_handler,
                                       solution,
                                       *exact_solution,
                                       difference_per_cell,
                                       QGauss<dim>(3),
                                       VectorTools::L2_norm, 
                                       &velocity_mask);
    const double L2_error_velocity = difference_per_cell.l2_norm();
    
    VectorTools::integrate_difference (dof_handler,
										solution,
										*exact_solution,
										difference_per_cell,
										QGauss<dim>(3),
										VectorTools::H1_seminorm,
										&velocity_mask);
	const double H1_error_velocity = difference_per_cell.l2_norm();

	const QTrapez<1> q_trapez;
	const QIterated<dim> q_iterated (q_trapez, 5);
	VectorTools::integrate_difference (dof_handler,
										solution,
										*exact_solution,
										difference_per_cell,
										q_iterated,
										VectorTools::Linfty_norm,
										&velocity_mask);
	const double Linfty_error_velocity = difference_per_cell.linfty_norm();
	
	//Calculate pressure errors*************************************************
	VectorTools::integrate_difference (dof_handler,
                                       solution,
                                       *exact_solution,
                                       difference_per_cell,
                                       QGauss<dim>(3),
                                       VectorTools::L2_norm, 
                                       &pressure_mask);
    const double L2_error_pressure = difference_per_cell.l2_norm();
    
    VectorTools::integrate_difference (dof_handler,
										solution,
										*exact_solution,
										difference_per_cell,
										QGauss<dim>(3),
										VectorTools::H1_seminorm,
										&pressure_mask);
	const double H1_error_pressure = difference_per_cell.l2_norm();

	VectorTools::integrate_difference (dof_handler,
										solution,
										*exact_solution,
										difference_per_cell,
										q_iterated,
										VectorTools::Linfty_norm,
										&pressure_mask);
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
void DarcyFlow<dim>::run()
{
	read_inputs(); 
	prm.print_parameters (std::cout, ParameterHandler::Text);

	pressure_induced = prm.get_bool("pressure induced");
	
	prm.enter_subsection("Output options");
	verbose = prm.get_bool("verbose");
	prm.leave_subsection();
	
	prm.enter_subsection("Run options");	
	bool convergence_study = prm.get_bool("convergence study");
	bool adaptive_meshing = prm.get_bool("adaptive meshing");
	prm.leave_subsection();
	
	if (convergence_study)
	{
		prm.enter_subsection("Run options");
		int nc = prm.get_integer("number of convergence cycles");
		prm.leave_subsection();		
		
		prm.enter_subsection("Geometry");
		double dx = prm.get_double("length")/(double) (prm.get_integer("nx")-1);	
		double dy = prm.get_double("initial height")/(double) (prm.get_integer("ny")-1);
		prm.leave_subsection();	  
		  
		for (int cycle = 0; cycle <= nc; cycle++)
		{		
			
			if (cycle!=0)
			{
				mesh.refine_global(1);
				dx        /= 2;
				dy        /= 2;
			}
			
			printf("Cycle %i, dx = %f, dy = %f\n", cycle, dx, dy);	
			
			setup_geometry(cycle);			
			assemble();
			solve();
			
			output_results(cycle);
			calculate_error(cycle);			
		}
		
		print_errors();
	}
	
	if (adaptive_meshing)
	{
		prm.enter_subsection("Run options");
		int nc = prm.get_integer("number of adaptive cycles");
		prm.leave_subsection();
		for (int cycle = 0; cycle <= nc; cycle++)
		{		
				
			printf("Adpative cycle %i\n", cycle);
			
			if (cycle!=0)
			{
				refine_grid();
				
			}

			setup_geometry(cycle);														
			assemble();
			solve();
			
			output_results(cycle);
		}
		
		move_mesh();
		assemble();
		solve();
		output_results(nc+1);
	}
}

template<int dim>
void DarcyFlow<dim>::print_errors()
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
	
	DarcyFlow<dim> stokes_problem;
	/*ExactSolutionForcingFunction<dim>  ff;
	ExactSolutionBoundaryValues<dim>   bv;
	ExactSolution<dim>            ex;*/
	
	MyZeroFunction<dim>  ff;
	MyZeroFunction<dim>   bv;
	ExactSolutionPressureDriven<dim>            ex;
	
	stokes_problem.forcing_function   = &ff;
	stokes_problem.boundary_values    = &bv;
	stokes_problem.exact_solution     = &ex;
    stokes_problem.run();
}
