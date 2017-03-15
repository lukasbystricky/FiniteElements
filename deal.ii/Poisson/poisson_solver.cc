#include "myfunctions.h"

/***************************************************************************************************
 * Poisson solver
 * 
 * Author: Lukas Bystricky, lb13f@my.fsu.edu
 * 
 * Solves Poisson equation
 * 
 * Supports adaptive meshing and convergence studies.
 * 
 * Boundary conditions, and exact solution are specified in "myfunctions.h". 
 * Exact solution is used only if errors are to be calculated.
 **************************************************************************************************/
template<int dim>
class PoissonSolver
{
	public:
		PoissonSolver(int degree, int mapping_order);
		void run();
			
		Function<dim> *forcing_function;
        Function<dim> *boundary_function;
        Function<dim> *exact_solution;
		
	private:
        void read_inputs();
		void setup_geometry(int cycle);
		void assemble_system();
		void solve();
		void refine_grid();
		void output_results(int cycle);
        void calculate_error(int cycle);  
        void print_errors();      
        
		struct InputFile
        {
            int          refinement_steps, p_degree, mapping_order;
            bool         print_error, verbose;
            std::string  output_name;
        };
        
        MappingQ<dim>          mapping;
        ParameterHandler 	   prm;   
        InputFile              input;        
		Triangulation<dim>     mesh;
		FE_Q<dim>              fe;
		DoFHandler<dim>        dof_handler;
		
		SparsityPattern        sparsity_pattern;
		SparseMatrix<double>   stiffness_matrix;
		
		Vector<double>         rhs;
		Vector<double>         solution;

		ConvergenceTable	   convergence_table;
};

template<int dim>
PoissonSolver<dim>::PoissonSolver(int degree, int mapping_order):
	/***********************************************************************************************
	 *  Lagrange elements
	 **********************************************************************************************/
    fe(degree),
    dof_handler(mesh),
    mapping(mapping_order)
{}

/***************************************************************************************************
 * Read and store inputs from user parameter file "input.in"
 **************************************************************************************************/
template<int dim>
void PoissonSolver<dim>::read_inputs()
{
	/***********************************************************************************************
	 * Set up paramter file patterns
	 **********************************************************************************************/
	
	prm.enter_subsection("Problem setup");
	{
		prm.declare_entry("basis degree", "2", Patterns::Integer(0,5),
						"degree of basis functions");
						
		prm.declare_entry("mapping order", "2", Patterns::Integer(0,5),
						"degree of mapping to boundary");
	}
	prm.leave_subsection();
	
	prm.enter_subsection("Run options");
	{							
		prm.declare_entry("refinement steps", "5", Patterns::Integer(0,20),
						"refinement steps");
	}
	prm.leave_subsection();  
													
	prm.enter_subsection("Output options");
	{
	  prm.declare_entry("output name", "poisson", Patterns::Anything(), 
						"name of the output file (without extension)");
	  
	  prm.declare_entry("print errors", "true", Patterns::Bool(),
						"print errors to console.");
						
	  prm.declare_entry("verbose", "true", Patterns::Bool(),
						"verbose output to console");
	}
	prm.leave_subsection();	
	
	/***********************************************************************************************
	 * Read and store paramters
	 **********************************************************************************************/
	 prm.read_input("input.in");
	
	 prm.enter_subsection("Problem setup");
	 input.p_degree = prm.get_double("basis degree");
	 input.mapping_order = prm.get_double("mapping order");
	 prm.leave_subsection();
	 
	 prm.enter_subsection("Run options");
	 input.refinement_steps = prm.get_double("refinement steps");
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
void PoissonSolver<dim>::setup_geometry (int cycle)
{  	
	if (input.verbose)
	{
		printf("Setting up geometry...");
	}
	Timer timer;
	timer.start ();
		
	if (cycle == 0)
	{
		Point<2> center (0,0);
		double radius = 1;
		GridGenerator::hyper_ball(mesh, center, radius);
		static const HyperBallBoundary<dim> boundary_description(center, radius);
		mesh.set_boundary(0, boundary_description);
	}
	
	//write mesh to file
	GridOut grid_out;
	GridOutFlags::Gnuplot gnuplot_flags(false, 30);
	grid_out.set_flags(gnuplot_flags);
	std::string filename = "mesh";
	filename+=('0'+input.mapping_order);
	filename+=".dat";
	
	std::ofstream gnuplot_file (filename.c_str());
	grid_out.write_gnuplot (mesh, gnuplot_file, &mapping);
  
    dof_handler.distribute_dofs(fe);  
    
    if (input.verbose)
	{
		std::printf("Number of active cells:%d\n", mesh.n_active_cells());
		std::printf("Number of degrees of freedom:%d\n", dof_handler.n_dofs()); 
	}
		
	CompressedSparsityPattern c_sparsity(dof_handler.n_dofs(), dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern (dof_handler, c_sparsity);
	sparsity_pattern.copy_from(c_sparsity);
	
	stiffness_matrix.reinit (sparsity_pattern);
	solution.reinit (dof_handler.n_dofs());
	rhs.reinit (dof_handler.n_dofs());	  
	
	timer.stop();								  
    if (input.verbose)
	{
		printf("done (%es)\n", timer());
	}     
}

template<int dim>
void PoissonSolver<dim>::assemble_system()
{
	if (input.verbose)
	{
		printf("Assembling matrix...");
	}
	Timer timer;
	timer.start ();
	
	MatrixCreator::create_laplace_matrix(mapping, dof_handler,
							QGauss<dim>(fe.degree+1), stiffness_matrix);
    
    VectorTools::create_right_hand_side(mapping, dof_handler,
					QGauss<dim>(fe.degree+1), *forcing_function, rhs);
    
    // apply boundary conditions
    std::map<types::global_dof_index,double> boundary_values;
	VectorTools::interpolate_boundary_values (mapping, dof_handler, 0, 
							*boundary_function, boundary_values);
							
	MatrixTools::apply_boundary_values (boundary_values, stiffness_matrix,
									  solution, rhs);
	
	timer.stop();								  
    if (input.verbose)
	{
		printf("done (%es)\n", timer());
	}                        
}

template<int dim>
void PoissonSolver<dim>::solve()
{
	if (input.verbose)
	{
		printf("Solving linear system... ");
	}
	Timer timer;
	timer.start ();
	

	SolverControl solver_control (10000, 1e-12);
    SolverCG<>    solver (solver_control);
	solver.solve (stiffness_matrix, solution, rhs,
                PreconditionIdentity());
    
    timer.stop ();
    if (input.verbose)
	{
		printf("done (%gs)\n",timer());
	}    
}

template<int dim>
void PoissonSolver<dim>::output_results(int cycle)
{	
	std::vector<std::string> solution_names;
	solution_names.push_back("u");
			
	DataOut<2> data_out;
	data_out.attach_dof_handler(dof_handler);
	data_out.add_data_vector(solution, solution_names);
			
	data_out.build_patches(1);
	std::ostringstream filename;
	filename << input.output_name << cycle << ".gpl";
	std::ofstream output (filename.str().c_str());
	data_out.write_gnuplot (output);
}

template<int dim>
void PoissonSolver<dim>::calculate_error(int cycle)
{

	/**Calculate errors*******************************************************************/
    Vector<float> difference_per_cell (mesh.n_active_cells());
    VectorTools::integrate_difference (mapping, dof_handler, solution, *exact_solution, difference_per_cell,
                                       QGauss<dim>(3), VectorTools::L2_norm);
    const double L2_error = difference_per_cell.l2_norm();
    
    VectorTools::integrate_difference (mapping, dof_handler, solution, *exact_solution, difference_per_cell,
										QGauss<dim>(3),	VectorTools::H1_seminorm);
	const double H1_error = difference_per_cell.l2_norm();

	const QTrapez<1> q_trapez;
	const QIterated<dim> q_iterated (q_trapez, 5);
	VectorTools::integrate_difference (mapping, dof_handler,	solution, *exact_solution, difference_per_cell,
										q_iterated,	VectorTools::Linfty_norm);
	const double Linfty_error = difference_per_cell.linfty_norm();
	

	convergence_table.add_value("cycle", cycle);
	convergence_table.add_value("cells", mesh.n_active_cells());
	convergence_table.add_value("h", GridTools::maximal_cell_diameter(mesh));
	convergence_table.add_value("dofs", dof_handler.n_dofs());
	convergence_table.add_value("L2", L2_error);
	convergence_table.add_value("H1", H1_error);
	convergence_table.add_value("Linfty", Linfty_error);
}

template<int dim>
void PoissonSolver<dim>::run()
{
	read_inputs(); 
	  
	for (int cycle = 0; cycle <= input.refinement_steps; cycle++)
	{						
		if (cycle!=0)
		{
			if (input.verbose)
			{
				printf("Refining mesh... ");
			}
			Timer timer;
			timer.start ();
	
			mesh.refine_global(1);
			
			timer.stop();								  
			if (input.verbose)
			{
				printf("done (%es)\n", timer());
			} 
		}				
		
		setup_geometry(cycle);
		
		if (input.verbose)
		{
			printf("Uniform refinement cycle %i\nh = %f\n", cycle, 
						GridTools::maximal_cell_diameter(mesh) );
		}
		
		assemble_system();
		solve();
		output_results(cycle);
		calculate_error(cycle);
	}		
	
	if (input.print_error)
	{
		print_errors();	
	}
}

template<int dim>
void PoissonSolver<dim>::print_errors()
{
	convergence_table.set_precision("h", 3);
	convergence_table.set_precision("L2", 3);
	convergence_table.set_precision("H1", 3);
	convergence_table.set_precision("Linfty", 3);
	convergence_table.set_scientific("h", true);
	convergence_table.set_scientific("L2", true);
	convergence_table.set_scientific("H1", true);
	convergence_table.set_scientific("Linfty", true);
	
	convergence_table.evaluate_convergence_rates("L2", ConvergenceTable::reduction_rate_log2);
	convergence_table.evaluate_convergence_rates("H1", ConvergenceTable::reduction_rate_log2);

	
	printf("\n");
	printf("Errors\n");
	convergence_table.write_text(std::cout);
}

int main ()
{
	const int dim = 2;
	
	// read in degree of basis function
	ParameterHandler prm;  
	
	prm.enter_subsection("Problem setup");
	{
		prm.declare_entry("basis degree", "2", Patterns::Integer(0,5),
						"degree of basis functions");
		prm.declare_entry("mapping order", "2", Patterns::Integer(0,5),
						"degree of mapping to boundary");
	}
	
	prm.read_input("input.in");
	prm.enter_subsection("Problem setup");
	int basis_degree = prm.get_double("basis degree");
	int map_degree = prm.get_double("mapping order");
	prm.leave_subsection();
	 
	 
	printf("basis degree %i, mapping degree %i\n", basis_degree, map_degree);
	PoissonSolver<dim> poisson_problem(basis_degree, map_degree);
	
	ExactSolutionForcingFunction<dim>  ff;
	ExactSolution<dim>           ex;
	BoundaryValues<dim>   		 bv;
	
	poisson_problem.forcing_function   = &ff;
	poisson_problem.boundary_function  = &bv;
	poisson_problem.exact_solution     = &ex;
    poisson_problem.run();
}
