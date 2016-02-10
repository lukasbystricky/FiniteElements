#include "myfunctions.h"

template<int dim>
class DarcyFlow
{
	public:
		DarcyFlow(const std::vector<const FiniteElement<dim> *> fes, 
				  const std::vector<unsigned int> multiplicities);
		void run();
			
		Function<dim> *forcing_function;
        Function<dim> *boundary_values;	
		
	private:
        void read_inputs();
		void setup_geometry();
		void assemble();
		void solve();
		void move_mesh();
		void calculate_displacement_x();
		void calculate_displacement_y();
		double calculate_permiability(double porosity);
		double calculate_porosity_1d_integral(Point<dim, double> vertex);
		void calculate_initial_porosity();
		void calculate_stress();
		void calculate_gamma();
		double calculate_c_1(Point<dim> p,
							 double delta_y);
		Vector<double> get_gamma(const double y,
								 const Vector<double> &c_1) const;
		void calculate_derivatives(Point<dim> p,
								   int component, 
								   double &partial_x, 
								   double &partial_y, 
								   int dof_number);
		void output_results(int iter);     
        
		ParameterHandler prm;     
		Triangulation<dim>     mesh;
		FESystem<dim>          fe;
		DoFHandler<dim>        dof_handler;
		
		ConstraintMatrix       constraints;
		SparsityPattern        sparsity_pattern;
		SparseMatrix<double>   system_matrix;
		Vector<double>         system_rhs;
		Vector<double>         solution;	
		Vector<double>         previous_solution;
		
		double average_porosity;		

		std::vector<bool> top_boundary_dofs;
		std::vector<bool> bottom_boundary_dofs;
		std::vector<bool> inlet_boundary_dofs;
		std::vector<bool> outlet_boundary_dofs;
		
		bool pressure_induced;
		bool verbose;
		
		int nx, ny, iter;
};

template<int dim>
DarcyFlow<dim>::DarcyFlow(const std::vector<const FiniteElement<dim> *> fes, 
						  const std::vector<unsigned int> multiplicities):
    fe(fes, multiplicities),
	dof_handler(mesh)
{}

template<int dim>
void DarcyFlow<dim>::read_inputs()
{
	prm.declare_entry("mu", "0.1", Patterns::Double(1e-7, 10000), 
					"dynamic viscosity of the fluid (in Pa s)");
					
	prm.declare_entry("Young's modulus", "1e5", Patterns::Double(1e-7, 1e10), 
					"Young's modulus of the foam");
	
	prm.declare_entry("Poisson's ratio", "0.1", Patterns::Double(1e-7, 1), 
					"Poisson's ratio of the foam");
					
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

	  prm.declare_entry("length", "100", Patterns::Double(0.1, 100000), 
						"length of the foam (in mm)");
						
	  prm.declare_entry("nx", "250", Patterns::Integer(2, 10000),
						"number of mesh points along foam length");

	  prm.declare_entry("ny", "10", Patterns::Integer(2, 10000),
						"number of mesh points along foam height");
	}
	prm.leave_subsection();

	prm.enter_subsection("Constants");
	{
	  prm.declare_entry("permiability coefficent a", "1", Patterns::Double(0, 100), 
						"coefficient a in k = a*phi^b");
	  
	  prm.declare_entry("permiability coefficent b", "2", Patterns::Double(0, 100), 
						"coefficient b in k = a*phi^b");
						
	  prm.declare_entry("min inner pore radius", "1e-6", Patterns::Double(1e-8, 1), 
						"minimum void space radius (mm)");
						
	  prm.declare_entry("max inner pore radius", "1e-4", Patterns::Double(1e-8, 1), 
						"maximum void space radius (mm)");
						
	  prm.declare_entry("min outer pore radius", "1.1e-4", Patterns::Double(1e-8, 1), 
						"minimum partial void space radius (mm)");
						
	  prm.declare_entry("max outer pore radius", "5e-4", Patterns::Double(1e-8, 1), 
						"maximum partial void space radius (mm)");
						
	  prm.declare_entry("number of pores", "1", Patterns::Integer(0, 100000), 
						"number of pores");
						
	  prm.declare_entry("max porosity", "0.95", Patterns::Double(0.5, 0.99999), 
						"maximum porosity (must be < 1)");
		
	  prm.declare_entry("min porosity", "0.01", Patterns::Double(1e-8, 0.49999), 
						"minimum porosity (must be > 0)");				
						
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
	
	prm.read_input("input.in");
	
	prm.enter_subsection("Geometry");
	nx = prm.get_integer("nx");
	ny = prm.get_integer("ny");
	prm.leave_subsection();
}

template<int dim>
void DarcyFlow<dim>::setup_geometry()
{  			
	if (iter == 0)
	{
		/*GridIn<dim> grid_in;
		grid_in.attach_triangulation(mesh);
		std::ifstream input_file("mesh.msh");
		grid_in.read_msh(input_file);	*/	
		
		prm.enter_subsection("Geometry");
		
		std::vector<unsigned int> number_elements(2);
		number_elements[0] = nx-1;
		number_elements[1] = ny-1;
  
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
	FEValuesExtractors::Scalar porosity(dim+1);
	FEValuesExtractors::Scalar stress_xx(dim+2);
	FEValuesExtractors::Scalar stress_yy(dim+3);
	FEValuesExtractors::Scalar stress_xy(dim+4);
	FEValuesExtractors::Scalar strain_xx(dim+5);
	FEValuesExtractors::Scalar strain_yy(dim+6);
	FEValuesExtractors::Scalar strain_xy(dim+7);
	FEValuesExtractors::Scalar disp_x(dim+8);
	FEValuesExtractors::Scalar disp_y(dim+9);
	FEValuesExtractors::Scalar gamma(dim+10);
	FEValuesExtractors::Scalar stiffness(dim+11);
	
    constraints.clear (); 				
    DoFTools::make_hanging_node_constraints (dof_handler, constraints);	       
    std::vector<bool> boundary_dofs(dof_handler.n_dofs(), false);
    std::vector<bool> porosity_dofs(dof_handler.n_dofs(), false);
    std::vector<bool> stress_xx_dofs(dof_handler.n_dofs(), false);
    std::vector<bool> stress_yy_dofs(dof_handler.n_dofs(), false);
    std::vector<bool> stress_xy_dofs(dof_handler.n_dofs(), false);
    std::vector<bool> strain_xx_dofs(dof_handler.n_dofs(), false);
    std::vector<bool> strain_yy_dofs(dof_handler.n_dofs(), false);
    std::vector<bool> strain_xy_dofs(dof_handler.n_dofs(), false);
    std::vector<bool> disp_x_dofs(dof_handler.n_dofs(), false);
    std::vector<bool> disp_y_dofs(dof_handler.n_dofs(), false);
    std::vector<bool> gamma_dofs(dof_handler.n_dofs(), false);
    std::vector<bool> stiffness_dofs(dof_handler.n_dofs(), false);
    
    system_rhs.reinit(dof_handler.n_dofs());   
    
	prm.enter_subsection("Geometry");
	double domain_height = prm.get_double("initial height");
	double domain_length = prm.get_double("length");
	prm.leave_subsection();
		
    typename DoFHandler<dim>::active_cell_iterator  
			cell = dof_handler.begin_active(), endc = dof_handler.end();
    
    //set boundary labels:
    // (1) : top
    // (2) : bottom
    // (3) : inlet
    // (4) : outlet
	for (; cell!=endc; ++cell)
	{
		for (int f=0; f<GeometryInfo<dim>::faces_per_cell; f++)
		{
			if (cell->face(f)->at_boundary())
			{
				if (fabs(cell->face(f)->center()[0]) < 1e-8)
				{
					cell->face(f)->set_boundary_indicator(3);
				}
				else if (fabs(cell->face(f)->center()[0] - domain_length) < 1e-8)	
				{
					cell->face(f)->set_boundary_indicator(4);
				}		
				else if (fabs(cell->face(f)->center()[1]) < 1e-8)
				{
					cell->face(f)->set_boundary_indicator(2);
				}
				else
				{
					cell->face(f)->set_boundary_indicator(1);
				}		
			}
		}
	}		
	
	//apply no-through boundary conditions on top and bottom
	//for RT elements the degrees of freedom on the faces are the normal
	//velocities, so we just apply 0 boundary conditions to these vectors
							
	const FEValuesExtractors::Vector VectorField(0); 
	ComponentMask vector_field_mask = fe.component_mask(VectorField);
	DoFTools::make_zero_boundary_constraints(dof_handler, 1, constraints, vector_field_mask);
	DoFTools::make_zero_boundary_constraints(dof_handler, 2, constraints, vector_field_mask);

	
	//mark boundary dofs based on initial position						
	if (iter == 0)
	{
		top_boundary_dofs.reserve(dof_handler.n_dofs());
		bottom_boundary_dofs.reserve(dof_handler.n_dofs());
		inlet_boundary_dofs.reserve(dof_handler.n_dofs());
		outlet_boundary_dofs.reserve(dof_handler.n_dofs());
		
		for (int i = 0; i < dof_handler.n_dofs(); i++)
		{
			top_boundary_dofs[i] = false;
			bottom_boundary_dofs[i] = false;
			inlet_boundary_dofs[i] = false;
			outlet_boundary_dofs[i] = false;
		}
		
		//add top left corner dof to top_boundary_dofs, by default this dof is on boundary 3
		cell = dof_handler.begin_active(), endc = dof_handler.end();
		std::vector<bool> vertex_touched (mesh.n_vertices(), false);
				
		for (; cell!=endc; ++cell)
		{
			for (int v=0; v<GeometryInfo<dim>::vertices_per_cell; v++)
			{
				if (vertex_touched[cell->vertex_index(v)] == false)
				{
					vertex_touched[cell->vertex_index(v)] == true;
					Point<dim> vertex = cell->vertex(v);			
					int dof_number = cell->vertex_dof_index(v, 0);
					
					if (fabs(vertex[1] - domain_height) < 1e-8)
					{					
						top_boundary_dofs[dof_number] = true;
					}
					
					if (fabs(vertex[1]) < 1e-8)
					{
						bottom_boundary_dofs[dof_number] = true;
					}
					
					if (fabs(vertex[0] - domain_length) < 1e-8)
					{
						outlet_boundary_dofs[dof_number] = true;
					}
					
					if (fabs(vertex[0]) < 1e-8)
					{
						inlet_boundary_dofs[dof_number] = true;
					}
					
				}
			}
		}
					
		//add porosity, stress/strain and displacement constraints
		solution.reinit(dof_handler.n_dofs());    
		calculate_initial_porosity();
	}
	
	DoFTools::extract_dofs(dof_handler, fe.component_mask(porosity), 
					porosity_dofs);
	DoFTools::extract_dofs(dof_handler, fe.component_mask(stress_xx), 
					stress_xx_dofs);
	DoFTools::extract_dofs(dof_handler, fe.component_mask(stress_yy), 
					stress_yy_dofs);
	DoFTools::extract_dofs(dof_handler, fe.component_mask(stress_xy), 
					stress_xy_dofs);
	DoFTools::extract_dofs(dof_handler, fe.component_mask(strain_xx), 
					strain_xx_dofs);
	DoFTools::extract_dofs(dof_handler, fe.component_mask(strain_yy), 
					strain_yy_dofs);
	DoFTools::extract_dofs(dof_handler, fe.component_mask(strain_xy), 
					strain_xy_dofs);				
	DoFTools::extract_dofs(dof_handler, fe.component_mask(disp_x), 
					disp_x_dofs);
	DoFTools::extract_dofs(dof_handler, fe.component_mask(disp_y), 
					disp_y_dofs);
	DoFTools::extract_dofs(dof_handler, fe.component_mask(gamma), 
					gamma_dofs);
	DoFTools::extract_dofs(dof_handler, fe.component_mask(stiffness), 
					stiffness_dofs);
											
	for (unsigned int i=0; i<dof_handler.n_dofs(); i++)
	{
		if (porosity_dofs[i])
		{
			constraints.add_line(i);	
			constraints.set_inhomogeneity(i, solution[i]);
		}
		else if (stress_xx_dofs[i] || stress_yy_dofs[i] || stress_xy_dofs[i] ||
				 strain_xx_dofs[i] || strain_yy_dofs[i] || strain_xy_dofs[i] ||
				 disp_x_dofs[i] || disp_y_dofs[i] || gamma_dofs[i] || stiffness_dofs[i])
		{
			constraints.add_line(i);
		}
	}	
				
    constraints.close();
            
	//calculate sparsity pattern 
	CompressedSparsityPattern c_sparsity(dof_handler.n_dofs(), dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern (dof_handler, c_sparsity);
    
    constraints.condense(c_sparsity);
	sparsity_pattern.copy_from(c_sparsity);
    system_matrix.reinit(sparsity_pattern);  
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
	double						permiability;
	double 						total_porosity = 0;
	double                      total_permiability = 0;
    int                         total_quad_points = 0;    
	
	prm.enter_subsection("Geometry");
	double domain_length = prm.get_double("length");
	prm.leave_subsection();
	
	std::vector<Tensor<2,dim> >  grad_phi_u (dofs_per_cell);
	std::vector<double>          div_phi_u (dofs_per_cell);
	std::vector<double>          phi_p (dofs_per_cell);
	std::vector<Tensor<1,dim> >  phi_u (dofs_per_cell);
	
	std::vector<double>   porosity_values(n_q_points);
	std::vector<Vector<double> > rhs_values (n_q_points, Vector<double>(dim+12)); 
	
	FullMatrix<double>          cell_matrix (dofs_per_cell, dofs_per_cell);
    Vector<double>              cell_rhs (dofs_per_cell);
    
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);    
	const FEValuesExtractors::Vector velocities (0);
	const FEValuesExtractors::Scalar pressure (dim);
	const FEValuesExtractors::Scalar porosity (dim + 1);
	
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
        
        fe_values[porosity].get_function_values(solution, porosity_values);
        cell_matrix = 0;
        cell_rhs  = 0;
        
        //calculate cell contribution to system         
        for (int q = 0; q < n_q_points; q++)
        {		
			for (int k=0; k<dofs_per_cell; k++)
			{
				div_phi_u[k]  = fe_values[velocities].divergence (k, q);
				phi_p[k]      = fe_values[pressure].value (k, q);
				phi_u[k]      = fe_values[velocities].value (k, q);				
			}
			
			permiability = calculate_permiability(porosity_values[q]);
			total_porosity += porosity_values[q];
			total_permiability += permiability;
			total_quad_points++;
			
			for (int i = 0; i < dofs_per_cell; i++)
			{		
				for (int j = 0; j < dofs_per_cell; j++)
				{
						
					cell_matrix(i,j) += 
							((mu/permiability)*phi_u[j]*phi_u[i]										
							- phi_p[j]*div_phi_u[i]
							- phi_p[i]*div_phi_u[j])
							*fe_values.JxW(q);									
				}				
								
				cell_rhs[i] += 
							(phi_u[i][0]*rhs_values[q][0] + phi_u[i][1]*rhs_values[q][1])//no longer dimension independent
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
							cell_rhs[i] += prm.get_double("r inlet")
									*phi_u[i][0]*fe_face_values.JxW(q_boundary);
						}
						else if (fabs(x - domain_length) < 1e-8)
						{
							cell_rhs[i] -= prm.get_double("r outlet")
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
		
		average_porosity = total_porosity/total_quad_points;
			
		printf("done (%gs)\n", timer()); 
		printf("average porosity = %f\n", average_porosity);
		printf("average permiability = %e\n", total_permiability/total_quad_points);
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
void DarcyFlow<dim>::move_mesh()
{
	Timer timer;
	if (verbose)
	{
		timer.start();
		printf("Moving mesh...");
	}
	
	std::vector<bool> vertex_touched(mesh.n_vertices(), false);
	typename DoFHandler<dim>::active_cell_iterator  
			cell = dof_handler.begin_active(), endc = dof_handler.end();
	
	double total_displacement = 0;
	int total_points = 0;
			
	for (; cell!=endc; ++cell)
	{
		for (int v=0; v<GeometryInfo<dim>::vertices_per_cell; v++)
		{			
			if (vertex_touched[cell->vertex_index(v)] == false)
			{				
				vertex_touched[cell->vertex_index(v)] = true;	
				Point<dim> vertex = cell->vertex(v);	
				double y_max = vertex[1];
						
				int dof_number = cell->vertex_dof_index(v, 0);					
				//dof numbers at nodes go:
				//0 - dim-1 : velocities
				//dim : pressure
				//dim+1 : porosity
				//dim+2 : stress_xx
				//dim+3 : stress_yy
				//dim+4 : stress_xy
				//dim+5 : strain_xx
				//dim+6 : strain_yy
				//dim+7 : strain_xy
				//dim+8 : disp_x
				//dim+9 : disp_y
				
				Point<dim> new_vertex;
				new_vertex[0] = vertex[0] + solution[dof_number + 8];
				new_vertex[1] = vertex[1] + solution[dof_number + 9];
				
				if (top_boundary_dofs[dof_number])
				{
					total_displacement += pow(solution[dof_number + 8], 2)
										+ pow(solution[dof_number + 9], 2);
					total_points++;
				}
				
				cell->vertex(v) = new_vertex;
			}
		}
	}			
	
	if (verbose)
	{
		timer.stop ();
		printf("done (%gs)\n",timer());
		printf("average displacement of top boundary = %e\n", sqrt(total_displacement)/total_points);
	}
}

template<int dim>
void DarcyFlow<dim>::calculate_displacement_x()
{
	Timer timer;
	double wtime_0;
	if (verbose)
	{
		timer.start();
		//wtime_0 = omp_get_wtime();
		printf("Calculating x displacement...");
	}
	
	prm.enter_subsection("Geometry");
	int nx = prm.get_integer("nx");
	int ny = prm.get_integer("ny");
	double length = prm.get_double("length");
	double initial_height = prm.get_double("initial height");		
	prm.leave_subsection();

	double hx = length/(nx - 1);
	double hy = initial_height/(ny - 1);
	
	//3 point Gauss quadrature
	double quad_pointsa[3] = {-sqrt(3.0/5.0), 0.0, sqrt(3.0/5.0)};
	double weightsa[3] = {5.0/9.0, 8.0/9.0, 5.0/9.0};
	std::vector<double> quad_points(&quad_pointsa[0], &quad_pointsa[0]+3);
	std::vector<double> weights(&weightsa[0], &weightsa[0]+3);
	
	std::vector<double> mapped_quad_points(3);
	std::vector<double> mapped_weights(3);
	
	std::vector<bool> vertex_touched(mesh.n_vertices(), false);
	typename DoFHandler<dim>::active_cell_iterator  
			cell, endc = dof_handler.end();				
	
	//#pragma omp parallel private(cell) shared(vertex_touched)
	for (cell = dof_handler.begin_active(); cell!=endc; ++cell)
	{
		for (int v=0; v<GeometryInfo<dim>::vertices_per_cell; v++)
		{			
			Point<dim> vertex = cell->vertex(v);	
			if (vertex_touched[cell->vertex_index(v)] == false)
			{				
				vertex_touched[cell->vertex_index(v)] = true;	
				
				double x_max = vertex[0];
				double y_max = vertex[1];
				
				int n_intervals_x = floor(x_max/hx);
				int n_intervals_y = floor(y_max/hy); 
						
				int dof_number = cell->vertex_dof_index(v, 0);					
				//dof numbers at nodes go:
				//0 - dim-1 : velocities
				//dim : pressure
				//dim+1 : porosity
				//dim+2 : stress_xx
				//dim+3 : stress_yy
				//dim+4 : stress_xy
				//dim+5 : strain_xx
				//dim+6 : strain_yy
				//dim+7 : strain_xy
				//dim+8 : disp_x
				//dim+9 : disp_y
				//dim+10 : gamma
				//dim+11: local_stiffness					
				
				for (int k = 0; k < n_intervals_x; k++)
				{
					double a = k*x_max/(n_intervals_x);
					double b = (k + 1)*x_max/(n_intervals_x);
					
					//map points and weights to domain [y_min_element, y_max_element]
					for (int i = 0; i < quad_points.size(); i++)
					{
						mapped_quad_points[i] = ((b - a)/2.0)*quad_points[i] + (a + b)/2.0;
						mapped_weights[i] = ((b - a)/2.0)*weights[i];
					}
					
					Vector<double> solution_values(dim+12); 
					
					for (int i = 0; i < quad_points.size(); i++)
					{
						Point<dim, double> p;
						p[0] = mapped_quad_points[i];
						p[1] = vertex[1];
			
						VectorTools::point_value(dof_handler, solution, p, solution_values);
						solution[dof_number + 7] += solution_values[dim + 5]*mapped_weights[i];	
					}				
				}
				
				double delta_y = y_max/n_intervals_y;
				solution[dof_number + 7] += calculate_c_1(vertex, delta_y);			
			}
		}
	}
	
	if (verbose)
	{
		timer.stop ();
		printf("done (%gs)\n",timer());
		//printf("done (%gs)\n",omp_get_wtime() - wtime_0);
	}	
}

template<int dim>
void DarcyFlow<dim>::calculate_displacement_y()
{
	Timer timer;
	//double wtime_0;
	if (verbose)
	{
		timer.start();
		//wtime_0 = omp_get_wtime();
		printf("Calculating y displacement...");
	}
	
	prm.enter_subsection("Geometry");
	int ny = prm.get_integer("ny");
	double initial_height = prm.get_double("initial height");
	prm.leave_subsection();
	
	double hy = initial_height/(ny - 1);
	
	//3 point Gauss quadrature
	double quad_pointsa[3] = {-sqrt(3.0/5.0), 0.0, sqrt(3.0/5.0)};
	double weightsa[3] = {5.0/9.0, 8.0/9.0, 5.0/9.0};
	std::vector<double> quad_points(&quad_pointsa[0], &quad_pointsa[0]+3);
	std::vector<double> weights(&weightsa[0], &weightsa[0]+3);
	
	std::vector<double> mapped_quad_points(3);
	std::vector<double> mapped_weights(3);

	std::vector<bool> vertex_touched(mesh.n_vertices(), false);
	typename DoFHandler<dim>::active_cell_iterator  
			cell, endc = dof_handler.end();				
	
	//#pragma omp parallel private(cell) shared(vertex_touched)
	for (cell = dof_handler.begin_active(); cell!=endc; ++cell)
	{
		for (int v=0; v<GeometryInfo<dim>::vertices_per_cell; v++)
		{			
			Point<dim> vertex = cell->vertex(v);	
			if (vertex_touched[cell->vertex_index(v)] == false && vertex[1] > 0)
			{				
				vertex_touched[cell->vertex_index(v)] = true;	
				
				double y_max = vertex[1];
				int n_intervals = floor(y_max/hy);
						
				int dof_number = cell->vertex_dof_index(v, 0);					
				//dof numbers at nodes go:
				//0 - dim-1 : velocities
				//dim : pressure
				//dim+1 : porosity
				//dim+2 : stress_xx
				//dim+3 : stress_yy
				//dim+4 : stress_xy
				//dim+5 : strain_xx
				//dim+6 : strain_yy
				//dim+7 : strain_xy
				//dim+8 : disp_x
				//dim+9 : disp_y
				//dim+10 : gamma
				//dim+11: local_stiffness					
				
				for (int k = 0; k < n_intervals; k++)
				{
					double a = k*y_max/(n_intervals);
					double b = (k + 1)*y_max/(n_intervals);
					
					//map points and weights to domain [y_min_element, y_max_element]
					for (int i = 0; i < quad_points.size(); i++)
					{
						mapped_quad_points[i] = ((b - a)/2.0)*quad_points[i] + (a + b)/2.0;
						mapped_weights[i] = ((b - a)/2.0)*weights[i];
					}
					
					Vector<double> solution_values(dim+12); 
					
					for (int i = 0; i < quad_points.size(); i++)
					{
						Point<dim, double> p;
						p[0] = vertex[0];
						p[1] = mapped_quad_points[i];
			
						VectorTools::point_value(dof_handler, solution, p, solution_values);
						solution[dof_number + 9] += solution_values[dim + 6]*mapped_weights[i];	
					}				
				}				
			}
		}
	}
	
	if (verbose)
	{
		timer.stop ();
		printf("done (%gs)\n",timer());
		//printf("done (%gs)\n",omp_get_wtime() - wtime_0);
	}
}

template<int dim>
double DarcyFlow<dim>::calculate_porosity_1d_integral(Point<dim> vertex)
{
	prm.enter_subsection("Geometry");
	int ny = prm.get_integer("ny");
	prm.leave_subsection();

	//3 point Gauss quadrature
	double quad_pointsa[3] = {-sqrt(3.0/5.0), 0.0, sqrt(3.0/5.0)};
	double weightsa[3] = {5.0/9.0, 8.0/9.0, 5.0/9.0};
	std::vector<double> quad_points(&quad_pointsa[0], &quad_pointsa[0]+3);
	std::vector<double> weights(&weightsa[0], &weightsa[0]+3);
	
	std::vector<double> mapped_quad_points(3);
	std::vector<double> mapped_weights(3);
	
	double y_max = vertex[1];
	double int_porosity = 0;
	
	for (int k = 0; k < ny - 1; k++)
	{
		double a = k*y_max/(ny - 1);
		double b = (k + 1)*y_max/(ny - 1);
		
		//map points and weights to domain [y_min_element, y_max_element]
		for (int i = 0; i < quad_points.size(); i++)
		{
			mapped_quad_points[i] = ((b - a)/2.0)*quad_points[i] + (a + b)/2.0;
			mapped_weights[i] = ((b - a)/2.0)*weights[i];
		}	

		Vector<double> solution_values(dim+12); 
		
		for (int i = 0; i < quad_points.size(); i++)
		{
			Point<dim, double> p;
			p[0] = vertex[0];
			p[1] = mapped_quad_points[i];
			
			VectorTools::point_value(dof_handler, solution, p, solution_values);
			
			int_porosity += solution_values[3]*mapped_weights[i];
		}
	}
	
	return int_porosity;
}

template<int dim>
void DarcyFlow<dim>::calculate_initial_porosity()
{	
	//srand(time(NULL));
	srand(123458);
			
	prm.enter_subsection("Geometry");
	double initial_height = prm.get_double("initial height");
	double length = prm.get_double("length");
	prm.leave_subsection();
	
	prm.enter_subsection("Constants");
	double r_in_min = prm.get_double("min inner pore radius");						
	double r_in_max = prm.get_double("max inner pore radius");
	double r_out_min = prm.get_double("min outer pore radius");						
	double r_out_max = prm.get_double("max outer pore radius");
	double phi_min = prm.get_double("min porosity");
	double phi_max = prm.get_double("max porosity");					
	int N =  prm.get_integer("number of pores");
	prm.leave_subsection();
	
	double x[N], y[N], r_in[N], r_out[N];	
	double min_porosity_level = 0;
	
	while (min_porosity_level < 0.1)
	{
		min_porosity_level = 1;		
			
		for (int i = 0; i < N; i++)
		{
			x[i] = length*rand()/((double) RAND_MAX);
			y[i] = initial_height*rand()/((double) RAND_MAX);
			r_in[i] = r_in_min + (r_in_max - r_in_min)*rand()/((double) RAND_MAX);
			r_out[i] = r_out_min + (r_out_max - r_out_min)*rand()/((double) RAND_MAX);		
		}
		
		std::vector<bool> vertex_touched (mesh.n_vertices(), false);		
		typename DoFHandler<dim>::active_cell_iterator  cell, endc = dof_handler.end();				

		//#pragma omp parallel private(cell) shared(vertex_touched)
		for (cell = dof_handler.begin_active(); cell!=endc; ++cell)
		{
			for (int v=0; v<GeometryInfo<dim>::vertices_per_cell; v++)
			{
				if (vertex_touched[cell->vertex_index(v)] == false)
				{				
					vertex_touched[cell->vertex_index(v)] = true;	
					Point<dim> vertex = cell->vertex(v);			
					int dof_number = cell->vertex_dof_index(v, 0);
					
					/*if ((vertex[1] < 0.1*initial_height) || (vertex[1] > 0.9*initial_height))
					{
						 solution[dof_number + dim + 1] = phi_min;
					}
					else*/
					{
						double tmp = 0;
						for (int i = 0; i < N; i++)
						{
							double dist = pow(vertex[0] - x[i],2) + pow(vertex[1] - y[i],2);
							//tmp += exp(-dist/pow(0.05,2));
							
							tmp += std::max((dist > r_in[i]) 
										? (dist - r_out[i])/(r_in[i] - r_out[i]) : 1.0, 0.0);
						}
						
						solution[dof_number] = std::min(std::max(tmp, phi_min), phi_max);
					}
				}
			}
		}
		
		//check that total (integrated) porosity is above some minimum for all vertical cross 
		//sections If not, generate new initial porosity field
		cell = dof_handler.begin_active(), endc = dof_handler.end();
		
		for (int i = 0; i < mesh.n_vertices(); i++)
		{
			vertex_touched[i] = false;
		}
		
		for (; cell!=endc; ++cell)
		{
			for (int v=0; v<GeometryInfo<dim>::vertices_per_cell; v++)
			{
				if (vertex_touched[cell->vertex_index(v)] == false)
				{
					Point<dim> vertex = cell->vertex(v);
					vertex_touched[cell->vertex_index(v)] = true;
					
					//if (cell->at_boundary(v))
					//{
						int dof_number = cell->vertex_dof_index(v, 0);
						
						if (top_boundary_dofs[dof_number])
						{
							double int_porosity = calculate_porosity_1d_integral(vertex);
							
							min_porosity_level = (int_porosity < min_porosity_level) 
													? int_porosity : min_porosity_level;
						}
					//}
				}
			}
		}
		
		printf("minimum porosity = %e\n", min_porosity_level);
	}
}

template<int dim>
double DarcyFlow<dim>::calculate_permiability(double porosity)
{
	prm.enter_subsection("Constants");
	double a = prm.get_double("permiability coefficent a");
	double b = prm.get_double("permiability coefficent b");
	prm.leave_subsection();
	
	return a*pow(porosity, b);
}

template<int dim>
void DarcyFlow<dim>::calculate_stress()
{
	Timer timer;
	//double wtime_0;
	if (verbose)
	{
		timer.start();
		//wtime_0 = omp_get_wtime();
		printf("Calculating stress...");
	}
	
	const double EPSILON = 1e-8;
	const double p_0 = 0.980655e5; //initial pressure, 1 atm = 0.980655e5 Pa
	double 	mu = prm.get_double("mu");
	double 	E = prm.get_double("Young's modulus");
	double 	nu = prm.get_double("Poisson's ratio");
	
	QGauss<dim> midpoint_quad(1);
	
	std::vector<Tensor< 2, dim> > velocity_gradient(midpoint_quad.size());
	std::vector<double> pressure_value(midpoint_quad.size());
	
	const FEValuesExtractors::Vector velocities(0);	
	const FEValuesExtractors::Scalar pressure(dim);	
	        
	std::vector<bool> vertex_touched(mesh.n_vertices(), false);
	typename DoFHandler<dim>::active_cell_iterator  
			cell, endc = dof_handler.end();		
				
	FEValues<dim> fe_values(fe, midpoint_quad, update_values | update_gradients | update_quadrature_points);
	
	//#pragma omp parallel private(cell, velocity_gradient) shared(vertex_touched)
	{		
		for (cell = dof_handler.begin_active(); cell!=endc; ++cell)
		{		
			fe_values.reinit(cell);
			fe_values[velocities].get_function_gradients(solution, velocity_gradient);
			fe_values[pressure].get_function_values(solution, pressure_value);
			
			for (int v=0; v<GeometryInfo<dim>::vertices_per_cell; v++)
			{		
				int dof_number = cell->vertex_dof_index(v, 0);
				const double BASE_PRESSURE = 0.980655e5;//initial pressure, 1 atm = 0.980655e5 Pa
			
				//double pressure = solution[dof_number] - BASE_PRESSURE;	
				double pressure = pressure_value[0] - BASE_PRESSURE;	
				
				//find scaling factor:
				int num_cells = 4;
				
				//corner vertex attached to 1 cell
				if ((inlet_boundary_dofs[dof_number] && top_boundary_dofs[dof_number])
					|| (inlet_boundary_dofs[dof_number] && bottom_boundary_dofs[dof_number])
					|| (outlet_boundary_dofs[dof_number] && top_boundary_dofs[dof_number])
					|| (outlet_boundary_dofs[dof_number] && bottom_boundary_dofs[dof_number]))
				{
					num_cells = 1;
				}
				//boundary vertex attached to 2 cells
				else if (inlet_boundary_dofs[dof_number] || bottom_boundary_dofs[dof_number]
							||outlet_boundary_dofs[dof_number] || top_boundary_dofs[dof_number])
				{
					num_cells = 2;
				}				
				
				solution[dof_number + 1] += (-pressure + 2*mu*velocity_gradient[0][0][0])
												   /num_cells;
				solution[dof_number + 2] += (-pressure + 2*mu*velocity_gradient[0][1][1])
												   /num_cells;
				solution[dof_number + 3] += (mu*(velocity_gradient[0][0][1] 
												  + velocity_gradient[0][1][0]))
												  /num_cells;				
			}
		}
	}
	
	FEValuesExtractors::Scalar strain_xx(dim+5);
	std::vector<bool> strain_xx_dofs(dof_handler.n_dofs(), false);
	DoFTools::extract_dofs(dof_handler, fe.component_mask(strain_xx), strain_xx_dofs);
	
	for (int i = 0; i < dof_handler.n_dofs(); i++)
	{
		//dof numbers at nodes go:
		//0 : porosity
		//1 : stress_xx
		//2 : stress_yy
		//3 : stress_xy
		//4 : strain_xx
		//5 : strain_yy
		//6 : strain_xy
		//7 : disp_x
		//8 : disp_y
		//9 : gamma
		//10: local_stiffness	
		
		if (strain_xx_dofs[i])
		{
			double porosity = solution[i - 4];
			double E_local = E*(1 - porosity)/(1 - average_porosity);
			
			double old_stress_xx;
			double old_stress_yy;
			double old_stress_xy;
				
			if (iter == 0)
				{
					old_stress_xx = 0;
					old_stress_yy = 0;
					old_stress_xy = 0;
				}
				else
				{
					old_stress_xx = previous_solution[i - 3];
					old_stress_yy = previous_solution[i - 2];
					old_stress_xy = previous_solution[i - 1];
				}
				
			solution[i] = ((solution[i - 3] - nu*solution[i - 2]) - (old_stress_xx -nu*old_stress_yy))/E_local;
			solution[i + 1] = ((solution[i - 2] - nu*solution[i - 3]) - (old_stress_yy - nu*old_stress_xx))/E_local;
			solution[i + 2] = ((1 + nu)*solution[i - 1] - (1 + nu)*old_stress_xy)/E_local;	
			
			solution[i + 6] = E_local;		
		}
	}
	
	if (verbose)
	{
		timer.stop ();
		printf("done (%gs)\n",timer());
		//printf("done (%gs)\n",omp_get_wtime() - wtime_0);
	}
}

template<int dim>
void DarcyFlow<dim>::calculate_gamma()
{	
	Timer timer;
	double wtime_0;
	if (verbose)
	{
		//timer.start();
		wtime_0 = omp_get_wtime();
		printf("Calculating gamma...");
	}
	
	prm.enter_subsection("Geometry");
	double length = prm.get_double("length");
	double initial_height = prm.get_double("initial height");	
	prm.leave_subsection();
	
	double hx = length/(nx - 1);
	double hy = initial_height/(ny - 1);
	
	//3 point Gauss quadrature
	double quad_pointsa[3] = {-sqrt(3.0/5.0), 0.0, sqrt(3.0/5.0)};
	double weightsa[3] = {5.0/9.0, 8.0/9.0, 5.0/9.0};
	std::vector<double> quad_points(&quad_pointsa[0], &quad_pointsa[0]+3);
	std::vector<double> weights(&weightsa[0], &weightsa[0]+3);
	
	std::vector<bool> vertex_touched (mesh.n_vertices(), false);
	typename DoFHandler<dim>::active_cell_iterator  
			cell, endc = dof_handler.end();				

	//#pragma omp parallel private(cell) shared(vertex_touched)
	for (cell = dof_handler.begin_active(); cell!=endc; ++cell)
	{
		for (int v=0; v<GeometryInfo<dim>::vertices_per_cell; v++)
		{			
			if (vertex_touched[cell->vertex_index(v)] == false)
			{				
				vertex_touched[cell->vertex_index(v)] = true;	
				Point<dim> vertex = cell->vertex(v);			
				int dof_number = cell->vertex_dof_index(v, 0);
				
				double x_max = vertex[0];
				double y_max = vertex[1];
				double int_y = 0;
				double int_x = 0;
				
				//dof numbers at nodes go:
				//0 - dim-1 : velocities
				//dim : pressure
				//dim+1 : porosity
				//dim+2 : stress_xx
				//dim+3 : stress_yy
				//dim+4 : stress_xy
				//dim+5 : strain_xx
				//dim+6 : strain_yy
				//dim+7 : strain_xy
				//dim+8 : disp_x
				//dim+9 : disp_y
				//dim+10: gamma	
				//dim+11: local_stiffness					
				
				if (inlet_boundary_dofs[dof_number] && bottom_boundary_dofs[dof_number])
				{
					int_x = 0;
					int_y = 0;
				}
				else if (inlet_boundary_dofs[dof_number])
				{
					int_x = 0;
					
					//calculate y integral of partial(strain_yy)/partial(x)
					int n_intervals_y = floor(y_max/hy);
					for (int k = 0; k < n_intervals_y; k++)
					{
						double a = k*y_max/(n_intervals_y);
						double b = (k + 1)*y_max/(n_intervals_y);
						
						//map points and weights to domain [a,b]
						std::vector<double> mapped_quad_points(3);
						std::vector<double> mapped_weights(3);
						
						for (int i = 0; i < quad_points.size(); i++)
						{
							mapped_quad_points[i] = ((b - a)/2.0)*quad_points[i] + (a + b)/2.0;
							mapped_weights[i] = ((b - a)/2.0)*weights[i];
						}	

						double partial_x, partial_y;						
						for (int i = 0; i < quad_points.size(); i++)
						{
							Point<dim, double> p;
							p[0] = vertex[0];
							p[1] = mapped_quad_points[i];
							
							calculate_derivatives(p, dim+6, partial_x, partial_y, dof_number);
							
							int_y += partial_x*mapped_weights[i];									
						}
					}
				}
				else if (bottom_boundary_dofs[dof_number])
				{
					int_y = 0;
					
					//calculate x integral of partial(strain_xx)/partial(y)
					int n_intervals_x = floor(x_max/hx);
					for (int k = 0; k < n_intervals_x; k++)
					{
						double a = k*x_max/(n_intervals_x);
						double b = (k + 1)*x_max/(n_intervals_x);
						
						//map points and weights to domain [a,b]
						std::vector<double> mapped_quad_points(3);
						std::vector<double> mapped_weights(3);
						
						for (int i = 0; i < quad_points.size(); i++)
						{
							mapped_quad_points[i] = ((b - a)/2.0)*quad_points[i] + (a + b)/2.0;
							mapped_weights[i] = ((b - a)/2.0)*weights[i];
						}	
						
						double partial_x, partial_y;						
						for (int i = 0; i < quad_points.size(); i++)
						{
							Point<dim, double> p;
							p[0] = mapped_quad_points[i];
							p[1] = vertex[1];							
							
							calculate_derivatives(p, dim+5, partial_x, partial_y, dof_number);
							
							int_x += partial_y*mapped_weights[i];									
						}
					}
				}
				else
				{				
					//calculate y integral of partial(strain_yy)/partial(x)
					int n_intervals_y = floor(y_max/hy);
					for (int k = 0; k < n_intervals_y; k++)
					{
						double a = k*y_max/(n_intervals_y);
						double b = (k + 1)*y_max/(n_intervals_y);
						
						//map points and weights to domain [a,b]
						std::vector<double> mapped_quad_points(3);
						std::vector<double> mapped_weights(3);
						
						for (int i = 0; i < quad_points.size(); i++)
						{
							mapped_quad_points[i] = ((b - a)/2.0)*quad_points[i] + (a + b)/2.0;
							mapped_weights[i] = ((b - a)/2.0)*weights[i];
						}	

						double partial_x, partial_y;						
						for (int i = 0; i < quad_points.size(); i++)
						{
							Point<dim, double> p;
							p[0] = vertex[0];
							p[1] = mapped_quad_points[i];
							
							calculate_derivatives(p, dim+6, partial_x, partial_y, dof_number);
							
							int_y += partial_x*mapped_weights[i];									
						}
					}
					
					//calculate x integral of partial(strain_xx)/partial(y)
					int n_intervals_x = floor(x_max/hx);
					for (int k = 0; k < n_intervals_x; k++)
					{
						double a = k*x_max/(n_intervals_x);
						double b = (k + 1)*x_max/(n_intervals_x);
						
						//map points and weights to domain [a,b]
						std::vector<double> mapped_quad_points(3);
						std::vector<double> mapped_weights(3);
						
						for (int i = 0; i < quad_points.size(); i++)
						{
							mapped_quad_points[i] = ((b - a)/2.0)*quad_points[i] + (a + b)/2.0;
							mapped_weights[i] = ((b - a)/2.0)*weights[i];
						}	
						
						double partial_x, partial_y;						
						for (int i = 0; i < quad_points.size(); i++)
						{
							Point<dim, double> p;
							p[0] = mapped_quad_points[i];
							p[1] = vertex[1];							
							
							calculate_derivatives(p, dim+5, partial_x, partial_y, dof_number);
							
							int_x += partial_y*mapped_weights[i];									
						}
					}
				}
				
				//gamma = 2*strain_xy - (int_x + int_y)
				solution[dof_number + 9] = 2*solution[dof_number + 6] - (int_x + int_y);				
			}
		}
	}
	
	if (verbose)
	{
		//timer.stop ();
		//printf("done (%gs)\n",timer());
		printf("done (%gs)\n", omp_get_wtime() - wtime_0);
	}
}

template<int dim>
double DarcyFlow<dim>::calculate_c_1(Point<dim> p, double delta_y)
{
	//Taken largely from tutorial 52
	
	double y = 0;
	double y_max = p[1];
	Vector<double> c_1;
	const double coarsen_param = 1.2;
	const double refine_param = 0.8;
	const double min_delta = 1e-8;
	const double max_delta = 10*delta_y;
	const double refine_tol = 1e-1;
	const double coarsen_tol = 1e-5;
	
	c_1.reinit(1);
	c_1 = 0.0;
	TimeStepping::EmbeddedExplicitRungeKutta<Vector<double> > embedded_explicit_runge_kutta(
					TimeStepping::HEUN_EULER, coarsen_param, refine_param, min_delta, max_delta, 
					refine_tol,	coarsen_tol);
					
	while (y < y_max)
	{
		if (y + delta_y > y_max)
		{
			delta_y = y_max - y;
		}
		
		y = embedded_explicit_runge_kutta.evolve_one_time_step(
							std_cxx11::bind(&DarcyFlow::get_gamma,this,std_cxx11::_1,std_cxx11::_2),
							y, delta_y, c_1);
							
		delta_y = embedded_explicit_runge_kutta.get_status().delta_t_guess;
	}

	return c_1[0];
}

template<int dim>
Vector<double> DarcyFlow<dim>::get_gamma(const double y, const Vector<double> &c_1) const
{
	Point<dim> p;
	p[0] = 0;//gamma should be independent of x
	p[1] = y;
	
	Vector<double> gamma;
	gamma.reinit(1);
	
	Vector<double> solution_values(dim+12); 
	VectorTools::point_value(dof_handler, solution, p, solution_values);
	
	gamma[0] = solution_values[dim+10];
	
	return gamma;
}

template<int dim>
void DarcyFlow<dim>::calculate_derivatives(Point<dim> p, int component, 
			double &partial_x, double &partial_y, int dof_number)
{
	const double EPSILON = 1e-8;
	
	Vector<double> solution_values_left(dim+12); 	
	Vector<double> solution_values_right(dim+12); 
	Vector<double> solution_values_bottom(dim+12); 
	Vector<double> solution_values_top(dim+12); 
	Vector<double> solution_values_self(dim+12);
					
	Point<dim> p_left, p_right, p_bottom, p_top;	

	//several boundary cases to condsider:
	//1. top left corner
	if (inlet_boundary_dofs[dof_number] && top_boundary_dofs[dof_number])
	{
		p_right[0] = p[0] + EPSILON/2.0;
		p_right[1] = p[1];
		
		p_bottom[0] = p[0];
		p_bottom[1] = p[1] - EPSILON/2.0;
		
		VectorTools::point_value(dof_handler, solution, p_right, solution_values_right);
		VectorTools::point_value(dof_handler, solution, p_bottom, solution_values_bottom);
		VectorTools::point_value(dof_handler, solution, p, solution_values_self);
		
		partial_x = (solution_values_right[component] 
					- solution_values_self[component])/(EPSILON/2.0);
		partial_y = (solution_values_self[component] 
					- solution_values_bottom[component])/(EPSILON/2.0);
	}
	//2. top right corner
	else if (outlet_boundary_dofs[dof_number] && top_boundary_dofs[dof_number])
	{
		p_left[0] = p[0] - EPSILON/2.0;
		p_left[1] = p[1];
		
		p_bottom[0] = p[0];
		p_bottom[1] = p[1] - EPSILON/2.0;
		
		VectorTools::point_value(dof_handler, solution, p_left, solution_values_left);
		VectorTools::point_value(dof_handler, solution, p_bottom, solution_values_bottom);
		VectorTools::point_value(dof_handler, solution, p, solution_values_self);
		
		partial_x = (solution_values_self[component] 
					- solution_values_left[component])/(EPSILON/2.0);
		partial_y = (solution_values_self[component] 
					- solution_values_bottom[component])/(EPSILON/2.0);
	}
	//3. bottom right corner
	else if (outlet_boundary_dofs[dof_number] && bottom_boundary_dofs[dof_number])
	{
		p_left[0] = p[0] - EPSILON/2.0;
		p_left[1] = p[1];
		
		p_top[0] = p[0];
		p_top[1] = p[1] + EPSILON/2.0;
		
		VectorTools::point_value(dof_handler, solution, p_left, solution_values_left);
		VectorTools::point_value(dof_handler, solution, p_top, solution_values_top);
		VectorTools::point_value(dof_handler, solution, p, solution_values_self);
		
		partial_x = (solution_values_self[component] 
					- solution_values_left[component])/(EPSILON/2.0);
		partial_y = (solution_values_top[component] 
					- solution_values_self[component])/(EPSILON/2.0);
	}
	//4. bottom left corner
	else if (inlet_boundary_dofs[dof_number] && bottom_boundary_dofs[dof_number])
	{
		p_right[0] = p[0] + EPSILON/2.0;
		p_right[1] = p[1];
		
		p_top[0] = p[0];
		p_top[1] = p[1] + EPSILON/2.0;
		
		VectorTools::point_value(dof_handler, solution, p_right, solution_values_right);
		VectorTools::point_value(dof_handler, solution, p_top, solution_values_top);
		VectorTools::point_value(dof_handler, solution, p, solution_values_self);
		
		partial_x = (solution_values_right[component] 
					- solution[dof_number + component])/(EPSILON/2.0);
		partial_y = (solution_values_top[component] 
					- solution[dof_number + component])/(EPSILON/2.0);
	}
	//5. top bounday 
	else if (top_boundary_dofs[dof_number])
	{
		p_right[0] = p[0] + EPSILON;
		p_right[1] = p[1];
		
		p_left[0] = p[0] - EPSILON;
		p_left[1] = p[1];
		
		p_bottom[0] = p[0];
		p_bottom[1] = p[1] - EPSILON/2.0;
		
		VectorTools::point_value(dof_handler, solution, p_right, solution_values_right);
		VectorTools::point_value(dof_handler, solution, p_left, solution_values_left);
		VectorTools::point_value(dof_handler, solution, p_bottom, solution_values_bottom);
		VectorTools::point_value(dof_handler, solution, p, solution_values_self);
		
		partial_x = (solution_values_right[component] 
					- solution_values_left[component])/(2*EPSILON);
		partial_y = (solution_values_self[component] 
					- solution_values_bottom[component])/(EPSILON/2.0);
	}
	//6. outlet boundary
	else if (outlet_boundary_dofs[dof_number])
	{						
		p_left[0] = p[0] - EPSILON/2.0;
		p_left[1] = p[1];
		
		p_top[0] = p[0];
		p_top[1] = p[1] + EPSILON;
		
		p_bottom[0] = p[0];
		p_bottom[1] = p[1] - EPSILON;
		
		VectorTools::point_value(dof_handler, solution, p_bottom, solution_values_bottom);
		VectorTools::point_value(dof_handler, solution, p_left, solution_values_left);
		VectorTools::point_value(dof_handler, solution, p_top, solution_values_top);
		VectorTools::point_value(dof_handler, solution, p, solution_values_self);
		
		partial_x = (solution_values_self[component] 
					- solution_values_left[component])/(EPSILON/2.0);
		partial_y = (solution_values_top[component] 
					- solution_values_bottom[component])/(2*EPSILON);
	}
	//7. bottom boundary
	else if (bottom_boundary_dofs[dof_number])
	{		
		p_right[0] = p[0] + EPSILON;
		p_right[1] = p[1];
		
		p_left[0] = p[0] - EPSILON;
		p_left[1] = p[1];
		
		p_top[0] = p[0];
		p_top[1] = p[1] + EPSILON/2.0;
		
		VectorTools::point_value(dof_handler, solution, p_right, solution_values_right);
		VectorTools::point_value(dof_handler, solution, p_left, solution_values_left);
		VectorTools::point_value(dof_handler, solution, p_top, solution_values_top);
		VectorTools::point_value(dof_handler, solution, p, solution_values_self);
		
		partial_x = (solution_values_right[component] 
					- solution_values_left[component])/(2*EPSILON);
		partial_y = (solution_values_top[component] 
					- solution_values_self[component])/(EPSILON/2.0);
	}
	//8. inlet boundary
	else if (inlet_boundary_dofs[dof_number])
	{						
		p_right[0] = p[0] + EPSILON/2.0;
		p_right[1] = p[1];
		
		p_top[0] = p[0];
		p_top[1] = p[1] + EPSILON;
		
		p_bottom[0] = p[0];
		p_bottom[1] = p[1] - EPSILON;
		
		VectorTools::point_value(dof_handler, solution, p_bottom, solution_values_bottom);
		VectorTools::point_value(dof_handler, solution, p_right, solution_values_right);
		VectorTools::point_value(dof_handler, solution, p_top, solution_values_top);
		VectorTools::point_value(dof_handler, solution, p, solution_values_self);
		
		partial_x = (solution_values_right[component] 
					- solution_values_self[component])/(EPSILON/2.0);
		partial_y = (solution_values_top[component] 
					- solution_values_bottom[component])/(2*EPSILON);
	}
	//9. interior p
	else
	{					
		p_left[0] = p[0] - EPSILON;
		p_left[1] = p[1];						
		
		p_right[0] = p[0] + EPSILON;
		p_right[1] = p[1];						
		
		p_bottom[0] = p[0];
		p_bottom[1] = p[1] - EPSILON;						
		
		p_top[0] = p[0];
		p_top[1] = p[1] + EPSILON;			
				
		VectorTools::point_value(dof_handler, solution, p_left, solution_values_left);
		VectorTools::point_value(dof_handler, solution, p_right, solution_values_right);
		VectorTools::point_value(dof_handler, solution, p_bottom, solution_values_bottom);
		VectorTools::point_value(dof_handler, solution, p_top, solution_values_top);
		
		partial_x = (solution_values_right[component] 
					- solution_values_left[component])/(2*EPSILON);
		partial_y = (solution_values_top[component] 
					- solution_values_bottom[component])/(2*EPSILON);
	}
}

template<int dim>
void DarcyFlow<dim>::output_results(int iter)
{	
	std::vector<std::string> solution_names;
	solution_names.push_back("u");
	solution_names.push_back("v");
	solution_names.push_back("pressure");
	solution_names.push_back("porosity");
	solution_names.push_back("stress_xx");
	solution_names.push_back("stress_yy");
	solution_names.push_back("stress_xy");
	solution_names.push_back("strain_xx");
	solution_names.push_back("strain_yy");
	solution_names.push_back("strain_xy");
	solution_names.push_back("disp_x");
	solution_names.push_back("disp_y");
	solution_names.push_back("gamma");
	solution_names.push_back("local_stiffness");
	
	DataOut<2> data_out;
	data_out.attach_dof_handler(dof_handler);
	data_out.add_data_vector(solution, solution_names);
			
	data_out.build_patches();
	std::ostringstream filename;
	
	prm.enter_subsection("Output options");
	filename << prm.get("output name") << iter << ".vtk";
	prm.leave_subsection();
	
	std::ofstream output (filename.str().c_str());
	data_out.write_vtk (output);
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
	
	iter = 0;
	for (; iter < 3; iter++)
	{
		printf("iter %i\n", iter);
		
		//mesh.clear();			
		
		/*if (iter!=0)
		{
			nx*=2;
			ny*=2;
		}*/
				
		setup_geometry();		
												
		assemble();
		solve();
		
		calculate_stress();	
		//calculate_gamma();	
		//calculate_displacement_x();
		calculate_displacement_y();
		
		output_results(iter);		
		previous_solution = solution;
		
		move_mesh();			
		printf("\n");
				
		//printf("residual = %e for iter %i\n", mean_displacement, iter);
		
		/*if (mean_displacement < 1e-7)
		{
			break;
		}*/
	}
	
	output_results(iter);	
}

int main ()
{
	const int dim = 2;	
	//omp_set_num_threads(1);
	
	std::vector<const FiniteElement<dim> *> fes;
	std::vector<unsigned int> multiplicities;	
	
	FE_RaviartThomas<dim> velocities = FE_RaviartThomas<dim>(1);
	FE_DGQ<dim> pressure = FE_DGQ<dim>(1);
	FE_Q<dim> porosity = FE_Q<dim>(1);
	FE_Q<dim> stress_xx = FE_Q<dim>(1);
	FE_Q<dim> stress_yy = FE_Q<dim>(1);
	FE_Q<dim> stress_xy = FE_Q<dim>(1);
	
	FE_Q<dim> strain_xx = FE_Q<dim>(1);
	FE_Q<dim> strain_yy = FE_Q<dim>(1);
	FE_Q<dim> strain_xy = FE_Q<dim>(1);
	
	FE_Q<dim> disp_x = FE_Q<dim>(1);
	FE_Q<dim> disp_y = FE_Q<dim>(1);
	
	FE_Q<dim> gamma = FE_Q<dim>(1);
	FE_Q<dim> local_stiffness = FE_Q<dim>(1);
	
	fes.push_back(&velocities);
	fes.push_back(&pressure);
	fes.push_back(&porosity);
	fes.push_back(&stress_xx);
	fes.push_back(&stress_yy);
	fes.push_back(&stress_xy);
	
	fes.push_back(&strain_xx);
	fes.push_back(&strain_yy);
	fes.push_back(&strain_xy);
	
	fes.push_back(&disp_x);
	fes.push_back(&disp_y);
	
	fes.push_back(&gamma);
	fes.push_back(&local_stiffness);
	
	multiplicities.push_back(1);
	multiplicities.push_back(1);
	multiplicities.push_back(1);
	multiplicities.push_back(1);
	multiplicities.push_back(1);
	multiplicities.push_back(1);	
	multiplicities.push_back(1);
	multiplicities.push_back(1);
	multiplicities.push_back(1);
	multiplicities.push_back(1);
	multiplicities.push_back(1);
	multiplicities.push_back(1);
	multiplicities.push_back(1);
	
	DarcyFlow<dim> darcy(fes, multiplicities);
	
	MyZeroFunction<dim>  ff;
	MyZeroFunction<dim>  bv;
	
	darcy.forcing_function   = &ff;
	darcy.boundary_values    = &bv;
    darcy.run();
}
