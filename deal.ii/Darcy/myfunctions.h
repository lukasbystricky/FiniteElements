#ifndef MYFUNCTIONS_H
#define MYFUNCTIONS_H

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/base/point.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/convergence_table.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/solution_transfer.h>


#include <fstream>
#include <cmath>

using namespace dealii;

/*******************************************************************************
 * Define zero function
 ******************************************************************************/
template<int dim>
class MyZeroFunction : public Function<dim>
{
    public:
        MyZeroFunction() : Function<dim>(3) {};       
        virtual void vector_value(const Point<dim> &p, 
				Vector<double> &values) const;
};

template<int dim>
void MyZeroFunction<dim>::vector_value(const Point<dim> &p, 
                Vector<double> &values) const
{
	values(0) = 0;
	values(1) = 0;
	values(2) = 0;
}


/*******************************************************************************
 * Define boundary conditions, forcing function and 
 * exact solution for convergence study
 ******************************************************************************/
template<int dim>
class ExactSolutionBoundaryValues : public Function<dim>
{
    public:
        ExactSolutionBoundaryValues() : Function<dim>(3) {}
        virtual void vector_value(const Point<dim> &p, 
                Vector<double> &values) const;
};

template<int dim>
void ExactSolutionBoundaryValues<dim>::vector_value(const Point<dim> &p, 
                Vector<double> &values) const
{  	
	const double x = p[0];
	const double y = p[1];
		
	values(0) = -cos(M_PI*x)/M_PI;
	values(1) = -y*sin(M_PI*x);
	values(2) = 0; //no need to specify boundary conditions for the pressure
}

template<int dim>
class ExactSolutionForcingFunction : public Function<dim>
{
    public:
        ExactSolutionForcingFunction() : Function<dim>(3) {};       
        virtual void vector_value(const Point<dim> &p, 
				Vector<double> &values) const;
};

template<int dim>
void ExactSolutionForcingFunction<dim>::vector_value(const Point<dim> &p, 
                Vector<double> &values) const
{
	const double x = p[0];
	const double y = p[1];
	
	values(0) = cos(M_PI*x)*(M_PI + 1/M_PI) + exp(-x)*y;
	values(1) = sin(M_PI*x)*y*(pow(M_PI,2) + 1) - exp(-x);
	values(2) = 0; //incompressible flow
}

template<int dim>
class ExactSolution : public Function<dim>
{
    public:
        ExactSolution() : Function<dim>(3) {}
        virtual void vector_value(const Point<dim> &p, 
                Vector<double> &values) const;
        virtual void vector_gradient(const Point<dim> &p, 
                std::vector<Tensor<1,dim> > &gradients) const;
};

template<int dim>
void ExactSolution<dim>::vector_value(const Point<dim> &p, 
                Vector<double> &values) const
{
	const double x = p[0];
	const double y = p[1];
		
	values(0) = -cos(M_PI*x)/M_PI;
	values(1) = -y*sin(M_PI*x);
	values(2) = exp(-x)*y; 
}

template<int dim>
void ExactSolution<dim>::vector_gradient(const Point<dim> &p, 
                std::vector<Tensor<1,dim> > & gradients) const
{
	const double x = p[0];
	const double y = p[1];
	
    gradients[0][0] = sin(M_PI*x);
    gradients[0][1] = 0;
    
    gradients[1][0] = -y*M_PI*cos(M_PI*x);
    gradients[1][1] = -sin(M_PI*x);
    
    gradients[2][0] = -y*exp(-x);
    gradients[2][1] = exp(-x);
}

template<int dim>
class ExactSolutionPressureDriven : public Function<dim>
{
    public:
        ExactSolutionPressureDriven() : Function<dim>(3) {}
        virtual void vector_value(const Point<dim> &p, 
                Vector<double> &values) const;
        virtual void vector_gradient(const Point<dim> &p, 
                std::vector<Tensor<1,dim> > &gradients) const;
};

template<int dim>
void ExactSolutionPressureDriven <dim>::vector_value(const Point<dim> &p, 
                Vector<double> &values) const
{
	const double x = p[0];
	const double y = p[1];
	
	values(0) = 1 - (double) cosh(y - 0.5)/((double) cosh(0.5));
	values(1) = 0;
	values(2) = 1-x; 
}

template<int dim>
void ExactSolutionPressureDriven <dim>::vector_gradient(const Point<dim> &p, 
                std::vector<Tensor<1,dim> > & gradients) const
{
	const double x = p[0];
	const double y = p[1];
	
    gradients[0][0] = 0;
    gradients[0][1] = -(double) sinh(y - 0.5)/((double) cosh(0.5));
    
    gradients[1][0] = 0;
    gradients[1][1] = 0;
    
    gradients[2][0] = -1;
    gradients[2][1] = 0;
}

#endif
