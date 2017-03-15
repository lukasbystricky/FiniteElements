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
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/base/point.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/convergence_table.h>
#include <deal.II/base/parameter_handler.h>
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
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/grid/tria_boundary_lib.h>

#include <fstream>
#include <cmath>


using namespace dealii;

/*******************************************************************************
 * Define boundary conditions, initial conditions, forcing function and 
 * exact solution for convergence study
 ******************************************************************************/
template<int dim>
class ExactSolution : public Function<dim>
{
    public:
        ExactSolution() : Function<dim>() {}
        virtual double value(const Point<dim> &p, const unsigned int component = 0) const;
        virtual Tensor<1,dim>  gradient(const Point<dim> &p, const unsigned int component = 0) const;
};

template<int dim>
double ExactSolution<dim>::value(const Point<dim> &point, 
                const unsigned int component) const
{  	
	double x = point[0];
	double y = point[1]; 
	
	double r = sqrt(pow(x,2)+pow(y,2));
	return exp(-pow(r,2));
} 
 
template<int dim>
Tensor<1,dim> ExactSolution<dim>::gradient(const Point<dim> &point, 
				const unsigned int component) const
{
	double x = point[0];
	double y = point[1];
	double r = sqrt(pow(x,2)+pow(y,2));
	double u = exp(-pow(r,2));
	
	Tensor<1,dim> gradient;
	
	gradient[0] = -2*x*u;
	gradient[1] = -2*y*u;

	return gradient;

}


template<int dim>
class ExactSolutionForcingFunction : public Function<dim>
{
    public:
        ExactSolutionForcingFunction() : Function<dim>() {}       
        virtual double value(const Point<dim> &p, const unsigned int = 0) const;
};

template<int dim>
double ExactSolutionForcingFunction<dim>::value(const Point<dim> &point, 
                const unsigned int component) const
{
	double x = point[0];
	double y = point[1];
	double r = sqrt(pow(x,2)+pow(y,2));
	double u = exp(-pow(r,2));
	
	return -4*pow(x,2)*u - 4*pow(y,2)*u + 4*u;

}

template<int dim>
class BoundaryValues : public Function<dim>
{
    public:
        BoundaryValues() : Function<dim>() {}       
        virtual double value(const Point<dim> &p, const unsigned int = 0) const;
};

template<int dim>
double BoundaryValues<dim>::value(const Point<dim> &point, 
                const unsigned int component) const
{
	double r = 1;
	
	return exp(-pow(r,2));

}
 
#endif
