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
        virtual void vector_value(const Point<dim> &p, Vector<double> &values) const;
};

template<int dim>
void MyZeroFunction<dim>::vector_value(const Point<dim> &p, Vector<double> &values) const
{
	values(0) = 0;
	values(1) = 0;
	values(2) = 0;
}


/*******************************************************************************
 * Define boundary conditions, initial conditions, forcing function and 
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

	values(0) = cos(M_PI*p[0]);
	values(1) = p[1]*M_PI*sin(M_PI*p[0]);
	values(2) = 0; 
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
	values(0) = cos(M_PI*p[0])*pow(M_PI,2) - 1000*(p[0] - 0.5)*exp(-500*(pow(p[0]-0.5,2) + pow(p[1]-0.5,2)));
	values(1) = sin(M_PI*p[0])*p[1]*pow(M_PI,3) - 1000*(p[1] - 0.5)*exp(-500*(pow(p[0]-0.5,2) + pow(p[1]-0.5,2)));
	values(2) = 0;
}

template<int dim>
class ExactSolution : public Function<dim>
{
    public:
        ExactSolution() : Function<dim>(3) {}
        virtual void vector_value(const Point<dim> &p, Vector<double> &values) const;
        virtual void vector_gradient(const Point<dim> &p, 
							std::vector<Tensor<1,dim> > &gradients) const;
};

template<int dim>
void ExactSolution<dim>::vector_value(const Point<dim> &p, Vector<double> &values) const
{
    values(0) = cos(M_PI*p[0]);
	values(1) = p[1]*M_PI*sin(M_PI*p[0]);
	values(2) = exp(-500*(pow(p[0]-0.5,2) + pow(p[1]-0.5,2))); 
}

template<int dim>
void ExactSolution<dim>::vector_gradient(const Point<dim> &p, 
				std::vector<Tensor<1,dim> > & gradients) const
{
    gradients[0][0] = -M_PI*sin(M_PI*p[0]);
    gradients[0][1] = 0;
    
    gradients[1][0] = pow(M_PI,2)*p[1]*cos(M_PI*p[0]);
    gradients[1][1] = M_PI*sin(M_PI*p[0]);
    
    gradients[2][0] = -1000*(p[0] - 0.5)*exp(-500*(pow(p[0]-0.5,2) + pow(p[1]-0.5,2)));
    gradients[2][1] = -1000*(p[1] - 0.5)*exp(-500*(pow(p[0]-0.5,2) + pow(p[1]-0.5,2)));
}

/*******************************************************************************
 * Define boundary conditions and forcing functions for lid driven cavity
 ******************************************************************************/
template<int dim>
class DrivenCavityBoundaryValues : public Function<dim>
{
    public:
        DrivenCavityBoundaryValues() : Function<dim>(3) {}
        virtual void vector_value(const Point<dim> &p, 
                Vector<double> &values) const;
};

template<int dim>
void DrivenCavityBoundaryValues<dim>::vector_value(const Point<dim> &p, 
                Vector<double> &values) const
{  
	values(0) = (std::fabs(p[1]-1.0) < 1.0e-14 ? 1 : 0);
	values(1) = 0;
	values(2) = 0;     
}

/*******************************************************************************
 * Define boundary conditions, initial conditions, forcing function for a 
 * test case for adaptive meshing
 ******************************************************************************/
template<int dim>
class AdaptiveBoundaryValues : public Function<dim>
{
    public:
        AdaptiveBoundaryValues() : Function<dim>(3) {}
        virtual void vector_value(const Point<dim> &p, 
                Vector<double> &values) const;
};

/*******************************************************************************
 * Define exact solution for slit flow
 ******************************************************************************/
template<int dim>
class ExactSolutionSlitFlow : public Function<dim>
{
    public:
        ExactSolutionSlitFlow() : Function<dim>(3) {}
        virtual void vector_value(const Point<dim> &p, 
                Vector<double> &values) const;
        virtual void vector_gradient(const Point<dim> &p, 
				std::vector<Tensor<1,dim> > & gradients) const;
};

template<int dim>
void ExactSolutionSlitFlow<dim>::vector_value(const Point<dim> &p, 
                Vector<double> &values) const
{  	

	values(0) = -5*p[1]*(p[1] - 1)/20;
	values(1) = 0;
	values(2) = 10 - 5*p[0]/10; 
 }
 
 template<int dim>
void ExactSolutionSlitFlow<dim>::vector_gradient(const Point<dim> &p, 
				std::vector<Tensor<1,dim> > & gradients) const
{
    gradients[0][0] = 0;
    gradients[0][1] = 5*(1 - 2*p[1])/20;
    
    gradients[1][0] = 0;
    gradients[1][1] = 0;
    
    gradients[2][0] = -1.0/2;
    gradients[2][1] = 0;
} 

 
#endif
