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

template<int dim>
class StepBoundaryValues : public Function<dim>
{
    public:
        StepBoundaryValues() : Function<dim>(3) {}
        virtual void vector_value(const Point<dim> &p, 
                Vector<double> &values) const;
};

template<int dim>
void StepBoundaryValues<dim>::vector_value(const Point<dim> &p, 
                Vector<double> &values) const
{  
	values(0) = ((std::fabs(p[0]) < 1.0e-14) || (fabs(p[0]-5.0) < 1.0e-14 ))? p[1]*(p[1]-1) : 0;
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
	double x = p[0];
	double y = p[1];
	
	values(0) = cos(M_PI*x);
	values(1) = y*M_PI*sin(M_PI*x);
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
	const double x = p[0];
	const double y = p[1];
	
	const double RE = 1;
	
	values(0) = (1/RE)*pow(M_PI,2)*cos(M_PI*x) - M_PI*cos(M_PI*x)*sin(M_PI*x) + 2*pow(y,2)*x;
	values(1) = (1/RE)*pow(M_PI,3)*y*sin(M_PI*x) + y*pow(M_PI,2) + 2*pow(x,2)*y;
	values(2) = 0;
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
	double x = p[0];
	double y = p[1];
	
	values(0) = cos(M_PI*x);
	values(1) = y*M_PI*sin(M_PI*x);
	values(2) = pow(x*y,2); 
}

template<int dim>
void ExactSolution<dim>::vector_gradient(const Point<dim> &p, 
                std::vector<Tensor<1,dim> > & gradients) const
{
	double x = p[0];
	double y = p[1];
	
    gradients[0][0] = -M_PI*sin(M_PI*x);
    gradients[0][1] = 0;
    
    gradients[1][0] = pow(M_PI,2)*y*cos(M_PI*x);
    gradients[1][1] = M_PI*sin(M_PI*x);
    
    gradients[2][0] = 2*pow(y,2)*x;
    gradients[2][1] = 2*pow(x,2)*y;
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

template<int dim>
void AdaptiveBoundaryValues<dim>::vector_value(const Point<dim> &p, 
                Vector<double> &values) const
{  	
	const double t = this->get_time();
	
	if (std::fabs(p[1]) < 1.0e-14) //bottom boundary
	{
		values(0) = sin(2*M_PI*t + M_PI/2);
		values(1) = 0;
	}
	else if (std::fabs(p[0] - 1) < 1.0e-14) //right boundary
	{	
		values(0) = 0;
		values(1) = 0;
	}
	else if (std::fabs(p[1] - 1) < 1.0e-14) //top boundary
	{
		values(0) = sin(2*M_PI*t);
		values(1) = 0;
	}
	else //left boundary
	{
		values(0) = 0;
		values(1) = 0;
	}
	
	values(2) = 0; 
 }
 
template<int dim>
class AdaptiveInitialCondition : public Function<dim>
{
    public:
        AdaptiveInitialCondition() : Function<dim>(3) {};       
        virtual void vector_value(const Point<dim> &p, 
				Vector<double> &values) const;
};

template<int dim>
void AdaptiveInitialCondition<dim>::vector_value(const Point<dim> &p, 
                Vector<double> &values) const
{
	values(0) = 0;
	values(1) = 0;
	values(2) = 0;
}

template<int dim>
class AdaptiveForcingFunction : public Function<dim>
{
    public:
        AdaptiveForcingFunction() : Function<dim>(3) {};       
        virtual void vector_value(const Point<dim> &p, 
				Vector<double> &values) const;
};

template<int dim>
void AdaptiveForcingFunction<dim>::vector_value(const Point<dim> &p, 
                Vector<double> &values) const
{
	const double t = this->get_time();
	
	values(0) = 0;
	values(1) = 0;
	values(2) = 0;
}

/*******************************************************************************
 * Define boundary conditions, initial conditions, forcing function for a 
 * simple test case for adaptive meshing
 ******************************************************************************/
template<int dim>
class SimpleAdaptiveBoundaryValues : public Function<dim>
{
    public:
        SimpleAdaptiveBoundaryValues() : Function<dim>(3) {}
        virtual void vector_value(const Point<dim> &p, 
                Vector<double> &values) const;
};

template<int dim>
void SimpleAdaptiveBoundaryValues<dim>::vector_value(const Point<dim> &p, 
                Vector<double> &values) const
{  	
	const double t = this->get_time();
	
	if (std::fabs(p[1]) < 1.0e-14) //bottom boundary
	{
		values(0) = 0;
		values(1) = 0;
	}
	else if (std::fabs(p[0] - 1) < 1.0e-14) //right boundary
	{	
		values(0) = p[1];
		values(1) = 0;
	}
	else if (std::fabs(p[1] - 1) < 1.0e-14) //top boundary
	{
		values(0) = 1;
		values(1) = 0;
	}
	else //left boundary
	{
		values(0) = p[1];
		values(1) = 0;
	}
	
	values(2) = 0; 
 }
 
template<int dim>
class SimpleAdaptiveInitialCondition : public Function<dim>
{
    public:
        SimpleAdaptiveInitialCondition() : Function<dim>(3) {};       
        virtual void vector_value(const Point<dim> &p, 
				Vector<double> &values) const;
};

template<int dim>
void SimpleAdaptiveInitialCondition<dim>::vector_value(const Point<dim> &p, 
                Vector<double> &values) const
{
	values(0) = p[1];
	values(1) = 0;
	values(2) = 0;
}

template<int dim>
class SimpleAdaptiveForcingFunction : public Function<dim>
{
    public:
        SimpleAdaptiveForcingFunction() : Function<dim>(3) {};       
        virtual void vector_value(const Point<dim> &p, 
				Vector<double> &values) const;
};

template<int dim>
void SimpleAdaptiveForcingFunction<dim>::vector_value(const Point<dim> &p, 
                Vector<double> &values) const
{
	const double t = this->get_time();
	
	values(0) = 0;
	values(1) = 0;
	values(2) = 0;
}

#endif
