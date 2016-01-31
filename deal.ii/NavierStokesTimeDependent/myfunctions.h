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
	const double t = this->get_time();
	
	double x = p[0];
	double y = p[1];
	
	values(0) = exp(-t)*cos(M_PI*x);
	values(1) = exp(-t)*y*M_PI*sin(M_PI*x);
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
	const double t = this->get_time();
	
	const double nu = 1;
	const double x = p[0];
	const double y = p[1];
	
	values(0) = -exp(-t)*cos(M_PI*x) + nu*exp(-t)*pow(M_PI,2)*cos(M_PI*x) - exp(-2*t)*M_PI*cos(M_PI*x)*sin(M_PI*x) + 2*pow(y*t,2)*x;
	values(1) = -exp(-t)*y*M_PI*sin(M_PI*x) + nu*exp(-t)*pow(M_PI,3)*y*sin(M_PI*x) + exp(-2*t)*y*pow(M_PI,2) + 2*pow(x*t,2)*y;
	values(2) = 0;
}


template<int dim>
class ExactSolutionInitialConditions : public Function<dim>
{
    public:
        ExactSolutionInitialConditions() : Function<dim>(3) {}
        virtual void vector_value(const Point<dim> &p, 
                Vector<double> &values) const;
};

template<int dim>
void ExactSolutionInitialConditions<dim>::vector_value(const Point<dim> &p, 
                Vector<double> &values) const
{
	double x = p[0];
	double y = p[1];
	
	values(0) = cos(M_PI*x);
	values(1) = y*M_PI*sin(M_PI*x);
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
	const double t = this->get_time();
	
	double x = p[0];
	double y = p[1];
	
	values(0) = exp(-t)*cos(M_PI*x);
	values(1) = exp(-t)*y*M_PI*sin(M_PI*x);
	values(2) = pow(t*x*y,2); 
}

template<int dim>
void ExactSolution<dim>::vector_gradient(const Point<dim> &p, 
                std::vector<Tensor<1,dim> > & gradients) const
{
	const double t = this->get_time();
	
	double x = p[0];
	double y = p[1];
	
    gradients[0][0] = -exp(-t)*M_PI*sin(M_PI*x);
    gradients[0][1] = 0;
    
    gradients[1][0] = exp(-t)*pow(M_PI,2)*y*cos(M_PI*x);
    gradients[1][1] = exp(-t)*M_PI*sin(M_PI*x);
    
    gradients[2][0] = pow(t,2)*2*pow(y,2)*x;
    gradients[2][1] = pow(t,2)*2*pow(x,2)*y;
}
 

/******************************************************************************
 * Define inflow boundary conditions
 ******************************************************************************/
 
 template<int dim>
class InflowBoundaryValues : public Function<dim>
{
    public:
        InflowBoundaryValues() : Function<dim>(3) {}
        virtual void vector_value(const Point<dim> &p, 
                Vector<double> &values) const;
};

template<int dim>
void InflowBoundaryValues<dim>::vector_value(const Point<dim> &p, 
                Vector<double> &values) const
{  	
	double v_max = 1.5;
	double h = 0.41;
	double L = 2.2;
	
	if ((std::fabs(p[0]) < 1.0e-14) || (std::fabs(p[0] - L) < 1.0e-14))
	{
		values(0) = (4*v_max*p[1]*(h - p[1]))/pow(h,2);
		values(1) = 0;
	}
	else
	{
		values(0) = 0;
		values(1) = 0;
	}
	
	values(2) = 0;
 }
 
template<int dim>
class PeriodicBoundaryValues : public Function<dim>
{
    public:
        PeriodicBoundaryValues() : Function<dim>(3) {}
        virtual void vector_value(const Point<dim> &p, 
                Vector<double> &values) const;
};

template<int dim>
void PeriodicBoundaryValues<dim>::vector_value(const Point<dim> &p, 
                Vector<double> &values) const
{  	
	const double t = this->get_time();
	double h = 0.41;
	double L = 2.2;
	
	if ((std::fabs(p[0]) < 1.0e-14) || (std::fabs(p[0] - L) < 1.0e-14))
	{
		values(0) = (6*sin(M_PI*t/8)*p[1]*(h - p[1]))/pow(h,2);
		values(1) = 0;
	}
	else
	{
		values(0) = 0;
		values(1) = 0;
	}
	
	values(2) = 0;
 }
 
 
 /******************************************************************************
 * Define rotating circle forcing function
 ******************************************************************************/
 
 template<int dim>
class RotatingCircleForcingFunction : public Function<dim>
{
    public:
        RotatingCircleForcingFunction() : Function<dim>(3) {}
        virtual void vector_value(const Point<dim> &p, 
                Vector<double> &values) const;
};

template<int dim>
void RotatingCircleForcingFunction<dim>::vector_value(const Point<dim> &p, 
                Vector<double> &values) const
{ 	
	double x = p[0];
	double y = p[1];
	
	values(0) = -4*y*(1 - pow(x,2) - pow(y,2));
	values(1) = 4*x*(1 - pow(x,2) - pow(y,2));
	values(2) = 0;
 }
 

#endif
