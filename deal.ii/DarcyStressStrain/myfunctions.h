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
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/base/point.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/convergence_table.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/time_stepping.h>
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
#include <deal.II/fe/fe_dgq.h>
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
#include <time.h>
#include <omp.h>

using namespace dealii;

/*******************************************************************************
 * Define zero function
 ******************************************************************************/
template<int dim>
class MyZeroFunction : public Function<dim>
{
    public:
        MyZeroFunction() : Function<dim>(dim+12) {};       
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
	values(3) = 0;
	values(4) = 0;
	values(5) = 0;
	values(6) = 0;
	values(7) = 0;
	values(8) = 0;
	values(9) = 0;
	values(9) = 0;
	values(10) = 0;
	values(11) = 0;
	values(12) = 0;
}

#endif
