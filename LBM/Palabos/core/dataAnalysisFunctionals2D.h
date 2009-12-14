/* This file is part of the Palabos library.
 * Copyright (C) 2009 Jonas Latt
 * E-mail contact: jonas@lbmethod.org
 * The most recent release of Palabos can be downloaded at 
 * <http://www.lbmethod.org/palabos/>
 *
 * The library Palabos is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * The library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/** \file
 * Helper functions for domain initialization -- header file.
 */
#ifndef DATA_ANALYSIS_FUNCTIONALS_2D_H
#define DATA_ANALYSIS_FUNCTIONALS_2D_H

#include "core/globalDefs.h"
#include "core/blockLatticeBase2D.h"
#include "core/dataFieldBase2D.h"
#include "atomicBlock/dataProcessorWrapper2D.h"
#include "atomicBlock/reductiveDataProcessorWrapper2D.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"

namespace plb {


/* *************** PART I ******************************************** */
/* *************** Analysis of the block-lattice ********************* */
/* ******************************************************************* */

/* *************** Reductive Data Functionals for BlockLattice ******* */

template<typename T, template<typename U> class Descriptor> 
class BoxSumRhoBarFunctional2D : public ReductiveBoxProcessingFunctional2D_L<T,Descriptor>
{
public:
    BoxSumRhoBarFunctional2D();
    virtual void process(Box2D domain, BlockLattice2D<T,Descriptor>& lattice);
    virtual BoxSumRhoBarFunctional2D<T,Descriptor>* clone() const;
    T getSumRhoBar() const;
private:
    plint sumRhoBarId;
};

template<typename T, template<typename U> class Descriptor> 
class BoxSumEnergyFunctional2D : public ReductiveBoxProcessingFunctional2D_L<T,Descriptor>
{
public:
    BoxSumEnergyFunctional2D();
    virtual void process(Box2D domain, BlockLattice2D<T,Descriptor>& lattice);
    virtual BoxSumEnergyFunctional2D<T,Descriptor>* clone() const;
    T getSumEnergy() const;
private:
    plint sumEnergyId;
};

template<typename T, template<typename U> class Descriptor, class BoolMask> 
class CountLatticeElementsFunctional2D : public ReductiveBoxProcessingFunctional2D_L<T,Descriptor>
{
public:
    CountLatticeElementsFunctional2D(BoolMask boolMask_);
    virtual void process(Box2D domain, BlockLattice2D<T,Descriptor>& lattice);
    virtual CountLatticeElementsFunctional2D<T,Descriptor,BoolMask>* clone() const;
    plint getCount() const;
private:
    plint countId;
    BoolMask boolMask;
};


/* *************** Data Functionals for BlockLattice ***************** */

template<typename T, template<typename U> class Descriptor> 
class ExtractLatticeSubDomainFunctional2D : public BoxProcessingFunctional2D_LL<T,Descriptor,Descriptor>
{
public:
    virtual void process(Box2D domain, BlockLattice2D<T,Descriptor>& lattice1,
                                       BlockLattice2D<T,Descriptor>& lattice2);
    virtual ExtractLatticeSubDomainFunctional2D<T,Descriptor>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, template<typename U> class Descriptor> 
class BoxDensityFunctional2D : public BoxProcessingFunctional2D_LS<T,Descriptor>
{
public:
    virtual void process(Box2D domain, BlockLattice2D<T,Descriptor>& lattice, ScalarField2D<T>& tensorField);
    virtual BoxDensityFunctional2D<T,Descriptor>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, template<typename U> class Descriptor> 
class BoxRhoBarFunctional2D : public BoxProcessingFunctional2D_LS<T,Descriptor>
{
public:
    virtual void process(Box2D domain, BlockLattice2D<T,Descriptor>& lattice, ScalarField2D<T>& tensorField);
    virtual BoxRhoBarFunctional2D<T,Descriptor>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, template<typename U> class Descriptor> 
class BoxKineticEnergyFunctional2D : public BoxProcessingFunctional2D_LS<T,Descriptor>
{
public:
    virtual void process(Box2D domain, BlockLattice2D<T,Descriptor>& lattice, ScalarField2D<T>& tensorField);
    virtual BoxKineticEnergyFunctional2D<T,Descriptor>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, template<typename U> class Descriptor> 
class BoxVelocityNormFunctional2D : public BoxProcessingFunctional2D_LS<T,Descriptor>
{
public:
    virtual void process(Box2D domain, BlockLattice2D<T,Descriptor>& lattice, ScalarField2D<T>& tensorField);
    virtual BoxVelocityNormFunctional2D<T,Descriptor>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, template<typename U> class Descriptor> 
class BoxVelocityComponentFunctional2D : public BoxProcessingFunctional2D_LS<T,Descriptor>
{
public:
    BoxVelocityComponentFunctional2D(int iComponent_);
    virtual void process(Box2D domain, BlockLattice2D<T,Descriptor>& lattice, ScalarField2D<T>& scalarField);
    virtual BoxVelocityComponentFunctional2D<T,Descriptor>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    int iComponent;
};

template<typename T, template<typename U> class Descriptor> 
class BoxVelocityFunctional2D : public BoxProcessingFunctional2D_LT<T,Descriptor, Descriptor<T>::d>
{
public:
    virtual void process(Box2D domain, BlockLattice2D<T,Descriptor>& lattice,
                                       TensorField2D<T, Descriptor<T>::d>& tensorField);
    virtual BoxVelocityFunctional2D<T,Descriptor>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, template<typename U> class Descriptor> 
class BoxDeviatoricStressFunctional2D :
    public BoxProcessingFunctional2D_LT<T,Descriptor, SymmetricTensor<T,Descriptor>::n>
{
public:
    virtual void process(Box2D domain, BlockLattice2D<T,Descriptor>& lattice,
                                       TensorField2D<T, SymmetricTensor<T,Descriptor>::n>& PiNeq);
    virtual BoxDeviatoricStressFunctional2D<T,Descriptor>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, template<typename U> class Descriptor> 
class BoxStrainRateFromStressFunctional2D :
    public BoxProcessingFunctional2D_LT<T,Descriptor, SymmetricTensor<T,Descriptor>::n>
{
public:
    virtual void process(Box2D domain, BlockLattice2D<T,Descriptor>& lattice,
                                       TensorField2D<T, SymmetricTensor<T,Descriptor>::n>& S);
    virtual BoxStrainRateFromStressFunctional2D<T,Descriptor>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, template<typename U> class Descriptor> 
class BoxPopulationFunctional2D : public BoxProcessingFunctional2D_LS<T,Descriptor>
{
public:
    BoxPopulationFunctional2D(plint iComponent_);
    virtual void process(Box2D domain, BlockLattice2D<T,Descriptor>& lattice,
                                       ScalarField2D<T>& population);
    virtual BoxPopulationFunctional2D<T,Descriptor>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    plint iComponent;
};


/* *************** PART II ******************************************* */
/* *************** Analysis of the scalar-field ********************** */
/* ******************************************************************* */

/* *************** Reductive Data Functionals for ScalarField ******** */

template<typename T>
class BoxScalarSumFunctional2D : public ReductiveBoxProcessingFunctional2D_S<T>
{
public:
    BoxScalarSumFunctional2D();
    virtual void process(Box2D domain, ScalarField2D<T>& scalarField);
    virtual BoxScalarSumFunctional2D<T>* clone() const;
    T getSumScalar() const;
private:
    plint sumScalarId;
};

template<typename T>
class BoxScalarMinFunctional2D : public ReductiveBoxProcessingFunctional2D_S<T>
{
public:
    BoxScalarMinFunctional2D();
    virtual void process(Box2D domain, ScalarField2D<T>& scalarField);
    virtual BoxScalarMinFunctional2D<T>* clone() const;
    T getMinScalar() const;
private:
    plint maxScalarId;
};

template<typename T>
class BoxScalarMaxFunctional2D : public ReductiveBoxProcessingFunctional2D_S<T>
{
public:
    BoxScalarMaxFunctional2D();
    virtual void process(Box2D domain, ScalarField2D<T>& scalarField);
    virtual BoxScalarMaxFunctional2D<T>* clone() const;
    T getMaxScalar() const;
private:
    plint maxScalarId;
};

template<typename T>
class BoundedBoxScalarSumFunctional2D : public BoundedReductiveBoxProcessingFunctional2D_S<T>
{
public:
    BoundedBoxScalarSumFunctional2D();
    virtual void processBulk(Box2D domain, ScalarField2D<T>& scalarField);
    virtual void processEdge( int direction, int orientation, Box2D domain,
                              ScalarField2D<T>& scalarField );
    virtual void processCorner( int normalX, int normalY, Box2D domain,
                                ScalarField2D<T>& scalarField );
    virtual BoundedBoxScalarSumFunctional2D<T>* clone() const;
    T getSumScalar() const;
private:
    plint sumScalarId;
};

template<typename T, class BoolMask> 
class CountScalarElementsFunctional2D : public ReductiveBoxProcessingFunctional2D_S<T>
{
public:
    CountScalarElementsFunctional2D(BoolMask boolMask_);
    virtual void process(Box2D domain, ScalarField2D<T>& field);
    virtual CountScalarElementsFunctional2D<T,BoolMask>* clone() const;
    plint getCount() const;
private:
    plint countId;
    BoolMask boolMask;
};

/* *************** Data Functionals for scalar-fields **************** */

template<typename T>
class ExtractScalarSubDomainFunctional2D : public BoxProcessingFunctional2D_SS<T>
{
public:
    virtual void process(Box2D domain, ScalarField2D<T>& field1, ScalarField2D<T>& field2);
    virtual ExtractScalarSubDomainFunctional2D<T>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, class Function>
class ApplyScalarFunctional2D : public BoxProcessingFunctional2D_S<T> {
public:
    ApplyScalarFunctional2D(Function f_);
    virtual void process(Box2D domain, ScalarField2D<T>& field);
    virtual ApplyScalarFunctional2D<T,Function>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    Function f;
};

template<typename T, class EvalFunction>
class EvaluateScalarFunctional2D : public BoxProcessingFunctional2D_SS<T> {
public:
    EvaluateScalarFunctional2D(EvalFunction f_);
    virtual void process(Box2D domain, ScalarField2D<T>& field,
                                       ScalarField2D<T>& result);
    virtual EvaluateScalarFunctional2D<T,EvalFunction>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    EvalFunction f;
};

template<typename T>
class A_plus_alpha_functional2D : public BoxProcessingFunctional2D_SS<T>
{
public:
    A_plus_alpha_functional2D(T alpha_);
    virtual void process(Box2D domain, ScalarField2D<T>& A,
                                       ScalarField2D<T>& result);
    virtual A_plus_alpha_functional2D<T>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    T alpha;
};

template<typename T>
class A_minus_alpha_functional2D : public BoxProcessingFunctional2D_SS<T>
{
public:
    A_minus_alpha_functional2D(T alpha_);
    virtual void process(Box2D domain, ScalarField2D<T>& A,
                                       ScalarField2D<T>& result);
    virtual A_minus_alpha_functional2D<T>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    T alpha;
};

template<typename T>
class Alpha_minus_A_functional2D : public BoxProcessingFunctional2D_SS<T>
{
public:
    Alpha_minus_A_functional2D(T alpha_);
    virtual void process(Box2D domain, ScalarField2D<T>& A,
                                       ScalarField2D<T>& result);
    virtual Alpha_minus_A_functional2D<T>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    T alpha;
};

template<typename T>
class A_times_alpha_functional2D : public BoxProcessingFunctional2D_SS<T>
{
public:
    A_times_alpha_functional2D(T alpha_);
    virtual void process(Box2D domain, ScalarField2D<T>& A,
                                       ScalarField2D<T>& result);
    virtual A_times_alpha_functional2D<T>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    T alpha;
};

template<typename T>
class A_dividedBy_alpha_functional2D : public BoxProcessingFunctional2D_SS<T>
{
public:
    A_dividedBy_alpha_functional2D(T alpha_);
    virtual void process(Box2D domain, ScalarField2D<T>& A,
                                       ScalarField2D<T>& result);
    virtual A_dividedBy_alpha_functional2D<T>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    T alpha;
};

template<typename T>
class Alpha_dividedBy_A_functional2D : public BoxProcessingFunctional2D_SS<T>
{
public:
    Alpha_dividedBy_A_functional2D(T alpha_);
    virtual void process(Box2D domain, ScalarField2D<T>& A,
                                       ScalarField2D<T>& result);
    virtual Alpha_dividedBy_A_functional2D<T>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    T alpha;
};

template<typename T>
class A_plus_alpha_inplace_functional2D : public BoxProcessingFunctional2D_S<T>
{
public:
    A_plus_alpha_inplace_functional2D(T alpha_);
    virtual void process(Box2D domain, ScalarField2D<T>& A);
    virtual A_plus_alpha_inplace_functional2D<T>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    T alpha;
};

template<typename T>
class A_minus_alpha_inplace_functional2D : public BoxProcessingFunctional2D_S<T>
{
public:
    A_minus_alpha_inplace_functional2D(T alpha_);
    virtual void process(Box2D domain, ScalarField2D<T>& A);
    virtual A_minus_alpha_inplace_functional2D<T>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    T alpha;
};

template<typename T>
class A_times_alpha_inplace_functional2D : public BoxProcessingFunctional2D_S<T>
{
public:
    A_times_alpha_inplace_functional2D(T alpha_);
    virtual void process(Box2D domain, ScalarField2D<T>& A);
    virtual A_times_alpha_inplace_functional2D<T>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    T alpha;
};

template<typename T>
class A_dividedBy_alpha_inplace_functional2D : public BoxProcessingFunctional2D_S<T>
{
public:
    A_dividedBy_alpha_inplace_functional2D(T alpha_);
    virtual void process(Box2D domain, ScalarField2D<T>& A);
    virtual A_dividedBy_alpha_inplace_functional2D<T>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    T alpha;
};


template<typename T>
class A_plus_B_functional2D : public ScalarFieldBoxProcessingFunctional2D<T>
{
public:
    virtual void process(Box2D domain, std::vector<ScalarField2D<T>*> scalarFields);
    virtual A_plus_B_functional2D<T>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class A_minus_B_functional2D : public ScalarFieldBoxProcessingFunctional2D<T>
{
public:
    virtual void process(Box2D domain, std::vector<ScalarField2D<T>*> scalarFields);
    virtual A_minus_B_functional2D<T>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class A_times_B_functional2D : public ScalarFieldBoxProcessingFunctional2D<T>
{
public:
    virtual void process(Box2D domain, std::vector<ScalarField2D<T>*> scalarFields);
    virtual A_times_B_functional2D<T>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class A_dividedBy_B_functional2D : public ScalarFieldBoxProcessingFunctional2D<T>
{
public:
    virtual void process(Box2D domain, std::vector<ScalarField2D<T>*> scalarFields);
    virtual A_dividedBy_B_functional2D<T>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class A_plus_B_inplace_functional2D : public BoxProcessingFunctional2D_SS<T>
{
public:
    virtual void process(Box2D domain, ScalarField2D<T>& A, ScalarField2D<T>& B);
    virtual A_plus_B_inplace_functional2D<T>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class A_minus_B_inplace_functional2D : public BoxProcessingFunctional2D_SS<T>
{
public:
    virtual void process(Box2D domain, ScalarField2D<T>& A, ScalarField2D<T>& B);
    virtual A_minus_B_inplace_functional2D<T>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class A_times_B_inplace_functional2D : public BoxProcessingFunctional2D_SS<T>
{
public:
    virtual void process(Box2D domain, ScalarField2D<T>& A, ScalarField2D<T>& B);
    virtual A_times_B_inplace_functional2D<T>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class A_dividedBy_B_inplace_functional2D : public BoxProcessingFunctional2D_SS<T>
{
public:
    virtual void process(Box2D domain, ScalarField2D<T>& A, ScalarField2D<T>& B);
    virtual A_dividedBy_B_inplace_functional2D<T>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

/* *************** PART III ****************************************** */
/* *************** Analysis of the tensor-field ********************** */
/* ******************************************************************* */

template<typename T, int nDim, class BoolMask> 
class CountTensorElementsFunctional2D : public ReductiveBoxProcessingFunctional2D_T<T,nDim>
{
public:
    CountTensorElementsFunctional2D(BoolMask boolMask_);
    virtual void process(Box2D domain, TensorField2D<T,nDim>& field);
    virtual CountTensorElementsFunctional2D<T,nDim,BoolMask>* clone() const;
    plint getCount() const;
private:
    plint countId;
    BoolMask boolMask;
};

template<typename T, int nDim>
class ExtractTensorSubDomainFunctional2D : public BoxProcessingFunctional2D_TT<T,nDim,nDim>
{
public:
    virtual void process(Box2D domain, TensorField2D<T,nDim>& field1, TensorField2D<T,nDim>& field2);
    virtual ExtractTensorSubDomainFunctional2D<T,nDim>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, int nDim>
class ExtractTensorComponentFunctional2D : public BoxProcessingFunctional2D_ST<T,nDim>
{
public:
    ExtractTensorComponentFunctional2D(int iComponent_);
    virtual void process(Box2D domain, ScalarField2D<T>& scalarField,
                                       TensorField2D<T,nDim>& tensorField);
    virtual ExtractTensorComponentFunctional2D<T,nDim>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    int iComponent;
};

template<typename T, int nDim>
class ComputeNormFunctional2D : public BoxProcessingFunctional2D_ST<T,nDim>
{
public:
    virtual void process(Box2D domain, ScalarField2D<T>& scalarField,
                                       TensorField2D<T,nDim>& tensorField);
    virtual ComputeNormFunctional2D<T,nDim>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, int nDim>
class ComputeNormSqrFunctional2D : public BoxProcessingFunctional2D_ST<T,nDim>
{
public:
    virtual void process(Box2D domain, ScalarField2D<T>& scalarField,
                                       TensorField2D<T,nDim>& tensorField);
    virtual ComputeNormSqrFunctional2D<T,nDim>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class ComputeSymmetricTensorNormFunctional2D :
    public BoxProcessingFunctional2D_ST<T,3>
{
public:
    virtual void process(Box2D domain, ScalarField2D<T>& scalarField,
                                       TensorField2D<T,3>& tensorField);
    virtual ComputeSymmetricTensorNormFunctional2D<T>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class ComputeSymmetricTensorNormSqrFunctional2D :
    public BoxProcessingFunctional2D_ST<T,3>
{
public:
    virtual void process(Box2D domain, ScalarField2D<T>& scalarField,
                                       TensorField2D<T,3>& tensorField);
    virtual ComputeSymmetricTensorNormSqrFunctional2D<T>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class ComputeSymmetricTensorTraceFunctional2D :
    public BoxProcessingFunctional2D_ST<T,3>
{
public:
    virtual void process(Box2D domain, ScalarField2D<T>& scalarField,
                                       TensorField2D<T,3>& tensorField);
    virtual ComputeSymmetricTensorTraceFunctional2D<T>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, int nDim>
class BoxBulkVorticityFunctional2D : public BoxProcessingFunctional2D_ST<T,nDim>
{
public:
    virtual void process(Box2D domain, ScalarField2D<T>& vorticity,
                                       TensorField2D<T,nDim>& velocity);
    virtual BoxBulkVorticityFunctional2D<T,nDim>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, int nDim>
class BoxVorticityFunctional2D : public BoundedBoxProcessingFunctional2D_ST<T,nDim>
{
public:
    virtual void processBulk( Box2D domain, ScalarField2D<T>& vorticity,
                                            TensorField2D<T,nDim>& velocity );
    virtual void processEdge( int direction, int orientation, Box2D domain,
                              ScalarField2D<T>& vorticity,
                              TensorField2D<T,nDim>& velocity );
    virtual void processCorner( int normalX, int normalY, Box2D domain,
                                ScalarField2D<T>& vorticity,
                                TensorField2D<T,nDim>& velocity );
    virtual BoxVorticityFunctional2D<T,nDim>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, int nDim>
class BoxBulkStrainRateFunctional2D :
    public BoxProcessingFunctional2D_TT<T,nDim,SymmetricTensorImpl<T,nDim>::n>
{
public:
    virtual void process(Box2D domain, TensorField2D<T,nDim>& velocity,
                                       TensorField2D<T,SymmetricTensorImpl<T,nDim>::n>& S);
    virtual BoxBulkStrainRateFunctional2D<T,nDim>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, int nDim>
class BoxStrainRateFunctional2D :
    public BoundedBoxProcessingFunctional2D_TT<T,nDim,SymmetricTensorImpl<T,nDim>::n>
{
public:
    virtual void processBulk( Box2D domain, TensorField2D<T,nDim>& velocity,
                                            TensorField2D<T,SymmetricTensorImpl<T,nDim>::n>& S );
    virtual void processEdge( int direction, int orientation, Box2D domain,
                              TensorField2D<T,nDim>& velocity,
                              TensorField2D<T,SymmetricTensorImpl<T,nDim>::n>& S );
    virtual void processCorner( int normalX, int normalY, Box2D domain,
                                TensorField2D<T,nDim>& velocity,
                                TensorField2D<T,SymmetricTensorImpl<T,nDim>::n>& S );
    virtual BoxStrainRateFunctional2D<T,nDim>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
};


template<typename T, int nDim>
class Tensor_A_plus_B_functional2D : public TensorFieldBoxProcessingFunctional2D<T,nDim>
{
public:
    virtual void process(Box2D domain, std::vector<TensorField2D<T,nDim>*> tensorFields);
    virtual Tensor_A_plus_B_functional2D<T,nDim>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, int nDim>
class Tensor_A_minus_B_functional2D : public TensorFieldBoxProcessingFunctional2D<T,nDim>
{
public:
    virtual void process(Box2D domain, std::vector<TensorField2D<T,nDim>*> tensorFields);
    virtual Tensor_A_minus_B_functional2D<T,nDim>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, int nDim>
class Tensor_A_times_B_functional2D : public TensorFieldBoxProcessingFunctional2D<T,nDim>
{
public:
    virtual void process(Box2D domain, std::vector<TensorField2D<T,nDim>*> tensorFields);
    virtual Tensor_A_times_B_functional2D<T,nDim>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, int nDim>
class Tensor_A_dividedBy_B_functional2D : public TensorFieldBoxProcessingFunctional2D<T,nDim>
{
public:
    virtual void process(Box2D domain, std::vector<TensorField2D<T,nDim>*> tensorFields);
    virtual Tensor_A_dividedBy_B_functional2D<T,nDim>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, int nDim>
class Tensor_A_plus_B_inplace_functional2D : public BoxProcessingFunctional2D_TT<T,nDim,nDim>
{
public:
    virtual void process(Box2D domain, TensorField2D<T,nDim>& A, TensorField2D<T,nDim>& B);
    virtual Tensor_A_plus_B_inplace_functional2D<T,nDim>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, int nDim>
class Tensor_A_minus_B_inplace_functional2D : public BoxProcessingFunctional2D_TT<T,nDim,nDim>
{
public:
    virtual void process(Box2D domain, TensorField2D<T,nDim>& A, TensorField2D<T,nDim>& B);
    virtual Tensor_A_minus_B_inplace_functional2D<T,nDim>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, int nDim>
class Tensor_A_times_B_inplace_functional2D : public BoxProcessingFunctional2D_TT<T,nDim,nDim>
{
public:
    virtual void process(Box2D domain, TensorField2D<T,nDim>& A, TensorField2D<T,nDim>& B);
    virtual Tensor_A_times_B_inplace_functional2D<T,nDim>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, int nDim>
class Tensor_A_dividedBy_B_inplace_functional2D : public BoxProcessingFunctional2D_TT<T,nDim,nDim>
{
public:
    virtual void process(Box2D domain, TensorField2D<T,nDim>& A, TensorField2D<T,nDim>& B);
    virtual Tensor_A_dividedBy_B_inplace_functional2D<T,nDim>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

}  // namespace plb

#endif  // DATA_ANALYSIS_FUNCTIONALS_2D_H

// Explicitly include generic algorithms which are never precompiled (not even in precompiled version)
#include "core/dataAnalysisGenerics2D.h"
