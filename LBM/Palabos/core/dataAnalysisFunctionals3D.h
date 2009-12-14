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
#ifndef DATA_ANALYSIS_FUNCTIONALS_3D_H
#define DATA_ANALYSIS_FUNCTIONALS_3D_H

#include "core/globalDefs.h"
#include "core/blockLatticeBase3D.h"
#include "core/dataFieldBase3D.h"
#include "atomicBlock/dataProcessorWrapper3D.h"
#include "atomicBlock/reductiveDataProcessorWrapper3D.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"

namespace plb {

/* *************** PART I ******************************************** */
/* *************** Analysis of the block-lattice ********************* */
/* ******************************************************************* */

/* *************** Reductive Data Functionals for BlockLattice ******* */

template<typename T, template<typename U> class Descriptor> 
class BoxSumRhoBarFunctional3D : public ReductiveBoxProcessingFunctional3D_L<T,Descriptor>
{
public:
    BoxSumRhoBarFunctional3D();
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice);
    virtual BoxSumRhoBarFunctional3D<T,Descriptor>* clone() const;
    T getSumRhoBar() const;
private:
    plint sumRhoBarId;
};

template<typename T, template<typename U> class Descriptor> 
class BoxSumEnergyFunctional3D : public ReductiveBoxProcessingFunctional3D_L<T,Descriptor>
{
public:
    BoxSumEnergyFunctional3D();
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice);
    virtual BoxSumEnergyFunctional3D<T,Descriptor>* clone() const;
    T getSumEnergy() const;
private:
    plint sumEnergyId;
};

template<typename T, template<typename U> class Descriptor, class BoolMask> 
class CountLatticeElementsFunctional3D : public ReductiveBoxProcessingFunctional3D_L<T,Descriptor>
{
public:
    CountLatticeElementsFunctional3D(BoolMask boolMask_);
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice);
    virtual CountLatticeElementsFunctional3D<T,Descriptor,BoolMask>* clone() const;
    plint getCount() const;
private:
    plint countId;
    BoolMask boolMask;
};


/* *************** Data Functionals for BlockLattice ***************** */

template<typename T, template<typename U> class Descriptor> 
class ExtractLatticeSubDomainFunctional3D : public BoxProcessingFunctional3D_LL<T,Descriptor,Descriptor>
{
public:
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice1,
                                       BlockLattice3D<T,Descriptor>& lattice2);
    virtual ExtractLatticeSubDomainFunctional3D<T,Descriptor>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, template<typename U> class Descriptor> 
class BoxDensityFunctional3D : public BoxProcessingFunctional3D_LS<T,Descriptor>
{
public:
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice, ScalarField3D<T>& tensorField);
    virtual BoxDensityFunctional3D<T,Descriptor>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, template<typename U> class Descriptor> 
class BoxRhoBarFunctional3D : public BoxProcessingFunctional3D_LS<T,Descriptor>
{
public:
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice, ScalarField3D<T>& tensorField);
    virtual BoxRhoBarFunctional3D<T,Descriptor>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, template<typename U> class Descriptor> 
class BoxKineticEnergyFunctional3D : public BoxProcessingFunctional3D_LS<T,Descriptor>
{
public:
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice, ScalarField3D<T>& tensorField);
    virtual BoxKineticEnergyFunctional3D<T,Descriptor>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, template<typename U> class Descriptor> 
class BoxVelocityNormFunctional3D : public BoxProcessingFunctional3D_LS<T,Descriptor>
{
public:
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice, ScalarField3D<T>& tensorField);
    virtual BoxVelocityNormFunctional3D<T,Descriptor>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, template<typename U> class Descriptor> 
class BoxVelocityComponentFunctional3D : public BoxProcessingFunctional3D_LS<T,Descriptor>
{
public:
    BoxVelocityComponentFunctional3D(int iComponent);
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice, ScalarField3D<T>& scalarField);
    virtual BoxVelocityComponentFunctional3D<T,Descriptor>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    int iComponent;
};

template<typename T, template<typename U> class Descriptor> 
class BoxVelocityFunctional3D : public BoxProcessingFunctional3D_LT<T,Descriptor, Descriptor<T>::d>
{
public:
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
                                       TensorField3D<T, Descriptor<T>::d>& tensorField);
    virtual BoxVelocityFunctional3D<T,Descriptor>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, template<typename U> class Descriptor> 
class BoxDeviatoricStressFunctional3D :
    public BoxProcessingFunctional3D_LT<T,Descriptor, SymmetricTensor<T,Descriptor>::n>
{
public:
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
                                       TensorField3D<T, SymmetricTensor<T,Descriptor>::n>& PiNeq);
    virtual BoxDeviatoricStressFunctional3D<T,Descriptor>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, template<typename U> class Descriptor> 
class BoxStrainRateFromStressFunctional3D :
    public BoxProcessingFunctional3D_LT<T,Descriptor, SymmetricTensor<T,Descriptor>::n>
{
public:
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
                                       TensorField3D<T, SymmetricTensor<T,Descriptor>::n>& S);
    virtual BoxStrainRateFromStressFunctional3D<T,Descriptor>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, template<typename U> class Descriptor> 
class BoxPopulationFunctional3D : public BoxProcessingFunctional3D_LS<T,Descriptor>
{
public:
    BoxPopulationFunctional3D(plint iComponent_);
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
                                       ScalarField3D<T>& population);
    virtual BoxPopulationFunctional3D<T,Descriptor>* clone() const;
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
class BoxScalarSumFunctional3D : public ReductiveBoxProcessingFunctional3D_S<T>
{
public:
    BoxScalarSumFunctional3D();
    virtual void process(Box3D domain, ScalarField3D<T>& scalarField);
    virtual BoxScalarSumFunctional3D<T>* clone() const;
    T getSumScalar() const;
private:
    plint sumScalarId;
};

template<typename T>
class BoxScalarMinFunctional3D : public ReductiveBoxProcessingFunctional3D_S<T>
{
public:
    BoxScalarMinFunctional3D();
    virtual void process(Box3D domain, ScalarField3D<T>& scalarField);
    virtual BoxScalarMinFunctional3D<T>* clone() const;
    T getMinScalar() const;
private:
    plint maxScalarId;
};

template<typename T>
class BoxScalarMaxFunctional3D : public ReductiveBoxProcessingFunctional3D_S<T>
{
public:
    BoxScalarMaxFunctional3D();
    virtual void process(Box3D domain, ScalarField3D<T>& scalarField);
    virtual BoxScalarMaxFunctional3D<T>* clone() const;
    T getMaxScalar() const;
private:
    plint maxScalarId;
};

template<typename T>
class BoundedBoxScalarSumFunctional3D : public BoundedReductiveBoxProcessingFunctional3D_S<T>
{
public:
    BoundedBoxScalarSumFunctional3D();
    virtual void processBulk(Box3D domain, ScalarField3D<T>& scalarField);
    virtual void processPlane( int direction, int orientation, Box3D domain,
                               ScalarField3D<T>& scalarField );
    virtual void processEdge( int plane, int normal1, int normal2,
                              Box3D domain, ScalarField3D<T>& scalarField );
    virtual void processCorner( int normalX, int normalY, int normalZ,
                                Box3D domain, ScalarField3D<T>& scalarField );
    virtual BoundedBoxScalarSumFunctional3D<T>* clone() const;
    T getSumScalar() const;
private:
    plint sumScalarId;
};


template<typename T, class BoolMask> 
class CountScalarElementsFunctional3D : public ReductiveBoxProcessingFunctional3D_S<T>
{
public:
    CountScalarElementsFunctional3D(BoolMask boolMask_);
    virtual void process(Box3D domain, ScalarField3D<T>& field);
    virtual CountScalarElementsFunctional3D<T,BoolMask>* clone() const;
    plint getCount() const;
private:
    plint countId;
    BoolMask boolMask;
};


/* *************** Data Functionals for scalar-fields **************** */

template<typename T>
class ExtractScalarSubDomainFunctional3D : public BoxProcessingFunctional3D_SS<T>
{
public:
    virtual void process(Box3D domain, ScalarField3D<T>& field1, ScalarField3D<T>& field2);
    virtual ExtractScalarSubDomainFunctional3D<T>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, class Function>
class ApplyScalarFunctional3D : public BoxProcessingFunctional3D_S<T> {
public:
    ApplyScalarFunctional3D(Function f_);
    virtual void process(Box3D domain, ScalarField3D<T>& field);
    virtual ApplyScalarFunctional3D<T,Function>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    Function f;
};

template<typename T, class EvalFunction>
class EvaluateScalarFunctional3D : public BoxProcessingFunctional3D_SS<T> {
public:
    EvaluateScalarFunctional3D(EvalFunction f_);
    virtual void process(Box3D domain, ScalarField3D<T>& field,
                                       ScalarField3D<T>& result);
    virtual EvaluateScalarFunctional3D<T,EvalFunction>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    EvalFunction f;
};

template<typename T>
class A_plus_alpha_functional3D : public BoxProcessingFunctional3D_SS<T>
{
public:
    A_plus_alpha_functional3D(T alpha_);
    virtual void process(Box3D domain, ScalarField3D<T>& A,
                                       ScalarField3D<T>& result);
    virtual A_plus_alpha_functional3D<T>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    T alpha;
};

template<typename T>
class A_minus_alpha_functional3D : public BoxProcessingFunctional3D_SS<T>
{
public:
    A_minus_alpha_functional3D(T alpha_);
    virtual void process(Box3D domain, ScalarField3D<T>& A,
                                       ScalarField3D<T>& result);
    virtual A_minus_alpha_functional3D<T>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    T alpha;
};

template<typename T>
class Alpha_minus_A_functional3D : public BoxProcessingFunctional3D_SS<T>
{
public:
    Alpha_minus_A_functional3D(T alpha_);
    virtual void process(Box3D domain, ScalarField3D<T>& A,
                                       ScalarField3D<T>& result);
    virtual Alpha_minus_A_functional3D<T>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    T alpha;
};

template<typename T>
class A_times_alpha_functional3D : public BoxProcessingFunctional3D_SS<T>
{
public:
    A_times_alpha_functional3D(T alpha_);
    virtual void process(Box3D domain, ScalarField3D<T>& A,
                                       ScalarField3D<T>& result);
    virtual A_times_alpha_functional3D<T>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    T alpha;
};

template<typename T>
class A_dividedBy_alpha_functional3D : public BoxProcessingFunctional3D_SS<T>
{
public:
    A_dividedBy_alpha_functional3D(T alpha_);
    virtual void process(Box3D domain, ScalarField3D<T>& A,
                                       ScalarField3D<T>& result);
    virtual A_dividedBy_alpha_functional3D<T>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    T alpha;
};

template<typename T>
class Alpha_dividedBy_A_functional3D : public BoxProcessingFunctional3D_SS<T>
{
public:
    Alpha_dividedBy_A_functional3D(T alpha_);
    virtual void process(Box3D domain, ScalarField3D<T>& A,
                                       ScalarField3D<T>& result);
    virtual Alpha_dividedBy_A_functional3D<T>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    T alpha;
};

template<typename T>
class A_plus_alpha_inplace_functional3D : public BoxProcessingFunctional3D_S<T>
{
public:
    A_plus_alpha_inplace_functional3D(T alpha_);
    virtual void process(Box3D domain, ScalarField3D<T>& A);
    virtual A_plus_alpha_inplace_functional3D<T>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    T alpha;
};

template<typename T>
class A_minus_alpha_inplace_functional3D : public BoxProcessingFunctional3D_S<T>
{
public:
    A_minus_alpha_inplace_functional3D(T alpha_);
    virtual void process(Box3D domain, ScalarField3D<T>& A);
    virtual A_minus_alpha_inplace_functional3D<T>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    T alpha;
};

template<typename T>
class A_times_alpha_inplace_functional3D : public BoxProcessingFunctional3D_S<T>
{
public:
    A_times_alpha_inplace_functional3D(T alpha_);
    virtual void process(Box3D domain, ScalarField3D<T>& A);
    virtual A_times_alpha_inplace_functional3D<T>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    T alpha;
};

template<typename T>
class A_dividedBy_alpha_inplace_functional3D : public BoxProcessingFunctional3D_S<T>
{
public:
    A_dividedBy_alpha_inplace_functional3D(T alpha_);
    virtual void process(Box3D domain, ScalarField3D<T>& A);
    virtual A_dividedBy_alpha_inplace_functional3D<T>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    T alpha;
};


template<typename T>
class A_plus_B_functional3D : public ScalarFieldBoxProcessingFunctional3D<T>
{
public:
    virtual void process(Box3D domain, std::vector<ScalarField3D<T>*> scalarFields);
    virtual A_plus_B_functional3D<T>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class A_minus_B_functional3D : public ScalarFieldBoxProcessingFunctional3D<T>
{
public:
    virtual void process(Box3D domain, std::vector<ScalarField3D<T>*> scalarFields);
    virtual A_minus_B_functional3D<T>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class A_times_B_functional3D : public ScalarFieldBoxProcessingFunctional3D<T>
{
public:
    virtual void process(Box3D domain, std::vector<ScalarField3D<T>*> scalarFields);
    virtual A_times_B_functional3D<T>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class A_dividedBy_B_functional3D : public ScalarFieldBoxProcessingFunctional3D<T>
{
public:
    virtual void process(Box3D domain, std::vector<ScalarField3D<T>*> scalarFields);
    virtual A_dividedBy_B_functional3D<T>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class A_plus_B_inplace_functional3D : public BoxProcessingFunctional3D_SS<T>
{
public:
    virtual void process(Box3D domain, ScalarField3D<T>& A, ScalarField3D<T>& B);
    virtual A_plus_B_inplace_functional3D<T>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class A_minus_B_inplace_functional3D : public BoxProcessingFunctional3D_SS<T>
{
public:
    virtual void process(Box3D domain, ScalarField3D<T>& A, ScalarField3D<T>& B);
    virtual A_minus_B_inplace_functional3D<T>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class A_times_B_inplace_functional3D : public BoxProcessingFunctional3D_SS<T>
{
public:
    virtual void process(Box3D domain, ScalarField3D<T>& A, ScalarField3D<T>& B);
    virtual A_times_B_inplace_functional3D<T>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class A_dividedBy_B_inplace_functional3D : public BoxProcessingFunctional3D_SS<T>
{
public:
    virtual void process(Box3D domain, ScalarField3D<T>& A, ScalarField3D<T>& B);
    virtual A_dividedBy_B_inplace_functional3D<T>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
};


/* *************** PART III ****************************************** */
/* *************** Analysis of the tensor-field ********************** */
/* ******************************************************************* */

template<typename T, int nDim, class BoolMask> 
class CountTensorElementsFunctional3D : public ReductiveBoxProcessingFunctional3D_T<T,nDim>
{
public:
    CountTensorElementsFunctional3D(BoolMask boolMask_);
    virtual void process(Box3D domain, TensorField3D<T,nDim>& field);
    virtual CountTensorElementsFunctional3D<T,nDim,BoolMask>* clone() const;
    plint getCount() const;
private:
    plint countId;
    BoolMask boolMask;
};

template<typename T, int nDim>
class ExtractTensorSubDomainFunctional3D : public BoxProcessingFunctional3D_TT<T,nDim,nDim>
{
public:
    virtual void process(Box3D domain, TensorField3D<T,nDim>& field1, TensorField3D<T,nDim>& field2);
    virtual ExtractTensorSubDomainFunctional3D<T,nDim>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, int nDim>
class ExtractTensorComponentFunctional3D : public BoxProcessingFunctional3D_ST<T,nDim>
{
public:
    ExtractTensorComponentFunctional3D(int iComponent_);
    virtual void process(Box3D domain, ScalarField3D<T>& scalarField,
                                       TensorField3D<T,nDim>& tensorField );
    virtual ExtractTensorComponentFunctional3D<T,nDim>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    int iComponent;
};

template<typename T, int nDim>
class ComputeNormFunctional3D : public BoxProcessingFunctional3D_ST<T,nDim>
{
public:
    virtual void process(Box3D domain, ScalarField3D<T>& scalarField,
                                       TensorField3D<T,nDim>& tensorField);
    virtual ComputeNormFunctional3D<T,nDim>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, int nDim>
class ComputeNormSqrFunctional3D : public BoxProcessingFunctional3D_ST<T,nDim>
{
public:
    virtual void process(Box3D domain, ScalarField3D<T>& scalarField,
                                       TensorField3D<T,nDim>& tensorField);
    virtual ComputeNormSqrFunctional3D<T,nDim>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class ComputeSymmetricTensorNormFunctional3D :
    public BoxProcessingFunctional3D_ST<T,6>
{
public:
    virtual void process(Box3D domain, ScalarField3D<T>& scalarField,
                                       TensorField3D<T,6>& tensorField);
    virtual ComputeSymmetricTensorNormFunctional3D<T>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class ComputeSymmetricTensorNormSqrFunctional3D :
    public BoxProcessingFunctional3D_ST<T,6>
{
public:
    virtual void process(Box3D domain, ScalarField3D<T>& scalarField,
                                       TensorField3D<T,6>& tensorField);
    virtual ComputeSymmetricTensorNormSqrFunctional3D<T>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class ComputeSymmetricTensorTraceFunctional3D :
    public BoxProcessingFunctional3D_ST<T,6>
{
public:
    virtual void process(Box3D domain, ScalarField3D<T>& scalarField,
                                       TensorField3D<T,6>& tensorField);
    virtual ComputeSymmetricTensorTraceFunctional3D<T>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, int nDim>
class BoxBulkVorticityFunctional3D : public BoxProcessingFunctional3D_TT<T,nDim,nDim>
{
public:
    virtual void process(Box3D domain, TensorField3D<T,nDim>& velocity,
                                       TensorField3D<T,nDim>& vorticity);
    virtual BoxBulkVorticityFunctional3D<T,nDim>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, int nDim>
class BoxVorticityFunctional3D : public BoundedBoxProcessingFunctional3D_TT<T,nDim,nDim>
{
public:
    virtual void processBulk( Box3D domain, TensorField3D<T,nDim>& velocity,
                                            TensorField3D<T,nDim>& vorticity );
    virtual void processPlane( int direction, int orientation, Box3D domain,
                               TensorField3D<T,nDim>& velocity,
                               TensorField3D<T,nDim>& vorticity );
    virtual void processEdge( int plane, int normal1, int normal2, Box3D domain,
                              TensorField3D<T,nDim>& velocity,
                              TensorField3D<T,nDim>& vorticity );
    virtual void processCorner( int normalX, int normalY, int normalZ, Box3D domain,
                                TensorField3D<T,nDim>& velocity,
                                TensorField3D<T,nDim>& vorticity );
    virtual BoxVorticityFunctional3D<T,nDim>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, int nDim>
class BoxBulkStrainRateFunctional3D :
    public BoxProcessingFunctional3D_TT<T,nDim,SymmetricTensorImpl<T,nDim>::n>
{
public:
    virtual void process(Box3D domain, TensorField3D<T,nDim>& velocity,
                                       TensorField3D<T,SymmetricTensorImpl<T,nDim>::n>& S);
    virtual BoxBulkStrainRateFunctional3D<T,nDim>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, int nDim>
class BoxStrainRateFunctional3D :
    public BoundedBoxProcessingFunctional3D_TT<T,nDim,SymmetricTensorImpl<T,nDim>::n>
{
public:
    virtual void processBulk( Box3D domain, TensorField3D<T,nDim>& velocity,
                                            TensorField3D<T,SymmetricTensorImpl<T,nDim>::n>& S );
    virtual void processPlane( int direction, int orientation, Box3D domain,
                               TensorField3D<T,nDim>& velocity,
                               TensorField3D<T,SymmetricTensorImpl<T,nDim>::n>& S );
    virtual void processEdge( int plane, int normal1, int normal2, Box3D domain,
                              TensorField3D<T,nDim>& velocity,
                              TensorField3D<T,SymmetricTensorImpl<T,nDim>::n>& S );
    virtual void processCorner( int normalX, int normalY, int normalZ, Box3D domain,
                                TensorField3D<T,nDim>& velocity,
                                TensorField3D<T,SymmetricTensorImpl<T,nDim>::n>& S );
    virtual BoxStrainRateFunctional3D<T,nDim>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
};


template<typename T, int nDim>
class Tensor_A_plus_B_functional3D : public TensorFieldBoxProcessingFunctional3D<T,nDim>
{
public:
    virtual void process(Box3D domain, std::vector<TensorField3D<T,nDim>*> tensorFields);
    virtual Tensor_A_plus_B_functional3D<T,nDim>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, int nDim>
class Tensor_A_minus_B_functional3D : public TensorFieldBoxProcessingFunctional3D<T,nDim>
{
public:
    virtual void process(Box3D domain, std::vector<TensorField3D<T,nDim>*> tensorFields);
    virtual Tensor_A_minus_B_functional3D<T,nDim>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, int nDim>
class Tensor_A_times_B_functional3D : public TensorFieldBoxProcessingFunctional3D<T,nDim>
{
public:
    virtual void process(Box3D domain, std::vector<TensorField3D<T,nDim>*> tensorFields);
    virtual Tensor_A_times_B_functional3D<T,nDim>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, int nDim>
class Tensor_A_dividedBy_B_functional3D : public TensorFieldBoxProcessingFunctional3D<T,nDim>
{
public:
    virtual void process(Box3D domain, std::vector<TensorField3D<T,nDim>*> tensorFields);
    virtual Tensor_A_dividedBy_B_functional3D<T,nDim>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, int nDim>
class Tensor_A_plus_B_inplace_functional3D : public BoxProcessingFunctional3D_TT<T,nDim,nDim>
{
public:
    virtual void process(Box3D domain, TensorField3D<T,nDim>& A, TensorField3D<T,nDim>& B);
    virtual Tensor_A_plus_B_inplace_functional3D<T,nDim>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, int nDim>
class Tensor_A_minus_B_inplace_functional3D : public BoxProcessingFunctional3D_TT<T,nDim,nDim>
{
public:
    virtual void process(Box3D domain, TensorField3D<T,nDim>& A, TensorField3D<T,nDim>& B);
    virtual Tensor_A_minus_B_inplace_functional3D<T,nDim>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, int nDim>
class Tensor_A_times_B_inplace_functional3D : public BoxProcessingFunctional3D_TT<T,nDim,nDim>
{
public:
    virtual void process(Box3D domain, TensorField3D<T,nDim>& A, TensorField3D<T,nDim>& B);
    virtual Tensor_A_times_B_inplace_functional3D<T,nDim>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, int nDim>
class Tensor_A_dividedBy_B_inplace_functional3D : public BoxProcessingFunctional3D_TT<T,nDim,nDim>
{
public:
    virtual void process(Box3D domain, TensorField3D<T,nDim>& A, TensorField3D<T,nDim>& B);
    virtual Tensor_A_dividedBy_B_inplace_functional3D<T,nDim>* clone() const;
    virtual void getModificationPattern(std::vector<bool>& isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

}  // namespace plb

#endif  // DATA_ANALYSIS_FUNCTIONALS_3D_H

// Explicitly include generic algorithms which are never precompiled (not even in precompiled version)
#include "core/dataAnalysisGenerics3D.h"
