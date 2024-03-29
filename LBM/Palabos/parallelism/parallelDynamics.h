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
 * Parallel dynamics object -- header file.
 */

#ifndef PARALLEL_DYNAMICS_H
#define PARALLEL_DYNAMICS_H

#include "core/globalDefs.h"
#include "core/dynamics.h"
#include "core/cell.h"


namespace plb {

#ifdef PLB_MPI_PARALLEL

template<typename T, template<typename U> class Descriptor>
class ParallelDynamics : public Dynamics<T,Descriptor> {
public:
    ParallelDynamics(std::vector<Cell<T,Descriptor>*>& baseCells_, bool hasBulkCell_);
    virtual Dynamics<T,Descriptor>* clone() const;
    virtual void collide(Cell<T,Descriptor>& cell,
                         BlockStatistics<T>& statistics_);
    virtual T computeEquilibrium(plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
                                 T jSqr, T thetaBar=T()) const;
    virtual void regularize(Cell<T,Descriptor>& cell, T rhoBar, Array<T,Descriptor<T>::d> const& j,
                            T jSqr, Array<T,SymmetricTensor<T,Descriptor>::n> const& PiNeq, T thetaBar=T() ) const;
    virtual T computeDensity(Cell<T,Descriptor> const& cell) const;
    virtual T computePressure(Cell<T,Descriptor> const& cell) const;
    virtual void computeVelocity( Cell<T,Descriptor> const& cell,
                                  Array<T,Descriptor<T>::d>& u ) const;
    virtual T computeTemperature(Cell<T,Descriptor> const& cell) const;
    virtual void computeDeviatoricStress (
        Cell<T,Descriptor> const& cell, Array<T,SymmetricTensor<T,Descriptor>::n>& PiNeq ) const;
    virtual void computeHeatFlux( Cell<T,Descriptor> const& cell,
                                  Array<T,Descriptor<T>::d>& q ) const;
    virtual void computeMoment( Cell<T,Descriptor> const& cell,
                                plint momentId, T* moment ) const;
    virtual T getOmega() const;
    virtual void setOmega(T omega_);
    virtual T getParameter(plint whichParameter) const;
    virtual void setParameter(plint whichParameter, T value);
    virtual plint numDecomposedVariables(plint order) const;
    virtual void decompose(Cell<T,Descriptor> const& cell, std::vector<T>& rawData, plint order) const;
    virtual void recompose(Cell<T,Descriptor>& cell, std::vector<T> const& rawData, plint order) const;
    virtual void rescale(std::vector<T>& rawData, T xDxInv, T xDt, plint order) const;
    virtual void getPopulations(Cell<T,Descriptor> const& cell, Array<T,Descriptor<T>::q>& f) const;
    virtual void getExternalField (
            Cell<T,Descriptor> const& cell, plint pos, plint size, T* ext ) const;
    virtual void setPopulations(Cell<T,Descriptor>& cell, Array<T,Descriptor<T>::q> const& f);
    virtual void setExternalField (
            Cell<T,Descriptor>& cell, plint pos, plint size, const T* ext);
    virtual void defineDensity(Cell<T,Descriptor>& cell, T density);
    virtual void defineVelocity(Cell<T,Descriptor>& cell, Array<T,Descriptor<T>::d> const& u);
    virtual void defineTemperature(Cell<T,Descriptor>& cell, T temperature);
    virtual void defineHeatFlux(Cell<T,Descriptor>& cell, Array<T,Descriptor<T>::d> const& q);
    virtual void defineDeviatoricStress(Cell<T,Descriptor>& cell,
                                        Array<T,SymmetricTensor<T,Descriptor>::n> const& PiNeq);
    virtual void defineMoment(Cell<T,Descriptor>& cell, plint momentId, T const* value);
    virtual T computeRhoBar(Cell<T,Descriptor> const& cell) const;
    virtual void computeRhoBarJ(Cell<T,Descriptor> const& cell,
                                T& rhoBar, Array<T,Descriptor<T>::d>& j) const;
    virtual void computeRhoBarJPiNeq(Cell<T,Descriptor> const& cell,
                                     T& rhoBar, Array<T,Descriptor<T>::d>& j,
                                     Array<T,SymmetricTensor<T,Descriptor>::n>& PiNeq) const;
    virtual T computeEbar(Cell<T,Descriptor> const& cell) const;
private:
    std::vector<Cell<T,Descriptor>*>& baseCells;
    bool hasBulkCell;
};

template<typename T, template<typename U> class Descriptor>
class ConstParallelDynamics : public Dynamics<T,Descriptor> {
public:
    ConstParallelDynamics(std::vector<Cell<T,Descriptor> const*>& baseCells_, bool hasBulkCell_);
    virtual Dynamics<T,Descriptor>* clone() const;
    virtual void collide(Cell<T,Descriptor>& cell,
                         BlockStatistics<T>& statistics_);
    virtual T computeEquilibrium(plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
                                 T jSqr, T thetaBar=T()) const;
    virtual void regularize(Cell<T,Descriptor>& cell, T rhoBar, Array<T,Descriptor<T>::d> const& j,
                            T jSqr, Array<T,SymmetricTensor<T,Descriptor>::n> const& PiNeq, T thetaBar=T() ) const;
    virtual T computeDensity(Cell<T,Descriptor> const& cell) const;
    virtual T computePressure(Cell<T,Descriptor> const& cell) const;
    virtual void computeVelocity( Cell<T,Descriptor> const& cell,
                                  Array<T,Descriptor<T>::d>& u ) const;
    virtual T computeTemperature(Cell<T,Descriptor> const& cell) const;
    virtual void computeDeviatoricStress (
        Cell<T,Descriptor> const& cell, Array<T,SymmetricTensor<T,Descriptor>::n>& PiNeq ) const;
    virtual void computeHeatFlux( Cell<T,Descriptor> const& cell,
                                  Array<T,Descriptor<T>::d>& q ) const;
    virtual void computeMoment( Cell<T,Descriptor> const& cell,
                                plint momentId, T* moment ) const;
    virtual T getOmega() const;
    virtual void setOmega(T omega_);
    virtual T getParameter(plint whichParameter) const;
    virtual void setParameter(plint whichParameter, T value);
    virtual plint numDecomposedVariables(plint order) const;
    virtual void decompose(Cell<T,Descriptor> const& cell, std::vector<T>& rawData, plint order) const;
    virtual void recompose(Cell<T,Descriptor>& cell, std::vector<T> const& rawData, plint order) const;
    virtual void rescale(std::vector<T>& rawData, T xDxInv, T xDt, plint order) const;
    virtual void getPopulations(Cell<T,Descriptor> const& cell, Array<T,Descriptor<T>::q>& f) const;
    virtual void getExternalField (
            Cell<T,Descriptor> const& cell, plint pos, plint size, T* ext ) const;
    virtual void setPopulations(Cell<T,Descriptor>& cell, Array<T,Descriptor<T>::q> const& f);
    virtual void setExternalField (
            Cell<T,Descriptor>& cell, plint pos, plint size, const T* ext);
    virtual void defineDensity(Cell<T,Descriptor>& cell, T density);
    virtual void defineVelocity(Cell<T,Descriptor>& cell, Array<T,Descriptor<T>::d> const& u);
    virtual void defineTemperature(Cell<T,Descriptor>& cell, T temperature);
    virtual void defineHeatFlux(Cell<T,Descriptor>& cell, Array<T,Descriptor<T>::d> const& q);
    virtual void defineDeviatoricStress(Cell<T,Descriptor>& cell,
                                        Array<T,SymmetricTensor<T,Descriptor>::n> const& PiNeq);
    virtual void defineMoment(Cell<T,Descriptor>& cell, plint momentId, T const* value);
    virtual T computeRhoBar(Cell<T,Descriptor> const& cell) const;
    virtual void computeRhoBarJ(Cell<T,Descriptor> const& cell,
                                T& rhoBar, Array<T,Descriptor<T>::d>& j) const;
    virtual void computeRhoBarJPiNeq(Cell<T,Descriptor> const& cell,
                                     T& rhoBar, Array<T,Descriptor<T>::d>& j,
                                     Array<T,SymmetricTensor<T,Descriptor>::n>& PiNeq) const;
    virtual T computeEbar(Cell<T,Descriptor> const& cell) const;
private:
    std::vector<Cell<T,Descriptor> const*>& baseCells;
    bool hasBulkCell;
};


#endif // PLB_MPI_PARALLEL

}  // namespace plb

#endif // defined PARALLEL_DYNAMICS_H
