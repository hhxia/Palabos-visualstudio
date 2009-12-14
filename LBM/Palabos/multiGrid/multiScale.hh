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
 * Multi scale manager -- generic implementation.
 */
#ifndef MULTI_SCALE_HH
#define MULTI_SCALE_HH

#include "multiGrid/multiScale.h"
#include "multiBlock/multiBlockManagement2D.h"
#include "multiBlock/multiBlockManagement3D.h"


namespace plb {


///////////////// Class PowerTwoMultiScaleManager //////////////////////////

template<typename T>
plint PowerTwoMultiScaleManager<T>::twoToTheLevel(plint nLevel) {
    plint powerOfTwo = 1;
    for (plint iL=0; iL<nLevel; ++iL) {
        powerOfTwo *= 2;
    }
    return powerOfTwo;
}

template<typename T>
Box2D PowerTwoMultiScaleManager<T>::scaleBox(Box2D box, plint nLevel) const
{
    if (nLevel>0) {
        return box.multiply(twoToTheLevel(nLevel));
    }
    else if (nLevel<0) {
        // If the fine-grid box does not fit with the coarse-grid box,
        //   let's shrink the coarse-grid box by one fine cell. This
        //   makes sure we'll never exceed the fine-grid box
        return box.divideAndFitSmaller(twoToTheLevel(-nLevel));
    }
    else {
        return box;
    }
}

template<typename T>
Box3D PowerTwoMultiScaleManager<T>::scaleBox(Box3D box, plint nLevel) const
{
    if (nLevel>0) {
        return box.multiply(twoToTheLevel(nLevel));
    }
    else if (nLevel<0) {
        // If the fine-grid box does not fit with the coarse-grid box,
        //   let's shrink the coarse-grid box by one fine cell. This
        //   makes sure we'll never exceed the fine-grid box
        return box.divideAndFitSmaller(twoToTheLevel(-nLevel));
    }
    else {
        return box;
    }
}

template<typename T>
MultiBlockManagement2D PowerTwoMultiScaleManager<T>::scaleMultiBlockManagement (
        MultiBlockManagement2D const& multiBlockManagement, plint nLevel ) const
{
    MultiBlockDistribution2D const& multiBlockDistribution = multiBlockManagement.getMultiBlockDistribution();
    MultiBlockDistribution2D rescaledBlockDistribution(scaleBox(multiBlockDistribution.getBoundingBox(), nLevel));
    for (plint iBlock=0; iBlock<multiBlockDistribution.getNumBlocks(); ++iBlock) {
        BlockParameters2D const& parameters = multiBlockDistribution.getBlockParameters(iBlock);
        rescaledBlockDistribution.addBlock (
                scaleBox(parameters.getBulk(), nLevel),
                parameters.getEnvelopeWidth(),
                parameters.getProcId() );
    }
    return MultiBlockManagement2D (
               rescaledBlockDistribution,
               multiBlockManagement.getThreadAttribution().clone(),
               multiBlockManagement.getRefinementLevel() );
}

template<typename T>
MultiBlockManagement3D PowerTwoMultiScaleManager<T>::scaleMultiBlockManagement (
        MultiBlockManagement3D const& multiBlockManagement, plint nLevel ) const
{
    MultiBlockDistribution3D const& multiBlockDistribution = multiBlockManagement.getMultiBlockDistribution();
    MultiBlockDistribution3D rescaledBlockDistribution(scaleBox(multiBlockDistribution.getBoundingBox(), nLevel));
    for (plint iBlock=0; iBlock<multiBlockDistribution.getNumBlocks(); ++iBlock) {
        BlockParameters3D const& parameters = multiBlockDistribution.getBlockParameters(iBlock);
        rescaledBlockDistribution.addBlock (
                scaleBox(parameters.getBulk(), nLevel),
                parameters.getEnvelopeWidth(),
                parameters.getProcId() );
    }
    return MultiBlockManagement3D (
               rescaledBlockDistribution,
               multiBlockManagement.getThreadAttribution().clone(),
               multiBlockManagement.getRefinementLevel() );
}


///////////////// Class ConvectiveMultiScaleManager //////////////////////////

template<typename T>
void ConvectiveMultiScaleManager<T>::scaleVelocity(Array<T,2>& u, plint nLevel) const
{
    // Velocity is scale-invariant in convective scaling.
}

template<typename T>
void ConvectiveMultiScaleManager<T>::scaleVelocity(Array<T,3>& u, plint nLevel) const
{
    // Velocity is scale-invariant in convective scaling.
}

template<typename T>
T ConvectiveMultiScaleManager<T>::scaleDeltaX(plint nLevel) const {
    if (nLevel>0) {
        return (T)1 / (T)this->twoToTheLevel(nLevel);
    }
    else if (nLevel<0) {
        return (T)this->twoToTheLevel(-nLevel);
    }
    else {
        return (T)1;
    }
}

template<typename T>
T ConvectiveMultiScaleManager<T>::scaleDeltaT(plint nLevel) const {
    // dt scale is equal to dx scale by definition in convective scaling.
    return scaleDeltaX(nLevel);
}

template<typename T>
ConvectiveMultiScaleManager<T>* ConvectiveMultiScaleManager<T>::clone() const {
    return new ConvectiveMultiScaleManager(*this);
}


namespace global {

    template<typename T>
    DefaultMultiScaleManager<T>::DefaultMultiScaleManager() {
        defaultManager = new ConvectiveMultiScaleManager<T>();
    }

    template<typename T>
    DefaultMultiScaleManager<T>::~DefaultMultiScaleManager() {
        delete defaultManager;
    }

    template<typename T>
    void DefaultMultiScaleManager<T>::set(MultiScaleManager<T>* newManager) {
        delete defaultManager;
        defaultManager = newManager;
    }

    template<typename T>
    MultiScaleManager<T> const& DefaultMultiScaleManager<T>::get() const {
        return *defaultManager;
    }


    template<typename T>
    DefaultMultiScaleManager<T>& accessDefaultMultiScaleManager() {
        static DefaultMultiScaleManager<T> defaultMultiScaleManager;
        return defaultMultiScaleManager;
    }

    template<typename T>
    MultiScaleManager<T> const& getDefaultMultiScaleManager() {
        return accessDefaultMultiScaleManager<T>().get();
    }

    template<typename T>
    void setDefaultMultiScaleManager(MultiScaleManager<T>* newMultiScaleManager) {
        accessDefaultMultiScaleManager<T>().set(newMultiScaleManager);
    }

}  // namespace global

}  // namespace plb

#endif
