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
 * Multi scale manager -- header file.
 */
#ifndef MULTI_SCALE_H
#define MULTI_SCALE_H

#include "core/globalDefs.h"
#include "core/array.h"

namespace plb {

class MultiBlockManagement2D;
class Box2D;

class MultiBlockManagement3D;
class Box3D;


template<typename T>
struct MultiScaleManager {
    virtual ~MultiScaleManager() { }
    virtual MultiScaleManager* clone() const =0;

    virtual Box2D scaleBox(Box2D box, plint nLevel) const =0;
    virtual Box3D scaleBox(Box3D box, plint nLevel) const =0;

    virtual MultiBlockManagement2D scaleMultiBlockManagement (
            MultiBlockManagement2D const& multiBlockManagement, plint nLevel ) const =0;
    virtual MultiBlockManagement3D scaleMultiBlockManagement (
            MultiBlockManagement3D const& multiBlockManagement, plint nLevel ) const =0;

    virtual void scaleVelocity(Array<T,2>& u, plint nLevel) const =0;
    virtual void scaleVelocity(Array<T,3>& u, plint nLevel) const =0;

    virtual T scaleDeltaX(plint nLevel) const =0;
    virtual T scaleDeltaT(plint nLevel) const =0;
};

template<typename T>
class PowerTwoMultiScaleManager : public MultiScaleManager<T> {
public:
    virtual Box2D scaleBox(Box2D box, plint nLevel) const;
    virtual Box3D scaleBox(Box3D box, plint nLevel) const;

    virtual MultiBlockManagement2D scaleMultiBlockManagement (
            MultiBlockManagement2D const& multiBlockManagement, plint nLevel ) const;
    virtual MultiBlockManagement3D scaleMultiBlockManagement (
            MultiBlockManagement3D const& multiBlockManagement, plint nLevel ) const;
public:
    static plint twoToTheLevel(plint nLevel);
};

template<typename T>
class ConvectiveMultiScaleManager : public PowerTwoMultiScaleManager<T> {
public:
    virtual ConvectiveMultiScaleManager<T>* clone() const;
    virtual void scaleVelocity(Array<T,2>& u, plint nLevel) const;
    virtual void scaleVelocity(Array<T,3>& u, plint nLevel) const;
    virtual T scaleDeltaX(plint nLevel) const;
    virtual T scaleDeltaT(plint nLevel) const;
};

namespace global {

    template<typename T>
    class DefaultMultiScaleManager {
    private:
        DefaultMultiScaleManager();
        ~DefaultMultiScaleManager();
        void set(MultiScaleManager<T>* newManager);
        MultiScaleManager<T> const& get() const;
    private:
        MultiScaleManager<T>* defaultManager;

    template<typename U>
    friend DefaultMultiScaleManager<U>& accessDefaultMultiScaleManager();
    template<typename U>
    friend MultiScaleManager<U> const& getDefaultMultiScaleManager();
    template<typename U>
    friend void setDefaultMultiScaleManager(MultiScaleManager<U>* newMultiScaleManager);

    };

    template<typename T>
    DefaultMultiScaleManager<T>& accessDefaultMultiScaleManager();

    template<typename T>
    MultiScaleManager<T> const& getDefaultMultiScaleManager();

    template<typename T>
    void setDefaultMultiScaleManager(MultiScaleManager<T>* newMultiScaleManager);

}  // namespace global

}  // namespace plb

#endif
