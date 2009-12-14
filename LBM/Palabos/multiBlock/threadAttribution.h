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
 * Attribution of block id to thread for parallel execution -- header file.
 */

#ifndef THREAD_ATTRIBUTION_H
#define THREAD_ATTRIBUTION_H

#include "core/globalDefs.h"
#include "core/globalDefs.h"

namespace plb {

struct ThreadAttribution {
    virtual ~ThreadAttribution() { }
    virtual bool isLocal(plint id) const =0;
    virtual bool allBlocksAreLocal() const =0;
    virtual ThreadAttribution* clone() const =0;
};

struct SerialThreadAttribution : public ThreadAttribution {
    virtual bool isLocal(plint id) const;
    virtual bool allBlocksAreLocal() const;
    virtual ThreadAttribution* clone() const;
};

struct OneToOneThreadAttribution : public ThreadAttribution {
    virtual bool isLocal(plint id) const;
    virtual bool allBlocksAreLocal() const;
    virtual ThreadAttribution* clone() const;
};

}  // namespace plb

#endif
