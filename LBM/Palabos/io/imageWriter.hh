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

#ifndef IMAGE_WRITER_HH
#define IMAGE_WRITER_HH

#include "parallelism/mpiManager.h"
#include "core/globalDefs.h"
#include "io/imageWriter.h"
#include "io/colormaps.h"
#include "atomicBlock/dataField2D.h"
#include "atomicBlock/dataField3D.h"
#include "core/dataAnalysis2D.h"
#include "core/runTimeDiagnostics.h"
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <sstream>
#include <string>
#include <vector>

namespace plb {

////////// class ImageWriter ////////////////////////////////////////

template<typename T>
ImageWriter<T>::ImageWriter(std::string const& map)
    : colorRange(1024),
      numColors(1024),
      colorMap( mapGenerators::generateMap(map) )
{ }

template<typename T>
ImageWriter<T>::ImageWriter(std::string const& map, plint colorRange_, plint numColors_)
    : colorRange(colorRange_),
      numColors(numColors_),
      colorMap( mapGenerators::generateMap(map) )
{ }

template<typename T>
void ImageWriter<T>::setMap(std::string const& map, plint colorRange_, plint numColors_) {
    colorRange = colorRange_;
    numColors  = numColors_;
    colorMap   = mapGenerators::generateMap(map);
}

template<typename T>
void ImageWriter<T>::writePpm (
        std::string const& fName,
        ScalarFieldBase2D<T>& field,
        T minVal, T maxVal) const
{
    ScalarField2D<T> localField(field.getNx(), field.getNy());
    copySerializedBlock(field, localField);
    writePpmImplementation(fName, localField, minVal, maxVal);
}

template<typename T>
void ImageWriter<T>::writeGif(std::string const& fName,
                              ScalarFieldBase2D<T>& field,
                              T minVal, T maxVal) const
{
    writePpm(fName, field, minVal, maxVal);
    imageMagickPpmToGif(fName);
}

template<typename T>
void ImageWriter<T>::writeGif(std::string const& fName,
                              ScalarFieldBase2D<T>& field,
                              T minVal, T maxVal,
                              plint sizeX, plint sizeY) const
{
    writePpm(fName, field, minVal, maxVal);
    imageMagickResize(fName, sizeX, sizeY);
}

template<typename T>
void ImageWriter<T>::writeScaledGif(std::string const& fName,
                                    ScalarFieldBase2D<T>& field) const
{
    writeGif(fName, field, T(), T());
}

template<typename T>
void ImageWriter<T>::writeScaledGif(std::string const& fName,
                                    ScalarFieldBase2D<T>& field,
                                    plint sizeX, plint sizeY) const
{
    writeGif(fName, field, T(), T(), sizeX, sizeY);
}

template<typename T>
void ImageWriter<T>::writeScaledPpm(std::string const& fName,
                                    ScalarFieldBase2D<T>& field) const
{
    writePpm(fName, field, T(), T());
}


template<typename T>
void ImageWriter<T>::writePpm (
        std::string const& fName,
        ScalarFieldBase3D<T>& field,
        T minVal, T maxVal) const
{
    plint nx=0, ny=0;
    if (field.getNx()==1) {
        nx = field.getNy();
        ny = field.getNz();
    }
    else if (field.getNy()==1) {
        nx = field.getNx();
        ny = field.getNz();
    }
    else if (field.getNz()==1) {
        nx = field.getNx();
        ny = field.getNy();
    }
    else {
        return;
    }

    ScalarField2D<T> localField(nx,ny);
    serializerToUnSerializer(
            field.getBlockSerializer(field.getBoundingBox(), IndexOrdering::forward),
            localField.getBlockUnSerializer(localField.getBoundingBox(), IndexOrdering::forward) );
    writePpmImplementation(fName, localField, minVal, maxVal);
}

template<typename T>
void ImageWriter<T>::writeGif(std::string const& fName,
                              ScalarFieldBase3D<T>& field,
                              T minVal, T maxVal) const
{
    writePpm(fName, field, minVal, maxVal);
    imageMagickPpmToGif(fName);
}

template<typename T>
void ImageWriter<T>::writeGif(std::string const& fName,
                              ScalarFieldBase3D<T>& field,
                              T minVal, T maxVal,
                              plint sizeX, plint sizeY) const
{
    writePpm(fName, field, minVal, maxVal);
    imageMagickResize(fName, sizeX, sizeY);
}

template<typename T>
void ImageWriter<T>::writeScaledGif(std::string const& fName,
                                    ScalarFieldBase3D<T>& field) const
{
    writeGif(fName, field, T(), T());
}

template<typename T>
void ImageWriter<T>::writeScaledGif(std::string const& fName,
                                    ScalarFieldBase3D<T>& field,
                                    plint sizeX, plint sizeY) const
{
    writeGif(fName, field, T(), T(), sizeX, sizeY);
}

template<typename T>
void ImageWriter<T>::writeScaledPpm(std::string const& fName,
                                    ScalarFieldBase3D<T>& field) const
{
    writePpm(fName, field, T(), T());
}


template<typename T>
void ImageWriter<T>::writePpmImplementation (
        std::string const& fName,
        ScalarField2D<T>& localField,
        T minVal, T maxVal) const
{
    if (global::mpi().isMainProcessor()) {
        if (std::fabs(minVal-maxVal)<1.e-12) {
            minVal = computeMin(localField);
            maxVal = computeMax(localField);
        }
        std::string fullName = global::directories().getImageOutDir() + fName+".ppm";
        std::ofstream fout(fullName.c_str());
        fout << "P3\n";
        fout << localField.getNx() << " " << localField.getNy() << "\n";
        fout << (colorRange-1) << "\n";

        for (plint iY=localField.getNy()-1; iY>=0; --iY) {
            for (plint iX=0; iX<localField.getNx(); ++iX) {
                T outputValue = T();
                if (! (minVal==maxVal) ) {
                    outputValue = ( (double) (localField.get(iX,iY)-minVal) /
                                    (double) (maxVal-minVal) *
                                    (double) (numColors-1) / (double) numColors );
                }
                if (outputValue <   0.) outputValue = 0.;
                if (outputValue >=  1.) outputValue = (double) (numColors-1) / (double) numColors;
                rgb color = colorMap.get(outputValue);
                fout << (int) (color.r*(colorRange-1)) << " "
                     << (int) (color.g*(colorRange-1)) << " "
                     << (int) (color.b*(colorRange-1)) << "\n";
            }
        }
    }
}

template<typename T>
void ImageWriter<T>::imageMagickPpmToGif(std::string const& fName) const {
    if (global::mpi().isMainProcessor()) {
#ifdef PLB_USE_POSIX
        std::string convCommand =
            std::string("convert ") +
            global::directories().getImageOutDir() + fName + ".ppm " +
            global::directories().getImageOutDir() + fName + ".gif ";
            
        std::string rmCommand =
            std::string("/bin/rm ") +
            global::directories().getImageOutDir() + fName + ".ppm";

        plint errorConv = system(convCommand.c_str());
        plbMainProcWarning(errorConv!=0, "Error in using ImageMagick convert command.");
        plint errorRm = system(rmCommand.c_str());
        plbMainProcWarning(errorRm!=0, "Error in removing temporary ppm file.");
#endif  // PLB_USE_POSIX
    }
}

template<typename T>
void ImageWriter<T>::imageMagickResize( std::string const& fName,
                                        plint sizeX, plint sizeY) const
{
#ifdef PLB_USE_POSIX
    if (global::mpi().isMainProcessor()) {
        std::stringstream imStream;
        imStream << "convert -resize "
                 << sizeX << "x" << sizeY << " "
                 << global::directories().getImageOutDir()
                 << fName << ".ppm "
                 << global::directories().getImageOutDir()
                 << fName << ".gif";
        plint errorConv = system(imStream.str().c_str());
        plbMainProcWarning(errorConv!=0, "Error in using ImageMagick convert command.");

        std::string rmCommand =
            std::string("/bin/rm ") +
            global::directories().getImageOutDir() + fName + ".ppm";
        plint errorRm = system(rmCommand.c_str());
        plbMainProcWarning(errorRm!=0, "Error in removing temporary ppm file.");
    }
#endif  // PLB_USE_POSIX
}

}  // namespace plb

#endif
