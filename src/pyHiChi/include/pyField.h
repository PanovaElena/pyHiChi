#pragma once
#include "pyFieldInterface.h"
#include "ScalarField.h"
#include "Mapping.h"
#include "Fdtd.h"
#include "Psatd.h"
#include "Pstd.h"
#include "Mapping.h"

#include "pybind11/pybind11.h"
#include "pybind11/numpy.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace pfc
{
    // wrapper over ScalarField class object
    class pyScalarField {
    public:

        pyScalarField(ScalarField<FP>* scalarField) :scalarField(scalarField) {}

        FP* getData() const {
            return scalarField->getData();
        }

        Int3 getSize() const {
            return scalarField->getSize();
        }

        // read-only accessors
        FP get(const Int3& index) const {
            return (*scalarField)(index);
        }

        FP get(int i, int j, int k) const {
            return (*scalarField)(i, j, k);
        }

    private:

        ScalarField<FP>* scalarField;
    };


    // Base class for all pyFields
    class pyFieldBase {
    public:

        virtual FP getEx(const FP3& coords) const = 0;
        virtual FP getEy(const FP3& coords) const = 0;
        virtual FP getEz(const FP3& coords) const = 0;

            OMP_FOR_COLLAPSE()
            for (int i = 0; i < nx; i++)
                for (int j = 0; j < ny; j++)
                    for (int chunk = 0; chunk < nChunks + 1; chunk++) {
                        FP3 cEx[chunkSize], cEy[chunkSize], cEz[chunkSize];
                        FP3 cBx[chunkSize], cBy[chunkSize], cBz[chunkSize];
                        int kLast = chunk == nChunks ? chunkRem : chunkSize;
#pragma ivdep
                        for (int k = 0; k < kLast; k++) {
                            cEx[k] = derived->convertCoords(fieldEntity->ExPosition(i, j, chunk * chunkSize),
                                fieldEntity->timeShiftE);
                            cEy[k] = derived->convertCoords(fieldEntity->EyPosition(i, j, chunk * chunkSize),
                                fieldEntity->timeShiftE);
                            cEz[k] = derived->convertCoords(fieldEntity->EzPosition(i, j, chunk * chunkSize),
                                fieldEntity->timeShiftE);

                            cBx[k] = derived->convertCoords(fieldEntity->BxPosition(i, j, chunk * chunkSize),
                                fieldEntity->timeShiftB);
                            cBy[k] = derived->convertCoords(fieldEntity->ByPosition(i, j, chunk * chunkSize),
                                fieldEntity->timeShiftB);
                            cBz[k] = derived->convertCoords(fieldEntity->BzPosition(i, j, chunk * chunkSize),
                                fieldEntity->timeShiftB);
                        }
OMP_SIMD()
                        for (int k = 0; k < kLast; k++) {
                            fieldEntity->Ex(i, j, k) = fieldConf->getE(cEx[k].x, cEx[k].y, cEx[k].z).x;
                            fieldEntity->Ey(i, j, k) = fieldConf->getE(cEy[k].x, cEy[k].y, cEy[k].z).y;
                            fieldEntity->Ez(i, j, k) = fieldConf->getE(cEz[k].x, cEz[k].y, cEz[k].z).z;

                            fieldEntity->Bx(i, j, k) = fieldConf->getB(cBx[k].x, cBx[k].y, cBx[k].z).x;
                            fieldEntity->By(i, j, k) = fieldConf->getB(cBy[k].x, cBy[k].y, cBy[k].z).y;
                            fieldEntity->Bz(i, j, k) = fieldConf->getB(cBz[k].x, cBz[k].y, cBz[k].z).z;
                        }
                    }
        }

    };

    // Interface for collocated grids
    template <class TGrid, class TFieldSolver, class TDerived>
    class pyStraggeredFieldInterface<TGrid, TFieldSolver, TDerived, false>
    {
    public:

        template <class FieldConfigurationType>
        void setFieldConfiguration(const FieldConfigurationType* fieldConf) {
            TDerived* derived = static_cast<TDerived*>(this);
            pyFieldEntity<TGrid, TFieldSolver>* fieldEntity = derived->getFieldEntity();
            const int chunkSize = 32;
            const int nChunks = fieldEntity->numCells.z / chunkSize;
            const int chunkRem = fieldEntity->numCells.z % chunkSize;
            const int nx = fieldEntity->numCells.x, ny = fieldEntity->numCells.y;

            OMP_FOR_COLLAPSE()
            for (int i = 0; i < nx; i++)
                for (int j = 0; j < ny; j++)
                    for (int chunk = 0; chunk < nChunks + 1; chunk++) {
                        FP3 coords[chunkSize];
                        int kLast = chunk == nChunks ? chunkRem : chunkSize;
                        FP3 startPosition = fieldEntity->ExPosition(i, j, chunk * chunkSize);
#pragma ivdep
                        for (int k = 0; k < kLast; k++) {
                            FP3 position(startPosition.x, startPosition.y, startPosition.z + k * fieldEntity->steps.z);
                            coords[k] = derived->convertCoords(position);
                        }
OMP_SIMD()
                        for (int k = 0; k < kLast; k++) {
                            FP3 E, B;
                            fieldConf->getEB(coords[k].x, coords[k].y, coords[k].z, &E, &B);

                            fieldEntity->Ex(i, j, k + chunk * chunkSize) = E.x;
                            fieldEntity->Ey(i, j, k + chunk * chunkSize) = E.y;
                            fieldEntity->Ez(i, j, k + chunk * chunkSize) = E.z;

                            fieldEntity->Bx(i, j, k + chunk * chunkSize) = B.x;
                            fieldEntity->By(i, j, k + chunk * chunkSize) = B.y;
                            fieldEntity->Bz(i, j, k + chunk * chunkSize) = B.z;
                        }
                    }
        }

        // functions to get cross sections
        py::array_t<FP> getSlice1d(
            CoordinateEnum crossAxis1, FP pos1,
            CoordinateEnum crossAxis2, FP pos2,
            CoordinateEnum axis, FP minCoord, FP maxCoord, size_t size,
            FP(pyFieldBase::*getFieldValue)(const FP3&) const)
        {
            py::array_t<FP> res({ size });
            auto accRes = res.mutable_unchecked<1>();
            FP step = (maxCoord - minCoord) / (FP)size;
            OMP_FOR()
            for (py::ssize_t i = 0; i < size; i++) {
                FP3 coords;
                coords[(int)crossAxis1] = pos1;
                coords[(int)crossAxis1] = pos2;
                coords[(int)axis] = minCoord + step * i;
                accRes(i) = (this->*getFieldValue)(coords);
            }
            return res;
        }

        py::array_t<FP> getSlice2d(
            CoordinateEnum crossAxis, FP pos,
            CoordinateEnum axis1, FP minCoord1, FP maxCoord1, size_t size1,
            CoordinateEnum axis2, FP minCoord2, FP maxCoord2, size_t size2,
            FP(pyFieldBase::* getFieldValue)(const FP3&) const)
        {
            TDerived* derived = static_cast<TDerived*>(this);
            pyFieldEntity<TGrid, TFieldSolver>* fieldEntity = derived->getFieldEntity();
            void(*fValueField)(FP, FP, FP, FP*) = (void(*)(FP, FP, FP, FP*))_fValueField;
            const int chunkSize = 32;
            const int nChunks = fieldEntity->numCells.z / chunkSize;
            const int chunkRem = fieldEntity->numCells.z % chunkSize;
            const int nx = fieldEntity->numCells.x, ny = fieldEntity->numCells.y;

            OMP_FOR_COLLAPSE()
            for (int i = 0; i < nx; i++)
                for (int j = 0; j < ny; j++)
                    for (int chunk = 0; chunk < nChunks + 1; chunk++) {
                        FP3 coords[chunkSize];
                        int kLast = chunk == nChunks ? chunkRem : chunkSize;
                        FP3 startPosition = fieldEntity->ExPosition(i, j, chunk * chunkSize);
#pragma ivdep
                        for (int k = 0; k < kLast; k++) {
                            FP3 position(startPosition.x, startPosition.y,
                                startPosition.z + k * fieldEntity->steps.z);
                            coords[k] = derived->convertCoords(position);
                        }
OMP_SIMD()
                        for (int k = 0; k < kLast; k++) {
                            ValueField field(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
                            fValueField(coords[k].x, coords[k].y, coords[k].z, &(field.E.x));

                            fieldEntity->Ex(i, j, k + chunk * chunkSize) = field.E.x;
                            fieldEntity->Ey(i, j, k + chunk * chunkSize) = field.E.y;
                            fieldEntity->Ez(i, j, k + chunk * chunkSize) = field.E.z;

                            fieldEntity->Bx(i, j, k + chunk * chunkSize) = field.B.x;
                            fieldEntity->By(i, j, k + chunk * chunkSize) = field.B.y;
                            fieldEntity->Bz(i, j, k + chunk * chunkSize) = field.B.z;
                        }
                    }
        }

        py::array_t<FP> getSlice3d(
            CoordinateEnum axis1, FP minCoord1, FP maxCoord1, size_t size1,
            CoordinateEnum axis2, FP minCoord2, FP maxCoord2, size_t size2,
            CoordinateEnum axis3, FP minCoord3, FP maxCoord3, size_t size3,
            FP(pyFieldBase::*getFieldValue)(const FP3&) const)
        {
            py::array_t<FP> res({ size1, size2, size3 });
            auto accRes = res.mutable_unchecked<3>();
            FP step1 = (maxCoord1 - minCoord1) / (FP)size1,
                step2 = (maxCoord2 - minCoord2) / (FP)size2,
                step3 = (maxCoord3 - minCoord3) / (FP)size3;
            OMP_FOR()
            for (py::ssize_t i = 0; i < size1; i++)
                for (py::ssize_t j = 0; j < size2; j++)
                    for (py::ssize_t k = 0; k < size3; k++) {
                        FP3 coords;
                        coords[(int)axis1] = minCoord1 + step1 * i;
                        coords[(int)axis2] = minCoord2 + step2 * j;
                        coords[(int)axis3] = minCoord3 + step3 * k;
                        accRes(i, j, k) = (this->*getFieldValue)(coords);
                    }
            return res;
        }

        void applyFunction(int64_t _func)
        {
            TDerived* derived = static_cast<TDerived*>(this);
            pyFieldEntity<TGrid, TFieldSolver>* fieldEntity = derived->getFieldEntity();
            void(*fValueField)(FP, FP, FP, FP*) = (void(*)(FP, FP, FP, FP*))_func;
            const int chunkSize = 32;
            const int nChunks = fieldEntity->numCells.z / chunkSize;
            const int chunkRem = fieldEntity->numCells.z % chunkSize;
            const int nx = fieldEntity->numCells.x, ny = fieldEntity->numCells.y;

            OMP_FOR_COLLAPSE()
            for (int i = 0; i < nx; i++)
                for (int j = 0; j < ny; j++)
                    for (int chunk = 0; chunk < nChunks + 1; chunk++) {
                        FP3 coords[chunkSize];
                        int kLast = chunk == nChunks ? chunkRem : chunkSize;
                        FP3 startPosition = fieldEntity->ExPosition(i, j, chunk * chunkSize);
#pragma ivdep
                        for (int k = 0; k < kLast; k++) {
                            FP3 position(startPosition.x, startPosition.y,
                                startPosition.z + k * fieldEntity->steps.z);
                            coords[k] = derived->convertCoords(position);
                        }
OMP_SIMD()
                        for (int k = 0; k < kLast; k++) {
                            int zIndex = k + chunk * chunkSize;
                            ValueField field(fieldEntity->Ex(i, j, zIndex),
                                fieldEntity->Ey(i, j, zIndex),
                                fieldEntity->Ez(i, j, zIndex),
                                fieldEntity->Bx(i, j, zIndex),
                                fieldEntity->By(i, j, zIndex),
                                fieldEntity->Bz(i, j, zIndex));

                            fValueField(coords[k].x, coords[k].y, coords[k].z, &(field.E.x));

                            fieldEntity->Ex(i, j, zIndex) = field.E.x;
                            fieldEntity->Ey(i, j, zIndex) = field.E.y;
                            fieldEntity->Ez(i, j, zIndex) = field.E.z;
                            
                            fieldEntity->Bx(i, j, zIndex) = field.B.x;
                            fieldEntity->By(i, j, zIndex) = field.B.y;
                            fieldEntity->Bz(i, j, zIndex) = field.B.z;
                        }
                    }
        }
    };


    template<class TFieldSolver>
    class pyMappedField;

    // Simple pyField class
    template<class TFieldSolver>
    class pyField : public pyFieldInterface<TFieldSolver, pyField<TFieldSolver>>,
        public pyFieldBase
    {
        using BaseInterface = pyFieldInterface<TFieldSolver, pyField<TFieldSolver>>;
    
    public:
    
        // analytical field constructor
        pyField(FP dt) :
            grid(new typename TFieldSolver::GridType()),
            fieldSolver(new TFieldSolver(grid.get(), dt))
        {}

        // numerical field constructor
        pyField(const Int3 & numInternalCells,
            const FP3 & minCoords, const FP3 & steps, FP dt) :
            grid(new typename TFieldSolver::GridType(numInternalCells,
                minCoords, steps, numInternalCells)),
            fieldSolver(new TFieldSolver(grid.get(), dt))
        {}

        void setExyz(int64_t _fEx, int64_t _fEy, int64_t _fEz)
        {
            TDerived* derived = static_cast<TDerived*>(this);
            pyFieldEntity<TGrid, TFieldSolver>* fieldEntity = derived->getFieldEntity();
            FP(*fEx)(FP, FP, FP) = (FP(*)(FP, FP, FP))_fEx;
            FP(*fEy)(FP, FP, FP) = (FP(*)(FP, FP, FP))_fEy;
            FP(*fEz)(FP, FP, FP) = (FP(*)(FP, FP, FP))_fEz;
            OMP_FOR_COLLAPSE()
            for (int i = 0; i < fieldEntity->numCells.x; i++)
                for (int j = 0; j < fieldEntity->numCells.y; j++)
                    for (int k = 0; k < fieldEntity->numCells.z; k++)
                    {
                        FP3 cEx, cEy, cEz;
                        cEx = derived->convertCoords(fieldEntity->ExPosition(i, j, k), fieldEntity->timeShiftE);
                        cEy = derived->convertCoords(fieldEntity->EyPosition(i, j, k), fieldEntity->timeShiftE);
                        cEz = derived->convertCoords(fieldEntity->EzPosition(i, j, k), fieldEntity->timeShiftE);
                        fieldEntity->Ex(i, j, k) = fEx(cEx.x, cEx.y, cEx.z);
                        fieldEntity->Ey(i, j, k) = fEy(cEy.x, cEy.y, cEy.z);
                        fieldEntity->Ez(i, j, k) = fEz(cEz.x, cEz.y, cEz.z);
                    }
        }

        void setExyzt(int64_t _fEx, int64_t _fEy, int64_t _fEz, FP t)
        {
            TDerived* derived = static_cast<TDerived*>(this);
            pyFieldEntity<TGrid, TFieldSolver>* fieldEntity = derived->getFieldEntity();
            FP(*fEx)(FP, FP, FP, FP) = (FP(*)(FP, FP, FP, FP))_fEx;
            FP(*fEy)(FP, FP, FP, FP) = (FP(*)(FP, FP, FP, FP))_fEy;
            FP(*fEz)(FP, FP, FP, FP) = (FP(*)(FP, FP, FP, FP))_fEz;
            OMP_FOR_COLLAPSE()
            for (int i = 0; i < fieldEntity->numCells.x; i++)
                for (int j = 0; j < fieldEntity->numCells.y; j++)
                    for (int k = 0; k < fieldEntity->numCells.z; k++)
                    {
                        FP3 cEx, cEy, cEz;
                        cEx = derived->convertCoords(fieldEntity->ExPosition(i, j, k), fieldEntity->timeShiftE);
                        cEy = derived->convertCoords(fieldEntity->EyPosition(i, j, k), fieldEntity->timeShiftE);
                        cEz = derived->convertCoords(fieldEntity->EzPosition(i, j, k), fieldEntity->timeShiftE);
                        fieldEntity->Ex(i, j, k) = fEx(cEx.x, cEx.y, cEx.z, t + fieldEntity->timeShiftE);
                        fieldEntity->Ey(i, j, k) = fEy(cEy.x, cEy.y, cEy.z, t + fieldEntity->timeShiftE);
                        fieldEntity->Ez(i, j, k) = fEz(cEz.x, cEz.y, cEz.z, t + fieldEntity->timeShiftE);
                    }
        }

        void setE(int64_t _fE)
        {
            TDerived* derived = static_cast<TDerived*>(this);
            pyFieldEntity<TGrid, TFieldSolver>* fieldEntity = derived->getFieldEntity();
            FP3(*fE)(FP, FP, FP) = (FP3(*)(FP, FP, FP))_fE;
            OMP_FOR_COLLAPSE()
            for (int i = 0; i < fieldEntity->numCells.x; i++)
                for (int j = 0; j < fieldEntity->numCells.y; j++)
                    for (int k = 0; k < fieldEntity->numCells.z; k++)
                    {
                        FP3 cEx, cEy, cEz;
                        cEx = derived->convertCoords(fieldEntity->ExPosition(i, j, k), fieldEntity->timeShiftE);
                        cEy = derived->convertCoords(fieldEntity->EyPosition(i, j, k), fieldEntity->timeShiftE);
                        cEz = derived->convertCoords(fieldEntity->EzPosition(i, j, k), fieldEntity->timeShiftE);
                        fieldEntity->Ex(i, j, k) = fE(cEx.x, cEx.y, cEx.z).x;
                        fieldEntity->Ey(i, j, k) = fE(cEy.x, cEy.y, cEy.z).y;
                        fieldEntity->Ez(i, j, k) = fE(cEz.x, cEz.y, cEz.z).z;
                    }
        }

        void pySetBxyz(py::function fBx, py::function fBy, py::function fBz)
        {
            TDerived* derived = static_cast<TDerived*>(this);
            pyFieldEntity<TGrid, TFieldSolver>* fieldEntity = derived->getFieldEntity();
            for (int i = 0; i < fieldEntity->numCells.x; i++)
                for (int j = 0; j < fieldEntity->numCells.y; j++)
                    for (int k = 0; k < fieldEntity->numCells.z; k++)
                    {
                        FP3 cBx, cBy, cBz;
                        cBx = derived->convertCoords(fieldEntity->BxPosition(i, j, k), fieldEntity->timeShiftB);
                        cBy = derived->convertCoords(fieldEntity->ByPosition(i, j, k), fieldEntity->timeShiftB);
                        cBz = derived->convertCoords(fieldEntity->BzPosition(i, j, k), fieldEntity->timeShiftB);
                        fieldEntity->Bx(i, j, k) = fBx("x"_a = cBx.x, "y"_a = cBx.y, "z"_a = cBx.z).template cast<FP>();
                        fieldEntity->By(i, j, k) = fBy("x"_a = cBy.x, "y"_a = cBy.y, "z"_a = cBy.z).template cast<FP>();
                        fieldEntity->Bz(i, j, k) = fBz("x"_a = cBz.x, "y"_a = cBz.y, "z"_a = cBz.z).template cast<FP>();
                    }
        }
        FP3 getJ(const FP3& coords) const override {
            return BaseInterface::getJ(coords);
        }

        void setBxyz(int64_t _fBx, int64_t _fBy, int64_t _fBz)
        {
            TDerived* derived = static_cast<TDerived*>(this);
            pyFieldEntity<TGrid, TFieldSolver>* fieldEntity = derived->getFieldEntity();
            FP(*fBx)(FP, FP, FP) = (FP(*)(FP, FP, FP))_fBx;
            FP(*fBy)(FP, FP, FP) = (FP(*)(FP, FP, FP))_fBy;
            FP(*fBz)(FP, FP, FP) = (FP(*)(FP, FP, FP))_fBz;
            OMP_FOR_COLLAPSE()
            for (int i = 0; i < fieldEntity->numCells.x; i++)
                for (int j = 0; j < fieldEntity->numCells.y; j++)
                    for (int k = 0; k < fieldEntity->numCells.z; k++)
                    {
                        FP3 cBx, cBy, cBz;
                        cBx = derived->convertCoords(fieldEntity->BxPosition(i, j, k), fieldEntity->timeShiftB);
                        cBy = derived->convertCoords(fieldEntity->ByPosition(i, j, k), fieldEntity->timeShiftB);
                        cBz = derived->convertCoords(fieldEntity->BzPosition(i, j, k), fieldEntity->timeShiftB);
                        fieldEntity->Bx(i, j, k) = fBx(cBx.x, cBx.y, cBx.z);
                        fieldEntity->By(i, j, k) = fBy(cBy.x, cBy.y, cBy.z);
                        fieldEntity->Bz(i, j, k) = fBz(cBz.x, cBz.y, cBz.z);
                    }
        }

        void setBxyzt(int64_t _fBx, int64_t _fBy, int64_t _fBz, FP t)
        {
            TDerived* derived = static_cast<TDerived*>(this);
            pyFieldEntity<TGrid, TFieldSolver>* fieldEntity = derived->getFieldEntity();
            FP(*fBx)(FP, FP, FP, FP) = (FP(*)(FP, FP, FP, FP))_fBx;
            FP(*fBy)(FP, FP, FP, FP) = (FP(*)(FP, FP, FP, FP))_fBy;
            FP(*fBz)(FP, FP, FP, FP) = (FP(*)(FP, FP, FP, FP))_fBz;
            OMP_FOR_COLLAPSE()
            for (int i = 0; i < fieldEntity->numCells.x; i++)
                for (int j = 0; j < fieldEntity->numCells.y; j++)
                    for (int k = 0; k < fieldEntity->numCells.z; k++)
                    {
                        FP3 cBx, cBy, cBz;
                        cBx = derived->convertCoords(fieldEntity->BxPosition(i, j, k), fieldEntity->timeShiftB);
                        cBy = derived->convertCoords(fieldEntity->ByPosition(i, j, k), fieldEntity->timeShiftB);
                        cBz = derived->convertCoords(fieldEntity->BzPosition(i, j, k), fieldEntity->timeShiftB);
                        fieldEntity->Bx(i, j, k) = fBx(cBx.x, cBx.y, cBx.z, t + fieldEntity->timeShiftB);
                        fieldEntity->By(i, j, k) = fBy(cBy.x, cBy.y, cBy.z, t + fieldEntity->timeShiftB);
                        fieldEntity->Bz(i, j, k) = fBz(cBz.x, cBz.y, cBz.z, t + fieldEntity->timeShiftB);
                    }
        }

        void setB(int64_t _fB)
        {
            TDerived* derived = static_cast<TDerived*>(this);
            pyFieldEntity<TGrid, TFieldSolver>* fieldEntity = derived->getFieldEntity();
            FP3(*fB)(FP, FP, FP) = (FP3(*)(FP, FP, FP))_fB;
            OMP_FOR_COLLAPSE()
            for (int i = 0; i < fieldEntity->numCells.x; i++)
                for (int j = 0; j < fieldEntity->numCells.y; j++)
                    for (int k = 0; k < fieldEntity->numCells.z; k++)
                    {
                        FP3 cBx, cBy, cBz;
                        cBx = derived->convertCoords(fieldEntity->BxPosition(i, j, k), fieldEntity->timeShiftB);
                        cBy = derived->convertCoords(fieldEntity->ByPosition(i, j, k), fieldEntity->timeShiftB);
                        cBz = derived->convertCoords(fieldEntity->BzPosition(i, j, k), fieldEntity->timeShiftB);
                        fieldEntity->Bx(i, j, k) = fB(cBx.x, cBx.y, cBx.z).x;
                        fieldEntity->By(i, j, k) = fB(cBy.x, cBy.y, cBy.z).y;
                        fieldEntity->Bz(i, j, k) = fB(cBz.x, cBz.y, cBz.z).z;
                    }
        }

        void pySetJxyz(py::function fJx, py::function fJy, py::function fJz)
        {
            TDerived* derived = static_cast<TDerived*>(this);
            pyFieldEntity<TGrid, TFieldSolver>* fieldEntity = derived->getFieldEntity();
            for (int i = 0; i < fieldEntity->numCells.x; i++)
                for (int j = 0; j < fieldEntity->numCells.y; j++)
                    for (int k = 0; k < fieldEntity->numCells.z; k++)
                    {
                        FP3 cJx, cJy, cJz;
                        cJx = derived->convertCoords(fieldEntity->JxPosition(i, j, k), fieldEntity->timeShiftJ);
                        cJy = derived->convertCoords(fieldEntity->JyPosition(i, j, k), fieldEntity->timeShiftJ);
                        cJz = derived->convertCoords(fieldEntity->JzPosition(i, j, k), fieldEntity->timeShiftJ);
                        fieldEntity->Jx(i, j, k) = fJx("x"_a = cJx.x, "y"_a = cJx.y, "z"_a = cJx.z).template cast<FP>();
                        fieldEntity->Jy(i, j, k) = fJy("x"_a = cJy.x, "y"_a = cJy.y, "z"_a = cJy.z).template cast<FP>();
                        fieldEntity->Jz(i, j, k) = fJz("x"_a = cJz.x, "y"_a = cJz.y, "z"_a = cJz.z).template cast<FP>();
                    }
        }

        void pySetJ(py::function fJ)
        {
            TDerived* derived = static_cast<TDerived*>(this);
            pyFieldEntity<TGrid, TFieldSolver>* fieldEntity = derived->getFieldEntity();
            for (int i = 0; i < fieldEntity->numCells.x; i++)
                for (int j = 0; j < fieldEntity->numCells.y; j++)
                    for (int k = 0; k < fieldEntity->numCells.z; k++)
                    {
                        FP3 cJx, cJy, cJz;
                        cJx = derived->convertCoords(fieldEntity->JxPosition(i, j, k), fieldEntity->timeShiftJ);
                        cJy = derived->convertCoords(fieldEntity->JyPosition(i, j, k), fieldEntity->timeShiftJ);
                        cJz = derived->convertCoords(fieldEntity->JzPosition(i, j, k), fieldEntity->timeShiftJ);
                        fieldEntity->Jx(i, j, k) = fJ("x"_a = cJx.x, "y"_a = cJx.y, "z"_a = cJx.z).template cast<FP3>().x;
                        fieldEntity->Jy(i, j, k) = fJ("x"_a = cJy.x, "y"_a = cJy.y, "z"_a = cJy.z).template cast<FP3>().y;
                        fieldEntity->Jz(i, j, k) = fJ("x"_a = cJz.x, "y"_a = cJz.y, "z"_a = cJz.z).template cast<FP3>().z;
                    }
        }

        void setJxyz(int64_t _fJx, int64_t _fJy, int64_t _fJz)
        {
            TDerived* derived = static_cast<TDerived*>(this);
            pyFieldEntity<TGrid, TFieldSolver>* fieldEntity = derived->getFieldEntity();
            FP(*fJx)(FP, FP, FP) = (FP(*)(FP, FP, FP))_fJx;
            FP(*fJy)(FP, FP, FP) = (FP(*)(FP, FP, FP))_fJy;
            FP(*fJz)(FP, FP, FP) = (FP(*)(FP, FP, FP))_fJz;
            OMP_FOR_COLLAPSE()
            for (int i = 0; i < fieldEntity->numCells.x; i++)
                for (int j = 0; j < fieldEntity->numCells.y; j++)
                    for (int k = 0; k < fieldEntity->numCells.z; k++)
                    {
                        FP3 cJx, cJy, cJz;
                        cJx = derived->convertCoords(fieldEntity->JxPosition(i, j, k), fieldEntity->timeShiftJ);
                        cJy = derived->convertCoords(fieldEntity->JyPosition(i, j, k), fieldEntity->timeShiftJ);
                        cJz = derived->convertCoords(fieldEntity->JzPosition(i, j, k), fieldEntity->timeShiftJ);
                        fieldEntity->Jx(i, j, k) = fJx(cJx.x, cJx.y, cJx.z);
                        fieldEntity->Jy(i, j, k) = fJy(cJy.x, cJy.y, cJy.z);
                        fieldEntity->Jz(i, j, k) = fJz(cJz.x, cJz.y, cJz.z);
                    }
        }

        void setJxyzt(int64_t _fJx, int64_t _fJy, int64_t _fJz, FP t)
        {
            TDerived* derived = static_cast<TDerived*>(this);
            pyFieldEntity<TGrid, TFieldSolver>* fieldEntity = derived->getFieldEntity();
            FP(*fJx)(FP, FP, FP, FP) = (FP(*)(FP, FP, FP, FP))_fJx;
            FP(*fJy)(FP, FP, FP, FP) = (FP(*)(FP, FP, FP, FP))_fJy;
            FP(*fJz)(FP, FP, FP, FP) = (FP(*)(FP, FP, FP, FP))_fJz;
            OMP_FOR_COLLAPSE()
            for (int i = 0; i < fieldEntity->numCells.x; i++)
                for (int j = 0; j < fieldEntity->numCells.y; j++)
                    for (int k = 0; k < fieldEntity->numCells.z; k++)
                    {
                        FP3 cJx, cJy, cJz;
                        cJx = derived->convertCoords(fieldEntity->JxPosition(i, j, k), fieldEntity->timeShiftJ);
                        cJy = derived->convertCoords(fieldEntity->JyPosition(i, j, k), fieldEntity->timeShiftJ);
                        cJz = derived->convertCoords(fieldEntity->JzPosition(i, j, k), fieldEntity->timeShiftJ);
                        fieldEntity->Jx(i, j, k) = fJx(cJx.x, cJx.y, cJx.z, t + fieldEntity->timeShiftJ);
                        fieldEntity->Jy(i, j, k) = fJy(cJy.x, cJy.y, cJy.z, t + fieldEntity->timeShiftJ);
                        fieldEntity->Jz(i, j, k) = fJz(cJz.x, cJz.y, cJz.z, t + fieldEntity->timeShiftJ);
                    }
        }

        void setJ(int64_t _fJ)
        {
            TDerived* derived = static_cast<TDerived*>(this);
            pyFieldEntity<TGrid, TFieldSolver>* fieldEntity = derived->getFieldEntity();
            FP3(*fJ)(FP, FP, FP) = (FP3(*)(FP, FP, FP))_fJ;
            OMP_FOR_COLLAPSE()
            for (int i = 0; i < fieldEntity->numCells.x; i++)
                for (int j = 0; j < fieldEntity->numCells.y; j++)
                    for (int k = 0; k < fieldEntity->numCells.z; k++)
                    {
                        FP3 cJx, cJy, cJz;
                        cJx = derived->convertCoords(fieldEntity->JxPosition(i, j, k), fieldEntity->timeShiftJ);
                        cJy = derived->convertCoords(fieldEntity->JyPosition(i, j, k), fieldEntity->timeShiftJ);
                        cJz = derived->convertCoords(fieldEntity->JzPosition(i, j, k), fieldEntity->timeShiftJ);
                        fieldEntity->Jx(i, j, k) = fJ(cJx.x, cJx.y, cJx.z).x;
                        fieldEntity->Jy(i, j, k) = fJ(cJy.x, cJy.y, cJy.z).y;
                        fieldEntity->Jz(i, j, k) = fJ(cJz.x, cJz.y, cJz.z).z;
                    }
        }
        FP getBz(const FP3& coords) const override {
            return BaseInterface::getBz(coords);
        }

        FP getJx(const FP3& coords) const override {
            return BaseInterface::getJx(coords);
        }
        FP getJy(const FP3& coords) const override {
            return BaseInterface::getJy(coords);
        }
        FP getJz(const FP3& coords) const override {
            return BaseInterface::getJz(coords);
        }
    
        void updateFields() override {
            return BaseInterface::updateFields();
        }
        
        TGrid* getGrid() {
            return static_cast<TGrid*>(static_cast<TDerived*>(this)->getFieldEntity());
        }

        std::shared_ptr<pyField<TFieldSolver>> zoom(const FP3& minCoord,
            const FP3& zoomedGridSize, const FP3& zoomedGridStep) const {
            std::shared_ptr<pyField<TFieldSolver>> zoomedField;
            zoomedField.reset(new pyField<TFieldSolver>(
                (Int3)zoomedGridSize, minCoord, zoomedGridStep,
                this->getFieldSolver()->getTimeStep())
            );
            this->getGrid()->copyValues(zoomedField->getGrid());
            return zoomedField;
        }

        std::shared_ptr<pyScalarField> getExArray() {
            return std::make_shared<pyScalarField>(&(this->getGrid()->Ex));
        }
        std::shared_ptr<pyScalarField> getEyArray() {
            return std::make_shared<pyScalarField>(&(this->getGrid()->Ey));
        }
        std::shared_ptr<pyScalarField> getEzArray() {
            return std::make_shared<pyScalarField>(&(this->getGrid()->Ez));
        }

        std::shared_ptr<pyScalarField> getBxArray() {
            return std::make_shared<pyScalarField>(&(this->getGrid()->Bx));
        }
        std::shared_ptr<pyScalarField> getByArray() {
            return std::make_shared<pyScalarField>(&(this->getGrid()->By));
        }
        std::shared_ptr<pyScalarField> getBzArray() {
            return std::make_shared<pyScalarField>(&(this->getGrid()->Bz));
        }

        std::shared_ptr<pyScalarField> getJxArray() {
            return std::make_shared<pyScalarField>(&(this->getGrid()->Jx));
        }
        std::shared_ptr<pyScalarField> getJyArray() {
            return std::make_shared<pyScalarField>(&(this->getGrid()->Jy));
        }
    };

    
    template <class TGrid, class TFieldSolver, class TDerived, bool>
    class pyPoissonFieldSolverInterface {};

    template <class TGrid, class TFieldSolver, class TDerived>
    class pyPoissonFieldSolverInterface<TGrid, TFieldSolver, TDerived, true> {
    public:
        void convertFieldsPoissonEquation() {
            static_cast<TDerived*>(this)->getFieldEntity()->convertFieldsPoissonEquation();
        }
    
    private:
   
        std::unique_ptr<typename TFieldSolver::GridType> grid;
        std::unique_ptr<TFieldSolver> fieldSolver;
    
    };

    typedef pyField<FDTD> pyYeeField;
    typedef pyField<PSTD> pyPSTDField;
    typedef pyField<PSATD> pyPSATDField;
    typedef pyField<PSATDPoisson> pyPSATDPoissonField;
    typedef pyField<PSATDTimeStaggered> pyPSATDTimeStaggeredField;
    typedef pyField<PSATDTimeStaggeredPoisson> pyPSATDTimeStaggeredPoissonField;

    typedef pyField<AnalyticalFieldSolver> pyAnalyticalField;


    // pyField class that supports mappings
    // wrapper over pyField class object or pyMappedField class object
    template<class TFieldSolver>
    class pyMappedField : public pyFieldInterface<TFieldSolver,
        pyMappedField<TFieldSolver>>, public pyFieldBase
    {
        using BaseInterface = pyFieldInterface<TFieldSolver, pyMappedField<TFieldSolver>>;

    // Interface for all fields
    template<class TGrid, class TFieldSolver, class TDerived>
    class pyFieldInterface:
        public pyGridFieldInterface<TGrid, TFieldSolver, TDerived,
        std::is_same<TFieldSolver, NoFieldSolver>::value>,
        public pyFieldSolverInterface<TGrid, TFieldSolver, TDerived>
    {
    public:

        pyMappedField(const std::shared_ptr<pyFieldBase>& other,
            const std::shared_ptr<Mapping>& mapping) :
            pyWrappedField(other), mapping(mapping)
        {}

        inline TFieldSolver* getFieldSolver() const {
            std::shared_ptr<pyField<TFieldSolver>> pyFieldPointer =
                std::dynamic_pointer_cast<pyField<TFieldSolver>>(pyWrappedField);
            if (pyFieldPointer) return pyFieldPointer->getFieldSolver();

            std::shared_ptr<pyMappedField<TFieldSolver>> pyMappedFieldPointer =
                std::dynamic_pointer_cast<pyMappedField<TFieldSolver>>(pyWrappedField);
            if (pyMappedFieldPointer) return pyMappedFieldPointer->getFieldSolver();

            return nullptr;
        }
        inline typename TFieldSolver::GridType* getGrid() const {
            std::shared_ptr<pyField<TFieldSolver>> pyFieldPointer =
                std::dynamic_pointer_cast<pyField<TFieldSolver>>(pyWrappedField);
            if (pyFieldPointer) return pyFieldPointer->getGrid();

            std::shared_ptr<pyMappedField<TFieldSolver>> pyMappedFieldPointer =
                std::dynamic_pointer_cast<pyMappedField<TFieldSolver>>(pyWrappedField);
            if (pyMappedFieldPointer) return pyMappedFieldPointer->getGrid();

            return nullptr;
        }

        inline FP3 convertCoords(const FP3& coords) const {
            bool status = true;
            return getDirectCoords(coords, getFieldSolver()->getTime(), &status);
        }

        std::shared_ptr<pyFieldBase> applyMapping(
            const std::shared_ptr<pyFieldBase>& self,
            const std::shared_ptr<Mapping>& mapping) const override {
            return std::static_pointer_cast<pyFieldBase>(
                std::make_shared<pyMappedField<TFieldSolver>>(
                    std::static_pointer_cast<pyFieldBase>(self), mapping
                    )
                );
        }

        FP3 getE(const FP3& coords) const override {
            return this->getFieldComp3(coords, &BaseInterface::getE);
        }
        FP3 getB(const FP3& coords) const override {
            return this->getFieldComp3(coords, &BaseInterface::getB);
        }
        FP3 getJ(const FP3& coords) const override {
            return this->getFieldComp3(coords, &BaseInterface::getJ);
        }

        FP getEx(const FP3& coords) const override {
            return this->getFieldComp(coords, &BaseInterface::getEx);
        }
        FP getEy(const FP3& coords) const override {
            return this->getFieldComp(coords, &BaseInterface::getEy);
        }
        FP getEz(const FP3& coords) const override {
            return this->getFieldComp(coords, &BaseInterface::getEz);
        }

        FP getBx(const FP3& coords) const override {
            return this->getFieldComp(coords, &BaseInterface::getBx);
        }
        FP getBy(const FP3& coords) const override {
            return this->getFieldComp(coords, &BaseInterface::getBy);
        }
        FP getBz(const FP3& coords) const override {
            return this->getFieldComp(coords, &BaseInterface::getBz);
        }

        FP getJx(const FP3& coords) const override {
            return this->getFieldComp(coords, &BaseInterface::getJx);
        }
        FP getJy(const FP3& coords) const override {
            return this->getFieldComp(coords, &BaseInterface::getJy);
        }
        FP getJz(const FP3& coords) const override {
            return this->getFieldComp(coords, &BaseInterface::getJz);
        }

        void updateFields() override {
            return BaseInterface::updateFields();
        }

        void advance(FP dt) override {
            return BaseInterface::advance(dt);
        }

    protected:

        inline FP3 getDirectCoords(const FP3& coords, FP time, bool* status) const {
            FP3 coords_ = coords;
            *status = true;
            std::shared_ptr<pyMappedField<TFieldSolver>> pyMappedFieldPointer =
                std::dynamic_pointer_cast<pyMappedField<TFieldSolver>>(pyWrappedField);
            if (pyMappedFieldPointer)
                coords_ = pyMappedFieldPointer->getDirectCoords(coords_, time, status);
            bool status2 = true;
            if (mapping) coords_ = mapping->getDirectCoords(coords_, time, &status2);
            *status = *status && status2;
            return coords_;
        }

        inline FP3 getInverseCoords(const FP3& coords, FP time, bool* status) const {
            FP3 coords_ = coords;
            *status = true;
            std::shared_ptr<pyMappedField<TFieldSolver>> pyMappedFieldPointer =
                std::dynamic_pointer_cast<pyMappedField<TFieldSolver>>(pyWrappedField);
            if (pyMappedFieldPointer)
                coords_ = pyMappedFieldPointer->getInverseCoords(coords_, time, status);
            bool status2 = true;
            if (mapping) coords_ = mapping->getInverseCoords(coords_, time, &status2);
            *status = *status && status2;
            return coords_;
        }

    private:

        std::shared_ptr<pyFieldBase> pyWrappedField;
        std::shared_ptr<Mapping> mapping;

        FP getFieldComp(const FP3& coords,
            FP(BaseInterface::* getFieldValue)(const FP3&) const) const
        {
            bool status = true;
            FP time = getFieldSolver()->getTime();
            FP3 inverseCoords = getInverseCoords(coords, time, &status);
            if (!status) return 0.0;
            return (this->*getFieldValue)(inverseCoords);
        }

        FP3 getFieldComp3(const FP3& coords,
            FP3(BaseInterface::* getFieldValue)(const FP3&) const) const
        {
            bool status = true;
            FP time = getFieldSolver()->getTime();
            FP3 inverseCoords = getInverseCoords(coords, time, &status);
            if (!status) return FP3(0.0, 0.0, 0.0);
            return (this->*getFieldValue)(inverseCoords);
        }

    };

    typedef pyMappedField<FDTD> pyMappedYeeField;
    typedef pyMappedField<PSTD> pyMappedPSTDField;
    typedef pyMappedField<PSATD> pyMappedPSATDField;
    typedef pyMappedField<PSATDPoisson> pyMappedPSATDPoissonField;
    typedef pyMappedField<PSATDTimeStaggered> pyMappedPSATDTimeStaggeredField;
    typedef pyMappedField<PSATDTimeStaggeredPoisson> pyMappedPSATDTimeStaggeredPoissonField;

    typedef pyMappedField<AnalyticalFieldSolver> pyMappedAnalyticalField;


    // Object returned when summing fields
    class pySumField : public pyFieldBase
    {
    public:

        pySumField(const std::shared_ptr<pyFieldBase>& pyWrappedField1,
            const std::shared_ptr<pyFieldBase>& pyWrappedField2) :
            pyWrappedField1(pyWrappedField1), pyWrappedField2(pyWrappedField2)
        {}

        pySumField(const std::shared_ptr<pySumField>& other,
            const std::shared_ptr<Mapping>& mapping) :
            pyWrappedField1(other->pyWrappedField1->applyMapping(other->pyWrappedField1, mapping)),
            pyWrappedField2(other->pyWrappedField2->applyMapping(other->pyWrappedField2, mapping))
        {}

        std::shared_ptr<pyFieldBase> applyMapping(
            const std::shared_ptr<pyFieldBase>& self,
            const std::shared_ptr<Mapping>& mapping) const override {
            return std::static_pointer_cast<pyFieldBase>(
                std::make_shared<pySumField>(
                    std::static_pointer_cast<pySumField>(self), mapping
                    )
                );
        }

        FP3 getE(const FP3& coords) const override {
            return pyWrappedField1->getE(coords) + pyWrappedField2->getE(coords);
        }
        FP3 getB(const FP3& coords) const override {
            return pyWrappedField1->getB(coords) + pyWrappedField2->getB(coords);
        }
        FP3 getJ(const FP3& coords) const override {
            return pyWrappedField1->getJ(coords) + pyWrappedField2->getJ(coords);
        }

        FP getEx(const FP3& coords) const override {
            return pyWrappedField1->getEx(coords) + pyWrappedField2->getEx(coords);
        }
        FP getEy(const FP3& coords) const override {
            return pyWrappedField1->getEy(coords) + pyWrappedField2->getEy(coords);
        }
        FP getEz(const FP3& coords) const override {
            return pyWrappedField1->getEz(coords) + pyWrappedField2->getEz(coords);
        }

        FP getBx(const FP3& coords) const override {
            return pyWrappedField1->getBx(coords) + pyWrappedField2->getBx(coords);
        }
        FP getBy(const FP3& coords) const override {
            return pyWrappedField1->getBy(coords) + pyWrappedField2->getBy(coords);
        }
        FP getBz(const FP3& coords) const override {
            return pyWrappedField1->getBz(coords) + pyWrappedField2->getBz(coords);
        }

        FP getJx(const FP3& coords) const override {
            return pyWrappedField1->getJx(coords) + pyWrappedField2->getJx(coords);
        }
        FP getJy(const FP3& coords) const override {
            return pyWrappedField1->getJy(coords) + pyWrappedField2->getJy(coords);
        }
        FP getJz(const FP3& coords) const override {
            return pyWrappedField1->getJz(coords) + pyWrappedField2->getJz(coords);
        }

        void updateFields() override {
            pyWrappedField1->updateFields();
            pyWrappedField2->updateFields();
        }

        void advance(FP dt) override {
            pyWrappedField1->advance(dt);
            pyWrappedField2->advance(dt);
        }

    private:

        std::shared_ptr<pyFieldBase> pyWrappedField1;
        std::shared_ptr<pyFieldBase> pyWrappedField2;
    };


    // Object returned when multiplying fields by factor
    class pyMulField : public pyFieldBase {
    public:

        pyMulField(const std::shared_ptr<pyFieldBase>& pyWrappedField, FP factor) :
            pyWrappedField(pyWrappedField),
            factor(factor)
        {}

        pyMulField(const std::shared_ptr<pyMulField>& other,
            const std::shared_ptr<Mapping>& mapping) :
            pyWrappedField(other->pyWrappedField->applyMapping(other->pyWrappedField, mapping))
        {}

        std::shared_ptr<pyFieldBase> applyMapping(
            const std::shared_ptr<pyFieldBase>& self,
            const std::shared_ptr<Mapping>& mapping) const override {
            return std::static_pointer_cast<pyFieldBase>(
                std::make_shared<pyMulField>(
                    std::static_pointer_cast<pyMulField>(self), mapping
                    )
                );
        }

        FP3 getE(const FP3& coords) const override {
            return pyWrappedField->getE(coords) * factor;
        }
        FP3 getB(const FP3& coords) const override {
            return pyWrappedField->getB(coords) * factor;
        }
        FP3 getJ(const FP3& coords) const override {
            return pyWrappedField->getJ(coords) * factor;
        }

        FP getEx(const FP3& coords) const override {
            return pyWrappedField->getEx(coords) * factor;
        }
        FP getEy(const FP3& coords) const override {
            return pyWrappedField->getEy(coords) * factor;
        }
        FP getEz(const FP3& coords) const override {
            return pyWrappedField->getEz(coords) * factor;
        }

        FP getBx(const FP3& coords) const override {
            return pyWrappedField->getBx(coords) * factor;
        }
        FP getBy(const FP3& coords) const override {
            return pyWrappedField->getBy(coords) * factor;
        }
        FP getBz(const FP3& coords) const override {
            return pyWrappedField->getBz(coords) * factor;
        }

        FP getJx(const FP3& coords) const override {
            return pyWrappedField->getJx(coords) * factor;
        }
        FP getJy(const FP3& coords) const override {
            return pyWrappedField->getJy(coords) * factor;
        }
        FP getJz(const FP3& coords) const override {
            return pyWrappedField->getJz(coords) * factor;
        }

        void updateFields() override {
            pyWrappedField->updateFields();
        }

        void advance(FP dt) override {
            pyWrappedField->advance(dt);
        }

    private:

        FP factor = 1.0;
        std::shared_ptr<pyFieldBase> pyWrappedField;
    };


    template<class TGrid, class TFieldSolver>
    class pyParticleFieldInterface {
    public:
        pyParticleFieldInterface(pyField<TGrid, TFieldSolver>* _field) {
            field = _field;
        }

    protected:
        pyField<TGrid, TFieldSolver>* field;
    };


    template<class TGrid, class TFieldSolver>
    class pyParticleGenerator : public ParticleGenerator, public pyParticleFieldInterface<TGrid, TFieldSolver> {

    public:
        typedef FP WeightType;
        typedef ParticleTypes TypeIndexType;

        pyParticleGenerator(pyField<TGrid, TFieldSolver>* _field) : ParticleGenerator(), pyParticleFieldInterface<TGrid, TFieldSolver>(_field) {}

        template <class TParticleArray>
        void operator()(TParticleArray* particleArray,
            int64_t particleDensity,
            int64_t initialTemperature,
            FP init_mx, FP init_my, FP inti_mz,
            WeightType weight = 1.0,
            TypeIndexType typeIndex = ParticleTypes::Electron)
        {
            ParticleGenerator::operator()(particleArray, static_cast<TGrid*>(this->field->getFieldEntity()), (FP(*)(FP, FP, FP))particleDensity,
                (FP(*)(FP, FP, FP))initialTemperature, [](FP init_mx, FP init_my, FP init_mz) -> FP3 {return FP3(init_mx, init_my, init_mz); }, weight, typeIndex);
        }
    };

    template<class TGrid, class TFieldSolver>
    class pyCurrentDepositionCIC : public CurrentDepositionCIC<TGrid>, public pyParticleFieldInterface<TGrid, TFieldSolver> {
    public:
        pyCurrentDepositionCIC(pyField<TGrid, TFieldSolver>* _field, double _dt) : CurrentDepositionCIC<TGrid>(_dt), pyParticleFieldInterface<TGrid, TFieldSolver>(_field) {}

        template<class TParticleArray>
        void operator()(TParticleArray* particleArray) {
            CurrentDepositionCIC<TGrid>::operator()(static_cast<TGrid*>(this->field->getFieldEntity()), particleArray);
        }

        FP3 getJ(const Int3& idx) {
            return this->field->getFieldEntity()->getJ(idx);
        }
    };

    template<class TGrid, class TFieldSolver>
    class pyPeriodicalParticleBoundaryConditions : public PeriodicalParticleBoundaryConditions, public pyParticleFieldInterface<TGrid, TFieldSolver> {
    public:
        pyPeriodicalParticleBoundaryConditions(pyField<TGrid, TFieldSolver>* _field) : PeriodicalParticleBoundaryConditions(), pyParticleFieldInterface<TGrid, TFieldSolver>(_field) {}

        template<class TParticleArray>
        void update(TParticleArray* particleArray) {
            std::cout << typeid(static_cast<TGrid*>(this->field->getFieldEntity())).name() << std::endl;
            PeriodicalParticleBoundaryConditions::updateParticlePosition(static_cast<TGrid*>(this->field->getFieldEntity()), particleArray);
        }

    };

    template<class TGrid, class TFieldSolver>
    class pyInterpolation : public pyParticleFieldInterface<TGrid, TFieldSolver> {
    public:
        pyInterpolation(pyField<TGrid, TFieldSolver>* _field) : pyParticleFieldInterface<TGrid, TFieldSolver>(_field) {}

        FP getExCIC(const FP3& coords) {
            return this->field->getFieldEntity()->getExCIC(coords);
        }
        FP getEyCIC(const FP3& coords) {
            return this->field->getFieldEntity()->getEyCIC(coords);
        }
        FP getEzCIC(const FP3& coords) {
            return this->field->getFieldEntity()->getEzCIC(coords);
        }
        FP getBxCIC(const FP3& coords) {
            return this->field->getFieldEntity()->getBxCIC(coords);
        }
        FP getByCIC(const FP3& coords) {
            return this->field->getFieldEntity()->getByCIC(coords);
        }
        FP getBzCIC(const FP3& coords) {
            return this->field->getFieldEntity()->getBzCIC(coords);
        }
        FP getJxCIC(const FP3& coords) {
            return this->field->getFieldEntity()->getJxCIC(coords);
        }
        FP getJyCIC(const FP3& coords) {
            return this->field->getFieldEntity()->getJyCIC(coords);
        }
        FP getJzCIC(const FP3& coords) {
            return this->field->getFieldEntity()->getJzCIC(coords);
        }
    };


    template<class TGrid, class TFieldSolver, GridTypes gridType>
    class pyPeriodicalCurrentBC : public PeriodicalCurrentBoundaryConditions<gridType> {
    public:
        pyPeriodicalCurrentBC(pyField<TGrid, TFieldSolver>* _field) : PeriodicalCurrentBoundaryConditions<gridType>(static_cast<RealFieldSolver<gridType>*>(_field->getFieldEntity())) {}
        void update() {
            this->fieldSolver->updateDims();
            PeriodicalCurrentBoundaryConditions<gridType>::updateCurrentBoundaries();
        }
    };
}