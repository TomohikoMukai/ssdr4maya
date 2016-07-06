#include <vector>
#include <map>
#include <string>
#include <boost/python.hpp>

#include "SSDR.h"

SSDR::Output ssdrOutput;

double build(int numMinBones, int numMaxInfluences, int numMaxIterations, boost::python::list bindVertices, boost::python::list animVertices, int numVertices, int numFrames)
{
    SSDR::Parameter ssdrParam;
    ssdrParam.numMaxInfluences = numMaxInfluences;
    ssdrParam.numMaxIterations = numMaxIterations;
    ssdrParam.numMinBones = numMinBones;

    SSDR::Input ssdrInput;
    ssdrInput.numVertices = numVertices;
    ssdrInput.numExamples = numFrames;
    ssdrInput.bindModel.resize(ssdrInput.numVertices);
    ssdrInput.sample.resize(ssdrInput.numVertices * numFrames);
    for (int vid = 0; vid < ssdrInput.numVertices; ++vid)
    {
        ssdrInput.bindModel[vid].x = boost::python::extract<float>(bindVertices[vid * 3 + 0]);
        ssdrInput.bindModel[vid].y = boost::python::extract<float>(bindVertices[vid * 3 + 1]);
        ssdrInput.bindModel[vid].z = boost::python::extract<float>(bindVertices[vid * 3 + 2]);
    }
    for (int f = 0; f < numFrames; ++f)
    {
        for (int vid = 0; vid < ssdrInput.numVertices; ++vid)
        {
            ssdrInput.sample[f * numVertices + vid].x = boost::python::extract<float>(animVertices[(f * numVertices + vid) * 3 + 0]);
            ssdrInput.sample[f * numVertices + vid].y = boost::python::extract<float>(animVertices[(f * numVertices + vid) * 3 + 1]);
            ssdrInput.sample[f * numVertices + vid].z = boost::python::extract<float>(animVertices[(f * numVertices + vid) * 3 + 2]);
        }
    }

    double sqe = SSDR::Decompose(ssdrOutput, ssdrInput, ssdrParam);
    return std::sqrt(sqe / ssdrInput.numVertices);
}

boost::python::list getSkinningWeight()
{
    boost::python::list skinningWeight;
    for (auto it = ssdrOutput.weight.begin(); it != ssdrOutput.weight.end(); ++it)
    {
        skinningWeight.append(*it);
    }
    return skinningWeight;
}

boost::python::list getSkinningIndex()
{
    boost::python::list skinningIndex;
    for (auto it = ssdrOutput.index.begin(); it != ssdrOutput.index.end(); ++it)
    {
        skinningIndex.append(*it);
    }
    return skinningIndex;
}

int getNumBones()
{
    return ssdrOutput.numBones;
}

boost::python::list getBoneTranslation(int boneIdx, int frame)
{
    boost::python::list retval;
    retval.append(ssdrOutput.boneTrans[frame * ssdrOutput.numBones + boneIdx].Translation().x);
    retval.append(ssdrOutput.boneTrans[frame * ssdrOutput.numBones + boneIdx].Translation().y);
    retval.append(ssdrOutput.boneTrans[frame * ssdrOutput.numBones + boneIdx].Translation().z);
    return retval;
}

boost::python::list getBoneRotation(int boneIdx, int frame)
{
    boost::python::list retval;
    retval.append(ssdrOutput.boneTrans[frame * ssdrOutput.numBones + boneIdx].Rotation().x);
    retval.append(ssdrOutput.boneTrans[frame * ssdrOutput.numBones + boneIdx].Rotation().y);
    retval.append(ssdrOutput.boneTrans[frame * ssdrOutput.numBones + boneIdx].Rotation().z);
    retval.append(ssdrOutput.boneTrans[frame * ssdrOutput.numBones + boneIdx].Rotation().w);
    return retval;
}

BOOST_PYTHON_MODULE(ssdr)
{
    boost::python::def("build", &build);
    boost::python::def("getSkinningWeight", &getSkinningWeight);
    boost::python::def("getSkinningIndex", &getSkinningIndex);
    boost::python::def("getNumBones", &getNumBones);
    boost::python::def("getBoneTranslation", &getBoneTranslation);
    boost::python::def("getBoneRotation", &getBoneRotation);
}
