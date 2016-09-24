#include "SSDR.h"
#include <limits>
#include <algorithm>
#include <Eigen/Core>
#include <Eigen/Eigen>
#include "QuadProg.h"
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>

using namespace Eigen;

namespace SSDR {

    double ComputeApproximationErrorSq(const Output& output, const Input& input, const Parameter& param)
    {
        std::vector<double> errsq(input.numExamples);
        double rsqsum = 0;
        for (int s = 0; s < input.numExamples; ++s)
        {
            const int numVertices = input.numVertices;
            const int numIndices = param.numMaxInfluences;
            const int numBones = output.numBones;
            const int numExamples = input.numExamples;

            for (int s = 0; s < numExamples; ++s)
            {
                for (int v = 0; v < numVertices; ++v)
                {
                    MVector residual = input.sample[s * numVertices + v];
                    const MPoint& p = input.bindModel[v];
                    for (int i = 0; i < numIndices; ++i)
                    {
                        const int b = output.index[v * numIndices + i];
                        const double w = output.weight[v * numIndices + i];
                        const MTransformationMatrix& rt = output.boneTrans[s * numBones + b];
                        residual -= w * MVector(p * rt.asMatrix());
                    }
                    rsqsum += residual * residual;
                }
            }
        }
        return rsqsum;
    }

    class WeightMapUpdator
    {
    private:
        Output* output;
        const Input* input;
        const Parameter* param;
        const MatrixXd* cem;
        const MatrixXd* cim;
        const VectorXd* cev;
        const VectorXd* civ;
        const MatrixXd* scem;
        const MatrixXd* scim;
        const VectorXd* sciv;
    public:
        WeightMapUpdator(Output* output_, const Input* input_, const Parameter* param_,
            const MatrixXd* cem_, const MatrixXd* cim_, const VectorXd* cev_, const VectorXd* civ_,
            const MatrixXd* scem_, const MatrixXd* scim_, const VectorXd* sciv_)
            : output(output_), input(input_), param(param_),
            cem(cem_), cim(cim_), cev(cev_), civ(civ_),
            scem(scem_), scim(scim_), sciv(sciv_)
        {
        }
        void operator ()(const tbb::blocked_range<int>& range) const
        {
            const int numVertices = input->numVertices;
            const int numExamples = input->numExamples;
            const int numIndices = param->numMaxInfluences;
            const int numBones = output->numBones;

            MatrixXd gm = MatrixXd::Zero(numBones, numBones), sgm = MatrixXd::Zero(numIndices, numIndices);
            VectorXd gv = VectorXd::Zero(numBones), sgv = VectorXd::Zero(numIndices);

            VectorXd weight = VectorXd::Zero(numBones), w0, sweight = VectorXd::Zero(numIndices);
            MatrixXd basis = MatrixXd::Zero(numBones, numExamples * 3), sbasis = MatrixXd::Zero(numIndices, numExamples * 3);
            VectorXd targetVertex = VectorXd::Zero(numExamples * 3);

            for (int v = range.begin(); v != range.end(); ++v)
            {
                const MPoint& restVertex = input->bindModel[v];
                for (int s = 0; s < numExamples; ++s)
                {
                    for (int b = 0; b < numBones; ++b)
                    {
                        const MTransformationMatrix& rt = output->boneTrans[s * numBones + b];
                        MPoint tv = restVertex * rt.asMatrix();
                        basis(b, s * 3 + 0) = tv.x;
                        basis(b, s * 3 + 1) = tv.y;
                        basis(b, s * 3 + 2) = tv.z;
                    }
                }
                for (int s = 0; s < numExamples; ++s)
                {
                    targetVertex[s * 3 + 0] = input->sample[s * numVertices + v].x;
                    targetVertex[s * 3 + 1] = input->sample[s * numVertices + v].y;
                    targetVertex[s * 3 + 2] = input->sample[s * numVertices + v].z;
                }
                gm = basis * basis.transpose();
                gv = -basis * targetVertex;

                double qperr = SolveQP(gm, gv, *cem, *cev, *cim, *civ, weight);
                assert(qperr != std::numeric_limits<double>::infinity());

                double weightSum = 0;
                for (int i = 0; i < numIndices; ++i)
                {
                    double maxw = -std::numeric_limits<double>::max();
                    int bestbone = -1;
                    for (int b = 0; b < numBones; ++b)
                    {
                        if (weight[b] > maxw)
                        {
                            maxw = weight[b];
                            bestbone = b;
                        }
                    }
                    if (maxw <= 0)
                    {
                        break;
                    }

                    output->index[v * numIndices + i] = bestbone;
                    output->weight[v * numIndices + i] = maxw;
                    weightSum += maxw;
                    weight[bestbone] = 0;
                }
                if (weightSum < 1.0f)
                {
                    for (int j = 0; j < numExamples * 3; ++j)
                    {
                        for (int i = 0; i < numIndices; ++i)
                        {
                            sbasis(i, j) = basis(output->index[v * numIndices + i], j);
                        }
                    }
                    sgm = sbasis * sbasis.transpose();
                    sgv = -sbasis * targetVertex;
                    qperr = SolveQP(sgm, sgv, *scem, *cev, *scim, *sciv, sweight);
                    if (qperr != std::numeric_limits<double>::infinity())
                    {
                        weightSum = 0;
                        for (int i = 0; i < numIndices; ++i)
                        {
                            weightSum += sweight[i];
                            output->weight[v * numIndices + i] = sweight[i];
                        }
                    }
                }
                for (int i = 0; i < numIndices; ++i)
                {
                    output->weight[v * numIndices + i] /= weightSum;
                }
            }
        }
    };
    void UpdateWeightMap(Output& output, const Input& input, const Parameter& param)
    {
        const int numBones = output.numBones;
        const int numIndices = param.numMaxInfluences;

        MatrixXd cem = MatrixXd::Zero(1, numBones);
        MatrixXd scem = MatrixXd::Zero(1, numIndices);
        VectorXd cev = VectorXd::Zero(1);
        MatrixXd cim = MatrixXd::Zero(numBones, numBones);
        MatrixXd scim = MatrixXd::Zero(numIndices, numIndices);
        VectorXd civ = VectorXd::Zero(numBones);
        VectorXd sciv = VectorXd::Zero(numIndices);
        for (int b = 0; b < numBones; ++b)
        {
            cem(0, b) = -1.0;
            cim(b, b) = 1.0;
            civ(b) = 0;
        }
        for (int i = 0; i < numIndices; ++i)
        {
            scem(0, i) = -1.0;
            scim(i, i) = 1.0;
            sciv(i) = 0;
        }
        cev(0) = 1.0;

        tbb::parallel_for(tbb::blocked_range<int>(0, input.numVertices),
            WeightMapUpdator(&output, &input, &param,
            &cem, &cim, &cev, &civ, &scem,
            &scim, &sciv));
    }

    MTransformationMatrix CalcPointsAlignment(size_t numPoints, std::vector<MPoint>::const_iterator ps, std::vector<MPoint>::const_iterator pd)
    {
        MTransformationMatrix transform;

        MPoint cs = MVector::zero, cd = MVector::zero;
        std::vector<MPoint>::const_iterator sit = ps;
        std::vector<MPoint>::const_iterator dit = pd;
        for (size_t i = 0; i < numPoints; ++i, ++sit, ++dit)
        {
            cs += *sit;
            cd += *dit;
        }
        cs = cs / numPoints;
        cd = cd / numPoints;

        if (numPoints < 3)
        {
            transform.setTranslation(cd - cs, MSpace::kTransform);
            return transform;
        }

        Matrix<double, 4, 4> moment;
        double sxx = 0, sxy = 0, sxz = 0, syx = 0, syy = 0, syz = 0, szx = 0, szy = 0, szz = 0;
        sit = ps;
        dit = pd;
        for (size_t i = 0; i < numPoints; ++i, ++sit, ++dit)
        {
            sxx += (sit->x - cs.x) * (dit->x - cd.x);
            sxy += (sit->x - cs.x) * (dit->y - cd.y);
            sxz += (sit->x - cs.x) * (dit->z - cd.z);
            syx += (sit->y - cs.y) * (dit->x - cd.x);
            syy += (sit->y - cs.y) * (dit->y - cd.y);
            syz += (sit->y - cs.y) * (dit->z - cd.z);
            szx += (sit->z - cs.z) * (dit->x - cd.x);
            szy += (sit->z - cs.z) * (dit->y - cd.y);
            szz += (sit->z - cs.z) * (dit->z - cd.z);
        }
        moment(0, 0) = sxx + syy + szz;
        moment(0, 1) = syz - szy;        moment(1, 0) = moment(0, 1);
        moment(0, 2) = szx - sxz;        moment(2, 0) = moment(0, 2);
        moment(0, 3) = sxy - syx;        moment(3, 0) = moment(0, 3);
        moment(1, 1) = sxx - syy - szz;
        moment(1, 2) = sxy + syx;        moment(2, 1) = moment(1, 2);
        moment(1, 3) = szx + sxz;        moment(3, 1) = moment(1, 3);
        moment(2, 2) = -sxx + syy - szz;
        moment(2, 3) = syz + szy;        moment(3, 2) = moment(2, 3);
        moment(3, 3) = -sxx - syy + szz;

        if (moment.norm() > 0)
        {
            EigenSolver<Matrix<double, 4, 4>> es(moment);
            int maxi = 0;
            for (int i = 1; i < 4; ++i)
            {
                if (es.eigenvalues()(maxi).real() < es.eigenvalues()(i).real())
                {
                    maxi = i;
                }
            }
            transform.setRotationQuaternion(
                es.eigenvectors()(1, maxi).real(),
                es.eigenvectors()(2, maxi).real(),
                es.eigenvectors()(3, maxi).real(),
                es.eigenvectors()(0, maxi).real());
        }

        MPoint cs0 = cs * transform.asMatrix();
        transform.setTranslation(cd - cs0, MSpace::kTransform);
        return transform;
    }

    void ComputeExamplePoints(std::vector<MPoint>& example, int sid, int bone, const Output& output, const Input& input, const Parameter& param)
    {
        const int numVertices = input.numVertices;
        const int numIndices = param.numMaxInfluences;
        const int numBones = output.numBones;
        for (int v = 0; v < numVertices; ++v)
        {
            example[v] = input.sample[sid * numVertices + v];
            const MPoint& s = input.bindModel[v];
            for (int i = 0; i < numIndices; ++i)
            {
                const int b = output.index[v * numIndices + i];
                if (b != bone)
                {
                    const double w = output.weight[v * numIndices + i];
                    const MTransformationMatrix& at = output.boneTrans[sid * numBones + b];
                    example[v] -= w * (s * at.asMatrix());
                }
            }
        }
    }

    void SubtractCentroid(std::vector<MPoint>& model, std::vector<MPoint>& example, MPoint& corModel, MPoint& corExample, const VectorXd& weight, const Output& output, const Input& input)
    {
        const int numVertices = input.numVertices;

        double wsqsum = 0;
        corModel = MVector::zero;
        corExample = MVector::zero;
        for (int v = 0; v < numVertices; ++v)
        {
            const double w = weight[v];
            corModel += w * w * input.bindModel[v];
            corExample += w * example[v];
            wsqsum += w * w;
        }
        corModel = corModel / wsqsum;
        corExample = corExample / wsqsum;
        for (int v = 0; v < numVertices; ++v)
        {
            model[v] = weight[v] * (input.bindModel[v] - corModel);
            example[v] -= weight[v] * corExample;
        }
    }

    class BoneTransformUpdator
    {
    private:
        Output* output;
        const Input* input;
        const Parameter* param;
        const VectorXd* weight;
        int bone;
    public:
        BoneTransformUpdator(Output* output_, const Input* input_, const Parameter* param_, const VectorXd* weight_)
            : bone(0), output(output_), input(input_), param(param_), weight(weight_)
        {
        }
        void ChangeBone(int b)
        {
            bone = b;
        }
        void operator ()(const tbb::blocked_range<int>& range) const
        {
            std::vector<MPoint> model(input->numVertices), example(input->numVertices);
            for (int s = range.begin(); s != range.end(); ++s)
            {
                ComputeExamplePoints(example, s, bone, *output, *input, *param);

                MPoint corModel(0, 0, 0), corExample(0, 0, 0);
                SubtractCentroid(model, example, corModel, corExample, *weight, *output, *input);

                MTransformationMatrix transform = CalcPointsAlignment(model.size(), model.begin(), example.begin());
                MVector d = corExample - corModel * transform.asMatrix();
                transform.setTranslation(d + transform.getTranslation(MSpace::kTransform), MSpace::kTransform);
                output->boneTrans[s * output->numBones + bone] = transform;
            }
        }
    };
    void UpdateBoneTransform(Output& output, const Input& input, const Parameter& param)
    {
        const int numVertices = input.numVertices;
        const int numExamples = input.numExamples;
        const int numIndices = param.numMaxInfluences;
        const int numBones = output.numBones;

        VectorXd boneWeight = VectorXd::Zero(numVertices);
        BoneTransformUpdator transformUpdator(&output, &input, &param, &boneWeight);
        tbb::blocked_range<int> blockedRange(0, numExamples);
        for (int bone = 0; bone < numBones; ++bone)
        {
            for (int v = 0; v < numVertices; ++v)
            {
                boneWeight[v] = 0;
                for (int i = 0; i < numIndices; ++i)
                {
                    if (output.index[v * numIndices + i] == bone)
                    {
                        boneWeight[v] = output.weight[v * numIndices + i];
                        break;
                    }
                }
            }
            transformUpdator.ChangeBone(bone);
            tbb::parallel_for(blockedRange, transformUpdator);
        }
    }

    void UpdateBoneTransform(std::vector<MTransformationMatrix>& boneTrans, int numBones, const Output& output, const Input& input, const Parameter& param)
    {
        const int numVertices = input.numVertices;
        const int numExamples = input.numExamples;
        const int numIndices = param.numMaxInfluences;

        std::vector<int> numBoneVertices(numBones, 0);
        for (int v = 0; v < numVertices; ++v)
        {
            ++numBoneVertices[output.index[v * numIndices + 0]];
        }
        std::vector<int> boneVertexId(numBones, 0);
        for (int i = 1; i < numBones; ++i)
        {
            boneVertexId[i] = boneVertexId[i - 1] + numBoneVertices[i - 1];
        }
        std::vector<MPoint> skin(numVertices, MPoint(0, 0, 0));
        std::vector<MPoint> anim(numVertices * numExamples, MPoint(0, 0, 0));
        for (int v = 0; v < numVertices; ++v)
        {
            const int bs = output.index[v * numIndices + 0];
            const int bd = boneVertexId[bs];
            skin[bd] = input.bindModel[v];
            for (int s = 0; s < numExamples; ++s)
            {
                anim[s * numVertices + bd] = input.sample[s * numVertices + v];
            }
            ++boneVertexId[bs];
        }
        boneVertexId[0] = 0;
        for (int i = 1; i < numBones; ++i)
        {
            boneVertexId[i] = boneVertexId[i - 1] + numBoneVertices[i - 1];
        }
        for (int b = 0; b < numBones; ++b)
        {
            for (int s = 0; s < numExamples; ++s)
            {
                boneTrans[s * numBones + b] = CalcPointsAlignment(numBoneVertices[b], skin.begin() + boneVertexId[b], anim.begin() + s * numVertices + boneVertexId[b]);
            }
        }
    }

    int BindVertexToBone(Output& output, std::vector<MTransformationMatrix>& boneTrans, const Input& input, const Parameter& param)
    {
        const int numVertices = input.numVertices;
        const int numExamples = input.numExamples;
        const int numIndices = param.numMaxInfluences;
        int numBones = static_cast<int>(boneTrans.size() / numExamples);

        std::vector<int> numBoneVertices(numBones, 0);
        std::vector<double> vertexError(numVertices, 0);

        for (int v = 0; v < numVertices; ++v)
        {
            int bestBone = 0;
            double minErr = std::numeric_limits<double>::max();
            const MPoint& bindModelPos = input.bindModel[v];
            for (int b = 0; b < numBones; ++b)
            {
                double errsq = 0;
                for (int s = 0; s < numExamples; ++s)
                {
                    MMatrix am = boneTrans[s * numBones + b].asMatrix();
                    MVector diff = input.sample[s * numVertices + v] - bindModelPos * am;
                    errsq += diff * diff;
                }
                if (errsq < minErr)
                {
                    bestBone = b;
                    minErr = errsq;
                }
            }
            ++numBoneVertices[bestBone];
            output.index[v * numIndices + 0] = bestBone;
            vertexError[v] = minErr;
        }

        std::vector<int>::iterator smallestBoneSize = std::min_element(numBoneVertices.begin(), numBoneVertices.end());
        while (*smallestBoneSize <= 0)
        {
            const int smallestBone = static_cast<int>(smallestBoneSize - numBoneVertices.begin());
            numBoneVertices.erase(numBoneVertices.begin() + smallestBone);
            for (int s = input.numExamples - 1; s >= 0; --s)
            {
                boneTrans.erase(boneTrans.begin() + s * numBones + smallestBone);
            }
            for (int v = 0; v < numVertices; ++v)
            {
                int b = output.index[v * numIndices + 0];
                if (b >= smallestBone)
                {
                    --output.index[v * numIndices + 0];
                }
            }
            --numBones;
            smallestBoneSize = std::min_element(numBoneVertices.begin(), numBoneVertices.end());
        }
        return static_cast<int>(numBoneVertices.size());
    }

    int ClusterInitialBones(Output& output, const Input& input, const Parameter& param)
    {
        const int numVertices = input.numVertices;
        const int numExamples = input.numExamples;
        const int numIndices = param.numMaxInfluences;

        std::fill(output.index.begin(), output.index.end(), 0);
        std::fill(output.weight.begin(), output.weight.end(), 0.0);
        for (int v = 0; v < numVertices; ++v)
        {
            output.weight[v * numIndices + 0] = 1.0;
        }

        int numClusters = 1;
        std::vector<MTransformationMatrix> boneTrans(numExamples);
        UpdateBoneTransform(boneTrans, numClusters, output, input, param);

        while (numClusters < param.numMinBones)
        {
            std::vector<MPoint> clusterCenter(numClusters, MPoint(0, 0, 0));
            std::vector<int> numBoneVertices(numClusters, 0);
            for (int v = 0; v < numVertices; ++v)
            {
                const int c = output.index[v * numIndices + 0];
                clusterCenter[c] += input.bindModel[v];
                ++numBoneVertices[c];
            }
            for (int c = 0; c < numClusters; ++c)
            {
                clusterCenter[c] = clusterCenter[c] / numBoneVertices[c];
            }

            std::vector<double> maxClusterError(numClusters, -std::numeric_limits<double>::max());
            std::vector<int> mostDistantVertex(numClusters, -1);
            for (int v = 0; v < numVertices; ++v)
            {
                const int c = output.index[v * numIndices + 0];
                double sumApproxErrorSq = 0;
                for (int s = 0; s < numExamples; ++s)
                {
                    MMatrix bm = boneTrans[s * numClusters + c].asMatrix();
                    MVector diff = input.sample[s * numVertices + v] - input.bindModel[v] * bm;
                    sumApproxErrorSq += diff * diff;
                }
                MVector d = input.bindModel[v] - clusterCenter[c];
                double errSq = sumApproxErrorSq * (d * d);
                if (errSq > maxClusterError[c])
                {
                    maxClusterError[c] = errSq;
                    mostDistantVertex[c] = v;
                }
            }
            int numPrevClusters = numClusters;
            for (int c = 0; c < numPrevClusters; ++c)
            {
                output.index[mostDistantVertex[c] * numIndices + 0] = numClusters++;
                --numBoneVertices[c];
                numBoneVertices.push_back(1);
            }
            boneTrans.resize(numExamples * numClusters);

            UpdateBoneTransform(boneTrans, numClusters, output, input, param);
            numClusters = BindVertexToBone(output, boneTrans, input, param);
        }
        return numClusters;
    }

#pragma region Decompose
    double Decompose(Output& output, const Input& input, const Parameter& param)
    {
        const int numVertices = input.numVertices;
        const int numExamples = input.numExamples;
        const int numIndices = param.numMaxInfluences;

        output.index.assign(numVertices * numIndices, 0);
        output.weight.assign(numVertices * numIndices, 0.0f);

        output.numBones = ClusterInitialBones(output, input, param);
        output.boneTrans.assign(numExamples * output.numBones, MTransformationMatrix::identity);
        UpdateBoneTransform(output.boneTrans, output.numBones, output, input, param);

        for (int loop = 0; loop < param.numMaxIterations; ++loop)
        {             
            UpdateWeightMap(output, input, param);
            UpdateBoneTransform(output, input, param);
        }
        return ComputeApproximationErrorSq(output, input, param);
    }
#pragma endregion

} //namespace SSDR
