#pragma once

#include <vector>
#include <DirectXMath.h>
#include "RigidTransform.h"

namespace SSDR
{
    // 入力データ構造体
    struct Input
    {
        //! 頂点数
        int numVertices;
        //! 例示データ数
        int numExamples;
        //! バインド頂点座標（頂点数）
        std::vector<DirectX::XMFLOAT3A> bindModel;
        //! 例示形状頂点座標 (例示データ数 x 頂点数）
        std::vector<DirectX::XMFLOAT3A> sample;

        Input() : numVertices(0), numExamples(0) {}
        ~Input() {}
    };

    // 出力データ構造体
    struct Output
    {
        //! ボーン数
        int numBones;
        //! スキニングウェイト（頂点数 x インデクス数）
        std::vector<double> weight;
        //! インデクス（頂点数 x インデクス数）
        std::vector<int> index;
        //! スキニング行列（例示データ数 x 骨数）
        std::vector<RigidTransform> boneTrans;
    };

    // 計算パラメータ構造体
    struct Parameter
    {
        //! 最小ボーン数
        int numMinBones;
        //! 頂点毎にバインドされる最大ボーン数
        int numMaxInfluences;
        //! 最大反復回数
        int numMaxIterations;
    };

    extern double Decompose(Output& output, const Input& input, const Parameter& param);
    extern double ComputeApproximationErrorSq(const Output& output, const Input& input, const Parameter& param);
}
