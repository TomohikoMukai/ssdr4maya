# SSDR4Maya
本コードは、論文"Smooth Skinning Decomposition with Rigid Bones"のMaya2016用プラグイン実装のサンプルです。

## 実行方法＆使用例
binフォルダにあるビルド済みパッケージ一式を、Mayaのプラグインパス（MAYA_PLUG_IN_PATH）が通っているフォルダに置きます。

[![SSDR4Maya](http://img.youtube.com/vi/ZPKKR24gGbg/0.jpg)](http://www.youtube.com/watch?v=ZPKKR24gGbg)

## ビルドと実行方法
拡張ライブラリ ssdr.pyd は Visual Studio 2013 Professional プロジェクトとして作成しています。ビルドには、外部ライブラリとして [Eigen](http://eigen.tuxfamily.org/ "Eigen")、 [QuadProg++](http://quadprog.sourceforge.net/ "QuadProg++")、[Boost](http://www.boost.org/ "Boost") 、および[Maya 2016.3 Developer Kit](https://apps.autodesk.com/MAYA/ja/Detail/Index?id=6303159649350432165&appLang=en&os=Win64 "MayaDevKit)が必要です。なお、ビルドおよび実行テストには Eigen 3.2.8、QuadProg++ 1.2.1、およびBoost 1.6.1 を用いました。

ビルド手順は次の通りです。

1. Eigenのインストールフォルダにインクルードパスを通す。
2. Boostのインストールフォルダにインクルードパスを通す。
3. //MAYA_LOCATION/include および //MAYA_LOCATION/include/python2.7 フォルダにインクルードパスを通す。
4. QuadProg++をダウンロードし、下記4つのファイルをssdrフォルダにコピーする。
 * QuadProg++.hh
 * QuadProg++.cc
 * Array.hh
 * Array.cc
5. Visual Studio 2013上でビルド＆実行

## 計算パラメータの調整
SSDRの計算パラメータは、mlSsdrBuilder.py内、ssdrBuildCmdクラスの冒頭にまとめられています。
* numMaxInfluences： 各頂点当たりに割り当てられる最大ボーン数
* numMinBones： 頂点アニメーション近似に用いる最小ボーン数
* numMaxIterations： 最大反復回数

これら3つのパラメータの変更することで、それにともなう計算結果の変化を確認できると思います。現状では、最小ボーン数に大きな値を与えると計算が破綻することを確認しています。

## 開発＆テスト環境
* Windows 10 Pro
* Maya 2016 SP6
* Visual Studio 2013 Update 5
* Maya 2016.3 Developer Kit
* Eigen 3.2.8
* QuadProg++ 1.2.1
* Boost 1.6.1

## 参考文献

1. Binh Huy Le and Zhigang Deng, Smooth Skinning Decomposition with Rigid Bones, ACM Transactions on Graphics, 31 (6), (2012), 199:1-199:10.
2. Binh Huy Le and Zhigang Deng, Robust and Accurate Skeletal Rigging from Mesh Sequences, ACM Transactions on Graphics, 33 (4), (2014), 84:1-84:10.
3. 向井 智彦, スキニング分解, Computer Graphics Gems JP 2015 7章：スキニング分解（ボーンデジタル）, 2015.

## 変更履歴
1. 2016/07/06 初版公開
