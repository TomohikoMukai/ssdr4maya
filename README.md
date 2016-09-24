# SSDR4Maya
本コードは、Binh Huy Le氏とZhigang Deng氏による論文[Smooth Skinning Decomposition with Rigid Bones](http://graphics.cs.uh.edu/ble/papers/2012sa-ssdr/ "SSDR paper")のMaya2016用プラグイン実装のサンプルです。

 1. Binh Huy Le and Zhigang Deng, Smooth Skinning Decomposition with Rigid Bones, ACM Transactions on Graphics, 31 (6), (2012), 199:1-199:10.
 2. Binh Huy Le and Zhigang Deng, Robust and Accurate Skeletal Rigging from Mesh Sequences, ACM Transactions on Graphics, 33 (4), (2014), 84:1-84:10.
 3. 向井 智彦, スキニング分解, Computer Graphics Gems JP 2015 7章：スキニング分解（ボーンデジタル）, 2015.

## インストール方法
* binフォルダにあるビルド済みパッケージ一式を、Mayaのプラグインパス（MAYA_PLUG_IN_PATH）が通っているフォルダに置きます。

## 使用方法
このプラグインは、各フレームのシェイプをスキン＋ボーン姿勢で近似します。隣り合うフレーム同士は必ずしも滑らかに変化する必要はありません。言い換えれば、バラバラのポーズが1フレームずつ記録されているようなシーケンスでも処理可能です。

手元ではnClothシミュレーションへのボーンアニメーションへのベイク、[Mesh Data from 
Deformation Transfer for Triangle Meshes](https://people.csail.mit.edu/sumner/research/deftransfer/data.html "MeshData@CSAIL")の公開データ、およびMaya Muscleなどの少数のデータでのみテストしています。

###利用手順
1. アニメーション開始時間と終了時間を指定します。
 - 指定した時間範囲のみが処理されます。
 - アニメーションが設定されていない余分な範囲も選択されていると、計算時間が長くなったり、計算が不安定になるなどの不具合が生じます。
2. 処理対象となるシェイプを選択します。
3. メニューの[MukaiLab]->[SSDR]より処理を開始します。
4. 処理が終わったら、コマンドラインに近似誤差（RMSE）と使用しているボーン数（#Bones）が表示されます。
 - 最小ボーン数や最大インフルーエンス数などの計算パラメータは、mlSsdrBuilder.py を直接編集することで変更できます。
5. 変換後のスキンとボーンは「SsdrResult」グループにまとめられます。
 - 変換前のシェイプと同じ位置に表示されています。
 - 全てのボーンは、バインドポーズでは必ずワールド座標系の原点に配置されます。原点中心の動き（変位）がシェイプに作用するイメージです。

利用イメージは下記のYouTubeビデオもご参照下さい。

[![SSDR4Maya](http://img.youtube.com/vi/ZPKKR24gGbg/0.jpg)](http://www.youtube.com/watch?v=ZPKKR24gGbg)

### 計算パラメータの調整
SSDRの計算パラメータは、mlSsdrBuilder.py内、ssdrBuildCmdクラスの冒頭にまとめられています。

- numMaxInfluences： 各頂点当たりに割り当てられる最大ボーン数
- numMinBones： 頂点アニメーション近似に用いる最小ボーン数
- numMaxIterations： 最大反復回数

これら3つのパラメータの変更することで、それにともなう計算結果の変化を確認できると思います。現状では、最小ボーン数に大きな値を与えると計算が破綻することを確認しています。

## ビルドと実行方法
拡張ライブラリ ssdr.pyd は Visual Studio 2012 Professional Update 4 プロジェクトとして作成しています（※ソリューションファイルはVS2013にて作成）。ビルドには、外部ライブラリとして [Eigen](http://eigen.tuxfamily.org/ "Eigen")、 [QuadProg++](http://quadprog.sourceforge.net/ "QuadProg++")、[Boost](http://www.boost.org/ "Boost") 、および[Maya 2016.3 Developer Kit](https://apps.autodesk.com/MAYA/ja/Detail/Index?id=6303159649350432165&appLang=en&os=Win64 "MayaDevKit")が必要です。なお、ビルドおよび実行テストには Eigen 3.2.8、QuadProg++ 1.2.1、およびBoost 1.55.0 を用いました。

###ビルド手順

1. Eigenのインストールフォルダにインクルードパスを通す。
2. Boostのインストールフォルダにインクルードパスを通す。
3. //MAYA_LOCATION/include および //MAYA_LOCATION/include/python2.7 フォルダにインクルードパスを通す。
4. QuadProg++をダウンロードし、下記4つのファイルをssdrフォルダにコピーする。
 * QuadProg++.hh
 * QuadProg++.cc
 * Array.hh
 * Array.cc
5. Visual Studio 上でビルド＆実行

### 開発＆テスト環境
* Windows 10 Pro
* Maya 2016 SP6
* Visual Studio 2012 Update 4
* Maya 2016.3 Developer Kit
* Eigen 3.2.8
* QuadProg++ 1.2.1
* Boost 1.55.0

## 変更履歴
2016/09/24 アルファ版公開
　- Maya2016標準の開発環境に準拠
    * Visual Studio 2012 Update 4
    * boost 1.55.0
　- ssdr.pyd モジュール内の数学ライブラリをOpenMayaに置き換え
　- メニュー生成関連のpythonコードの修正

2016/07/06 初版公開
