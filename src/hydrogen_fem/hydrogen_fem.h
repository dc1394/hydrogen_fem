/*! \file hydrogen_fem.h
    \brief FEMで水素原子に対するSchrödinger方程式を解くクラスの宣言

    Copyright © 2019 @dc1394 All Rights Reserved.
    This software is released under the BSD 2-Clause License.    
*/

#ifndef _HYDROGEN_FEM_H_
#define _HYDROGEN_FEM_H_

#pragma once

#include <cstdint>      // for std::int32_t
#include <Eigen/Core>   // for Eigen::MatrixXd

namespace hydrogen_fem {
    //! A class.
    /*!
        FEMで水素原子に対するSchrödinger方程式を解くクラス
    */
    class Hydrogen_FEM final {
        // #region コンストラクタ・デストラクタ

    public:
        //! A constructor.
        /*!
            唯一のコンストラクタ
        */
        Hydrogen_FEM();

        //! A destructor.
        /*!
            デストラクタ
        */
        ~Hydrogen_FEM() = default;

        // #region publicメンバ関数

        //! A public member function.
        /*!
            実際に計算を行い、水素原子の基底状態のエネルギー固有値を返す
            \return 水素原子の基底状態のエネルギー固有値
        */
        double do_run();

        //! A public member function.
        /*!
            計算結果をファイルに出力する
        */
        void save_result() const;
        
        // #endregion publicメンバ関数

        // #region privateメンバ関数

    private:
        //! A private member function (const).
        /*!
            左辺の要素マトリックスを計算する
            \param n n番目の要素のn
            \param p 左辺の要素マトリックスの行
            \param q 左辺の要素マトリックスの列
            \return 左辺の要素マトリックスの要素
        */
        double get_he_matrix_element(std::int32_t n, std::int32_t p, std::int32_t q) const;

        //! A private member function (const).
        /*!
            右辺の要素マトリックスを計算する
            \param n n番目の要素のn
            \param p 右辺の要素マトリックスの行
            \param q 右辺の要素マトリックスの列
            \return 右辺の要素マトリックスの要素
        */
        double get_ue_matrix_element(std::int32_t n, std::int32_t p, std::int32_t q) const;
                
        //! A private member function.
        /*!
            左辺の全体マトリックスを生成する
        */
        void make_hg_matrix();

        //! A private member function.
        /*!
            右辺の全体マトリックスを生成する
        */
        void make_ug_matrix();
        
        // #endregion privateメンバ関数

        // #region メンバ変数

    public:
        //! A public member variable (constant expression).
        /*!
            出力ファイル名
        */
        static auto constexpr RESULTFILENAME = "result.csv";

    private:
        //! A private member variable (constant expression).
        /*!
            Δh
        */
        static auto constexpr DH = 0.02;

        //! A private member variable (constant expression).
        /*!
            全体マトリックスの行列数
        */
        static auto constexpr N = 1000;
                        
        //! A private member variable.
        /*!
            固有ベクトル
        */
        Eigen::VectorXd c_;

        //! A private member variable.
        /*!
            左辺の全体マトリックス
        */
        Eigen::MatrixXd hg_;

        //! A private member variable.
        /*!
            右辺の全体マトリックス
        */
        Eigen::MatrixXd ug_;

        // #region 禁止されたコンストラクタ・メンバ関数

    public:
        //! A public copy constructor (deleted).
        /*!
            コピーコンストラクタ（禁止）
            \param dummy コピー元のオブジェクト（未使用）
        */
        Hydrogen_FEM(Hydrogen_FEM const & dummy) = delete;

        //! A public member function (deleted).
        /*!
            operator=()の宣言（禁止）
            \param dummy コピー元のオブジェクト（未使用）
            \return コピー元のオブジェクト
        */
        Hydrogen_FEM & operator=(Hydrogen_FEM const & dummy) = delete;

        // #endregion 禁止されたコンストラクタ・メンバ関数
    };
        
}

#endif  // _HYDROGEN_FEM_H_
