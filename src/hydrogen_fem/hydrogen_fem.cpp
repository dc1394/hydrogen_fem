﻿/*! \file hydrogen_fem.cpp
    \brief FEMで水素原子に対するSchrödinger方程式を解くクラスの実装

    Copyright © 2019 @dc1394 All Rights Reserved.
    (but this is originally adapted by sunsetyuhi for fem1d_poisson.py from https://github.com/sunsetyuhi/fem_py/blob/master/fem1d_poisson )
    This software is released under the BSD 2-Clause License.
*/

#include "hydrogen_fem.h"
#include <cmath>                // for std::exp
#include <cstdio>               // for FILE, std::fclose, std::fopen, std::fprintf 
#include <memory>               // for std::unique_ptr
#include <boost/assert.hpp>     // for BOOST_ASSERT
#include <Eigen/Eigenvalues>    // for Eigen::GeneralizedSelfAdjointEigenSolver

namespace hydrogen_fem {
    // #region コンストラクタ

    Hydrogen_FEM::Hydrogen_FEM()
        :   hg_(Eigen::MatrixXd::Zero(NODE_TOTAL, NODE_TOTAL)),
            length_(ELE_TOTAL),
            mat_A_ele_(boost::extents[ELE_TOTAL][2][2]),
            mat_B_ele_(boost::extents[ELE_TOTAL][2][2]),
            nod_num_seg_(boost::extents[ELE_TOTAL][2]),
            node_r_ele_(boost::extents[ELE_TOTAL][2]),
            ug_(Eigen::MatrixXd::Zero(NODE_TOTAL, NODE_TOTAL))
    {
    }

    // #endregion コンストラクタ

    // #region publicメンバ関数 

    double Hydrogen_FEM::do_run()
    {
        // 入力データの生成
        make_input_data();

        // 要素行列の生成
        make_element_matrix();

        // 全体行列を生成
        make_global_matrix();
                
        // 一般化固有値問題を解く
        Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> es(hg_, ug_);

        // Eを取得
        auto const e = es.eigenvalues()[0];

        // 固有ベクトルを取得
        phi_ = es.eigenvectors().col(0);

        return e;
    }

    void Hydrogen_FEM::save_result() const
    {
        std::unique_ptr<FILE, decltype(&std::fclose)> fp(std::fopen(Hydrogen_FEM::RESULT_FILENAME, "w"), std::fclose);
        
        for (auto i = 0; i < ELE_TOTAL; i++) {
            auto const r = static_cast<double>(i) * length_[i];
            // 厳密な結果と比較
            std::fprintf(fp.get(), "%.14f, %.14f, %.14f\n", r, -phi_[i], 2.0 * std::exp(-r));
        }
    }
        
    // #endregion publicメンバ関数

    // #region privateメンバ関数

    double Hydrogen_FEM::get_A_matrix_element(std::int32_t e, double le, std::int32_t p, std::int32_t q) const
    {
        auto const ed = static_cast<double>(e);
        switch (p) {
        case 0:
            switch (q) {
            case 0:
                return 0.5 * le * (ed * ed + ed + 1.0 / 3.0) - le * le * (ed / 3.0 + 1.0 / 12.0);

            case 1:
                return -0.5 * le * (ed * ed + ed + 1.0 / 3.0) - le * le * (ed / 6.0 + 1.0 / 12.0);

            default:
                BOOST_ASSERT(!"heの添字が2以上！");
                return 0.0;
            }

        case 1:
            switch (q) {
            case 0:
                return -0.5 * le * (ed * ed + ed + 1.0 / 3.0) - le * le * (ed / 6.0 + 1.0 / 12.0);

            case 1:
                return 0.5 * le * (ed * ed + ed + 1.0 / 3.0) - le * le * (ed / 3.0 + 0.25);

            default:
                BOOST_ASSERT(!"heの添字が2以上！");
                return 0.0;
            }

        default:
            BOOST_ASSERT(!"heの添字が2以上！");
            return 0.0;
        }
    }
    
    double Hydrogen_FEM::get_B_matrix_element(std::int32_t e, double le, std::int32_t p, std::int32_t q) const
    {
        auto const ed = static_cast<double>(e);
        switch (p) {
        case 0:
            switch (q) {
            case 0:
                return le * le * le * (ed * ed / 3.0 + ed / 6.0 + 1.0 / 30.0);

            case 1:
                return le * le * le * (ed * ed / 6.0 + ed / 6.0 + 0.05);

            default:
                BOOST_ASSERT(!"ueの添字が2以上！");
                return 0.0;
            }

        case 1:
            switch (q) {
            case 0:
                return le * le * le * (ed * ed / 6.0 + ed / 6.0 + 0.05);

            case 1:
                return le * le * le * (ed * ed / 3.0 + ed / 2.0 + 0.2);

            default:
                BOOST_ASSERT(!"ueの添字が2以上！");
                return 0.0;
            }

        default:
            BOOST_ASSERT(!"ueの添字が2以上！");
            return 0.0;
        }
    }

    void Hydrogen_FEM::make_element_matrix()
    {
        // 各線分要素の長さを計算
        for (auto e = 0; e < ELE_TOTAL; e++) {
            length_[e] = std::fabs(node_r_ele_[e][1] - node_r_ele_[e][0]);
        }

        // 要素行列の各成分を計算
        for (auto e = 0; e < ELE_TOTAL; e++) {
            auto const le = length_[e];
            for (auto i = 0; i < 2; i++) {
                for (auto j = 0; j < 2; j++) {
                    mat_A_ele_[e][i][j] = get_A_matrix_element(e, le, i, j);
                    mat_B_ele_[e][i][j] = get_B_matrix_element(e, le, i, j);
                }
            }
        }
    }

    void Hydrogen_FEM::make_input_data()
    {
        std::valarray<double> node_r_glo(NODE_TOTAL);
        auto const dr = (R_MAX - R_MIN) / static_cast<double>(ELE_TOTAL);

        for (auto i = 0; i <= ELE_TOTAL; i++) {
            node_r_glo[i] = R_MIN + static_cast<double>(i) * dr;
        }

        for (auto e = 0; e < ELE_TOTAL; e++) {
            nod_num_seg_[e][0] = e;
            nod_num_seg_[e][1] = e + 1;
        }
        
        for (auto e = 0; e < ELE_TOTAL; e++) {
            for (auto i = 0; i < 2; i++) {
                node_r_ele_[e][i] = node_r_glo[nod_num_seg_[e][i]];
            }
        }
    }
    
    void Hydrogen_FEM::make_global_matrix()
    {
        for (auto e = 0; e < ELE_TOTAL; e++) {
            for (auto i = 0; i < 2; i++) {
                for (auto j = 0; j < 2; j++) {
                    hg_(nod_num_seg_[e][i], nod_num_seg_[e][j]) += mat_A_ele_[e][i][j];
                    ug_(nod_num_seg_[e][i], nod_num_seg_[e][j]) += mat_B_ele_[e][i][j];
                }
            }
        }
    }


    // #endregion privateメンバ関数
}
