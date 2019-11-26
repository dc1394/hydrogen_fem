/*! \file hydrogen_fem.cpp
    \brief FEMで水素原子に対するSchrödinger方程式を解くクラスの実装

    Copyright © 2019 @dc1394 All Rights Reserved.
    This software is released under the BSD 2-Clause License.
*/

#include "hydrogen_fem.h"
#include <cmath>                // for std::exp
#include <cstdio>               // for FILE, std::fclose, std::fopen, std::fprintf 
#include <memory>               // for std::unique_ptr
#include <boost/assert.hpp>     // for BOOST_ASSERT
#include <Eigen/Eigenvalues>    // for Eigen::GeneralizedSelfAdjointEigenSolver
#include <iostream>

namespace hydrogen_fem {
    // #region コンストラクタ

    Hydrogen_FEM::Hydrogen_FEM()
        :   hg_(Eigen::MatrixXd::Zero(NODE_TOTAL, NODE_TOTAL)),
            length_(ELE_TOTAL),
            mat_A_ele_(boost::extents[ELE_TOTAL][2][2]),
            mat_B_ele_(boost::extents[ELE_TOTAL][2][2]),
            node_num_glo_in_seg_ele_(boost::extents[ELE_TOTAL][2]),
            node_x_ele_(boost::extents[ELE_TOTAL][2]),
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
        c_ = es.eigenvectors().col(0);

        return e;
    }

    void Hydrogen_FEM::save_result() const
    {
        //std::unique_ptr<FILE, decltype(&std::fclose)> fp(std::fopen(Hydrogen_FEM::RESULTFILENAME, "w"), std::fclose);

        for (auto i = 0; i < ELE_TOTAL; i++) {
            auto const r = static_cast<double>(i) * length_[i];
        //    fprintf(fp.get(), "%.14f, %.14f, %.14f\n", r, -c_[i], 2.0 * std::exp(-r));
        }
    }
        
    // #endregion publicメンバ関数

    // #region privateメンバ関数

    double Hydrogen_FEM::get_A_matrix_element(double dh, std::int32_t n, std::int32_t p, std::int32_t q) const
    {
        auto const nd = static_cast<double>(n);
        switch (p) {
        case 0:
            switch (q) {
            case 0:
                return 0.5 * dh * (nd * nd + nd + 1.0 / 3.0) - dh * dh * (nd / 3.0 + 1.0 / 12.0);

            case 1:
                return -0.5 * dh * (nd * nd + nd + 1.0 / 3.0) - dh * dh * (nd / 6.0 + 1.0 / 12.0);

            default:
                BOOST_ASSERT(!"heの添字が2以上！");
                return 0.0;
            }

        case 1:
            switch (q) {
            case 0:
                return -0.5 * dh * (nd * nd + nd + 1.0 / 3.0) - dh * dh * (nd / 6.0 + 1.0 / 12.0);

            case 1:
                return 0.5 * dh * (nd * nd + nd + 1.0 / 3.0) - dh * dh * (nd / 3.0 + 0.25);

            default:
                BOOST_ASSERT(!"heの添字が2以上！");
                return 0.0;
            }

        default:
            BOOST_ASSERT(!"heの添字が2以上！");
            return 0.0;
        }
    }
    
    double Hydrogen_FEM::get_B_matrix_element(double dh, std::int32_t n, std::int32_t p, std::int32_t q) const
    {
        auto const nd = static_cast<double>(n);
        switch (p) {
        case 0:
            switch (q) {
            case 0:
                return dh * dh * dh * (nd * nd / 3.0 + nd / 6.0 + 1.0 / 30.0);

            case 1:
                return dh * dh * dh * (nd * nd / 6.0 + nd / 6.0 + 1.0 / 20.0);

            default:
                BOOST_ASSERT(!"ueの添字が2以上！");
                return 0.0;
            }

        case 1:
            switch (q) {
            case 0:
                return dh * dh * dh * (nd * nd / 6.0 + nd / 6.0 + 1.0 / 20.0);

            case 1:
                return dh * dh * dh * (nd * nd / 3.0 + nd / 2.0 + 1.0 / 5.0);

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
            length_[e] = std::fabs(node_x_ele_[e][1] - node_x_ele_[e][0]);
        }

        // 要素行列の各成分を計算
        for (auto e = 0; e < ELE_TOTAL; e++) {
            auto const dh = length_[e];
            for (auto i = 0; i < 2; i++) {
                for (auto j = 0; j < 2; j++) {
                    mat_A_ele_[e][i][j] = get_A_matrix_element(dh, e, i, j);
                    mat_B_ele_[e][i][j] = get_B_matrix_element(dh, e, i, j);
                }
            }
        }
    }

    void Hydrogen_FEM::make_input_data()
    {
        std::valarray<double> node_x_glo(NODE_TOTAL);
        auto const dr = (R_MAX - R_MIN) / static_cast<double>(ELE_TOTAL);

        for (auto i = 0; i <= ELE_TOTAL; i++) {
            node_x_glo[i] = R_MIN + static_cast<double>(i) * dr;
        }

        for (auto e = 0; e < ELE_TOTAL; e++) {
            node_num_glo_in_seg_ele_[e][0] = e;
            node_num_glo_in_seg_ele_[e][1] = e + 1;
        }
        
        for (auto e = 0; e < ELE_TOTAL; e++) {
            for (auto i = 0; i < 2; i++) {
                node_x_ele_[e][i] = node_x_glo[node_num_glo_in_seg_ele_[e][i]];
            }
        }
    }
    
    void Hydrogen_FEM::make_global_matrix()
    {
        for (auto e = 0; e < ELE_TOTAL; e++) {
            for (auto i = 0; i < 2; i++) {
                for (auto j = 0; j < 2; j++) {
                    hg_(node_num_glo_in_seg_ele_[e][i], node_num_glo_in_seg_ele_[e][j]) += mat_A_ele_[e][i][j];
                    ug_(node_num_glo_in_seg_ele_[e][i], node_num_glo_in_seg_ele_[e][j]) += mat_B_ele_[e][i][j];
                }
            }
        }
    }


    // #endregion privateメンバ関数
}
