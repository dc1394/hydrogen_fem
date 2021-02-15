/*! \file hydrogen_fem.cpp
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
            node_num_seg_(boost::extents[ELE_TOTAL][2]),
            node_r_ele_(boost::extents[ELE_TOTAL][2]),
            node_r_glo_(NODE_TOTAL),
            ug_(Eigen::MatrixXd::Zero(NODE_TOTAL, NODE_TOTAL))
    {
    }

    // #endregion コンストラクタ

    // #region publicメンバ関数 

    double Hydrogen_FEM::do_run()
    {
        // 各種データの生成
        make_data();

        // 要素行列の生成
        make_element_matrix();

        // 全体行列を生成
        make_global_matrix();
        
        // 境界条件処理を行う
        boundary_conditions();

        // 一般化固有値問題を解く
        Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> es(hg_, ug_);

        // エネルギー固有値を取得
        eigenval_ = es.eigenvalues();

        // 基底状態の固有関数（波動関数）（波動関数）を取得
        phi_ = es.eigenvectors().col(0);

        // 固有ベクトル（波動関数）のN要素目を追加
        phi_.resize(NODE_TOTAL);

        // 固有ベクトル（波動関数）を規格化
        normalize();

        return eigenval_[0];
    }

    void Hydrogen_FEM::save_result() const
    {
        std::unique_ptr<FILE, decltype(&std::fclose)> fp_eigenfunc(std::fopen(Hydrogen_FEM::EIGENFUNC_FILENAME, "w"), std::fclose);
        
        for (auto i = 0; i < NODE_TOTAL; i++) {
            auto const r = static_cast<double>(i) * length_[i];
            // 厳密な結果と比較
            std::fprintf(fp_eigenfunc.get(), "%.14f, %.14f, %.14f\n", r, phi_[i], 2.0 * std::exp(-r));
        }

        std::unique_ptr<FILE, decltype(&std::fclose)> fp_eigenval(std::fopen(Hydrogen_FEM::EIGENVAL_FILENAME, "w"), std::fclose);
        
        for (auto i = 0; i < NODE_TOTAL - 1; i++) {
            // 厳密な結果と比較
            std::fprintf(fp_eigenval.get(), "%d, %.14f, %.14f\n", i + 1, eigenval_[i], - 0.5 / static_cast<double>((i + 1) * (i + 1)));
        }

    }
        
    // #endregion publicメンバ関数

    // #region privateメンバ関数

    void Hydrogen_FEM::boundary_conditions()
    {
        // 左辺の全体行列のN + 1行とN + 1列を削る
        hg_.conservativeResize(hg_.rows() - 1, hg_.cols() - 1);

        // 右辺の全体行列のN + 1行とN + 1列を削る
        ug_.conservativeResize(ug_.rows() - 1, ug_.cols() - 1);
    }

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

    void Hydrogen_FEM::make_data()
    {
        // Global節点のx座標を定義(R_MIN～R_MAX）
        auto const dr = (R_MAX - R_MIN) / static_cast<double>(ELE_TOTAL);
        for (auto i = 0; i <= ELE_TOTAL; i++) {
            // 計算領域を等分割
            node_r_glo_[i] = R_MIN + static_cast<double>(i) * dr;
        }

        for (auto e = 0; e < ELE_TOTAL; e++) {
            node_num_seg_[e][0] = e;
            node_num_seg_[e][1] = e + 1;
        }
        
        for (auto e = 0; e < ELE_TOTAL; e++) {
            for (auto i = 0; i < 2; i++) {
                node_r_ele_[e][i] = node_r_glo_[node_num_seg_[e][i]];
            }
        }
    }
    
    void Hydrogen_FEM::make_global_matrix()
    {
        for (auto e = 0; e < ELE_TOTAL; e++) {
            for (auto i = 0; i < 2; i++) {
                for (auto j = 0; j < 2; j++) {
                    hg_(node_num_seg_[e][i], node_num_seg_[e][j]) += mat_A_ele_[e][i][j];
                    ug_(node_num_seg_[e][i], node_num_seg_[e][j]) += mat_B_ele_[e][i][j];
                }
            }
        }
    }

    void Hydrogen_FEM::normalize()
    {
        auto sum = 0.0;
        auto const size = phi_.size();
        auto const max = size - 2;

        // Simpsonの公式によって数値積分する
        for (auto i = 0; i < max; i += 2) {
            auto const f0 = phi_[i] * phi_[i] * node_r_glo_[i] * node_r_glo_[i];
            auto const f1 = phi_[i + 1] * phi_[i + 1] * node_r_glo_[i + 1] * node_r_glo_[i + 1];
            auto const f2 = phi_[i + 2] * phi_[i + 2] * node_r_glo_[i + 2] * node_r_glo_[i + 2];
            sum += (f0 + 4.0 * f1 + f2);
        }
        
        auto const a_1 = 1.0 / std::sqrt(sum * length_[0] / 3.0);

        for (auto i = 0; i < size; i++) {
            phi_[i] *= -a_1;
        }
    }

    // #endregion privateメンバ関数
}
