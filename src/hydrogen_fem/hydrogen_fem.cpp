/*! \file hydrogen_fem.cpp
    \brief FEMで水素原子に対するSchrödinger方程式を解くクラスの実装

	Copyright © 2019 @dc1394 All Rights Reserved.
     This software is released under the BSD 2-Clause License.
*/

#include "hydrogen_fem.h"
#include <cstdio>				// for FILE, std::fclose, std::fopen, std::fprintf 
#include <memory>				// for std::unique_ptr
#include <boost/assert.hpp>		// for BOOST_ASSERT
#include <Eigen/Eigenvalues>	// for Eigen::GeneralizedSelfAdjointEigenSolver

namespace hydrogen_fem {
    // #region コンストラクタ

    Hydrogen_FEM::Hydrogen_FEM()
		:	hg_(Eigen::MatrixXd::Zero(N, N)),
			ug_(Eigen::MatrixXd::Zero(N, N))
    {
    }

    // #endregion コンストラクタ

    // #region publicメンバ関数 

    double Hydrogen_FEM::do_run()
    {
        // 左辺の全体マトリックスを生成
		make_hg_matrix();

		// 右辺の全体マトリックスを生成
		make_ug_matrix();

        // 一般化固有値問題を解く
        Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> es(hg_, ug_);

        // Eを取得
        auto const e = es.eigenvalues()[0];

        // 固有ベクトルを取得
        c_ = es.eigenvectors().col(0);

        // SCF計算が収束しなかった
        return e;
    }

	void Hydrogen_FEM::save_result() const
	{
		std::unique_ptr<FILE, decltype(&std::fclose)> fp(std::fopen(Hydrogen_FEM::RESULTFILENAME, "w"), std::fclose);

    	for (auto i = 0; i < N; i++) {
			auto const r = static_cast<double>(i) * DH;
			fprintf(fp.get(), "%.14f, %.14f\n", r, c_[i]);
    	}
	}
		
    // #endregion publicメンバ関数

    // #region privateメンバ関数

	double Hydrogen_FEM::get_he_matrix_element(std::int32_t n, std::int32_t p, std::int32_t q) const
	{
		auto const nd = static_cast<double>(n);
		switch (p) {
		case 0:
			switch (q) {
			case 0:
				return 0.5 * DH * (nd * nd + nd + 1.0 / 3.0) - DH * DH * (nd / 3.0 + 1.0 / 12.0);

			case 1:
				return -0.5 * DH * (nd * nd + nd + 1.0 / 3.0) - DH * DH * (nd / 6.0 + 1.0 / 12.0);

			default:
				BOOST_ASSERT(!"heの添字が2以上！");
				return 0.0;
			}

		case 1:
			switch (q) {
			case 0:
				return -0.5 * DH * (nd * nd + nd + 1.0 / 3.0) - DH * DH * (nd / 6.0 + 1.0 / 12.0);

			case 1:
				return 0.5 * DH * (nd * nd + nd + 1.0 / 3.0) - DH * DH * (nd / 3.0 + 0.25);

			default:
				BOOST_ASSERT(!"heの添字が2以上！");
				return 0.0;
			}

		default:
			BOOST_ASSERT(!"heの添字が2以上！");
			return 0.0;
		}
	}
	
    double Hydrogen_FEM::get_ue_matrix_element(std::int32_t n, std::int32_t p, std::int32_t q) const
    {
		auto const nd = static_cast<double>(n);
		switch (p) {
		case 0:
			switch (q) {
			case 0:
				return DH * DH * DH * (nd * nd / 3.0 + nd / 6.0 + 1.0 / 30.0);

			case 1:
				return DH * DH * DH * (nd * nd / 6.0 + nd / 6.0 + 1.0 / 20.0);

			default:
				BOOST_ASSERT(!"ueの添字が2以上！");
				return 0.0;
			}

		case 1:
			switch (q) {
			case 0:
				return DH * DH * DH * (nd * nd / 6.0 + nd / 6.0 + 1.0 / 20.0);

			case 1:
				return DH * DH * DH * (nd * nd / 3.0 + nd / 2.0 + 1.0 / 5.0);

			default:
				BOOST_ASSERT(!"ueの添字が2以上！");
				return 0.0;
			}

		default:
			BOOST_ASSERT(!"ueの添字が2以上！");
			return 0.0;
		}
    }

    void Hydrogen_FEM::make_hg_matrix()
    {
		for (auto i = 0; i < N; i++) {
			if (i == 0)
			{
				hg_(0, 0) = get_he_matrix_element(0, 0, 0);

			}
			else if (i == N - 1)
			{
				hg_(N - 1, N - 1) = get_he_matrix_element(N - 1, 1, 1);
			}
			else
			{
				hg_(i, i) = get_he_matrix_element(i - 1, 1, 1) + get_he_matrix_element(i, 0, 0);
			}

			if (i != N - 1)
			{
				hg_(i, i + 1) = get_he_matrix_element(i, 0, 1);
				hg_(i + 1, i) = get_he_matrix_element(i, 0, 1);
			}
		}
    }

    void Hydrogen_FEM::make_ug_matrix()
    {
		for (auto i = 0; i < N; i++) {
			if (i == 0)
			{
				ug_(0, 0) = get_ue_matrix_element(0, 0, 0);
			}
			else if (i == N - 1)
			{
				ug_(N - 1, N - 1) = get_ue_matrix_element(N - 1, 1, 1);
			}
			else
			{
				ug_(i, i) = get_ue_matrix_element(i - 1, 1, 1) + get_ue_matrix_element(i, 0, 0);
			}
			
			if (i != N - 1)
			{
				ug_(i, i + 1) = get_ue_matrix_element(i, 0, 1);
				ug_(i + 1, i) = get_ue_matrix_element(i, 0, 1);
			}
		}
    }

    // #endregion privateメンバ関数
}
