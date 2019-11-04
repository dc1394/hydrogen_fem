/*! \file hydrogen_fem_main.cpp
    \brief FEMで水素原子に対するSchrödinger方程式を解く

    Copyright © 2019 @dc1394 All Rights Reserved.
	This software is released under the BSD 2-Clause License.
*/

#include "hydrogen_fem.h"
#include <iostream>             // for std::cout
#include <boost/format.hpp>     // for boost::format

int main()
{
    hydrogen_fem::Hydrogen_FEM hl;
    std::cout << boost::format("計算が終わりました: 基底状態のエネルギー固有値E = %.14f (Hartree)\n") % hl.do_run();
	hl.save_result();

	std::cout << "計算結果を" << hydrogen_fem::Hydrogen_FEM::RESULTFILENAME << "に書き込みました" << std::endl;

	return 0;
}
