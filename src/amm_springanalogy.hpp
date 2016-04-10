/** \brief Spring-analogy based mesh movement
 * \author Aditya Kashi
 */

#ifndef __AMM_SPRINGANALOGY_H

#define __AMM_SPRINGANALOGY_H

#ifndef __AMM_BASE_H
#include <amm_base.hpp>
#endif

template <int ndim> class SpringAnalogyMeshMovement : public MeshMove
{
public:
	
