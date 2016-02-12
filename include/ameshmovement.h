#ifndef __AMESHMOVEMENT_H
#define __AMESHMOVEMENT_H 1

namespace amc {

/// Abstract class for mesh movement
class MeshMove
{
public:
	virtual void move() = 0;
};

} // end namespace

#endif