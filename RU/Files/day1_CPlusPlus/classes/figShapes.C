#include "figShapes.H"

Foam::Shape::Shape()
: refCount()
{
}

Foam::Shape::Shape (const Shape& s)
: refCount()
{
}

Foam::Shape::~Shape()
{
}

Foam::scalar Foam::Shape::area() const
{
    return 0.0;
}

Foam::Circle::Circle()
:
Shape(),
x_(0.0),
y_(0.0),
r_(0.0)
{
}

Foam::Circle::Circle(scalar x, scalar y, scalar r)
:
Shape(),
x_(x),
y_(y),
r_(r)
{
}

Foam::Circle::Circle(const Circle& c)
:
Shape (c),
x_(c.x_),
y_(c.y_),
r_(c.r_)
{
}

Foam::Circle::~Circle()
{
}

Foam::scalar Foam::Circle::area() const
{
    return M_PI*r_*r_;
}

Foam::scalar Foam::Circle::centerX() const
{
    return x_;
}

Foam::scalar Foam::Circle::centerY() const
{
    return y_;
}

void Foam::Circle::centerX(scalar x)
{
    x_ = x;
}

void Foam::Circle::centerY(scalar y)
{
    y_ = y;
}
