#ifndef figShapes_H
#define figShapes_H

#include <scalar.H>
#include <label.H>
#include <Switch.H>
#include <word.H>
#include <tmp.H>

namespace Foam
{

class Shape : public refCount
{
public:
    Shape();
    Shape(const Shape&);
    virtual ~Shape();
    virtual scalar area() const;
    virtual scalar centerX() const = 0;
    virtual scalar centerY() const = 0;
};

class Circle : public Shape
{
private:
    scalar x_;
    scalar y_;
    scalar r_;
public:
    Circle();
    Circle(scalar rx, scalar ry, scalar r);
    Circle(const Circle&);
    virtual ~Circle();
    virtual scalar area() const;
    virtual scalar centerX() const;
    virtual scalar centerY() const;
    virtual void   centerX(scalar x);
    virtual void   centerY(scalar y);
};

}

#endif

