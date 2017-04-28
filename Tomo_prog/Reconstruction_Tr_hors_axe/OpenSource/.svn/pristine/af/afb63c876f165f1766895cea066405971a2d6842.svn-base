BOOL
sign_of(double n)
{
  return (n >= 0.0);
}


int
integerize(BOOL b)
{
  return (int)b;
}


// assuming there is one, return the index of the element of sign
// different to the other
int
singular_index(BOOL tab[3])
{
  int sum = (integerize(tab[0]) + integerize(tab[1]) + integerize(tab[2]));  

  if (sum == 1) // singular is 1
    {
      if (tab[0]) return 0;
      else if (tab[1]) return 1;
      else return 2;
    }
  else // sum == 2 and singular is 0
    {
      if (! tab[0]) return 0;
      else if (! tab[1]) return 1;
      else return 2;
    }
}


template <class T>
BOOL
Triangle<T>::intersects_plane_p(const Plane<T> &pl)
{
  plane[0] = sign_of(this -> A.distance_signed(pl));
  plane[1] = sign_of(this -> B.distance_signed(pl));
  plane[2] = sign_of(this -> C.distance_signed(pl));

  this -> plane_set = TRUE;

  int sum = (integerize(plane[0]) + integerize(plane[1]) + integerize(plane[2]));
  // intersection only occurs when at least one point is not on the
  // same plane side from the other
  return (!(sum == 0 || sum == 3));
}


template <class T>
Segment<T>
Triangle<T>::intersection_plane(const Plane<T> &pl) 
{
  if (!this -> plane_set)
    {
      BOOL inter = intersects_plane_p(pl);
      assert(inter);
    }

  int singular = singular_index(this -> plane);


  // intersect. points
  Point<T> first;  Point<T> last;
  switch(singular)
    {
    case 0: //A
      first = this -> AB.intersection_plane(pl);
      last = this -> AC.intersection_plane(pl);
      break;
    case 1: //B
      first = this -> BC.intersection_plane(pl);
      last = this -> AB.intersection_plane(pl);
      break;
    case 2: //C
      first = this -> AC.intersection_plane(pl);
      last = this -> BC.intersection_plane(pl);
      break;
    default:
      cerr << "undefined error";
    }


  this -> plane_set = FALSE;
  return Segment<T>(first, last);
}


template <class T>
Point<T> 
Triangle<T>::barycenter() const
{
  Point<T> b(0, 0, 0);
  b += A;
  b += B;
  b += C;
  b /= 3;
  return b;
}


template <class T>
void
Triangle<T>::set_barycenter(Point<T> &b) const
{
  b.set(0, 0, 0);
  b += A;
  b += B;
  b += C;
  b /= 3;
}


template <class T>
double 
Triangle<T>::area() const
{
  double a, h;
  a = BC.length();
  h = this -> height();
  return 0.5f * a * h;
}


template <class T>
double 
Triangle<T>::height() const
{
  double a, b, c, exp;
  a = BC.length();
  b = AC.length();
  c = AB.length();
  exp = (a*a + c*c - b*b)/(2.0f * a);
  return sqrt(c*c - exp*exp);
}


template <class T>
ostream& 
operator<<(ostream &o, const Triangle<T> &t) 
{
  o << "A: " << t.A << "B: " << t.B << "C: " << t.C;

  return o;
}
