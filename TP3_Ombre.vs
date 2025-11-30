const float M_PI = 3.14159265358979323846;

// Rotation matrix around z axis
// a : Angle
mat3 Rz(float a)
{
  float sa=sin(a);float ca=cos(a);
  return mat3(ca,sa,0.,-sa,ca,0.,0.,0.,1.);
}

struct Ray {
vec3 ro; // Ray origin
vec3 rd; // Direction
};

// Compute a point on the ray
// ray : Ray
// t   : depth
vec3 Point(Ray ray,float t)
{
  return ray.ro+t*ray.rd;
}

// Compute the ray
//     m : Mouse position
//     p : Pixel
Ray CreateRay(vec2 m,vec2 p)
{
  float a=3.*3.14*m.x; 
  float le=3.5;
  
  // Origin
  vec3 ro=vec3(45.,0.,10.)*Rz(a);
  
  // Target point
  vec3 ta=vec3(0.,0.,3.);
  
  // Orthonormal frame
  vec3 w=normalize(ta-ro);
  vec3 u=normalize(cross(w,vec3(0.,0.,1.)));
  vec3 v=normalize(cross(u,w));
  vec3 rd=normalize(p.x*u+p.y*v+le*w);
  return Ray(ro,rd);
}

// Primitives

// Sphere
// p : point
// c : center of skeleton
// r : radius
float Sphere(vec3 p,vec3 c,float r)
{
  return length(p-c)-r;
}
//Capsule
float capsule (vec3 p,vec3 a,vec3 b,float r){ 
   vec3 u=(b-a)/length (b-a);
   float la= dot(p-a,u);
   if(la  <0.0) 
   {
   return length(p-a)-r;
   }
   float lb= dot (p-b, u);
   if (lb>0.0){return length (p-b)-r;}
   return (sqrt((length (p-a))*(length (p-a))-la *la)-r);
}
//DEMI Plan
float Demiplan (vec3 p, vec3 c, vec3 n) {
  float norme= length(n);
  vec3 normalise=n/norme;
  return dot ((p-c), normalise);
}
//Cube
float Cube(vec3 p, vec3 c, float L) {
    vec3 normales[6] = vec3[6](
        vec3( 1.0, 0.0, 0.0), vec3(-1.0, 0.0, 0.0),
        vec3( 0.0, 1.0, 0.0), vec3( 0.0,-1.0, 0.0),
        vec3( 0.0, 0.0, 1.0), vec3( 0.0, 0.0,-1.0)
    );
    
    vec3 centres[6] = vec3[6](
        c + vec3( L/2.0, 0.0, 0.0), c + vec3(-L/2.0, 0.0, 0.0),
        c + vec3( 0.0, L/2.0, 0.0), c + vec3( 0.0,-L/2.0, 0.0),
        c + vec3( 0.0, 0.0, L/2.0), c + vec3( 0.0, 0.0,-L/2.0)
    );

    float d = -1e9; // très petit au départ
    for (int i = 0; i < 6; i++) {
        d = max(d, Demiplan(p, centres[i], normales[i]));
    }
    return d;
}
// Plane
// p : point
// c : center of skeleton
// n : Normal
float Plane(vec3 p,vec3 c,vec3 n)
{
  return dot(p-c,n);
}

// Operators

// Union
// a,b : field function of left and right sub-trees
float Union(float a,float b)
{
  return min(a,b);
}

// Union, extension to four sub-trees
// a,b,c : field function of left and right sub-trees
float Union(float a,float b,float c,float d)
{
  return min(min(a,b),min(c,d));
}

// Difference
// a,b : field function of left and right sub-trees
float Difference(float a,float b)
{
  return max(a,-b);
}

// Potential field of the object (Du prof)
// p : point
float object(vec3 p)
{

  float v=  Sphere(p,vec3(3.,0.,6.),4.);

  v=Union(v,
    Sphere(p,vec3(0.,2.,3.),3.));
  
  float d=Union(
      Sphere(p,vec3(5.,1.,5.),3.),
    Sphere(p,vec3(-1.,1.,6.),3.),
    Sphere(p,vec3(-1.,2.,2.),2.),
    Sphere(p,vec3(1.,1.,6.),2.));
  v=Difference(v,d);
  
  v=Union(v,
    Plane(p,vec3(0.,0.,-1.),vec3(0.,0.0,1.0))
  );
  
  return v;
}
//My Oject
/*float object(vec3 p)
{

  float v=  Cube(p,vec3(3.,0.,6.),6.);

  v=Union(v,
    Cube(p,vec3(0.,5.,3.),4.));
  
  
  return v;
}*/
// Analysis of the scalar field

const int Steps=200;// Number of steps
const float Epsilon=.01;// Marching epsilon

// Object normal
// p : point
vec3 ObjectNormal(vec3 p)
{
  const float eps=.001;
  vec3 n;
  float v=object(p);
  n.x=object(vec3(p.x+eps,p.y,p.z))-v;
  n.y=object(vec3(p.x,p.y+eps,p.z))-v;
  n.z=object(vec3(p.x,p.y,p.z+eps))-v;
  return normalize(n);
}

// Trace ray using ray marching
// ray : Ray 
// e : Maximum distance
// h : hit
// s : Number of steps
float SphereTrace(Ray ray,float e,out bool h,out int s)
{
  h=false;
  
  // Start at the origin
  float t=0.;
  
  for(int i=0;i<Steps;i++)
  {
    s=i;
    vec3 p=Point(ray,t);
    float v=object(p);
    // Hit object
    if(v<0.)
    {
      h=true;
      break;
    }
    // Move along ray
    t+=max(Epsilon,v);
    // Escape marched too far away
    if(t>e)
    {
      break;
    }
  }
  return t;
}

// Lighting

// Background color
// d : Ray direction
vec3 Background(vec3 d)
{
  return mix(vec3(.45,.55,.99),vec3(.65,.69,.99),d.z*.5+.5);
}
/*vec3 Background(vec3 d)
{
  return mix(vec3(.1,.6,.9),vec3(.10,.69,.9),d.y+1.);
}*/

// Shading and lighting
// Ombre dure + mouvement 
float Shadow(vec3 p,vec3 n)
{
  // Point light
 // const vec3 lp=vec3(10.,10.,30.);
  vec3 lp=vec3(sin(iTime)*10.,cos(iTime)*-20.,30.);
  // Light direction to point light
  p=vec3(p.x+n.x*0.11,p.y+n.y*0.01,p.z+n.z*0.01);
  vec3 l=normalize(lp-p);
  // nouvelle ombre 
  Ray ray= Ray(p,l);
  
  
  // Hit and number of steps
  bool hit;
  int s;
  
   // Trace ray
 float t=SphereTrace(ray,75.,hit,s);
 float shadow =1.;

  if(hit)
  {
    return 0.;
  }else{return 1.;}
}

//Ombre Douce
float SoftSegment(vec3 p,vec3 n)
{
  // Point light
  //const vec3 lp=vec3(5.,10.,30.);
  p = p + n*0.01;
  const int nbligth = 20;
  vec3 l;
  int j = 0;
  
  vec3 a = vec3(-5.,10.,30.);
  vec3 b = vec3(5.,10.,30.);
  vec3 tabLum[nbligth];
  
  int compteur =0;
  float shadow=0.;
  
  //for(float i=a.x; i<b.x; i=i+10./100.)
  //for(float i=a.x; i<b.x; i=i+(length()/float(nbligth)))

  
  for(float i=a.x; i<b.x; i=i+10./float(nbligth))
  {
   
      tabLum[j] = vec3(i, a.y, a.z);
      
      
      l=normalize(tabLum[j]-p);
      
       
      Ray ray = Ray(p, l);
  
      // Hit and number of steps
   bool hit;
   int s;
      // Trace ray
      float t=SphereTrace(ray,75.,hit,s);
 
      if(hit)
      {
        compteur++;
      }
      
       
      j++;
    
  }
  if(compteur==0)
     shadow=1.;
  else
   return shadow = 1.0 - float(compteur) / float(j);}


//Le cube avec l'aide d'un amie de la promo
float SoftSquare(vec3 p, int n)
{
  float nbhit;
 
  vec3 pl = vec3(5.,5.,30.);
  vec3 pointORG = pl;
  Ray ray;
  bool hit;
  int s;
  float e = 100.;
 
  for(int i=0; i<n; i++)
  {
    pl = vec3(pl.x + 1., pl.y, pl.z);
 
    for(int j=0; j<n; j++)
    {
      pl = vec3(pl.x,pl.y+1.,pl.z);
      ray.ro = p + normalize(pl - p)*0.1;  
      ray.rd = normalize(pl-p);
 
      SphereTrace(ray, e, hit, s);
 
      if(!hit)
      {nbhit++;}
    }
    pl = vec3(pl.x, pointORG.y, pl.z); 
  }
 
  float res;
  res = nbhit/float(n);
  return res;
}
//FIBO du cours
vec3 Fibonacci(int i, int n)
{
    float phi = acos(1.0 - 2.0 * float(i) / float(n));
    float theta = M_PI * (1.0 + sqrt(5.0)) * float(i);
    
    float x = cos(theta) * sin(phi);
    float y = sin(theta) * sin(phi);
    float z = cos(phi);
    

    return vec3(x, y, z);
}
//Ombre douce avec fibo
float SoftShadow(vec3 p, vec3 n)
{
    p = p + n * 0.01;

    const int nblight = 5;
    vec3 tabLum[nblight];

    vec3 lp = vec3(5., 10., 30.); 
    float radius = 5.0;               

    float compteur = 0.;

    bool hit;
    int s;

    for(int j = 0; j < nblight; j++)
    {
        vec3 dir = Fibonacci(j, nblight);   
        tabLum[j] = lp + dir * radius;  
    }

    for(int j = 0; j < nblight; j++)
    {
        vec3 l = normalize(tabLum[j] - p);

        Ray ray = Ray(p, l);

        float t = SphereTrace(ray, 75., hit, s);

        if(hit)
            compteur++;
    }

    return 1.0 -compteur / float(nblight);
}
// Shading and lighting
// p : Point
// n : Normal at point
// e : Eye direction

vec3 shade(vec3 p,vec3 n,vec3 e)
{
  // Point light
  const vec3 lp=vec3(5.,10.,30.);
  /*vec3 lp = vec3(
        0.0 * cos(iTime),
        10.0,
        0.0 + 5.0 * sin(iTime)
    );*/
  // Light direction to point light
  vec3 l=normalize(lp-p);
  
  // Ambient color
  vec3 ambient=.2+.2*Background(n);
  
  // Shadow computation
  //float shadow=1.0;
  //float shadow=Shadow(p,n);
  //float shadow=SoftSegment(p,n);
  //float shadow=SoftSquare(p,1);
  float shadow=SoftShadow(p,n);
  //float shadow=SoftShadow(p,n)+ Shadow(p,n);

  // Phong diffuse
  vec3 diffuse=.35*clamp(dot(n,l),0.,1.)*vec3(.0, .6, 0.);
  
  // Specular
  vec3 r=reflect(e,n);
  vec3 specular=.15*pow(clamp(dot(r,l),0.,1.),35.)*vec3(.0,.2,0.);
  vec3 c=ambient+shadow*(diffuse+specular);
  return c;
}

// Image
void mainImage(out vec4 color,in vec2 pxy)  
{
  // Pixel
  vec2 pixel=(-iResolution.xy+2.*pxy)/iResolution.y;

  // Mouse
  vec2 m=iMouse.xy/iResolution.xy;
  
  // Camera
  Ray ray=CreateRay(m,pixel);
  
  
  // Hit and number of steps
  bool hit;
  int s;
  
   // Trace ray
 float t=SphereTrace(ray,75.,hit,s);
  
  // Shade background
  vec3 rgb=Background(ray.rd);
  
  if(hit)
  {
    // Position
    vec3 p=Point(ray,t);
    
    // Compute normal
    vec3 n=ObjectNormal(p);
    
    // Shade object with light
    rgb=shade(p,n,ray.rd);
  }
    
  color=vec4(rgb,1.);
}