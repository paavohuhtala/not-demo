
uniform float time;
uniform vec2 resolution;

in vec2 v_uv;
out vec4 frag;

const float PI = 3.1415926538;

const int MAX_MARCHING_STEPS = 512;
const float MIN_DIST = 0.0;
const float MAX_DIST = 100.0;
const float EPSILON = 0.0001;

mat4 rotationMatrix(vec3 axis, float angle) {
    axis = normalize(axis);
    float s = sin(angle);
    float c = cos(angle);
    float oc = 1.0 - c;
    
    return mat4(oc * axis.x * axis.x + c,           oc * axis.x * axis.y - axis.z * s,  oc * axis.z * axis.x + axis.y * s,  0.0,
                oc * axis.x * axis.y + axis.z * s,  oc * axis.y * axis.y + c,           oc * axis.y * axis.z - axis.x * s,  0.0,
                oc * axis.z * axis.x - axis.y * s,  oc * axis.y * axis.z + axis.x * s,  oc * axis.z * axis.z + c,           0.0,
                0.0,                                0.0,                                0.0,                                1.0);
}

vec3 rotate(vec3 v, vec3 axis, float angle) {
	mat4 m = rotationMatrix(axis, angle);
	return (m * vec4(v, 1.0)).xyz;
}

float opSmoothUnion( float d1, float d2, float k ) {
    float h = clamp( 0.5 + 0.5*(d2-d1)/k, 0.0, 1.0 );
    return mix( d2, d1, h ) - k*h*(1.0-h);
}

vec4 opUnionColor(vec4 obj1, vec4 obj2) {
  if (obj2.w < obj1.w) return obj2;
  return obj1;
}

float opSmoothSubtraction( float d1, float d2, float k ) {
    float h = clamp( 0.5 - 0.5*(d2+d1)/k, 0.0, 1.0 );
    return mix( d2, -d1, h ) + k*h*(1.0-h);
}

float opSmoothIntersection( float d1, float d2, float k ) {
    float h = clamp( 0.5 - 0.5*(d2-d1)/k, 0.0, 1.0 );
    return mix( d2, d1, h ) + k*h*(1.0-h);
}

float sphereSDF(vec3 p, float r) {
    return length(p) - r;
}

float boxSDF( vec3 p, vec3 b )
{
  vec3 q = abs(p) - b;
  return length(max(q,0.0)) + min(max(q.x,max(q.y,q.z)),0.0);
}

float sdRoundCone( vec3 p, float r1, float r2, float h )
{
  // sampling independent computations (only depend on shape)
  float b = (r1-r2)/h;
  float a = sqrt(1.0-b*b);

  // sampling dependant computations
  vec2 q = vec2( length(p.xz), p.y );
  float k = dot(q,vec2(-b,a));
  if( k<0.0 ) return length(q) - r1;
  if( k>a*h ) return length(q-vec2(0.0,h)) - r2;
  return dot(q, vec2(a,b) ) - r1;
}

float sdCutHollowSphere( vec3 p, float r, float h, float t )
{
  // sampling independent computations (only depend on shape)
  float w = sqrt(r*r-h*h);
  
  // sampling dependant computations
  vec2 q = vec2( length(p.xz), p.y );
  return ((h*q.x<w*q.y) ? length(q-vec2(w,h)) : 
                          abs(length(q)-r) ) - t;
}

float planeSDF(vec3 p, vec3 n, float h)
{
  return dot(p,n) + h;
}

vec3 opRep(vec3 p, vec3 c)
{
    vec3 q = mod(p+0.5*c,c)-0.5*c;
    return q;
}

vec3 opRepLim(vec3 p, float c, in vec3 l )
{
    vec3 q = p-c*clamp(round(p/c),-l,l);
    return q;
}

vec3 opCheapBend(vec3 p)
{
    const float k = 10.0; // or some other amount
    float c = cos(k*p.x);
    float s = sin(k*p.x);
    mat2  m = mat2(c,-s,s,c);
    vec3  q = vec3(m*p.xy,p.z);

    return q;
}

vec4 mushroom(vec3 p) {
    // p = opRepLim(p, 4.0, vec3(10.0, 0.0, 1.0));

    p = rotate(p, vec3(0.0, 1.0, 0.0), time);

    vec3 stemP = rotate(p, vec3(0.0, 0.0, 1.0), 0.06);
    float depth = sdRoundCone(stemP, 0.9, 0.7, 2.0);

    vec3 capP = rotate(p + vec3(0.0, -2.0, 0.0), normalize(vec3(1.0, 0.1, 0.0)), PI);

    depth = opSmoothUnion(depth, sdCutHollowSphere(capP, 1.5, 0.1, 0.1 + sin(time) * 0.02), 0.9);
    return vec4(vec3(0.501, 0.164, 0.078), depth);
}

vec4 scene(vec3 p) {
    // p += vec3(sin(time) * 2.0, sin(time * 0.3) * 1.5 + 2.0, 0.0);

    float h = 1.0 * sin(p.x * 0.1) * sin(p.y * 0.2) + 0.7;

    float planeDepth = planeSDF(p, vec3(0.0, 1.0, 0.0), h);
    vec4 c = vec4(vec3(0.2, 0.4, 0.2), planeDepth);

    c = opUnionColor(c, mushroom(p));

    return c;
}

vec3 calcNormal(vec3 p)
{
    const float h = 0.0001; // replace by an appropriate value
    const vec2 k = vec2(1, -1);
    return normalize(
        k.xyy * scene(p + k.xyy*h).w + 
        k.yyx * scene(p + k.yyx*h).w + 
        k.yxy * scene(p + k.yxy*h).w + 
        k.xxx * scene(p + k.xxx*h).w
    );
}

vec4 shortestDistanceToSurface(vec3 eye, vec3 marchingDirection, float start, float end) {
    float depth = start;
    for (int i = 0; i < MAX_MARCHING_STEPS; i++) {
        vec4 distColor = scene(eye + depth * marchingDirection);
        float dist = distColor.w;
        if (dist < EPSILON) {
			return vec4(distColor.rgb, depth);
        }
        depth += dist;
        if (depth >= end) {
            return vec4(0, 0, 0, end);
        }
    }
    return vec4(vec3(0), end);
}

vec3 phongContribForLight(vec3 k_d, vec3 k_s, float alpha, vec3 p, vec3 eye,
                          vec3 lightPos, vec3 lightIntensity) {
    vec3 N = calcNormal(p);

    vec3 L = normalize(lightPos - p);
    vec3 V = normalize(eye - p);
    vec3 R = normalize(reflect(-L, N));
    
    float dotLN = dot(L, N);
    float dotRV = dot(R, V);
    
    if (dotLN < 0.0) {
        // Light not visible from this point on the surface
        return vec3(0.0, 0.0, 0.0);
    } 
    
    if (dotRV < 0.0) {
        // Light reflection in opposite direction as viewer, apply only diffuse
        // component
        return lightIntensity * (k_d * dotLN);
    }
    return lightIntensity * (k_d * dotLN + k_s * pow(dotRV, alpha));
}

float softshadow( in vec3 ro, in vec3 rd, float mint, float maxt, float k )
{
    float res = 1.0;
    float ph = 1e20;
    for( float t=mint; t<maxt; )
    {
        float h = scene(ro + rd*t).w;
        if( h<0.001 )
            return 0.0;
        float y = h*h/(2.0*ph);
        float d = sqrt(h*h-y*y);
        res = min( res, k*d/max(0.0,t-y) );
        ph = h;
        t += h;
    }
    return res;
}
vec3 phongIllumination(vec3 k_a, vec3 k_d, vec3 k_s, float alpha, vec3 p, vec3 eye) {
    vec3 light1Pos = vec3(2, 10, 20);
    vec3 light2Pos = vec3(2, 5, 4);
    vec3 light1Intensity = vec3(0.9);

    vec3 color = clamp(k_a * softshadow(p, light2Pos, 0.011, 10.0, 0.2), k_a * 0.3, vec3(1.0));
    
    /*color += phongContribForLight(
        k_d,
        k_s,
        alpha,
        p,
        eye,
        light1Pos,
        light1Intensity
    );*/

    color += phongContribForLight(
        k_d,
        k_s,
        alpha,
        p,
        eye,
        light2Pos,
        light1Intensity
    );

    return color;
}


// With better camera setup
vec3 rayDirection(float fieldOfView, vec2 size, vec2 fragCoord) {
    vec2 xy = fragCoord - size / 2.0; // translate screenspace to centre
    float theta = radians(fieldOfView) / 2.0;
    float halfHeight = size.y / 2.0;
    float z = halfHeight / tan(theta);
    return normalize(vec3(xy, -z));
}

void main() {
    vec3 dir = rayDirection(65.0, resolution, gl_FragCoord.xy);
    vec3 eye = vec3(0.0, 1.0, 8.0);

    vec3 backColor = mix(
        vec3(0.780, 0.972, 1),
        vec3(0.619, 0.898, 0.882),
        smoothstep(
            0.0,
            0.7,
            v_uv.y
        )
    );

    vec4 distColor = shortestDistanceToSurface(eye, dir, MIN_DIST, MAX_DIST);
    float dist = distColor.w;

    if (dist > MAX_DIST - EPSILON) {
        // Didn't hit anything
        frag = vec4(backColor, 0.0);
		return;
    }

    // The closest point on the surface to the eyepoint along the view ray
    vec3 p = eye + dist * dir;
    
    vec3 K_a = vec3(0.839, 0.980, 1) * 0.1;
    vec3 K_d = distColor.rgb;
    vec3 K_s = vec3(0.1, 0.1, 0.1);
    float shininess = 1.0;
    
    vec3 color = phongIllumination(K_a, K_d, K_s, shininess, p, eye);
    
    frag = vec4(color, 1.0);
}
