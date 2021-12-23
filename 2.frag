
// Analytic bilinear surface ray intersection. Also shows how such a surface is contained
// within a tetrahedron defined by the same 4 points as the surface itself.
//
// This can be useful for when rasterising the bounding surface as a pre-step when mixing
// raytracing with rasterisation.
//
// This example also shows the use of ray differentials, using surface differentials
// obtained directly from the intersection point which are then used to transform the ray
// differentials in to texture space.
//

#define SHOW_BOUNDING_POLYTOPE	1
#define SHOW_BOUNDING_WIREFRAME	1


// Deformed cube geometry
vec3 vertices[8];
int indices[6 * 4];

// Ray intersection test with bilinear surface defined by 4 points
bool traceBilinearPatch(vec3 ro, vec3 rd, vec3 pa, vec3 pb, vec3 pc, vec3 pd,
                        inout vec3 outs, inout vec3 outt, inout vec3 outuvi)
{
    vec3 va = pc - pa;
    vec3 vb = (pd - pc) - (pb - pa);
    vec3 vd = ro - pa;
    vec3 vc = pa - pb;

    // Note that the coefficients are in reverse order here. Solving this
    // equation gives 1 / x and avoids a divide-by-zero case for coplanar controlpoints
    // by switching that case to c = 0 instead of a = 0, thus putting the denominator
    // in the solution as -b instead of 0.
    // Thanks to IQ for this trick!
    
    float c = dot(cross(vb, vc), rd);
    float b = dot(cross(va, vc) + cross(vb, vd), rd);
    float a = dot(cross(va, vd), rd);

    float desc = b * b - 4. * a * c;

    if(desc < 0.0)
        return false;

    // Put PA at the origin of the coordinate system

    ro -=pa;

    pc = va;

    pd -= pa;
    pb -= pa;
    pa -= pa;

    float i;
    float u, v;

    // Solve for U at each intersection point, which are two quadratics

    float u0 = (2. * a) / (-b - sqrt(desc));
    float u1 = (2. * a) / (-b + sqrt(desc));

    vec3 pu0 = pb * u0;
    vec3 pu1 = pb * u1;

    // Solve for V at each intersection point, geometrically

    vec3 vv0 = mix(pc, pd, u0) - pu0;
    vec3 m20 = ro - pu0 - vv0 * dot(vv0, ro - pu0) / dot(vv0, vv0);

    vec3 vv1 = mix(pc, pd, u1) - pu1;
    vec3 m21 = ro - pu1 - vv1 * dot(vv1, ro - pu1) / dot(vv1, vv1);


    float v0, v1;

    {
        vec3 n = cross(va + vb * u0, vd + vc * u0);
        vec3 m = cross(n, rd);

        float da = dot(pu0 - ro, m);
        float db = dot(pc + (pd - pc) * u0 - ro, m);

        v0 = da / (da - db);
    }

    {
        vec3 n = cross(va + vb * u1, vd + vc * u1);
        vec3 m = cross(n, rd);

        float da = dot(pu1 - ro, m);
        float db = dot(pc + (pd - pc) * u1 - ro, m);

        v1 = da / (da - db);
    }

    // Solve for the ray intersection distance at each intersection point

    float da20 = dot(ro - pu0, m20);
    float db20 = dot((ro + rd) - pu0, m20);
    float i0 = da20 / (da20 - db20);

    float da21 = dot(ro - pu1, m21);
    float db21 = dot((ro + rd) - pu1, m21);
    float i1 = da21 / (da21 - db21);

    // Resolve which valid intersection point is nearest to ray origin

    if(u0 < 0. || u0 > 1. || i0 < 0. || v0 < 0. || v0 > 1.)
    {
        u = u1;
        v = v1;
        i = i1;
    }
    else if(u1 < 0. || u1 > 1. || i1 < 0. || v1 < 0. || v1 > 1.)
    {
        u = u0;
        v = v0;
        i = i0;
    }
    else
    {
        u = mix(u0, u1, step(i1, i0));
        v = mix(v0, v1, step(i1, i0));
        i = min(i0, i1);
    }

    if(u < 0. || u > 1. || i < 0. || v < 0. || v > 1.)
        return false;

    outuvi = vec3(u, v, i);

    // Surface differentials in worldspace
    outs = pb * u - mix(pc, pd, u);
    outt = pc * v - mix(pb, pd, v);

    return true;
}

// Ray intersection test with tetrahedron
vec2 traceTetrahedron(vec3 ro, vec3 rd, vec3 pa, vec3 pb, vec3 pc, vec3 pd, inout vec3 outn)
{
    vec3 tn0 = cross(pb - pa, pc - pa);
    vec3 tn1 = cross(pb - pc, pd - pc);
    vec3 tn2 = cross(pd - pa, pb - pa);
    vec3 tn3 = cross(pc - pa, pd - pa);

    if(dot(pd - pc, tn0) > 0.0)
    {
        tn0 = -tn0;
        tn1 = -tn1;
        tn2 = -tn2;
        tn3 = -tn3;
    }

    float td0 = dot(rd, tn0);
    float td1 = dot(rd, tn1);
    float td2 = dot(rd, tn2);
    float td3 = dot(rd, tn3);

    float tt0 = dot(pa - ro, tn0) / td0;
    float tt1 = dot(pc - ro, tn1) / td1;
    float tt2 = dot(pa - ro, tn2) / td2;
    float tt3 = dot(pa - ro, tn3) / td3;

    float tmin = -1e4, tmax = +1e4;

    if(td0 > 0.0)
    {
        tmax = min(tmax, tt0);
    }
    else
    {
        if(tt0 > tmin)
        {
            tmin = tt0;
            outn = tn0;
        }
    }

    if(td1 > 0.0)
    {
        tmax = min(tmax, tt1);
    }
    else
    {
        if(tt1 > tmin)
        {
            tmin = tt1;
            outn = tn1;
        }
    }

    if(td2 > 0.0)
    {
        tmax = min(tmax, tt2);
    }
    else
    {
        if(tt2 > tmin)
        {
            tmin = tt2;
            outn = tn2;
        }
    }

    if(td3 > 0.0)
    {
        tmax = min(tmax, tt3);
    }
    else
    {
        if(tt3 > tmin)
        {
            tmin = tt3;
            outn = tn3;
        }
    }

    return vec2(tmin, tmax);
}

// Ray intersection test with deformed cube composed of six bilinear patches
bool traceDeformedCube(vec3 ro, vec3 rd, inout vec3 outs, inout vec3 outt, inout vec3 outuvi)
{
    vec3 closest_uvi = vec3(1e4), closest_s = vec3(0), closest_t = vec3(0);

    vec3 uvi, norm, s, t;

    // Intersection test against sides of deformed cube

    for(int i = 0; i < 6; ++i)
    {
        vec3 pa = vertices[indices[i * 4 + 0]];
        vec3 pb = vertices[indices[i * 4 + 1]];
        vec3 pc = vertices[indices[i * 4 + 2]];
        vec3 pd = vertices[indices[i * 4 + 3]];

        if(traceBilinearPatch(ro, rd, pa, pb, pc, pd, s, t, uvi))
        {
            if(uvi.z > 0.0 && uvi.z < closest_uvi.z)
            {
                closest_uvi = uvi;
                closest_s = s;
                closest_t = t;
            }
        }
    }

    float u = closest_uvi.x;
    float v = closest_uvi.y;
    float i = closest_uvi.z;

    if(u < 0. || u > 1. || i < 0. || v < 0. || v > 1.)
        return false;

    outs = closest_s;
    outt = closest_t;
    outuvi = closest_uvi;

    return true;
}

vec3 closestPointsOnLines(vec3 p0, vec3 v0, vec3 p1, vec3 v1)
{
    return inverse(mat3(v0, -v1, cross(v1, v0))) * (p1 - p0);
}

float lineMask(vec3 ro, vec3 rd, vec3 pa, vec3 pb, float r, float maxt)
{
    vec3 t = closestPointsOnLines(ro, rd, pa, pb - pa);

    vec3 lp = mix(pa, pb, clamp(t.y, 0., 1.));

    return step(distance(ro + rd * t.x, lp), r) * step(t.x, maxt);
}

float tetrahedronWireframeMask(vec3 ro, vec3 rd, vec3 pa, vec3 pb, vec3 pc, vec3 pd, float r, float maxt)
{
    float mask = 0.;

    mask = max(mask, lineMask(ro, rd, pa, pb, r, maxt));
    mask = max(mask, lineMask(ro, rd, pb, pd, r, maxt));
    mask = max(mask, lineMask(ro, rd, pd, pc, r, maxt));
    mask = max(mask, lineMask(ro, rd, pc, pa, r, maxt));
    mask = max(mask, lineMask(ro, rd, pa, pd, r, maxt));
    mask = max(mask, lineMask(ro, rd, pc, pb, r, maxt));

    return mask;
}

float deformedCubeWireframeMask(vec3 ro, vec3 rd, float r, float maxt)
{
    float mask = 0.;

    for(int i = 0; i < 6; ++i)
    {
        vec3 pa = vertices[indices[i * 4 + 0]];
        vec3 pb = vertices[indices[i * 4 + 1]];
        vec3 pc = vertices[indices[i * 4 + 2]];
        vec3 pd = vertices[indices[i * 4 + 3]];

        mask = max(mask, tetrahedronWireframeMask(ro, rd, pa, pb, pc, pd, r, maxt));
    }

    return mask;
}

float traceDeformedCubeBounds(vec3 ro, vec3 rd, inout vec3 outn)
{
    vec3 closest_uvi = vec3(1e4), closest_s = vec3(0), closest_t = vec3(0);

    float min_i = 1e4;

    // Intersection test against tetrahedral bounds of sides of deformed cube

    for(int i = 0; i < 6; ++i)
    {
        vec3 pa = vertices[indices[i * 4 + 0]];
        vec3 pb = vertices[indices[i * 4 + 1]];
        vec3 pc = vertices[indices[i * 4 + 2]];
        vec3 pd = vertices[indices[i * 4 + 3]];

        vec3 n;

        vec2 is = traceTetrahedron(ro, rd, pa, pb, pc, pd, n);

        if(is.x < is.y && is.x < min_i)
        {
            min_i = is.x;
            outn = n;
        }
    }

    return min_i;
}

// Box-filtered grid, from http://iquilezles.org/www/articles/filterableprocedurals/filterableprocedurals.htm
float grid( in vec2 p, in vec2 dpdx, in vec2 dpdy )
{
    const float N = 10.0; // grid ratio
    vec2 w = max(abs(dpdx), abs(dpdy));
    vec2 a = p + 0.5*w;                        
    vec2 b = p - 0.5*w;           
    vec2 i = (floor(a)+min(fract(a)*N,1.0)-
              floor(b)-min(fract(b)*N,1.0))/(N*w);
    return (1.0-i.x)*(1.0-i.y);
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    float time = iTime;
    vec2 uv = fragCoord / iResolution.xy * 2. - 1.;
    uv.x *= iResolution.x / iResolution.y;

    vec3 ro = vec3(0., -.1, 4.);

    // Set up primary ray including differentials

    vec3 rd = normalize(vec3(uv.xy, -1.6));
    vec3 rdx = rd + dFdx(rd);
    vec3 rdy = rd - dFdy(rd);

    // Base vertices and indices of cube

    vertices[0] = vec3(-1, -1, -1);
    vertices[1] = vec3(-1, -1, +1);
    vertices[2] = vec3(-1, +1, -1);
    vertices[3] = vec3(-1, +1, +1);
    vertices[4] = vec3(+1, -1, -1);
    vertices[5] = vec3(+1, -1, +1);
    vertices[6] = vec3(+1, +1, -1);
    vertices[7] = vec3(+1, +1, +1);

    indices[0 * 4 + 0] = 0;
    indices[0 * 4 + 1] = 1;
    indices[0 * 4 + 2] = 2;
    indices[0 * 4 + 3] = 3;

    indices[1 * 4 + 0] = 7;
    indices[1 * 4 + 1] = 6;
    indices[1 * 4 + 2] = 5;
    indices[1 * 4 + 3] = 4;

    indices[2 * 4 + 0] = 0;
    indices[2 * 4 + 1] = 1;
    indices[2 * 4 + 2] = 4;
    indices[2 * 4 + 3] = 5;

    indices[3 * 4 + 0] = 2;
    indices[3 * 4 + 1] = 3;
    indices[3 * 4 + 2] = 6;
    indices[3 * 4 + 3] = 7;

    indices[4 * 4 + 0] = 0;
    indices[4 * 4 + 1] = 2;
    indices[4 * 4 + 2] = 4;
    indices[4 * 4 + 3] = 6;

    indices[5 * 4 + 0] = 1;
    indices[5 * 4 + 1] = 3;
    indices[5 * 4 + 2] = 5;
    indices[5 * 4 + 3] = 7;

    // Vertex deformation, twisting etc.

    for(int i = 0; i < 8; ++i)
    {
        vec3 p = vertices[i];

        float an = cos(time + p.y) / 2.;
        p.xz = mat2(cos(an), sin(an), sin(an), -cos(an)) * p.xz;

        an = cos(time*2.+p.x+2.);
        p.yz = mat2(cos(an),sin(an),sin(an),-cos(an)) * p.yz;

        an = time + 5.;
        p.xz = mat2(cos(an), sin(an), sin(an), -cos(an)) * p.xz;

        an = time / 3.;
        p.yz = mat2(cos(an), sin(an), sin(an), -cos(an)) * p.yz;

        vertices[i] = p;
    }

    // The convex bounds are traced here as a small speedup
    vec3 boundsn;
    float bounds_i = traceDeformedCubeBounds(ro, rd, boundsn);

    bool hit = false;
    vec3 closest_uvi = vec3(1e4), closest_s = vec3(0), closest_t = vec3(0);

    if(bounds_i > 0. && bounds_i < 1e3)
    {   
        hit = traceDeformedCube(ro, rd, closest_s, closest_t, closest_uvi);
    }

    int textureID = 0;
    float shadow = 1.;

    // Light direction
    vec3 l = normalize(vec3(6, 7, 0));

    // Intersection with floor
    float floor_i = (-2. - ro.y) / rd.y;

    if(floor_i > 0.0 && floor_i < closest_uvi.z)
    {
        // Ray hit the floor
        vec3 rp = ro + rd * floor_i;

        closest_s = vec3(1, 0, 0) * .2;
        closest_t = vec3(0, 0, 1) * .2;

        closest_uvi.x = dot(rp, closest_s);
        closest_uvi.y = dot(rp, closest_t);
        closest_uvi.z = floor_i;

        // vec3 dummy_s, dummy_t, dummy_uvi;
        //if(traceDeformedCube(rp, l, dummy_s, dummy_t, dummy_uvi))
        //{
        //   shadow = .5;
        //}

        textureID = 1;
    }
    else if(hit && (closest_uvi.z < floor_i || floor_i < 0.0))
    {
        textureID = 2;
    }

    if(textureID == 0)
    {
        // Background
        fragColor.rgb = vec3(.1);
    }
    else
    {
        float u = closest_uvi.x;
        float v = closest_uvi.y;
        float i = closest_uvi.z;

        vec3 closest_norm = normalize(cross(closest_s, closest_t));

        // Ensure that the normal is forward-facing

        if(dot(rd, closest_norm) > 0.)
            closest_norm = -closest_norm;

        // Use ray differentials to get intersection points for neighbouring pixels
        // and transform them in to texture space for texture sampling.

        vec3 rp = ro + rd * i;
        vec3 rpx = ro + rdx * dot(rp - ro, closest_norm) / dot(rdx, closest_norm);
        vec3 rpy = ro + rdy * dot(rp - ro, closest_norm) / dot(rdy, closest_norm);

        vec2 duvx = vec2(dot(rpx - rp, closest_s), dot(rpx - rp, closest_t));
        vec2 duvy = vec2(dot(rpy - rp, closest_s), dot(rpy - rp, closest_t));

        // Texturing

        vec3 c;

        if(textureID == 1)
        {
            c = textureGrad(iChannel2, vec2(u, v), duvx / 2., duvy / 2.).rgb;
            c *= mix(.5, 1., smoothstep(0., 2., length(rp.xz)));
        }
        else if(textureID == 2)
        {
            c = textureGrad(iChannel0, vec2(u, v), duvx / 2., duvy / 2.).rgb;
            c *= mix(.5, 1., smoothstep(0., .1, 1. - 2. * max(abs(u - .5), abs(v - .5))));
            float gridscale = 8.0;
            vec3 gridc = vec3(grid(vec2(u, v) * gridscale + .05, duvx / 2. * gridscale, duvy / 2. * gridscale));
            c = mix(gridc, c, smoothstep(.2, .26, fract(time / 5. - .7)) - smoothstep(.7, .76, fract(time / 5. - .7)));
        }

        // Shading

        vec3 r = reflect(rd, closest_norm);
        float fres = mix(.01, .8, pow(clamp(1. - dot(-rd, closest_norm), 0., 1.), 2.));

        vec3 spec = texture(iChannel1, r).rgb * .5;

        spec += pow(max(0., dot(closest_norm, normalize(l + normalize(-rd)))), 64.);

        vec3 diff = c * max(0., .5 + .5 * dot(l, closest_norm));

        fragColor.rgb = mix(diff, spec, fres) * shadow;
    }


    #if SHOW_BOUNDING_POLYTOPE
    if(bounds_i > 0. && bounds_i < 1e3)
    {
        boundsn = normalize(boundsn);
        vec3 bounds_col = vec3(.5 + .5 * dot(boundsn, l));
        fragColor.rgb = mix(fragColor.rgb, bounds_col, .5);
    }
    #endif

    #if SHOW_BOUNDING_WIREFRAME
    float wireframe_mask = deformedCubeWireframeMask(ro, rd, 0.01, closest_uvi.z);
    fragColor.rgb = mix(fragColor.rgb, vec3(1), wireframe_mask * .5);
    #endif

    // Gamma correction

    fragColor.rgb = pow(fragColor.rgb, vec3(1. / 2.2));
}
