
var fragmentShaderCode3 = `#version 300 es

#ifdef GL_FRAGMENT_PRECISION_HIGH
	precision highp float;
#else
	precision mediump float;
#endif

uniform float width;
uniform float height;
uniform float time;

out vec4 FragmentColor;

float sdSphere( vec3 p, float s )
{
    return length(p)-s;
}

vec2 opU( vec2 d1, vec2 d2 )
{
	return (d1.x<d2.x) ? d1 : d2;
}

vec2 map( in vec3 pos )
{
    vec2 res = vec2( 1e10, 0.0 );

    res = opU( res, vec2( sdSphere(    pos-vec3(-2.0,0.25, 0.0), 0.25 ), 26.9 ) );

    return res;
}

vec2 iBox( in vec3 ro, in vec3 rd, in vec3 rad ) 
{
    vec3 m = 1.0/rd;
    vec3 n = m*ro;
    vec3 k = abs(m)*rad;
    vec3 t1 = -n - k;
    vec3 t2 = -n + k;
	return vec2( max( max( t1.x, t1.y ), t1.z ),
	             min( min( t2.x, t2.y ), t2.z ) );
}

vec2 raycast( in vec3 ro, in vec3 rd )
{
    vec2 res = vec2(-1.0,-1.0);

    float tmin = 1.0;
    float tmax = 20.0;

    float tp1 = (0.0-ro.y)/rd.y;
    if( tp1>0.0 )
    {
        tmax = min( tmax, tp1 );
        res = vec2( tp1, 1.0 );
    }
    
    vec2 tb = iBox( ro-vec3(0.0,0.4,-0.5), rd, vec3(2.5,0.41,3.0) );
    if( tb.x<tb.y && tb.y>0.0 && tb.x<tmax)
    {
        tmin = max(tb.x,tmin);
        tmax = min(tb.y,tmax);

        float t = tmin;
        for( int i=0; i<70; i++ )
        {
            vec2 h = map( ro+rd*t );
            if( abs(h.x)<(0.0001*t) )
            { 
                res = vec2(t,h.y); 
                break;
            }
            t += h.x;
			if(t >= tmax)
			{
				break;
			}
        }
    }
    
    return res;
}

float calcSoftshadow( in vec3 ro, in vec3 rd, in float mint, in float tmax )
{
    // bounding volume
    float tp = (0.8-ro.y)/rd.y; if( tp>0.0 ) tmax = min( tmax, tp );

    float res = 1.0;
    float t = mint;
    for( int i=0; i<24; i++ )
    {
		float h = map( ro + rd*t ).x;
        float s = clamp(8.0*h/t,0.0,1.0);
        res = min( res, s*s*(3.0-2.0*s) );
        t += clamp( h, 0.02, 0.2 );
        if( res<0.004 || t>tmax ) break;
    }
    return clamp( res, 0.0, 1.0 );
}

vec3 calcNormal( in vec3 pos )
{
    vec3 n = vec3(0.0);
    for( int i=0; i<4; i++ )
    {
        vec3 e = 0.5773*(2.0*vec3((((i+3)>>1)&1),((i>>1)&1),(i&1))-1.0);
        n += e*map(pos+0.0005*e).x;
    }
    return normalize(n);
}

// float calcAO( in vec3 pos, in vec3 nor )
// {
// 	float occ = 0.0;
//     float sca = 1.0;
//     for( int i=ZERO; i<5; i++ )
//     {
//         float h = 0.01 + 0.12*float(i)/4.0;
//         float d = map( pos + h*nor ).x;
//         occ += (h-d)*sca;
//         sca *= 0.95;
//         if( occ>0.35 ) break;
//     }
//     return clamp( 1.0 - 3.0*occ, 0.0, 1.0 ) * (0.5+0.5*nor.y);
// }

// float checkersGradBox( in vec2 p, in vec2 dpdx, in vec2 dpdy )
// {
//     vec2 w = abs(dpdx)+abs(dpdy) + 0.001;
//     vec2 i = 2.0*(abs(fract((p-0.5*w)*0.5)-0.5)-abs(fract((p+0.5*w)*0.5)-0.5))/w;
//     return 0.5 - 0.5*i.x*i.y;                  
// }

// vec3 render( in vec3 ro, in vec3 rd, in vec3 rdx, in vec3 rdy )
// { 
//     // background
//     vec3 col = vec3(0.7, 0.7, 0.9) - max(rd.y,0.0)*0.3;
    
//     // raycast scene
//     vec2 res = raycast(ro,rd);
//     float t = res.x;
// 	float m = res.y;
//     if( m>-0.5 )
//     {
//         vec3 pos = ro + t*rd;
//         vec3 nor = (m<1.5) ? vec3(0.0,1.0,0.0) : calcNormal( pos );
//         vec3 ref = reflect( rd, nor );
        
//         // material        
//         col = 0.2 + 0.2*sin( m*2.0 + vec3(0.0,1.0,2.0) );
//         float ks = 1.0;
        
//         if( m<1.5 )
//         {
//             // project pixel footprint into the plane
//             vec3 dpdx = ro.y*(rd/rd.y-rdx/rdx.y);
//             vec3 dpdy = ro.y*(rd/rd.y-rdy/rdy.y);

//             float f = checkersGradBox( 3.0*pos.xz, 3.0*dpdx.xz, 3.0*dpdy.xz );
//             col = 0.15 + f*vec3(0.05);
//             ks = 0.4;
//         }

//         // lighting
//         float occ = calcAO( pos, nor );
        
// 		vec3 lin = vec3(0.0);

//         // sun
//         {
//             vec3  lig = normalize( vec3(-0.5, 0.4, -0.6) );
//             vec3  hal = normalize( lig-rd );
//             float dif = clamp( dot( nor, lig ), 0.0, 1.0 );
//           //if( dif>0.0001 )
//         	      dif *= calcSoftshadow( pos, lig, 0.02, 2.5 );
// 			float spe = pow( clamp( dot( nor, hal ), 0.0, 1.0 ),16.0);
//                   spe *= dif;
//                   spe *= 0.04+0.96*pow(clamp(1.0-dot(hal,lig),0.0,1.0),5.0);
//                 //spe *= 0.04+0.96*pow(clamp(1.0-sqrt(0.5*(1.0-dot(rd,lig))),0.0,1.0),5.0);
//             lin += col*2.20*dif*vec3(1.30,1.00,0.70);
//             lin +=     5.00*spe*vec3(1.30,1.00,0.70)*ks;
//         }
//         // sky
//         {
//             float dif = sqrt(clamp( 0.5+0.5*nor.y, 0.0, 1.0 ));
//                   dif *= occ;
//             float spe = smoothstep( -0.2, 0.2, ref.y );
//                   spe *= dif;
//                   spe *= 0.04+0.96*pow(clamp(1.0+dot(nor,rd),0.0,1.0), 5.0 );
//           //if( spe>0.001 )
//                   spe *= calcSoftshadow( pos, ref, 0.02, 2.5 );
//             lin += col*0.60*dif*vec3(0.40,0.60,1.15);
//             lin +=     2.00*spe*vec3(0.40,0.60,1.30)*ks;
//         }
//         // back
//         {
//         	float dif = clamp( dot( nor, normalize(vec3(0.5,0.0,0.6))), 0.0, 1.0 )*clamp( 1.0-pos.y,0.0,1.0);
//                   dif *= occ;
//         	lin += col*0.55*dif*vec3(0.25,0.25,0.25);
//         }
//         // sss
//         {
//             float dif = pow(clamp(1.0+dot(nor,rd),0.0,1.0),2.0);
//                   dif *= occ;
//         	lin += col*0.25*dif*vec3(1.00,1.00,1.00);
//         }
        
// 		col = lin;

//         col = mix( col, vec3(0.7,0.7,0.9), 1.0-exp( -0.0001*t*t*t ) );
//     }

// 	return vec3( clamp(col,0.0,1.0) );
// }

// mat3 setCamera( in vec3 ro, in vec3 ta, float cr )
// {
// 	vec3 cw = normalize(ta-ro);
// 	vec3 cp = vec3(sin(cr), cos(cr),0.0);
// 	vec3 cu = normalize( cross(cw,cp) );
// 	vec3 cv =          ( cross(cu,cw) );
//     return mat3( cu, cv, cw );
// }

void main(void)
{
    // vec2 iResolution = vec2(width, height);
    // vec3 ta = vec3( 0.5, -0.5, -0.6 );
    // vec3 ro = ta + vec3( 4.5*cos(0.1*time), 1.3, 4.5*sin(0.1*time) );
    // mat3 ca = setCamera( ro, ta, 0.0 );

    // vec3 tot = vec3(0.0);
    // vec2 p = (2.0*gl_FragCoord-iResolution.xy)/iResolution.y;

    // const float fl = 2.5;
    
    // vec3 rd = ca * normalize( vec3(p,fl) );

    // vec2 px = (2.0*(gl_FragCoord+vec2(1.0,0.0))-iResolution.xy)/iResolution.y;
    // vec2 py = (2.0*(gl_FragCoord+vec2(0.0,1.0))-iResolution.xy)/iResolution.y;
    // vec3 rdx = ca * normalize( vec3(px,fl) );
    // vec3 rdy = ca * normalize( vec3(py,fl) );

    // vec3 col = render( ro, rd, rdx, rdy );

    // col = pow( col, vec3(0.4545) );

    // tot += col;
    
    // FragmentColor = vec4( tot, 1.0 );
	FragmentColor = vec4(sin(time / 1000.0) * gl_FragCoord.x / width, gl_FragCoord.y / height, gl_FragCoord.z, 1.0);
}

`;
