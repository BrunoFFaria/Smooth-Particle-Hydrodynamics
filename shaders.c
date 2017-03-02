//
// Created by bruno on 2/14/17.
//

#define STRINGIFY(A) #A

/* vertex shader */
const char * vertex_shader = STRINGIFY(
    uniform float point_radius;
    uniform float point_scale;
    uniform float density_scale;
    uniform float density_offset;
    uniform float max_z;
    uniform float min_z;
    vec3 vertex_position;
    varying vec3 pos_eye;
    void main(){

        pos_eye = vec3(gl_ModelViewMatrix * vec4(gl_Vertex.xyz, 1.0));
        float dist = length(pos_eye);

        gl_PointSize = point_radius * (point_scale / dist);
        gl_TexCoord[0] = gl_MultiTexCoord0;
        gl_Position = gl_ModelViewProjectionMatrix * vec4(gl_Vertex.xyz, 1.0);

        gl_FrontColor = gl_Color;
    }
);
/*
const char * sphere_shader = STRINGIFY(
     uniform vec3 light_dir = vec3(0.577, 0.577, 0.577);

     void main(){
         vec3 N;

         N.xy = gl_TexCoord[0].xy * vec2(2.0, -2.0) + vec2(-1.0, 1.0);
         float mag = dot(N.xy, N.xy);
         if(mag > 1.0f) discard;
         N.z = sqrt(1-mag);
         float diffuse = max(0.0, dot(light_dir, N));
         gl_FragColor = gl_Color * diffuse;
     }
);
*/

const char * sphere_shader = STRINGIFY(
    uniform float point_radius;
    uniform float near;
    uniform float far;
    varying vec3 pos_eye;
    uniform vec3 light_dir = vec3(0.577, 0.577, 0.577);
    void main(){

        vec3 N;
        N.xy = gl_TexCoord[0].xy * vec2(2.0, -2.0) + vec2(-1.0, 1.0);
        float mag = dot(N.xy,N.xy);
        if(mag > 1.0f) discard;
        N.z = sqrt(1-mag);

        // point on surface on sphere in eye space
        vec4 sphere_pos_eye = vec4(pos_eye+N*point_radius, 1.0);
        vec4 clip_space_pos = gl_ProjectionMatrix * sphere_pos_eye;
        float norm_depth = clip_space_pos.z / clip_space_pos.w;

        //gl_FragDepth = 1.0f - (0.5f*(far-near)*norm_depth + 0.5f * (far+near));
        float diffuse = max(0.0, dot(light_dir, N));
        gl_FragColor = gl_Color*diffuse;

    }
);
